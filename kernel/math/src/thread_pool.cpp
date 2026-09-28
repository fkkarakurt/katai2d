#include <katai/math/thread_pool.hpp>

#include <atomic>
#include <condition_variable>
#include <exception>
#include <mutex>
#include <thread>
#include <vector>

namespace katai::math {
namespace {

// Are we inside a pool worker? (nested parallel_for → run serially, avoid deadlock)
thread_local bool t_in_pool_worker = false;

class Pool {
public:
    static Pool& instance() {
        static Pool pool;
        return pool;
    }

    int width() const { return nthreads_; }

    // The pool serves ONE job at a time. A second thread asking while a job is running used to
    // write its own job over the first one's -- the function, the range and the count of
    // workers still out -- so the workers ran pieces of the wrong job, the count never came back
    // to zero, and two jobs solved concurrently in one process hung or crashed (measured: 8 copies
    // of test_jobs in parallel, most rounds hung and some segfaulted). Whoever finds it busy runs
    // the whole range on its own thread instead; every caller writes only to per-index storage,
    // so the answer is the same either way.
    bool try_acquire() { return busy_.try_lock(); }
    void release() { busy_.unlock(); }

    // Split [0,n) into nthreads_ fixed contiguous blocks; block 0 on the caller,
    // 1..nthreads_-1 on the workers. On return all blocks have finished.
    void run(int n, const std::function<void(int, int)>& fn) {
        {
            std::lock_guard<std::mutex> lk(m_);
            job_fn_ = &fn;
            job_n_ = n;
            job_error_ = nullptr;
            pending_ = nthreads_ - 1;
            ++generation_;
        }
        cv_work_.notify_all();

        run_chunk(0);  // the calling thread processes block 0

        std::unique_lock<std::mutex> lk(m_);
        cv_done_.wait(lk, [this] { return pending_ == 0; });
        job_fn_ = nullptr;
        if (job_error_) {
            std::exception_ptr e = job_error_;
            job_error_ = nullptr;
            lk.unlock();
            std::rethrow_exception(e);
        }
    }

    void run_chunk(int t) {
        // Fixed contiguous blocking: block boundaries depend only on (n, nthreads_) →
        // deterministic from call to call.
        const int n = job_n_;
        const long long b = (long long)n * t / nthreads_;
        const long long e = (long long)n * (t + 1) / nthreads_;
        if (b >= e) return;
        try {
            (*job_fn_)((int)b, (int)e);
        } catch (...) {
            std::lock_guard<std::mutex> lk(m_);
            if (!job_error_) job_error_ = std::current_exception();
        }
    }

private:
    Pool() {
        const unsigned hc = std::thread::hardware_concurrency();
        nthreads_ = hc > 1 ? (int)hc : 1;
        workers_.reserve(nthreads_ - 1);
        for (int t = 1; t < nthreads_; ++t)
            workers_.emplace_back([this, t] { worker_loop(t); });
    }

    ~Pool() {
        {
            std::lock_guard<std::mutex> lk(m_);
            shutdown_ = true;
            ++generation_;
        }
        cv_work_.notify_all();
        for (auto& w : workers_) w.join();
    }

    void worker_loop(int t) {
        t_in_pool_worker = true;
        unsigned long long seen = 0;
        for (;;) {
            {
                std::unique_lock<std::mutex> lk(m_);
                cv_work_.wait(lk, [&] { return shutdown_ || generation_ != seen; });
                if (shutdown_) return;
                seen = generation_;
            }
            run_chunk(t);
            {
                std::lock_guard<std::mutex> lk(m_);
                if (--pending_ == 0) cv_done_.notify_one();
            }
        }
    }

    int nthreads_ = 1;
    std::mutex busy_;
    std::vector<std::thread> workers_;
    std::mutex m_;
    std::condition_variable cv_work_, cv_done_;
    unsigned long long generation_ = 0;
    bool shutdown_ = false;
    int pending_ = 0;
    const std::function<void(int, int)>* job_fn_ = nullptr;
    int job_n_ = 0;
    std::exception_ptr job_error_ = nullptr;
};

} // namespace

void parallel_for(int n, const std::function<void(int, int)>& fn) {
    if (n <= 0) return;
    // Small job or a call from inside a worker (nested): straight serial — the pool never engages.
    if (n < 64 || t_in_pool_worker || Pool::instance().width() == 1) {
        fn(0, n);
        return;
    }
    Pool& pool = Pool::instance();
    if (!pool.try_acquire()) {
        fn(0, n);
        return;
    }
    struct Release {
        Pool& p;
        ~Release() { p.release(); }
    } release{pool};
    pool.run(n, fn);
}

int parallel_threads() {
    return Pool::instance().width();
}

} // namespace katai::math
