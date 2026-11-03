"""A run can be cancelled from Python while it runs: Job.run() lets go of the interpreter while it
solves, so another thread can ask it to stop, and it stops at its next poll.

The Prandtl strip footing of the corpus, loaded to collapse, runs for most of a minute. It is
started on a worker thread and asked to stop, from the main thread, a second after its loading
phase has begun: the request returns at once (it does not wait for the interpreter the run used to
hold to the end), the run returns within a few seconds of it, cancelled, with the initial phase
kept and the phase under way dropped. The request waits for the loading phase to begin, not for a
fixed time: on a loaded machine the initial phase took longer than a second, and a request that
came in it left no phase to keep.
"""
import os
import sys
import threading
import time

import katai

failures = 0


def check(ok, what):
    global failures
    print(("ok:   " if ok else "FAIL: ") + what)
    if not ok:
        failures += 1


corpus = sys.argv[1]
pr, _ = katai.load_project(os.path.join(corpus, "kv-fnd-010-prandtl-strip-footing.k2d"))
job = katai.Job(pr)
returned = {}
loading = threading.Event()
job.set_on_phase(lambda i, n, name: loading.set() if i >= 1 else None)


def run():
    job.run()
    returned["t"] = time.perf_counter()


worker = threading.Thread(target=run)
worker.start()
loading.wait(120.0)
time.sleep(1.0)
t_ask = time.perf_counter()
job.request_cancel()
t_asked = time.perf_counter()
worker.join(120.0)
alive = worker.is_alive()
latency = returned.get("t", float("inf")) - t_ask
print("   the request took %.3f s; the run returned %.3f s after it: %s" % (t_asked - t_ask, latency, job.message()))
check(not alive, "the run returned")
check(t_asked - t_ask < 0.5, "asking does not wait for the run")
check(latency < 5.0, "the run stops within seconds of the request")
check(job.state() == katai.JobState.Cancelled and len(job.results()) == 1 and job.results()[0].ok,
      "cancelled, with the initial phase kept and the phase under way dropped")
sys.exit(1 if failures else 0)
