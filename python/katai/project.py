"""The engineer-facing project builder.

Design rule (maintainer directive, 2026-08-02): easy to LEARN, easy to APPLY.
A script should read like the engineering problem -- materials, regions, water,
loads, phases -- with keyword arguments, plain names and defaults that match
what an engineer would assume. Units are FIXED and stated everywhere: kN, m, day.

A tiny slope-stability analysis, end to end::

    import katai

    prj = katai.Project("My slope", mesh_size=3.0)
    soil = prj.materials.mohr_coulomb("Clayey sand", E=1e5, nu=0.3,
                                      c=3.0, phi=19.6, gamma=20.2)
    prj.geometry.polygon([(20, 20), (70, 20), (70, 35), (50, 35),
                          (30, 25), (20, 25)], material=soil,
                         fix=["full", "horizontal", "free", "free",
                              "free", "horizontal"])
    prj.initial(procedure="safety")          # phi-c reduction of the gravity state
    job = prj.run()
    print("FoS", job.results()[-1].fos)

Everything built here becomes a plain schema ``Project`` (``prj.build()``), so a
script, a saved ``.k2d`` and the GUI describe the same run -- the acceptance test
pins that a DSL-built corpus case is BYTE-identical to its checked-in file.
"""

from . import _core

# Plain-name vocabularies. Keys are what an engineer types; values are the schema's.
_FIXITY = {
    "free": _core.BCType.Free,
    "normal": _core.BCType.NormallyFixed,
    "horizontal": _core.BCType.HorizontallyFixed,
    "vertical": _core.BCType.VerticallyFixed,
    "full": _core.BCType.FullyFixed,
}
_DRAINAGE = {
    "drained": _core.Drainage.Drained,
    "undrained_a": _core.Drainage.Undrained,
    "undrained_b": _core.Drainage.UndrainedB,
    "undrained_c": _core.Drainage.UndrainedC,
    "non_porous": _core.Drainage.NonPorous,
}
_PROCEDURE = {
    "k0": _core.InitialProcedure.K0Procedure,
    "gravity": _core.InitialProcedure.GravityLoading,
    "safety": _core.InitialProcedure.Safety,
}
_FLOW = {
    "closed": _core.FlowBCType.Closed,
    "seepage": _core.FlowBCType.Seepage,
}
_WAVE = {
    "harmonic": _core.SeismicWave.Harmonic,
    "ricker": _core.SeismicWave.Ricker,
    "record": _core.SeismicWave.Record,
}
# Rectangle edge order as the geometry helper emits it (counter-clockwise).
_RECT_EDGES = ("bottom", "right", "top", "left")


def _pick(table, value, what):
    try:
        return table[value]
    except KeyError:
        raise ValueError(
            f"unknown {what} '{value}'; choose one of {', '.join(sorted(table))}") from None


class _Handle:
    """An index into one of the project's lists, safe to pass around by name."""

    def __init__(self, index, name):
        self.index = index
        self.name = name

    def __repr__(self):
        return f"<{type(self).__name__} {self.index}: {self.name!r}>"


class MaterialHandle(_Handle):
    pass


class RegionHandle(_Handle):
    pass


class LoadHandle(_Handle):
    pass


class DispHandle(_Handle):
    pass


class _Materials:
    """``prj.materials`` -- one constructor per constitutive model, named as in the
    registry. Common ground parameters: ``gamma`` [kN/m3] above water (``gamma_sat``
    defaults to it), ``k`` [m/day] isotropic permeability (or ``kx``/``ky``),
    ``drainage`` in {"drained", "undrained_a", "undrained_b", "undrained_c",
    "non_porous"} -- "undrained_c" is a TOTAL stress analysis: E, nu and c are the
    undrained values and no pore pressure is generated.

    For an undrained material the pore fluid's stiffness Kw/n follows the equivalent
    undrained Poisson ratio ``nu_u`` (default 0.495, as in PLAXIS), or Skempton's B if
    ``und_mode=1, skempton_B=...`` is given instead. Both are per material."""

    def __init__(self, prj):
        self._prj = prj

    def _add(self, name, model, *, E=None, nu=None, c=None, phi=None, psi=None,
             gamma=None, gamma_sat=None, k=None, kx=None, ky=None,
             drainage="drained", **extra):
        m = _core.Material()
        m.name = name
        m.model = model
        m.drainage = _pick(_DRAINAGE, drainage, "drainage")
        if E is not None: m.E = E
        if nu is not None: m.nu = nu
        if c is not None: m.c = c
        if phi is not None: m.phi = phi
        if psi is not None: m.psi = psi
        if gamma is not None:
            m.gamma_unsat = gamma
            m.gamma_sat = gamma if gamma_sat is None else gamma_sat
        elif gamma_sat is not None:
            m.gamma_sat = gamma_sat
        if k is not None:
            m.kx = k
            m.ky = k
        if kx is not None: m.kx = kx
        if ky is not None: m.ky = ky
        for field, value in extra.items():
            if not hasattr(m, field):
                raise ValueError(f"Material has no parameter '{field}'")
            setattr(m, field, value)
        self._prj._materials.append(m)
        return MaterialHandle(len(self._prj._materials) - 1, name)

    def linear_elastic(self, name, *, E, nu, **kw):
        """Linear elastic. E [kN/m2], nu [-]."""
        return self._add(name, _core.SoilModel.LinearElastic, E=E, nu=nu, **kw)

    def mohr_coulomb(self, name, *, E, nu, c, phi, psi=0.0, **kw):
        """Mohr-Coulomb. E [kN/m2], c [kN/m2], phi/psi [deg]."""
        return self._add(name, _core.SoilModel.MohrCoulomb,
                         E=E, nu=nu, c=c, phi=phi, psi=psi, **kw)

    def hardening_soil(self, name, *, E50, Eoed, Eur, c, phi, psi=0.0, m=0.5, **kw):
        """Hardening Soil. Reference stiffnesses [kN/m2] at p_ref (default 100 kPa)."""
        return self._add(name, _core.SoilModel.HardeningSoil, c=c, phi=phi, psi=psi,
                         E50ref=E50, Eoedref=Eoed, Eurref=Eur, m=m, **kw)

    def hs_small(self, name, *, E50, Eoed, Eur, G0, gamma07, c, phi, psi=0.0, **kw):
        """HS-small: Hardening Soil + small-strain overlay (G0 [kN/m2], gamma_0.7 [-])."""
        return self._add(name, _core.SoilModel.HSsmall, c=c, phi=phi, psi=psi,
                         E50ref=E50, Eoedref=Eoed, Eurref=Eur, G0ref=G0,
                         gamma07=gamma07, **kw)

    def soft_soil(self, name, *, lam_star, kap_star, c, phi, psi=0.0, **kw):
        """Soft Soil. Modified compression/swelling indices lambda*/kappa* [-]."""
        return self._add(name, _core.SoilModel.SoftSoil, c=c, phi=phi, psi=psi,
                         lam_star=lam_star, kap_star=kap_star, **kw)

    def soft_soil_creep(self, name, *, lam_star, kap_star, mu_star, c, phi, psi=0.0, **kw):
        """Soft Soil Creep: Soft Soil + modified creep index mu* [-]."""
        return self._add(name, _core.SoilModel.SoftSoilCreep, c=c, phi=phi, psi=psi,
                         lam_star=lam_star, kap_star=kap_star, mu_star=mu_star, **kw)


class _Geometry:
    """``prj.geometry`` -- soil regions. Fixities use plain names:
    "free", "normal", "horizontal", "vertical", "full"."""

    def __init__(self, prj):
        self._prj = prj

    def polygon(self, points, *, material, name="Soil", fix=None, flow=None, head=None):
        """A soil region from ``points`` [(x, y), ...] in metres, counter-clockwise.

        ``fix``: one fixity name per edge (edge i runs from point i to point i+1).
        ``flow``: optional, one per edge -- "closed"/"seepage", ("head", value_m),
        or ("flux", q_m_per_day) with inflow positive.
        """
        poly = _core.SoilPolygon()
        poly.name = name
        poly.material = material.index
        poly.x = [float(p[0]) for p in points]
        poly.y = [float(p[1]) for p in points]
        if fix is not None:
            if len(fix) != len(points):
                raise ValueError(
                    f"fix has {len(fix)} entries for {len(points)} edges")
            poly.edge_bc = [int(_pick(_FIXITY, f, "fixity")) for f in fix]
        if flow is not None:
            if len(flow) != len(points):
                raise ValueError(
                    f"flow has {len(flow)} entries for {len(points)} edges")
            kinds, heads, fluxes = [], [], []
            any_flux = False
            for f in flow:
                if f in ("head", "flux"):
                    raise ValueError(
                        f'"{f}" needs a value: use ("{f}", '
                        + ("12.0) for the head in metres" if f == "head"
                           else "0.02) for the flux in m/day"))
                if isinstance(f, tuple) and f[0] == "head":
                    kinds.append(int(_core.FlowBCType.Head))
                    heads.append(float(f[1]))
                    fluxes.append(0.0)
                elif isinstance(f, tuple) and f[0] == "flux":
                    kinds.append(int(_core.FlowBCType.Flux))
                    heads.append(0.0)
                    fluxes.append(float(f[1]))
                    any_flux = True
                else:
                    kinds.append(int(_pick(_FLOW, f, "flow boundary")))
                    heads.append(0.0)
                    fluxes.append(0.0)
            poly.edge_flow = kinds
            poly.edge_head = heads
            # Only written when a flux edge is actually used: an all-zero array would
            # say the same thing in more bytes, and a model with no flux stays exactly
            # the file it has always been.
            if any_flux:
                poly.edge_flux = fluxes
        self._prj._polygons.append(poly)
        return RegionHandle(len(self._prj._polygons) - 1, name)

    def rectangle(self, x0, y0, x1, y1, *, material, name="Soil",
                  bottom="full", right="horizontal", top="free", left="horizontal",
                  flow=None):
        """A rectangular region [m] with per-side fixities (engineering defaults:
        fixed base, roller sides, free surface). ``flow``: optional dict with keys
        bottom/right/top/left, values "closed"/"seepage"/("head", value_m)/
        ("flux", q_m_per_day); unmentioned sides are closed."""
        points = [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]
        fix = [bottom, right, top, left]
        flow_list = None
        if flow is not None:
            unknown = set(flow) - set(_RECT_EDGES)
            if unknown:
                raise ValueError(f"unknown side(s) {sorted(unknown)}; use {_RECT_EDGES}")
            flow_list = [flow.get(side, "closed") for side in _RECT_EDGES]
        return self.polygon(points, material=material, name=name, fix=fix,
                            flow=flow_list)


class _Water:
    """``prj.water`` -- the phreatic level. Absent = dry."""

    def __init__(self, prj):
        self._prj = prj

    def table(self, y):
        """A horizontal water table at elevation ``y`` [m], across the model."""
        self._prj._water = ([None, None], [float(y), float(y)])  # spans, filled at build

    def phreatic_line(self, points):
        """A (possibly sloped) phreatic polyline [(x, y), ...] in metres."""
        self._prj._water = ([float(p[0]) for p in points],
                            [float(p[1]) for p in points])


class _Loads:
    """``prj.loads`` -- external loads [kN/m]. Handles feed phase activation."""

    def __init__(self, prj):
        self._prj = prj

    def line_load(self, a, b, *, qx=0.0, qy=0.0, name="Load"):
        """A distributed load from point ``a`` to ``b`` [m], qx/qy [kN/m/m]
        (downward = negative qy)."""
        L = _core.Load()
        L.kind = _core.LoadKind.Distributed
        L.name = name
        L.x1, L.y1, L.x2, L.y2 = a[0], a[1], b[0], b[1]
        L.qx1, L.qy1, L.qx2, L.qy2 = qx, qy, qx, qy
        self._prj._loads.append(L)
        return LoadHandle(len(self._prj._loads) - 1, name)

    def point_load(self, at, *, qx=0.0, qy=0.0, name="Load"):
        """A point load at ``at`` [m] -- in plane strain this is a LINE load per
        metre out of plane, qx/qy [kN/m]. Only the first-point components are
        set; the unused second point keeps its schema defaults, exactly as the
        GUI leaves them."""
        L = _core.Load()
        L.kind = _core.LoadKind.Point
        L.name = name
        L.x1, L.y1 = at[0], at[1]
        L.qx1, L.qy1 = qx, qy
        self._prj._loads.append(L)
        return LoadHandle(len(self._prj._loads) - 1, name)


class _Displacements:
    """``prj.displacements`` -- line prescribed displacements (PLAXIS style)."""

    def __init__(self, prj):
        self._prj = prj

    def line(self, a, b, *, ux=None, uy=None, name="Displacement"):
        """Impose displacement components on every mesh node along the line
        ``a``-``b`` [m]. ``ux``/``uy`` [m]: a number prescribes that component
        (ramped 0 -> value over the phase; 0 makes a rigid support line), None
        leaves it free. The classic use is a displacement-controlled rigid
        footing -- its force is then read from ``result.reaction``. Static
        (plastic) phases only; activate/deactivate by handle like loads."""
        if ux is None and uy is None:
            raise ValueError("prescribe at least one component (ux= and/or uy=)")
        D = _core.PrescribedDisp()
        D.name = name
        D.x1, D.y1, D.x2, D.y2 = a[0], a[1], b[0], b[1]
        D.set_ux = ux is not None
        D.ux = 0.0 if ux is None else float(ux)
        D.set_uy = uy is not None
        D.uy = 0.0 if uy is None else float(uy)
        self._prj._disps.append(D)
        return DispHandle(len(self._prj._disps) - 1, name)


class _PhaseBuilder:
    """One staged-construction phase; created via ``prj.phases.*``."""

    def __init__(self, prj, phase):
        self._prj = prj
        self.phase = phase


class _Phases:
    """``prj.phases`` -- staged construction after the initial phase, in order.

    Activation is INHERITED from the previous phase (PLAXIS-style); pass
    ``activate=[...]`` / ``deactivate=[...]`` with region or load handles to
    change only what the phase changes."""

    def __init__(self, prj):
        self._prj = prj

    def _add(self, name, ptype, *, activate=(), deactivate=(), duration=None,
             steps=None, tolerance=None, load_steps=None, max_iterations=None,
             apply_fraction=None, ignore_undrained=None):
        ph = _core.Phase()
        ph.name = name
        ph.type = ptype
        if duration is not None: ph.duration = duration
        if steps is not None: ph.time_steps = steps
        # apply_fraction is PLAXIS's Sum-Mstage: 0.5 applies half the staged change
        # and leaves the rest for a later phase.
        if apply_fraction is not None: ph.sum_mstage = apply_fraction
        if ignore_undrained is not None: ph.ignore_undrained = ignore_undrained
        # Numerical controls; unset = the program chooses by material class. They are
        # written into the .k2d, so a script that pins them publishes a run someone
        # else can reproduce exactly.
        if tolerance is not None: ph.tolerance = tolerance
        if load_steps is not None: ph.load_steps = load_steps
        if max_iterations is not None: ph.max_iterations = max_iterations
        self._prj._phase_toggles.append((ph, list(activate), list(deactivate)))
        return _PhaseBuilder(self._prj, ph)

    def plastic(self, name, **kw):
        """A plastic (staged construction) phase: fill, excavate, install, load.

        ``apply_fraction=0.5`` applies half of the stage (PLAXIS Sum-Mstage) and
        leaves the rest; ``ignore_undrained=True`` solves undrained materials as
        drained for this phase."""
        return self._add(name, _core.PhaseType.Plastic, **kw)

    def consolidation(self, name, *, duration, steps, **kw):
        """Time-dependent Biot consolidation. duration [day]."""
        return self._add(name, _core.PhaseType.Consolidation,
                         duration=duration, steps=steps, **kw)

    def safety(self, name="Safety", **kw):
        """phi-c reduction of the current state -> factor of safety."""
        return self._add(name, _core.PhaseType.Safety, **kw)

    def dynamic(self, name, *, duration, steps, wave="harmonic", amplitude=1.0,
                frequency=1.0, damping=0.05, rayleigh=None, record=None,
                record_dt=None, compliant_base=False, free_field=False,
                nonlinear=False, **kw):
        """A seismic time-history phase. duration [s]; steps = Newmark steps
        over it (dt = duration/steps).

        wave: "harmonic" (base acceleration A sin(2 pi f t); amplitude [m/s2],
        frequency [Hz]), "ricker" (pulse), or "record" (pass record=[a_g ...]
        [m/s2] sampled at record_dt [s]; amplitude scales it, 1 = as recorded;
        the record is stored IN the project file, so the run is reproducible
        without loose side files). damping is the Rayleigh target ratio at the
        two frequencies rayleigh=(f1, f2) [Hz].

        compliant_base: the bottom absorbs outgoing waves (Lysmer dashpots)
        and the input is applied there as the upward wave. free_field: Lysmer
        free-field lateral boundaries instead of free/roller sides. nonlinear:
        full Newton on the internal force each step (plasticity during
        shaking); opt-in, much slower."""
        b = self._add(name, _core.PhaseType.Dynamic,
                      duration=duration, steps=steps, **kw)
        ph = b.phase
        ph.seismic_wave = _pick(_WAVE, wave, "wave")
        ph.seismic_amp = amplitude
        if wave == "record":
            if record is None or record_dt is None:
                raise ValueError('wave="record" needs record=[...] and record_dt')
            ph.accel_record = [float(a) for a in record]
            ph.record_dt = record_dt
        else:
            ph.seismic_freq = frequency
        ph.damping_ratio = damping
        if rayleigh is not None:
            ph.rayleigh_f1, ph.rayleigh_f2 = rayleigh
        ph.seismic_compliant_base = compliant_base
        ph.seismic_free_field = free_field
        ph.dynamic_nonlinear = nonlinear
        return b


class Project:
    """A KATAI 2D analysis, built the way the problem is described.

    Parameters: ``kind`` "plane_strain" (default) or "axisymmetric";
    ``element`` "tri6" (default) or "tri15"; ``mesh_size`` target element edge [m];
    ``auto_refine`` refines near structures and loads.

    Build with ``prj.materials`` / ``prj.geometry`` / ``prj.water`` /
    ``prj.loads`` / ``prj.phases``, choose ``prj.initial(...)``, then
    ``prj.run()`` (or ``prj.save("case.k2d")`` for the GUI and the CLI --
    the three front ends run the same job).
    """

    def __init__(self, name, *, kind="plane_strain", element="tri6",
                 mesh_size=2.0, auto_refine=True):
        if kind not in ("plane_strain", "axisymmetric"):
            raise ValueError("kind must be 'plane_strain' or 'axisymmetric'")
        if element not in ("tri6", "tri15"):
            raise ValueError("element must be 'tri6' or 'tri15'")
        self.name = name
        self._kind = kind
        self._element = element
        self._mesh_size = mesh_size
        self._auto_refine = auto_refine
        self._materials = []
        self._polygons = []
        self._loads = []
        self._disps = []
        self._water = None
        self._procedure = "k0"
        self._initial_exclude = []
        self._initial_numerics = {}
        self._phase_toggles = []
        self.materials = _Materials(self)
        self.geometry = _Geometry(self)
        self.water = _Water(self)
        self.loads = _Loads(self)
        self.displacements = _Displacements(self)
        self.phases = _Phases(self)

    def initial(self, *, procedure="k0", exclude=(), tolerance=None,
                load_steps=None, max_iterations=None, ignore_undrained=None):
        """How the in-situ state is established: "k0" (geostatic, default),
        "gravity" (switch-on self-weight) or "safety" (single-phase slope FoS).
        ``exclude``: loads or regions NOT present in the initial phase.

        ``tolerance`` / ``load_steps`` / ``max_iterations``: numerical controls
        for this phase; unset = chosen by material class. ``ignore_undrained``:
        establish the state with undrained materials solved as drained."""
        self._procedure = procedure
        self._initial_exclude = list(exclude)
        self._initial_numerics = {"tolerance": tolerance, "load_steps": load_steps,
                                  "max_iterations": max_iterations,
                                  "ignore_undrained": ignore_undrained}
        return self

    # ------------------------------------------------------------------ build --
    def _extent(self):
        xs = [x for p in self._polygons for x in p.x]
        ys = [y for p in self._polygons for y in p.y]
        if not xs:
            raise ValueError("the project has no soil regions; add geometry first")
        return min(xs), max(xs), min(ys), max(ys)

    def build(self):
        """The plain schema Project this analysis means -- what ``.k2d`` stores and
        every front end runs. Model extents come from the geometry."""
        pr = _core.Project()
        pr.name = self.name
        pr.axisymmetric = self._kind == "axisymmetric"
        pr.initial_procedure = _pick(_PROCEDURE, self._procedure, "initial procedure")
        mesh = _core.MeshSettings()
        mesh.elem_size = self._mesh_size
        mesh.order = 6 if self._element == "tri6" else 15
        mesh.auto_refine = self._auto_refine
        pr.mesh = mesh
        pr.x_min, pr.x_max, pr.y_min, pr.y_max = self._extent()
        if self._water is not None:
            wx, wy = self._water
            if wx[0] is None:                      # horizontal table: span the model
                wx = [pr.x_min, pr.x_max]
            pr.has_water = True
            pr.wx = wx
            pr.wy = wy
        else:
            pr.has_water = False
        pr.materials = self._materials
        pr.polygons = self._polygons
        pr.loads = self._loads
        pr.disps = self._disps

        # The initial phase: everything active except the exclusions. Vectors are
        # materialized ONLY for a touched class, exactly like the GUI's files.
        initial = _core.Phase()
        excl_loads = [h.index for h in self._initial_exclude
                      if isinstance(h, LoadHandle)]
        excl_regions = [h.index for h in self._initial_exclude
                        if isinstance(h, RegionHandle)]
        excl_disps = [h.index for h in self._initial_exclude
                      if isinstance(h, DispHandle)]
        load_state = [1] * len(self._loads)
        poly_state = [1] * len(self._polygons)
        disp_state = [1] * len(self._disps)
        if excl_loads:
            for i in excl_loads:
                load_state[i] = 0
            initial.load_active = load_state
        if excl_regions:
            for i in excl_regions:
                poly_state[i] = 0
            initial.poly_active = poly_state
        if excl_disps:
            for i in excl_disps:
                disp_state[i] = 0
            initial.disp_active = disp_state
        for field, value in self._initial_numerics.items():
            if value is not None:
                setattr(initial, field, value)
        pr.initial = initial

        # Staged phases: inherit the previous phase's effective state; write only
        # the vectors a toggle touches.
        phases = []
        for ph, activate, deactivate in self._phase_toggles:
            touched_loads = touched_polys = touched_disps = False
            for h in activate + deactivate:
                on = 1 if h in activate else 0
                if isinstance(h, LoadHandle):
                    load_state[h.index] = on
                    touched_loads = True
                elif isinstance(h, RegionHandle):
                    poly_state[h.index] = on
                    touched_polys = True
                elif isinstance(h, DispHandle):
                    disp_state[h.index] = on
                    touched_disps = True
                else:
                    raise ValueError(f"cannot toggle {h!r}")
            if touched_loads:
                ph.load_active = list(load_state)
            if touched_polys:
                ph.poly_active = list(poly_state)
            if touched_disps:
                ph.disp_active = list(disp_state)
            phases.append(ph)
        pr.phases = phases
        return pr

    def validate(self):
        """The input-contract report for this project (field-path issues)."""
        return _core.validate_project(self.build())

    def save(self, path):
        """Write the analysis as a ``.k2d`` -- the same file the GUI opens and
        ``katai solve`` runs."""
        _core.save_project(self.build(), path)

    def run(self, on_phase=None):
        """Run every phase headlessly; returns the finished Job. An invalid
        project raises :class:`katai.Refusal` with the field-path report."""
        from . import run as _run
        return _run(self.build(), on_phase=on_phase)
