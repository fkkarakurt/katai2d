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
_BARRIER = {
    "permeable": 0,          # the default: the flow net runs through the line
    "impermeable": 1,        # a screen; the two sides' pore pressures are separated
    "semi_permeable": 2,     # passes q = dh / R, R = d/k [day]
}
_CONNECTION = {"hinged": 0, "free": 1}   # embedded-beam top attachment
# Design approaches, by the name the code itself uses. "none" is the default: the
# phase is a characteristic-value analysis and no partial factors are applied.
_DESIGN = {
    "none": _core.DesignApproach.None_,   # the binding spells it None_; the engineer types "none"
    "ec7_da1_c1": _core.DesignApproach.EC7_DA1_C1,
    "ec7_da1_c2": _core.DesignApproach.EC7_DA1_C2,
    "ec7_da2": _core.DesignApproach.EC7_DA2,
    "ec7_da3": _core.DesignApproach.EC7_DA3,
    "tbdy2018_static": _core.DesignApproach.TBDY2018_Static,
    "tbdy2018_seismic": _core.DesignApproach.TBDY2018_Seismic,
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


class StructHandle(_Handle):
    """A structural element (wall, anchor, geogrid, pile, interface)."""


class DispHandle(_Handle):
    pass


class HydroHandle(_Handle):
    pass


class _Materials:
    """``prj.materials`` -- one constructor per constitutive model, named as in the
    registry. Common ground parameters: ``gamma`` [kN/m3] above water (``gamma_sat``
    defaults to it), ``k`` [m/day] isotropic permeability (or ``kx``/``ky``),
    ``drainage`` in {"drained", "undrained_a", "undrained_b", "undrained_c",
    "non_porous"} -- "undrained_c" is a TOTAL stress analysis: E, nu and c are the
    undrained values and no pore pressure is generated.

    For an undrained material the pore fluid's stiffness Kw/n follows the equivalent
    undrained Poisson ratio ``nu_u`` (0.495 by default), or Skempton's B if
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


class _Structures:
    """``prj.structures`` -- retaining walls, anchors, geogrids, piles, interfaces.

    Every element is one line from ``a`` to ``b`` [m] plus the property set it is
    made of, created in a single call, and every one returns a handle that phases
    activate and deactivate exactly like a load or a region -- which is how a wall
    gets installed before the excavation that needs it::

        wall = prj.structures.plate((0, 0), (0, -12), EA=1.2e7, EI=1.0e5)
        prj.phases.plastic("Install wall", activate=[wall])
        prj.phases.plastic("Excavate", deactivate=[pit])

    Stiffnesses are PER METRE out of plane, as everywhere else in a plane-strain
    analysis. Anchors and piles are the exception and say so: they are discrete
    members, so they carry a ``spacing`` [m] and the engine divides by it.
    """

    def __init__(self, prj):
        self._prj = prj

    def _add(self, kind, a, b, name, material_index, *, barrier=None,
             resistance=0.0, interfaces=None, coarseness=1.0, conn=None):
        S = _core.StructElement()
        S.kind = kind
        S.name = name
        S.x1, S.y1, S.x2, S.y2 = a[0], a[1], b[0], b[1]
        S.material = material_index
        S.coarseness = coarseness
        if interfaces in ("both", "positive"): S.iface_pos = True
        if interfaces in ("both", "negative"): S.iface_neg = True
        if interfaces not in (None, "none", "both", "positive", "negative"):
            raise ValueError("interfaces must be 'positive', 'negative', 'both' or 'none'")
        if barrier is not None:
            S.flow_barrier = _pick(_BARRIER, barrier, "flow barrier")
            S.hydraulic_resistance = resistance
            if barrier == "semi_permeable" and not resistance > 0.0:
                raise ValueError("a semi-permeable line needs resistance=d/k [day] > 0")
        if conn is not None:
            S.conn = _pick(_CONNECTION, conn, "connection")
        self._prj._structs.append(S)
        return StructHandle(len(self._prj._structs) - 1, name)

    def plate(self, a, b, *, EA, EI, w=0.0, nu=0.0, Mp=0.0, Np=0.0,
              elastoplastic=None, name="Plate", interfaces=None, barrier=None,
              resistance=0.0, coarseness=1.0):
        """A wall, lining or slab: EA [kN/m], EI [kN m2/m], self-weight w [kN/m/m].

        ``Mp`` / ``Np`` (0 = elastic) make it elastoplastic -- the plastic moment
        and axial force it cannot exceed. ``interfaces='both'`` puts a soil-
        structure interface on each side, which is what makes a wall able to slip
        rather than glue the soil to it. ``barrier='impermeable'`` makes it a
        groundwater screen ('semi_permeable' needs ``resistance`` = d/k [day])."""
        M = _core.PlateMaterial()
        M.name, M.EA, M.EI, M.w, M.nu = name, EA, EI, w, nu
        M.Mp, M.Np = Mp, Np
        # A capacity the engine never reads is worse than no capacity: the plate's caps are
        # only applied when the elastoplastic flag is set, so giving Mp or Np turns it on.
        # Pass elastoplastic=False explicitly to keep a stated capacity DORMANT.
        M.elastoplastic = (Mp > 0.0 or Np > 0.0) if elastoplastic is None else elastoplastic
        self._prj._plates.append(M)
        return self._add(_core.StructKind.Plate, a, b, name, len(self._prj._plates) - 1,
                         barrier=barrier, resistance=resistance, interfaces=interfaces,
                         coarseness=coarseness)

    def anchor(self, a, b, *, EA, spacing=1.0, prestress=0.0, Fmax_tens=0.0,
               Fmax_comp=0.0, elastoplastic=None, name="Anchor", coarseness=1.0):
        """A node-to-node anchor or strut: EA [kN], out-of-plane ``spacing`` [m].

        ``prestress`` [kN] is the lock-off force, tension-positive: the member is
        tensioned against the wall when it is installed, and that force holds the
        excavation before anything moves. It is an elastic spring from there on,
        so the force follows the wall rather than staying put. ``Fmax_tens`` /
        ``Fmax_comp`` (0 = unlimited) cap it."""
        M = _core.AnchorMaterial()
        M.name, M.EA, M.Lspacing, M.prestress = name, EA, spacing, prestress
        M.Fmax_tens, M.Fmax_comp = Fmax_tens, Fmax_comp
        M.elastoplastic = ((Fmax_tens > 0.0 or Fmax_comp > 0.0)
                           if elastoplastic is None else elastoplastic)
        self._prj._anchors.append(M)
        return self._add(_core.StructKind.Anchor, a, b, name,
                         len(self._prj._anchors) - 1, coarseness=coarseness)

    def geogrid(self, a, b, *, EA, Np=0.0, elastoplastic=None, name="Geogrid",
                coarseness=1.0):
        """A tension-only reinforcement sheet: EA [kN/m], tension cap ``Np``
        [kN/m] (0 = unlimited). It carries no compression and no bending -- the
        force goes to zero rather than negative, which is the whole point of it.

        Giving ``Np`` ARMS the cap (the engine applies it only when the element is
        marked elastoplastic, so a bare cap would silently do nothing); pass
        ``elastoplastic=False`` to state a capacity and deliberately leave it off."""
        M = _core.GeogridMaterial()
        M.name, M.EA, M.Np = name, EA, Np
        # The engine reads the cap only when the flag is set (driver: Np applies if
        # elastoplastic && Np > 0), so a bare Np= would otherwise be a silent no-op.
        M.elastoplastic = (Np > 0.0) if elastoplastic is None else elastoplastic
        self._prj._geogrids.append(M)
        return self._add(_core.StructKind.Geogrid, a, b, name,
                         len(self._prj._geogrids) - 1, coarseness=coarseness)

    def pile(self, a, b, *, E, diameter, spacing, gamma=24.0, Tskin_max=0.0,
             Fmax_base=0.0, connection="hinged", name="Pile", coarseness=1.0):
        """An embedded beam (pile row): E [kN/m2], ``diameter`` [m], out-of-plane
        ``spacing`` [m], unit weight ``gamma`` [kN/m3].

        ``Tskin_max`` [kN/m] and ``Fmax_base`` [kN] are the shaft and base
        capacities -- CAREFUL, 0 means UNLIMITED here, not zero capacity.
        ``connection`` is how the pile top attaches to the soil: 'hinged' (the
        default) means the top moves with the soil, which is what
        makes the pile loadable at its head at all; 'free' leaves it coupled only
        through the shaft interface, as a grout body is."""
        M = _core.EmbeddedBeamMaterial()
        M.name, M.E, M.gamma, M.diameter, M.Lspacing = name, E, gamma, diameter, spacing
        M.Tskin_max, M.Fmax_base = Tskin_max, Fmax_base
        self._prj._embedded.append(M)
        return self._add(_core.StructKind.EmbeddedBeam, a, b, name,
                         len(self._prj._embedded) - 1, conn=connection,
                         coarseness=coarseness)

    def interface(self, a, b, *, name="Interface", barrier=None, resistance=0.0,
                  coarseness=1.0):
        """A bare soil-structure interface line: it takes the strength of the soil
        beside it, reduced by that soil's ``Rinter``. Use this for a slip surface
        that has no structure on it; a wall's own interfaces come from
        ``plate(..., interfaces='both')``."""
        return self._add(_core.StructKind.Interface, a, b, name, -1,
                         barrier=barrier, resistance=resistance, coarseness=coarseness)


class _Displacements:
    """``prj.displacements`` -- prescribed displacements along a line."""

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


class _Dewatering:
    """``prj.dewatering`` -- wells and drains drawn inside the model, the
    hydraulic conditions. Both are lines, both switch on and off per phase
    like loads, and they say two different things to the ground: a well
    prescribes how much water is taken out, a drain prescribes how low the
    head is held."""

    def __init__(self, prj):
        self._prj = prj

    def well(self, a, b, *, q, h_min, name="Well", infiltration=False):
        """Pump ``q`` [m3/day per m out of plane] out along the line ``a``-``b``
        [m] (``infiltration=True`` puts water in instead). ``h_min`` [m] is the
        lowest head the well can draw the ground down to: once it is reached the
        well extracts only what the ground can supply, and the run says so."""
        H = _core.HydroLine()
        H.name = name
        H.kind = _core.HydroKind.Well
        H.behaviour = 1 if infiltration else 0
        H.x1, H.y1, H.x2, H.y2 = a[0], a[1], b[0], b[1]
        H.q = float(q)
        H.h_min = float(h_min)
        self._prj._hydros.append(H)
        return HydroHandle(len(self._prj._hydros) - 1, name)

    def drain(self, a, b, *, head, name="Drain", vacuum=False):
        """Hold the head at ``head`` [m] along the line ``a``-``b`` [m]. A normal
        drain only takes water away -- ground already drier than the drain is
        left alone; ``vacuum=True`` holds the head in both directions (vacuum
        consolidation). In a consolidation phase a drain sets the excess pore
        pressure to zero."""
        H = _core.HydroLine()
        H.name = name
        H.kind = _core.HydroKind.Drain
        H.behaviour = 1 if vacuum else 0
        H.x1, H.y1, H.x2, H.y2 = a[0], a[1], b[0], b[1]
        H.head = float(head)
        self._prj._hydros.append(H)
        return HydroHandle(len(self._prj._hydros) - 1, name)


class _PhaseBuilder:
    """One staged-construction phase; created via ``prj.phases.*``."""

    def __init__(self, prj, phase):
        self._prj = prj
        self.phase = phase


class _Phases:
    """``prj.phases`` -- staged construction after the initial phase, in order.

    Activation is INHERITED from the previous phase; pass
    ``activate=[...]`` / ``deactivate=[...]`` with region or load handles to
    change only what the phase changes."""

    def __init__(self, prj):
        self._prj = prj

    def _add(self, name, ptype, *, activate=(), deactivate=(), duration=None,
             steps=None, tolerance=None, load_steps=None, max_iterations=None,
             apply_fraction=None, ignore_undrained=None, reset_small_strain=None,
             design=None, water=None):
        ph = _core.Phase()
        ph.name = name
        ph.type = ptype
        if duration is not None: ph.duration = duration
        if steps is not None: ph.time_steps = steps
        # apply_fraction is the staged-change fraction (Sum-Mstage): 0.5 applies half
        # the staged change and leaves the rest for a later phase.
        if apply_fraction is not None: ph.sum_mstage = apply_fraction
        if ignore_undrained is not None: ph.ignore_undrained = ignore_undrained
        # Reset small strain: start this phase with the HS-small strain history
        # cleared, so the soil meets it at G0. Use it when the phases before this one built a
        # state rather than continued a loading path -- a surcharge placed and removed to leave
        # an overconsolidation behind leaves a strain history too, and that history is an
        # artefact of the modelling, not something the soil would still remember.
        if reset_small_strain is not None: ph.reset_small_strain = reset_small_strain
        # Design code: the partial factors of EC7 or TBDY 2018 applied to THIS phase, so a
        # characteristic run and its design check live in one file and are read side by side.
        # Unset = characteristic values, which is what an unfactored analysis means.
        if design is not None: ph.design_approach = _pick(_DESIGN, design, "design approach")
        # A phase-level water table: staged dewatering without touching the project's own.
        # water=y lowers/raises a horizontal table; water=[(x, y), ...] gives a phreatic line.
        # The schema's phase water is a POLYLINE and needs >= 2 points, so a horizontal table
        # is a two-point line spanning the model -- and the model's extent is not known until
        # the geometry is complete. It is left empty here and filled in build(), exactly as
        # prj.water.table() is; writing a one-point wy with an empty wx produced a project the
        # validator refused ("0 x-value(s) but 1 y-value(s)"), which made the documented
        # scalar form unusable.
        if water is not None:
            ph.water_override = True
            if isinstance(water, (int, float)):
                ph.wx = []                                    # spans, filled at build
                ph.wy = [float(water), float(water)]
            else:
                pts = list(water)
                ph.wx = [float(p[0]) for p in pts]
                ph.wy = [float(p[1]) for p in pts]
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

        ``apply_fraction=0.5`` applies half of the stage (Sum-Mstage) and
        leaves the rest; ``ignore_undrained=True`` solves undrained materials as
        drained for this phase; ``reset_small_strain=True`` clears the HS-small
        strain history first, so the soil meets this phase at G0 instead of the
        stiffness the earlier phases degraded it to."""
        return self._add(name, _core.PhaseType.Plastic, **kw)

    def consolidation(self, name, *, duration, steps, **kw):
        """Time-dependent Biot consolidation. duration [day]."""
        return self._add(name, _core.PhaseType.Consolidation,
                         duration=duration, steps=steps, **kw)

    def safety(self, name="Safety", **kw):
        """phi-c reduction of the current state -> factor of safety."""
        return self._add(name, _core.PhaseType.Safety, **kw)

    def transient_flow(self, name, *, duration, steps, **kw):
        """Time-dependent GROUNDWATER FLOW alone: heads and the free surface move,
        the skeleton does not. duration [day], ``steps`` time steps over it. Use it
        for a drawdown or a filling reservoir when the deformation is not the
        question; for both at once use :meth:`fully_coupled`."""
        return self._add(name, _core.PhaseType.TransientFlow,
                         duration=duration, steps=steps, **kw)

    def fully_coupled(self, name, *, duration, steps, **kw):
        """Fully coupled flow-deformation: pore pressure and displacement solved
        TOGETHER over time. Differs from :meth:`consolidation` in that the
        hydraulic boundary conditions may change during the phase -- the coupling
        runs both ways rather than dissipating a pressure that is already there."""
        return self._add(name, _core.PhaseType.FullyCoupled,
                         duration=duration, steps=steps, **kw)

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
        self._structs = []
        self._plates = []
        self._anchors = []
        self._geogrids = []
        self._embedded = []
        self._disps = []
        self._hydros = []
        self._water = None
        self._procedure = "k0"
        self._initial_exclude = []
        self._initial_numerics = {}
        self._phase_toggles = []
        self.materials = _Materials(self)
        self.geometry = _Geometry(self)
        self.water = _Water(self)
        self.loads = _Loads(self)
        self.structures = _Structures(self)
        self.displacements = _Displacements(self)
        self.dewatering = _Dewatering(self)
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
        # A wall toed below the soil block, or an anchor reaching outside it, is still part of
        # the model: the extent has to contain it or the mesher never sees the element.
        xs += [v for S in self._structs for v in (S.x1, S.x2)]
        ys += [v for S in self._structs for v in (S.y1, S.y2)]
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
        pr.structs = self._structs
        pr.plates = self._plates
        pr.anchors = self._anchors
        pr.geogrids = self._geogrids
        pr.embedded = self._embedded
        pr.disps = self._disps
        pr.hydros = self._hydros

        # The initial phase: everything active except the exclusions. Vectors are
        # materialized ONLY for a touched class, exactly like the GUI's files.
        initial = _core.Phase()
        excl_loads = [h.index for h in self._initial_exclude
                      if isinstance(h, LoadHandle)]
        excl_regions = [h.index for h in self._initial_exclude
                        if isinstance(h, RegionHandle)]
        excl_disps = [h.index for h in self._initial_exclude
                      if isinstance(h, DispHandle)]
        excl_hydros = [h.index for h in self._initial_exclude
                       if isinstance(h, HydroHandle)]
        excl_structs = [h.index for h in self._initial_exclude
                        if isinstance(h, StructHandle)]
        load_state = [1] * len(self._loads)
        poly_state = [1] * len(self._polygons)
        disp_state = [1] * len(self._disps)
        hydro_state = [1] * len(self._hydros)
        struct_state = [1] * len(self._structs)
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
        if excl_hydros:
            for i in excl_hydros:
                hydro_state[i] = 0
            initial.hydro_active = hydro_state
        if excl_structs:
            for i in excl_structs:
                struct_state[i] = 0
            initial.struct_active = struct_state
        for field, value in self._initial_numerics.items():
            if value is not None:
                setattr(initial, field, value)
        pr.initial = initial

        # Staged phases: inherit the previous phase's effective state; write only
        # the vectors a toggle touches.
        phases = []
        for ph, activate, deactivate in self._phase_toggles:
            touched_loads = touched_polys = touched_disps = False
            touched_hydros = touched_structs = False
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
                elif isinstance(h, HydroHandle):
                    hydro_state[h.index] = on
                    touched_hydros = True
                elif isinstance(h, StructHandle):
                    struct_state[h.index] = on
                    touched_structs = True
                else:
                    raise ValueError(f"cannot toggle {h!r}")
            if touched_loads:
                ph.load_active = list(load_state)
            if touched_polys:
                ph.poly_active = list(poly_state)
            if touched_disps:
                ph.disp_active = list(disp_state)
            if touched_hydros:
                ph.hydro_active = list(hydro_state)
            if touched_structs:
                ph.struct_active = list(struct_state)
            # A phase's horizontal water table spans the model, like the project's own.
            if ph.water_override and not ph.wx:
                ph.wx = [pr.x_min, pr.x_max]
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
