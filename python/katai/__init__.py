"""KATAI 2D -- scripting over the published facade.

This package is the engineer-facing layer over ``katai._core`` (the 1:1 nanobind
binding of the C++ facade). Units are the schema's own and FIXED: kN, m, day --
values pass through verbatim, and every parameter's docstring states its unit.

First increment (contract-first): the raw ``_core`` surface re-exported, plus the
``Refusal`` exception and the ``run()`` convenience. The engineer-facing
construction DSL (``Project.geometry.polygon(...)``, ``materials.mohr_coulomb(...)``,
phase chaining) lands on top of this in the next increment; the checked-in input
corpus is its acceptance suite.

The scripting boundary (roadmap section 7.3) is inherited from the facade: the
input contract is fully scriptable, the implementation is closed -- no stiffness
matrices, no element state, no solver iterations.
"""

import os as _os

# Windows: Python 3.8+ ignores PATH when resolving extension-module DLL
# dependencies (the dynamic MKL runtime here). KATAI_DLL_PATH lists the
# directories to trust; a future wheel bundles its runtime instead and this
# hook stays a harmless no-op.
if hasattr(_os, "add_dll_directory"):
    for _d in _os.environ.get("KATAI_DLL_PATH", "").split(_os.pathsep):
        if _d and _os.path.isdir(_d):
            _os.add_dll_directory(_d)

from . import _core
from ._core import (  # noqa: F401 -- the sanctioned vocabulary, re-exported by name
    AnchorMaterial, BCType, DesignApproach, Drainage, EmbeddedBeamMaterial,
    FlowBCType, GeogridMaterial, InitialProcedure, Issue, Job, JobState, Load,
    LoadKind, Material, MeshSettings, Phase, PhaseType, PlateMaterial, PrescribedDisp,
    SeismicWave, Severity, SoilModel, SoilPolygon, StructElement, StructKind,
    ValidationReport, backend_name, fnv1a64, load_project, project_from_json,
    project_to_json, save_project, validate_project,
    PROJECT_FILE_VERSION, RESULTS_FILE_VERSION,
)

__version__ = _core.__version__   # the ONE identity (<katai/api/version.hpp>)

# `katai.Project` is the ENGINEER-FACING builder (easy to learn, easy to apply --
# the maintainer's rule). The raw schema struct stays available as
# `katai._core.Project`; `load_project` returns that raw form.
from .project import Project  # noqa: F401,E402


class Refusal(RuntimeError):
    """The engine or the input contract refused the run, honestly.

    The message carries the refusal text byte-for-byte as the engine states it;
    ``report`` (when the refusal came from validation) holds the field-path issues.
    """

    def __init__(self, message, report=None):
        super().__init__(message)
        self.report = report


def run(project, on_phase=None):
    """Run every phase of ``project`` headlessly; return the finished Job.

    The same layer-2 work the GUI's solve button and ``katai solve`` submit: the
    input contract executes first (an invalid project raises :class:`Refusal` with
    the field-path report), the mesh comes from the project's own settings, and
    every phase goes through the engine's phase strategies.

    ``on_phase(index, total, name)`` fires before each phase solve.
    """
    job = Job(project)
    if on_phase is not None:
        job.set_on_phase(on_phase)
    if not job.run():
        rep = job.report()
        raise Refusal(job.message(), report=None if rep.ok() else rep)
    return job
