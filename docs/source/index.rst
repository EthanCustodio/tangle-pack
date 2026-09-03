.. tanglepack documentation master file, created by
   sphinx-quickstart on Wed Jul 30 12:04:06 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

tanglepack documentation
========================

``tanglepack`` computes and visualises heteroclinic/homoclinic tangles — the stable
and unstable manifolds of saddle fixed points in 2D area-preserving maps.

The library is organised into three subpackages:

* :mod:`tanglepack.numerics` — the numerical engine (dynamical system, fixed points,
  manifold growth, intersection detection, bridges). Entry point:
  :class:`~tanglepack.numerics.TangleWorkbench.TangleWorkbench`.
* :mod:`tanglepack.topology` — the topological view of a computed tangle
  (:class:`~tanglepack.topology.Trellis.Trellis`, strong pips, pseudoneighbors,
  the planar :class:`~tanglepack.topology.Arrangement.Arrangement` and its regions).
* :mod:`tanglepack.loom` — cross-layer algorithms that read topology results and act on
  the numerical layer (:class:`~tanglepack.loom.TangleSession.TangleSession`, resonance
  zones, blasting).

:mod:`tanglepack.examples` holds the Hénon fixtures used by the scripts, notebooks and
tests. The API reference below is generated with ``sphinx-apidoc``; regenerate it with
``sphinx-apidoc -o docs/source src/tanglepack -f --no-toc`` after adding a module.


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   tanglepack
