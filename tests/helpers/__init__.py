"""Shared, non-collected helpers of the tanglepack test suite.

* :mod:`helpers.invariants` -- manifold-level checks (cdist order, spikes,
  iterate law, one-to-one, area along a chain).
* :mod:`helpers.laws` -- the physical-law checks over a built
  :class:`cases.Case`, each returning how many items it checked.
* :mod:`helpers.law_tier` -- the law tier's plumbing: shared-build products,
  the fingerprint guard's snapshot, the ``law_test`` factory.
* :mod:`helpers.fakes` -- synthetic fixed points, partitions, trellises and
  dual graphs for the walk and naming kernels.
* :mod:`helpers.names` -- letter-free spellers of classes, words, brackets and
  transition matrices.
* :mod:`helpers.logs` -- ``assert_logged`` (level and logger only).
"""
