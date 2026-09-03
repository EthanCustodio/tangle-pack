import logging

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)-8s %(name)s: %(message)s",
    datefmt="%H:%M:%S",
)

import tanglepack, numpy as np
from tanglepack.examples import (
    HENON_P3,
    henon_map as _henon_map_factory,
    henon_map_inverse as _henon_map_inverse_factory,
    henon_jacobian as _henon_jacobian_factory,
)
import matplotlib.pyplot as plt


henon_map = _henon_map_factory(*HENON_P3)
henon_map_inverse = _henon_map_inverse_factory(*HENON_P3)
henon_jacobian = _henon_jacobian_factory(*HENON_P3)


# ── Numeric phase ──────────────────────────────────────────────────────────
wb = tanglepack.TangleWorkbench(henon_map, henon_map_inverse, henon_jacobian)

wb._man_machine.area_cutoff = 1e-8

# Period-3 inner fixed point
fp3 = wb.construct_fixed_point([[0, 1], [-1, 0], [-1, 1]])
wb.orient_eigenvectors(
    fp3, {"unstable": np.array([0, -1]), "stable": np.array([-1, -1])}
)
wb.initialize_both_manifolds(fp3)
wb.grow_n_times(fp3, "unstable", num_iterations=14)
wb.grow_n_times(fp3, "stable", num_iterations=14)


# Outer period-1 fixed point
wb._man_machine.area_cutoff = 1e-4
fp1 = wb.construct_fixed_point([4, -4])
wb.orient_eigenvectors(fp1, {"unstable": np.array([-1, 0]), "stable": np.array([0, 1])})
wb.initialize_both_manifolds(fp1)
wb.grow_n_times(fp1, "unstable", num_iterations=7)
wb.grow_n_times(fp1, "stable", num_iterations=7)

# ── Intersection phase (period-3 homoclinic intersections) ─────────────────
wb.compute_intersections([fp3, fp1])
wb.trim_stable_manifolds(fp3)
wb.trim_stable_manifolds(fp1)

bridges = wb.create_bridges(fp3)
bridges1 = wb.create_bridges(fp1)

new_links = wb.infer_iterate_table()
print(f"Recorded {new_links} iterate relationships")

# ── Topological phase ──────────────────────────────────────────────────────
registry = wb.intersection_registry
print(f"Total intersections: {len(registry)}")

first_id = registry.by_unstable_cdist[0]
p = registry[first_id]
print(
    f"Intersection {first_id}: coords={p.coords}, u_cdist={p.unstable_cdist:.4f}, s_cdist={p.stable_cdist:.4f}"
)
print(f"  Fixed point: {p.manifold_a_key[0].coordinates[0].ravel()}")
print(f"  Branch: {p.manifold_a_key[3]}")

# Iterate chains
fwd = registry.iterate_table.forward_chain(first_id)
print(f"Forward chain from {first_id}: {fwd}")

# Orderings
u_order = registry.by_unstable_cdist
s_order = registry.by_stable_cdist

# --- Query interface ---

# Which intersections will map into the cdist range [5, 10] on the stable manifold?
sources = registry.on_interval(5.0, 10.0, stability="stable")
print(f"{len(sources)} intersections map into s-cdist [5, 10]")

# All intersections currently sitting in a cdist range (no iteration)
in_range = registry.on_cdist_range(0.0, 1.0, stability="stable")
print(f"{len(in_range)} intersections have s-cdist in [0, 1]")

# All intersections involving the period-3 fixed point
from_fp3 = registry.from_fixed_point(fp3)

# Only branch-0 intersections on the stable side
branch0_stable = registry.from_branch(0, stability="stable")

# Custom filter: only intersections whose iterate is also registered
has_iterate = registry.filter(lambda ix: (ix.id, 1) in registry.iterate_table)

# --- Live graph ---
G = wb.build_intersection_graph()
print(f"Graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

# Adjacency edges along W^u
adj_u = [(u, v) for u, v, d in G.edges(data=True) if d.get("stability") == "unstable"]
# Iterate edges
iter_edges = [
    (u, v, d["n"]) for u, v, d in G.edges(data=True) if d.get("type") == "iterate"
]

# ── Tangle plot ────────────────────────────────────────────────────────────
plt.figure()
wb.plot_tangle(fp3, "unstable", color="b")
wb.plot_tangle(fp3, "stable", color="r")
wb.plot_tangle(fp1, "unstable", color="b")
wb.plot_tangle(fp1, "stable", color="r")
wb.plot_intersections(fp3)
wb.plot_intersections(fp1)
wb.plot_all_bridges()

plt.xlim([-1, 1])
plt.ylim([0.28, 3.1])
plt.title("k=2, b=1 Period 3 Nested Tangle")
plt.tight_layout()

# --- Dense array exports ---
F = registry.as_forward_array(max_depth=5)  # shape (N, 5), -1 = unknown
B = registry.as_backward_array(max_depth=5)

plt.show()
