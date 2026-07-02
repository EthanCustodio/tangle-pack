import tanglepack, numpy as np
import matplotlib.pyplot as plt


def henon_map(point):
    k, b = 10, 1
    x, y = point
    return np.array([y - k + x**2, -b * x])


def henon_map_inverse(point):
    k, b = 10, 1
    x, y = point
    return np.array([-y / b, x + k - (y**2) / (b**2)])


# ── Numeric phase ──────────────────────────────────────────────────────────
wb = tanglepack.TangleWorkbench(henon_map, henon_map_inverse)
fp = wb.construct_fixed_point([4, -4])
wb.orient_eigenvectors(fp, {"unstable": np.array([-1, 0]), "stable": np.array([0, 1])})
wb.initialize_both_manifolds(fp)
wb.grow_n_times(fp, "unstable", num_iterations=8)
wb.grow_until_turnaround(fp, "stable")

wb.compute_intersections(fp)  # also populates registry
wb.trim_stable_manifolds(fp)
bridges = wb.create_bridges(fp)

# for _ in range(3):
#     wb.iterate_all_bridges()
new_bridges = wb.iterate_bridge(bridges[2])
# wb.iterate_bridge(new_bridges[0])


new_links = wb.infer_iterate_table()
print(f"Recorded {new_links} iterate relationships")

# ── Topological phase ──────────────────────────────────────────────────────
registry = wb.intersection_registry
print(f"Total intersections: {len(registry)}")

# Lookup by ID
p = registry[3]
print(
    f"Intersection 3: coords={p.coords}, u_cdist={p.unstable_cdist:.4f}, s_cdist={p.stable_cdist:.4f}"
)
print(f"  Fixed point: {p.manifold_a_key[0].coordinates[0].ravel()}")
print(f"  Branch: {p.manifold_a_key[3]}")

# Iterate chains
fwd = registry.iterate_table.forward_chain(3)
print(f"Forward chain from 3: {fwd}")

# Orderings
u_order = registry.by_unstable_cdist
s_order = registry.by_stable_cdist

# --- Query interface ---

# Which intersections will map into the cdist range [5, 10] on the unstable manifold?
sources = registry.on_interval(5.0, 10.0, stability="stable")
print(f"{len(sources)} intersections map into s-cdist [5, 10]")

# All intersections currently sitting in a cdist range (no iteration)
in_range = registry.on_cdist_range(2.0, 6.0, stability="stable")
print(f"{len(in_range)} intersections have s-cdist in [2, 8]")

# All intersections involving the specific fixed point
from_fp = registry.from_fixed_point(fp)

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
# plt.figure(figsize=(8, 8))
# wb.plot_tangle(fp, "unstable", color="b")
# wb.plot_tangle(fp, "stable", color="r")
# wb.plot_intersections(fp)
# plt.title("Henon Tangle")
# plt.axis("equal")
# plt.tight_layout()
# plt.show()

plt.figure()
wb.plot_tangle(fp, "stable", color="r")
wb.plot_intersections(fp)
wb.plot_all_bridges()

plt.xlim([-15, 15])
plt.ylim([-15, 15])

# ── Intersection graph ─────────────────────────────────────────────────────
wb.visualize_intersection_graph(G, label_mode="all")

# --- Dense array exports ---
F = registry.as_forward_array(max_depth=5)  # shape (N, 5), -1 = unknown
B = registry.as_backward_array(max_depth=5)
