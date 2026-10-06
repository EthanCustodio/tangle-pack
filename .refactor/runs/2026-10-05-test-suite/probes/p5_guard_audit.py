"""P5: regression-guard audit of the deletion candidates (plan Phases 2 and 7).

For each candidate test, greps the nine regression-guard memory notes for the
test name and for subject keywords, and prints every hit (file:line, trimmed)
for a human verdict. The verdicts are recorded in docs/test_suite_refactor_ledger.md §0.
"""

from __future__ import annotations

import re
from pathlib import Path

MEMORY = Path.home() / ".claude/projects/-home-dezhu-Development-TanglePack-tangle-pack/memory"
NOTES = [
    "codebase-audit-2026-09.md", "pseudoneighbor-collision-fragility.md", "pseudoneighbor-partition.md",
    "numerics-test-suite.md", "regions-refactor-2026-09.md", "topology-layer.md",
    "cdist-strict-monotonicity-fix.md", "hole-side-own-blast-2026-10-02.md", "audit-polish-2026-07.md",
]

# candidate -> subject keywords (regex, case-insensitive)
CANDIDATES: dict[str, list[str]] = {
    "visualizations/": [r"viz_", r"visualizations"],
    "test_workbench_split.py": [r"workbench.split", r"patchab"],
    "test_manifold_machine.py": [r"pre-1\.8", r"grow loop"],
    "test_deleted_genealogy_attributes_are_gone": [r"genealogy"],
    "test_grown_until_intersection_is_gone": [r"grown_until", r"grow_until_intersection"],
    "test_advance_key_returns_the_next_orbit_point": [r"advance_key"],
    "test_per_step_beta_matches_the_expression_it_replaces": [r"per_step_beta", r"k.th root", r"k_value.root"],
    "test_k10_intersection_graph_is_unchanged": [r"8/18", r"graph count"],
    "test_branch_cycle_reproduces_the_topology_branch_orders": [r"branch_cycle"],
    "test_bridge_cutting_pin_k10/_period_3": [r"bridge.cutting", r"cut at the exact"],
    "test_tangle_keeps_only_index_state": [r"index state"],
    "test_string_dispatched_iter_method_is_gone": [r"string.dispatch"],
    "test_stability_alias_is_defined_once_in_numerics": [r"Stability alias", r"alias"],
    "test_no_test_or_script_defines_its_own_henon": [r"own henon", r"henon_map\("],
    "test_small_tangle_crossing_count_is_unchanged": [r"crossing count"],
    "test_fixed_point_takes_no_branch_count": [r"num_branches", r"branch count"],
    "test_get_lambda_u_* / test_on_interval_*": [r"_get_lambda_u", r"on_interval", r"lambda_u"],
    "test_curvature_area_batch_matches_scalar": [r"_curvature_area"],
    "test_invalidate_trellises_is_a_deprecated_no_op": [r"invalidate_trellises"],
    "test_boundary_arc_dataclass_is_gone / test_shapely_is_not_a_dependency": [r"BoundaryArc", r"shapely"],
    "test_k10/p3_zone_area(s)_and_containment_are_unchanged": [r"zone area", r"_segments_of", r"tail fix"],
    "test_initialization_unstable/stable": [r"3 init", r"three init", r"init points"],
    "test_cdists_are_positive_and_increasing": [r"positive and increasing"],
    "test_image_chain_is_the_bridge_class_function": [r"_image_chain", r"image_chain"],
    "test_kinds_are_distinct": [r"kind"],
    "test_names_are_deterministic_across_builds": [r"deterministic", r"not deterministic"],
    "plotting layout tests (4)": [r"fanout helper", r"hole style", r"drawing code"],
    "test_henon_partition_still_built_after_wiring": [r"wiring"],
    "test_image_cdist (3: additive / scaled-element / every-bounded)": [r"image_cdist", r"additive"],
    "test_image_of_element_lands_on_the_advanced_branch": [r"image_of_element"],
    "test_session_result_matches_gathered_partitions_helper / test_p3_no_class_mixes_fixed_points": [r"_gathered_partitions", r"mixes"],
    "test_describe_reports_the_image_evidence": [r"describe"],
    "test_henon_reference_holes_hug_the_stable_manifold / _stable_arc_midpoint x2": [r"_stable_arc_midpoint", r"hug", r"chord midpoint", r"chord"],
    "test_forward_pairs_beyond_the_fundamental_segment_get_no_hole[1..k-1]": [r"forward.hole", r"\+1\.\.\+\(k", r"k_value"],
    "summary/repr/describe string pins (3)": [r"summary\(\)"],
    "P7: test_trim_stable_manifolds_cuts_just_past_the_outermost_crossing": [r"trim_stable", r"outermost crossing"],
    "P7: test_forward_unstable_branch_cycle_matches_orbit_order": [r"forward_unstable_branch_cycle"],
    "P7: mixed-sign rejection x2": [r"mixed", r"set_k_value", r"disagree"],
    "P7: test_faces_closed_raises_at_the_cap_for_the_anchor": [r"faces_closed", r"anchor.*cap"],
    "P7: test_refined_cdist_is_mean_of_neighbours": [r"mean of", r"refined cdist", r"midpoint cdist"],
    "P7: _region_key test": [r"_region_key", r"singleton bridge"],
    "P7: test_sparse_kept_endpoints_are_degree_three_without_unstable_stubs": [r"degree", r"half.edge"],
    "P7: test_the_snap_does_not_fire_on_the_default_path": [r"snap"],
}


def main() -> None:
    """Print the hits per candidate."""
    texts = {n: (MEMORY / n).read_text().splitlines() for n in NOTES if (MEMORY / n).exists()}
    missing = [n for n in NOTES if n not in texts]
    if missing:
        print("MISSING NOTES:", missing)
    for cand, kws in CANDIDATES.items():
        names = re.findall(r"test_[A-Za-z0-9_]+", cand)
        pats = [re.compile(re.escape(n)) for n in names] + [re.compile(k, re.I) for k in kws]
        hits = []
        for note, lines in texts.items():
            for i, line in enumerate(lines, 1):
                if any(p.search(line) for p in pats):
                    hits.append(f"{note}:{i}: {line.strip()[:170]}")
        print(f"\n## {cand}: {len(hits)} hit(s)")
        for h in hits[:8]:
            print("   ", h)


if __name__ == "__main__":
    main()
