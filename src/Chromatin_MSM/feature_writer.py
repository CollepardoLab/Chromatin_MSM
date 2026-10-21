#!/usr/bin/env python3
"""
multi_feature_writer.py — PLUMED collective variable feature library

Feature types (--features):
  distance   Pairwise distances between the 4 nucleosome core centers (6 pairs)
  angle      Core triple-angles (1-2-3, 2-3-4) with sin/cos transformations
  torsion    Torsion angle across 4 cores with sin/cos transformations
  bend       DNA linker bending angles (3 linkers) with sin/cos transformations
  contacts   Histone-tail / DNA-phosphate / core coordination numbers
             plus per-tail and per-nucleosome phosphate COMs
  rg         Radius of gyration (non-histone-tail atoms)
  rmsd       RMSD to reference structures (DNA and core_all, from --ref-dir)
  seg_dist   Pairwise core-segment distances with exp(-x) transformations

Usage examples
--------------
  # All features, linker length 15
  python multi_feature_writer.py --pdb ref/15/15rep.pdb --linker-length 15

  # Only distances and angles
  python multi_feature_writer.py --pdb ref/15/15rep.pdb \\
      --features distance angle torsion

  # Contacts only, with custom cutoffs
  python multi_feature_writer.py --pdb ref/15/15rep.pdb \\
      --linker-length 15 --features contacts rg \\
      --htail-phos-cf 3.8 --htail-nucle-cf 3.6

  # Custom tail ranges from JSON file
  python multi_feature_writer.py --pdb ref/15/15rep.pdb \\
      --tail-ranges tail_ranges.json --output multi_feature.dat
"""

import argparse
import json
import os
import sys
import mdtraj as md
from itertools import combinations


ALL_FEATURES = ["distance", "angle", "torsion", "bend", "contacts", "rg", "rmsd", "seg_dist"]

# Default histone tail residue index ranges for each chain (A-D), based on the reference structure. These can be overridden by a user-provided JSON file.
DEFAULT_TAIL_RANGES = {
    'A': [(1,40),(136,160),(238,257),(351,364),(366,390),
          (488,527),(623,647),(725,744),(838,851),(853,877)],
    'B': [(975,1013),(1110,1134),(1212,1231),(1324,1337),(1339,1362),
          (1462,1507),(1596,1619),(1698,1716),(1811,1824),(1827,1850)],
    'C': [(1949,1987),(2084,2108),(2186,2205),(2299,2312),(2314,2337),
          (2435,2473),(2570,2593),(2672,2690),(2786,2799),(2801,2824)],
    'D': [(2923,2962),(3058,3082),(3160,3179),(3273,3286),(3288,3312),
          (3409,3447),(3544,3567),(3646,3664),(3760,3773),(3775,3799)],
}

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--pdb", required=True,
        help="Reference PDB file path",
    )
    p.add_argument(
        "--output", default="mdp_features.dat",
        help="Output PLUMED .dat file (default: mdp_features.dat, matching what "
             "adaptive.py's --feature-file expects)",
    )
    p.add_argument(
        "--features", nargs="+", default=["all"],
        choices=ALL_FEATURES + ["all"],
        metavar="FEATURE",
        help=(
            "Feature types to include (space-separated). "
            f"Choices: all {ALL_FEATURES}. Default: all."
        ),
    )
    p.add_argument(
        "--ref-dir", default=None,
        help="Directory holding the RMSD references DNA.pdb and core_all.pdb, "
             "as PLUMED will see it from its working directory "
             "(default: the directory containing --pdb, e.g. ref/25_A)",
    )
    p.add_argument(
        "--linker-length", type=int,
        help="Linker DNA length",
    )
    p.add_argument(
        "--n-segments-core", type=int, default=15,
        help="Number of segments per nucleosome core (default: 15)",
    )
    p.add_argument(
        "--segment-stride", type=int, default=5,
        help="Atom stride within each core-segment CENTER (default: 5)",
    )
    p.add_argument(
        "--htail-phos-cf", type=float, default=3.6,
        help="Coordination cutoff: histone tail vs phosphate (default: 3.6)",
    )
    p.add_argument(
        "--htail-nucle-cf", type=float, default=3.6,
        help="Coordination cutoff: histone tail vs core segment (default: 3.6)",
    )
    p.add_argument(
        "--tail-ranges", default=None,
        help=(
            "JSON file defining tail residue ranges: "
            '{"A": [[start,end], ...], "B": [...], ...}. '
            "Uses built-in defaults when omitted."
        ),
    )
    return p.parse_args()


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def resolve_features(feature_list):
    """Expand 'all' shorthand and return a set of feature names."""
    if "all" in feature_list:
        return set(ALL_FEATURES)
    return set(feature_list)


def load_tail_ranges(path):
    with open(path) as fh:
        raw = json.load(fh)
    return {k: [tuple(pair) for pair in v] for k, v in raw.items()}


# ---------------------------------------------------------------------------
# Main writer
# ---------------------------------------------------------------------------

def write_features(
    feature_file,
    pdb_file,
    features,
    tail_ranges,
    linker_length,
    n_segments_core,
    segment_stride,
    htail_phos_cf,
    htail_nucle_cf,
    ref_dir=None,
):
    # RMSD references live alongside the reference topology (ref/<L>[_A|_AT]/).
    if ref_dir is None:
        ref_dir = os.path.dirname(pdb_file) or "."
    traj = md.load(pdb_file)
    top = traj.topology
    chain_labels = ['A', 'B', 'C', 'D']

    # ---- Atom index collections ----------------------------------------
    core_indices = [
        top.select(f'chainid == {i} and element Pt').tolist()
        for i in range(top.n_chains - 1)
    ]
    if len(core_indices) < 4:
        raise ValueError("Topology must contain at least 4 nucleosome cores (Pt chains).")
    core_indices = core_indices[:4]

    DNA_indices = top.select('element Au').tolist()
    n_dna = len(DNA_indices)

    # Histone tail: Cu atoms in range, filtered to pR/pK for contacts
    filtered_tail_selections = {}
    for chain, ranges in tail_ranges.items():
        for i, (start, end) in enumerate(ranges):
            label = f"tail{i+1}_{chain}"
            inds = top.select(f"element Cu and resid {start} to {end}").tolist()
            filtered_tail_selections[label] = [
                idx for idx in inds
                if top.atom(idx).residue.name in ('pR', 'pK')
            ]

    # Dyad positions and per-nucleosome DNA windows
    first_dyad = linker_length // 2 + 73
    dyad_positions = [first_dyad + i * (147 + linker_length) for i in range(4)]

    dna_groups = {
        chain: DNA_indices[dyad_positions[idx] - 63 : dyad_positions[idx] + 64]
        for idx, chain in enumerate(chain_labels)
    }
    # Au-indexed phosphate beads (offset by n_dna in the CG topology)
    phos_groups = {
        chain: [idx + n_dna for idx in inds]
        for chain, inds in dna_groups.items()
    }
    # Actual phosphorus atoms (element P) used for coordination contacts
    phos_P_indices = top.select('element P').tolist()

    # ---- Decide which prerequisite blocks are required ------------------
    need_core_centers  = bool(features & {"distance", "angle", "torsion"})
    need_core_segments = bool(features & {"contacts", "seg_dist"})
    need_core_all      = "contacts" in features
    need_tail_ct       = "contacts" in features
    need_phos_groups   = "contacts" in features
    need_bend_groups   = "bend" in features

    # ---- Validate index bounds for bend feature -------------------------
    if need_bend_groups:
        max_bend_idx = max(
            dyad_positions[idx] + 63 + linker_length + 20
            for idx in range(len(chain_labels) - 1)
        )
        if max_bend_idx >= n_dna:
            print(
                f"WARNING: bend endpoint index {max_bend_idx} exceeds "
                f"DNA_indices length {n_dna}. Check --linker-length.",
                file=sys.stderr,
            )

    with open(feature_file, 'w') as f:

        # ==================================================================
        # File header
        # ==================================================================
        f.write(f"# PLUMED input: tetranucleosomes CG model\n")
        f.write(f"# Features   : {', '.join(sorted(features))}\n")
        f.write(f"# Linker     : {linker_length} base pairs\n")
        f.write(f"# Core segs  : {n_segments_core}  stride={segment_stride}\n\n")
        f.write(f"MOLINFO STRUCTURE={pdb_file}\n\n")

        # ==================================================================
        # BLOCK 1 — Shared GROUP / CENTER prerequisites
        # (defined once; referenced by multiple feature blocks below)
        # ==================================================================

        # --- Core center virtual atoms (distance / angle / torsion) ------
        if need_core_centers:
            f.write("# Core center virtual atoms\n")
            for i, inds in enumerate(core_indices, start=1):
                f.write(f"core{i}_center: CENTER ATOMS={','.join(map(str, inds))}\n")
            f.write("\n")

        # --- Core segment virtual atoms (contacts / seg_dist) ------------
        core_segment_names = {}
        if need_core_segments:
            f.write("# Core segment virtual atoms\n")
            for core_idx, core in enumerate(core_indices, start=1):
                core_name = f"core{core_idx}"
                core_segment_names[core_name] = []
                seg_length = max(1, len(core) // n_segments_core)
                for seg in range(n_segments_core):
                    seg_start = core[seg * seg_length]
                    seg_end = (
                        core[-1] if seg == n_segments_core - 1
                        else core[(seg + 1) * seg_length - 1]
                    )
                    seg_name = f"{core_name}_seg{seg}"
                    core_segment_names[core_name].append(seg_name)
                    f.write(f"{seg_name}: CENTER ATOMS={seg_start}-{seg_end}:{segment_stride}\n")
            f.write("\n")

        # --- Core overall centers (contacts only) ------------------------
        if need_core_all:
            f.write("# Core overall centers (mean of segment virtual atoms)\n")
            for core_name, segs in core_segment_names.items():
                f.write(f"{core_name}_all: CENTER ATOMS={','.join(segs)}\n")
            f.write("\n")

        # --- Histone tail contact groups (contacts only) -----------------
        if need_tail_ct:
            f.write("# Histone tail contact groups (pR/pK residues)\n")
            for label, inds in filtered_tail_selections.items():
                if inds:
                    f.write(f"{label}_ct: GROUP ATOMS={','.join(map(str, inds))}\n")
            f.write("\n")

        # --- DNA phosphate bead groups and actual P atoms (contacts) -----
        if need_phos_groups:
            f.write("# Per-nucleosome DNA phosphate bead groups (Au-offset in CG)\n")
            for chain, inds in phos_groups.items():
                if inds:
                    f.write(f"phos{chain}: GROUP ATOMS={','.join(map(str, inds))}\n")
            f.write("\n")
            if phos_P_indices:
                f.write("# All phosphorus (element P) atoms for coordination\n")
                f.write(f"phos_all: GROUP ATOMS={','.join(map(str, phos_P_indices))}\n")
                f.write("\n")

        # --- DNA linker bending endpoint groups (bend only) --------------
        if need_bend_groups:
            f.write("# DNA linker bending: endpoint atom groups\n")
            for idx, chain in enumerate(chain_labels[:-1]):
                p1 = DNA_indices[dyad_positions[idx] + 63]
                p2 = DNA_indices[dyad_positions[idx] + 63 + linker_length + 20]
                f.write(f"dna{chain}_start: GROUP ATOMS={p1}\n")
                f.write(f"dna{chain}_end:   GROUP ATOMS={p2}\n")
            f.write("\n")
            f.write("# DNA linker midpoint virtual atoms\n")
            for chain in chain_labels[:-1]:
                f.write(f"dna{chain}_center: CENTER ATOMS=dna{chain}_start,dna{chain}_end\n")
            f.write("\n")

        # ==================================================================
        # BLOCK 2 — Feature CVs
        # ==================================================================

        # ------------------------------------------------------------------
        # distance: pairwise distances between core centers
        # ------------------------------------------------------------------
        if "distance" in features:
            f.write("# ---- Core-center pairwise distances ----\n")
            for i, j in combinations(range(1, 5), 2):
                f.write(f"dist_core{i}_{j}: DISTANCE ATOMS=core{i}_center,core{j}_center\n")
            f.write("\n")

        # ------------------------------------------------------------------
        # angle: triple angles (1-2-3, 2-3-4) with sin/cos
        # ------------------------------------------------------------------
        if "angle" in features:
            f.write("# ---- Core triple angles ----\n")
            f.write("angle_123: ANGLE ATOMS=core1_center,core2_center,core3_center\n")
            f.write("angle_234: ANGLE ATOMS=core2_center,core3_center,core4_center\n")
            f.write("sin_angle_123: CUSTOM ARG=angle_123 FUNC=sin(x) PERIODIC=NO\n")
            f.write("cos_angle_123: CUSTOM ARG=angle_123 FUNC=cos(x) PERIODIC=NO\n")
            f.write("sin_angle_234: CUSTOM ARG=angle_234 FUNC=sin(x) PERIODIC=NO\n")
            f.write("cos_angle_234: CUSTOM ARG=angle_234 FUNC=cos(x) PERIODIC=NO\n")
            f.write("\n")

        # ------------------------------------------------------------------
        # torsion: torsion angle across 4 cores with sin/cos
        # ------------------------------------------------------------------
        if "torsion" in features:
            f.write("# ---- Core torsion angle ----\n")
            f.write("psi_nucl: TORSION ATOMS=core1_center,core2_center,core3_center,core4_center\n")
            f.write("sin_psi_nucl: CUSTOM ARG=psi_nucl FUNC=sin(x) PERIODIC=NO\n")
            f.write("cos_psi_nucl: CUSTOM ARG=psi_nucl FUNC=cos(x) PERIODIC=NO\n")
            f.write("\n")

        # ------------------------------------------------------------------
        # bend: DNA linker bending angles with sin/cos
        # ------------------------------------------------------------------
        if "bend" in features:
            f.write("# ---- DNA linker bending angles ----\n")
            for chain in chain_labels[:-1]:
                f.write(f"bend_{chain}: ANGLE ATOMS=dna{chain}_start,dna{chain}_center,dna{chain}_end\n")
                f.write(f"sin_bend_{chain}: CUSTOM ARG=bend_{chain} FUNC=sin(x) PERIODIC=NO\n")
                f.write(f"cos_bend_{chain}: CUSTOM ARG=bend_{chain} FUNC=cos(x) PERIODIC=NO\n")
            f.write("\n")

        # ------------------------------------------------------------------
        # contacts: histone-tail / DNA / core coordination + COMs
        # ------------------------------------------------------------------
        if "contacts" in features:
            f.write("# ---- Histone-tail COMs ----\n")
            for label, inds in filtered_tail_selections.items():
                if inds:
                    f.write(f"{label}_ct_com: CENTER ATOMS={label}_ct\n")
            f.write("\n")

            f.write("# ---- Per-nucleosome DNA phosphate COMs ----\n")
            for chain in phos_groups:
                f.write(f"phos{chain}_com: CENTER ATOMS=phos{chain}\n")
            f.write("\n")

            f.write("# ---- Histone-tail coordination numbers ----\n")
            for label, inds in filtered_tail_selections.items():
                if not inds:
                    continue
                f.write(
                    f"ct_{label}_phos: COORDINATION "
                    f"GROUPA={label}_ct GROUPB=phos_all R_0={htail_phos_cf}\n"
                )
                for i in range(1, 5):
                    f.write(
                        f"ct_{label}_core{i}: COORDINATION "
                        f"GROUPA={label}_ct GROUPB=core{i}_all R_0={htail_nucle_cf}\n"
                    )
            f.write("\n")

        # ------------------------------------------------------------------
        # rg: radius of gyration (non-Cu atoms (histone-tails))
        # ------------------------------------------------------------------
        if "rg" in features:
            f.write("# ---- Radius of gyration ----\n")
            rg_indices = top.select('not element Cu').tolist()
            if rg_indices:
                f.write(f"Rg: GYRATION ATOMS={','.join(map(str, rg_indices))}\n")
            f.write("\n")

        # ------------------------------------------------------------------
        # rmsd: RMSD to reference structures
        # ------------------------------------------------------------------
        if "rmsd" in features:
            f.write("# ---- RMSD to reference structures ----\n")
            if DNA_indices:
                f.write(f"RMSD_DNA: RMSD REFERENCE={os.path.join(ref_dir, 'DNA.pdb')} TYPE=OPTIMAL\n")
            f.write(f"RMSD_core_all: RMSD REFERENCE={os.path.join(ref_dir, 'core_all.pdb')} TYPE=OPTIMAL\n")
            f.write("\n")

        # ------------------------------------------------------------------
        # seg_dist: pairwise core-segment distances with exp(-x)
        # ------------------------------------------------------------------
        if "seg_dist" in features:
            f.write("# ---- Core-segment pairwise distances with exp(-x) ----\n")
            core_pairs = [
                ("core1", "core2"), ("core1", "core3"), ("core1", "core4"),
                ("core2", "core3"), ("core2", "core4"), ("core3", "core4"),
            ]
            for c1, c2 in core_pairs:
                for seg1 in core_segment_names[c1]:
                    for seg2 in core_segment_names[c2]:
                        d_label = (
                            f"d{c1[-1]}{c2[-1]}"
                            f"_{seg1.split('_')[-1]}"
                            f"_{seg2.split('_')[-1]}"
                        )
                        f.write(f"{d_label}: DISTANCE ATOMS={seg1},{seg2}\n")
                        f.write(f"e{d_label}: CUSTOM ARG={d_label} FUNC=exp(-x) PERIODIC=NO\n")
            f.write("\n")

        # ==================================================================
        # Output directive
        # ==================================================================
        f.write("PRINT ARG=* STRIDE=1 FILE=COLVAR\n")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    args = parse_args()
    features = resolve_features(args.features)

    tail_ranges = (
        load_tail_ranges(args.tail_ranges)
        if args.tail_ranges
        else DEFAULT_TAIL_RANGES
    )

    write_features(
        feature_file=args.output,
        pdb_file=args.pdb,
        features=features,
        tail_ranges=tail_ranges,
        linker_length=args.linker_length,
        n_segments_core=args.n_segments_core,
        segment_stride=args.segment_stride,
        htail_phos_cf=args.htail_phos_cf,
        htail_nucle_cf=args.htail_nucle_cf,
        ref_dir=args.ref_dir,
    )
    print(f"Written : {args.output}")
    print(f"Features: {', '.join(sorted(features))}")
    print(f"Linker  : {args.linker_length} bp  |  "
          f"Core segments: {args.n_segments_core} (stride {args.segment_stride})")


if __name__ == "__main__":
    main()
