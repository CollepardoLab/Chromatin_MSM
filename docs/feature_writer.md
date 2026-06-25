# multi_feature_writer.py

PLUMED collective variable (CV) input generator for tetranucleosomes, built on a coarse-grained (CG) model where bead types are encoded as element labels:

| Element | Bead type |
|---------|-----------|
| `Au`    | DNA backbone / phosphate beads |
| `Pt`    | Nucleosome core beads |
| `Cu`    | Histone tail beads |
| `P`     | Phosphorus atoms (for coordination contacts) |

The script reads a reference PDB, extracts atom indices, and writes a `multi_feature.dat` file ready for PLUMED

---

## Requirements

- Python 3.7+
- [MDTraj](https://mdtraj.org/) (`pip install mdtraj`)
- PLUMED 2.x (to run the generated `.dat` file)

---

## Quick start

```bash
# All features, linker length 15
python multi_feature_writer.py --pdb ref/15rep.pdb --linker-length 15

# Selected features only
python multi_feature_writer.py --pdb ref/15rep.pdb --linker-length 15 \
    --features distance angle torsion

# Custom cutoffs and output name
python multi_feature_writer.py --pdb ref/15rep.pdb \
    --features contacts rg \
    --htail-phos-cf 3.8 --htail-nucle-cf 3.6 \
    --output contact_rg.dat
```

---

## CLI options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--pdb` | path | *(required)* | Reference PDB file |
| `--output` | path | `multi_feature.dat` | Output PLUMED `.dat` file |
| `--features` | list | `all` | Feature types to include (see below) |
| `--linker-length` | int | `15` | Linker DNA length |
| `--n-segments-core` | int | `15` | Number of segments per nucleosome core |
| `--segment-stride` | int | `5` | Atom stride within each core-segment `CENTER` |
| `--htail-phos-cf` | float | `3.6` | Coordination cutoff: histone tail vs phosphate (nm) |
| `--htail-nucle-cf` | float | `3.6` | Coordination cutoff: histone tail vs core (nm) |
| `--tail-ranges` | path | built-in | JSON file with histone tail residue ranges (see below) |

---

## Feature types

Pass one or more names to `--features` (space-separated), or omit for all.

### `distance`
Pairwise distances between the 4 nucleosome core centers (6 pairs: 1-2, 1-3, 1-4, 2-3, 2-4, 3-4).

```
dist_core1_2, dist_core1_3, ..., dist_core3_4
```

### `angle`
Triple angles along the nucleosome chain, with sin/cos encoding to avoid periodicity issues.

```
angle_123, angle_234
sin_angle_123, cos_angle_123, sin_angle_234, cos_angle_234
```

### `torsion`
Torsion angle defined by all 4 core centers, with sin/cos encoding.

```
psi_nucl, sin_psi_nucl, cos_psi_nucl
```

### `bend`
DNA linker bending angle for each of the 3 inter-nucleosome linkers, with sin/cos encoding. Uses atoms at ±63 beads from each dyad plus a midpoint virtual atom.

```
bend_A, bend_B, bend_C
sin_bend_A, cos_bend_A, ...
```

### `contacts`
Histone-tail / DNA-phosphate / core coordination numbers via PLUMED `COORDINATION`, plus center-of-mass (COM) positions for pK, pR amino acids in each tail group and per-nucleosome phosphate group.

```
ct_tail1_A_phos, ct_tail1_A_core1, ...   (coordination numbers)
tail1_A_ct_com, ...                       (tail COMs)
phosA_com, phosB_com, ...                 (phosphate COMs)
```

### `rg`
Radius of gyration computed over all non-histone-tail beads (`not element Cu`).

```
Rg
```

### `rmsd`
RMSD to reference structures. Requires `./ref/DNA.pdb` and `./ref/core_all.pdb` relative to the working directory where PLUMED is run.

```
RMSD_DNA, RMSD_core_all
```

### `seg_dist`
Each nucleosome core is divided into `--n-segments-core` segments. All pairwise inter-core segment distances are computed with an `exp(-x)` transformation for smooth encoding (6 core pairs × segments²).

```
d12_seg0_seg0, ed12_seg0_seg0, ...
```

---

## Histone tail ranges

By default the script uses built-in ranges for a standard 4-nucleosome system. To use a different system, supply a JSON file:

```bash
python multi_feature_writer.py --pdb ref/rep.pdb --tail-ranges tail_ranges.json
```

**JSON format** — keys are histone chain labels (`A`–`D`), values are lists of `[start, end]` residue index pairs:

```json
{
  "A": [[1, 40], [136, 160], [238, 257]],
  "B": [[975, 1013], [1110, 1134]],
  "C": [[1949, 1987]],
  "D": [[2923, 2962]]
}
```

Only residues with names `pR` or `pK` within these ranges are included in contact group definitions.

---

## Output file structure

The generated `.dat` file is organised into two sections:

1. **Shared prerequisites** — `GROUP` and `CENTER` virtual atoms defined once, in dependency order. Only blocks required by the selected features are written.
2. **Feature CVs** — one clearly labelled section per feature type.

A `PRINT ARG=* STRIDE=1 FILE=COLVAR` directive is appended at the end.

### Example output (distance + torsion)

```
MOLINFO STRUCTURE=ref/29rep.pdb

# Core center virtual atoms
core1_center: CENTER ATOMS=0,1,2,...
...

# ---- Core-center pairwise distances ----
dist_core1_2: DISTANCE ATOMS=core1_center,core2_center
...

# ---- Core torsion angle ----
psi_nucl: TORSION ATOMS=core1_center,core2_center,core3_center,core4_center
sin_psi_nucl: CUSTOM ARG=psi_nucl FUNC=sin(x) PERIODIC=NO
cos_psi_nucl: CUSTOM ARG=psi_nucl FUNC=cos(x) PERIODIC=NO

PRINT ARG=* STRIDE=1 FILE=COLVAR
```

---

## Running with PLUMED

```bash
plumed driver --plumed multi_feature.dat --mf_xtc traj.xtc
```

Or add `multi_feature.dat` as the PLUMED input in your GROMACS `.mdp`:

```
plumed-active = yes
plumed-dat    = multi_feature.dat
```

---

## Notes

- PLUMED is sensitive to **duplicate label definitions** — this script ensures every `GROUP`/`CENTER` is written exactly once regardless of which features are selected.
- The `rmsd` feature expects reference PDB files at `./ref/DNA.pdb` and `./ref/core_all.pdb` relative to where PLUMED is run, not where this script is run.
- `exp(-x)` is used for the `seg_dist` exponential encoding.
