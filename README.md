# 👑 Welcome to the RoyalMD Gold Edition
✨ Click on the image to watch the video in 4K:<br>
<a href="https://youtu.be/3tgJmtr9DHs"><img src="https://img.youtube.com/vi/3tgJmtr9DHs/maxresdefault.jpg" alt="Watch the video" width="800"></a>

## 💫 QUICK LAUNCH
```bash
python ./RoyalMD_Gold.py ./test_systems/CarbonicAnhydraseII.pdb --metalsurf \
       --interface-root ./INTERFACE_FF_1_5
```
This builds a gold surface, places the protein on top of it, and then runs the simulation. The metal surface is harmonically restrained during minimization, heating, and production, following the official AMBER tutorial:
[Setting Up A Protein System at the FCC Metal Surface](https://ambermd.org/tutorials/advanced/tutorial27/pro_metal.php).

The metal is modeled using the INTERFACE force field:
([INTERFACE-MD, v1.5](https://bionanostructures.com/interface-md/) · [Heinz *et al.*, *J. Phys.
Chem. C* **2008**, *112*, 17281](https://doi.org/10.1021/jp801931d)). It provides parameters for neutral fcc metals using 12-6 Lennard-Jones potentials. 
Download the force-field files once, then point `--interface-root` to the directory where you saved them:

```bash
curl -O https://bionanostructures.com/wp-content/uploads/2016/02/interface_ff_1_5.zip
unzip interface_ff_1_5.zip
```

For a complete list of command-line options:

```bash
python ./RoyalMD_Gold.py --help
```

## 🔭 Overview
**RoyalMD Gold Edition** is a lightweight molecular dynamics pipeline built on
**OpenMM**, designed to run **protein–metal interface** simulations on portable hardware. It
automates on the fly every step of surface modeling: from PDB fixing to building an fcc metal
slab, aligning the protein above it, solvation, minimization, multi-step equilibration and
production runs, with a single command-line interface. The metal surface follows the
**AMBER tutorial 27** route using the **INTERFACE force field** (Heinz et al.), but it is
built natively in OpenMM and runs solvated instead of in vacuo. It automatically detects your
hardware and uses GPU if available.

## 👤 Author

Original code & test: **Gleb Novikov**

Technical support: **Claude Opus**

## 🌟 How to cite

If you use **RoyalMD Gold Edition** in your projects, please cite it via the **Cite this repository**
(GitHub reads [CITATION.cff](CITATION.cff)), and cite the **INTERFACE force field**
it depends on: Heinz, Vaia, Farmer & Naik, *J. Phys. Chem. C* **2008**, *112*, 17281.


## ✨ Features
- Automatic **PDB fixing** (missing atoms, residues, protonation via `pdbfixer`)
- **fcc metal slab construction** from the INTERFACE model database (8 metals, 3 facets)
- **Automatic slab sizing**: replication factors derived from the protein footprint + padding
- **Principal-axis alignment** that lays the flattest face down and points the active site up
- **Solvation and ion placement** with box vectors locked to the metal lattice
- **Layer-selective restraints** — freeze the bulk, let the surface breathe
- Energy **minimization** and multi-step **NPT equilibration**
- Production MD runs with benchmark-ready logging
- Cross-platform optimization with GPU detection

---

## 🪄 Requirements
- Python 3.12
- OpenMM
- pdbfixer
- mdtraj (for saving MD trajectories in the NetCDF format)
- **AmberTools** (only `PropPDB` is used, to replicate the unit cell)
- **INTERFACE force field v1.5** — the metal model database

```bash
conda create -n MDsims python=3.12
conda activate MDsims
conda install -c conda-forge openmm pdbfixer mdtraj ambertools
conda install -c conda-forge "numpy<2.0" # <<< if you have any numpy-related errors"
```

Download and unpack the INTERFACE package (free, no registration):
```bash
curl -O https://bionanostructures.com/wp-content/uploads/2016/02/interface_ff_1_5.zip
unzip interface_ff_1_5.zip          # -> INTERFACE_FF_1_5/
```
Only `INTERFACE_FF_1_5/MODEL_DATABASE/METALS/` is needed. If it lives elsewhere, point
`--interface-root` at it.

## ⚜️ Usage Examples

```bash
# Carbonic anhydrase II on gold {111} — defaults, the whole slab restrained:
python ./RoyalMD_Gold.py ./test_systems/CarbonicAnhydraseII.pdb --metalsurf

# Let the top 2 metal layers breathe from NPT onwards (recommended):
python ./RoyalMD_Gold.py ./test_systems/CarbonicAnhydraseII.pdb --metalsurf \
       --metal-relax-top 2

# Silver instead of gold, {100} facet, thicker slab:
python ./RoyalMD_Gold.py ./test_systems/CarbonicAnhydraseII.pdb --metalsurf \
       --metal ag --metal-facet 100 --metal-cells-z 3

# Platinum, closer start, 50 ns production:
python ./RoyalMD_Gold.py ./test_systems/CarbonicAnhydraseII.pdb --metalsurf \
       --metal pt --metal-separation 2.5 --production-time 50000

# INTERFACE package somewhere else:
python ./RoyalMD_Gold.py ./test_systems/CarbonicAnhydraseII.pdb --metalsurf \
       --interface-root ~/software/INTERFACE_FF_1_5

# NO metal at all — plain protein MD in a dodecahedron:
python ./RoyalMD_Gold.py ./test_systems/CarbonicAnhydraseII.pdb
```

## ⚙️ Configuration

Every parameter is a command-line flag — nothing needs editing inside the script.

```bash
# --- files ---
--pre-sim-pdb       solvated.pdb   # solvated system written before the run
--output-nc         production.nc  # production trajectory

# --- integrator / timing ---
--timestep          0.002          # ps
--equil-time        1000           # ps of NPT equilibration (after 100 ps of heating)
--production-time   50000          # ps of production
--report-interval   5000           # save a frame every N steps
--log-freq          10             # print every N% of production

# --- conditions ---
--temperature       310            # K
--pressure          1.0            # bar
--ionic-strength    0.15           # mol/L
--pH                7.0            # used for protonation

# --- box and solvation ---
--box-type          dodecahedron   # dodecahedron | octahedron | cube  (ignored with --metalsurf)
--box-padding       1.2            # nm; ALSO sets the slab XY margin

# --- force fields ---
--ff-protein        amber14        # see the table below
--ff-water          tip3p

# --- cleanup ---
--strip-residues    SO4 HOH EDO GOL LIG LIH   # this build simulates the macromolecule only
```

## 🏆 Metal surface options

```bash
--metalsurf                        # master switch: build the slab and immobilize the solute
--interface-root  INTERFACE_FF_1_5 # path to the unpacked INTERFACE package
--metal           au               # ag | al | au | cu | ni | pb | pd | pt
--metal-facet     111              # 100 | 110 | 111
--metal-cells-z   2                # slab thickness in UNIT CELLS (fcc{111}: 3 layers/cell -> 6)
--metal-separation 3.0             # A between the solute's lowest atom and the top metal plane
--metal-headroom  25.0             # A of water above the solute
--metal-under     2.0              # A below the slab before the box floor
--metal-orient    (on)             # principal-axis alignment; --no-metal-orient to disable
--metal-relax-top 0                # release the top N layers once NPT begins
--metal-k-min     83680.0          # metal restraint during minimization (kJ/mol/nm²)
--metal-k-heat    41840.0          # ... during heating
--metal-k-prod    16736.0          # ... during NPT + production
```

### 🔬 How the surface is built

1. **Unit cell** — `<metal>_cell_P1_<facet>.car` is read from the INTERFACE model database and
   converted to PDB with **every metal atom as its own residue**, residue name identical to the
   atom name (`Au0`, `Ag0`, …). This is what keeps the slab loadable against the Amber-format
   INTERFACE parameters.
2. **Replication** — `PropPDB -ix -iy -iz` tiles the cell into a slab. **`-ix` and `-iy` are
   derived, never guessed**:

   ```
   ix = ceil((footprint_x + 2 × box-padding) / cell_a)      box_x = ix × cell_a
   ```

   The box edge is therefore an *exact integer multiple* of the lattice — the reason the slab
   stays continuous across the periodic boundary. `-iz` is yours to choose (`--metal-cells-z`),
   because thickness is a modelling decision, not a fit.
3. **Alignment** — the solute's inertia tensor is diagonalised; the axis with the **smallest
   extent** becomes the surface normal (flattest face down, shortest box). The sign is chosen by
   ray-casting from the centre of mass: the **least-obstructed direction** — the mouth of any
   active-site funnel — is pointed *away* from the metal, so the enzyme stays accessible.
4. **Solvation** — explicit `boxVectors`, not `padding`. Padding would pick an arbitrary box and
   tear the lattice.

### ⚛️ Force field for metals

INTERFACE fcc metals are **neutral, 12-6 Lennard-Jones only** — no charges, no bonded terms
(Heinz, Vaia, Farmer & Naik, *J. Phys. Chem. C* **2008**, *112*, 17281). Values are read
straight from the package and converted to OpenMM units:

| Metal | ε (kcal/mol) | Rmin/2 (Å) |
|---|---|---|
| Ag | 4.56 | 1.4775 |
| Al | 4.02 | 1.4625 |
| Au | 5.29 | 1.4755 |
| Cu | 4.72 | 1.3080 |
| Ni | 5.65 | 1.2760 |
| Pb | 2.93 | 1.7825 |
| Pd | 6.15 | 1.4095 |
| Pt | 7.80 | 1.4225 |

Since a 12-6 LJ has the same functional form in CHARMM and AMBER, and both use
Lorentz–Berthelot mixing, the table transfers 1:1. **Use it with AMBER** (`amber14`, `amber19`,
`amber99sb` are all verified). CHARMM36 currently refuses to merge because its 1-4 scaling factors
differ from the generated metal XML.

### 🔒 Restraints and the barostat

The metal gets **its own restraint force**, independent of the protein's. The protein's `k`
ramps 1000 → 0 across NPT as usual; the metal's does **not** — it is held all the way through
production, following AMBER tutorial 27 (`ntr=1`, `restraint_wt` 100 → 50 → 20 kcal/mol/Å²).

With `--metal-relax-top N` the top N layers are split off onto a third force that is driven to
zero the moment NPT begins:

| Phase | Bottom layers | Top N layers |
|---|---|---|
| Minimization | 83 680 | 83 680 — frozen |
| NVT heating | 41 840 | 41 840 — frozen |
| **NPT** | 16 736 | **0 — released** |
| **Production** | 16 736 | **0 — free** |

Why bother: uniform restraints give a *flat* vibrational profile, but real surfaces vibrate
~1.5–2× more than the bulk. Releasing the top 2 layers of a 6-layer Au{111} slab reproduces
that gradient (measured 0.092 Å anchored vs 0.133 Å at the surface) while the lower layers
still hold the lattice at its correct spacing.

> ⚠️ `--metal-relax-top 6` frees the **whole** slab. It stays a crystal, but with nothing
> anchored it drifts and the lattice sits ~3.5 % compressed, because x/y are pinned to the box.
> `2`–`4` is the better physics.

A **z-only anisotropic barostat** is used whenever a slab is present. Isotropic scaling would
rescale x/y — and therefore the metal's lattice constant itself — fighting the restraints and
moving the surface away from the spacing the INTERFACE parameters are defined for.

## 🌀 Supported force fields and water models

```bash
# Select with --ff-protein and --ff-water
  ff_map = {
        'amber19': 'amber19-all.xml', # with OPC (4 point water model)
        'amber14': 'amber14-all.xml', # with TIP3P (3 point water model)
        'amber99sb': 'amber99sb.xml', # with TIP3P water
        'amber99sbildn': 'amber99sbildn.xml', # with TIP3P water
        'amber03': 'amber03.xml', # with TIP3P water
        'charmm36': 'charmm36.xml' # with "water" model | ⚠️ does not work with --metalsurf
    }
```

## ✨ Test system

```bash
./test_systems/CarbonicAnhydraseII.pdb
```
**Human carbonic anhydrase II** (PDB **[3KS3](https://www.rcsb.org/structure/3KS3)**, X-ray at
**0.90 Å** — one of the highest-resolution protein structures available). A single 257-residue
chain dominated by a large, twisted **antiparallel β-sheet** that packs into a rigid, compact
fold, with only short helices around the edges. It also carries one Zn²⁺ ion, which is kept for the
tutorial and treated purely electrostatically (via a non-bonded model), plus one glycerol
and 481 crystallographic waters. The latter two are removed by the default `--strip-residues`.

## 🔮 Output

```
MD_<name>_<ddmmyy>/
├── solvated.pdb        # the built system, before minimization
├── minimized.cif       # after the three-step minimization
├── heating.nc          # NVT heating trajectory
├── equilibration.nc    # 10-stage NPT equilibration
├── equilibrated.cif    # final equilibrated snapshot
├── production.nc       # production trajectory
└── metal_param/        # generated unit cell, slab PDB and metal force field XML
```

## 📜 License

The code in this repository is released under the **MIT License** — see [LICENSE](LICENSE).
Copyright © TheVisualHub.

The INTERFACE force field v1.5 is not included in this repository. 
It is distributed under its own terms.
