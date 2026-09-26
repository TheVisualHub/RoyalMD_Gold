# 👑 Royal MD — GOLD EDITION: MOLECULAR DYNAMICS on your portable hardware
# Created by Gleb Novikov
# last update 19/09/2026:
# v1.6p:
# - NEW --metal-relax-top N: the top N metal layers stay frozen through minimization
#   and NVT, then are RELEASED when NPT begins and stay free in production, so the
#   surface can breathe against the solute while the lower layers anchor the lattice.
#   Implemented as a third CustomExternalForce with its own global parameter
#   (k_metal_top -> 0 at NPT). Layers are detected by clustering metal z-coordinates,
#   so it works for any metal, facet and thickness. Default 0 = v1.5 behaviour.
#   Recovers the surface-vs-bulk vibration difference that uniform restraints suppress.
#
# (c) The VisualHub 2026
import os
import sys
import glob
import math
import time
import shutil
import random
import argparse
import subprocess
import numpy as np
from openmm.app import *
from openmm import *
from openmm.unit import *
from pdbfixer import PDBFixer
from datetime import datetime

# Registry of supported force fields (selectable via --ff-protein).
# internal; the selectable value is exposed through --ff-protein / --ff-water.
FF_MAP = {
    'amber19': 'amber19-all.xml',       # with OPC water
    'amber14': 'amber14-all.xml',       # with TIP3P water
    'amber99sb': 'amber99sb.xml',
    'amber99sbildn': 'amber99sbildn.xml',
    'amber03': 'amber03.xml',
    'charmm36': 'charmm36.xml'
}

# --- INTERFACE FF v1.5 — neutral FCC metals, 12-6 Lennard-Jones ---------------
# Heinz, Vaia, Farmer & Naik, J. Phys. Chem. C 2008, 112, 17281.
# Values are transcribed VERBATIM from the package file
#     INTERFACE_FF_1_5/FORCE_FIELDS/charmm27_interface_v1_5.prm
# whose NONBONDED format is:  <type>  0.0  -epsilon(kcal/mol)  Rmin/2(A)
# CHARMM and AMBER share this convention, so the numbers transfer 1:1. This was
# verified by rebuilding AMBER tutorial 27's Ag system from these values: all 40
# parameter-bearing prmtop blocks matched the published reference byte-for-byte.
# q = 0 for every metal: the INTERFACE fcc metals are strictly neutral.
METAL_LJ = {
    #  symbol : (epsilon kcal/mol, Rmin/2 A, atomic mass)
    'ag': (4.56, 1.4775, 107.8682),
    'al': (4.02, 1.4625,  26.9815385),
    'au': (5.29, 1.4755, 196.96657),
    'cu': (4.72, 1.3080,  63.546),
    'ni': (5.65, 1.2760,  58.6934),
    'pb': (2.93, 1.7825, 207.2),
    'pd': (6.15, 1.4095, 106.42),
    'pt': (7.80, 1.4225, 195.084),
}
KCAL_TO_KJ = 4.184

def print_intro_banner(add_slogan_separator=True):
    """
    Prints an inverted tetrahedral brand logo.
    Geometry: Inverted sp3 tripod (top-heavy) with wide-angled base.
    Palette: Champagne Gold Gradient.
    """
    # Champagne Gold Palette
    GOLD_D  = "\033[38;5;172m" # Deep Gold (Bonds/Springs)
    GOLD_M  = "\033[38;5;220m" # Bright Gold (Atom Borders)
    GOLD_L  = "\033[38;5;221m" # Champagne Gold (Atom Centers)
    
    BOLD    = "\033[1m"
    RESET   = "\033[0m"

    # Inverted geometry: One top atom, spreading into three wide legs (no center foot)
    logo = rf"""
{GOLD_M}{BOLD}
               Welcome to the Royal MD!

                 {GOLD_L}▟██▙{RESET}    {GOLD_L}▟██▙{RESET}    {GOLD_L}▟██▙{RESET}
                 {GOLD_L}▜██▛{RESET}    {GOLD_L}▜██▛{RESET}    {GOLD_L}▜██▛{RESET}
                  {GOLD_D}▜█▙{RESET}     {GOLD_D}█{RESET}     {GOLD_D}▟█▛{RESET}
                   {GOLD_M}▜██▙{GOLD_D}  {GOLD_M}▟█▙{GOLD_D}  {GOLD_M}▟██▛{RESET}
                    {GOLD_M}▜█▙{GOLD_D}▀ {GOLD_M}▜█▛{GOLD_D} ▀{GOLD_M}▜█▛{RESET}
                          {GOLD_D}█{RESET}
                          {GOLD_D}█{RESET}
                        {GOLD_L}▟██▙{RESET}
                        {GOLD_L}▜██▛{RESET}
{RESET}
    """
    print(logo)
    time.sleep(0.5)
    slogan = "GOLD EDITION"
    separators = ["🥂", "👑", "🔱", "✨"]
    available_seps = separators.copy()
    print("          ", end="")
    time.sleep(1.0)
    if add_slogan_separator:
        # two random separators before and two after - never repeated,
        # so every launch draws a different crest
        picks = random.sample(available_seps, 4)
        pre, post = picks[:2], picks[2:]
        line = f"{BOLD}{GOLD_M}{' '.join(pre)}  {slogan}  {' '.join(post)}{RESET}"
    else:
        line = f"{BOLD}{GOLD_M}    {slogan}    {RESET}"
    # Close the line (reset + newline) BEFORE the pause: otherwise the cursor
    # block parks right next to the final emoji and shows as a stray square.
    print(line, flush=True)
    time.sleep(2.0)
    print()  # blank line, as before

# --- STEP 0: CLEAN OLD DATA ---
def clean_old_files(target_pdb, activate: bool = False):
    if not activate:
        return
    print(f"🧹Cleaning Old Data:")
    patterns = ["*.dcd", "*.nc", target_pdb, "*.cif"]
    for pattern in patterns:
        for file in glob.glob(pattern):
            try:
                if os.path.isfile(file):
                    os.remove(file)
                print(f"Removed: {file}")
            except Exception as e:
                pass

def clean_old_dir(work_dir, activate: bool = False):
    if not activate:
        return
    if os.path.exists(work_dir):
        print(f"🧹 Cleaning Old Data: Removing {work_dir}")
        shutil.rmtree(work_dir)
    os.makedirs(work_dir, exist_ok=True)

# --- STEP 1: FIX PROTEIN TOPOLOGY ---
def fix_protein_topology(input_pdb, env_pH, strip_residues):
    print(f"--- Step 1: Cleaning Protein Topology ---")
    fixer = PDBFixer(filename=input_pdb)

    # 1. STRIP UNWANTED RESIDUES FIRST (stops the NoneType error)
    modeller = Modeller(fixer.topology, fixer.positions)
    residues_to_strip = [r for r in modeller.topology.residues() if r.name in strip_residues]
    
    if residues_to_strip:
        removed_names = {r.name for r in residues_to_strip}
        print(f"🧹 Pre-cleaning: Removing {len(residues_to_strip)} residues {removed_names}")
        modeller.delete(residues_to_strip)
        
        # Synchronize fixer so it never "sees" the stripped residues
        fixer.topology = modeller.topology
        fixer.positions = modeller.positions
        # 👑 THE FIX: Clear the "memory" of the old gaps
        fixer.missingResidues = {}
    else:
        print("✨ No blacklisted residues found to strip.")

    # 2. Convert MSE to MET
    mse_count = 0
    cys_count = 0
    for residue in fixer.topology.residues():
        if residue.name == 'MSE':
            residue.name = 'MET'
            mse_count += 1
            for atom in residue.atoms():
                if atom.element.symbol == 'Se':
                    atom.element = Element.getBySymbol('S')
                    atom.name = 'SD' 
        elif residue.name in ['CYM', 'CYX']:
            residue.name = 'CYS'
            cys_count += 1


    if mse_count > 0:
        print(f"🧬 Converted {mse_count} 'MSE' residues to 'MET'.")
    if cys_count > 0:
        print(f"🧪 Standardized {cys_count} 'CYM/CYX' residues to 'CYS'.")

    # 3. Analyze and Repair (Only on clean protein)
    fixer.findMissingResidues()
    # added 29/12: quick GAP diagnostics (currently shows only gaps in the last chain ;-))
    if fixer.missingResidues:
        # Use [ ] to create a list explicitly, satisfying OpenMM's sum() 
        total_missing = sum([len(residues) for residues in fixer.missingResidues.values()])
        print(f"💫 GAPS DETECTED: Total of {total_missing} missing residues found.")
        
        chains = list(fixer.topology.chains())
        num_chains = len(chains)
        
        # Sort keys to ensure Chain A comes before Chain B
        sorted_keys = sorted(fixer.missingResidues.keys())
        
        for key in sorted_keys:
            chain_idx, res_idx = key
            residues = fixer.missingResidues[key]
            
            if chain_idx < num_chains:
                chain_id = chains[chain_idx].id
            else:
                chain_id = f"Unknown({chain_idx})"
            
            gap_size = len(residues)
            real_resid = res_idx +1
            print(f"  - Chain {chain_id}: Gap at residue {real_resid} ({gap_size} residues: {', '.join(residues)})")
    else:
        print("📜 NO GAPS DETECTED! The STRUCTURE is OKAY 👌")
    
    unique_residues_kept = {r.name for r in fixer.topology.residues()}
    print(f"⛩️ Residues kept: {sorted(list(unique_residues_kept))}")

    print(f"⚛️ Filling missing atoms and hydrogens at pH '{env_pH}'...")
    fixer.findMissingAtoms()
    fixer.addMissingAtoms() # Will not crash now
    fixer.addMissingHydrogens(env_pH)

    return fixer


# =============================================================================
#  STEP 0.7 (GOLD EDITION): METAL SURFACE CONSTRUCTION
#  Follows AMBER tutorial 27 (Pengfei Li): INTERFACE unit cell -> PDB with every
#  metal atom its own residue (resname == atomname) -> PropPDB replication.
#  The VMD hand-placement step is replaced by programmatic principal-axis
#  alignment, and the in vacuo sander run by solvated OpenMM.
# =============================================================================

def locate_interface_cell(interface_root, metal, facet):
    """
    Find <metal>_cell_P1_<facet>.car inside an extracted INTERFACE_FF_1_5 tree.
    Accepts either the package root or the METALS directory itself.
    """
    metal, facet = metal.lower(), str(facet)
    candidates = [
        os.path.join(interface_root, 'MODEL_DATABASE', 'METALS'),
        os.path.join(interface_root, 'INTERFACE_FF_1_5', 'MODEL_DATABASE', 'METALS'),
        interface_root,
    ]
    fname = f"{metal}_cell_P1_{facet}.car"
    for d in candidates:
        p = os.path.join(d, fname)
        if os.path.isfile(p):
            return p
    raise FileNotFoundError(
        f"❌ Could not find '{fname}'.\n"
        f"   Looked under: {', '.join(candidates)}\n"
        f"   Point --interface-root at your extracted INTERFACE_FF_1_5 folder.\n"
        f"   The package (v1.5) is at https://bionanostructures.com/interface-md/")


def car_to_pdb(car_path, pdb_path, metal):
    """
    Convert an INTERFACE Materials-Studio .car cell to an Amber-compatible PDB.

    Tutorial convention: EVERY metal atom becomes its own residue, and the
    residue name is identical to the atom name (e.g. 'Au0'). That is what keeps
    the structure loadable against the Amber-format INTERFACE parameters.
    Returns (n_atoms, (a, b, c)) with the cell edges in Angstrom.
    """
    name = f"{metal.capitalize()}0"           # Au -> Au0
    cell, atoms = None, []
    with open(car_path) as fh:
        for line in fh:
            f = line.split()
            if line.startswith('PBC ') and len(f) >= 7:
                cell = [float(x) for x in f[1:7]]
            elif f and f[0].lower().startswith(metal.lower()) and len(f) >= 9 and f[4] == 'XXXX':
                atoms.append((float(f[1]), float(f[2]), float(f[3])))
    if cell is None or not atoms:
        raise ValueError(f"❌ Could not parse cell/atoms from {car_path}")
    a, b, c, al, be, ga = cell
    with open(pdb_path, 'w') as fh:
        fh.write(f"CRYST1{a:9.3f}{b:9.3f}{c:9.3f}{al:7.2f}{be:7.2f}{ga:7.2f} P 1           1\n")
        for i, (x, y, z) in enumerate(atoms, 1):
            fh.write(f"ATOM  {i:5d}  {name:<3s} {name:<3s} A{i:4d}    "
                     f"{x:8.3f}{y:8.3f}{z:8.3f}{0.0:6.2f}{0.0:6.2f}"
                     f"{'':10s}{metal.upper():>2s}\n")
        fh.write("END\n")
    return len(atoms), (a, b, c)


def run_proppdb(cell_pdb, slab_pdb, ix, iy, iz):
    """Replicate the unit cell into a slab with AmberTools' PropPDB."""
    if shutil.which('PropPDB') is None:
        raise RuntimeError("❌ PropPDB not found on PATH (AmberTools required for --metalsurf).")
    subprocess.run(['PropPDB', '-p', cell_pdb, '-o', slab_pdb,
                    '-ix', str(ix), '-iy', str(iy), '-iz', str(iz)], check=True)
    coords = []
    with open(slab_pdb) as fh:
        for l in fh:
            if l.startswith(('ATOM', 'HETATM')):
                coords.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    return np.array(coords)


def write_metal_forcefield(xml_path, metal):
    """
    Emit an OpenMM ForceField XML for the INTERFACE metal, converting the
    distributed CHARMM/AMBER-convention values into OpenMM units:
        sigma   = 2*(Rmin/2) / 2^(1/6)   [A -> nm]
        epsilon = eps_kcal * 4.184       [kcal/mol -> kJ/mol]
    The metal is neutral and has no bonded terms, so one <Atom> is the whole model.
    """
    eps_kcal, rmin_half, mass = METAL_LJ[metal.lower()]
    name = f"{metal.capitalize()}0"
    sigma_nm = (2.0 * rmin_half / (2.0 ** (1.0 / 6.0))) / 10.0
    eps_kj = eps_kcal * KCAL_TO_KJ
    os.makedirs(os.path.dirname(xml_path) or '.', exist_ok=True)
    with open(xml_path, 'w') as fh:
        fh.write(f"""<ForceField>
 <!-- INTERFACE FF v1.5 neutral fcc metal: {name}
      Heinz, Vaia, Farmer, Naik, J. Phys. Chem. C 2008, 112, 17281 (12-6 LJ).
      Package value (charmm27_interface_v1_5.prm): eps = {eps_kcal} kcal/mol,
      Rmin/2 = {rmin_half} A, q = 0.  Converted only, never refitted. -->
 <AtomTypes>
  <Type name="{name}" class="{name}" element="{metal.capitalize()}" mass="{mass}"/>
 </AtomTypes>
 <Residues>
  <Residue name="{name}">
   <Atom name="{name}" type="{name}" charge="0.0"/>
  </Residue>
 </Residues>
 <NonbondedForce coulomb14scale="0.8333333333333334" lj14scale="0.5">
  <UseAttributeFromResidue name="charge"/>
  <Atom type="{name}" sigma="{sigma_nm:.8f}" epsilon="{eps_kj:.6f}"/>
 </NonbondedForce>
</ForceField>
""")
    print(f"🔬 {name} LJ from INTERFACE: eps = {eps_kcal} kcal/mol -> {eps_kj:.5f} kJ/mol, "
          f"Rmin/2 = {rmin_half} A -> sigma = {sigma_nm*10:.5f} A, q = 0")
    return xml_path, name


def principal_axis_rotation(positions_A, masses, elements, orient=True):
    """
    Rotation that lays the solute on the surface: the principal axis with the
    SMALLEST extent is mapped to +z, so the flattest aspect faces the metal and
    the box stays short in z.

    The sign is chosen so the largest solvent-exposed pocket (the direction from
    the solute centre with the least atom density, i.e. the mouth of any active
    site) points AWAY from the surface. For an enzyme this is what keeps the
    catalytic site accessible after immobilization.

    Returns (R, com); apply as  new = (old - com) @ R.T
    """
    heavy = np.array([e != 'H' for e in elements])
    P, m = positions_A[heavy], masses[heavy]
    com = (P * m[:, None]).sum(0) / m.sum()
    if not orient:
        return np.eye(3), com
    X = P - com
    I = np.zeros((3, 3))
    for xi, mi in zip(X, m):
        I += mi * (np.dot(xi, xi) * np.eye(3) - np.outer(xi, xi))
    ev, evec = np.linalg.eigh(I)
    evec = evec[:, np.argsort(ev)]                  # col0 longest ... col2 shortest extent
    if np.linalg.det(evec) < 0:
        evec[:, 2] *= -1
    Q = X @ evec
    ext = Q.max(0) - Q.min(0)
    k = int(np.argmin(ext))                         # shortest extent -> surface normal

    # least-obstructed direction from the centre = solvent-exposed pocket
    n = 2000
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    th = np.pi * (1 + 5 ** 0.5) * i
    D = np.c_[np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)]
    d = np.linalg.norm(Q, axis=1)
    sel = d < 25.0
    if sel.sum() > 10:
        cosang = (Q[sel] @ D.T) / d[sel][:, None]
        pocket = D[np.argmin((cosang > np.cos(np.radians(30))).sum(0))]
        s = 1.0 if pocket[k] > 0 else -1.0          # pocket should end up pointing +z
    else:
        s = 1.0
    e3 = s * evec[:, k]
    rest = [j for j in range(3) if j != k]
    e1 = evec[:, rest[0]]
    e1 = e1 - np.dot(e1, e3) * e3
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(e3, e1)
    R = np.vstack([e1, e2, e3])
    if np.linalg.det(R) < 0:
        R[1] *= -1
    return R, com


def build_metal_slab(modeller, metal_cfg, work_dir):
    """
    Orient the solute, size the slab from its footprint + --box-padding, replicate
    it with PropPDB, drop it under the solute and return the periodic box vectors.

    The XY replication is DERIVED, never guessed: the box edge must be an exact
    integer multiple of the metal cell edge or the slab will not be continuous
    across the periodic boundary, so
            ix = ceil((footprint_x + 2*padding) / cell_a)
    and the box vector is then exactly ix*cell_a. Z thickness is manual
    (--metal-cells-z) because it is a physical modelling choice, not a fit.
    """
    print(f"🏗️  Step 0.7: Building {metal_cfg['metal'].capitalize()}{{{metal_cfg['facet']}}} surface ...")
    param_dir = os.path.join(work_dir, 'metal_param')
    os.makedirs(param_dir, exist_ok=True)

    car = locate_interface_cell(metal_cfg['interface_root'], metal_cfg['metal'], metal_cfg['facet'])
    cell_pdb = os.path.join(param_dir, f"{metal_cfg['metal']}_cell_P1_{metal_cfg['facet']}.pdb")
    n_cell, (ca, cb, cc) = car_to_pdb(car, cell_pdb, metal_cfg['metal'])
    print(f"   unit cell: {os.path.basename(car)} -> {n_cell} atoms, "
          f"{ca:.4f} x {cb:.4f} x {cc:.4f} A")

    pos_A = np.array(modeller.positions.value_in_unit(angstrom))
    elements = [a.element.symbol if a.element is not None else 'X'
                for a in modeller.topology.atoms()]
    masses = np.array([a.element.mass.value_in_unit(dalton) if a.element is not None else 12.0
                       for a in modeller.topology.atoms()])
    R, com = principal_axis_rotation(pos_A, masses, elements, orient=metal_cfg['orient'])
    P = (pos_A - com) @ R.T
    fx, fy = P[:, 0].ptp(), P[:, 1].ptp()
    pad_A = metal_cfg['padding'] * 10.0
    ix = max(1, int(math.ceil((fx + 2 * pad_A) / ca)))
    iy = max(1, int(math.ceil((fy + 2 * pad_A) / cb)))
    iz = metal_cfg['cells_z']
    LX, LY = ix * ca, iy * cb
    print(f"   solute footprint {fx:.1f} x {fy:.1f} A + 2 x {pad_A:.1f} A padding")
    print(f"   -> PropPDB -ix {ix} -iy {iy} -iz {iz}   (box locked to lattice: "
          f"{LX:.4f} x {LY:.4f} A)")

    slab_pdb = os.path.join(param_dir, f"{metal_cfg['metal']}_slab.pdb")
    S = run_proppdb(cell_pdb, slab_pdb, ix, iy, iz)
    n_layers = len(np.unique(np.round(S[:, 2], 3)))
    print(f"   slab: {len(S)} atoms, {n_layers} atomic layers, "
          f"{S[:,0].ptp():.2f} x {S[:,1].ptp():.2f} x {S[:,2].ptp():.2f} A")

    # centre slab in XY on the box, place its TOP face --metal-separation below the solute
    S[:, 0] -= LX / 2.0
    S[:, 1] -= LY / 2.0
    S[:, 2] += (P[:, 2].min() - metal_cfg['separation']) - S[:, 2].max()

    # shift everything into a box starting at the origin
    zlo = S[:, 2].min() - metal_cfg['under']
    zhi = P[:, 2].max() + metal_cfg['headroom']
    LZ = zhi - zlo
    shift = np.array([LX / 2.0, LY / 2.0, -zlo])
    P += shift
    S += shift

    for i, p in enumerate(P):
        modeller.positions[i] = Vec3(*p) * angstrom
    chain = modeller.topology.addChain('S')
    el = Element.getBySymbol(metal_cfg['metal'].capitalize())
    resname = f"{metal_cfg['metal'].capitalize()}0"
    for i, xyz in enumerate(S, 1):
        r = modeller.topology.addResidue(resname, chain, str(i))
        modeller.topology.addAtom(resname, el, r)
        modeller.positions.append(Vec3(*xyz) * angstrom)

    box = (Vec3(LX, 0, 0), Vec3(0, LY, 0), Vec3(0, 0, LZ)) * angstrom
    modeller.topology.setPeriodicBoxVectors(box)
    print(f"   separation {metal_cfg['separation']:.1f} A | box {LX:.3f} x {LY:.3f} x {LZ:.3f} A")
    return box, resname, len(S)

# --- STEP 2: MERGE AND SOLVATE -----------------------------------------------
def solvate_system(protein_fixer, ff_protein, ff_water, ff_map, box_type,
                   box_padding, ionic_strength, pre_sim_pdb, input_pdb,
                   metal_cfg=None, work_dir='.'):
    print(f"-- Step 2: Solvating 🫧 --")

    # 1. Start with the fixed protein
    modeller = Modeller(protein_fixer.topology, protein_fixer.positions)

    # 3. Load Forcefields
    protein_xml = ff_map[ff_protein]
    try:
        water_xml = f"{ff_protein}/{ff_water}.xml"
        # added: test the ForceField creation the valided path
        _checkFF = ForceField(protein_xml, water_xml)
        ff_args = [protein_xml, water_xml]
    except:
        water_xml = f"{ff_water}.xml"
        ff_args = [protein_xml, water_xml]

    # --- GOLD EDITION: metal XML joins the SAME ForceField, no second pipeline ---
    metal_resname, n_metal, box_vectors = None, 0, None
    if metal_cfg is not None:
        metal_xml, metal_resname = write_metal_forcefield(
            os.path.join(work_dir, 'metal_param', 'metal.xml'), metal_cfg['metal'])
        ff_args.append(metal_xml)

    print(f"🥂 Successfully loading: {ff_args}")
    forcefield = ForceField(*ff_args)

    # The slab is built AFTER the solute topology is final, so everything rotates
    # together with the protein.
    if metal_cfg is not None:
        box_vectors, metal_resname, n_metal = build_metal_slab(modeller, metal_cfg, work_dir)

    # --- FIX FOR 4-POINT WATER (like OPC with amber19) ---
    # .. check if we need the 4-point 'tip4pew' generator
    solvent_model = ff_water
    if 'opc' in ff_water.lower() or 'tip4p' in ff_water.lower():
        print(f"💧 4-point water detected ({ff_water}). Setting geometry to 'tip4pew'.")
        solvent_model = 'tip4pew'
    elif 'tip3p' in ff_water.lower() or ff_water.lower() == 'water':
        # CHARMM 'water' or standard 'tip3p' both use tip3p geometry
        print(f"💧 3-point water detected ({ff_water}). Setting geometry to 'tip3p'.")
        solvent_model = 'tip3p'
    elif 'spce' in ff_water.lower():
        print(f"💧 SPC/E water detected ({ff_water}). Setting geometry to 'spce'.")
        solvent_model = 'spce'
    else:
        print(f"😅 Unknown water '{ff_water}', defaulting to 'tip3p' geometry.")
        solvent_model = 'tip3p'

    # 4. Solvate
    print("💦 Adding solvent...")
    res_templates = {}
    if metal_cfg is not None:
        # amber14 ships MONATOMIC ION templates for Ag, Al, Cu, Ni, Pb, Pd and Pt.
        # A lone slab atom matches both that ion and our <Metal>0 template, and
        # OpenMM refuses the ambiguity ("Multiple non-identical matching templates").
        # Gold is the only one of the eight with no amber14 counterpart, so this bites
        # every metal except Au. Naming the template explicitly resolves it.
        res_templates = {r: metal_resname for r in modeller.topology.residues()
                         if r.name == metal_resname}
        # Explicit box vectors, NOT padding/boxShape: the XY edges must stay exact
        # integer multiples of the metal lattice or the slab tears at the boundary.
        modeller.addSolvent(forcefield,
                            boxVectors=box_vectors,
                            model=solvent_model,
                            ionicStrength=ionic_strength*molar,
                            residueTemplates=res_templates)
    else:
        modeller.addSolvent(forcefield,
                            padding=box_padding*nanometer,
                            model=solvent_model,
                            boxShape=box_type,
                            ionicStrength=ionic_strength*molar)

    print(f"Total System Size: {modeller.topology.getNumAtoms()} atoms.")
    if metal_cfg is not None:
        counts = {}
        for r in modeller.topology.residues():
            counts[r.name] = counts.get(r.name, 0) + 1
        nw = counts.get('HOH', 0)
        print(f"🏅 Metal atoms ({metal_resname}): {n_metal} | 💧 waters: {nw} | "
              f"Na+: {counts.get('NA',0)} | Cl-: {counts.get('CL',0)}")
        bv = modeller.topology.getPeriodicBoxVectors().value_in_unit(angstrom)
        print(f"📦 Box vectors (A): {bv[0][0]:.4f} x {bv[1][1]:.4f} x {bv[2][2]:.4f}")

    # addSolvent rebuilds the topology, so re-key the map onto the NEW residues
    # for the later createSystem call.
    if metal_cfg is not None:
        res_templates = {r: metal_resname for r in modeller.topology.residues()
                         if r.name == metal_resname}

    with open(pre_sim_pdb, 'w') as f:
        PDBFile.writeFile(modeller.topology, modeller.positions, f)

    return modeller, forcefield, metal_resname, res_templates
# --- MODULAR ADD-ON 2: REFRESH RESTRAINT REFERENCES (anti-NaN) ---
def refresh_restraint_references(simulation, forces):
    """
    Reset every restraint's reference coords (x0,y0,z0) to the CURRENT
    minimized positions, BEFORE NVT heating turns k back up to k_max.

    Why this is needed here: the references are captured from the raw
    pre-minimization (solvated) coordinates. Minimization step 3 runs with
    k=0.0, so restrained atoms relax freely and drift off those references.
    Heating then sets k=k_max against the STALE refs -> F=k*d yanks every
    displaced atom hard -> energy spike -> "Particle coordinate is NaN" on
    the first step. Re-referencing to the minimized structure makes the
    heating restraints hold atoms where they actually are.
    """
    state = simulation.context.getState(getPositions=True)
    pos = state.getPositions().value_in_unit(nanometer)
    if not isinstance(forces, (list, tuple)):
        forces = [forces]
    for force in forces:
        for i in range(force.getNumParticles()):
            atom_index, _params = force.getParticleParameters(i)
            p = pos[atom_index]
            force.setParticleParameters(i, atom_index, [p[0], p[1], p[2]])
        force.updateParametersInContext(simulation.context)

# --- STEPS 3 - 5: MINIMIZE, EQUILIBRATE AND RUN ---
def setup_and_run(modeller, forcefield, temperature, pressure, timestep, equil_time, production_time, report_interval, output_nc, log_freq, metal_resname=None, metal_cfg=None,
                  res_templates=None):
    print("--- Step 3: Setup & Minimize ⚡ ---")

    # Force conversion to numbers
    timestep = float(timestep)
    report_interval = int(report_interval)
    # production run steps:
    total_steps = int(production_time / timestep)

    # Create the System, Integrator, and Simulation
    system = forcefield.createSystem(modeller.topology,
                                    residueTemplates=(res_templates or {}), 
                                    nonbondedMethod=PME, 
                                    nonbondedCutoff=1.0*nanometer, 
                                    constraints=HBonds,
                                    rigidWater=True)

    if metal_resname is not None:
        # Anisotropic, Z-ONLY. Isotropic scaling would rescale X/Y and therefore the
        # metal lattice constant itself, fighting the position restraints and moving
        # the surface away from the spacing the INTERFACE parameters are defined for.
        barostat = MonteCarloAnisotropicBarostat(
            Vec3(pressure, pressure, pressure)*bar, temperature*kelvin,
            False, False, True)
        print("🧿 Anisotropic barostat (Z only): metal lattice in X/Y is never rescaled.")
    else:
        barostat = MonteCarloBarostat(pressure*bar, temperature*kelvin)
    system.addForce(barostat)

    # 3. Multi-step NPT equilibration with Harmonic Restraints
    equil_time = int(equil_time)
    num_stages = 10
    time_per_stage = equil_time / num_stages # for print report only
    k_max = 1000.0

    #formula = "0.5 * k * ((x-x0)^2 + (y-y0)^2 + (z-z0)^2)" # this may work only for vacum simulations
    formula = "0.5 * k * periodicdistance(x, y, z, x0, y0, z0)^2" # this works with PBC
    restraint_force = CustomExternalForce(formula)
    restraint_force.addGlobalParameter("k", k_max) 
    restraint_force.addPerParticleParameter("x0")
    restraint_force.addPerParticleParameter("y0")
    restraint_force.addPerParticleParameter("z0")

    # --- GOLD EDITION: the metal gets its OWN force with its OWN global parameter.
    # RoyalMD ramps 'k' to zero for production; the slab must stay restrained the
    # whole way through or the cluster falls apart (AMBER tutorial 27 keeps ntr=1
    # with restraint_wt=20 kcal/mol/A^2 even in production).
    metal_force = None
    if metal_resname is not None:
        metal_force = CustomExternalForce(formula.replace(" k ", " k_metal "))
        metal_force.addGlobalParameter("k_metal", metal_cfg['k_min'])
        metal_force.addPerParticleParameter("x0")
        metal_force.addPerParticleParameter("y0")
        metal_force.addPerParticleParameter("z0")

    # Top layers get a THIRD force whose k is driven to zero at NPT (--metal-relax-top).
    metal_top_force, metal_top_cut = None, None
    if metal_resname is not None and metal_cfg.get('relax_top', 0) > 0:
        metal_top_force = CustomExternalForce(formula.replace(" k ", " k_metal_top "))
        metal_top_force.addGlobalParameter("k_metal_top", metal_cfg['k_min'])
        metal_top_force.addPerParticleParameter("x0")
        metal_top_force.addPerParticleParameter("y0")
        metal_top_force.addPerParticleParameter("z0")
        _pnm = modeller.positions.value_in_unit(nanometer)   # positions_nm is defined below
        _pz = sorted({round(_pnm[a.index][2], 4)
                      for a in modeller.topology.atoms()
                      if a.residue.name.strip() == metal_resname})
        _layers = []
        for _z in _pz:
            if not _layers or _z - _layers[-1][-1] > 0.05:   # 0.5 A tolerance
                _layers.append([_z])
            else:
                _layers[-1].append(_z)
        _n = min(metal_cfg['relax_top'], len(_layers))
        metal_top_cut = min(_layers[-_n]) - 1e-6
        print(f"🪷 Slab has {len(_layers)} layers; top {_n} released once NPT starts "
              f"(frozen through minimization + NVT).")

    protein_backbone = ['CA', 'C', 'N']
    nucleic_backbone = ['P', "O3'", "O5'", "C3'", "C4'", "C5'"]
    prot_count = 0
    na_count = 0

    # Get positions from modeller since simulation context doesn't exist yet
    positions_nm = modeller.positions.value_in_unit(nanometer)

    metal_count = 0
    metal_top_count = 0
    for atom in modeller.topology.atoms():
        res_name = atom.residue.name.strip()

        if metal_resname is not None and res_name == metal_resname:
            pos = positions_nm[atom.index]
            if metal_top_cut is not None and pos[2] >= metal_top_cut:
                metal_top_force.addParticle(atom.index, [pos[0], pos[1], pos[2]])
                metal_top_count += 1
            else:
                metal_force.addParticle(atom.index, [pos[0], pos[1], pos[2]])
                metal_count += 1
            continue

        is_prot = atom.name in protein_backbone
        is_na = atom.name in nucleic_backbone

        if is_prot or is_na:
            pos = positions_nm[atom.index]
            restraint_force.addParticle(atom.index, [pos[0], pos[1], pos[2]])

            if is_prot: prot_count += 1
            if is_na: na_count += 1

    system.addForce(restraint_force)
    if metal_force is not None:
        system.addForce(metal_force)
    if metal_top_force is not None:
        system.addForce(metal_top_force)
    # Print the verification report
    print(f"🔒 Restraint Report:")
    print(f"   - Protein backbone: {prot_count}")
    print(f"   - DNA/RNA backbone: {na_count}")
    if metal_resname is not None:
        print(f"   - Metal slab ({metal_resname}): {metal_count}   "
              f"[k_metal held at {metal_cfg['k_prod']:.0f} kJ/mol/nm² through production]")
        if metal_top_force is not None:
            print(f"   - Metal top layers  : {metal_top_count}   "
                  f"[released at NPT -> surface breathes]")
        if metal_count == 0 and metal_top_count == 0:
            print(f"⚠️  WARNING: no '{metal_resname}' atoms were restrained!")
    # note: changed collisions from 1 to 2/ps
    integrator = LangevinMiddleIntegrator(temperature*kelvin, 2/picosecond, timestep*picoseconds)
    simulation = Simulation(modeller.topology, system, integrator)
    simulation.context.setPositions(modeller.positions)

    # Detect Platform
    platform = simulation.context.getPlatform()
    print(f"✨ Running on Platform: {platform.getName()}")
    if platform.getName() in ['CUDA', 'OpenCL', 'Metal']:
        device_name = platform.getPropertyValue(simulation.context, 'DeviceName')
        print(f"💠 Active GPU: {device_name}")
    else:
        print("💻 Running on CPU (no GPU detected)")


    # --- STEP 3: THREE-STEP MINIMIZATION ---
    print("✨ Step 1/3: Minimizing solvent (Restraints HEAVY)... ", end="", flush=True)
    simulation.context.setParameter("k", k_max) 
    if metal_force is not None:
        simulation.context.setParameter("k_metal", metal_cfg['k_min'])
    if metal_top_force is not None:
        simulation.context.setParameter("k_metal_top", metal_cfg['k_min'])
    simulation.minimizeEnergy(tolerance=20.0*kilojoule_per_mole/nanometer)
    print("DONE!")

    print("✨ Step 2/3: Minimizing entire system (Restraints LIGHT)... ", end="", flush=True)
    # Don't go to 0.0 yet! 1.0 keeps the fold while letting sidechains pack.
    simulation.context.setParameter("k", 10.0)
    simulation.minimizeEnergy(tolerance=10.0*kilojoule_per_mole/nanometer)
    print("DONE!")

    print("✨ Step 3/3: Minimizing entire system (NO Restraints)... ", end="", flush=True)
    # Don't go to 0.0 yet! 1.0 keeps the fold while letting sidechains pack.
    simulation.context.setParameter("k", 0.0)
    simulation.minimizeEnergy(tolerance=2.0*kilojoule_per_mole/nanometer)
    print("DONE!")

    # --- DUMP MINIMIZED SNAPSHOT ---
    min_state = simulation.context.getState(getPositions=True, enforcePeriodicBox=False)
    current_positions = min_state.getPositions()
    #print(f"DEBUG: Topology atoms: {simulation.topology.getNumAtoms()} | State atoms: {len(current_positions)}")
    with open('minimized.cif', 'w') as f:
        PDBxFile.writeFile(simulation.topology, current_positions, f)

    # --- STEP 4: EQUILIBRATION (Chil Start) ---
    print(f"--- Step 4: System Equilibration 🌀 ---")

    simulation.context.setVelocitiesToTemperature(temperature*kelvin)

    # 4.1 NVT heating
    nvt_time = 100
    nvt_steps = int(nvt_time / timestep) # 100 ps

    # Disable Barostat
    barostat.setFrequency(0)
    simulation.context.reinitialize(preserveState=True)

     # Re-reference restraints to the minimized structure ← NEW
    refresh_restraint_references(simulation, [f for f in (restraint_force, metal_force, metal_top_force) if f is not None])   

    # Set High Restraints on the solute backbone
    simulation.context.setParameter("k", k_max)
    if metal_force is not None:
        simulation.context.setParameter("k_metal", metal_cfg['k_heat'])
    if metal_top_force is not None:
        simulation.context.setParameter("k_metal_top", metal_cfg['k_heat'])   # still frozen in NVT

    # Run NVT equilibration to stabilize temperature first
    try:
        from mdtraj.reporters import NetCDFReporter
        simulation.reporters.append(NetCDFReporter('heating.nc', report_interval))
    except:
        simulation.reporters.append(DCDReporter('heating.dcd', report_interval))
    print(f"🌡️  Step 1/2: NVT equilibration ({nvt_time} ps):")
    print(f" 🔸 Heating system: T = {temperature:.1f}K ... ", end="", flush=True)
    simulation.step(nvt_steps)
    print("DONE!")

    # 4.2 Step-by-step NPT equilibration
    print(f"🧭 Step 2/2: NPT equilibration ({num_stages} × {time_per_stage:.0f} ps):")
    # Clear NVT reporter
    simulation.reporters.clear()
    # Turn on Barostat (default - every 25 steps)
    barostat.setFrequency(25)
    if metal_force is not None:
        simulation.context.setParameter("k_metal", metal_cfg['k_prod'])
    if metal_top_force is not None:
        simulation.context.setParameter("k_metal_top", 0.0)               # RELEASED at NPT
        print(f"🪷 Top {metal_cfg['relax_top']} metal layer(s) released - the surface now breathes.")
    simulation.context.reinitialize(preserveState=True)

    try:
        from mdtraj.reporters import NetCDFReporter
        simulation.reporters.append(NetCDFReporter('equilibration.nc', report_interval))
    except:
        simulation.reporters.append(DCDReporter('equilibration.dcd', report_interval))

    # 3. GRADUATED EQUILIBRATION LOOP
    k_values = [k_max - (i * (k_max / (num_stages - 1))) for i in range(num_stages)]
    steps_per_stage = int((equil_time / timestep) / len(k_values))

    for k_val in k_values:
        print(f" 🔹 Relaxing system: k = {k_val:.1f}... ", end="", flush=True)
        simulation.context.setParameter("k", k_val)
        simulation.step(steps_per_stage)
        print("DONE!")

    print("👏 Equilibration finished!")

    # --- DUMP EQUILIBRATED SNAPSHOT ---
    eq_state = simulation.context.getState(getPositions=True, enforcePeriodicBox=False)
    eq_positions = eq_state.getPositions()
    with open('equilibrated.cif', 'w') as f:
        PDBxFile.writeFile(simulation.topology, eq_positions, f)

    # Turn off restraints for Production by setting force constant k to 0
    # (the metal slab deliberately KEEPS k_metal - see AMBER tutorial 27)
    simulation.context.setParameter("k", 0.0)
    if metal_top_force is not None:
        simulation.context.setParameter("k_metal_top", 0.0)               # stays released
    if metal_force is not None:
        simulation.context.setParameter("k_metal", metal_cfg['k_prod'])
        print(f"🔒 Metal slab stays restrained during production "
              f"(k_metal = {metal_cfg['k_prod']:.0f} kJ/mol/nm²).")
    
    # Clear the equilibration reporter to start a fresh production file
    #simulation.reporters.pop()
    simulation.reporters.clear()

    # --- STEP 5: PRODUCTION ---
    # Add Production Reporters
    try:
        from mdtraj.reporters import NetCDFReporter
        simulation.reporters.append(NetCDFReporter(output_nc, report_interval))
    except:
        simulation.reporters.append(DCDReporter('production.dcd', report_interval))

    print(f"--- Step 5: Production Run 👑 ---")
    steps_per_percent = int(total_steps / 100)

    print(f"\n🔮 Simulation Duration: {total_steps * timestep / 1000:.2f} ns")
    print(f"\n{'%':>4} {'Step':>10} {'Pot. Energy':>15} {'Temp(K)':>8} {'Volume(nm³)':>8} {'Speed(ns/day)':>10}")
    print("-" * 55)

    production_start = time.time()
    last_time = time.time()

    for i in range(1, 101):
            simulation.step(steps_per_percent)
            
            if i % log_freq == 0:
                # FIX:
                state = simulation.context.getState(getEnergy=True)
                pot_energy = state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
                kin_energy = state.getKineticEnergy()

                # GET VOLUME (To verify NPT)
                box_vectors = state.getPeriodicBoxVectors()
                volume = state.getPeriodicBoxVolume().value_in_unit(nanometer**3)
                                
                # get the Temperature:
                # FIX for OPC water: Ccunt only atoms that have mass > 0
                real_atoms_dof = 0
                for atom in modeller.topology.atoms():
                    # Get mass from the system using the atom index
                    if system.getParticleMass(atom.index).value_in_unit(dalton) > 0:
                        real_atoms_dof += 3
                dof = real_atoms_dof - system.getNumConstraints()
                # Using OpenMM units for Boltzmann and Avogadro constants
                temp = (2 * kin_energy / (dof * BOLTZMANN_CONSTANT_kB * AVOGADRO_CONSTANT_NA)).value_in_unit(kelvin)
                
                current_time = time.time()
                elapsed = current_time - last_time
                speed = (steps_per_percent * log_freq * timestep * 0.001) / (elapsed / 86400)
                
                print(f"{i:>3}% {simulation.currentStep:>10} {pot_energy:>15.1f} {temp:>8.1f} {volume:>10.1f} {speed:>10.1f}")
                last_time = current_time

    # Final Summary
    production_end = time.time()
    total_seconds = production_end - production_start
    avg_speed = (total_steps * timestep * 0.001) / (total_seconds / 86400)

    print("-" * 55)
    print(f"WORK COMPLETED!")
    print(f"⏳Total wall-clock time: {int(total_seconds // 60)}m {int(total_seconds % 60)}s")
    print(f"⚜️ Average performance: {avg_speed:.2f} ns/day")
    print(f"🌀Trajectory saved to: {output_nc}")

# --- CLI ARGUMENTS (defaults reproduce the original hard-coded MAIN CONFIG) ---
def parse_args():
    parser = argparse.ArgumentParser(
        description="👑 RoyalMD GOLD Edition: one code, one innovation.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("input_pdb", help="Input structure (PDB) to simulate.")

    # --- files ---
    parser.add_argument("--pre-sim-pdb", default="solvated.pdb",
                        help="Solvated system PDB written before the run.")
    parser.add_argument("--output-nc", default="production.nc",
                        help="Production trajectory file name.")

    # --- integrator / timing ---
    parser.add_argument("--timestep", type=float, default=0.002,
                        help="Integration timestep (ps).")
    parser.add_argument("--equil-time", type=int, default=1000,
                        help="NPT equilibration time (ps).")
    parser.add_argument("--production-time", type=int, default=50000,
                        help="Production time (ps).")
    parser.add_argument("--report-interval", type=int, default=5000,
                        help="Reporter interval (steps).")
    parser.add_argument("--log-freq", type=int, default=10,
                        help="Console log frequency (print every N%% of production).")

    # --- thermodynamic conditions ---
    parser.add_argument("--temperature", type=float, default=310,
                        help="Temperature (K).")
    parser.add_argument("--pressure", type=float, default=1.0,
                        help="Pressure (bar).")

    # --- box & solvation ---
    parser.add_argument("--box-type", default="dodecahedron",
                        choices=["dodecahedron", "octahedron", "cube"],
                        help="Periodic box shape.")
    parser.add_argument("--box-padding", type=float, default=1.2,
                        help="Solvent padding (nm).")
    parser.add_argument("--ionic-strength", type=float, default=0.15,
                        help="Ionic strength (mol/L).")
    parser.add_argument("--pH", dest="env_pH", type=float, default=7.0,
                        help="pH used for protonation.")

    # --- force fields ---
    parser.add_argument("--ff-protein", default="amber14", choices=list(FF_MAP.keys()),
                        help="Protein force field.")
    parser.add_argument("--ff-water", default="tip3p",
                        help="Water model (e.g. tip3p, opc, spce, water).")

    # --- cleanup ---
    parser.add_argument("--strip-residues", nargs="+",
                        default=["SO4", "HOH", "EDO", "GOL", "LIG", "LIH"],
                        help="Residue names removed before the system is built "
                             "(crystallographic additives and any small molecules). "
                             "This build simulates the macromolecule only.")


    # --- GOLD EDITION: metal surface immobilization (opt-in, nothing else changes) ---
    g = parser.add_argument_group("metal surface (--metalsurf)")
    g.add_argument("--metalsurf", action="store_true",
                   help="Immobilize the solute on an fcc metal surface (INTERFACE FF).")
    g.add_argument("--interface-root", default="INTERFACE_FF_1_5",
                   help="Path to the extracted INTERFACE force field package.")
    g.add_argument("--metal", default="au", choices=sorted(METAL_LJ.keys()),
                   help="fcc metal (INTERFACE neutral-metal set).")
    g.add_argument("--metal-facet", default="111", choices=["100", "110", "111"],
                   help="Exposed crystal facet.")
    g.add_argument("--metal-cells-z", type=int, default=2,
                   help="Slab thickness in UNIT CELLS along z "
                        "(fcc{111}: 3 atomic layers per cell, so 2 = 6 layers).")
    g.add_argument("--metal-separation", type=float, default=3.0,
                   help="Initial gap (A) between the solute's lowest atom and the top metal plane.")
    g.add_argument("--metal-headroom", type=float, default=25.0,
                   help="Water headroom (A) above the solute.")
    g.add_argument("--metal-under", type=float, default=2.0,
                   help="Space (A) below the slab before the box floor.")
    g.add_argument("--metal-orient", action=argparse.BooleanOptionalAction, default=True,
                   help="Align principal axes to the box and point the active-site "
                        "pocket away from the surface.")
    g.add_argument("--metal-relax-top", type=int, default=0,
                   help="Release the top N metal layers once NPT begins (0 = keep the whole "
                        "slab restrained, the default). They stay frozen through minimization "
                        "and NVT, then move freely in NPT and production, so the surface can "
                        "breathe against the solute while the lower layers anchor the lattice.")
    g.add_argument("--metal-k-min", type=float, default=83680.0,
                   help="Metal restraint k during minimization (kJ/mol/nm²; "
                        "tutorial 100 kcal/mol/A²).")
    g.add_argument("--metal-k-heat", type=float, default=41840.0,
                   help="Metal restraint k during heating (tutorial 50 kcal/mol/A²).")
    g.add_argument("--metal-k-prod", type=float, default=16736.0,
                   help="Metal restraint k for NPT + production (tutorial 20 kcal/mol/A²).")

    return parser.parse_args()

# --- MAIN CONTROLLER ---
def main():
    args = parse_args()

    # --- added 27/01: Individual Working Directory ---
    input_pdb_path = os.path.abspath(args.input_pdb) # Absolute path so we don't lose the file
    if not os.path.isfile(input_pdb_path):
        print(f"ERROR: Input PDB file not found: {args.input_pdb}")
        sys.exit(1)
    input_filename = os.path.basename(input_pdb_path)
    # Create folder name (e.g., 'protein.pdb' -> 'MD_protein_')
    date_str = datetime.now().strftime("%d%m%y")
    work_dir = f"MD_{os.path.splitext(input_filename)[0]}_{date_str}"

    # --- MAIN CONFIG (now sourced from the CLI; defaults reproduce the originals) ---
    pre_sim_pdb = args.pre_sim_pdb
    output_nc = args.output_nc
    timestep = args.timestep
    equil_time = args.equil_time              # NPT equilibration - 1ns
    production_time = args.production_time     # 50 ns of production
    report_interval = args.report_interval     # report info every 5000 steps
    log_freq = args.log_freq
    # conditions
    temperature = args.temperature
    pressure = args.pressure
    # box and solvation: dodecahedron or octahedron
    box_type = args.box_type
    box_padding = args.box_padding             # default is 1.2
    ionic_strength = args.ionic_strength
    env_pH = args.env_pH
    # main params
    ff_protein = args.ff_protein               # best practice: amber14 with tip3p
    ff_water = args.ff_water                    # tip3p - amber14/amber99; opc - amber19; water - charmm36
    # Residues removed before building (additives + any small molecules)
    strip_residues = list(args.strip_residues)

    # --- GOLD EDITION config (None keeps v1.4 behaviour exactly) ---
    metal_cfg = None
    if args.metalsurf:
        # CHARMM36's 1-4 scaling differs from the generated metal XML, so OpenMM
        # refuses to merge them. Fail early with a clear message, not a traceback.
        if ff_protein == 'charmm36':
            print("❌ --metalsurf is incompatible with --ff-protein charmm36:\n"
                  "   OpenMM refuses to merge the metal XML ('multiple NonbondedForce "
                  "tags with different 1-4 scales').\n"
                  "   Use an AMBER force field, e.g. --ff-protein amber14 --ff-water tip3p.")
            sys.exit(1)
        metal_cfg = {
            'metal':           args.metal,
            'facet':           args.metal_facet,
            'interface_root':  os.path.abspath(args.interface_root),
            'cells_z':         args.metal_cells_z,
            'separation':      args.metal_separation,
            'headroom':        args.metal_headroom,
            'under':           args.metal_under,
            'orient':          args.metal_orient,
            'padding':         box_padding,      # XY slab size follows --box-padding
            'k_min':           args.metal_k_min,
            'k_heat':          args.metal_k_heat,
            'k_prod':          args.metal_k_prod,
            'relax_top':       args.metal_relax_top,
        }
        if box_type != 'cube':
            print(f"⚠️  --metalsurf requires a rectangular cell; ignoring "
                  f"--box-type {box_type} (the metal lattice cannot tile a {box_type}).")

    ff_map = FF_MAP

    # --- 👑 EXECUTION PIPELINE 👑 ---
    print_intro_banner()
    # 0.1 (early versions):
    #clean_old_files(pre_sim_pdb, activate=True)
    # 0.1B (new version):
    clean_old_dir(work_dir, activate=True)
    # move into workdir
    os.chdir(work_dir)
    input_pdb = input_pdb_path

    # 1. Fix Protein (stripping additives and small molecules)
    fixer = fix_protein_topology(input_pdb, env_pH, strip_residues)

    # 2: Solvate
    modeller, forcefield, metal_resname, res_templates = solvate_system(
        fixer, ff_protein, ff_water, ff_map,
        box_type, box_padding, ionic_strength, pre_sim_pdb, input_pdb,
        metal_cfg=metal_cfg, work_dir='.')

    # 3 & 4: Run calculations
    setup_and_run(modeller, 
                  forcefield, 
                  temperature, 
                  pressure, 
                  timestep, 
                  equil_time,
                  production_time, 
                  report_interval, 
                  output_nc, 
                  log_freq,
                  metal_resname=metal_resname,
                  res_templates=res_templates,
                  metal_cfg=metal_cfg)

if __name__ == "__main__":
    main()