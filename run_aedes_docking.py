import os
import numpy as np
import pandas as pd
from Bio import PDB
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy
from vina import Vina

# ==========================================
# DIRECTORY SETUP
# ==========================================
DATA_DIR = "./aedes_data"
PDB_DIR = os.path.join(DATA_DIR, "structures")
PDBQT_DIR = os.path.join(DATA_DIR, "pdbqt")
RESULTS_DIR = os.path.join(DATA_DIR, "docking_results")

os.makedirs(PDBQT_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)

RECEPTOR_PDB = os.path.join(PDB_DIR, "AaegOr4_A0A1S7UE45.pdb")
RECEPTOR_PDBQT = os.path.join(PDBQT_DIR, "AaegOr4.pdbqt")

# Volatile ligands to dock into AaegOr4
LIGANDS_SMILES = {
    "Sulcatone": "CC(=CCCC(=O)C)C",
    "1-Octen-3-ol": "CCCCCC(C=C)O",
    "Lactic_acid": "CC(C(=O)O)O",
    "DEET": "CCN(CC)C(=O)C1=CC=CC(=C1)C"
}

# ==========================================
# 1. STRICT FIXED-WIDTH PDBQT RECEPTOR WRITER
# ==========================================
def clean_receptor_to_pdbqt(input_pdb: str, output_pdbqt: str):
    """Converts PDB to strict fixed-column PDBQT format required by Vina."""
    print(f"[*] Preparing receptor PDBQT for {input_pdb}...")
    with open(input_pdb, "r") as infile, open(output_pdbqt, "w") as outfile:
        for line in infile:
            if line.startswith(("ATOM  ", "HETATM")):
                # Extract strict fields according to PDBQT specification
                atom_num = line[6:11]
                atom_name = line[12:16]
                res_name = line[17:20]
                chain_id = line[21] if len(line) > 21 else "A"
                res_num = line[22:26]
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                
                # Determine element symbol
                elem = line[76:78].strip() if len(line) >= 78 else atom_name.strip()[0]
                if not elem or elem.isdigit():
                    elem = atom_name.strip()[0]
                
                # Format strict fixed-width PDBQT line
                # Columns: ATOM, serial, name, resName, chain, resSeq, X, Y, Z, occ, temp, charge, type
                pdbqt_line = f"ATOM  {atom_num:>5} {atom_name:<4} {res_name:>3} {chain_id}{res_num:>4}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 00.00    +0.000 {elem:<2}\n"
                outfile.write(pdbqt_line)
            elif line.startswith("END"):
                outfile.write("END\n")
    print(f"[+] Saved Receptor PDBQT to {output_pdbqt}")

# ==========================================
# 2. GENERATE LIGAND 3D CONFORMERS & PDBQT
# ==========================================
def ligand_smiles_to_pdbqt(name: str, smiles: str) -> str:
    """Generates 3D structure from SMILES using RDKit and exports PDBQT with Meeko."""
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    
    # Generate 3D coordinates via ETKDG algorithm
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    AllChem.MMFFOptimizeMolecule(mol)
    
    # Prepare ligand setup using Meeko
    preparator = MoleculePreparation()
    setup_objects = preparator.prepare(mol)
    mol_setup = setup_objects[0] if isinstance(setup_objects, list) else setup_objects
    
    # Write PDBQT string
    pdbqt_string = PDBQTWriterLegacy.write_string(mol_setup)[0]
    
    out_path = os.path.join(PDBQT_DIR, f"{name}.pdbqt")
    with open(out_path, "w") as f:
        f.write(pdbqt_string)
        
    print(f"[+] Prepared 3D Ligand PDBQT: {name}")
    return out_path

# ==========================================
# 3. CALCULATE BINDING POCKET CENTER OF MASS
# ==========================================
def get_receptor_center_and_size(pdb_path: str):
    """Calculates geometric center and box dimensions of receptor protein."""
    parser = PDB.PDBParser(QUIET=True)
    struct = parser.get_structure("receptor", pdb_path)
    coords = [atom.get_coord() for atom in struct.get_atoms() if atom.element != 'H']
    
    coords = np.array(coords)
    center = coords.mean(axis=0)
    box_size = coords.max(axis=0) - coords.min(axis=0) + 10.0
    
    return center.tolist(), box_size.tolist()

# ==========================================
# 4. RUN AUTODOCK VINA VIRTUAL SCREENING
# ==========================================
def run_docking_pipeline():
    clean_receptor_to_pdbqt(RECEPTOR_PDB, RECEPTOR_PDBQT)
    center, box_size = get_receptor_center_and_size(RECEPTOR_PDB)
    
    print(f"\n[*] Receptor Grid Center (X, Y, Z): {[round(c, 2) for c in center]}")
    print(f"[*] Search Box Size (X, Y, Z): {[round(s, 2) for s in box_size]}")
    
    docking_results = []
    
    for name, smiles in LIGANDS_SMILES.items():
        ligand_pdbqt = ligand_smiles_to_pdbqt(name, smiles)
        out_docked_pdbqt = os.path.join(RESULTS_DIR, f"{name}_docked.pdbqt")
        
        # Initialize AutoDock Vina Engine
        v = Vina(sf_name='vina')
        v.set_receptor(RECEPTOR_PDBQT)
        v.set_ligand_from_file(ligand_pdbqt)
        v.compute_vina_maps(center=center, box_size=[30.0, 30.0, 30.0])
        
        # Perform Virtual Screening (Exhaustiveness = 8)
        print(f"[*] Docking {name} into AaegOr4 receptor pocket...")
        v.dock(exhaustiveness=8, n_poses=3)
        v.write_poses(out_docked_pdbqt, n_poses=1, overwrite=True)
        
        # Extract Top Binding Energy (kcal/mol)
        energies = v.energies(n_poses=1)
        top_affinity = round(energies[0][0], 2)
        
        docking_results.append({
            "Ligand": name,
            "Binding_Affinity_kcal_mol": top_affinity,
            "Status": "Strong Binding" if top_affinity <= -5.0 else "Moderate/Weak"
        })
        
    df_res = pd.DataFrame(docking_results).sort_values(by="Binding_Affinity_kcal_mol")
    
    print("\n" + "="*60)
    print("AUTODOCK VINA MOLECULAR DOCKING RESULTS (AaegOr4)")
    print("="*60)
    print(df_res.to_string(index=False))
    
    # Save CSV Report
    df_res.to_csv(os.path.join(RESULTS_DIR, "binding_affinities.csv"), index=False)
    print(f"\n[+] Saved docking report to {RESULTS_DIR}/binding_affinities.csv")

if __name__ == "__main__":
    run_docking_pipeline()