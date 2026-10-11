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

WILD_PDB = os.path.join(PDB_DIR, "AaegOr4_A0A1S7UE45.pdb")
MUTANT_PDB = os.path.join(PDB_DIR, "AaegOr4_Mutant_Y150A.pdb")

WILD_PDBQT = os.path.join(PDBQT_DIR, "AaegOr4_Wild.pdbqt")
MUTANT_PDBQT = os.path.join(PDBQT_DIR, "AaegOr4_Mutant_Y150A.pdbqt")

# Control ligands to measure resistance impact
TEST_LIGANDS = {
    "DEET": "CCN(CC)C(=O)C1=CC=CC(=C1)C",     # Synthetic Repellent
    "Sulcatone": "CC(=CCCC(=O)C)C"              # Primary Human Skin Volatile
}

# ==========================================
# 1. IN SILICO MUTAGENESIS ENGINE (Y150A)
# ==========================================
def mutate_tyrosine_to_alanine(input_pdb: str, output_pdb: str, target_resnum: int = 150):
    """
    Parses PDB and replaces Tyrosine (TYR) at target residue position with Alanine (ALA)
    by retaining backbone atoms (N, CA, C, O, CB) and removing side-chain heavy atoms.
    """
    print(f"[*] Performing in silico mutagenesis on {input_pdb}...")
    print(f"[*] Targeting Residue position: {target_resnum} (TYR -> ALA)...")
    
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("receptor", input_pdb)
    
    mutated_count = 0
    
    for model in structure:
        for chain in model:
            for residue in chain:
                # Check residue sequence number
                if residue.id[1] == target_resnum and residue.resname in ["TYR", "TRP", "PHE"]:
                    old_res = residue.resname
                    residue.resname = "ALA"  # Rename residue to Alanine
                    
                    # Remove side chain atoms beyond Beta-Carbon (CB)
                    atoms_to_remove = [atom.id for atom in residue if atom.id not in ["N", "CA", "C", "O", "CB"]]
                    for atom_id in atoms_to_remove:
                        residue.detach_child(atom_id)
                        
                    mutated_count += 1
                    print(f"[+] Successfully mutated Residue {old_res}{target_resnum} -> ALA{target_resnum}")

    io = PDB.PDBIO()
    io.set_structure(structure)
    io.save(output_pdb)
    print(f"[+] Saved Mutant Protein PDB to {output_pdb}")
    return output_pdb

# ==========================================
# 2. STRICT FIXED-WIDTH PDBQT WRITER
# ==========================================
def clean_receptor_to_pdbqt(input_pdb: str, output_pdbqt: str):
    """Converts PDB to strict fixed-column PDBQT format required by Vina."""
    with open(input_pdb, "r") as infile, open(output_pdbqt, "w") as outfile:
        for line in infile:
            if line.startswith(("ATOM  ", "HETATM")):
                atom_num = line[6:11]
                atom_name = line[12:16]
                res_name = line[17:20]
                chain_id = line[21] if len(line) > 21 else "A"
                res_num = line[22:26]
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                
                elem = line[76:78].strip() if len(line) >= 78 else atom_name.strip()[0]
                if not elem or elem.isdigit():
                    elem = atom_name.strip()[0]
                
                pdbqt_line = f"ATOM  {atom_num:>5} {atom_name:<4} {res_name:>3} {chain_id}{res_num:>4}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 00.00    +0.000 {elem:<2}\n"
                outfile.write(pdbqt_line)
            elif line.startswith("END"):
                outfile.write("END\n")

# ==========================================
# 3. GENERATE LIGAND PDBQT
# ==========================================
def ligand_smiles_to_pdbqt(name: str, smiles: str) -> str:
    """Generates 3D structure from SMILES using RDKit and exports PDBQT with Meeko."""
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    AllChem.MMFFOptimizeMolecule(mol)
    
    preparator = MoleculePreparation()
    setup_objects = preparator.prepare(mol)
    mol_setup = setup_objects[0] if isinstance(setup_objects, list) else setup_objects
    
    pdbqt_string = PDBQTWriterLegacy.write_string(mol_setup)[0]
    out_path = os.path.join(PDBQT_DIR, f"{name}.pdbqt")
    
    with open(out_path, "w") as f:
        f.write(pdbqt_string)
        
    return out_path

# ==========================================
# 4. CALCULATE CENTER OF MASS
# ==========================================
def get_receptor_center(pdb_path: str):
    """Calculates geometric center of receptor protein."""
    parser = PDB.PDBParser(QUIET=True)
    struct = parser.get_structure("receptor", pdb_path)
    coords = [atom.get_coord() for atom in struct.get_atoms() if atom.element != 'H']
    return np.array(coords).mean(axis=0).tolist()

# ==========================================
# 5. DOCKING COMPARISON PIPELINE
# ==========================================
def run_mutational_analysis():
    print("="*60)
    print("STEP 1: IN SILICO MUTAGENESIS & RESISTANCE MODELING")
    print("="*60)
    
    # 1. Generate Mutant Structure
    mutate_tyrosine_to_alanine(WILD_PDB, MUTANT_PDB, target_resnum=150)
    
    # 2. Convert Receptors to PDBQT
    clean_receptor_to_pdbqt(WILD_PDB, WILD_PDBQT)
    clean_receptor_to_pdbqt(MUTANT_PDB, MUTANT_PDBQT)
    
    center = get_receptor_center(WILD_PDB)
    results = []
    
    receptors = [
        ("Wild-Type (AaegOr4)", WILD_PDBQT),
        ("Mutant (Y150A)", MUTANT_PDBQT)
    ]
    
    for lig_name, smiles in TEST_LIGANDS.items():
        lig_pdbqt = ligand_smiles_to_pdbqt(lig_name, smiles)
        energies_dict = {}
        
        for rec_label, rec_pdbqt in receptors:
            out_docked = os.path.join(RESULTS_DIR, f"{lig_name}_{rec_label.split()[0]}_docked.pdbqt")
            
            v = Vina(sf_name='vina')
            v.set_receptor(rec_pdbqt)
            v.set_ligand_from_file(lig_pdbqt)
            v.compute_vina_maps(center=center, box_size=[30.0, 30.0, 30.0])
            
            v.dock(exhaustiveness=8, n_poses=1)
            v.write_poses(out_docked, n_poses=1, overwrite=True)
            
            energy = round(v.energies(n_poses=1)[0][0], 2)
            energies_dict[rec_label] = energy
            
        wild_e = energies_dict["Wild-Type (AaegOr4)"]
        mut_e = energies_dict["Mutant (Y150A)"]
        delta_delta_g = round(mut_e - wild_e, 2)  # Positive value means loss of binding affinity
        
        results.append({
            "Ligand": lig_name,
            "Wild_Type_kcal_mol": wild_e,
            "Mutant_Y150A_kcal_mol": mut_e,
            "Delta_Delta_G_kcal_mol": delta_delta_g,
            "Impact": f"Binding Lost (+{delta_delta_g} kcal/mol)" if delta_delta_g > 0 else "No Resistance Effect"
        })

    df_res = pd.DataFrame(results)
    
    print("\n" + "="*60)
    print("GENETIC RESISTANCE ANALYSIS TABLE (Y150A MUTATION)")
    print("="*60)
    print(df_res.to_string(index=False))
    
    # Save CSV Report
    report_path = os.path.join(RESULTS_DIR, "genetic_resistance_Y150A.csv")
    df_res.to_csv(report_path, index=False)
    print(f"\n[+] Saved mutational resistance report to {report_path}")

if __name__ == "__main__":
    run_mutational_analysis()