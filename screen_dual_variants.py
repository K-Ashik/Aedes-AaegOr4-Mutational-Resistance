import os
import numpy as np
import pandas as pd
from Bio import PDB
import pubchempy as pcp
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

WILD_PDBQT = os.path.join(PDBQT_DIR, "AaegOr4_Wild.pdbqt")
MUTANT_PDBQT = os.path.join(PDBQT_DIR, "AaegOr4_Mutant_Y150A.pdbqt")
WILD_PDB = os.path.join(PDB_DIR, "AaegOr4_A0A1S7UE45.pdb")

# Natural Volatiles to Screen
NATURAL_COMPOUNDS = [
    "DEET",             # Control Benchmark
    "Citronellal",      # Citronella oil active component
    "Thymol",           # Thyme oil component
    "Carvacrol",        # Oregano/Thyme derivative
    "Eucalyptol",       # Eucalyptus oil (1,8-cineole)
    "Geraniol",         # Rose / Citronella oil constituent
    "Menthol",          # Mint oil derivative
    "Limonene",         # Citrus peel volatile
    "Eugenol",          # Clove oil constituent
    "Linalool",         # Lavender oil volatile
    "Alpha-Pinene",     # Pine resin terpene
    "Camphor",          # Camphor laurel terpene
    "Cinnamaldehyde",   # Cinnamon bark oil
    "Terpineol",        # Tea tree oil constituent
    "Pulegone"          # Pennyroyal oil terpene
]

# ==========================================
# 1. FETCH & CONVERT LIGANDS
# ==========================================
def name_to_pdbqt(name: str) -> str:
    """Retrieves 3D chemical structure from PubChem and exports PDBQT."""
    compounds = pcp.get_compounds(name, "name")
    if not compounds:
        return None
        
    smiles = getattr(compounds[0], "smiles", getattr(compounds[0], "isomeric_smiles", None))
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

def get_receptor_center(pdb_path: str):
    """Calculates geometric center of receptor protein."""
    parser = PDB.PDBParser(QUIET=True)
    struct = parser.get_structure("receptor", pdb_path)
    coords = [atom.get_coord() for atom in struct.get_atoms() if atom.element != 'H']
    return np.array(coords).mean(axis=0).tolist()

# ==========================================
# 2. DUAL VARIANT VIRTUAL SCREENING
# ==========================================
def run_dual_screening():
    print("="*60)
    print("STEP 2: SCREENING NATURAL VOLATILES AGAINST WILD & MUTANT AaegOr4")
    print("="*60)
    
    center = get_receptor_center(WILD_PDB)
    results = []
    
    receptors = [
        ("Wild_Type", WILD_PDBQT),
        ("Mutant_Y150A", MUTANT_PDBQT)
    ]
    
    for idx, name in enumerate(NATURAL_COMPOUNDS, 1):
        print(f"[{idx}/{len(NATURAL_COMPOUNDS)}] Screening {name}...")
        try:
            lig_pdbqt = name_to_pdbqt(name)
            if not lig_pdbqt:
                continue
                
            scores = {}
            for rec_label, rec_pdbqt in receptors:
                v = Vina(sf_name='vina')
                v.set_receptor(rec_pdbqt)
                v.set_ligand_from_file(lig_pdbqt)
                v.compute_vina_maps(center=center, box_size=[30.0, 30.0, 30.0])
                
                v.dock(exhaustiveness=8, n_poses=1)
                scores[rec_label] = round(v.energies(n_poses=1)[0][0], 2)
                
            w_score = scores["Wild_Type"]
            m_score = scores["Mutant_Y150A"]
            
            # Classification: Pan-Inhibitor if it maintains strong binding (<= -5.0) on BOTH
            if w_score <= -5.0 and m_score <= -5.0:
                classification = "Pan-Inhibitor (Resistance Proof!)"
            elif w_score <= -5.0:
                classification = "Wild-Type Selective"
            else:
                classification = "Moderate/Weak"
                
            results.append({
                "Compound": name,
                "Wild_Type_kcal_mol": w_score,
                "Mutant_Y150A_kcal_mol": m_score,
                "Delta_Delta_G": round(m_score - w_score, 2),
                "Classification": classification
            })
            print(f"[+] {name} -> Wild: {w_score} | Mutant: {m_score} | {classification}")
            
        except Exception as e:
            print(f"[-] Error with {name}: {e}")
            
    df_res = pd.DataFrame(results).sort_values(by="Mutant_Y150A_kcal_mol").reset_index(drop=True)
    
    print("\n" + "="*60)
    print("DUAL-VARIANT VIRTUAL SCREENING LEADERBOARD")
    print("="*60)
    print(df_res.to_string(index=False))
    
    output_csv = os.path.join(RESULTS_DIR, "dual_variant_screening_leaderboard.csv")
    df_res.to_csv(output_csv, index=False)
    print(f"\n[+] Saved dual screening leaderboard to {output_csv}")

if __name__ == "__main__":
    run_dual_screening()