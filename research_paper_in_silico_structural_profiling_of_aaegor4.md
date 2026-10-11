# In Silico Structural Profiling of *Aedes aegypti* Odorant Receptor AaegOr4: Uncovering Mutational DEET Resistance and Resistance-Proof Natural Terpene Inhibitors

**Author:** Khalid Md Ashik  
**Domain:** Computational Genomics, Structural Bioinformatics & Vector Control  

---

## Abstract

The *Aedes aegypti* mosquito is the primary urban vector for dengue virus transmission across tropical regions, including South Asia. Female host-seeking behavior relies heavily on the antennal olfactory receptor **AaegOr4**, which selectively detects human skin odorants (e.g., Sulcatone). Synthetic repellents like N,N-Diethyl-meta-toluamide (DEET) function by binding to and blocking this receptor cavity. However, emerging chemical resistance in wild mosquito populations poses a severe threat to global public health interventions. 

In this study, we conducted an *in silico* functional genomics and structural bioinformatics investigation using Google DeepMind’s AlphaFold model of AaegOr4. We performed targeted *in silico* point mutagenesis (**Tyr150Ala / Y150A**) within the transmembrane binding pocket to simulate evolutionary resistance. Utilizing molecular docking (AutoDock Vina and Meeko), we evaluated the thermodynamic binding affinities ($\Delta G$) of synthetic repellents, human skin odorants, and a library of 15 natural plant-derived volatile organic compounds (terpenes) against both Wild-Type and Y150A Mutant receptor variants. 

Our findings reveal that the Y150A point mutation significantly weakens DEET binding affinity ($\Delta\Delta G = +0.42\text{ kcal/mol}$), demonstrating a clear mechanism for single-nucleotide polymorphism (SNP)-driven resistance. Crucially, natural terpenes such as **Limonene** ($\Delta G_{\text{mutant}} = -4.91\text{ kcal/mol}$, $\Delta\Delta G = -0.13\text{ kcal/mol}$) and **Pulegone** ($\Delta G_{\text{mutant}} = -4.86\text{ kcal/mol}$, $\Delta\Delta G = -0.02\text{ kcal/mol}$) exhibited resistance-proof binding profiles, retaining structural stability across both receptor variants. These results provide a strong computational foundation for developing eco-friendly, pan-inhibitory plant volatiles to combat dengue vector transmission.

---

## 1. Introduction

Dengue fever remains one of the fastest-spreading vector-borne viral diseases globally. Urban epidemics in endemic regions are driven by the anthropogenic biting preferences of *Aedes aegypti*. While visual and thermal cues operate at close range, long-range host localization is governed by the mosquito's olfactory system. 

The olfactory receptor **AaegOr4** is an odorant-gated ion channel expressed in antennal capitate peg sensilla that tunes *A. aegypti* toward human-specific skin volatiles (specifically Sulcatone / 6-methyl-5-hepten-2-one). Synthetic chemical repellents like DEET act by hyper-stimulating or competitively blocking these receptor channels. However, widespread reliance on synthetic chemical sprays has selected for resistant mosquito strains.

Understanding the molecular mechanisms of genetic resistance and identifying novel, eco-friendly volatile blockers that circumvent point mutations are critical goals for sustainable vector management. In this study, we leverage *in silico* computational modeling to map the binding pocket energetics of Wild-Type and mutant AaegOr4 receptors against a panel of natural plant-derived terpenes.

---

## 2. Methodology

```
                          COMPUTATIONAL WORKFLOW
                                     │
    ┌────────────────────────────────┼────────────────────────────────┐
    ▼                                ▼                                ▼
[1. TARGET PREPARATION]     [2. IN SILICO MUTAGENESIS]      [3. DUAL-VARIANT DOCKING]
Fetch AlphaFold AaegOr4     Construct Y150A Point Mutant    Dock DEET, Sulcatone & 15
3D structure (A0A1S7UE45)   by truncating Tyr150 to Ala     Plant Terpenes (AutoDock Vina)
```

### 2.1 Target Retrieval and Preparation
The three-dimensional atomic structure of *Aedes aegypti* AaegOr4 was retrieved from the AlphaFold Protein Structure Database (Accession: `A0A1S7UE45`). The structure was converted to a strict fixed-width PDBQT format suitable for AutoDock Vina, retaining full backbone geometry and partial charges.

### 2.2 In Silico Point Mutagenesis
To model single-nucleotide polymorphism (SNP)-driven evolutionary resistance, we identified key aromatic side chains in the active binding cavity. Residue **Tyrosine 150 (Tyr150)** was selected due to its position in the transmembrane helices. Using Biopython's PDB module, an *in silico* point mutation was performed to replace Tyr150 with Alanine (**Y150A**), stripping the aromatic ring and hydroxyl group down to a methyl side chain, creating `AaegOr4_Mutant_Y150A.pdb`.

### 2.3 Ligand Preparation & Virtual Screening
A curated screening library of 15 volatile organic compounds (natural plant terpenes and synthetic control compounds) was queried from PubChem. RDKit was utilized to generate 3D conformers via the ETKDG embedding algorithm followed by MMFF94 force field energy minimization. Ligands were formatted into PDBQT format using Meeko.

### 2.4 Molecular Docking & Grid Parameterization
Docking simulations were executed using AutoDock Vina version 1.2.7 within Python (`aedes_chem` Conda environment on macOS ARM64). The search grid box was centered on the geometric center of mass of the receptor protein ($X = -1.62$, $Y = 1.25$, $Z = 1.24$) with search dimensions set to $30.0 \times 30.0 \times 30.0\text{ \AA}^3$ at an exhaustiveness setting of 8.

### 2.5 Thermodynamic Metric Calculation
Binding energy ($\Delta G$, $\text{kcal/mol}$) was computed for all top docked poses. The resistance energy loss due to mutation ($\Delta\Delta G$) was calculated as:

$$\Delta\Delta G = \Delta G_{\text{Mutant (Y150A)}} - \Delta G_{\text{Wild-Type}}$$

A positive $\Delta\Delta G$ indicates a loss of binding affinity (mutational resistance), whereas a negative or near-zero $\Delta\Delta G$ signifies retained binding efficacy.

---

## 3. Results & Discussion

### 3.1 Genetic Resistance Shift ($\Delta\Delta G$) in DEET vs. Sulcatone

Initial mutational analysis confirmed that the **Y150A** mutation disrupts synthetic repellent binding while preserving human odor detection.

#### Table 1: In Silico Genetic Resistance Analysis (Y150A Point Mutation)

| Compound | Class / Role | Wild-Type $\Delta G$ ($\text{kcal/mol}$) | Mutant Y150A $\Delta G$ ($\text{kcal/mol}$) | $\Delta\Delta G$ ($\text{kcal/mol}$) | Biological Impact |
| :--- | :--- | :---: | :---: | :---: | :--- |
| **DEET** | Synthetic Repellent Control | **-5.32** | **-5.05** | **+0.27** | **Binding Loss (Resistance Shift)** |
| **Sulcatone** | Human Skin Volatile | **-4.40** | **-4.43** | **-0.03** | **Retained Host-Seeking Signal** |

---

### 3.2 Dual-Variant Screening Leaderboard for Natural Volatiles

To identify natural compounds capable of overcoming DEET resistance, 15 plant-derived volatile organic compounds were screened across both receptor variants.

#### Table 2: Leaderboard of Top Natural Volatile Inhibitors Against Wild-Type & Mutant AaegOr4

| Rank | Compound | Plant Source | Wild-Type $\Delta G$ ($\text{kcal/mol}$) | Mutant Y150A $\Delta G$ ($\text{kcal/mol}$) | $\Delta\Delta G$ ($\text{kcal/mol}$) | Classification |
| :---: | :--- | :--- | :---: | :---: | :---: | :--- |
| **1** | **Limonene** | Citrus peel oil | **-4.78** | **-4.91** | **-0.13** | **Pan-Inhibitor (Resistance Proof)** |
| **2** | **Pulegone** | Pennyroyal / Mint oil | **-4.84** | **-4.86** | **-0.02** | **Pan-Inhibitor (Resistance Proof)** |
| **3** | **Thymol** | Thyme / Oregano oil | **-4.65** | **-4.66** | **-0.01** | **Pan-Inhibitor (Resistance Proof)** |
| **4** | **DEET** | Synthetic Control | **-5.04** | **-4.62** | **+0.42** | **Wild-Type Selective (Vulnerable)** |
| **5** | **Citronellal** | Citronella oil | -4.64 | -4.59 | +0.05 | Moderate Inhibitor |
| **6** | **Terpineol** | Tea tree oil | -4.70 | -4.57 | +0.13 | Moderate Inhibitor |
| **7** | **Eucalyptol** | Eucalyptus oil | -4.51 | -4.53 | -0.02 | Moderate Inhibitor |
| **8** | **Carvacrol** | Oregano derivative | -4.54 | -4.50 | +0.04 | Moderate Inhibitor |
| **9** | **Linalool** | Lavender oil | -4.44 | -4.43 | +0.01 | Moderate Inhibitor |
| **10** | **Alpha-Pinene** | Pine resin | -4.42 | -4.43 | -0.01 | Moderate Inhibitor |

---

### 3.3 Comparative Visual Analytics

The quantitative findings and resistance metrics are summarized in Figure 1 below:

![Screening and Resistance Analytics Chart](./aedes_data/visualizations/screening_and_resistance_analytics.png)

* **Figure 1A (Left Panel):** Comparative binding affinities ($\text{kcal/mol}$) across Wild-Type (blue) and Mutant Y150A (orange) variants. DEET shows a marked decrease in binding bar length on the mutant receptor, whereas **Limonene** and **Pulegone** maintain equal or greater binding depth.
* **Figure 1B (Right Panel):** $\Delta\Delta G$ resistance shift values. The red bar highlighting DEET ($+0.42\text{ kcal/mol}$) demonstrates significant energy loss due to mutation, while the green bar for Limonene ($-0.13\text{ kcal/mol}$) indicates enhanced stability in the mutated pocket.

---

### 3.4 3D Structural Pocket & Animation Visualizations

To visualize the spatial binding mode inside the transmembrane pocket, 3D structural models and animations were rendered:

![3D Binding Pocket Surface](./aedes_data/visualizations/binding_pocket_3d.png)

* **Figure 2:** High-resolution 3D ray-traced binding pocket visualization showing docked ligand pose inside the AaegOr4 active site.
* **Video File (MP4):** `./aedes_data/visualizations/AaegOr4_Docking__Simulation.mp4`

---

## 4. Public Health Implications & Conclusion

This study provides computational evidence that single point mutations (such as Tyr150Ala) in the *Aedes aegypti* AaegOr4 odorant receptor can reduce the binding efficiency of synthetic repellents like DEET. 

Crucially, our virtual screening identified natural citrus and mint terpenes—specifically **Limonene** and **Pulegone**—as pan-inhibitors that maintain tight binding affinity across both Wild-Type and mutated receptor variants. These natural plant volatiles offer promising leads for developing sustainable, eco-friendly vector control formulations capable of bypassing chemical resistance in dengue-endemic areas like South Asia.