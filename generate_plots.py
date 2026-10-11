import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Directory Setup
DATA_DIR = "./aedes_data"
RESULTS_DIR = os.path.join(DATA_DIR, "docking_results")
VIS_DIR = os.path.join(DATA_DIR, "visualizations")
os.makedirs(VIS_DIR, exist_ok=True)

# Load Leaderboard Data
csv_path = os.path.join(RESULTS_DIR, "dual_variant_screening_leaderboard.csv")
df = pd.read_csv(csv_path)

# Filter out Cinnamaldehyde docking anomaly (> 0 kcal/mol) for clean scale
df_clean = df[df["Wild_Type_kcal_mol"] < 0].copy()
df_top = df_clean.head(10)  # Top 10 natural compounds + DEET

# Set Scientific Styling
sns.set_theme(style="whitegrid", font="sans-serif")
fig, axes = plt.subplots(1, 2, figsize=(16, 7))

# ---------------------------------------------------------
# PANEL A: Binding Affinity Comparison (Top Volatiles)
# ---------------------------------------------------------
df_melted = df_top.melt(id_vars=["Compound"], value_vars=["Wild_Type_kcal_mol", "Mutant_Y150A_kcal_mol"],
                        var_name="Receptor_Variant", value_name="Binding_Affinity")

df_melted["Receptor_Variant"] = df_melted["Receptor_Variant"].replace({
    "Wild_Type_kcal_mol": "Wild-Type (AaegOr4)",
    "Mutant_Y150A_kcal_mol": "Mutant (Y150A)"
})

sns.barplot(data=df_melted, x="Binding_Affinity", y="Compound", hue="Receptor_Variant",
            palette=["#1f77b4", "#ff7f0e"], ax=axes[0])

axes[0].set_title("A: Binding Affinity Comparison (Top Volatiles)", fontsize=14, fontweight="bold")
axes[0].set_xlabel("Binding Affinity (kcal/mol) [More Negative = Stronger]", fontsize=12)
axes[0].set_ylabel("Volatile Organic Compound", fontsize=12)

# Zoom in on meaningful docking range (-6.0 to -3.5 kcal/mol)
axes[0].set_xlim(-5.5, -3.5)
axes[0].axvline(-5.0, color="red", linestyle="--", alpha=0.7, label="Strong Binding Threshold (-5.0)")

# Place Legend SAFELY outside the plot area at the top
axes[0].legend(loc="upper center", bbox_to_anchor=(0.5, 1.15), ncol=3, frameon=True)

# ---------------------------------------------------------
# PANEL B: Delta Delta G Resistance Shift
# ---------------------------------------------------------
colors = ["#d62728" if x > 0.2 else ("#2ca02c" if x < -0.05 else "#7f7f7f") for x in df_top["Delta_Delta_G"]]
sns.barplot(data=df_top, x="Delta_Delta_G", y="Compound", palette=colors, ax=axes[1])

axes[1].set_title("B: Resistance Energy Loss (ΔΔG = Mutant - Wild)", fontsize=14, fontweight="bold")
axes[1].set_xlabel("ΔΔG (kcal/mol) [Positive = Resistance / Loss of Affinity]", fontsize=12)
axes[1].set_ylabel("")
axes[1].set_xlim(-0.3, 0.6)
axes[1].axvline(0.0, color="black", linestyle="-", linewidth=1)

# Annotate DEET Resistance on Chart
deet_row = df_top[df_top["Compound"] == "DEET"]
if not deet_row.empty:
    deet_ddg = deet_row["Delta_Delta_G"].values[0]
    axes[1].text(deet_ddg + 0.02, df_top[df_top["Compound"] == "DEET"].index[0], 
                 f"DEET Resistance (+{deet_ddg})", va="center", fontweight="bold", color="#d62728")

plt.tight_layout()
out_plot = os.path.join(VIS_DIR, "screening_and_resistance_analytics.png")
plt.savefig(out_plot, dpi=300, bbox_inches="tight")
print(f"[+] Saved clean publication analytics chart to {out_plot}")