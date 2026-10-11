import os
import glob
import pymol
from pymol import cmd
from PIL import Image
import imageio

# ==========================================
# DIRECTORY SETUP
# ==========================================
DATA_DIR = "./aedes_data"
PDB_DIR = os.path.join(DATA_DIR, "structures")
RESULTS_DIR = os.path.join(DATA_DIR, "docking_results")
VIS_DIR = os.path.join(DATA_DIR, "visualizations")
FRAMES_DIR = os.path.join(VIS_DIR, "temp_frames")

os.makedirs(VIS_DIR, exist_ok=True)
os.makedirs(FRAMES_DIR, exist_ok=True)

RECEPTOR_PDB = os.path.join(PDB_DIR, "AaegOr4_A0A1S7UE45.pdb")
TOP_LIGAND_PDBQT = os.path.join(RESULTS_DIR, "Limonene_docked.pdbqt")
DEET_LIGAND_PDBQT = os.path.join(RESULTS_DIR, "DEET_docked.pdbqt")

def build_3d_visualization():
    print("="*60)
    print("STEP 3: GENERATING 3D RAY-TRACED FIGURES, MP4 VIDEO & GIF")
    print("="*60)
    
    # Initialize PyMOL in headless mode
    pymol.finish_launching(['pymol', '-cq'])
    cmd.reinitialize()
    
    # 1. Load Receptor
    print(f"[*] Loading receptor PDB: {RECEPTOR_PDB}...")
    cmd.load(RECEPTOR_PDB, "AaegOr4")
    
    # 2. Load Ligand
    ligand_to_load = TOP_LIGAND_PDBQT if os.path.exists(TOP_LIGAND_PDBQT) else DEET_LIGAND_PDBQT
    ligand_label = "Limonene" if os.path.exists(TOP_LIGAND_PDBQT) else "DEET"
    
    print(f"[*] Loading docked pose ({ligand_label}): {ligand_to_load}...")
    cmd.load(ligand_to_load, "docked_ligand")
        
    # 3. Styling Protein Structure (Cinematic View)
    cmd.hide("everything", "all")
    cmd.show("cartoon", "AaegOr4")
    cmd.color("cyan", "AaegOr4")
    cmd.set("cartoon_smooth_loops", 1)
    
    # 4. Styling Ligand & Active Binding Pocket
    cmd.show("sticks", "docked_ligand")
    cmd.color("yellow", "docked_ligand")
    
    # Highlight Key Active Site Residues
    cmd.select("pocket_res", "byres (AaegOr4 within 6.0 of docked_ligand)")
    cmd.show("sticks", "pocket_res")
    cmd.color("salmon", "pocket_res")
    
    # Center camera and adjust lighting/perspective
    cmd.zoom("docked_ligand", buffer=8.0)
    cmd.orient("docked_ligand")
    cmd.set("ray_shadows", 1)
    
    # 5. Save High-Resolution Static PNG Figure
    img_path = os.path.join(VIS_DIR, "binding_pocket_3d.png")
    cmd.png(img_path, width=1920, height=1080, dpi=300, ray=1)
    print(f"[+] Saved high-resolution 3D figure to {img_path}")
    
    # 6. Render 360-degree Smooth Rotation Frames (60 frames for smooth motion)
    print("[*] Rendering 60 high-definition rotation frames...")
    frame_files = []
    num_frames = 60
    
    for i in range(num_frames):
        cmd.turn("y", 360 / num_frames)
        frame_path = os.path.join(FRAMES_DIR, f"frame_{i:03d}.png")
        cmd.png(frame_path, width=1280, height=720, ray=0)
        frame_files.append(frame_path)
        
    # 7. Assemble MP4 Video via imageio & ffmpeg
    print("[*] Encoding MP4 video file...")
    mp4_path = os.path.join(VIS_DIR, "AaegOr4_docking_3D_simulation.mp4")
    
    writer = imageio.get_writer(mp4_path, fps=20, codec='libx264', quality=8)
    for frame_p in frame_files:
        writer.append_data(imageio.v3.imread(frame_p))
    writer.close()
    print(f"[+] Successfully generated 3D rotation MP4 VIDEO at {mp4_path}")
    
    # 8. Assemble Animated GIF
    print("[*] Assembling animated GIF...")
    gif_path = os.path.join(VIS_DIR, "AaegOr4_docking_spin.gif")
    images = [Image.open(f) for f in frame_files]
    images[0].save(gif_path, save_all=True, append_images=images[1:], duration=50, loop=0)
    print(f"[+] Successfully saved GIF at {gif_path}")

if __name__ == "__main__":
    try:
        build_3d_visualization()
    except Exception as e:
        print(f"[-] Visualization error: {e}")