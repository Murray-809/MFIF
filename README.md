# Cascade Alignment and Group-Based Fusion for Multi-Focus Images with Severe Geometric Misalignment

Official MATLAB implementation of the paper framework.

> ### 📢 Journal Submission Notice
> This repository is the official implementation for the manuscript currently submitted to **The Visual Computer**. 
> 
> **Important:** As requested by the journal, we remind readers that this code is directly associated with our submission. If you use this code or find our research helpful, please cite the relevant manuscript.

---

## 🌟 Key Contributions
Our framework addresses geometric distortions and data redundancy through three main innovations:
1. **Multi-Scale Cascade Alignment:** Robustly rectifies large geometric distortions to ensure spatial consistency even in blurred regions.
2. **Group-Based Refinement Fusion:** Progressively eliminates data redundancy and optimizes intermediate fused images via hierarchical group fusion.
3. **Local Cross-Correlation (LCC):** Accurately detects focused regions to generate reliable decision maps and preserve fine-grained edge details.

---

## 🛠 Prerequisites
* **Software:** MATLAB (Tested on R2022a or later)
* **Required Toolboxes:**  Image Processing Toolbox


---

## 🚀 Quick Start
1. **Clone/Download** this repository.
2. Open MATLAB and navigate to the project directory.
3. Run the demo script to see the fusion results:
   ```matlab
   run('align_demo.m') and run('fusion_demo.m')
