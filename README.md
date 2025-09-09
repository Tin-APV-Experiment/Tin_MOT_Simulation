# Tin MOT Simulation

This repository contains raw data generated from the **Optical Bloch Equation simulation software** developed by [T.K. Langin](https://iopscience.iop.org/article/10.1088/1367-2630/acc34d) for molecular laser cooling, and adapted by **G. Zheng** and **J. Wang** for studying:

- Laser slowing of tin (Sn) atoms  
- Laser cooling and magneto-optical trapping (MOT) of Sn atoms

More details on the Sn laser cooling and trapping simulation results and an outline of how such a cold Sn atom experiment can be used to study **atomic parity violation (APV)** and **isotope shift physics** can be found in our manuscript ([arXiv:2509.04635](https://arxiv.org/abs/2509.04635)).

The repository also includes MATLAB analysis code used to analyze the raw data. This includes:
- Plot acceleration heat maps and curves  
- Calculate slowing efficiency  
- Calculate MOT temperature, size, and average state populations 

---

## Repository Contents

### Jianwei’s Files
Contain OBE simulation raw data for the following stages:
- Red MOT  
- Compressed red MOT (cMOT)  
- Blue MOT  
- Conveyor belt blue MOT (CB MOT)  

Also contains Matlab analysis code to analyze the OBE raw data for the stages mentioned above, which creates relevant figures in our paper.
Additionally contains the Julia files which run the OBE solver for red MOT, cMOT, regular blue MOT, and CB MOT stages.

The **cMOT** and **CB MOT** results were used in our tin laser cooling and trapping paper ([arXiv:2509.04635](https://arxiv.org/abs/2509.04635)).


### Geoffrey’s Files
Contain OBE simulation raw data for the following stages:
- Red MOT  
- Compressed red MOT (cMOT)  
- White light slowing (WLS)

Also contains Matlab analysis code to analyze the OBE raw data for the stages mentioned above, which creates relevant figures in our paper.
Additionally contains the Julia files which run the OBE solver for red MOT, cMOT, and WLS stages.

Note that 2 independent runs were executed for red MOT and cMOT stages as a consistency check. The actual runs in each stage whose data we used for the paper are:

(Capture) Red MOT: 20250831_2359

Compressed red MOT step 1: 20250901_1046

Compressed red MOT step 2: 20250831_1046

WLS contains 8 separate runs (all in the folder) which are interpolated together to simulate a converging slowing beam with 300 mW of power.

The **red MOT**, **cMOT**, and **WLS** results were also used in our tin laser cooling and trapping paper ([arXiv:2509.04635](https://arxiv.org/abs/2509.04635)).

---

## Contact

For further questions about the material in this repository, please reach out to the **co-first authors** of the paper:  

- Geoffrey Zheng: [geoffreyz@uchicago.edu](mailto:geoffreyz@uchicago.edu)  
- Jianwei Wang: [jwang695@jh.edu](mailto:jwang695@jh.edu)  

