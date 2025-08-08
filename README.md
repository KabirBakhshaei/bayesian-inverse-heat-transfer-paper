# Bayesian Inverse Heat Transfer Paper (Data Assimilation Framework)

This repository provides all the necessary code and data to reproduce the results presented in our paper entitled **Optimized Bayesian Framework for Inverse Heat Transfer Problems Using Reduced Order Method**. It includes Gaussian and Multiquadric RBF-based reconstructions, visualization scripts, and reproducible containerized environments using **Docker** or **Singularity**.

---
## Repository Structure
```
.
├── Data_Assimilation_Gaussian_RBF/
│   ├── Files/
│   └── Results/
├── Data_Assimilation_Multiquadric_RBF/
│   ├── Files/
│   └── Results/
├── SupplementaryImages/
│   ├── Files/
│   └── Results/
├── requirements.txt
└── README.md
```

---

## Project Setup Instructions (Docker or Singularity)

This project uses a **pre-configured container** that includes:

- `OpenFOAM` for CFD simulations  
- `MUQ` for uncertainty quantification  
- `PyTorch`  
→ All are already installed and configured in the container.

---

##  Step-by-Step Instructions

### 1A. Clone this repository and move into the cloned repository
```bash
git clone https://github.com/KabirBakhshaei/bayesian-inverse-heat-transfer-paper
cd bayesian-inverse-heat-transfer-paper
```
### 1B. Save the absolute path of this folder into a shell variable (used for folder mounting)
```bash
path_files=$(pwd)
```

## If You're Using Docker
### 2A. Pull (download) the pre-built Docker image
The following code was copied from this [link](https://hub.docker.com/r/ithacafv/openfoam2412-muq2-pytorch).

```
docker pull ithacafv/openfoam2412-muq2-pytorch
```
### 3A. Create and start a Docker container based on the downloaded image and share the repo folder with that container
```
docker run -i -t -v $path_files:/data/paper_respository  ithacafv/openfoam2412-muq2-pytorch:latest bash
```
## If You're Using Singularity (e.g., SISSA Workstations)
### 2B. Load Singularity
```
module load singularity
```
### 3B. Pull the Docker image as a Singularity .sif file
```
singularity pull docker://ithacafv/openfoam2412-muq2-pytorch
```
### 4B. Run the container and mount the paper repo
```
singularity shell --bind $path_files:/data/paper_repository openfoam2412-muq2-pytorch_latest.sif
```

## Inside the Container (Linux Environment)
### 5A. Clone ITHACA-FV inside the container
```
cd /data/paper_respository
git clone --depth 1 https://github.com/ITHACA-FV/ITHACA-FV
```
### 5B.Compile ITHACA-FV with MUQ support
```
cd ITHACA-FV
git submodule update --init
source /usr/lib/openfoam/openfoam2412/etc/bashrc
# source etc/bashrc
./Allwmake -tauq
```
### 6. Navigate back to the mounted repo
```
cd ../paper_respository
```

### 7. Navigate to the Simulation Directory and Run the Solver 
```
cd Data_Assimilation_Multiquadric_RBF/Files/
```
Then load the required modules and compile/run the simulation:
```
source /usr/lib/openfoam/openfoam2412/etc/bashrc            # Load OpenFOAM environment (version 2212)
# module load muq   # Load MUQ module, already pre-installed and linked inside the Docker/Singularity image
source /data/paper_repository/ITHACA-FV/etc/bashrc # Then load ITHACA-FV environment
wclean
wmake             # Compile the solver
blockMesh
06enKFwDF_3dIHTP  # Run the solver 
```
For the following directory that uses **Gaussian RBF**, you need to follow some additional steps:

```
cd Data_Assimilation_Gaussian_RBF/Files/
```
Open the file
```ITHACA-FV-KF/src/ITHACA_FOMPROBLEMS/sequentialIHTP/sequentialIHTP.C```
inside the file, commant out the line marked ```heatFluxSpaceBasis[funcI][faceI] = Foam::sqrt(1 + (shapeParameter * radius) * (shapeParameter * radius));``` and uncommand the line marked with```heatFluxSpaceBasis[funcI][faceI] = Foam::exp(-1.0 * (shapeParameter * shapeParameter) * (radius * radius));```to switch the configuration to Gaussian RBF mode. 
This change toggles the reconstruction kernel used by the solver from multiquadric to Gaussian. Recompile ITHACA-FV
```
cd ~/ITHACA-FV-KF
of2212
source ~/ITHACA-FV-KF/etc/bashrc
./Allwmake -m -j 4
```
And run the simulation as done previously:

### 8. Generated Output Files and Folders After Simulation
**Folders:**
```
ITHACAoutput/
```
**Files:**
```
B.npy
Btemp.npy
condNumber.txt
condNumberAutoCovInverse.txt
condNumberCrossCov.txt
condNumberKalmanGain.txt
gTrueMatrix.npy
measNoiseCovTotal.txt
measurementsMat.npy
measurementsMatNoise.npy
parameterPriorCov.npy
parameterPriorMean.npy
parameterPriorMeanWithoutShifting.npy
radius_kb.npy
replay_pid2156.log
Temp.npy
Temp2.npy
thermocouplesCellsID_mat.txt
thermocoupleXValues.npy
thermocoupleYValues.npy
thermocoupleZValues.npy
xyz.npy
```

Once data is generated, postprocessing should be done outside Docker/Singularity using Python, and MATLAB.

## Postprocessing Guide

### Inside `Data_Assimilation_Multiquadric_RBF/Files/`

#### `plots.ipynb`

- **Reads the following as inputs:**
  - From `./ITHACAoutput/true/`:  
    - `trueTimeVec_mat.txt`  
    - `probe_true_mat.txt`
  - From `./ITHACAoutput/reconstruction/`:  
    - `probe_rec_mat.txt`  
    - `probeState_minConf_mat.txt`  
    - `probeState_maxConf_mat.txt`  
    - `gTrue_probe_mat.txt`  
    - `gRec_probe_mat.txt`  
    - `gRec_probeMaxConf_mat.txt`  
    - `gRec_probeMinConf_mat.txt`  
    - `parameterMean_mat.txt`  
    - `parameter_minConf_mat.txt`  
    - `parameter_maxConf_mat.txt`
  - From `./ITHACAoutput/projection/HeatFluxSpaceRBF/`:  
    - `heat_flux_space_basis_mat.txt`
  - From `./ITHACAoutput/projection/TrueHeatFlux/`:  
    - `HeatFluxTrue_mat.txt`
  - From current directory:  
    - `parameterPriorMean.npy`  
    - `parameterPriorMeanWithoutShifting.npy`  
    - `condNumberAutoCovInverse.txt`  
    - `condNumberCrossCov.txt`  
    - `condNumberKalmanGain.txt`  
    - `xyz.npy`  
    - `Temp.npy`

- **Outputs (to `../Results/`):**
  - Figures:  
    - `Figure 11B TrueAndReconstructedMeanTemperatureAtaProbe_0.91_0.02_0.55_OverTimeMultiquadric.png`  
    - `Figure 12B TrueAndReconstructedMeanHeatFluxAtProbe_0.91_0.0_0.55_over_time.png`  
    - `Figure 13B TrueAndReconstructedMeanHeatFluxAtTheHotSideWithConfidenceInterval.png`
  - Text files:  
    - `meanOfMeanRelativeError.txt`  
    - `meanRelativeErrorAtT0.txt`  
    - `meanRelativeErrorAtT0WithoutShifting.txt`  
    - `autoCovInverse_mean_std.txt`  
    - `crossCov_mean_std.txt`  
    - `kalmanGain_mean_std.txt`  
    - `xyz.txt`  
    - `Temp.txt`

#### `Surface3DAnimation.m`

- **Creates (to `../Results/`):**
  - `3D Combined surface plot.avi`
  - `3D Combined surface plot.gif`

#### `ContourFigurePaperMultiquadricRelative.m`

- **Creates (to `../Results/`):**
  - `Figure 14b SnapshotCountorsMultiquadric.png`  
  - `Figure 14C SnapshotCountorsTrue.png`  
  - `Figure 14e SnapshotCountorsMultiquadricRelative.png`

#### `RBFsThermocouples.m`

- **Creates (to `../Results/`):**
  - `Figure 2 Thermocouples_RBF_Centers.png`

---

### Inside `Data_Assimilation_Gaussian_RBF/Files/`

- Run:
  - `plots.ipynb`
  - `Surface3DAnimation.m`
  - `ContourFigurePaperGaussianRelative.m`

- Produces:
  - Figures analogous to the Multiquadric RBF setup:  
    `Figure 11A, 12A, 13A`–`Figure 14a, 14C, 14D`
  - `3D Combined surface plot.avi`
  - `3D Combined surface plot.gif`

---

### Inside `SupplementaryImages/Files/`

#### `RBFsComparison.m`

- **Reads (from local folder):**
  - `xyz.txt`
  - `TempG0_5.txt`, `TempG1.txt`, `TempG2.txt`, `TempG2_5.txt`
  - `TempM0_5.txt`, `TempM1.txt`, `TempM3.txt`, `TempM7_5.txt`

- **Creates (in `../Results/`):**
  - **Individual Surface Plots:**
    - `Figure 4a1`–`Figure 4a4` (Gaussian)  
    - `Figure 4b1`–`Figure 4b4` (Multiquadric)
  - **Combined Images:**
    - `Figure 4a Combined1.png` / `.pdf`  
    - `Figure 4b Combined2.png` / `.pdf`
  - **Face Centers:**
    - `Figure 4c XZCentersOfFaces.pdf`

---

## requirements.txt

The `requirements.txt` file specifies:

- Python version
- MATLAB version
- Required libraries
---

## Zenodo Archive

This repository is archived and citable via Zenodo.  
Click the badge below to access the DOI and download the versioned release:
In Progress..............

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)



