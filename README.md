# Bayesian Inverse Heat Transfer Paper (Data Assimilation Framework)

This repository provides all the necessary code and data to reproduce the results presented in our paper on Bayesian inverse heat transfer using data assimilation. It includes Gaussian and Multiquadric RBF-based reconstructions, visualization scripts, and reproducible containerized environments using **Docker** or **Singularity**.

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

### 1. Clone this repository
```bash
git clone https://github.com/KabirBakhshaei/bayesian-inverse-heat-transfer-paper
cd bayesian-inverse-heat-transfer-paper
path_files=$(pwd)
```

## If You're Using Docker
### 2A. Pull the Docker image
```
docker pull ithacafv/openfoam2106-muq2-pytorch
```
### 3A. Start a Docker container and mount the repo folder
```
docker run -it -v $path_files:/data/paper_respository ithacafv/openfoam2106-muq2-pytorch:latest bash
```
## If You're Using Singularity (e.g., SISSA Workstations)
### 2B. Load Singularity
```
module load singularity
```
### 3B. Pull the Docker image as a Singularity .sif file
```
singularity pull docker://ithacafv/openfoam2106-muq2-pytorch
```
### 4B. Run the container and mount the paper repo
```
singularity shell --bind $path_files:/data/paper_repository openfoam2106-muq2-pytorch_latest.sif
```
## Inside the Container (Linux Environment)
### 5. Clone and build ITHACA-FV with MUQ support
```
cd /data/paper_respository
git clone --depth 1 https://github.com/ITHACA-FV/ITHACA-FV
cd ITHACA-FV
git submodule update --init
source etc/bashrc
./Allwmake -tauq
```
### 6. Navigate back to the mounted repo
```
cd ../paper_respository
```
Once data is generated, postprocessing should be done outside Docker/Singularity using Python, and MATLAB.



