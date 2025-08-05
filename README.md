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

# If You're Using Docker
## 2A. Pull the Docker image

