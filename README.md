# Simulations and data analysis for *Causal Vaccine Effects on Post-infection Outcomes in the Naturally Infected*

This repository contains code for simulations and real data analysis accompanying the manuscript "Causal Vaccine Effects on Post-infection Outcomes in the Naturally Infected."

## Repository Structure

```
├── simulations/              # Simulation studies (Section 6)
│   ├── bash scripts          # Job submission scripts
│   ├── configuration files   # Simulation parameters
│   └── R/                    # Analysis and table generation scripts
│   └── results/              # Folder to hold results output
│   └── truth/                # Folder to hold truth output
│
└── real_data_analysis/       # PROVIDE data analysis
    ├── bash scripts          # Analysis job scripts
    ├── configuration files   # Analysis parameters
    ├── real_data_analysis.R  # Primary analysis
    └── sensitivity_analysis.R # Sensitivity analyses
```

## Simulations

### 6.1 Asymptotic Properties of Estimators

**Structure:**
- `run_simulation_1.sh` - Bash script to submit simulations for one-step estimators
- `config_sim_1.yml` - Configuration file containing settings used in simulation 1
- `run_simulation_1_bounds.sh` - Bash script to submit simulations for bounds
- `config_sim_1_bounds.yml` - Configuration file containing settings used in simulation 1 for bounds
- `R/` folder contains scripts for running analysis and generating tables (files named corresponding to "simulation_1" and "simulation_1_bounds")

**Outputs:**
- Supplement Section I.1, Tables 1-5

### 6.2 Comparing Power of Estimands in Realistic Setting

**Structure:**
- `run_simulation_contour.sh` - Bash script to submit simulations under PROVIDE-like data generating process
- `R/` folder contains scripts for running PROVIDE-like analysis and generating contour plot to compare power across settings (file names correspond to "contour")

**Outputs:**
- Main manuscript Figure 1

## Real Data Analysis

Analysis of data from the PROVIDE study 

### Primary Analysis

**Structure:**
- `real_data_analysis.R` - Rscript containing main analysis, covariate adjusted bounds, and table generation code

**Outputs:**
- Main manuscript Table 2 (primary results)
- Supplement Table 8 (covariate-adjusted bounds)

### Sensitivity Analysis

**Structure:**
- `sensitivity_analysis.R` - Sensitivity analysis script to read main results and make sensitivity analysis figure

**Outputs:**
- Supplement section J2, Figure 1

## Contact

For questions, submit an issue or email [allison.codi@emory.edu](mailto:allison.codi@emory.edu).

