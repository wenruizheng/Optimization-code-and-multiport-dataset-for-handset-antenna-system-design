# Optimization Code and Multiport Dataset for Handset Antenna System Design

## Overview
This repository provides MATLAB implementation of multi-objective optimization algorithms and a multiport dataset for **handset antenna system design**.  
The framework enables automated optimization and evaluation of antenna systems under realistic device constraints, supporting reproducibility and further research.

## Features
- Multi-objective optimization using NSGA-III and GA-based operators.  
- Support for multiport network evaluation with scattering matrix data.  
- Includes pre-computed datasets (`.mat` files) for quick testing and benchmarking.  
- Configurable number of ports and feeds for flexible handset antenna scenarios.  
- Example `main.m` script to reproduce optimization results.

## File Structure
multiportrdata_farfield/
│── main.m # Main script to run optimization
│── CalObj.m # Objective function calculation
│── GA.m # Genetic Algorithm operators
│── EnvironmentalSelection.m # NSGA-III environmental selection
│── TournamentSelection.m # Tournament selection operator
│── NDSort.m # Non-dominated sorting
│── UniformPoint.m # Reference point generation
│── funfun.m # Initial population generator
│── obj_key.m # Decode chromosome to topology
│── Matrix&Target.mat # Pre-computed matrix & targets
│── num_feeds.mat # Predefined configuration
│── synthesized_pattern_3d.mat # 3D synthesized pattern data


## Requirements
- MATLAB R2020a or later  


   git clone https://github.com/wenruizheng/Optimization-code-and-multiport-dataset-for-handset-antenna-system-design.git
   cd Optimization-code-and-multiport-dataset-for-handset-antenna-system-design
