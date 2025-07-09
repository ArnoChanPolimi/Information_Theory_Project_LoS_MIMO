# Information Theory Project: LoS MIMO Simulation

## Project Overview
This MATLAB project implements simulations based on the Line-of-Sight (LoS) MIMO channel model in information theory.  
The project encapsulates key functions such as channel modeling and water-filling algorithm to compute channel capacity, enabling easy reproduction and extension of the results.

## Directory Structure
```
Information_Theory_Project_LoS_MIMO/
├── main.m                            % Main script to run the simulation
├── utils/                            % Utility functions folder
│   ├── create_H_matrix.m            % Generates channel matrix H based on system parameters
│   ├── los_mimo_capacity.m          % Computes LoS MIMO channel capacity (equal power allocation)
│   ├── water_filling_capacity_bisect.m % Implements water-filling capacity optimization
│   ├── compute_singular_values.m    % Computes and displays singular values and related diagnostics
│   ├── draw.m                       % Custom plotting functions
│   └── save_figure_custom.m         % Customized function to save figures
├── results/                          % Optional folder to save generated figures and data
└── README.md                         % This documentation file
```
### Code Framework

* **main.m**: The main entry point of the project. It initializes parameters, calls functions in `utils/`, runs the simulation pipeline, and outputs results.
* **utils/**: Contains all auxiliary functions, including:
  * Channel matrix generation (`create_H_matrix.m`)
  * LoS MIMO capacity calculation with **traditional equal power allocation** (`los_mimo_capacity.m`)
  * Water-filling power allocation algorithm (`water_filling_capacity_bisect.m`)
  * Singular value analysis and diagnostics: (`compute_singular_values.m`)
  * Custom plotting functions (`draw.m`)
  * Custom figure saving function (`save_figure_custom.m`)
* The project follows a modular design for clarity and ease of maintenance and extension.

### Data Flow Overview
Parameter setup → Channel matrix generation → Capacity computation → Result output and visualization

## How to Run
1. Clone the repository:
```bash
git clone https://github.com/ArnoChanPolimi/Information_Theory_Project_LoS_MIMO.git
