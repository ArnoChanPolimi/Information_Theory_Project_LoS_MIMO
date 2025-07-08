# Information Theory Project: LoS MIMO Simulation

## Project Overview
This MATLAB project implements simulations based on the Line-of-Sight (LoS) MIMO channel model in information theory.  
The project encapsulates key functions such as channel modeling and water-filling algorithm to compute channel capacity, enabling easy reproduction and extension of the results.

## Directory Structure
Information_Theory_Project_LoS_MIMO/
├── main.m % Main script to run the simulation
├── utils/ % Utility functions folder including channel model, water-filling algorithm, etc.
│ ├── water_filling.m
│ ├── channel_model.m
│ └── ... % Other helper functions
├── results/ % (Optional) Folder to store simulation results and plots
└── README.md % This documentation file

## Code Framework
- **main.m**: The entry point of the project. It initializes parameters, calls functions from `utils`, runs the simulation pipeline, and outputs results.
- **utils/**: Contains all auxiliary functions, including:
  - Channel matrix generation (`channel_model.m`)
  - Water-filling power allocation algorithm (`water_filling.m`)
  - Other mathematical helper functions
- The project follows a modular design for clarity and ease of maintenance and extension.

### Data Flow Overview
Parameter setup → Channel matrix generation → Capacity computation → Result output and visualization

## How to Run
1. Clone the repository:
```bash
git clone https://github.com/ArnoChanPolimi/Information_Theory_Project_LoS_MIMO.git
