# Greedy Time-Optimal Collision Avoidance (GTO) CAMs

A complete software framework for computing the latest possible initiation time to perform a collision avoidance maneuver (CAMs).

## Overview

GTO CAM uses a greedy approach to minimize the initiation time of a maneuver before conjunction.

## System Requirements

- **MATLAB**: R2025b or newer
- **WSL (Ubuntu)** with the following:
  - DACE (Differential Algebra Computational Engine)
  - nlohmann-json3-dev
  - libeigen3-dev, libjsoncpp-dev, libdlib-dev

## Installation

### 1. Install Dependencies

Run the dependency checker script on Windows:

```bash
check_deps.bat
```

This script automatically:
- Detects MATLAB installation 
- Detects WSL installation
- Installs/verifies DACE in WSL
- Installs/verifies JSON library

### 2. Build the Project

Run the build script:

```bash
build_scp.bat
```

This compiles the C++ backend using CMake within WSL and generates the necessary binaries.

## Quick Start

1. **Configure the scenario** in `main.m` inside the `params` structure:
   
2. Define the conjunction using `generateInitShort.m` or `readCDM_CONGEN.m` or any custom function
   which outputs two structures `primary` and `secondary`.

3.  **Run the main script**:
   ```matlab
   main
   ```

4. **View results** via generated figures in postprocessing.

## High-Level Algorithm Flow

1. **Initialize**: Set up MATLAB path and load toolboxes
2. **Configure**: Assemble full configuration via `params` and `generateInitShort.m`
3. **Run Greedy optimization**: Run `!wsl ./build/bin/backSweep`
4. **Validate**: Check solution feasibility against a forward-propagation scheme
5. **Postprocess**: Generate visualizations and reports via `mainPostprocess()`

## Output

- **outSim**: structure with the relevant output of the optimization, such as control profile and minimum time 
- **Figures**: Trajectory plots, validation errors, control profiles, etc.
- **Logs**: Installation and build logs for troubleshooting

## Troubleshooting

- **WSL user detection failed**: Ensure Ubuntu WSL distribution is installed and accessible
- **DACE build errors**: Check `install_log.txt` for detailed error messages

## Authors

**Zeno Pavanello**  
E-mail: zeno.pavanello@polimi.it  
**Frank De Veld**
E-mail: frankdeveld@proton.me
Date: 2022–2026

## Compatibility

- MATLAB R2025b or newer
- Windows 10/11 with WSL2 (Ubuntu)
- GCC 64-bit (Linux backend)

## Credits
If you use this work, please cite the following pre-print:

- Z. Pavanello, F. D. Veld, and R. Armellin, ‘Time-Optimal Collision Avoidance Via a 
  Greedy Polynomial Backward Sweep’, 2026, arXiv:2606.01169. doi:10.48550/arXiv.2606.01169.
