# Graph-Vertex-Model
Graph Vertex Model for simulating tissue growth. 

Includes ET and TE transitions was well as cell divisions.

This script runs vertex-model tissue simulations using the GVM framework. 

It configures model parameters via command-line arguments, executes one or more stochastic simulations, and automatically recovers from failed runs by restarting from the last valid saved state. 

The script is intended to be used as a batch-simulation driver, for example in parameter sweeps or high-throughput numerical experiments.

Example Usage

python /path/to/src/test_gvm.py --nsim 1 --tmax 1000.0 --nmax 2000 --sigma 0.2 --drate 0.01 --tens 3.0 --lays 1 --g0b 1.0

Overview

The workflow of the script is:
- Parse simulation parameters from the command line.
- Load an initial tissue configuration (.vt3d file).
- Initialize a graph-based database for tracking tissue topology.
- Assign physical and numerical model parameters.
- Run the simulation and periodically write output.
- If the simulation fails: 
	- Reload the last valid output state
	- Restore all counters and simulation time
	- Resume the simulation automatically
	This ensures robust long-time simulations even when numerical instabilities occur.

Requirements

- Python ≥ 3.8 (Python 3.10 and 3.11 were used for the simulations)
- networkx
- argparse
A valid input file (.vt3d) - provided in the “/test” folder



Command-Line Arguments

All model parameters are provided via command-line flags:

Argument      |      Type      |      Description
- --nsim      |      int      |      Number of independent simulations
- --tmax      |      float      |      Maximum simulation time
- --nmax      |      int      |      Maximum allowed number of cells
- --sigma      |      float      |      Edge-force fluctuation amplitude (sigma)
- --drate      |      float      |      Cell division rate (k_div)
- --tens      |      float      |      Line tension between live and necrotic layers (Gamma_LNI)
- --lays      |      int      |      Number of live cell layers (lambda)
- --g0b      |      float      |      Boundary surface tension (Gamma)



Fixed Model Parameters

In addition to command-line arguments, the following parameters are set internally:

Parameter  Value  Description
- h      |      0.001      |      Time integration step
- dg0      |      1      |      Baseline line tension
- kM      |      1      |      Myosin turnover rate
- kV      |      100      |      Volume compressibility modulus
- outFreq      |      10.0      |      Output frequency

These can be modified directly in the script if needed.

Input and Output

Input

Initial tissue configuration: /test/input.vt3d

Output

Simulation output is written to: /test_output

Output files include saved tissue states that can be reused to restart simulations.

Automatic Failure Recovery

If a simulation step fails (e.g. due to numerical instability), the script converts the last valid output file back into an input file using output_to_input.py. 
It reloads the tissue from that file and restores: simulation time, file counters, event counters (ET, TE, division events), and continues the simulation from the last valid state. 
This loop continues until the simulation completes successfully.
