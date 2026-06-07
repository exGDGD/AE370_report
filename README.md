# AE370 Report

Code for the AE370 Miura optimization project.

## Requirements

- MATLAB

## Usage

Run the main MATLAB script from this folder:

```matlab
run_miura
```

Or run it from a terminal:

```powershell
matlab -batch run_miura
```

The script adds the local `matlab/` helper folder to the MATLAB path, performs
the coarse search, refines the best layout, and samples a local objective
landscape around the final best design.

Search ranges, optimizer limits, and output settings are configured near the
top of `run_miura.m`.

## Outputs

Generated results are written to `result/`, including:

- search ranking CSV files
- saved MAT result files
- best-design summaries
- convergence and landscape figures

The `result/` folder is ignored by git because it contains generated outputs.
