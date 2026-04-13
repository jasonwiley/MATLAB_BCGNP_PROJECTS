# MATLAB_BCGNP_PROJECTS

A collection of MATLAB tools for analyzing lamellae-forming block copolymer grafted nanoparticle (BCGNP) systems simulated using [HOOMD-blue](https://hoomd-blue.readthedocs.io/).

## Description

These tools are designed to process and analyze simulation output from HOOMD-blue, with a focus on characterizing lamellar microstructure in block copolymer grafted nanoparticle systems. Current tools include surface area analysis of rough vs. smooth lamellar interfaces.

## Requirements

- MATLAB R2023b
- Some tools require simulation data exported as `.mat` files from HOOMD-blue GSD output (see individual tool folders for details)

## Setup

1. Clone the repository:
```bash
   git clone https://github.com/jasonwiley/MATLAB_BCGNP_PROJECTS.git
```
2. Navigate to the tool of interest and refer to its folder for any data requirements.

## Usage

Each subfolder contains an independent analysis tool. Navigate to the relevant folder and run the main script:

```matlab
cd TOOL_FOLDER_NAME
run('main_script_name.m')
```

## Notes

- Tools that require simulation data expect a `.mat` file exported from a HOOMD-blue GSD output. These are noted in the relevant folders.
- Tools were developed for studying lamellar phase behavior in BCGNP systems.

## Author

Jason Wiley, M.S.  
Louisiana Tech University  
Department of Engineering Physics
