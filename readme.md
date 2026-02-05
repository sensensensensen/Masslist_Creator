# Masslist Creator

**Masslist Creator** is a Julia-based toolkit designed for the automated generation, calibration, and management of mass lists for Time-of-Flight Mass Spectrometry (ToF-MS) data.

It streamlines the workflow of processing raw ToF data (`.h5`), identifying peaks, and matching them with chemical formulas (including MCM species), making it an essential tool for atmospheric science and chemical analysis.

## 📂 Project Structure

The repository is organized as follows:

- **`src/`**: Contains the core source code.
  - `Masslist_Creator.jl`: Main module entry point.
  - `Masslist_final.jl`: Finalized logic for mass list generation.
- **`Test_Tof_data/`**: Contains raw TOF-MS data files (`.h5`) used for testing and validation.
- **`Test_masslist/`**: Includes reference mass lists and compound databases (e.g., MCM species).
- **`Literatures/`**: Related academic references and supplementary materials.
- **`test/`**: Unit tests and debugging scripts (Contains the fix for the current bug).

## ⚠️ Known Issues & Fixes

> **Important Note:** There is a known bug in the current processingProject.jl in Toftracer2.
>
> **Solution:** To run the program correctly, please use the test script located in **`test/GUI_MainTest.jl`**. This file contains the necessary patch to ensure correct processing.

> To change the threshold value, press enter when you've finished.

## 🚀 Features

- **Automated Peak Finding**: Algorithms tailored for TOF data to identify peaks from HDF5 raw files.
- **Mass List Generation**: Create calibrated mass lists (`.csv`) ready for analysis software (e.g., Tofware).
- **MCM Integration**: Built-in support for mapping Master Chemical Mechanism (MCM) compounds.
- **Batch Processing**: Capable of handling multiple datasets and time-series data.

## 🛠️ Installation & Usage

### Prerequisites
- **Julia 1.11.3
- HDF5 support (via `HDF5.jl`)

### Setup

1. **Clone the repository:**
   ```bash
   git clone [https://github.com/sensensensensen/Masslist_Creator.git](https://github.com/sensensensensen/Masslist_Creator.git)
   cd Masslist_Creator

using Pkg
Pkg.activate(".")
Pkg.instantiate()

include("test/GUI_MainTest.jl")
