# MassList Creator

**A high-resolution mass spectrometry analysis tool based on the methodology of Mickwitz et al. (2025, AMT).**

![Julia Version](https://img.shields.io/badge/Julia-v1.11%2B-blue.svg)
![Platform](https://img.shields.io/badge/Platform-Windows%20%7C%20Linux%20%7C%20macOS-lightgrey)

## 📖 Overview

**MassList Creator** is a Julia-based GUI application designed to automate the generation of mass lists for Time-of-Flight Chemical Ionization Mass Spectrometry (ToF-CIMS).

It processes HDF5 spectral data, performs peak fitting, generates chemical formula candidates based on configurable constraints, and assigns formulas using isotopic pattern matching. The output is optimized for creating calibration files compatible with **Toftracer2** or other Tof-processing tools.

## ✨ Key Features

* 🧪 **Smart Formula Generation**: Generates comprehensive chemical formula candidates based on user-defined rules and constraints.
* 📊 **Peak Fitting**: Fits detected peaks to theoretical profiles to determine precise peak centers and intensities.
* 🔍 **Isotopic Matching**: Computes expected isotopic distributions for accurate peak assignment and validation.
* 📈 **Interactive Visualization**: Built with `GLMakie` for high-performance, real-time spectral plotting and reactive analysis.
* ❄️ **Adduct Support**: Preconfigured for common adducts in `H+` and `NH4+` ionization modes.

---

## 💻 System Requirements

* **Operating System**: Windows 10/11, Linux, or macOS.
    * *Verified on Windows 11 (x86_64) / Intel i9-14900HX*.
* **Julia Version**: **v1.10** or higher.
    * *Developed and tested on **v1.11.3***.
* **Hardware**: Discrete GPU recommended for smooth plotting performance with `GLMakie`.

---

## 🛠️ Installation (Project Mode)

This tool is designed to run in a self-contained **Project Environment**. This ensures reproducibility and prevents conflicts with your global Julia packages.

### 1. Install Julia
Download and install Julia v1.11+ from the official website: [https://julialang.org/downloads/](https://julialang.org/downloads/)

### 2. Set Up the Environment
1.  Open a terminal (Command Prompt or PowerShell) and navigate to the folder containing `Masslist_final.jl`.
2.  Start Julia:
    ```
    julia
    ```
3.  Enter **Package Mode** by pressing `]` (the closing bracket key). The prompt will change to `(@v1.11) pkg>`.
4.  Run the following commands to activate the environment, install dependencies, and **precompile**:

    ```julia
    activate .
    instantiate
    precompile
    ```

    > **⚠️ Important**: The `precompile` step is crucial. `GLMakie` is a large library; precompiling ensures the GUI launches quickly. This may take a few minutes the first time.

5.  Press `Backspace` to exit Package Mode and return to the `julia>` prompt.

---

## 🚀 Usage

### Method: Running from the REPL (Recommended for Devs)
If you have just finished the installation steps above and are still in the Julia REPL:

```
julia
# 1. Ensure the environment is active
using Pkg; Pkg.activate(".")

# 2. Run the application
include("Masslist_final.jl")
```

## ⚙️ Workflow

    Load Data:
        Click Load HDF5. Select your .h5 file.
        Note: The file must contain MassAxis and AvgSpectrum (or SumSpecs) datasets.

    Configure Parameters:
        Set the Mass Range (e.g., 17 - 400 Da).
        Select Ion Mode (NH4+ or H+).
        Adjust Noise Threshold.

    Analyze:
        Click A. Fit Peaks: Identifies peaks in the spectrum using Gaussian or custom shapes.
        Click B. Assign Formulas: Matches fitted peaks to chemical formulas based on mass accuracy and isotopic patterns.

    Export:
        Click Export Masslist to save the results as a .csv file formatted for Tofware import.

## 🧪 Advanced Configuration

To customize the chemical formula generation rules (e.g., to allow for larger molecules or different element ratios), you must edit the FormulaConfig struct directly in Masslist_final.jl.

Locate the FormulaConfig struct (approx. line 60):
```
Julia

 Base.@kwdef struct FormulaConfig
     # --- 1. Atomic Count Limits (Max number of atoms) ---
     max_C::Int = 40
     max_H::Int = 80
     max_N::Int = 2      # Limit to 3 to reduce noise (sufficient for atmospheric nitrates)
     max_O::Int = 30
     max_S::Int = 2
     max_Si::Int = 2

     # --- 2. Organic Rules (Applies when C >= 1) ---
     # H/C Ratio Limits
     org_min_HC::Float64 = 0.3   # Min H/C ratio (Low value covers PAHs like Coronene)
     org_max_HC::Float64 = 2.5   # Max H/C ratio (High value covers saturated chains)
     
     # Other Ratio Limits
     org_max_OC::Float64 = 3.0   # Max O/C ratio (Allows highly oxidized molecules)
     org_max_NC::Float64 = 1.0   # Max N/C ratio (Kept at 1.0 to allow C1 nitrates)

     # --- 3. Advanced Filtering (Valence & DBE) ---
     max_DBE::Float64 = 10.0     # Max Double Bond Equivalent (Strict limit as requested)
     filter_radicals::Bool = false # If true, enforces integer DBE (removes radicals)

     # --- 4. Inorganic Rules (Applies when C = 0) ---
     inorg_max_H::Int = 12       # Safety cap for H count in inorganic clusters
     inorg_must_have_S::Bool = true # If true, Inorganic formulas MUST contain Sulfur (Filters HNO3)
 end
```
Modify these values and restart the application to apply changes.
## ❓ Troubleshooting

    Slow Startup: If the window takes a long time to appear or shows "Not Responding" initially, ensure you have run the precompile command in the installation steps.

    OpenGL/Plotting Errors: If the plot area is blank or crashes, it may be an issue with GPU drivers (common on some Windows laptops). Ensure your graphics drivers are up to date.

    "Package not found": Make sure you activated the environment using activate . (in Pkg mode).

## ⚠️ Known Issues & Fixes

> **Important Note 1:** There is a potential bug in the current processingProject.jl in Toftracer2 or Peakfitter (https://github.com/weikou/TOF-Tracer2/tree/dev; https://github.com/lstark-uibk/Manual_Pyeakfitter).
>
> **Solution:** To run the program correctly, please use the test script located in **`test/GUI_MainTest.jl`**. This file contains the necessary patch to ensure correct processing, but it is needed to change the pathway of the masslist files.

> **Important Note 2:** To change the threshold value, press enter when you've finished.

> **Important Note 3:** To start and create with a new masslist, it is better to restart the app, which may be updated later with new functions.

## 📄 Reference

The methodology described in:
Mickwitz et al. (2025), Atmospheric Measurement Techniques (AMT).
