# ThesisProject

## Overview

This project investigates stellar systems, focusing on the relationship between **orbital period (PERIOD)** and **equatorial rotational velocity (VBROAD)**.

Initial attempts to model a linear relationship between these variables revealed significant scatter in the data (see `Figures/scatter without solution.jpg`). Standard robust fitting techniques were insufficient to accurately determine the slope of the densest region in the distribution.

To address this, we developed a novel **Markov Chain Monte Carlo (MCMC)**-based algorithm that fits a **2D Gaussian distribution** to the data. The inclination of the fitted Gaussian provides an estimate of the underlying linear relationship between the variables.

A detailed description of the algorithm can be found in the following publication:
https://ui.adsabs.harvard.edu/abs/2025A%26A...701A.195H/abstract

---

## Project Structure

### Data

* **`EBs with vbroad on MS.csv`**
  Contains data for 1050 stellar systems, including:

  * Orbital period (`PERIOD`)
  * Rotational velocity (`VBROAD`)

---

### Core Modules

* **`global_variable.py`**
  Contains shared functions and global variables used across the project, including utilities for plotting and directory management.

* **`Gaussian_Function_Fitting.py`**
  Implements the core algorithm:

  * Accepts `PERIOD` and `VBROAD` as input vectors
  * Fits a 2D Gaussian using an MCMC approach
  * Extracts the slope from the Gaussian inclination
  * Generates diagnostic plots:

    * MCMC parameter chains
    * Corner plots
    * Fitted distribution vs. data

  Output figures are stored in:
  `Figures/MCMC plots`

  The final fitted result is shown in:
  `Figures/scatter with solution.pdf`

---

* **`Activate_Gaussian_Function.py`**

  * Applies filters to the dataset (e.g., radius, temperature constraints)
  * Calls the Gaussian fitting module
  * Generates additional plots used during the study (not all included in this repository)

---

### Figures

* **`Figures/MCMC figures`**
  Contains MCMC diagnostic plots:

  * Corner plots
  * Log chains
  * Parameter distributions

* **`Figures/SI hist.jpg`**
  Displays the distribution of the fitted linear parameters (slope and inclination).

---

## Key Idea

Instead of directly fitting a line to noisy, highly scattered data, this project models the **density distribution** of the data using a 2D Gaussian. The orientation of this distribution provides a more robust estimate of the underlying linear relationship.

---

## Requirements

* Python 3.x
* NumPy
* SciPy
* Matplotlib
* emcee (for MCMC)
* corner (for visualization)

---

## Usage

1. Prepare the dataset (`EBs with vbroad on MS.csv`)
2. Run:

   ```
   Activate_Gaussian_Function.py
   ```
3. Outputs (plots and fitted results) will be saved in the `Figures/` directory.

---

## Notes

* Some plots generated during the research phase are not included in this repository.
* File and folder names reflect the original research workflow and may be verbose.

---
