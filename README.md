# Specialized Parallel Pathways for Adaptive Control of Visual Object Pursuit

This repository accompanies the manuscript:

**Collie, M.F., Jin, C., Kellogg, E., Vanderbeck, Q.X., Hartman, A.K., Hol, S.L., & Wilson, R.I.**  
*Specialized parallel pathways for adaptive control of visual object pursuit*  
[bioRxiv, 2025.](https://doi.org/10.1101/2025.04.23.650240)

---

## Overview

This project explores the circuit-level mechanisms underlying adaptive visual object pursuit in Drosophila melanogaster. The study demonstrates that pursuit behavior relies on two functionally distinct and anatomically segregated feedback pathways through the anterior optic tubercle (AOTU): AOTU025 provides coarse, fixed-gain steering for objects located in the lateral visual field, while AOTU019 provides precise, flexible-gain steering for objects near the midline. By combining connectomic analysis, whole-cell electrophysiology, optogenetic activation, and closed-loop behavioral assays, the authors show that AOTU019 exhibits direction-selectivity, is modulated by arousal and forward locomotion, and is sufficient to drive turning and forward acceleration. Modeling work further reveals how direction-selective inhibition from AOTU019 and excitation from AOTU025 jointly generate smooth, adaptive steering behavior. The repository contains all code for reconstructing these analyses and simulations, including connectomic data parsing, neural data acquisition, biophysical modeling, and statistical analysis of pursuit behavior across conditions.
---

## Repository Contents

The repository is organized into the following directories:

### 1. `Connectomics_Code/`
Code for analyzing synaptic connectivity data from the FAFB/FlyWire Drosophila brain connectome. Includes:
- Extracting and quantifying LC10a → AOTU019/AOTU025 → DNa02 synaptic pathways.
- Estimating visual receptive fields based on lobula input distributions.
- Mapping retinotopic organization by projecting medial-lateral lobula synapse positions to anterior-posterior retinal coordinates.
- Analyzing synapse counts, cable morphology, and connectivity strength across parallel visual-motor channels.

### 2. `Data_Acquisition_Code/`
Code used for behavioral and electrophysiological data collection. Includes:
- Integration with the G4 modular LED display system for open- and closed-loop visual stimulation.
- Real-time spherical treadmill tracking using FicTrac to measure forward, lateral, and rotational locomotion.
- Control of optogenetic arousal stimulation (P1 activation) and synchronization with behavioral output.
- Patch-clamp electrophysiology setup interfacing with visual stimuli and locomotion feedback.

### 3. `Modeling_Code/`
Mathematical models of pursuit behavior, including:
- A dynamical network model combining contralateral inhibition (AOTU019) and ipsilateral excitation (AOTU025).
- Simulation of pursuit behavior under varying error, object motion direction, arousal state, and forward velocity.
- Direction selectivity and gain scaling modules to test how circuit specializations improve pursuit performance.
- Variants that include random object motion, step perturbations, and genetic manipulations (e.g., AOTU019 silencing).
- Tools for quantifying model performance metrics such as settling time and midline fixation accuracy.

### 4. `Data_Analysis_Code/`
Data analysis pipelines for:
- Spike detection, membrane voltage filtering, and conversion of locomotor signals into kinematic variables.
- Binning and aligning data to visual stimuli for calculating average responses, tuning curves, and latencies.
- Computation of direction selectivity indices (DSIs), pursuit index metrics, and right-minus-left asymmetry signals.
- Classification of behavioral epochs (e.g., Schmitt triggers on forward velocity).
- Analysis of current injection experiments and cross-correlation between neural activity and rotational velocity.
- Mixed-effects statistical models for within-animal comparisons and between-genotype group analyses.

---

## Requirements

- MATLAB 2023b+ for data acquisition and analysis
- R v4.4.2+ and RStudio v2024.12.0+ for connectome analyses
- Python 3.9+ for treadmill and display arena interfacing
- [FicTrac](https://github.com/rjdmoore/fictrac) for spherical treadmill tracking
- [Display_Tools](https://reiserlab.github.io/Modular-LED-Display/G4/) for modular LED display arena ("G4")

---

## Citation

If you use this code or data, please cite:

> Collie et al. (2025). *Specialized parallel pathways for adaptive control of visual object pursuit*. bioRxiv. https://doi.org/10.1101/2025.04.23.650240

---

## Contact

For questions or data/code requests, please contact:  
Rachel I. Wilson (rachel_wilson@hms.harvard.edu)

---

© 2025 The Authors. Distributed under a CC BY-NC-ND 4.0 License.
