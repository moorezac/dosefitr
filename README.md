# dosefitr

**dosefitr** is an R package for analysing longitudinal biological data, with a focus on growth dynamics and dose-response relationships. It provides a flexible and automated framework for non-linear regression, enabling robust and interpretable modelling of time-resolved experiments.

Originally developed for high-content live-cell imaging data, dosefitr is **platform-agnostic** and can be applied to any longitudinal dataset where growth, decay, or response trajectories are measured over time.

## **Key Features**

- **Flexible non-linear regression**
  - Built-in support for common models (sigmoidal, exponential, biphasic, etc.).
  - Easy extension to custom user-defined models.
  - Automatic model selection using BIC.
- **Longitudinal analysis**
  - Designed for time-resolved datasets rather than endpoint-only measurements.
  - Supports extraction of growth rates and doubling times.
- **Dose-response modelling**
  - Calculation of GR, NDR, and LGR metrics.
  - Improved comparability across models with different baseline growth rates.
- **Robust fitting**
  - Outlier detection (ROUT + robust regression).
  - Parameter constraints and shared parameters across groups.
- **Scalable workflows**
  - Works with tabular data from any source.
  - Compatible with high-throughput and batch-processing pipelines.

**Installation**

You can install incur via:

remotes::install_github("moorezac/dosefitr")

**Usage**

The principle goal of dosefitr is to provide a user-friendly pipeline to export, segment, and analyse IncuCyte-derived data.

However, the automated curve fitting functions are able to be used on any dataset, and includes option for shared parameters across groups, outlier detection, and upper/lower limits on parameters. See below for built-in options

<div align="center">
  <img src="images/example_curves.jpg" width="50%">
</div>

## Learn

See the vignettes.
