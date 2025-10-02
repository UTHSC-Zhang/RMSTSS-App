

# RMSTdesign: Interactive Power and Sample Size Calculator

This repository contains the source code for the `RMSTdesign` web application, an interactive tool for designing modern clinical trials.

-----

## Access the Live Application

The easiest way to use the calculator is to access the live, deployed version directly in your web browser. No installation is required.

### **[âž¡ï¸ Access the Live Application Here](https://arnab96.shinyapps.io/uthsc-app/)**

-----

## About the Application

This web application provides a user-friendly, point-and-click interface for performing complex power and sample size calculations for clinical trials that use the **Restricted Mean Survival Time (RMST)** as a primary endpoint.

It is designed for biostatisticians, clinical trialists, and researchers who need to design studies using advanced statistical methods without writing R code from scratch.

### Key Features

  * **Interactive Data Upload**: Upload your own pilot dataset in `.csv` format.
  * **Support for Advanced Models**: Implements a wide range of modern statistical methods:
      * Linear IPCW Models
      * Additive & Multiplicative Stratified Models
      * Semiparametric GAMs for non-linear effects
      * Dependent Censoring (Competing Risks) Models
  * **Intuitive Controls**: Set all analysis parameters, such as the truncation time ($\\tau$), target power, and sample sizes, through a graphical interface.
  * **Rich Outputs**: Instantly generate and view results, including:
      * Kaplan-Meier survival plots
      * Power vs. Sample Size curves
      * Summary tables of effect sizes and results
      * A downloadable analysis report

-----

## The `RMSTdesign` R Package

This application is powered by the **`RMSTdesign`** R package, which contains all the underlying statistical functions. For users who want to integrate these calculations into their own R scripts, see the full documentation, or contribute to the project, please visit the official package repository.

  * **Package Website & Source Code**: [**Click Here**](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/articles/RMSTpowerBoost.html)

-----

## Feedback and Issues

If you encounter a bug or have a suggestion for improving the application, please open an issue on this repository's "Issues" tab.
