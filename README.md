# Finite-Difference-Method-Perturbation-Theory-Study
This project implements the Finite Difference Method (FDM) to solve the time-independent Schrödinger equation for both the Double Well and Quantum Harmonic Oscillator (QHO) potentials. Additionally, First-Order Perturbation Theory is applied to the QHO to analyze the effect of a quartic perturbation on the ground state and the first excited state.

## Contents 
 - $\textbf{C++ Files:}$ Includes header, implementation, and test files.
 - $\textbf{Python Jupyter Notebook:}$ Used for graph creation and analysis.
 - $\textbf{PDF Document:}$ Outlines the theory behind the project and provides analysis and discussion of the results.

 ## Project Overview
This project uses an Object-Oriented Programming (OOP) approach to achieve the following:

 - $\textbf{Discretization of the x-axis:}$ Handles the problem in one dimension.
 - $\textbf{Creation of the Kinetic Matrix:}$ Based on the FDM method.
 - $\textbf{Choice of Potential:}$ Allows for the investigation of different potentials (two potentials are available at the time of writing).
 - $\textbf{Calculation of Eigenvectors & Eigenvalues:}$ Utilizes the Eigen library for these calculations.
 - $\textbf{Data Output:}$ Results are written to text files for further analysis.

 ## Analysis and Graphical Representation
The results obtained from the C++ implementation are analyzed and graphed using Python's Matplotlib library in a Jupyter Notebook. This approach is preferred for its flexibility and the wide variety of graphical libraries available in Python.

 ## Comparison with Analytical Predictions
The Eigenvectors and Eigenvalues calculated using the FDM are compared to analytical predictions to verify the accuracy of the numerical methods used.

 ## First-Order Perturbation Theory
A first-order perturbation is applied to the Quantum Harmonic Oscillator. The impact of this perturbation on the first two eigenvalues and eigenvectors is analyzed and compared with the results obtained from the FDM approach.

## Dependencies
$\textbf{C++:}$
 - Eigen library for linear algebra operations.
$\textbf{Python:}$
 - Matplotlib for plotting graphs.
 - Jupyter Notebook for interactive analysis.

## Conclusion
This project demonstrates the application of the Finite Difference Method and First-Order Perturbation Theory to solve and analyze quantum mechanical systems. The results are validated against analytical predictions, providing a comprehensive understanding of the methods and their accuracy.
