📘 Linear Algebra Solver (C++)
🚀 Overview

A console-based Linear Algebra Computational Solver developed in C++ that solves systems of linear equations using multiple numerical methods with step-by-step visualization.

✨ Features
Gaussian Elimination (with pivoting)
Gauss-Jordan Elimination
Cramer’s Rule (for square systems)
Determinant Calculation
Matrix Rank Analysis
Step-by-step solution tracing
Fraction and decimal input support (a/b format included)
Solution verification (AX = b check)
Handles:
Unique solutions
Infinite solutions
No solution cases
🧠 Concepts Used
Linear Algebra
Numerical Methods
Matrix Operations
Gaussian Elimination
System of Linear Equations
Algorithm Design
Input validation and parsing
🛠 Tech Stack
C++
STL (Vectors, Strings, I/O)
Chrono Library (performance timing)
📌 How It Works
User inputs system of equations
System builds augmented matrix [A | b]
User selects solving method:
Gaussian Elimination
Cramer’s Rule
Gauss-Jordan
Program computes solution step-by-step
Final result is verified using AX = b
📊 Sample Output
Solution:
x[0] = 2
x[1] = -1

Verification: AX = b ✔

⭐ Future Improvements
GUI version (Qt / Python frontend)
File saving/loading system
Graphical matrix visualization