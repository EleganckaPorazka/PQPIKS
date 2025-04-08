# PQPIK-S

## Introduction

This is a library implementing the Predictive Quadratic Programming Inverse Kinematics with Scaling (PQPIK-S) algorithm for solving the inverse kinematics of redundant robots, accounting for the joint constraints on position, velocity, and acceleration.

If you plan to use it in your paper, please cite the following:

- ### Ł. Woliński, M. Wojtyra, "Predictive inverse kinematics with trajectory scaling for redundant manipulators based on quadratic optimization", Mechanism and Machine Theory, Volume 209, 2025, 105988, ISSN 0094-114X, [DOI:10.1016/j.mechmachtheory.2025.105988](https://doi.org/10.1016/j.mechmachtheory.2025.105988).

- ### Ł. Woliński, M. Wojtyra, "An inverse kinematics solution with trajectory scaling for redundant manipulators", Mechanism and Machine Theory, Volume 191, 2024, 105493, ISSN 0094-114X, [DOI:10.1016/j.mechmachtheory.2023.105493](https://doi.org/10.1016/j.mechmachtheory.2023.105493).

- ### Ł. Woliński, M. Wojtyra, "A Novel QP-Based Kinematic Redundancy Resolution Method With Joint Constraints Satisfaction", IEEE Access, vol. 10, pp. 41023-41037, 2022, [DOI:10.1109/ACCESS.2022.3167403](https://doi.org/10.1109/ACCESS.2022.3167403).

The above papers describe the PQPIK-S method in detail.

## Dependencies

This project uses the [Eigen][eigen] and [qpOASES][qpoases] libraries.

[eigen]: http://eigen.tuxfamily.org/index.php?title=Main_Page
[qpoases]: https://github.com/coin-or/qpOASES

## Installation

First, download and install Eigen. Then, download and build qpOASES. Make sure that the compiled qpOASES library is in directory qpOASES_master/bin, and that qpOASES_master is in the same directory as PQPIKS (or modify CMakeLists.txt in PQPIKS).

## Example of installation

cd PQPIKS

mkdir build

cd build

cmake .. && make

## Example of usage

cd PQPIKS/examples

mkdir build

cd build

cmake .. && make

./example
