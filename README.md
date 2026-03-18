# Research on Edge Computing

[![Language](https://img.shields.io/badge/Language-MATLAB-orange.svg)](https://www.mathworks.com/products/matlab.html)

## 📖 Introduction

This repository focuses on **resource allocation and task offloading in Edge Computing networks**, providing theoretical studies and code demonstrations.

The project not only implements fundamental resource allocation algorithms but also places a strong emphasis on **algorithm optimization, computational acceleration**, and the **horizontal comparison of algorithm performance under various strategies**. It serves as an excellent reference and foundation for further research and development in the field of edge computing.

## ⚙️ Core Algorithms

All simulations and computational code in this repository are written entirely in **MATLAB**, featuring the implementation, improvement, and acceleration of the following core algorithms:

* **LODCO (Lyapunov Optimization-based Dynamic Computation Offloading)**: An algorithm designed to optimize resource scheduling while ensuring system stability.
* **Greedy Algorithm**: A fundamental local search algorithm used to quickly find sub-optimal solutions for resource allocation.
* **$\epsilon$-Greedy Algorithm**: Introduces an exploration mechanism to the traditional greedy algorithm, balancing "exploration" and "exploitation" to potentially escape local optima.

## 📂 Repository Structure

The project code is organized into different folders, each containing not just basic simulations, but also specific implementations for **algorithm acceleration and optimization**:

* `/LODCO/`: Contains the core code, simulation scripts, and **implementations of acceleration optimizations** for the classic LODCO algorithm.
* `/greedy_logco/`: Implementation of resource allocation algorithms utilizing a greedy strategy, **including specific optimization and acceleration mechanisms designed to improve execution efficiency**.
* `/eps_greedy_logco/`: Implementation of the $\epsilon$-greedy strategy, **with a strong focus on algorithmic acceleration schemes and performance enhancements under this approach**.
* `/eps_normal_diff/`: Code for differential comparison experiments under varying $\epsilon$ parameters, **as well as the validation and comparison of different acceleration methods**.

## 🚀 Research Focus & Current Limitations

* **Optimization & Acceleration**: A primary focus of this project is exploring and implementing various methods to accelerate the convergence and execution efficiency of the aforementioned algorithms.
* **Strategy Comparison**: By configuring different environments and strategies, it comprehensively compares the performance of the greedy algorithm and its variants in edge computing scenarios.
* **Known Issues**: Currently, some of the methods used to accelerate the algorithms introduce a certain bias in strategy selection. These optimization methods are still undergoing continuous testing and improvement. Discussions, Issues, and Pull Requests are highly welcome.

## 🛠️ Requirements

* **MATLAB** (A recent version is recommended to ensure matrix computation and simulation efficiency).



