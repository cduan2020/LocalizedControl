# Localized algorithms for controlling network systems
This repository contains MATLAB code for controllability analysis, driver placement, and optimal control design on large-scale dynamical networks using the concept of information neighborhood. The algorithms, underlying theory, and application example details will be described in an upcoming paper, <i>Prevalence and scalable control of localized networks</i> by C. Duan, T. Nishikawa, and A. E. Motter.

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.


The full text of the GNU General Public License can be found in the file "LICENSE.txt".


# Dependency

The control design computations require MATLAB Control System Toolbox. Some of the application examples leverage MATLAB Parallel Computing Toolbox. 
The power grid control example requires Matpower 6.0 (https://matpower.org/download/) for power flow calculation. 


# Usage

### Scripts for example applications:
(Note: The run times for these examples were measured on a desktop computer with Intel i7-7700K / 4 cores / 16GB RAM.)

* `driverplacement.m` runs the gradient-based greedy algorithm for driver placement (Algorithm 2) on a Barabási–Albert (BA) random network with Laplacian dynamics. (Run time: ~1 min)

* `run_kuramoto_control.m` runs the local control design algorithm (Algorithm 3) and dynamical simulations on the network of Kuramoto oscillators. The folder `Kuramoto` contains supporting code and data files. `newtonraphson.m` in this folder was obtained from: asmaa lakhdar (2021). newtonraphson (https://www.mathworks.com/matlabcentral/fileexchange/87884-newtonraphson), MATLAB Central File Exchange. (Run time: ~1 min)

* `run_power_control.m` runs the local control design algorithm (Algorithm 3) and dynamical simulations on a synthetic Texas power grid (`case_ACTIVSg2000.m`, obtained from https://electricgrids.engr.tamu.edu). The code requires Matpower 6.0 to be installed on the MATLAB search path. The folder `Power grid` contains supporting code and data files. (Run time: ~2 min)

* `run_epidemic_control.m` runs the local control design algorithm (Algorithm 3) and dynamical simulations on the epidemics spreading over the global airline transportation network. The folder `Epidemics` contains supporting code and data files. This script uses MATLAB Parallel Computing Toolbox. (Run time: ~22 min)

* `run_brain_control.m` runs the local control design algorithm (Algorithm 3) and dynamical simulations on the whole brain network. The folder `Brainnet` contains supporting code and data files. This script uses MATLAB Parallel Computing Toolbox. (Run time: ~3.5 min)

### Utility functions:

* `Wup = InfoDistanceUpperbdd(normA)` calculates the upper bound on the information distance matrix from the norm of the blocks of the system coupling matrix A, i.e., `normA(i,j) = norm(A_{ij})`.

* `[nlist, dlist] = ucs_geodesic_k(Wup,i,k)` calculates the information neighborhood of size `k` centered at node `i` (Algorithm 1). The output `nlist` is the node list for the information neighborhood, and `dlist` is the corresponding list of information distances.

* `[nlist, dlist] = ucs_geodesic_tau(Wup,i,tau)` calculates the information neighborhood of radius `tau` centered at node `i` (Algorithm 1). The output `nlist` is the node list for the information neighborhood, and `dlist` is the corresponding list of information distances.

* `k = DesignLocalContrl(A,B,Q,R,Cset_aug0,Cset_aug,i,n)` runs the algorithm to design a local optimal controller at node `i` by solving the projected Riccati equation (a subroutine of Algorithm 3).

* `averge_neighborhood_rate = average_k_rate(normA,k)` calculates the average `k`th neighbor reduction rate of the system coupling matrix A, i.e., `normA(i,j) = norm(A_{ij})`.

* `PriorityQueue.m` provides the priority queue data structure used in the UCS algorithm (Algorithm 1).

* `ba_net.m` generates BA random networks. This code was obtained from: Gregorio Alanis-Lobato (2021). CNM (https://www.mathworks.com/matlabcentral/fileexchange/45734-cnm), MATLAB Central File Exchange.


