# Controlling localized networks
This repository contains code for controllability analysis, driver placement, and optimal control design on large-scale dynamical networks using the concept of information neighborhood.

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.


The full text of the GNU General Public License can be found in the file "LICENSE.txt".


# Dependency


In the power grid control problem, Matpower 6.0 (https://matpower.org/download/) is used for power flow calculation.


# Usage

Wup = InfoDistanceUpperbdd(P): Caculate the upper bound on the information distance matrix from the system coupling matrix P.

[nlist, dlist] = ucs_geodesic_k(Wup,i,k): Calculate the information neighborhood of size k centered at node i. The output nlist is the node list for the information neighborhood, and dlist is the corresponding list of information distances.

[nlist, dlist] = ucs_geodesic_tau(Wup,i,tau): Calculate the information neighborhood of radius tau centered at node i. The output nlist is the node list for the information neighborhood, and dlist is the corresponding list of information distances.

driverplacement.m: An exmaple of applying the gradient-based greedy algorithm for driver placement (Algorithm 2) to a BA random network with Laplacian dynamics.

run_kuromoto_control.m: An exmaple of local control design and dynamical simulation for a Kuramoto oscillator network.

run_power_control.m: An exmaple of local control design and dynamical simulation for the dynamics of a synthetic Texas power grid.

run_epidemic_control.m: An exmaple of local control design and dynamical simulation for epidemic spreading dynamics over the global airline transportation network.

run_brain_control.m: An exmaple of local control design and dynamical simulation of brain network.


