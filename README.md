# LocalizedControl
Controllability analysis, driverplacement, and optimal control design on large-scale dynamical networks using the concept of information neighbourhood.

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.


The full text of the GNU General Public License can be found in the file "LICENSE.txt".


# Dependency

MATLAB 2017b or later.

In the power grid control problem, Matpower 6.0 https://matpower.org/download/ is used for power flow calculation.


# Usage

Wup=InfoDistanceUpperbdd(P) Caculate the upperbound of information distance matrix from the system coupling matrix P.

[nlist, dlist] = ucs_geodesic_k(Wup,i,k) calculate the information neighborhood of size k centered at node i. nlist is the node list of the information neighborhood and dlist is the corresponding list of information distances.

[nlist, dlist] = ucs_geodesic_tau(Wup,i,tau) calculate the information neighborhood of information radius tau centered at node i. nlist is the node list of the information neighborhood and dlist is the corresponding list of information distances.

driverplacement.m An exmaple of the gradient-based greedy algorithm for driverplacement on a BA random network with Lapacian dynamics.

run_kuromoto_control.m an exmaple of local control design and dynamical simulation of Kuramoto oscillator network.

run_power_control.m an exmaple of local control design and dynamical simulation of power grid.

run_epidemic_control.m an exmaple of local control design and dynamical simulation of epidemics over airline transportation network.

run_brain_control.m an exmaple of local control design and dynamical simulation of brain network.


