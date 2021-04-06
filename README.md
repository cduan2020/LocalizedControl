# Controlling localized networks
This repository contains code for controllability analysis, driver placement, and optimal control design on large-scale dynamical networks using the concept of information neighborhood.

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.


The full text of the GNU General Public License can be found in the file "LICENSE.txt".


# Dependency


In the power grid control problem, Matpower 6.0 (https://matpower.org/download/) is used for power flow calculation.


# Usage

`Wup = InfoDistanceUpperbdd(normA)`: Caculate the upper bound on the information distance matrix from the norm of the blocks of the system coupling matrix A, i.e., `normA(i,j) = norm(A_{ij})`.

`[nlist, dlist] = ucs_geodesic_k(Wup,i,k)`: Calculate the information neighborhood of size k centered at node i (Algorithm 1). The output nlist is the node list for the information neighborhood, and dlist is the corresponding list of information distances.

`[nlist, dlist] = ucs_geodesic_tau(Wup,i,tau)`: Calculate the information neighborhood of radius tau centered at node i (Algorithm 1). The output nlist is the node list for the information neighborhood, and dlist is the corresponding list of information distances.

`driverplacement.m`: An exmaple of applying the gradient-based greedy algorithm for driver placement (Algorithm 2) to a BA random network with Laplacian dynamics.

`run_kuromoto_control.m`: run the example for the control of Kuramoto oscillator.<br/>
`Kuromoto`: Folder containing supporting code and data files for local control design and dynamical simulation on a Kuramoto oscillator network.

`run_power_control.m`: run the example for the control of Texas power grid.<br/> 
`Power grid`: Folder containing supporting code and data files for local control design and dynamical simulation on the dynamics of a synthetic Texas power grid.

`run_epidemic_control.m`: run the example for the control of the epidemics spreading dynamics over the global airline transportation network.<br/> 
`Epidemics`: Folder containing supporting code and data files for local control design and dynamical simulation on the epidemic spreading dynamics over the global airline transportation network.

`run_brain_control.m`: run the example for the control of the whole brain network.<br/> 
`Brainnet`: Folder containing supporting code and data files for local control design and dynamical simulation on the whole brain network.

Other unitility functions：<br/>
`k=DesignLocalContrl(A,B,Q,R,Cset_aug0,Cset_aug,i,n)`: design local optimal controller at node i by solving the prejected Riccati equation (part of Algorithm 3).<br/>
`averge_neighborhood_rate = average_k_rate(normA,k)`: calculate the average kth neighbor reduction rate of the system coupling matrix A, i.e., `normA(i,j) = norm(A_{ij})`.<br/>
`PriorityQueue.m`: the priorigy queue data structure used in the USC algorithm (Algorithm 1).<br/>
`er_net.m`: generate ER random networks.<br/>
`ba_net.m`: generate BA random networks.<br/>

