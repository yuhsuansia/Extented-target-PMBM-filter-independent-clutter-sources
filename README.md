# Extented-target-PMBM-filter-independent-clutter-sources
This repository contains Matlab implementations of the following multi-target filters for extended targets, where the clutter is the union of Poisson-distributed clutter and a finite number of independent clutter sources:

* Poisson multi-Bernoulli mixture (PMBM) filter.
* Poisson multi-Bernoulli (PMB) filter.

These filters are described in

A. F. Garcia-Fernandez, J. L. Williams, L. Svensson and Y. Xia, "Poisson multi-Bernoulli mixture filter with general target-generated measurements and arbitrary clutter," in IEEE Transactions on Signal Processing, 2023 doi: 10.1109/TSP.2023.3278944.

https://arxiv.org/abs/2210.12983

The PMB approximation is described in 

Y. Xia, K. Granström, L. Svensson, M. Fatemi, A. F. García-Fernández and J. L. Williams, "Poisson multi-Bernoulli approximations for multiple extended object filtering," in IEEE Transactions on Aerospace and Electronic Systems, 2022 doi: 10.1109/TAES.2021.3111720.

The filters are evaluated using the generalised optimal subpattern-assignment (GOSPA) integrated with the Gaussian Wasserstein distance

A. S. Rahmathullah, A. F. Garcia-Fernandez, and L. Svensson, Generalized optimal sub-pattern assignment metric, in 20th International Conference on Information Fusion, 2017.

Video on GOSPA: https://www.youtube.com/watch?v=M79GTTytvCM

main.m runs the extended target PMBM/PMB filters with PMB clutter

NOTE: implementations assign2D.m and kBest2DAssign.m are copied from the GitHub repository https://github.com/USNavalResearchLaboratory/TrackerComponentLibrary.
