# CEMP

## Overview
Clustering Ensemble algorithm using optimized Multiobjective Particle (CEMP) in single-cell RNA-seq. CEMP is featured with several mechanisms:

- Multi-subspace projection method;

- Basic partition module in different subspaces;

- Transforming representation between clusters and particles;

- Clustering ensemble optimization.

## Main function
[CEMP](https://github.com/lqmmring/CEMP/blob/master/CEMP.m): Main CEMO algorithm consisting of the three steps.

[Mating selection of MPSO/D](https://github.com/lqmmring/CEMP/blob/master/MatingSelection.m): Mating selection in CEMP.

[NDSort](https://github.com/lqmmring/CEMP/blob/master/NDSort.m): Non-dominated sorting by efficient non-dominated sort.

[Operator](https://github.com/lqmmring/CEMP/blob/master/Operator.m):Particle swarm optimization in CEMP.

[TournamentSelection](https://github.com/lqmmring/CEMP/blob/master/TournamentSelection.m): Tournament selection processing.

[UniformPoint](https://github.com/lqmmring/CEMP/blob/master/UniformPoint.m): Generate a set of uniformly distributed points on the unit hyperplane.

[Classification](https://github.com/lqmmring/CEMP/blob/master/Classification.m): Classify solutions into sub-regions

## Requirement
MATLAB R2016b

## Data availability

[Data sets are provided in the directory.](https://github.com/ishspsy/project/tree/master/MPSSC)

## Example

% CEMP(path_data,index_of_similarity,CA,C1,C2,T_index,evaluation)

% [Results] = CEMP('Data_Treutlin.mat',1,1,1,1,2,50);

## Main compaired algorithm

[MPSS](https://github.com/ishspsy/project/tree/master/MPSSC)

[SinNLRR](https://github.com/zrq0123/SinNLRR)

[EMEP](https://github.com/lixt314/EMEP)

## Acknowlegement

[PlatEMO](https://github.com/BIMK/PlatEMO)

## Contact
cslqm@hit.edu.cn

## License
This project is licensed under the MIT License.

