# CEMP
Clustering Ensemble algorithm using optimized Multiobjective Particle (CEMP) in single-cell RNA-seq

### multi-subspace projection method 
It is used for mapping the original data to low-dimensional subspaces is applied in order to detect hidden data structure at both gene level and sample level. 
### basic partition module in different subspaces
It is utilized to generate clustering solutions. 
### transforming representation between clusters and particle
It is used to bridge the gap between the discrete clustering ensemble optimization problem and the continuous multiobjective optimization algorithm. 
### clustering ensemble optimization

## Requirement
Matlab

## Example
% > Matlab R2016b

% CEMP(path_data,index_of_similarity,CA,C1,C2,T_index)

% [Results] = CEMP('Data_Treutlin.mat',1,1,1,1,2);

## Acknowlege
[PlatEMO](https://github.com/BIMK/PlatEMO)

## Contact
cslqm@hit.edu.cn

