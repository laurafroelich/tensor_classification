# Purpose

This repository contains code for classifying tensors using different tensor projection and classification methods. Currently the implementations only support matrix observations, There will also be code for simulating data to allow comparisons between the methods.

# Methods
## Projection methods
The projection methods implemented are: Discriminant Analysis Tensor Representation [1] (DATEReig, with the generalised eigenvalue problem formulation), Constrained Multilinear Discriminant Analysis [2] (CMDA), Direct General Tensor Discriminant Analysis [2] (DGTDA), DATER (the same as DATEReig, but optimised using the standard  eigenvalue problem [2]). The repository also contains code for methods that  optimise the objective function rigorously on a manifold and allow for objective functions assuming the PARAFAC structure. 

## Direct classification methods

# References
[1] Shuicheng Yan, Dong Xu, Qiang Yang, Lei Zhang, Xiaoou Tang and Hong-Jiang Zhang, "Discriminant analysis with tensor representation," 2005 IEEE Computer Society Conference on Computer Vision and Pattern Recognition (CVPR'05), 2005, pp. 526-532 vol. 1.
doi: 10.1109/CVPR.2005.131

[2] Q. Li and D. Schonfeld, "Multilinear Discriminant Analysis for Higher-Order Tensor Data Classification," in IEEE Transactions on Pattern Analysis and Machine Intelligence, vol. 36, no. 12, pp. 2524-2537, Dec. 1 2014.
doi: 10.1109/TPAMI.2014.2342214

[3]
