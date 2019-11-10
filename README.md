# Purpose

This repository contains code for classifying tensors using different tensor projection and classification methods. Most of the implementations support tensors with arbitrarily many modes.

In particular, implementations of four novel methods, described in "Rigorous optimisation of multilinear discriminant analysis with Tucker and PARAFAC structures" [3], are included and work for arbitrarily sized tensors.

Currently, the implementations are in Matlab code. The long term plan is to translate the code to other languages such as Python, R, Scala, and possibly others, and also to support tensors with arbitrarily many modes for all methods. I would be more than happy to review pull requests towards this goal.



# Methods
## Projection methods that optimise all modes at once by leveraging manifold optimisation

## Mode-alternating projection methods
The implemented projection methods that alternate between modes during optimisation are: 

* Discriminant Analysis Tensor Representation [1] (DATEReig, with the generalised eigenvalue problem formulation)
* Constrained Multilinear Discriminant Analysis [2] (CMDA)
* Direct General Tensor Discriminant Analysis [2] (DGTDA)
* DATER (the same as DATEReig, but optimised using the standard  eigenvalue problem [2])
* Higher Order Discriminant Analysis (HODA)

## Direct classification methods

# Dependencies

## Matlab

* ManOpt for optimisation over manifolds in Matlab

* Tensor manipulation code by Morten Mørup, available at http://www.imm.dtu.dk/~mm/downloads/CPandTucker.zip

* immoptibox http://www2.imm.dtu.dk/projects/immoptibox/

* N-Way toolbox by Rasmus Bro & Claus A. Andersson, available at http://www.models.life.ku.dk/nwaytoolbox

* Code assumes Matalb R2018b or higher


## Python

# Usage

## Matlab

# Development status

## Matlab

Development will soon be stalled (end of 2019) as I will no longer have a Matlab license. Help from anyone with a license would be greatly appreciated.

## Python

More focus here from 2020.


# References
[1] Shuicheng Yan, Dong Xu, Qiang Yang, Lei Zhang, Xiaoou Tang and Hong-Jiang Zhang, "Discriminant analysis with tensor representation," 2005 IEEE Computer Society Conference on Computer Vision and Pattern Recognition (CVPR'05), 2005, pp. 526-532 vol. 1.
doi: 10.1109/CVPR.2005.131

[2] Q. Li and D. Schonfeld, "Multilinear Discriminant Analysis for Higher-Order Tensor Data Classification," in IEEE Transactions on Pattern Analysis and Machine Intelligence, vol. 36, no. 12, pp. 2524-2537, Dec. 1 2014.
doi: 10.1109/TPAMI.2014.2342214

[3] Frølich, L., Andersen, T. & Mørup, M. Rigorous optimisation of multilinear discriminant analysis with Tucker and PARAFAC structures. BMC Bioinformatics 19, 197 (2018) doi:10.1186/s12859-018-2188-0 (https://rdcu.be/bWwIL)
