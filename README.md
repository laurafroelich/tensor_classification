# Purpose

This repository contains code for classifying tensors using different tensor projection and classification methods. Most of the implementations support tensors with arbitrarily many modes.

In particular, implementations of four novel methods, described in "Rigorous optimisation of multilinear discriminant analysis with Tucker and PARAFAC structures" [3], are included and work for arbitrarily sized tensors.

Currently, the implementations are in Matlab code. The long term plan is to translate the code to other languages such as Python, R, Scala, and possibly others, and also to support tensors with arbitrarily many modes for all methods. I would be more than happy to review pull requests towards this goal.



# Methods


The projection methods find optimal projections that separate observations from different classes maximally. Hence the projection matrices found by the projection methods can be used to project each higher-order observation into a smaller space, and its values in that space can be used for classification.

The direct classification methods perform projection and classification in one step, by evaluating an expression similar to that in logistic regression, but using tensor multiplication to arrive at a scalar value in the exponents. For the direct classification methods, we have only implemented versions that work for matrix observations.

## Projection methods that optimise all modes at once by leveraging manifold optimisation

We propose four ways to optimise a projection over a manifold in [3]. Two of the methods assume the Tucker structure and two of the methods assume the PARAFAC structure for interactions between modes. Furthermore, two of the methods optimise a ratio of two traces whereas the other two methods optimise the trace of the ratio of two matrices. For further details see [3] or [5].

## Mode-alternating projection methods
The implemented projection methods that alternate between modes during optimisation are: 

* Discriminant Analysis Tensor Representation [1] (DATEReig, with the generalised eigenvalue problem formulation)
* Constrained Multilinear Discriminant Analysis [2] (CMDA)
* Direct General Tensor Discriminant Analysis [2] (DGTDA)
* DATER (the same as DATEReig, but optimised using the standard  eigenvalue problem [2])
* Higher Order Discriminant Analysis [4] (HODA)

## Direct classification methods

* Bilinear discriminant component analysis [6] (BDCA, PARAFAC structure)
* BDCA with Tucker structure, as described in [3] and [5]

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

[4] Phan A. H, Cichocki A. Tensor decompositions for feature extraction and classification of high dimensional datasets. Nonlinear Theory Appl IEICE. 2010; 1(1):37–68

[5] Decomposition and classification of electroencephalography data. / Frølich, Laura.
Kgs. Lyngby : Technical University of Denmark, 2016. 208 p. (DTU Compute PHD-2016; No. 408). (https://orbit.dtu.dk/en/publications/decomposition-and-classification-of-electroencephalography-data(35de6ca8-d5a6-467b-88d9-a0d97bbe4685)/export.html)

[6] M. Dyrholm, C. Christoforou, and L. C. Parra, “Bilinear discriminant component analysis,” The Journal of Machine Learning Research, vol. 8, pp. 1097–1111, 2007.