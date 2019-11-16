# Dependencies

* ManOpt for optimisation over manifolds in Matlab

* Tensor manipulation code by Morten Mørup, available at http://www.imm.dtu.dk/~mm/downloads/CPandTucker.zip

* immoptibox http://www2.imm.dtu.dk/projects/immoptibox/

* N-Way toolbox by Rasmus Bro & Claus A. Andersson, available at http://www.models.life.ku.dk/nwaytoolbox

* gridLegend_v1.4 - only used in plot_comparisons.m

* export_fig - only used in plot_comparisons.m

* Code assumes Matalb R2018b or higher


# Tensor input assumptions

## Observation mode
All methods assume that observations run along the first mode.

## Sizes allowed in implementations
* Arbitrarily sized input tensors
 * Optimisation over the Stiefel manifold with the Tucker structure enforced, and the trace of the ratio of two matrices as the objective function, ManTDA.m
 * Optimisation over the Stiefel manifold with the Tucker structure enforced, and the ratio of two matrix traces as the objective function, ManTDA_normsratio.m
 * Optimisation over the Stiefel manifold with the PARAFAC structure enforced, and the trace of the ratio of two matrices as the objective function, ManPDA.m
 * Discriminant Analysis Tensor Representation [1] (DATEReig, with the generalised eigenvalue problem formulation), DATEReig.m
 * Constrained Multilinear Discriminant Analysis [2] (CMDA), CMDA.m
 * Direct General Tensor Discriminant Analysis [2] (DGTDA), DGTDA.m
 * DATER (the same as DATEReig, but optimised using the standard  eigenvalue problem [2]), DATER.m
 * Higher Order Discriminant Analysis [4] (HODA), HODA.m

* Matrix observations only (that is, 3D input tensors)
  * Optimisation over the Stiefel manifold with the PARAFAC structure enforced, and the ratio of two matrix traces as the objective function,  ManPDA_normsratio.m
  * Bilinear discriminant component analysis [6] (BDCA, PARAFAC structure), bilinear_logreg.m
  * BDCA with Tucker structure, as described in [3] and [5], bilinear_logreg_tucker.m

# Usage

## Tests

To run the tests of the classification and projection methods, make sure the current working directory in Matlab is 'matlab'. Then open and run the file code/tests/script_to_run_classification_tests.m. This file will run all the tests defined in code/tests/test_classification_nd.m.

The other tests defined in code/tests can be run by opening the respective files and running the tests with the test run button.

## Simulation code

To run the simulation code from which figures were generated for the paper [3] (see README at repo root level), make sure the current working directory is 'matlab'. Then open and run code/run_comparisons.m.

To plot the results of the comparisons, make sure the current working directory is 'matlab'. Then open and run code/plot_comparisons.m (after run_comparisons completed successfully).

## Fit and predict

### Projection methods that optimise all modes at once by leveraging manifold optimisation
### Mode-alternating projection methods
### Direct classification methods


# Development status

Development will soon be stalled (end of 2019) as I will no longer have a Matlab license. Help from someone with a license would be greatly appreciated.

## Tasks

* Make fit and predict methods consistent between method types

* Rewrite implementation of get_diagonal_indices to make it more efficient, e.g. by vectorising and leveraging linear indices

* Generalise bilinear_logreg.m, bilinear_logreg_tucker.m, and ManPDA_normsratio.m to work for input tensors of arbitrary size.

* Generalise methods to handle more than two classes

* Add objective-oriented layer on top of basic implementations such that objects can be fitted and then used to predict (model objects with fit and predict methods)
