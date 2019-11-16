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

### Projection methods
To predict using a projection method, use the projection method to get the projection matrices. For example:

Us = ManTDA(Xs, classes, lower dims)

Where ´Xs´ contains the tensor observations, with observations running along the first mode. So if the observed tensors have n modes, the tensor ´Xs´ will have (n+1) modes. The vector ´classes´ contains class labels which must be given as the integers 1 and 2. The vector ´classes´ must be the same length as the first mode of ´Xs´ is long. The input ´lower_dims´ specifies the size of each observation in the projected space. Hence this vector must be of length n, i.e. the same number of modes as in the observed tensors.

Once optimal projection matrices, ´Us´, have been found, these can used to project observed tensors into a lower dimensional space. This can be done as follows, where Yi is the tensor to project:

nmodes = length(size(Yi))
for jmode = 1:nmodes 
    origdims = size(Yi);
    Yi = Us{jmode}'*matricizing(Yi, jmode);
    newdims  = origdims;
    newdims(jmode) = size(Us{jmode}, 2);
    Yi = unmatricizing(Yi, jmode, newdims);
end

Note that the above is both inefficient and unnecessarily complex. At the minimum, tmult (from CPandTucker by Morten Mørup should be used instead of multiplication and the matricizing and unmatricizing functions, also from CPandTucker).

Once training data has been projected into the smaller space, the scalar values in the tensors in that small space can be used to train a standard classification method such as logistic regression.

Test data can be projected using the projection matrices ´Us´ learned from the training data, resulting in scalar values that can be fed into the trained classifier to yield final classifications.


### Direct classification methods
Use bilinear_get_predictions, giving both training and test data as input, to train the model and get predictions for the test data. This function does not save the optimised parameters, so some code inspection and refactoring is necessary to save optimised parameters. 

# Development status

Development will soon be stalled (end of 2019) as I will no longer have a Matlab license. Help from someone with a license would be greatly appreciated.

## Tasks

* Make fit and predict methods consistent between method types

* Rewrite implementation of get_diagonal_indices to make it more efficient, e.g. by vectorising and leveraging linear indices

* Generalise bilinear_logreg.m, bilinear_logreg_tucker.m, and ManPDA_normsratio.m to work for input tensors of arbitrary size.

* Generalise methods to handle more than two classes

* Add objective-oriented layer on top of basic implementations such that objects can be fitted and then used to predict (model objects with fit and predict methods)

* Write efficient code to project all observations into a lower dimensional space

* Modularise fitting and predicting with the bilinear methods (bilinear_get_predictions.m)
