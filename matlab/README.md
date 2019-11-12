# Dependencies

* ManOpt for optimisation over manifolds in Matlab

* Tensor manipulation code by Morten Mørup, available at http://www.imm.dtu.dk/~mm/downloads/CPandTucker.zip

* immoptibox http://www2.imm.dtu.dk/projects/immoptibox/

* N-Way toolbox by Rasmus Bro & Claus A. Andersson, available at http://www.models.life.ku.dk/nwaytoolbox

* gridLegend_v1.4 - only used in plot_comparisons.m

* export_fig - only used in plot_comparisons.m

* Code assumes Matalb R2018b or higher


# Usage

## Tests

To run the tests of the classification and projection methods, make sure the current working directory in Matlab is 'matlab'. Then open and run the file code/tests/script_to_run_classification_tests.m. This file will run all the tests defined in code/tests/test_classification_nd.m.

The other tests defined in code/tests can be run by opening the respective files and running the tests with the test run button.

## Simulation code

To run the simulation code from which figures were generated for the paper [3] (see README at repo root level), make sure the current working directory is 'matlab'. Then open and run code/run_comparisons.m.

To plot the results of the comparisons, make sure the current working directory is 'matlab'. Then open and run code/plot_comparisons.m (after run_comparisons completed successfully).

## Fit and predict

# Development status

Development will soon be stalled (end of 2019) as I will no longer have a Matlab license. Help from someone with a license would be greatly appreciated.

