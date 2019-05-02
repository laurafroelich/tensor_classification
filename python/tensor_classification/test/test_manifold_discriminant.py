import pytest
from tensor_classification.classificationmethods import ManifoldDiscrimantAnalysisParafacTucker as pt
import numpy as np
from numpy import random as random
import scipy.stats as stats

"""
    Unit tests for python manifold discriminant analysis. 
"""


def generate_data(classes, observations_per_class, shape=None):
    """
    The purpose of this method is to simulate data for testing the manifold discriminant analysis tensor_classification.

    For each class i we simulate the same number of observations where the data is uniformly distributed around i.
       The method returns the simulated observations and the class for each observation.

    :param classes: Int, integer specifying number of classes to simulate observations for.
    :param observations_per_class: Int, integer specifying the number of observations per class.
    :param shape: List, list specifying the shape of each simulated observation.
    :return: observations: Tensor containing all simulated observations
    :return: output_classes: List of classes specifying the class for each simulated observation
    """
    if shape is None:
        shape = [5, 6]

    output_classes = np.zeros(classes*observations_per_class)
    observations = np.multiply.outer(output_classes, np.zeros(shape))

    for i in range(classes):
        observations[observations_per_class*i:observations_per_class*(i+1)] \
            = i + np.random.uniform(-0.5, 0.5, [observations_per_class]+shape)

        output_classes[observations_per_class*i:observations_per_class*(i+1)] = i

    return observations, output_classes

#def test_QtCheck():
#    assert False


def test_based_differences():
    model = pt.TuckerDiscriminantAnalysis()
    input_params = generate_data(5, 1000)
    input_data = input_params[0]
    input_classes = input_params[1]
    means = model.class_based_differences(input_data, input_classes)
    print(means)
    assert True


#def test_fail():
#    assert False


def test_generate_data():
    # generate_data(5, 31)
    pass


def test_fit_manifold():
    assert True
    # pass


def test_fit_pipeline():
    assert True
    # pass


def test_initialize_manifold_with_u():
    pass


def test_initialize_manifold_without_u():
    pass


def test_calculate_set_tolerances():
    modeller = pt.TuckerDiscriminantAnalysis()
    modeller.set_tolerances(0, 0)
    assert modeller.Fdifftol == 0, modeller.Udifftol == 0

def test_tucker_object_data_matrix():
    pass


def test_my_cost():
    pass
