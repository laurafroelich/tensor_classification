import pytest
from classificationmethods import ManifoldDiscrimantAnalysisParafacTucker as pt
import numpy as np
from numpy import random as random
import scipy.stats as stats

"""
    Unit tests for python manifold discriminant analysis. 
"""


def generate_data(classes, length_of_data_set, shape=None):
    end_padding = length_of_data_set % classes
    number_of_copies = length_of_data_set // classes
    output_classes = number_of_copies * list(range(classes)) + list(range(end_padding))
    if shape is None:
        shape = (5, 6)

    observations = np.multiply.outer(output_classes, np.zeros(shape))

    for i in range(classes):
        observations[i] = output_classes[i] + np.random.uniform(-0.5, 0.5, shape)

    return observations, output_classes

def test_QtCheck():
    assert False


def test_based_differences():
    assert False
    # model = pt.TuckerDiscriminantAnalysis()
    # input_params = generate_data(5, 1000)
    # input_data = input_params[0]
    # input_classes = input_params[1]
    # means = model.class_based_differences(input_data, input_classes)

def test_fail():
    assert False


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


def test_class_based_differences():
    modeller = pt.TuckerDiscriminantAnalysis()
    pass


def test_tucker_object_data_matrix():
    pass


def test_my_cost():
    pass
