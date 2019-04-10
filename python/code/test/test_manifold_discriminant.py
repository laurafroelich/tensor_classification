import pytest
from classificationmethods import ManifoldDiscrimantAnalysisParafacTucker as pt
import numpy as np
"""
    Unit tests for python manifold discriminant analysis. 
"""


def test_fit_manifold():
    assert True
    #pass


def test_fit_pipeline():
    assert True
    #pass


def test_initialize_manifold_with_u():
    pass


def test_initialize_manifold_without_u():
    pass

def test_calculate_set_tolerances():
    modeller = pt.TuckerDiscriminantAnalysis()
    modeller.set_tolerances(0, 0)
    assert modeller.Fdifftol == 0, modeller.Udifftol == 0



def test_class_based_differences():
    modeller=pt.TuckerDiscriminantAnalysis()
    pass

def test_tucker_object_data_matrix():
    pass


def test_my_cost():
    pass
