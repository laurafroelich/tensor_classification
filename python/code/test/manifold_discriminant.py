import pytest


class MyTestCase(unittest.TestCase):
    def test_something(self):
        self.assertEqual(True, False)


if __name__ == '__main__':
    unittest.main()

"""
    Unit tests for python manifold discriminant analysis. 
"""


def test_fit_manifold():
    pass


def test_fit_pipeline():
    pass


def test_initialize_manifold_with_u():
    pass


def test_initialize_manifold_without_u():
    pass


def test_calculate_set_tolerances():
    modeller=ManifoldDiscrimantAnalysis(np.zeros((2, 1)))
    modeller.SetTolerances(0, 0)
    assert modeller.FFdifftol==0, modeller.Udifftol==0



def test_class_based_differences():
    pass

def test_tucker_object_data_matrix():
    pass


def test_my_cost():
    pass