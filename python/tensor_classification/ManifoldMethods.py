from abc import ABC, abstractmethod
import numpy as np
from pymanopt import Problem
from pymanopt.manifolds import Product, Stiefel
import tensorflow as tf
from pymanopt.solvers import ConjugateGradient


class ManifoldDiscrimantAnalysis(ABC):
    """
    This is an abstract class serving to implement joint methods between PARAFAC and Tucker discriminant analysis.
    The class implements optimization, fitting and calculating class-based differences.
    Methods which have to be defined in either case but are not shared are declared as abstract.
    """
    def __init__(self, Rw=None, Rb=None):
        """
        The initialization method for the PARAFAC and Tucker models, creating an empty cost function, empty list of
            optimal rotations and a dictionary for saving Tensorflow functions defining 'Qt*'.

        :param Rw:  Matrix specifying pre-calculated'Rw' to avoid calculations.
        :param Rb:  Matrix specifying pre-calculated'Rw' to avoid calculations.
        """
        self.store ={'Rw': Rw, 'Rb': Rb, 'QtRb': None,
                     'QtRw':  None, 'QtBQ': None, 'QtWQ': None, 'QtWQinvQtBQ': None}
        self.f_diff_tol = 10 ** (-10)
        self.u_diff_tol = 10 ** (-12)
        self.MyCost = None
        self.rotations = None

    def set_tolerances(self, f_diff_tol=10 ** (-10), u_diff_tol=10 ** (-12)):
        self.f_diff_tol = f_diff_tol
        self.u_diff_tol = u_diff_tol

    # This is a virtual function that any inheriting subclass should realize
    @abstractmethod 
    def object_matrix_data(self, class_mean_diffs, observation_diffs, nis, k1, k2, ):
        return None

    @abstractmethod
    def my_cost(self):
        pass

    @staticmethod
    def class_based_differences(Xs, classes):
        """
        Calculate class based differences.

        The purpose of this method is, given a labelled input data set to output the means of each class,
        and the data set recentered around the class means for each point.

        :param Xs: numpy tensor containing the input data, indexed along the first mode.
        :param classes: The class for each input data point.

        :return nis : Numpy-array list of samples in each class
        :return cmeans_m_xmeans : Numpy-array containing means for each class minus the overall mean
        :return xi_m_cmeans : Numpy-array of the samples centered around the class means
        """
        classes = np.array(classes)
        Xs = np.array(Xs)
        nsamples = max(np.shape(Xs))
        nclasses = len(np.unique(classes))
        Xsum = np.sum(Xs, axis=0)
        Xmean = Xsum/nsamples
        shape = np.shape(Xs[0])

        Xsumsclasses = np.zeros([nclasses, *shape])
        Xmeansclasses = np.zeros([nclasses, *shape])
        nis = np.zeros(nclasses)
        xi_m_cmeans=np.zeros([nsamples]+list(shape))

        for i in range(nclasses):
            locations = np.nonzero(classes == i)
            nis[i] = max(np.shape(locations))
            Xsumsclasses[i] = np.sum(Xs[locations], axis=0)
            Xmeansclasses = Xsumsclasses[i]/nis[i]

        for i in range(nsamples):  # This should be vectorized
            xi_m_cmeans[i] = Xs[i] - Xmeansclasses[classes[i]]

        cmeans_m_xmeans = Xmeansclasses-Xmean

        return cmeans_m_xmeans, xi_m_cmeans, nis

    def optimize_on_manifold(self, options, optmeth):
        if optmeth not in ['bo13', 'wen12', 'ManOpt']:
            print("Chosen optimization method", optmeth, "has not been implemented, using 'ManOpt' ")
            optmeth = 'ManOpt'

        if optmeth == 'ManOpt':
            # This is hardcoding it to the two-dimensional case..
            manifold_one = Stiefel(np.shape(self.rotations[0])[0], np.shape(self.rotations[0])[1])
            manifold_two = Stiefel(np.shape(self.rotations[0])[0], np.shape(self.rotations[0])[1])
            manifold = Product((manifold_one, manifold_two))
            optimization_variable = tf.Variable(tf.placeholder(tf.float32))
            problem = Problem(manifold=manifold, cost=self.my_cost(), arg=optimization_variable)
            solver = ConjugateGradient(problem, optimization_variable, options)

            return solver

    def fit(self, Xs, classes, optmeth, options=None):
        lowerdims = np.shape(Xs[0])

        K1 = lowerdims[0] #This hardcodes the two-dimensional case.
        K2 = lowerdims[1]

        cmeans_m_xmeans, xi_m_cmeans, nis = self.class_based_differences(Xs, classes)

        self.object_matrix_data(cmeans_m_xmeans, xi_m_cmeans, nis, K1, K2)
        self.my_cost()
        self.optimize_on_manifold(options, optmeth)
        fitted_data = self.transform(Xs)
        return fitted_data

    def transform(self, Xs):
        """
        This method transforms the observed data with the optimized rotations from the manifold optimization step.

        The data is transformed by acting with the tensor product of the rotations on each observation.
           This method assumes two-dimensional observations as currently written.
           The ´transform´ method overrides the Pipeline transform method,
           and is a prerequisite for fitting into a pipeline.


        :params: Xs: Tensor containing the input data to be transformed.
        :return: transformed_data, Tensor containing the data transformed by the optimal transformations.
        """
        rotation1 = tf.Variable(tf.placeholder(tf.float32))
        rotation2 = tf.Variable(tf.placeholder(tf.float32))
        input_data = tf.Variable(tf.placeholder(tf.float32))
        product_rotation = tf.tensordot(tf.transpose(rotation1), tf.transpose(rotation2), axes=0)
        data_transformer = tf.tensordot(input_data, product_rotation, axes=[1, 2]) #Double check that this is correct
        with tf.Session as sess:
            transformed_data = sess.run(data_transformer,
                                        feed_dict={rotation1: self.rotations[0], rotation2: self.rotations[1],
                                                   input_data: Xs})

        return transformed_data


class TuckerDiscriminantAnalysis(ManifoldDiscrimantAnalysis):

    def Qt_initializer(self, Qt, K1, K2, N, M):

        rotation1 = tf.Variable(tf.placeholder(tf.float32))
        rotation2 = tf.Variable(tf.placeholder(tf.float32))
        if Qt in {'QtRw', 'QtRb'}:
            # Something might be horribly wrong here...
            if Qt == 'QtRw':
                dims = N
            else:
                dims = M
            Qt_mm = tf.tensordot(tf.tensordot(Qt[-2:], tf.linalg.transpose(rotation1), axes=(1, 0)),
                                 tf.linalg.transpose(rotation2), axes=(2, 0))
            # Not sure if this is the correct approach. Need MATLAB license to test original behavior
            Qt_temp=tf.reshape(tf.roll(Qt_mm, (1, 0, 2))(K2*K1, dims))
        if Qt in {'QtWQ', 'QtBQ'}:
            # #This is not the best way to do this.. But it's not the worst?
            Qt_temp=tf.matmul('QtR'+Qt[2].lower(), self.store['QtR'+Qt[2].lower()])
        if Qt in {'QtWQinvQtBQ'}:
            Qt_temp=tf.matmul(self.store['QtBQ'], tf.linalg.inv(self.store['QtWQ']))
        else:
            print("Wrong Qt specified")

        self.store[Qt] = Qt_temp

    def object_matrix_data(self, class_mean_diffs, observation_diffs, nis, k2, k1):
        if self.store['Rw'] is None or self.store['Rb'] is None:
            obsExample = class_mean_diffs[0]
            sizeObs = np.shape(obsExample)
            I = sizeObs[0]
            J = sizeObs[1]
            if self.store['Rw'] is None: #What if only one of them is missing? Shouldn't we still calculate the missing one? Split the tensor_classification below into two parts, one for each case. This has been done now, but not in the MATLAB tensor_classification
                nclasses = max(np.shape(class_mean_diffs)) #This is a straight-copy of the existing MATLAB tensor_classification, not sure if this is the intended behavior there either...
                classmeandiffstensor = np.reshape(class_mean_diffs, (I, J, nclasses), order="F")
                self.store['Rw']= classmeandiffstensor

            if self.store['Rb'] is None:
                nobs = max(np.shape(observation_diffs))
                observationdiffstensor = np.reshape(observation_diffs, (I, J, nobs), order="F")
                self.store['Rb'] = observationdiffstensor

        Rwsize = np.shape(self.store['Rw'])
        #nobs = Rwsize[0]
        datadims = np.shape(self.store['Rb'])
        #nclasses = len(datadims)
        N = datadims[0] #This is a point where the two-dimensionality is hardcoded..
        M = datadims[1]
        #We proceed to calculate all relevant matrices for the optimization step.
        for j in {'QtRw', 'QtRb', 'QtWQ', 'QtBQ', 'QtWQinvQtBQ'}:
            self.Qt_initializer(j, k2, k1, N, M)

    def my_cost(self):  # , x, store, classmeandiffs, observationdiffs, nis, K1, K2):
        self.MyCost = tf.linalg.trace(self.store['QtWQinvQtBQ'])


class PARAFACDiscriminantAnalysis(ManifoldDiscrimantAnalysis):
    def my_cost(self):
        pass

    def object_matrix_data(U, cmean_m_xmeans, xi_m_cmeans, nis, k1, k2):
        pass
