from abc import ABC, abstractmethod
import numpy as np
from scipy.stats import ortho_group  
from pymanopt import Problem
from pymanopt.manifolds import Product, Stiefel
import tensorflow as tf
from pymanopt.solvers import ConjugateGradient
from sklearn.pipeline import Pipeline


class ManifoldDiscrimantAnalysis(ABC, Pipeline):
    def __init__(self, classes=None, opts=None, usestoppingcrit=True, maxits=1000,
                 store={'Rw': None, 'Rb': None, 'QtRb': None, 'QtRw':  None, 'QtBQ': None, 'QtWQ': None, 'QtWQinvQtBQ': None},
                 optmeth='ManOpt'):
        self.classes = classes
        self.opts = opts
        self.usestoppingcrit = usestoppingcrit
        self.maxits = maxits
        self.store = store
        self.optmeth = optmeth
        if optmeth not in ['bo13', 'wen12', 'ManOpt']:
            self.optmeth = 'ManOpt'
            print("Chosen optimization method", optmeth, "has not been implemented")
        self.Fdifftol = Fdifftol = 10**(-10)
        self.Udifftol = Udifftol = 10**(-12)
        self.MyCost = None
        self.rotations = None


    def set_tolerances(self, Fdifftol=10**(-10), Udifftol=10**(-12)):
        self.Fdifftol = Fdifftol
        self.Udifftol = Udifftol
        
    @abstractmethod 
    def object_matrix_data(self, classmeandiffs, observationdiffs, nis, K1, K2, ): #This is a virtual function that any inheriting subclass should realize
        return(None)


    @abstractmethod
    def my_cost(self): #, x, store, classmeandiffs, observationdiffs, nis, K1, K2):
        pass
        """
        self.F = -np.linalg.trace(self.store['QtWQinvQtBQ'])
        Q = tf.Variable(tf.placeholder(tf.float32))
        """

        #problem = Problem(manifold = manifold, cost = MyCost, arg = Q)

        """
        With tensorflow:
            Q = tf.Variable(tf.placeholder(tf.float32))
            self.MyCost = tf.linalg.trace(g(Q)) Operations to construct my cost inserted herein - potential to save a lot of calculations and definitions in the code when porting?
                    
            problem = Problem(manifold = manifold, cost = MyCost, arg = Q)
            
        """
    @staticmethod
    def class_based_differences(Xs, classes):
        """
        Calculate class based differences.

        The purpose of this method is, given a labelled input data set to output the means of each class, and the data set recentered around the class means for each point.

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

        for i in range(nclasses):
            locations = np.nonzero(classes == i)
            nis[i] = max(np.shape(locations))
            Xsumsclasses[i] = np.sum(Xs[locations], axis=0)
            Xmeansclasses = Xsumsclasses[i]/nis[i]

        for i in range(nsamples): #This should be vectorized
            xi_m_cmeans = Xs[i] - Xmeansclasses[classes[i]]

        cmeans_m_xmeans = Xmeansclasses-Xmean
        return cmeans_m_xmeans, xi_m_cmeans, nis

    def optimize_on_manifold(self, options):
        if self.optmeth == 'ManOpt':
            manifold_one = Stiefel(np.shape(self.Us[0])[0], np.shape(self.Us[0])[1]) #This is hardcoding it to the two-dimensional case..
            manifold_two = Stiefel(np.shape(self.Us[0])[0], np.shape(self.Us[0])[1])
            manifold = Product(manifold_one, manifold_two)
            Q = tf.Variable(tf.placeholder(tf.float32))
            problem = Problem(manifold = manifold, cost = self.my_cost(self), arg = Q) #This assumes TensorFlow implementation... Else we have to implement the gradient and Hessian manually...
            solver = ConjugateGradient(problem, Q, options)
            return solver

    def fit(self, Xs, classes=None, rotations=None):
        lowerdims = np.shape(Xs[0])

        K1 = lowerdims[0] #This hardcodes the two-dimensional case.
        K2 = lowerdims[1]
        #if Us is None:
        #    for i in range(self.nmodes):
        #        Us[i] = ortho_group.rvs(dim = (np.shape(self.sizeX)[i])) #Contemplate changing this

        cmeans_m_xmeans, xi_m_cmeans, nis = self.ClassBasedDifferences(self, Xs, classes)

        self.object_matrix_data(self, cmeans_m_xmeans, xi_m_cmeans, nis, K1, K2)
        self.my_cost(self)
        self.optimize_on_manifold(self)
        self.transform(self, Xs)

    def _transform(self, Xs):
        rotation1 = tf.variable(tf.placeholder(tf.float32))
        rotation2 = tf.variable(tf.placeholder(tf.float32))
        input_data = tf.variable(tf.placeholder(tf.float32))
        product_rotation = tf.tensordot(tf.transpose(rotation1), tf.transpose(rotation2), axes=0)
        data_transformer = tf.tensordot(input_data, product_rotation, axes=[0, 1]) #Double check that this is correct
        with tf.Session as sess:
            transformed_data = sess.run(data_transformer,
                                        feed_dict={rotation1: self.rotations[0], rotation2: self.rotations[1], input_data: Xs})

        return transformed_data


class TuckerDiscriminantAnalysis(ManifoldDiscrimantAnalysis):
    def QtCheck(self, Qt):
        if self.store[Qt] is None:
            self.store[Qt] = self.QtCalculator(Qt)

    def Qt_initializer(self, Qt, K1, K2, N, M):

        rotation1 = tf.Variable(tf.placeholder(tf.float32))
        rotation2 = tf.Variable(tf.placeholder(tf.float32))
        if Qt in {'QtRw', 'QtRb'}:
            Qt_mm = tf.tensordot(tf.tensordot(self.Qt[-2:], tf.linalg.transpose(rotation1), axes=(1, 0)), tf.linalg.transpose(rotation2), axes = (2, 0)) #Something might be horribly wrong here...
            Qt_temp = tf.reshape(tf.roll(Qt_mm, (1, 0, 2))(K2*K1, N)) #Not sure if this is the correct approach. Need MATLAB license to test original behavior
        if Qt in {'QtWQ', 'QtBQ'}:
            Qt_temp = tf.matmul('QtR'+Qt[2].lower(), self.store['QtR'+Qt[2].lower()]) # #This is not the best way to do this.. But it's not the worst?
        if Qt in {'QtWQinvQtBQ'}:
            Qt_temp = tf.matmul(self.store['QtBQ'], tf.linalg.inv(self.store['QtWQ']))
        else:
            print("Wrong Qt specified")
        self.store[Qt] = Qt_temp


    def object_matrix_data(self, classmeandiffs, observationdiffs, nis, K2, K1):
        if self.store['Rw'] is None or self.store['Rb'] is None:
            obsExample = classmeandiffs[0]
            sizeObs = np.shape(obsExample)
            I = sizeObs[0]
            J = sizeObs[1]
            if self.store['Rw'] is None: #What if only one of them is missing? Shouldn't we still calculate the missing one? Split the code below into two parts, one for each case. This has been done now, but not in the MATLAB code
                nclasses = max(np.shape(classmeandiffs)) #This is a straight-copy of the existing MATLAB code, not sure if this is the intended behavior there either...
                classmeandiffstensor = np.reshape(classmeandiffs, (I, J, nclasses), order = "F")
                Rb = classmeandiffstensor
            if self.store['Rb'] is None:
                nobs = max(np.shape(observationdiffs))
                observationdiffstensor = np.reshape(observationdiffs, (I, J, nobs), order = "F")
                Rw = observationdiffstensor
            self.store['Rw'] = Rw
            self.store['Rb'] = Rb
#        Rw = store['Rw'] #This assignment SHOULD be superflous now...
#        Rw = store['Rb']
        Rwsize = np.shape(Rw)
        nobs = Rwsize[-1]
        datadims = np.shape(Rb)
        nclasses = len(datadims)
        N = datadims[0] #This is a point where the two-dimensionality is hardcoded..
        M = datadims[1]
        #We proceed to calculate all relevant matrices for the optimization step.
        for j in {'QtRw', 'QtRb', 'QtWQ', 'QtBQ', 'QtWQinvQtBQ'}:
            self.Qt_initializer(j, K2, K1, N,  M)

    def my_cost(self):  # , x, store, classmeandiffs, observationdiffs, nis, K1, K2):
        # self.object_matrix_data(self, x, store, classmeandiffs, observationdiffs, nis, K1, K2)

        self.MyCost = tf.linalg.trace(self.store['QtWQinvQtBQ'])


class PARAFACDiscriminantAnalysis(ManifoldDiscrimantAnalysis):
    def object_matrix_data(U, cmean_m_xmeans, xi_m_cmeans, nis, K1, K2):
        
        
        return()