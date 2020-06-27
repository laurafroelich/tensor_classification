import numpy as np


def differences(Xs, classes):
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
    Xmean = Xsum / nsamples
    shape = np.shape(Xs[0])

    Xsumsclasses = np.zeros([nclasses, *shape])
    Xmeansclasses = np.zeros([nclasses, *shape])
    nis = np.zeros(nclasses)
    xi_m_cmeans = np.zeros([nsamples] + list(shape))

    for i in range(nclasses):
        locations = np.nonzero(classes == i)
        nis[i] = max(np.shape(locations))
        Xsumsclasses[i] = np.sum(Xs[locations], axis=0)
        Xmeansclasses = Xsumsclasses[i] / nis[i]

    for i in range(nsamples):  # This should be vectorized
        xi_m_cmeans[i] = Xs[i] - Xmeansclasses[classes[i]]

    cmeans_m_xmeans = Xmeansclasses - Xmean

    return cmeans_m_xmeans, xi_m_cmeans, nis