
####################################################################################################
#                                                                                                  #
#   Script for the analysis of Free energy simulations using MBAR                                  #
#                                                                                                  #
#   author: Antonia Mey <antonia.mey@ed.ac.uk>                                                     #
#                                                                                                  #
####################################################################################################

try:
    import numpy as np
except ImportError:
    raise ImportError('Numpy is not installed. Please install numpy in order to use MBAR for your free energy analysis.')
try:
    from pymbar import MBAR
    from pymbar import timeseries
except ImportError:
    raise ImportError('pymbar is not installed. Please install pymbar in order to use MBAR for your free energy '
                      'analysis.')

class FreeEnergies(object):
    r"""This class contains all the different pmf informations"""
    def __init__(self, u_kln, N_k, lambda_array, gradients):
        r"""The data passed here is already subsampled"""
        print ("I am the constructor")
        self._u_kln = u_kln
        self._N_k = N_k
        self._lambda_array = lambda_array
        self._gradients_kn = gradients

        #initialise results containers
        self._deltaF_mbar = None
        self._deltaF_ti = None
        self._f_k = None
        self._pmf_ti = None

    def run_ti(self, cubic_spline=False):
        if cubic_spline:
            NotImplementedError("Cubic Spline TI has not been implemented yet")
        else:
            means = np.mean(self._gradients_kn, axis=1)
            self._pmf_ti = np.zeros(shape=(self._lambda_array.shape[0], 2))
            self._pmf_ti[:, 0] = self._lambda_array
            self._pmf_ti[:, 1] = means
            self._deltaF_ti = np.trapz(means, self._lambda_array)


    def run_mbar(self):
        MBAR_obj = MBAR(self._u_kln, self._N_k, verbose=True)
        self._f_k = MBAR_obj.f_k
        (deltaF_ij, dDeltaF_ij, theta_ij) = MBAR_obj.getFreeEnergyDifferences()
        self._deltaF_mbar = deltaF_ij[0,self._lambda_array.shape[0]-1]
        self._pmf_mbar = np.zeros(shape=(self._lambda_array.shape[0], 2))
        self._pmf_mbar[:, 0] = self._lambda_array
        self._pmf_mbar[:, 1] = self._f_k

    @property
    def pmf_ti(self):
        return self._pmf_ti

    @property
    def pmf_mbar(self):
        return self._pmf_mbar

    @property
    def deltaF_ti(self):
        return self._deltaF_ti

    @property
    def deltaF_mbar(self):
        return self._deltaF_mbar

class SubSample(object):
    r"""This class subsamples data based on the timeseries analysis or percentage of data ready for pmf use"""
    def __init__(self, gradients, energies, u_kln, N_k, percentage=100, subsample=timeseries):
        self._gradients_kn = gradients
        self._N_k = N_k
        self._energies_kn = energies
        self._u_kln = u_kln
        self._subsampled_u_kln = None
        self._subsampled_N_k_energies = None
        self._subsampled_N_k_gradients = None
        self._subsampled_grad_kn = None

        if N_k.shape[0]!=u_kln.shape[0]:
            RuntimeError("The number of thermodynamic states must be the same in u_kln and N_k!"
                         "u_kln has size %d and N_k has size %d" %(u_kln.shape[0], N_k.shape[0]))
        self.subsample_method = subsample
        self.percentage = percentage
        if percentage <0.0 or percentage>100.0:
            RuntimeError("You must provide a percentage between 0 and 100")


    def subsample_gradients(self):
        if self.subsample_method!=timeseries:
            print("We are only eliminating samples from the beginning of the data and are still working with highly"
                  " correlated data!")
            percentage_removal = self._N_k*(1-self.percentage/100.0)
            self._subsampled_grad_kn[k,:] = self._gradients_kn[:,percentage_removal:]
            self._subsampled_N_k_gradients = self._N_k-percentage_removal
        else:
            print("We are doing a timeseries analysis using the timeseries analysis module in pymbar and will subsample"
                  " according to that.")
            #first we compute statistical inefficiency
            g_k = np.zeros(shape=(self._gradients_kn.shape[0]))
            for i in xrange(g_k.shape[0]):
                g_k[i] = timeseries.statisticalInefficiency(self._gradients_kn[i,:])
            g = np.max(g_k)
            #now we need to figure out what the indices in the data are for subsampling
            indices_k = []
            for i in xrange(g_k.shape[0]):
                indices_k.append(timeseries.subsampleCorrelatedData(self._gradients_kn[i,:], g=g))
            self._subsampled_N_k_gradients = (np.ceil(self._N_k / g)).astype(int)
            N_max = np.max(self._subsampled_N_k_gradients)
            self._subsampled_grad_kn = np.zeros([self._gradients_kn.shape[0], N_max], np.float64)
            for k in range(self._gradients_kn.shape[0]):
                self._subsampled_grad_kn[k, :] = self._gradients_kn[k, indices_k[k]]

    def subsample_energies(self):
        if self.subsample_method!=timeseries:
            print("We are only eliminating samples from the beginning of the data and are still working with highly"
                  " correlated data!")
            percentage_removal = self._N_k*(1-self.percentage/100.0)
            self._subsampled_u_kln[k,:] = self._u_kln[:,:,percentage_removal:]
            self._subsampled_N_k_energies = self._N_k-percentage_removal
        else:
            print("We are doing a timeseries analysis using the timeseries analysis module in pymbar and will subsample"
                  " according to that.")
            #first we compute statistical inefficiency
            g_k = np.zeros(shape=(self._energies_kn.shape[0]))
            for i in xrange(g_k.shape[0]):
                g_k[i] = timeseries.statisticalInefficiency(self._energies_kn[i,:])
            g = np.max(g_k)
            print (g)
            #now we need to figure out what the indices in the data are for subsampling
            indices_k = []
            for i in xrange(g_k.shape[0]):
                indices_k.append(timeseries.subsampleCorrelatedData(self._energies_kn[i,:], g=g))
            self._subsampled_N_k_energies = (np.ceil(self._N_k / g)).astype(int)
            N_max = np.max(self._subsampled_N_k_energies)
            self._subsampled_u_kln = np.zeros([self._gradients_kn.shape[0],self._gradients_kn.shape[0], N_max], np.float64)
            for k in range(self._gradients_kn.shape[0]):
                self._subsampled_u_kln[k,:,:] = self._u_kln[k,:,indices_k[k]].transpose()

    @property
    def u_kln(self):
        return self._subsampled_u_kln
    @property
    def gradients_kn(self):
        return  self._subsampled_grad_kn
    @property
    def N_k_energies(self):
        return self._subsampled_N_k_energies
    @property
    def N_k_gradients(self):
        return self._subsampled_N_k_gradients