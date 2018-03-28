
####################################################################################################
#                                                                                                  #
#   Script for the analysis of Free energy simulations using MBAR                                  #
#                                                                                                  #
#   author: Antonia Mey <antonia.mey@ed.ac.uk>                                                     #
#                                                                                                  #
####################################################################################################
from Sire.Analysis import *
import Sire.Stream
import warnings
from Sire.Units import *

try:
    numpy = Sire.try_import("numpy")
except ImportError:
    raise ImportError('Numpy is not installed. Please install numpy in order to use MBAR for your free energy analysis.')
try:
    scipy = Sire.try_import("scipy")
except ImportError:
    raise ImportError('Numpy is not installed. Please install numpy in order to use MBAR for your free energy analysis.')
from pymbar import MBAR
from pymbar import timeseries
import matplotlib.pylab as plt
from ipywidgets import interact, interactive, fixed, interact_manual, Layout, Label
import ipywidgets as widgets
import bz2
import glob
import os
import sys
import warnings


class FreeEnergies(object):
    r"""This class contains all the different pmf information
    The constructor expects subsampled MBAR and TI compatible data.
    Parameters
    ----------

    u_kln : ndarray(shape=(therm_states, therm_states, nsamples), dtype=float)
        reduced perturbed energies used for MBAR estimates
    N_K : ndarray(shape=(therm_states), dtype=int)
        number of samples per thermodynamic state
    lambda_array : ndarray(shape=(therm_states), dtype=float)
        lambda thermodynamic values
    gradients_kn : ndarray(shape=(therm_state, nsamples), dtype=float)
        reduced gradients
    """

    def __init__(self, u_kln=None, N_k=None, lambda_array=None, gradients_kn=None):
        r"""The data passed here is already subsampled"""

        self._u_kln = u_kln
        self._N_k = N_k
        self._lambda_array = lambda_array
        self._gradients_kn = gradients_kn

        #initialise results containers
        self._deltaF_mbar = None
        self._deltaF_ti = None
        self._dDeltaF_mbar = None
        self._f_k = None
        self._pmf_ti = None
        self._overlap_matrix = None
        self._pairwise_F = None

    def run_ti(self, cubic_spline=False):
        r"""Runs Thermodynamic integration free energy estimate
        Parameters
        ----------

        cubic_spline : bool
            Use cubic spline estimation instead of trapezium rule.
        """
        means = numpy.mean(self._gradients_kn, axis=1)
        if cubic_spline:
            NotImplementedError("Cubic Spline TI has not been implemented yet")
        else:
            self._pmf_ti = numpy.zeros(shape=(self._lambda_array.shape[0], 2))
            self._pmf_ti[:, 0] = self._lambda_array
            for i in range(1, self._lambda_array.shape[0]):
                self._pmf_ti[i-1][1] = numpy.trapz(means[0:i], self._lambda_array[0:i])
            self._pmf_ti[-1][1] = numpy.trapz(means, self._lambda_array)
            self._deltaF_ti = numpy.trapz(means, self._lambda_array)


    def run_mbar(self, test_overlap = True):
        r"""Runs MBAR free energy estimate """
        MBAR_obj = MBAR(self._u_kln, self._N_k, verbose=True)
        self._f_k = MBAR_obj.f_k
        (deltaF_ij, dDeltaF_ij, theta_ij) = MBAR_obj.getFreeEnergyDifferences()
        self._deltaF_mbar = deltaF_ij[0, self._lambda_array.shape[0]-1]
        self._dDeltaF_mbar = dDeltaF_ij[0, self._lambda_array.shape[0]-1]
        self._pmf_mbar = numpy.zeros(shape=(self._lambda_array.shape[0], 3))
        self._pmf_mbar[:, 0] = self._lambda_array
        self._pmf_mbar[:, 1] = self._f_k
        self._pmf_mbar[:,2] = dDeltaF_ij[0]
        self._pairwise_F = numpy.zeros(shape=(self._lambda_array.shape[0]-1,4))
        self._pairwise_F[:,0] = self._lambda_array[:-1]
        self._pairwise_F[:,1] = self._lambda_array[1:]
        self._pairwise_F[:,2] = numpy.diag(deltaF_ij,1)
        self._pairwise_F[:,3] = numpy.diag(dDeltaF_ij,1)        


        ##testing data overlap:
        if test_overlap:
            overlap_matrix = MBAR_obj.computeOverlap()
            self._overlap_matrix = overlap_matrix[2]


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

    @property
    def errorF_mbar(self):
        return self._dDeltaF_mbar

    @property
    def overlap_matrix(self):
        return self._overlap_matrix

    @property
    def pairwise_F(self):
        return self._pairwise_F

class SubSample(object):
    r"""This class subsamples data based on the timeseries analysis or percentage of data ready for pmf use
    Parameters
    ----------
    gradients_kn : ndarray(shape=(therm_state, nsamples), dtype=float)
        reduced gradients
    energies : ndarray(shape=(therm_state, nsamples), trype=float)
        potential energies used to find statisitical inefficiency
    u_kln : ndarray(shape=(therm_states, therm_states, nsamples), dtype=float)
        reduced perturbed energies used for MBAR estimates
    N_K : ndarray(shape=(therm_states), dtype=int)
        number of samples per thermodynamic state
    lambda_array : ndarray(shape=(therm_states), dtype=float)
        lambda thermodynamic values
    percentage : int [0,100]
        percentage of the data that should be discarded from the beginning of the simulation
    subsample : string
        string idenfier for subsampling method, default='timeseries' from timeseries module in MBAR
    """

    def __init__(self, gradients_kn, energies, u_kln, N_k, percentage=100, subsample=True):
        self._gradients_kn = gradients_kn
        self._N_k = N_k
        self._energies_kn = energies
        self._u_kln = u_kln
        self._subsampled_u_kln = None
        self._subsampled_N_k_energies = None
        self._subsampled_N_k_gradients = None
        self._subsampled_grad_kn = None
        self._subsampled_energies_kn = None

        if N_k.shape[0]!=u_kln.shape[0]:
            RuntimeError("The number of thermodynamic states must be the same in u_kln and N_k!"
                         "u_kln has size %d and N_k has size %d" %(u_kln.shape[0], N_k.shape[0]))
        self.subsample = subsample
        self.percentage = percentage
        assert(percentage > 0.0 and percentage <=100.0), "You must provide a percentage between 0 and 100" 


    def subsample_gradients(self):
        r''' method to subsample gradients and get a better estiamte.
        '''
        if self.percentage == 100 and not self.subsample:
            warnings.warn("You are not subsampling your data according to the statistical inefficiency nor are "
                           "you discarding initial data. Please set percentage to another value than 100!")
        percentage_removal = (self._N_k*(1-self.percentage/100.0)).astype('int32')
        self._subsampled_N_k_gradients = self._N_k-percentage_removal
        N_max = int(numpy.max(self._subsampled_N_k_gradients))
        self._subsampled_grad_kn = numpy.zeros(shape=(self._N_k.shape[0], N_max))
        for p in range(percentage_removal.shape[0]):
            start = percentage_removal[p]
            finish = percentage_removal[p]+N_max
            self._subsampled_grad_kn[p,:] = self._gradients_kn[p,start:finish]
        if N_max <=50:
            warnings.warn("You have reduced your data to less than 50 samples, the results from these might not "
                           "be trustworthy. If you don't want to add more samples consider rerunning the analysis using the percentage option.")
        #if subsampling is percentage, then we are done here, otherwise we will now subsample according to timeseries

        if self.subsample:
            print("#Subsampling gradients according to statistical inefficiency")
            #first we compute statistical inefficiency
            self._gradients_kn = self._subsampled_grad_kn.copy()
            self._N_k = self._subsampled_N_k_gradients.copy()

            g_k = numpy.zeros(shape=(self._gradients_kn.shape[0]))
            self._subsampled_N_k_gradients = numpy.zeros(shape=(self._gradients_kn.shape[0]))
            for i in range(g_k.shape[0]):
                g_k[i] = timeseries.statisticalInefficiency(self._gradients_kn[i,:])
            g = int(numpy.max(g_k))
            #now we need to figure out what the indices in the data are for subsampling
            indices_k = []
            for i in range(g_k.shape[0]):
                indices_k.append(timeseries.subsampleCorrelatedData(self._gradients_kn[i,:], g=g))
                self._subsampled_N_k_gradients[i]=len(indices_k[i])
            N_max = int(numpy.max(self._subsampled_N_k_gradients))
            if N_max <=50:
                warnings.warn("You have reduced your data to less than 50 samples, the results from these might not "
                               "be trustworthy. If you don't want to add more samples consider rerunning the analysis using the percentage option.")
            self._subsampled_grad_kn = numpy.zeros([self._gradients_kn.shape[0], N_max], numpy.float64)
            for k in range(self._gradients_kn.shape[0]):
                self._subsampled_grad_kn[k, :] = self._gradients_kn[k, indices_k[k]]

    def subsample_energies(self):
        r''' This subsamples u_kln according to percentage, i.e. remove initial equilibration data and then can additionally subsample according to timeseries

        '''
        #removing percent
        if self.percentage == 100 and not self.subsample:
            warnings.warn("You are not subsampling your data according to the statistical inefficiency nor are "
                           "you discarding initial data. Please set percentage to another value than 100!")

        percentage_removal = (self._N_k*(1-self.percentage/100.0)).astype('int32')
        self._subsampled_N_k_energies = self._N_k-percentage_removal
        N_max = int(numpy.max(self._subsampled_N_k_energies))
        self._subsampled_u_kln = numpy.zeros(shape=(self._N_k.shape[0], self._N_k.shape[0], N_max))
        self._subsampled_energies_kn = numpy.zeros(shape=(self._N_k.shape[0], N_max))
        for k in range(0, self._N_k.shape[0]):
            self._subsampled_u_kln[k] = self._u_kln[k,:,percentage_removal[k]:percentage_removal[k]+N_max]
            self._subsampled_energies_kn[k] = self._energies_kn[k,percentage_removal[k]:percentage_removal[k]+N_max]
        if N_max <=50:
            warnings.warn("You have reduced your data to less than 50 samples, the results from these might not "
                           "be trustworthy. If you don't want to add more samples consider rerunning the analysis using the percentage option.")

        #Now we are doing some additional subsampling according to timeseries analysis
        if self.subsample:
            print("#Subsampling energies according to statistical inefficiency for pymbar")

            self._u_kln = self._subsampled_u_kln.copy()
            self._N_k = self._subsampled_N_k_energies.copy()
            self._energies_kn = self._subsampled_energies_kn.copy()
            #first we compute statistical inefficiency
            g_k = numpy.zeros(shape=(self._energies_kn.shape[0]))
            for i in range(g_k.shape[0]):
                g_k[i] = timeseries.statisticalInefficiency(self._energies_kn[i,percentage_removal[i]:])
            g = numpy.max(g_k)
            #now we need to figure out what the indices in the data are for subsampling
            indices_k = []
            self._subsampled_N_k_energies = numpy.zeros(shape=(self._energies_kn.shape[0]))
            for i in range(g_k.shape[0]):
                indices_k.append(timeseries.subsampleCorrelatedData(self._energies_kn[i,:], g=g))
                self._subsampled_N_k_energies[i]=len(indices_k[i])
            #self._subsampled_N_k_energies = (numpy.ceil(self._N_k / g)).astype(int)
            N_max = int(numpy.max(self._subsampled_N_k_energies))
            if N_max <=50:
                warnings.warn("You have reduced your data to less than 50 samples, the results from these might not "
                               "be trustworthy. If you don't want to add more samples consider rerunning the analysis using the percentage option.")
            self._subsampled_u_kln = numpy.zeros([self._gradients_kn.shape[0],self._gradients_kn.shape[0], N_max], numpy.float64)
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

class NotebookHelper(object):
    r"""

    """
    def __init__(self):
        self._out = None
        self._basedir = None
        self._outputdir = None
        self._free_energy_file = None
        self._discrad = None
        self._overlap = None
        self._subsample = None
        self._generate_notebook = None
        self._perturbation_list = None


    def initialise_notebook(self):
        ''' initialises the checkboxes for the jupyter notebook
        '''
        style = {'description_width': 'initial'}
        layout = Layout(flex='2 1 auto', width='auto')
        a = widgets.Text(placeholder='/path/to/simulation/data', description='Simulation Base Directory:',disabled=False,layout=layout, 
    style = style)
        b = widgets.Text(placeholder='/path/to/output/directory', description='Output Directory:',disabled=False,layout=layout, 
    style = style)
        c = widgets.Text(placeholder='Free_energy_summary.csv', description='Free energy file:',disabled=False,layout=layout, 
    style = style)
        d = widgets.BoundedIntText(value=0, min=0, max=1000000, step=1, description='Number of frames to discard',layout=layout, 
    style = style)
        e = widgets.Checkbox(False, description='Plot Overlap matrices',layout=layout, 
    style = style)
        f = widgets.Checkbox(False, description='Subsample data according to statistical inefficiency',layout=layout, 
    style = style)
        g = widgets.Checkbox(False, description='Generate Network analysis notebook',layout=layout, 
    style = style)
        ui = widgets.VBox([a, b, c, d, e, f, g], layout=Layout(display='flex',
                        flex_flow='column',
                        align_items='stretch',
                        border='solid',
                        width='100%'))
        def _func(a,b,c,d,e,f,g):
            print((a,b,c,d,e,f,g))

        self._out = widgets.interactive_output(_func, {'a': a, 'b': b, 'c': c, 'd': d, 'e': e, 'f': f, 'g': g})
        return ui

    def update(self):
        r"""updates the information passed from the widgets of the notebook cells. 

        """
        inputs = (self._out.get_state()['outputs'][0]['text'])
        bad_chars = '()"\n '
        inputs = ''.join(c for c in inputs if c not in bad_chars).split(',')
        self._basedir = inputs[0].strip("'")
        self._outputdir = inputs[1].strip("'")
        self._free_energy_file = inputs[2].strip("'")
        self._discard = int(inputs[3].strip("'"))
        self._overlap = self._get_boolean_state((inputs[4].strip("'")))
        self._subsample = self._get_boolean_state((inputs[5].strip("'")))
        self._generate_notebook = self._get_boolean_state((inputs[6].strip("'")))

        #Set the output directory to the current working directory if no path has been given
        if self._outputdir == '':
            self._outputdir = os.getcwd()
        else:
            #check if outputdir exists
            if not os.path.isdir(self._outputdir):
                os.makedirs(self._outputdir)

    def compute_free_energies(self, input_files, TI = False):
        r"""computes free energies
        Paratmerters:
        -------------
        input_files : FILES
            list of simulation.dat files for a given lambda
        TI : boolean
            decides whether to also compute TI free energies or just MBAR free energies
            Default: False

        Retruns:
        free_energies : FreeEnergy object
            object contains free energy differences 
        T : Float
            temperature at which simulations were run, as recorded in the simulation.dat files
        """

        T, lamdas = self._read_sim_parameters(input_files)
        data = self._read_data(input_files)
        free_energies = self._run_preprocessing(data, lamdas)
        free_energies.run_mbar()
        if TI:
            free_energies.run_ti()
        return free_energies, T

    def run_free_energy_analysis(self, perturbation, files_input_bound, files_input_free, TI= False, separator = '~'):
        r""" runs free energy analysis for ligand bound to protein and ligand in water
        Parameters
        ----------
        perturbation : string
            directory in which to find the perturbation. Usually it would take the form:
            compound1~compoound2. The separator '~' can be set separately. 
        files_input_bound : list
            list of path to files containing simfile.dat for a particular perturbation of the simulations where
            the ligand is bound to the protein during the alchemical perturbation.
        files_input_free : list
            list of path to files containing simfile.dat for a particular perturbation of the simulations, where
            the ligand is solvated in a water box during the alchemical perturbation.
        TI : boolean
            Default = False
            Indicates wheter to also evaluate Free energy differences using thermodynamic integration
        separator : string
            Default = '~'
            Separator used to indicate the perturbation, e.g. compound1-compound2, or compound1%compound2. The use of
            '~' as a separator is recommended. 

        Retruns:
        --------
        result : string
            preformatted string used to output free energies in the format: compound1,compound2,DDG,dDDG
            This preformatted output will support a further network analysis. 
        """
        #bound
        T, lamdas = self._read_sim_parameters(files_input_bound)
        data = self._read_data(files_input_bound)
        bound = self._run_preprocessing(data, lamdas)
        bound.run_mbar(self._overlap)
        if TI:
            bound.run_ti()
        if self._overlap:
            self._write_overlap_info(perturbation, bound.overlap_matrix, sim_type='bound')
        
        #free
        T_free, lamdas = self._read_sim_parameters(files_input_free)
        data = self._read_data(files_input_free)
        free = self._run_preprocessing(data, lamdas)
        free.run_mbar(self._overlap)
        if TI:
            free.run_ti()
        if self._overlap:
            self._write_overlap_info(perturbation, free.overlap_matrix, sim_type='free')
        
        #gathering results
        if TI:
            DDG = (bound.deltaF_ti* T * k_boltz)-(free.deltaF_ti* T_free * k_boltz)
            return ('%s,%s,%.2f' %(perturbation.split(separator)[0],perturbation.split(separator)[1],DDG))
        # mbar results
        else:
            DDG = (bound.deltaF_mbar* T * k_boltz)-(free.deltaF_mbar* T_free * k_boltz)
            dDDG = numpy.sqrt((bound.errorF_mbar * T * k_boltz)**2+(free.errorF_mbar * T_free * k_boltz)**2) 
            return ('%s,%s,%.2f,%.2f' %(perturbation.split(separator)[0],perturbation.split(separator)[1],DDG,dDDG))

    def write_free_energies(self, DDG_list):
        r"""
        writes energies either to file or std out
        Parameters:
        DDG_list : list of strings
            output from run__free_energy_analysis, preformatted free energy string
        """
        if self._free_energy_file is None or self._free_energy_file is '':
            fh = sys.stdout
        else:
            fname = os.path.join(self._outputdir,self._free_energy_file)
            fh = open(fname, 'w')
        for line in DDG_list:
            fh.write(line+'\n')
        fh.close()

    def _read_sim_parameters(self, input_files):
        r""" reading input files from a given file list of simfile.dat files
        Parameters:
        ----------
        input_files : list
            list containing all simfile.dat files to read from

        Returns:
        --------
        T : float
            simulation temperature as extracted from the simfile.dat
        lambda : numpy array
            array containing lambda values at which simulation was carried out

        """
        lamvals = None
        T_previous = None
        for f in input_files:
            #print ('working on input file %s' % f)
            # Compressed file
            if f.endswith('.bz2'):
                bz_file = bz2.BZ2File(f)
                for line in bz_file:
                    if line.startswith(b'#Alchemical'):
                        if lamvals is None:
                            lamvals = (line.split(b'(')[-1].split(b')')[0].split(b','))
                            lamvals = numpy.array([float(i) for i in lamvals])
                            lam = lamvals

                        else:
                            lam = (line.split(b'(')[-1].split(b')')[0].split(b','))
                            lam = numpy.array([float(i) for i in lam])
                    elif line.startswith(b'#Generating temperature'):
                        temp = line.split()[3]
                        unit = line.split()[4]
                        if unit == b'C':
                            T = float(temp)+273
                        else:
                            T = float(temp)
                        if T_previous is None:
                            T_previous = T
                            continue
                        else:
                            if T_previous !=T:
                                print ('Generating temperature is %.2f and does not match with other temperature given %.2f, '
                                       'please make sure all simulations are run at the same temperature' %(T, T_previous))
                                sys.exit(-1)
                            else:
                                T_previous = T
                        break
            # Normal file
            else:
                for line in open(f):
                    if line.startswith('#Alchemical'):
                        if lamvals is None:
                            lamvals = (line.split('(')[-1].split(')')[0].split(','))
                            lamvals = numpy.array([float(i) for i in lamvals])
                            lam = lamvals

                        else:
                            lam = (line.split('(')[-1].split(')')[0].split(','))
                            lam = numpy.array([float(i) for i in lam])
                    elif line.startswith('#Generating temperature'):
                        temp = line.split()[3]
                        unit = line.split()[4]
                        if unit == 'C':
                            T = float(temp)+273
                        else:
                            T = float(temp)
                        if T_previous is None:
                            T_previous = T
                            continue
                        else:
                            if T_previous !=T:
                                print ('Generating temperature is %.2f and does not match with other temperature given %.2f, '
                                       'please make sure all simulations are run at the same temperature' %(T, T_previous))
                                sys.exit(-1)
                            else:
                                T_previous = T
                        break
            if not numpy.array_equal(lam, lamvals):
                print ("Lambda arrays do not match! Make sure your input data is consistent")
                print (lam)
                print (lamvals)
                sys.exit(-1)
        return T_previous, lamvals

    def _read_data(self, input_files):
        r"""reads data from simfile.dat
        Parameters:
        -----------
        input_files : list of files
            list of files containing paths to simfile.dat

        Returns:
        --------
        data : numpy nd-array
            array containing all simfile.dat file data. 
        """
        data = []
        if self._discard is None:
            for f in input_files:
                data.append(numpy.loadtxt(f))
        else:
            print ('Loading data and discarding %d frames' %self._discard)
            for f in input_files:
                data.append(numpy.loadtxt(f, skiprows=self._discard))
        return data

    def _run_preprocessing(self, data, lamvals):
        r""" preprocess data files
        Parameters:
        -----------
        data : nd-array
            data from simfiles.dat
        lamvals : 1D numpy array
            list of lambda values from the simulation

        Returns:
        --------
        free_energy_obj = FreeEnergy object
            can be used to run MBAR and TI analysis

        """
        # N_k is the number of samples at generating thermodynamic state (lambda) k
        N_k = numpy.zeros(shape=lamvals.shape[0], dtype='int32')
        for k in range(0, lamvals.shape[0]):
            N_k[k] = data[k].shape[0]

        max_sample = int(max(N_k))
        grad_kn = numpy.zeros(shape=(lamvals.shape[0], max_sample))
        energies_kn = numpy.zeros(shape=(lamvals.shape[0], max_sample))

        for k in range(0, N_k.shape[0]):
            grad_kn[k, 0:N_k[k]] = data[k][:, 2]  # get the gradient information
            energies_kn[k, 0:N_k[k]] = data[k][:, 1]  # get the potential energies from the file.

        # Are the reduced perturbed potential energies generated at thermodynamic state k evaluated at state l, over
        # all n samples. This information is contained as is in the simulation file.
        u_kln = numpy.zeros(shape=(lamvals.shape[0], lamvals.shape[0], max_sample))
        for k in range(0, lamvals.shape[0]):
            u_kln[k, :, 0:N_k[k]] = data[k][:, 5:].transpose()

        # now we use the subsampling information to subsample the data.
        subsample_obj = SubSample(grad_kn, energies_kn, u_kln, N_k, subsample=self._subsample)
        subsample_obj.subsample_energies()
        subsample_obj.subsample_gradients()

        free_energy_obj = FreeEnergies(subsample_obj.u_kln, subsample_obj.N_k_energies, lamvals,
                                       subsample_obj.gradients_kn)
        return free_energy_obj

    def _get_boolean_state(self, s):
        r"""assess if a string is True or False
        Parameters:
        -----------
        s : string
            string name either True or False
        Reuturns:
        ---------
        b : boolean

        """
        if s == "True":
            return True
        elif s == "False":
            return False
        else:
            print ("The string you have supplied does not seem to contain boolean information")

    def _write_overlap_info(self,perturbation, M, sim_type=''):
        r"""computes and saves overlap matrix plots and info
        Parameters:
        -----------
        perturbation : string
            string identifying the perturbation for which the overlap is being computed
        M : 2D numpy array
            matrix containing the overlap information
        sim_type : string
            Identifier for bound type or free type simulation
        """
        diag_elements = numpy.array([numpy.diag(M, k=1), numpy.diag(M, k=-1)])
        fname = os.path.join(self._outputdir,perturbation)+'_'+sim_type+'_matrix.dat'
        fh = open(fname, 'wb')
        fh.write(bytes('#Overlap matrix\n', "UTF-8"))
        if numpy.min(diag_elements) < 0.03:
            fh.write(bytes('#Off diagonal elements of the overlap matrix are smaller than 0.03! Your free energy estimate is \n'
                        '#not reliable!\n', "UTF-8"))
        numpy.savetxt(fh, M, fmt='%.4f')
        fh.close()
        #now plotting the matrix
        fplot = os.path.join(self._outputdir,perturbation)+'_'+sim_type+'_matrix.png'
        # turns off interactive plotting
        plt.ioff()
        plt.matshow(M)
        plt.colorbar()
        plt.savefig(fplot, dpi=100)   

    @property
    def perturbation_list(self):
        r""""
        returns the list of pertubration directories
        """
        if self._perturbation_list == None:
            self._perturbation_list = os.listdir(self._basedir)
        return self._perturbation_list