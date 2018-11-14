import os, errno, sys
import numpy as np
import random
from scipy.constants import c
from ocelot.S2E_STFC import FEL_simulation_block
from ocelot.adaptors.genesis import generate_input, get_genesis_launcher, run_genesis, rematch_edist, edist2beam
from ocelot.gui.genesis_plot import fwhm3
from ocelot.common.math_op import *
from ocelot.cpbd.beam import Twiss
sys.path.append(os.path.abspath(__file__+'/../../../../'))
import SimulationFramework.Modules.read_beam_file as rbf
from SimulationFramework.Modules.optimiser import optimiser
opt = optimiser()
from SimulationFramework.Modules.constraints import constraintsClass
import matplotlib.pyplot as plt
import time
import csv
import multiprocessing
from scoop import futures
from deap import base, creator, tools, algorithms
import copy

# make folder for genesis output files
# os.system('mkdir genesis_out')
# specify location of output (usually pwd)
# outputdir = '/home/jkj62/Github/OnlineModel/SimulationFramework/Examples/CLARA/MOGA_Test/'

# Overwriting FEL simulation block to avoid plots etc. - returns Genesis output object
########################################################
class FEL_sim(FEL_simulation_block.FEL_simulation_block):
    def __init__(self,*initial_data,**kwargs):
       super(FEL_sim,self).__init__(*initial_data,**kwargs)

    def convert_HDF5_edist(self, inp, filename='CLA-FMS-APER-01.hdf5', center=True):
        from ocelot.adaptors.genesis import GenesisElectronDist
        inp_out = inp
        beam = rbf.beam()
        beam.read_HDF5_beam_file(filename)
        beam.rematchXPlane(beta=3.76, alpha=-0.189)
        # print 'beta_x = ', beam.beta_x
        beam.rematchYPlane(beta=1.45, alpha=0)
        # print 'beta_y = ', beam.beta_y
        edist = GenesisElectronDist()
        edist.x = beam.x
        edist.y = beam.y
        # edist.t = -(adist[:, 2] - np.mean(adist[:, 2])) / speed_of_light  # long position normalized to 0 and converted to time #HMCC add -
        edist.t = -1*beam.t
        edist.t = edist.t - np.mean(edist.t)
        edist.xp = beam.xp
        edist.yp = beam.yp
        p_tot = beam.p
        edist.g = beam.gamma
        edist.part_charge = 250e-12 / len(beam.x)

        if center:
            edist.x -= np.mean(edist.x)
            edist.y -= np.mean(edist.y)
            edist.xp -= np.mean(edist.xp)
            edist.yp -= np.mean(edist.yp)
        setattr(edist,'filePath',getattr(self,'file_pout')+'read_edist_astra')
        setattr(inp_out,'edist',edist)
        return inp_out, beam

    def GEN_simul_preproc(self, A_input, i_aft=0, dirno=1):
        if not self.file_pout.endswith('/'):
            self.file_pout=self.file_pout+'/'
        #print('\n\n\nis it working?\n\n')
        inp_arr = []
        A_bbeam = ['rxbeam','rybeam','emitx','emity','alphax','alphay','xbeam','ybeam','pxbeam','pybeam']
        A_simul = ['alignradf','npart','ncar','zsep','delz','dmpfld','fbess0','dgrid','rmax0','xkx','xky','iwityp',
                   'nptr','lbc','zrayl','zstop','zwaist','delgam','xlamd','nscr','nscz','curpeak',
                   'iphsty','nharm','curlen','nbins','gamma0','isravg','isrsig','eloss','version',
                   'multconv','imagl','convharm','idril','ffspec','ndcut','ibfield','nslice','ntail',
                   'ippart','ispart','ipradi','isradi']
        A_td = ['itdp','prad0','shotnoise']
        A_und = ['quadd', 'quadf','fl','dl','drl','nsec','nwig','aw0', 'awd']

        # print('++++ Output Path {0} ++++++'.format(self.file_pout))

        # Setting the number of noise realisations and scan (over quads or different parameters)

        if (self.i_scan ==0):
            s_scan = range(1)
            num = self.stat_run
            run_ids = xrange(0,num)
            # print('++++++++ No scan ++++++++++')
        elif (self.i_scan !=0):
            if (self.parameter in A_und):
                run_ids= range(1)
                if self.parameter !='aw0':
                    s_scan = range(int(self.init),int(self.end),int(np.ceil((self.end-self.init)/(self.n_scan))))
                else:
                    s_scan = np.linspace(self.init,self.end,self.n_scan)
                # print('++++ Quad scan, parameter  {0} ++++++'.format(self.parameter))
            elif (self.parameter=='xlamds'):
                run_ids= range(1)
                s_scan = np.linspace(self.init,self.end,self.n_scan)
                # print('++++ Quad scan, parameter  {0} ++++++'.format(self.parameter))
            else:
                s_scan = np.linspace(self.init,self.end,self.n_scan)
                num = self.stat_run
                run_ids = xrange(0,num)
                # print('++++ Number of noise realisations {0} ++++++'.format(num))

            # setting the undulator design( Magnetic Lattice)
        A_undulator = self.undulator_design(A_input)

            # Fill in the beam object
        A_beam = self.BeamDefinition(A_input)
        if (getattr(A_input,'itdp')==0):
            # print('++++ Steady State run +++++')
            i_tdp = False
        elif (getattr(A_input,'itdp')==1):
            i_tdp = True

            # Generate input object
        inp = generate_input(A_undulator['Undulator Parameters'],A_beam,itdp=i_tdp)

        # Overwrite the simulation attributes of the input object with the ones defined in the input file
        for key in A_input.__dict__:
            if (key in A_simul) or (key in A_und) or (key in A_td) or (key =='xlamds') or (key == 'f1st') or (key == 'nslice') or (key == 'ntail'):
                setattr(inp,key, getattr(A_input,key))
        for key in ['edist','beam','dfl']:
            if getattr(A_input,key)!=None:
                setattr(inp,key,getattr(A_input,key))

        # Set up some input parameters
        if getattr(inp,'itdp')==0:
            setattr(inp,'type','steady')
        else:
            setattr(inp,'type','tdp')
        setattr(inp, 'awd', float(getattr(inp, 'aw0')))

        # idump attribute
        if (getattr(self,'idump')) == 1:
            setattr(inp,'idump',1)
            setattr(inp,'idmpfld',1)

        # Existent dist or beam file (if exists)
        if (getattr(self,'i_edist') == 1) and (hasattr(self,'file_edist')):
            inp=self.GEN_existent_beam_dist_dpa_rad(inp,'edist')
        elif (getattr(self,'i_beam') == 1) and (hasattr(self,'file_beam')):
            #print inp.__dict__.keys()
            inp=self.GEN_existent_beam_dist_dpa_rad(inp,'beam')
            #print inp.__dict__.keys()
        elif  (getattr(self,'i_rad') == 1) and (hasattr(self,'file_rad')):
            inp=self.GEN_existent_beam_dist_dpa_rad(inp,'rad')
        elif  (getattr(self,'i_dpa') == 1) and (hasattr(self,'file_dpa')):
            inp=self.GEN_existent_beam_dist_dpa_rad(inp,'dpa')
        elif  (getattr(self,'i_dfl') == 1) and (hasattr(self,'file_dfl')):
            inp=self.GEN_existent_beam_dist_dpa_rad(inp,'dfl')
        else:
            pass
            # print('++++ No edist or beam or dpa or rad file available ++++++')
        # print inp.beam

        # Read HDF5 file.
        if hasattr(self,'i_HDF5') and getattr(self,'i_HDF5')==1 and hasattr(self,'HDF5_file'):
            inp, beam = self.convert_HDF5_edist(inp, getattr(self,'HDF5_file'))
            setattr(inp,'beam',None)
            energy = np.average(inp.edist.g)
            setattr(inp,'xlamds',float(inp.xlamd*(1.0+np.square(inp.aw0))/(2.0*np.square(energy))))
            print 'xlamds = ', inp.xlamds

        elif  (hasattr(self,'i_HDF5') and getattr(self,'i_HDF5')==1) and not (hasattr(self,'HDF5_file')):
            # print('Path of  the HDF5 file not provided')
            return
        else:
            pass
            # print('No need to HDF5 ASTRA file')

        if (getattr(self,'i_edist')==1) or (getattr(inp,'edist')!=None) or  (getattr(inp,'beam')!=None) :
            setattr(inp,'ntail',0)
        else:
            if (getattr(self,'i_edist')==0) and getattr(A_input,'ntail')!=0 :
                setattr(inp,'ntail',int(getattr(A_input,'ntail')))
            else:
                setattr(inp,'ntail',-int(np.floor(getattr(inp,'nslice')/2)))
        # Overwrite the simulation attributes if the user has new values for them defined in the input data structure
        if (hasattr(self, 'i_rewrite')) and (hasattr(self, 'par_rew')) and (getattr(self, 'i_rewrite') == 1):
            inp = self.GEN_rewrite_par(inp)
        else:
            pass

        # Running over noise realisations and/or scan parameters
        for n_par in s_scan:
            for run_id in run_ids:
                inp.runid = run_id
                inp.lout =  [1,1,1,1,1,0,1,1,1,1,1,1,1,1,0,0,0,0,0]
                inp.run_dir = getattr(self,'file_pout')

                try:
                    os.makedirs(inp.run_dir)
                except OSError as exc:
                    pass
                if self.i_scan==1 and inp.f1st==1:
                    inp= self.GEN_scan(n_par ,A_input,A_undulator,inp)
                elif self.i_scan==0 and inp.f1st==1:
                    inp.lat = A_undulator['Magnetic Lattice']
                    setattr(inp,'magin',1)
                else:
                    inp.lat =None
                    setattr(inp,'magin',0)
                # DJD 16/10/18 temporary hack to randomise seed
                ipseed = 56#np.random.randint(9999)
                setattr(inp,'ipseed', ipseed)

                inp_arr.append(inp)
                launcher=get_genesis_launcher(self.gen_launch)
                g = run_genesis(inp,launcher,i_aft=i_aft)
                setattr(g,'filePath',str(inp.run_dir))
        return g

########################################################

def evalBeamWithGenesis(dir):

    # Genesis run:
    data={'gen_file':'NEW_SHORT_NOM_TD_v7.in',
          'file_pout': dir,
          'file_beam':'beamfile.txt',
          'i_scan':0,
          'gen_launch':'/opt/OpenMPI-3.1.3/bin/mpiexec --timeout 300 -np 20 /opt/Genesis/bin/genesis2 < tmp.cmd 2>&1 > /dev/null',
          'stat_run':1,
          'idump':0,
          'i_edist':0,
          'i_beam':0,
          'i_rad':0,
          'i_dpa':0,
          'i_dfl':0,
          'i_HDF5':1,
          'HDF5_file': dir+'/'+'CLA-FMS-APER-01.hdf5',
          'i_match': 0}
    f = FEL_sim(data)
    A_inp=f.read_GEN_input_file()
    g = f.GEN_simul_preproc(A_inp, dirno=dir)

    # this section is copied from genesis_plot.py - to calculate bandwidth
    spectrum_lamdwidth_fwhm = np.zeros_like(g.z)
    spectrum_lamdwidth_std = np.zeros_like(g.z)
    for zz in range(g.nZ):
        spectrum_lamdwidth_fwhm[zz] = None
        spectrum_lamdwidth_std[zz] = None
        if np.sum(g.spec[:,zz])!=0:
            pos, width, arr = fwhm3(g.spec[:, zz])
            if width != None:
                if arr[0] == arr[-1]:
                    dlambda = abs(g.freq_lamd[pos] - g.freq_lamd[pos-1])
                else:
                    dlambda = abs( (g.freq_lamd[arr[0]] - g.freq_lamd[arr[-1]]) / (arr[0] - arr[-1]) )
                spectrum_lamdwidth_fwhm[zz] = dlambda * width / g.freq_lamd[pos]
                # spectrum_lamdwidth_fwhm[zz] = abs(g.freq_lamd[arr[0]] - g.freq_lamd[arr[-1]]) / g.freq_lamd[pos]  # the FWHM of spectral line (error when peakpos is at the edge of lamdscale)
            spectrum_lamdwidth_std[zz] = std_moment(g.freq_lamd, g.spec[:, zz]) / n_moment(g.freq_lamd, g.spec[:, zz], 0, 1)
    ######### end of section copied from genesis_plot.py

    return g.energy[-1], spectrum_lamdwidth_std[-1]

def saveState(args, fitness, *values):
    global csv_out
    args=list(args)
    for v in values:
        args.append(v)
    args.append(fitness)
    csv_out.writerow(args)
    # csv_out.flush()

global_best = 0
def runEleGenesis(inputargs, tmpdir, verbose=False, *args, **kwargs):
    global startingvalues, global_best
    try:
        fit = runEle.fitnessFunc(inputargs, tmpdir, *args, **kwargs)
        fitvalue = fit.calculateBeamParameters()
        sys.stdout = open(tmpdir+'/'+'std.out', 'w')
        sys.stderr = open(tmpdir+'/'+'std.err', 'w')
        e, b = evalBeamWithGenesis(tmpdir)
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        if e > global_best:
            global_best = e
            print 'New Best: ', 1e4*e, 1e2*b
        if verbose:
            print 'bandwidth = ', 1e2*b
            print 'pulse energy =', 1e4*e
        constraintsList = {
            'energy': {'type': 'greaterthan', 'value': abs(1e4*e), 'limit': 3, 'weight': 1},
            'bandwidth': {'type': 'lessthan', 'value': abs(1e2*b), 'limit': 1, 'weight': 1},
        }
        cons = constraintsClass()
        fitness = cons.constraints(constraintsList)
        saveState(inputargs, fitness, e, b)
        # print cons.constraintsList(constraintsList)

        return 1e4*e, 1e2*b
    except Exception as e:
        print 'Error! ', e
        return 0, 10

import runElegant as runEle
def optfunc(inputargs, verbose=False, dir=None, *args, **kwargs):
    process = multiprocessing.current_process()
    runno = process.pid
    # inputargs = [a*b for a,b in zip(startingvalues,inputargsmult)]
    # print 'inputargs = ', inputargs
    if dir == None:
        with runEle.TemporaryDirectory(dir=os.getcwd()) as tmpdir:
            return runEleGenesis(inputargs, tmpdir, verbose=verbose, *args, **kwargs)
    else:
        return runEleGenesis(inputargs, dir, verbose=verbose, *args, **kwargs)

def checkBounds(min, max):
    def decorator(func):
        def wrapper(*args, **kargs):
            offspring = func(*args, **kargs)
            for child in offspring:
                for i in xrange(len(child)):
                    if child[i] > max[i]:
                        child[i] = max[i]
                    elif child[i] < min[i]:
                        child[i] = min[i]
            return offspring
        return wrapper
    return decorator

startingvalues = best = [30000000.0, -23, 27000000.0, -8, 24000000.0, 184, 32000000.0, 45, -0.1185]
# best using genesis
# startingvalues = best = [2.5909474150098216e7,-17.061238835950345,2.6239817052311208e7,0.7941999561214287,1.802394629640912e7,173.6725806707975,2.8248813612251766e7,52.62317478360277,-0.17180888501850444]
#best from astra (no genesis)
# startingvalues = best = [27477812.476143554,-24.149079998274154,27185718.693607226,-8.078720301316636,24326877.45572499,188.9222526372106,30434366.8035302,45.17175610232019,-0.12770679325302878]
# print 'starting values = ', best
# optfunc(best, dir='testing', scaling=5, post_injector=True, verbose=True)
# exit()
# print best
startranges = [[0.95 * i, 1.05 * i] for i in best]
MIN = [0, -90, 0, -90, 0, 90, 0, -90, -0.2]
MAX = [32e6, 90, 32e6, 90, 45e6, 270, 32e6, 90, -0.1]

generateHasBeenCalled = False
def generate():
    global generateHasBeenCalled
    if not generateHasBeenCalled:
        generateHasBeenCalled = True
        return creator.Individual(list(best))
    else:
        return creator.Individual(random.uniform(a,b) for a,b in startranges)

# print generate()

creator.create("FitnessMin", base.Fitness, weights=(1.0, -1.0,))
creator.create("Individual", list, fitness=creator.FitnessMin)

toolbox = base.Toolbox()

# Attribute generator
toolbox.register("attr_bool", generate)

# Structure initializers
toolbox.register("Individual", generate)
toolbox.register("population", tools.initRepeat, list, toolbox.Individual)

if os.name == 'nt':
    toolbox.register("evaluate", optfunc, scaling=3, post_injector=True)
else:
    toolbox.register("evaluate", optfunc, scaling=5, post_injector=True)
toolbox.register("mate", tools.cxUniform, indpb=0.3)
toolbox.register("mutate", tools.mutGaussian, mu=0, sigma=[1e6,2,1e6,2,2e6,2,1e6,2,0.01], indpb=0.3)

toolbox.decorate("mate", checkBounds(MIN, MAX))
toolbox.decorate("mutate", checkBounds(MIN, MAX))

toolbox.register("select", tools.selNSGA2)


if __name__ == "__main__":
    global hof, global_best, csv_out
    random.seed(64)

    global_best = 0
    out = open('best_solutions_running_simplex_elegant_genesis.csv','a', buffering=0)
    csv_out=csv.writer(out)

    NGEN = 100
    MU = 24
    LAMBDA = 48
    CXPB = 0.7
    MUTPB = 0.2

    # Process Pool of 4 workers
    if not os.name == 'nt':
        pool = multiprocessing.Pool(processes=8)
    else:
        pool = multiprocessing.Pool(processes=2)
    toolbox.register("map", pool.map)

    if not os.name == 'nt':
        pop = toolbox.population(n=MU)
    else:
        pop = toolbox.population(n=MU)
    hof = tools.ParetoFront()
    stats = tools.Statistics(lambda ind: ind.fitness.values)
    stats.register("avg", np.mean, axis=0)
    stats.register("std", np.std, axis=0)
    stats.register("min", np.min, axis=0)
    stats.register("max", np.max, axis=0)

    algorithms.eaMuPlusLambda(pop, toolbox, MU, LAMBDA, CXPB, MUTPB, NGEN, stats,
                              halloffame=hof)

    pool.close()
    # print 'pop = ', pop
    # print logbook
    print hof

    # try:
    #     print 'best fitness = ', optfunc(hof[0], dir=os.getcwd()+'/CLARA_best_longitudinal_Genesis', scaling=5, overwrite=True, verbose=True, summary=True)
    #     with open('CLARA_best_longitudinal_Genesis/CLARA_best_longitudinal_Genesis_solutions.csv','wb') as out:
    #         csv_out=csv.writer(out)
    #         for row in hof:
    #             csv_out.writerow(row)
    #     with open('CLARA_best_longitudinal_Genesis/CLARA_best_longitudinal_Genesis_stats.csv','wb') as out:
    #         csv_out=csv.writer(out)
    #         for row in stats:
    #             csv_out.writerow(row)
    # except:
    with open('CLARA_best_longitudinal_Genesis_solutions_final.csv','wb') as out:
        hof_out=csv.writer(out)
        for row in hof:
            hof_out.writerow(row)
