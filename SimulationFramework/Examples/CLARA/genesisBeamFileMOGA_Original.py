import os, errno
import numpy as np
import random
from scipy.constants import c

from ocelot.S2E_STFC import FEL_simulation_block
from ocelot.adaptors.genesis import generate_input, get_genesis_launcher, run_genesis
from ocelot.gui.genesis_plot import fwhm3
from ocelot.common.math_op import *

import matplotlib.pyplot as plt
import time

from deap import base, creator, tools, algorithms
import copy

# make folder for genesis output files
os.system('mkdir genesis_out')
# specify location of output (usually pwd)
outputdir = '/scratch2b/djd63/MOGA/ForNeil/'

# Overwriting FEL simulation block to avoid plots etc. - returns Genesis output object
########################################################
class FEL_sim(FEL_simulation_block.FEL_simulation_block):
    def __init__(self,*initial_data,**kwargs):
       super(FEL_sim,self).__init__(*initial_data,**kwargs)

    def GEN_simul_preproc(self,A_input,i_aft=0):
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

        print('++++ Output Path {0} ++++++'.format(self.file_pout))

        # Setting the number of noise realisations and scan (over quads or different parameters)

        if (self.i_scan ==0):
            s_scan = range(1)
            num = self.stat_run
            run_ids = xrange(0,num)
            print('++++++++ No scan ++++++++++')
        elif (self.i_scan !=0):
            if (self.parameter in A_und):
                run_ids= range(1)
                if self.parameter !='aw0':
                    s_scan = range(int(self.init),int(self.end),int(np.ceil((self.end-self.init)/(self.n_scan))))
                else:
                    s_scan = np.linspace(self.init,self.end,self.n_scan)
                print('++++ Quad scan, parameter  {0} ++++++'.format(self.parameter))
            elif (self.parameter=='xlamds'):
                run_ids= range(1)
                s_scan = np.linspace(self.init,self.end,self.n_scan)
                print('++++ Quad scan, parameter  {0} ++++++'.format(self.parameter))
            else:
                s_scan = np.linspace(self.init,self.end,self.n_scan)
                num = self.stat_run
                run_ids = xrange(0,num)
                print('++++ Number of noise realisations {0} ++++++'.format(num))

            # setting the undulator design( Magnetic Lattice)
        A_undulator = self.undulator_design(A_input)

            # Fill in the beam object
        A_beam = self.BeamDefinition(A_input)
        if (getattr(A_input,'itdp')==0):
            print('++++ Steady State run +++++')
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
            print('++++ No edist or beam or dpa or rad file available ++++++')
        print inp.beam

        # Read ASTRA file.
        if hasattr(self,'i_astra') and getattr(self,'i_astra')==1 and hasattr(self,'astra_file'):
            inp=self.convert_ASTRA_edist(inp)
            setattr(inp,'beam',None)
        elif  (hasattr(self,'i_astra') and getattr(self,'i_astra')==1) and not (hasattr(self,'astra_file')):
            print('Path of  the ASTRA file not provided')
            return
        else:
            print('No need to read ASTRA file')

        # Rematch beam (if the flag has been set within the data dictionary)
        if (getattr(inp,'edist')!= None) and hasattr(self,'i_match') and (getattr(self,'i_match')==1):
            inp = self.rematch_edist(inp)

        # Setting up the time window of the distribution
        #if getattr(self,'i_beam')==0:
        #    setattr(inp,'nslice',getattr(A_input,'nslice'))
        #else:
        #    beamf = getattr(inp,'beam')
        #    if hasattr(beamf,'I'):
        #        setattr(inp,'curlen', getattr(beamf,tpulse) * speed_of_light / 1e15)
        #        setattr(inp,'nslice',8 * int(inp.curlen / inp.zsep / inp.xlamds))
        #    else:
        #        setattr(inp,'nslice',getattr(A_input,'nslice'))

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
                inp.lout =  [1,1,1,1,1,0,1,1,1,1,1,0,0,1,0,0,0,0,0]
                if ((self.stat_run==1) or (self.i_scan==1)):
                    setattr(inp,'ipseed',-1)
                else:
                    ipseed = np.random.randint(9999)
                    setattr(inp,'ipseed', ipseed)
                inp.run_dir = getattr(self,'file_pout')+'scan_'+str(n_par)+'/ip_seed_'+str(inp.ipseed)+'/'

                try:
                    os.makedirs(inp.run_dir)
                    if self.stat_run>1:
                        while ipseed in  [int(files[8:]) for files in os.listdir(getattr(self,'file_pout')+'scan_'+str(n_par)) if os.path.isdir(getattr(self,'file_pout')+'scan_'+str(n_par)+'/'+files)]:
                            ipseed = np.random.randint(9999)
                            shutil.rmtree(inp.run_dir)
                        else:
                            setattr(inp,'ipseed',ipseed)
                            inp.run_dir = getattr(self,'file_pout')+'scan_'+str(n_par)+'/ip_seed_'+str(inp.ipseed)+'/'
                            os.makedirs(inp.run_dir)
                except OSError as exc:
                    if (exc.errno == errno.EEXIST) and os.path.isdir(self.file_pout+'scan'+str(n_par)+'/run_'+str(run_id)):
                        pass
                    else:
                        raise
                if self.i_scan==1 and inp.f1st==1:
                    inp= self.GEN_scan(n_par ,A_input,A_undulator,inp)
                elif self.i_scan==0 and inp.f1st==1:
                    inp.lat = A_undulator['Magnetic Lattice']
                    setattr(inp,'magin',1)
                else:
                    inp.lat =None
                    setattr(inp,'magin',0)
                # DJD 16/10/18 temporary hack to randomise seed
                ipseed = np.random.randint(9999)
                setattr(inp,'ipseed', ipseed)

                inp_arr.append(inp)
                launcher=get_genesis_launcher(self.gen_launch)
                print('+++++ Starting simulation of noise realisation {0}'.format(run_id))
                g = run_genesis(inp,launcher,i_aft=i_aft)
                setattr(g,'filePath',str(inp.run_dir))
                #if (inp.itdp==1):
                #    plot_gen_out_all(handle=g, savefig=True, showfig=False, choice=(1, 1, 1, 1, 3.05, 0, 0, 0, 0, 0, 0, 0, 0), vartype_dfl=complex128, debug=1)
                inp.latticefile=None
                inp.outputfile=None
                inp.edistfile=None
                inp.beamfile=None
                inp.fieldfile=None
                inp.radfile=None
                #print g.__dict__.keys()

        #print('++++ postprocessing (results saved in {0}results++++'.format(self.file_pout))
        #self.post_processing(inp,s_scan)
        #plt.close("all")
        #print inp.out
        #return out
        return g

########################################################
# Generate s=array
def generate_sArray(nslice, win_len):
    step = win_len/nslice
    s = np.arange(0, win_len+step, step)
    #print s
    #print len(s)
    return s

def generateInitCur(nslice, av_current):
    # make array of random numbers between 0 and 1
    rands = np.random.uniform(0,1,nslice+1)
    return rands
    #print rands
    #print len(rands)

    # normalise so that every list of rands adds to 1
    sum_rands = sum(rands)
    #print sum_rands
    norm_rands = rands/sum_rands
    #print norm_rands
    sum_norm_rands = sum(norm_rands)
    #print sum_norm_rands
    current = norm_rands*av_current*nslice
    print current
    print current.tolist()
    time.sleep(5)
    return current.tolist()

def gaussian():
    # copied in elements for quick test
    a = [0.0031, 0.0060, 0.0111, 0.0198, 0.0340, 0.0561, 0.0889, 0.1353, 0.1979, 0.2780, 0.3753, 0.4868, 0.6065, 0.7261, 0.8353, 0.9231, 0.9802, 1.0000, 0.9802, 0.9231, 0.8353, 0.7261, 0.6065, 0.4868, 0.3753, 0.2780, 0.1979, 0.1353, 0.0889, 0.0561, 0.0340, 0.0198, 0.0111, 0.0060, 0.0031, 0.0015]
    #print type(a)
    #return np.array(a)
    return a

def write_beamfile(filename, zpos, curpeak):
    print '\nIn beamfile\n'
    print zpos
    print type(zpos), len(zpos)
    print curpeak
    print type(curpeak), len(curpeak)
    cur, = curpeak
    print cur
    #print curpeak[0]
    #print type(curpeak[0]), len(curpeak[0])
    #print zpos
    #print curpeak
    with open(filename, 'w') as beam_file:
        beam_file.write('? VERSION=1.0\n')
        beam_file.write('? SIZE = '+str(len(zpos))+'\n')
        beam_file.write('? COLUMNS ZPOS CURPEAK\n')
        for i in range(len(zpos)):
            #beam_file.write(str(zpos[i])+' '+str(curpeak[i])+'\n')
            beam_file.write('%.7E %.7E\n' % (zpos[i], cur[i]))


# Constants:
# Working with window length = 448um (xlambds=100e-9)*(zsep=8)*(nslice=560)
xlambds = 100e-9
zsep = 8
nslice = 560
window_length = xlambds*zsep*nslice
print window_length
# sub-divide into 560/16=35 slices
num_slices = nslice/16
print num_slices
charge = 250e-12
av_current = charge/(window_length/c)
print av_current

# generate s-array
s = generate_sArray(num_slices, window_length)
# generate current array
#I = generateInitCur(num_slices, av_current)
# write beamfile
#beamfile = 'beamfile.txt'
#write_beamfile(beamfile, s, I)

def evalBeamWithGenesis(ind, zpos=s):
    print 'in evalBeamWithGenesis\n\n\n'
    print ind

    # copy file to backup location with random file name (needs better identifier...)
    if os.path.exists('scan_0/ip_seed_-1/run.0.gout') == True:
        os.system('cp scan_0/ip_seed_-1/run.0.gout genesis_out/run%i.out' % random.randint(10000,99999))
    os.system('rm -R scan_0')

    # normalise so that every list of rands adds to 1, then write beamfile
    sum_rands = sum(ind)
    norm_rands = ind/sum_rands
    sum_norm_rands = sum(norm_rands)
    current = norm_rands*av_current*num_slices
    write_beamfile('beamfile.txt', zpos, current)

    # Genesis run:
    data={'gen_file':'NEW_SHORT_NOM_TD_v7.in',
          'file_pout':outputdir,
          'file_beam':'beamfile.txt',
          'i_scan':0,
          'gen_launch':'genesis2.0',
          'stat_run':1,
          'idump':0,
          'i_edist':0,
          'i_beam':1,
          'i_rad':0,
          'i_dpa':0,
          'i_dfl':0}
    #f = FEL_simulation_block.FEL_simulation_block(data)
    f = FEL_sim(data)
    A_inp=f.read_GEN_input_file()
    g = f.GEN_simul_preproc(A_inp)

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

    #return Emax
    #print spectrum_lamdwidth_fwhm[-1]
    print 'bandwidth = ', spectrum_lamdwidth_std[-1]
    print 'pulse energy =', g.energy[-1]
    return g.energy[-1], spectrum_lamdwidth_std[-1]

evalBeamWithGenesis([[0.0031, 0.006, 0.0111, 0.0198, 0.034, 0.0561, 0.0889, \
0.1353, 0.1979, 0.278, 0.3753, 0.4868, 0.6065, 0.7261, 0.8353, 0.9231, 0.9802, \
1.0, 0.9802, 0.9231, 0.8353, 0.7261, 0.6065, 0.4868, 0.3753, 0.278, 0.1979, \
0.1353, 0.0889, 0.0561, 0.034, 0.0198, 0.0111, 0.006, 0.0031, 0.0015]])
exit()
# create fitness class
# weights apply to maximise pulse energy and minimise bandwidth
creator.create("Fitness", base.Fitness, weights=(1.0, -1.0))

# create individual
creator.create("Individual", list, fitness=creator.Fitness)
#creator.create("Individual", generateInitCur(num_slices, av_current), fitness=creator.Fitness)

toolbox = base.Toolbox()
#toolbox.register("individual", creator.Individual, generateInitCur, num_slices, av_current)
#toolbox.register("genInitCur", generateInitCur, num_slices, av_current)
#toolbox.register("individual", creator.Individual, toolbox.genInitCur)
#toolbox.register("individual", creator.Individual)
#toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.genInitCur, 1)
#toolbox.register("attr_rand", random.uniform, 0, 1)
#toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.attr_rand, 36)
toolbox.register("gaussian", gaussian)
#toolbox.register("individual", creator.Individual, toolbox.gaussian )
toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.gaussian, 1)
toolbox.register("population", tools.initRepeat, list, toolbox.individual)

# crossover function interleaves even and odd elements of parents
def cx(ind1, ind2):
    print '\nin cx\n'
    #return ind1, ind2
    temp = copy.deepcopy(ind1)
    ind1[0][1::2]=ind2[0][1::2]
    ind2[0][1::2]=temp[0][1::2]
    return ind1, ind2

# mutation function changes a proportion of the elements
def mut(ind, nslice=num_slices, av_current=av_current):
    print '\n in mut\n'
    # can read from file if want to adjust during the run
    #with open('mutFlipBitProb.txt', 'r') as f:
    #    mutfbp = float(f.read())
    mutfbp = 0.05
    for i in range(len(ind[0])):
        #print i
        if random.random() < mutfbp:
            #print random.random()
            ind[0][i] =0.25*(3*ind[0][i]+random.random())#2*random.random()/len(ind)
            #print ind[i]
    return ind,

#toolbox.register("evalfn", evalBeamWithGenesis, toolbox.individual)
#toolbox.register("evaluate", toolbox.evalfn)
toolbox.register("evaluate", evalBeamWithGenesis)
toolbox.register("mate", cx)
toolbox.register("mutate", mut)
toolbox.register("select", tools.selNSGA2)

def main():
    print 'Starting main'
    NGEN=100#100#3#30
    MU=50#2#10
    LAMBDA=100#3#15
    CXPB=0.2
    MUTPB=0.7

    pop = toolbox.population(n=MU)
    print 'here', pop[0]
    print pop[1]
    #print 'fitness', pop[0].fitness
    #print pop[0][0]
    #time.sleep(5)
    hof = tools.ParetoFront()
    stats = tools.Statistics(lambda ind: ind.fitness.values)
    stats.register("avg", np.mean, axis=0)
    stats.register("std", np.std, axis=0)
    stats.register("min", np.min, axis=0)
    stats.register("max", np.max, axis=0)

    pop, stats, hof = algorithms.eaMuPlusLambdaPopPerGen(pop,toolbox,MU,LAMBDA,CXPB,MUTPB,NGEN,stats,halloffame=hof)
    return pop, stats, hof

main()
#pop, stats, hof = main()
#print pop
#print stats
#print hof
