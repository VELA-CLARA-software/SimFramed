import os, errno, sys
import numpy as np
import random
from scipy.constants import c
sys.path.append(os.path.abspath(__file__+'/../../../../../'))
from ocelot.S2E_STFC import FEL_simulation_block
from ocelot.adaptors.genesis import generate_input, get_genesis_launcher, run_genesis, rematch_edist, edist2beam
from ocelot.gui.genesis_plot import fwhm3
from ocelot.common.math_op import *
from ocelot.cpbd.beam import Twiss
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
import SimulationFramework.Examples.CLARA.Elegant.runElegant as runEle
from shutil import copyfile
import argparse

parser = argparse.ArgumentParser(description='Run Elegant + Genesis')
parser.add_argument('-g', '--gaussian', default=False)
parser.add_argument('-s', '--set', default=False)
parser.add_argument('-v', '--vbc', default=False)
parser.add_argument('-p', '--postinjector', default=True)

class FEL_sim(FEL_simulation_block.FEL_simulation_block):
    def __init__(self, initial_data, alphax=-0.189, betax=3.76, alphay=0, betay=1.44, **kwargs):
       print 'initial_data = ', (initial_data)
       super(FEL_sim,self).__init__((initial_data),**kwargs)
       self.alphax = alphax
       self.betax = betax
       self.alphay = alphay
       self.betay = betay

    def convert_HDF5_edist(self, inp, filename='CLA-FMS-APER-01.hdf5', center=True):
        from ocelot.adaptors.genesis import GenesisElectronDist
        inp_out = inp
        beam = rbf.beam()
        beam.read_HDF5_beam_file(filename)
        # print 'beta_x = ', beam.beta_x, ' alpha_x = ', beam.alpha_x
        # print 'beta_y = ', beam.beta_y, ' alpha_y = ', beam.alpha_y
        beam.rematchXPlanePeakISlice(beta=self.betax, alpha=self.alphax)
        beam.rematchYPlanePeakISlice(beta=self.betay, alpha=self.alphay)
        # print 'beta_x = ', beam.beta_x, ' alpha_x = ', beam.alpha_x
        # print 'beta_y = ', beam.beta_y, ' alpha_y = ', beam.alpha_y
        edist = GenesisElectronDist()
        edist.x = beam.x
        edist.y = beam.y
        # edist.t = -(adist[:, 2] - np.mean(adist[:, 2])) / speed_of_light  # long position normalized to 0 and converted to time #HMCC add -
        edist.t = beam.t
        edist.t = edist.t - np.mean(edist.t)# - 1.80e-13
        edist.xp = beam.xp
        edist.yp = beam.yp
        p_tot = beam.p
        edist.g = beam.gamma
        edist.part_charge = abs(self.startcharge)*1e-12 / len(beam.x)
        print 'self.startcharge = ', self.startcharge
        self.bunch_length = np.std(beam.t)

        if center:
            edist.x -= np.mean(edist.x)
            edist.y -= np.mean(edist.y)
            edist.xp -= np.mean(edist.xp)
            edist.yp -= np.mean(edist.yp)
        setattr(edist,'filePath',getattr(self,'file_pout')+'read_edist_astra')
        setattr(inp_out,'edist',edist)
        return inp_out, beam

    def GEN_simul_preproc(self, A_input, i_aft=0, dirno=1, startcharge=250):
        self.startcharge = startcharge
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
            gamma = np.average(inp.edist.g)
            # setattr(inp,'xlamds',float(inp.xlamd*(1.0+np.square(inp.aw0))/(2.0*np.square(gamma))))
            setattr(inp,'xlamds', float(0.022058051560136/(gamma**2))) ## DJD 15/02/2019
            print 'aw0 = ', inp.aw0
            print 'awd = ', inp.awd
            print 'xlamds = ', inp.xlamds

        elif  (hasattr(self,'i_HDF5') and getattr(self,'i_HDF5')==1) and not (hasattr(self,'HDF5_file')):
            # print('Path of  the HDF5 file not provided')
            return
        else:
            pass
            # print('No need to HDF5 ASTRA file')

        # if (getattr(self,'i_edist')==1) or (getattr(inp,'edist')!=None) or  (getattr(inp,'beam')!=None) :
        #     setattr(inp,'ntail',0)
        # else:
        #     if (getattr(self,'i_edist')==0) and getattr(A_input,'ntail')!=0 :
        #         setattr(inp,'ntail',int(getattr(A_input,'ntail')))
        #     else:
        #         setattr(inp,'ntail',-int(np.floor(getattr(inp,'nslice')/2)))
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
                except:
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
                # g.Lsat = g.z[np.argmax(g.bunching[np.argmax(g.p_int[:,-1])])]
                g.bunch_length = self.bunch_length
                g.Lsat = g.z[np.argmax(g.bunching[np.argmax(g.I)])]
        return g

########################################################

def find_saturation(power, z, n_smooth=20):
    p = np.diff(np.log10(power))

    u = np.convolve(p, np.ones(n_smooth) / float(n_smooth), mode='same')
    um = np.max(u)

    ii = 0

    for i in range(len(u)):
        if u[i] < 0.1 * um and z[i] > 5:
            ii = i
            print 'break! i = ', ii
            break

    #plt.plot(g.z[1:], u, lw=3)
    #plt.plot(g.z[ii+1], p[ii], 'rd')

    #plt.plot(g.z, power, lw=3)
    #plt.plot(z[ii+1], np.log10(power[ii]), 'rd')

    return z[ii+1], ii+1

def evalBeamWithGenesis(dir, run=True, startcharge=250):

    # Genesis run:
    data={'gen_file':'NEW_SHORT_NOM_TD_v7.in',
          'file_pout': dir,
          'file_beam':'beamfile.txt',
          'i_scan':0,
          'gen_launch':'',
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
    if run:
        data['gen_launch'] = '/opt/OpenMPI-3.1.3/bin/mpiexec --timeout 300 -np 44 /opt/Genesis/bin/genesis2 < tmp.cmd 2>&1 > /dev/null'
    else:
        data['gen_launch'] = ''
    f = FEL_sim(data)
    A_inp = f.read_GEN_input_file()
    g = f.GEN_simul_preproc(A_inp, dirno=dir, startcharge=startcharge)

    p = map(np.mean, zip(*g.p_int))
    g.Lsat, g.Lsatindex = find_saturation(p, g.z)

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
    brightness = g.energy / spectrum_lamdwidth_std

    for i in range(len(g.z)):
        if g.z[i] < 5:
            brightness[i] = 0
    max = np.argmax(brightness)
    maxe = np.argmax(g.energy)
    g.spectrum_lamdwidth_std = spectrum_lamdwidth_std
    return g.energy[max], g.spectrum_lamdwidth_std[max], g.z[max], g.energy[maxe], g.spectrum_lamdwidth_std[maxe], g.z[maxe], g.bunch_length, g

csv_out = ''

def saveState(args, id, *values):
    global csv_out
    args=list(args)
    for v in values:
        args.append(v)
    args.append(id)
    csv_out.writerow(args)
    # csv_out.flush()


def optfunc(inputargs, verbose=True, dir=None, savestate=True, runGenesis=True, runElegant=True, post_injector=True, *args, **kwargs):
    global global_best
    process = multiprocessing.current_process()
    runno = process.pid

    with runEle.TemporaryDirectory(dir=os.getcwd()) as tmpdir:
        if dir is not None:
            tmpdir = dir
            if not os.path.exists(tmpdir):
                os.makedirs(tmpdir)
        # try:
        if (not hasattr(inputargs, 'id')) or (hasattr(inputargs, 'id') and inputargs.id is None):
            idNumber = os.path.basename(tmpdir)
        else:
            idNumber = inputargs.id

        # print 'post_injector = ', post_injector

        if len(inputargs) == 9:
            inputargs = np.append(inputargs, 250)
        elif len(inputargs) == 15:
            inputargs = np.append(inputargs, 250)
            # print 'inputargs = ', inputargs
        sys.stdout = open(tmpdir+'/'+'std.out', 'w')
        sys.stderr = open(tmpdir+'/'+'std.err', 'w')
        if runElegant:
            if post_injector:
                fit = runEle.fitnessFunc(inputargs[:9], tmpdir, id=idNumber, startcharge=inputargs[9], post_injector=True, *args, **kwargs)
            else:
                print 'inputargs = ', inputargs
                fit = runEle.fitnessFunc(inputargs[:15], tmpdir, id=idNumber, startcharge=inputargs[15], post_injector=False, *args, **kwargs)
            fitvalue = fit.calculateBeamParameters()
        if post_injector:
            e, b, l, ee, be, le, bunchlength, g = evalBeamWithGenesis(tmpdir, run=runGenesis, startcharge=inputargs[9])
        else:
            e, b, l, ee, be, le, bunchlength, g = evalBeamWithGenesis(tmpdir, run=runGenesis, startcharge=inputargs[15])
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        if verbose:
            print 'bandwidth = ', 1e2*b, '  pulse energy =', 1e6*e, '  Sat. Length =', l
            print 'bandwidth E = ', 1e2*be, '  max pulse energy =', 1e6*ee, '  Sat. Length E =', le
        #### Save Output files to dir named after runno #####
        if not os.path.exists('./outData/'):
            os.makedirs('./outData/')
        dir = './outData/' + str(idNumber)
        if not os.path.exists(dir):
            os.makedirs(dir)
        copyfile(tmpdir+'/CLA-FMS-APER-01.hdf5',dir+'/CLA-FMS-APER-01.hdf5')
        copyfile(tmpdir+'/run.0.gout',dir+'/run.0.gout')
        copyfile(tmpdir+'/std.out',dir+'/std.out')
        if savestate:
            try:
                saveState(inputargs, idNumber, e, b, l, ee, be, le, bunchlength)
            except:
                pass
        return 1e4*e, 1e2*b, 1e4*ee, 1e2*be, l, g
        # except Exception as e:
            # print 'Error! ', e
            # return 0, 10, 0, 10


if __name__ == "__main__":
    args = parser.parse_args()

    injector_startingvalues = [-8.906156010951616,0.3420474160090586,2.0515744815221354e7,-16.281405933324855,0.05036027437405955,-0.0502414403704962]
    startingvalues = best = [1.879700695321506e7,-27.832442932602543,2.583534801151637e7,-1.2644316016467563,1.446505414414888e7,181.03888866546697,
                             3.1987431329329092e7,44.128256932519484,-0.1560424081528136]

    if args.postinjector == False or args.postinjector == 'False' or args.postinjector == 0:
        best = injector_startingvalues + startingvalues
        POST_INJECTOR = False
    else:
        best = startingvalues
        POST_INJECTOR = True
    print 'Post Injector = ', POST_INJECTOR, ' [', args.postinjector,']'

    global_best = 0
    from SimulationFramework.Framework import Framework
    framework = Framework('longitudinal_best', overwrite=False)
    framework.loadSettings('Lattices/clara400_v12_v3_elegant_jkj.def')
    parameters = []
    parameters.append(framework.getElement('CLA-L02-CAV', 'field_amplitude'))
    parameters.append(framework.getElement('CLA-L02-CAV', 'phase'))
    parameters.append(framework.getElement('CLA-L03-CAV', 'field_amplitude'))
    parameters.append(framework.getElement('CLA-L03-CAV', 'phase'))
    parameters.append(framework.getElement('CLA-L4H-CAV', 'field_amplitude'))
    parameters.append(framework.getElement('CLA-L4H-CAV', 'phase'))
    parameters.append(framework.getElement('CLA-L04-CAV', 'field_amplitude'))
    parameters.append(framework.getElement('CLA-L04-CAV', 'phase'))
    parameters.append(framework['bunch_compressor']['angle'])

    def gaussian_beam_test():
        ## This is for the Gaussian Beam Test!
        master_subdir = 'gaussianBeam'
        beam = rbf.beam()
        elegantbeamfilename = 'CLA-FMS-APER-01.sdds'
        beam.read_SDDS_beam_file(master_subdir + '/' + elegantbeamfilename, charge=250e-12)
        beam.beam['total_charge'] = 250e-12
        HDF5filename = elegantbeamfilename.replace('.sdds','.hdf5').replace('.SDDS','.hdf5').strip('\"')
        beam.write_HDF5_beam_file(master_subdir + '/' + HDF5filename, centered=False, sourcefilename=elegantbeamfilename)
        optfunc(best, dir=os.path.abspath('gaussianBeam'), scaling=6, post_injector=POST_INJECTOR, verbose=True, savestate=False, runGenesis=True, runElegant=False)

    # 1351 == 50uJ
    set253 = [2.3133382987191193e7,-19.372618068598406,3.0612137974426482e7,-3.772434513198637,2.198461283242174e7,173.7485023073137,3.13847229015289e7,42.038068437372736,-0.15496116784319712, 150]
    def npart_test():
        ## This is for the set Beam Test!
        for n in [253]:
            best = eval('set'+str(n))
            print n,' = ', best
            optfunc(best, dir=os.path.abspath('set'+str(n)+'_32k'), scaling=5, post_injector=POST_INJECTOR, verbose=True, savestate=False)
            # optfunc(best, dir=os.path.abspath('set'+str(n)+'_262k'), scaling=6, post_injector=True, verbose=True, savestate=False)

    def scan_vbc_angle():
        for vbc in np.arange(startingvalues[-1]-0.01, startingvalues[-1]+0.01,0.002):
            newbest = best
            newbest[-1] = vbc
            print 'VBC angle = ', vbc
            optfunc(newbest, dir='simplex/vbc_angle_'+str(vbc)+'_262k', scaling=6, post_injector=POST_INJECTOR, verbose=True, savestate=False)

    if args.gaussian or args.gaussian == 'True' or args.gaussian > 0:
        print 'Gaussian Beam'
        gaussian_beam_test()
    elif args.set or args.set == 'True' or args.set > 0:
        print 'Set Test'
        npart_test()
    elif args.vbc or args.vbc == 'True' or args.vbc > 0:
        print 'VBC Scan'
        scan_vbc_angle()
    exit()
