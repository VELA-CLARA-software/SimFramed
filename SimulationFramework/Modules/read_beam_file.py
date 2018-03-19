import os, time, csv
import numpy as np
import scipy.constants as constants
from scipy.spatial.distance import cdist
if os.name == 'nt':
    try:
        import sdds
    except:
        pass
import read_gdf_file as rgf

class beam(object):

    E0 = constants.m_e * constants.speed_of_light**2
    E0_eV = E0 / constants.elementary_charge
    q_over_c = (constants.e / constants.speed_of_light)

    def __init__(self):
        self.beam = {}

    def normalise_to_ref_particle(self, array, index=0,subtractmean=False):
        array[1:] = array[0] + array[1:]
        if subtractmean:
            array = array - np.mean(array)
        return array

    def reset_dicts(self):
        self.beam = {}
        self.twiss = {}
        self.slice = {}
        self._tbins = []
        self._pbins = []

    def read_sdds_file(self, fileName, charge=None):
        self.reset_dicts()
        self.sdds = sdds.SDDS(0)
        self.sdds.load(fileName)
        for col in range(len(self.sdds.columnName)):
            if len(self.sdds.columnData[col]) == 1:
                self.beam[self.sdds.columnName[col]] = np.array(self.sdds.columnData[col][0])
            else:
                self.beam[self.sdds.columnName[col]] = np.array(self.sdds.columnData[col])
        self.SDDSparameterNames = list()
        for param in self.sdds.columnName:
            if isinstance(self.beam[param][0], (float, long)):
                self.SDDSparameterNames.append(param)
        self.beam['cp'] = self.beam['p'] * self.E0_eV
        self.beam['cpz'] = np.sqrt(self.cp / (self.xp**2 + self.yp**2 + 1))
        self.beam['cpx'] = self.xp * self.beam['cpz']
        self.beam['cpy'] = self.yp * self.beam['cpz']
        self.beam['px'] = self.beam['cpx'] * self.q_over_c
        self.beam['py'] = self.beam['cpy'] * self.q_over_c
        self.beam['pz'] = self.beam['cpz'] * self.q_over_c
        self.beam['p'] = self.beam['cp'] * self.q_over_c
        self.beam['gamma'] = np.sqrt(1+(self.cp/self.E0_eV)**2)
        velocity_conversion = 1 / (constants.m_e * self.gamma)
        self.beam['vx'] = velocity_conversion * self.px
        self.beam['vy'] = velocity_conversion * self.py
        self.beam['vz'] = velocity_conversion * self.pz
        self.beam['Bx'] = self.vx / constants.speed_of_light
        self.beam['By'] = self.vy / constants.speed_of_light
        self.beam['Bz'] = self.vz / constants.speed_of_light
        self.beam['z'] = np.full(len(self.t), 0) #(-1 * self.Bz * constants.speed_of_light) * self.t
        if charge is None:
            self.beam['total_charge'] = 0
        else:
            self.beam['total_charge'] = charge
        self.beam['charge'] = []

    def set_beam_charge(self, charge):
        self.beam['total_charge'] = charge

    def read_astra_beam_file(self, file):
        starttime = time.time()
        self.reset_dicts()
        with open(file, 'r') as f:
            data = np.array([l for l in csv.reader(f, delimiter=' ',  quoting=csv.QUOTE_NONNUMERIC, skipinitialspace=True)])
        # datanp = np.loadtxt(file)
        # data = pd.read_csv(file, delim_whitespace=True, quoting=csv.QUOTE_NONNUMERIC).values
        self.interpret_astra_data(data)

    def read_hdf5_beam(self, data):
        self.reset_dicts()
        self.interpret_astra_data(data)

    def interpret_astra_data(self, data):
        x, y, z, cpx, cpy, cpz, clock, charge, index, status = np.transpose(data)
        z = self.normalise_to_ref_particle(z, subtractmean=False)
        cpz = self.normalise_to_ref_particle(cpz, subtractmean=False)
        clock = self.normalise_to_ref_particle(clock, subtractmean=True)
        cp = np.sqrt(cpx**2 + cpy**2 + cpz**2)
        self.beam['x'] = x
        self.beam['y'] = y
        self.beam['z'] = z
        self.beam['cpx'] = cpx
        self.beam['cpy'] = cpy
        self.beam['cpz'] = cpz
        self.beam['px'] = cpx * self.q_over_c
        self.beam['py'] = cpy * self.q_over_c
        self.beam['pz'] = cpz * self.q_over_c
        self.beam['cp'] = cp
        self.beam['p'] = cp * self.q_over_c
        self.beam['xp'] = np.arctan(self.px/self.pz)
        self.beam['yp'] = np.arctan(self.py/self.pz)
        self.beam['clock'] = clock
        self.beam['charge'] = 1e-9*charge
        self.beam['index'] = index
        self.beam['status'] = status
        self.beam['gamma'] = np.sqrt(1+(self.cp/self.E0_eV)**2)
        velocity_conversion = 1 / (constants.m_e * self.gamma)
        self.beam['vx'] = velocity_conversion * self.px
        self.beam['vy'] = velocity_conversion * self.py
        self.beam['vz'] = velocity_conversion * self.pz
        self.beam['Bx'] = self.vx / constants.speed_of_light
        self.beam['By'] = self.vy / constants.speed_of_light
        self.beam['Bz'] = self.vz / constants.speed_of_light
        self.beam['t'] = self.z / (-1 * self.Bz * constants.speed_of_light)#[time if status is -1 else 0 for time, status in zip(clock, status)]#
        self.beam['total_charge'] = np.sum(self.beam['charge'])

    def read_csrtrack_beam_file(self, file):
        self.reset_dicts()
        data = np.loadtxt(file, unpack=False)
        z, x, y, cpz, cpx, cpy, charge = np.transpose(data[1:])
        z = self.normalise_to_ref_particle(z, subtractmean=False)
        cpz = self.normalise_to_ref_particle(cpz, subtractmean=False)
        clock = self.normalise_to_ref_particle(clock, subtractmean=True)
        cp = np.sqrt(cpx**2 + cpy**2 + cpz**2)
        self.beam['x'] = x
        self.beam['y'] = y
        self.beam['z'] = z
        self.beam['cpx'] = cpx
        self.beam['cpy'] = cpy
        self.beam['cpz'] = cpz
        self.beam['px'] = cpx * self.q_over_c
        self.beam['py'] = cpy * self.q_over_c
        self.beam['pz'] = cpz * self.q_over_c
        self.beam['cp'] = cp
        self.beam['p'] = cp * self.q_over_c
        self.beam['xp'] = np.arctan(self.px/self.pz)
        self.beam['yp'] = np.arctan(self.py/self.pz)
        self.beam['clock'] = np.full(len(self.x), 0)
        self.beam['clock'][0] = data[0, 0] * 1e-9
        # self.beam['charge'] = 1e-9*charge
        self.beam['index'] = np.full(len(self.x), 5)
        self.beam['status'] = np.full(len(self.x), 1)
        self.beam['gamma'] = np.sqrt(1+(self.cp/self.E0_eV)**2)
        velocity_conversion = 1 / (constants.m_e * self.gamma)
        self.beam['vx'] = velocity_conversion * self.px
        self.beam['vy'] = velocity_conversion * self.py
        self.beam['vz'] = velocity_conversion * self.pz
        self.beam['Bx'] = self.vx / constants.speed_of_light
        self.beam['By'] = self.vy / constants.speed_of_light
        self.beam['Bz'] = self.vz / constants.speed_of_light
        self.beam['t'] = self.z / (-1 * self.Bz * constants.speed_of_light)# [time if status is -1 else 0 for time, status in zip(clock, self.beam['status'])]
        self.beam['charge'] = charge
        self.beam['total_charge'] = np.sum(self.beam['charge'])


    def convert_csrtrackfile_to_astrafile(self, csrtrackfile, astrafile):
        data = np.loadtxt(csrtrackfile, unpack=False)
        z, x, y, cpz, cpx, cpy, charge = np.transpose(data[1:])
        charge = -charge*1e9
        clock0 = (data[0, 0] / constants.speed_of_light) * 1e9
        # print 'clock0 = ', clock0
        clock = np.full(len(x), 0)
        clock[0] = clock0
        index = np.full(len(x), 1)
        status = np.full(len(x), 5)
        array = np.array([x, y, z, cpx, cpy, cpz, clock, charge, index, status]).transpose()
        np.savetxt(astrafile, array, fmt=('%.12e','%.12e','%.12e','%.12e','%.12e','%.12e','%.12e','%.12e','%.12e','%d'))

    def find_nearest_vector(self, nodes, node):
        return cdist([node], nodes).argmin()

    def rms(self, x, axis=None):
        return np.sqrt(np.mean(x**2, axis=axis))

    def create_ref_particle(self, array, index=0, subtractmean=False):
        array[1:] = array[0] + array[1:]
        if subtractmean:
            array = array - np.mean(array)
        return array

    def write_astra_beam_file(self, file, index=1, status=5, charge=None, normaliseZ=False):
        if not isinstance(index,(list, tuple, np.ndarray)):
            if len(self.beam['charge']) == len(self.x):
                chargevector = self.beam['charge']
            else:
                chargevector = np.full(len(self.x), self.charge/len(self.x))
        if not isinstance(index,(list, tuple, np.ndarray)):
            indexvector = np.full(len(self.x), index)
        if not isinstance(status,(list, tuple, np.ndarray)):
            statusvector = np.full(len(self.x), status)
        ''' if a particle is emitting from the cathode it's z value is 0 and it's clock value is finite, otherwise z is finite and clock is irrelevant (thus zero) '''
        zvector = [0 if status == -1 and t == 0 else z for status, z, t in zip(statusvector, self.z, self.t)]
        ''' if the clock value is finite, we calculate it from the z value, using Betaz '''
        clockvector = [z / (-1 * Bz * constants.speed_of_light) if status == -1 and t == 0 else t for status, z, t, Bz in zip(statusvector, self.z, self.t, self.Bz)]
        ''' this is the ASTRA array in all it's glory '''
        array = np.array([self.x, self.y, zvector, self.cpx, self.cpy, self.cpz, clockvector, chargevector, indexvector, statusvector]).transpose()
        ''' take the rms - if the rms is 0 set it to 1, so we don't get a divide by error '''
        rms_vector = [a if abs(a) > 0 else 1 for a in self.rms(array, axis=0)]
        ''' normalise the array '''
        norm_array = array / rms_vector
        ''' take the meen of the normalised array '''
        mean_vector = np.mean(norm_array, axis=0)
        ''' find the index of the vector that is closest to the mean - if you read in an ASTRA file, this should actually return the reference particle! '''
        nearest_idx = self.find_nearest_vector(norm_array, mean_vector)
        ref_particle = array[nearest_idx]
        ''' set the closest mean vector to be in position 0 in the array '''
        array = np.roll(array, -1*nearest_idx, axis=0)
        ''' normalise Z to the reference particle '''
        array[1:,2] = array[1:,2] - ref_particle[2]
        ''' should we leave Z as the reference value, set it to 0, or set it to be some offset? '''
        if not normaliseZ is False:
            array[0,2] = 0
        if not isinstance(normaliseZ,(bool)):
            array[0,2] += normaliseZ
        ''' normalise pz and the clock '''
        array[1:,5] = array[1:,5] - ref_particle[5]
        array[1:,6] = array[1:,6] - ref_particle[6]
        np.savetxt(file, array, fmt=('%.12e','%.12e','%.12e','%.12e','%.12e','%.12e','%.12e','%.12e','%d','%d'))

    def write_vsim_beam_file(self, file, normaliseT=False):
        if len(self.beam['charge']) == len(self.x):
            chargevector = self.beam['charge']
        else:
            chargevector = np.full(len(self.x), self.beam['total_charge']/len(self.x))
        if normaliseT:
            tvector = self.t - np.mean(self.t)
            zvector = self.z - np.mean(self.z)
        else:
            tvector = self.t
            zvector = self.z
        zvector = [t * (-1 * Bz * constants.speed_of_light) if z == 0 else z for z, t, Bz in zip(zvector, tvector, self.Bz)]
        ''' this is the VSIM array in all it's glory '''
        array = np.array([zvector, self.y, self.x, self.Bz*self.gamma*constants.speed_of_light, self.By*self.gamma*constants.speed_of_light, self.Bx*self.gamma*constants.speed_of_light]).transpose()
        ''' take the rms - if the rms is 0 set it to 1, so we don't get a divide by error '''
        np.savetxt(file, array, fmt=('%.12e','%.12e','%.12e','%.12e','%.12e','%.12e'))


    def read_gdf_beam_file_object(self, file):
        if isinstance(file, (str)):
            gdfbeam = rgf.read_gdf_file(file)
        elif isinstance(file, (rgf.read_gdf_file)):
            gdfbeam = file
        else:
            raise Exception('file is not str or gdf object!')
        return gdfbeam

    def read_gdf_beam_file(self, file, position=None, time=None, block=None, charge=None):
        self.reset_dicts()
        gdfbeamdata = None
        gdfbeam = self.read_gdf_beam_file_object(file)

        if position is not None and (time is not None or block is not None):
            print 'Assuming position over time!'
            gdfbeamdata = gdfbeam.get_position(position)
            if gdfbeamdata is not None:
                time = None
                block = None
            else:
                 position = None
        if position is None and time is not None and block is not None:
            print 'Assuming time over block!'
            gdfbeamdata = gdfbeam.get_time(time)
            if gdfbeamdata is not None:
                block = None
            else:
                 time = None
        if position is None and time is None and block is not None:
            gdfbeamdata = gdfbeam.get_grab(block)
            if gdfbeamdata is None:
                block = None
        if position is None and time is None and block is None:
            gdfbeamdata = gdfbeam.get_grab(0)
        self.beam['x'] = gdfbeamdata.x
        self.beam['y'] = gdfbeamdata.y
        self.beam['Bx'] = gdfbeamdata.Bx
        self.beam['By'] = gdfbeamdata.By
        self.beam['Bz'] = gdfbeamdata.Bz
        if hasattr(gdfbeamdata,'z'):
            self.beam['z'] = gdfbeamdata.z
            self.beam['t'] = np.full(len(self.z), 0)# self.z / (-1 * self.Bz * constants.speed_of_light)
        elif hasattr(gdfbeamdata,'t'):
            self.beam['t'] = gdfbeamdata.t
            self.beam['z'] = np.full(len(self.t), 0)#(-1 * self.Bz * constants.speed_of_light) * self.t
        self.beam['gamma'] = gdfbeamdata.G
        if hasattr(gdfbeamdata,'q'):
            self.beam['charge'] = gdfbeamdata.q
            self.beam['total_charge'] = np.sum(self.beam['charge'])
        else:
            if charge is None:
                self.beam['total_charge'] = 0
            else:
                self.beam['total_charge'] = charge
        self.beam['vx'] = self.Bx * constants.speed_of_light
        self.beam['vy'] = self.By * constants.speed_of_light
        self.beam['vz'] = self.Bz * constants.speed_of_light
        velocity_conversion = 1 / (constants.m_e * self.gamma)
        self.beam['px'] = self.beam['vx'] / velocity_conversion
        self.beam['py'] = self.beam['vy'] / velocity_conversion
        self.beam['pz'] = self.beam['vz'] / velocity_conversion
        self.beam['p'] = np.sqrt(self.px**2 + self.py**2 + self.pz**2)
        self.beam['cpx'] = self.px / self.q_over_c
        self.beam['cpy'] = self.py / self.q_over_c
        self.beam['cpz'] = self.pz / self.q_over_c
        self.beam['cp'] = self.p / self.q_over_c
        self.beam['xp'] = np.arctan(self.px/self.pz)
        self.beam['yp'] = np.arctan(self.py/self.pz)

    def write_HDF5_beam_file(self, filein, fileout):
        savescreens = [str(s) for s in screens]
        screenpositions = self.getScreenFiles()
        # print 'screenpositions = ', list(screenpositions.iteritems())
        if reference is not None:
            filename = '_'.join(map(str,[reference, self.settingsFile,self.globalSettings['charge'],self.globalSettings['npart']])) + '.hdf5'
        else:
            filename = '_'.join(map(str,[self.settingsFile,self.globalSettings['charge'],self.globalSettings['npart']])) + '.hdf5'
        # print filename
        f = h5py.File(filename, "w")
        inputgrp = f.create_group("Input")
        inputgrp['charge'] = self.globalSettings['charge']
        inputgrp['npart'] = self.globalSettings['npart']
        inputgrp['subdirectory'] = self.subdirectory
        xemitgrp = f.create_group("Xemit")
        yemitgrp = f.create_group("Yemit")
        zemitgrp = f.create_group("Zemit")
        screengrp = f.create_group("screens")
        if os.path.isfile(self.subdir+'/'+self.globalSettings['initial_distribution']):
            inputgrp.create_dataset('initial_distribution',data=numpy.loadtxt(self.subdir+'/'+self.globalSettings['initial_distribution']))
        for n, screendict in sorted(screenpositions.iteritems()):
            if os.path.isfile(self.subdir+'/'+n+'.in'):
                inputfile = file(self.subdir+'/'+n+'.in','r')
                inputfilecontents = inputfile.read()
                inputgrp.create_dataset(n, data=inputfilecontents)
            for emit, grp in {'X': xemitgrp,'Y': yemitgrp,'Z': zemitgrp}.iteritems():
                emitfile = self.subdir+'/'+n+'.'+emit+'emit.'+screendict['run']
                if os.path.isfile(emitfile):
                    grp.create_dataset(n, data=numpy.loadtxt(emitfile))
            for s in screendict['screenpositions']:
                screenfile = self.subdir+'/'+n+'.'+s+'.'+screendict['run']
                if screens == [] or s in savescreens:
                    # print 's = ', screenfile
                    screengrp.create_dataset(s, data=numpy.loadtxt(screenfile))
        csrtrackfiles = glob.glob(self.subdir+'/csrtrk.in')
        for csrfile in csrtrackfiles:
            inputfile = file(csrfile,'r')
            inputfilecontents = inputfile.read()
            inputgrp.create_dataset('csrtrack', data=inputfilecontents)
            screenfile = self.subdir+'/end.fmt2.astra'
            if os.path.isfile(screenfile):
                screengrp.create_dataset(screenfile, data=numpy.loadtxt(screenfile))










    def covariance(self, u, up):
        u2 = u - np.mean(u)
        up2 = up - np.mean(up)
        return np.mean(u2*up2) - np.mean(u2)*np.mean(up2)

    def emittance(self, x, xp, p=None):
        emittance = np.sqrt(self.covariance(x, x)*self.covariance(xp, xp) - self.covariance(x, xp)**2)
        if p is None:
            return emittance
        else:
            gamma = np.mean(p)/self.E0_eV
            return gamma*emittance

    @property
    def normalized_horizontal_emittance(self):
        return self.emittance(self.x, self.xp, self.cp)
    @property
    def normalized_vertical_emittance(self):
        return self.emittance(self.y, self.yp, self.cp)
    @property
    def horizontal_emittance(self):
        return self.emittance(self.x, self.xp)
    @property
    def vertical_emittance(self):
        return self.emittance(self.y, self.yp)
    @property
    def horizontal_emittance_90(self):
        emit = self.horizontal_emittance
        alpha = self.alpha_x
        beta = self.beta_x
        gamma = self.gamma_x
        emiti = gamma * self.x**2 + 2 * alpha * self.x * self.xp + beta * self.xp * self.xp
        return sorted(emiti)[int(0.9*len(emiti)-0.5)]
    @property
    def normalized_horizontal_emittance_90(self):
        emit = self.horizontal_emittance_90
        return np.mean(self.cp)/self.E0_eV * emit
    @property
    def vertical_emittance_90(self):
        emit = self.vertical_emittance
        alpha = self.alpha_y
        beta = self.beta_y
        gamma = self.gamma_y
        emiti = gamma * self.y**2 + 2 * alpha * self.y * self.yp + beta * self.yp * self.yp
        return sorted(emiti)[int(0.9*len(emiti)-0.5)]
    @property
    def normalized_vertical_emittance_90(self):
        emit = self.vertical_emittance_90
        return np.mean(self.cp)/self.E0_eV * emit

    @property
    def beta_x(self):
        self.twiss['beta_x'] = self.covariance(self.x,self.x) / self.horizontal_emittance
        return self.twiss['beta_x']
    @property
    def alpha_x(self):
        self.twiss['alpha_x'] = -1*self.covariance(self.x,self.xp) / self.horizontal_emittance
        return self.twiss['alpha_x']
    @property
    def gamma_x(self):
        self.twiss['gamma_x'] = self.covariance(self.xp,self.xp) / self.horizontal_emittance
        return self.twiss['gamma_x']
    @property
    def beta_y(self):
        self.twiss['beta_y'] = self.covariance(self.y,self.y) / self.vertical_emittance
        return self.twiss['beta_y']
    @property
    def alpha_y(self):
        self.twiss['alpha_y'] = -1*self.covariance(self.y,self.yp) / self.vertical_emittance
        return self.twiss['alpha_y']
    @property
    def gamma_y(self):
        self.twiss['gamma_y'] = self.covariance(self.yp,self.yp) / self.vertical_emittance
        return self.twiss['gamma_y']
    @property
    def twiss_analysis(self):
        return self.horizontal_emittance, self.alpha_x, self.beta_x, self.gamma_x, self.vertical_emittance, self.alpha_y, self.beta_y, self.gamma_y

    def eta_correlation(self, u):
        return self.covariance(u,self.p) / self.covariance(self.p, self.p)
    def eta_corrected(self, u):
        return u - self.eta_correlation(u)*self.p
    @property
    def horizontal_emittance_corrected(self):
        xc = self.eta_corrected(self.x)
        xpc = self.eta_corrected(self.xp)
        return self.emittance(xc, xpc)
    @property
    def vertical_emittance_corrected(self):
        yc = self.eta_corrected(self.y)
        ypc = self.eta_corrected(self.yp)
        return self.emittance(yc, ypc)
    @property
    def beta_x_corrected(self):
        xc = self.eta_corrected(self.x)
        self.twiss['beta_x'] = self.covariance(xc, xc) / self.horizontal_emittance_corrected
        return self.twiss['beta_x']
    @property
    def alpha_x_corrected(self):
        xc = self.eta_corrected(self.x)
        xpc = self.eta_corrected(self.xp)
        self.twiss['alpha_x'] = -1*self.covariance(xc, xpc) / self.horizontal_emittance_corrected
        return self.twiss['alpha_x']
    @property
    def gamma_x_corrected(self):
        xpc = self.eta_corrected(self.xp)
        self.twiss['gamma_x'] = self.covariance(xpc, xpc) / self.horizontal_emittance_corrected
        return self.twiss['gamma_x']
    @property
    def beta_y_corrected(self):
        yc = self.eta_corrected(self.y)
        self.twiss['beta_y'] = self.covariance(yc,yc) / self.vertical_emittance_corrected
        return self.twiss['beta_y']
    @property
    def alpha_y_corrected(self):
        yc = self.eta_corrected(self.y)
        ypc = self.eta_corrected(self.yp)
        self.twiss['alpha_y'] = -1*self.covariance(yc, ypc) / self.vertical_emittance_corrected
        return self.twiss['alpha_y']
    @property
    def gamma_y_corrected(self):
        ypc = self.eta_corrected(self.yp)
        self.twiss['gamma_y'] = self.covariance(ypc,ypc) / self.vertical_emittance_corrected
        return self.twiss['gamma_y']
    @property
    def twiss_analysis_corrected(self):
        return  self.horizontal_emittance_corrected, self.alpha_x_corrected, self.beta_x_corrected, self.gamma_x_corrected, \
                self.vertical_emittance_corrected, self.alpha_y_corrected, self.beta_y_corrected, self.gamma_y_corrected

    @property
    def slice_length(self):
        return self._slicelength

    @slice_length.setter
    def slice_length(self, slicelength):
        self._slicelength = slicelength

    def bin_time(self):
        if not hasattr(self,'slice'):
            self.slice = {}
        if not hasattr(self,'_slicelength'):
            self.slice_length = 0.1e-12
            print("Assuming slice length is 100 fs")
        twidth = (max(self.t) - min(self.t))
        if twidth == 0:
            t = self.z / (-1 * self.Bz * constants.speed_of_light)
            twidth = (max(t) - min(t))
        else:
            t = self.t
        if not self.slice_length > 0.0:
            self.slice_length = twidth / 20.0
        nbins = max([1,int(np.ceil(twidth / self.slice_length))])
        self._hist, binst =  np.histogram(t, bins=nbins)
        self.slice['t_Bins'] = binst
        self._t_binned = np.digitize(t, self.slice['t_Bins'])
        self._tbins = [[self._t_binned == i] for i in range(1, len(binst))]
        self._cpbins = [self.cp[tbin] for tbin in self._tbins]

    @property
    def slice_bins(self):
        if not hasattr(self,'slice'):
            self.bin_time()
        bins = self.slice['t_Bins']
        return (bins[:-1] + bins[1:]) / 2
        # return [t.mean() for t in ]
    @property
    def slice_momentum(self):
        if not hasattr(self,'_tbins') or not hasattr(self,'_cpbins'):
            self.bin_time()
        self.slice['Momentum'] = np.array([cpbin.mean() for cpbin in self._cpbins])
        return self.slice['Momentum']
    @property
    def slice_momentum_spread(self):
        if not hasattr(self,'_tbins') or not hasattr(self,'_cpbins'):
            self.bin_time()
        self.slice['Momentum_Spread'] = np.array([cpbin.std() for cpbin in self._cpbins])
        return self.slice['Momentum_Spread']
    @property
    def slice_relative_momentum_spread(self):
        if not hasattr(self,'_tbins') or not hasattr(self,'_cpbins'):
            self.bin_time()
        self.slice['Relative_Momentum_Spread'] = np.array([100*cpbin.std()/cpbin.mean() for cpbin in self._cpbins])
        return self.slice['Relative_Momentum_Spread']
    @property
    def slice_normalized_horizontal_emittance(self):
        if not hasattr(self,'_tbins') or not hasattr(self,'_cpbins'):
            self.bin_time()
        xbins = [self.x[tbin] for tbin in self._tbins]
        xpbins = [self.xp[tbin] for tbin in self._tbins]
        emitbins = zip(*[xbins, xpbins, self._cpbins])
        self.slice['Normalized_Horizontal_Emittance'] = np.array([self.emittance(xbin, xpbin, cpbin) for xbin, xpbin, cpbin in emitbins])
        return self.slice['Normalized_Horizontal_Emittance']
    @property
    def slice_normalized_vertical_emittance(self):
        if not hasattr(self,'_tbins') or not hasattr(self,'_cpbins'):
            self.bin_time()
        ybins = [self.y[tbin] for tbin in self._tbins]
        ypbins = [self.yp[tbin] for tbin in self._tbins]
        emitbins = zip(*[ybins, ypbins, self._cpbins])
        self.slice['Normalized_Vertical_Emittance'] = np.array([self.emittance(ybin, ypbin, cpbin) for ybin, ypbin, cpbin in emitbins])
        return self.slice['Normalized_Vertical_Emittance']
    @property
    def slice_peak_current(self):
        if not hasattr(self,'_hist'):
            self.bin_time()
        self.slice['Peak_Current'] = self.charge/(self.slice_length * len(self.t)) * self._hist
        return abs(self.slice['Peak_Current'])
    @property
    def slice_max_peak_current_slice(self):
        peakI = self.slice_peak_current
        self.slice['Max_Peak_Current_Slice'] = list(abs(peakI)).index(max(abs(peakI)))
        return self.slice['Max_Peak_Current_Slice']

    def sliceAnalysis(self):
        self.slice = {}
        self.bin_time()
        peakIPosition = self.slice_max_peak_current_slice
        return self.slice_peak_current[peakIPosition], \
            self.slice_relative_momentum_spread[peakIPosition], \
            self.slice_normalized_horizontal_emittance[peakIPosition], \
            self.slice_normalized_vertical_emittance[peakIPosition], \
            self.slice_momentum[peakIPosition]

    @property
    def chirp(self):
        self.bin_time()
        slice_current_centroid_indices = []
        slice_momentum_centroid = []
        peakIPosition = self.slice_max_peak_current_slice
        peakI = self.slice_peak_current[peakIPosition]
        slicemomentum = self.slice_momentum
        for index, slice_current in enumerate(self.slice_peak_current):
            if abs(peakI - slice_current) < (peakI * 0.75):
                slice_current_centroid_indices.append(index)
        for index in slice_current_centroid_indices:
            slice_momentum_centroid.append(slicemomentum[index])
        chirp = (1e-18 * (slice_momentum_centroid[-1] - slice_momentum_centroid[0]) / (len(slice_momentum_centroid) * self.slice_length))
        return chirp

    @property
    def x(self):
        return self.beam['x']
    @property
    def y(self):
        return self.beam['y']
    @property
    def z(self):
        return self.beam['z']
    @property
    def zn(self):
        return self.beam['z']-np.mean(self.beam['z'])
    @property
    def px(self):
        return self.beam['px']
    @property
    def py(self):
        return self.beam['py']
    @property
    def pz(self):
        return self.beam['pz']
    @property
    def cpx(self):
        return self.beam['cpx']
    @property
    def cpy(self):
        return self.beam['cpy']
    @property
    def cpz(self):
        return self.beam['cpz']
    @property
    def xp(self):
        return self.beam['xp']
    @property
    def yp(self):
        return self.beam['yp']
    @property
    def t(self):
        return self.beam['t']
    @property
    def p(self):
        return self.beam['p']
    @property
    def cp(self):
        return self.beam['cp']
    @property
    def gamma(self):
        return self.beam['gamma']
    @property
    def BetaGamma(self):
        return self.beam['cp']/self.E0_eV
    @property
    def vx(self):
        return self.beam['vx']
    @property
    def vy(self):
        return self.beam['vy']
    @property
    def vz(self):
        return self.beam['vz']
    @property
    def Bx(self):
        return self.beam['Bx']
    @property
    def By(self):
        return self.beam['By']
    @property
    def Bz(self):
        return self.beam['Bz']
    @property
    def charge(self):
        return self.beam['total_charge']
