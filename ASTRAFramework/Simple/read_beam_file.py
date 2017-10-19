import numpy as np
import scipy.constants as constants
import sdds
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
        self.beam['z'] = (-1 * self.Bz * constants.speed_of_light) * self.t
        if charge is None:
            self.beam['total_charge'] = 0
        else:
            self.beam['total_charge'] = charge

    def set_beam_charge(self, charge):
        self.beam['total_charge'] = charge

    def read_astra_beam_file(self, file):
        self.reset_dicts()
        x, y, z, px, py, pz, clock, charge, index, status = np.loadtxt(file, unpack=True)
        z = self.normalise_to_ref_particle(z, subtractmean=True)
        pz = self.normalise_to_ref_particle(pz, subtractmean=False)
        p = np.sqrt(px**2 + py**2 + pz**2)
        self.beam['x'] = x
        self.beam['y'] = y
        self.beam['z'] = z
        self.beam['cpx'] = px
        self.beam['cpy'] = py
        self.beam['cpz'] = pz
        self.beam['px'] = px * self.q_over_c
        self.beam['py'] = py * self.q_over_c
        self.beam['pz'] = pz * self.q_over_c
        self.beam['cp'] = p
        self.beam['p'] = p * self.q_over_c
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
        self.beam['t'] = self.z / (-1 * self.Bz * constants.speed_of_light)
        self.beam['total_charge'] = np.sum(self.beam['charge'])

    def read_gdf_beam_file(self, file, position=None, time=None, block=None, charge=None):
        self.reset_dicts()
        gdfbeamdata = None
        gdfbeam = rgf.read_gdf_file(file)
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
            self.beam['t'] = self.z / (-1 * self.Bz * constants.speed_of_light)
        elif hasattr(gdfbeamdata,'t'):
            self.beam['t'] = gdfbeamdata.t
            self.beam['z'] = (-1 * self.Bz * constants.speed_of_light) * self.t
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

    def rms(self, u):
        return np.mean(u)

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
        nbins = int(np.ceil(twidth / self.slice_length))
        self._hist, binst =  np.histogram(self.t, bins=nbins)
        self.slice['t_Bins'] = binst
        self._t_binned = np.digitize(self.t, self.slice['t_Bins'])
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
        if not hasattr(self,'_tbins'):
            self.bin_time()
        self.slice['Momentum'] = np.array([cpbin.mean() for cpbin in self._cpbins])
        return self.slice['Momentum']
    @property
    def slice_momentum_spread(self):
        if not hasattr(self,'_tbins'):
            self.bin_time()
        self.slice['Momentum_Spread'] = np.array([cpbin.std() for cpbin in self._cpbins])
        return self.slice['Momentum_Spread']
    @property
    def slice_relative_momentum_spread(self):
        if not hasattr(self,'_tbins'):
            self.bin_time()
        self.slice['Relative_Momentum_Spread'] = np.array([100*cpbin.std()/cpbin.mean() for cpbin in self._cpbins])
        return self.slice['Relative_Momentum_Spread']
    @property
    def slice_normalized_horizontal_emittance(self):
        if not hasattr(self,'_tbins'):
            self.bin_time()
        xbins = [self.x[tbin] for tbin in self._tbins]
        xpbins = [self.xp[tbin] for tbin in self._tbins]
        emitbins = zip(*[xbins, xpbins, self._cpbins])
        self.slice['Normalized_Horizontal_Emittance'] = np.array([self.emittance(xbin, xpbin, cpbin) for xbin, xpbin, cpbin in emitbins])
        return self.slice['Normalized_Horizontal_Emittance']
    @property
    def slice_normalized_vertical_emittance(self):
        if not hasattr(self,'_tbins'):
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
    def x(self):
        return self.beam['x']
    @property
    def y(self):
        return self.beam['y']
    @property
    def z(self):
        return self.beam['z']
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
