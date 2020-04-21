import os, math, h5py, sys
import numpy as np
from scipy import interpolate
import scipy.integrate as integrate
import scipy.constants as constants
import sdds
sys.path.append(os.path.dirname(__file__))
import read_gdf_file as rgf
import munch

class twiss(munch.Munch):

    E0 = constants.m_e * constants.speed_of_light**2
    E0_eV = E0 / constants.elementary_charge
    q_over_c = (constants.e / constants.speed_of_light)

    def __init__(self):
        super(twiss, self).__init__()
        self.reset_dicts()

    def find_nearest(self, array, value):
        idx = np.searchsorted(array, value, side="left")
        if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
            return idx-1
        else:
            return idx

    def reset_dicts(self):
        self.clear()
        self['z'] = []
        self['t'] = []
        self['kinetic_energy'] = []
        self['gamma'] = []
        self['cp'] = []
        self['cp_eV'] = []
        self['p'] = []
        self['enx'] = []
        self['ex'] = []
        self['eny'] = []
        self['ey'] = []
        self['enz'] = []
        self['ez'] = []
        self['beta_x'] = []
        self['gamma_x'] = []
        self['alpha_x'] = []
        self['beta_y'] = []
        self['gamma_y'] = []
        self['alpha_y'] = []
        self['beta_z'] = []
        self['gamma_z'] = []
        self['alpha_z'] = []
        self['sigma_x'] = []
        self['sigma_y'] = []
        self['sigma_z'] = []
        self['sigma_t'] = []
        self['sigma_p'] = []
        self['sigma_cp'] = []
        self['sigma_cp_eV'] = []
        self['mux'] = []
        self['muy'] = []
        self['eta_x'] = []
        self['eta_xp'] = []
        self['element_name'] = []
        self['x'] = []
        self['y'] = []
        self.elegant = {}

    def read_sdds_file(self, fileName, charge=None, ascii=False):
        # self.reset_dicts()
        self.sdds = sdds.SDDS(0)
        self.sdds.load(fileName)
        for col in range(len(self.sdds.columnName)):
            # print 'col = ', self.sdds.columnName[col]
            if len(self.sdds.columnData[col]) == 1:
                self.elegant[self.sdds.columnName[col]] = np.array(self.sdds.columnData[col][0])
            else:
                self.elegant[self.sdds.columnName[col]] = np.array(self.sdds.columnData[col])
        self.SDDSparameterNames = list()
        for i, param in enumerate(self.sdds.parameterName):
            # print 'param = ', param
            self.elegant[param] = self.sdds.parameterData[i]
            # if isinstance(self[param][0], (float, long)):
            #     self.SDDSparameterNames.append(param)

    def read_elegant_floor_file(self, filename, offset=[0,0,0], rotation=[0,0,0], reset=True):
        if reset:
            self.reset_dicts()
        self.read_sdds_file(filename)
        self['x'] = [np.round(x+offset[0], decimals = 6) for x in self.elegant['X']]
        self['y'] = [np.round(y+offset[1], decimals = 6) for y in self.elegant['Y']]
        self['z'] = [np.round(z+offset[2], decimals = 6) for z in self.elegant['Z']]
        self['theta'] = [np.round(theta+rotation[0], decimals = 6) for theta in self.elegant['theta']]
        self['phi'] = [np.round(phi+rotation[1], decimals = 6) for phi in self.elegant['phi']]
        self['psi'] = [np.round(psi+rotation[2], decimals = 6) for psi in self.elegant['psi']]
        xyz = list(zip(self['x'], self['y'], self['z']))
        thetaphipsi = list(zip(self['phi'], self['psi'], self['theta']))
        return list(zip(self.elegant['ElementName'], xyz[-1:] + xyz[:-1], xyz, thetaphipsi))[1:]

    def read_elegant_twiss_files(self, filename, startS=0, reset=True):
        if reset:
            self.reset_dicts()
        if isinstance(filename, (list, tuple)):
            for f in filename:
                self.read_elegant_twiss_files(f, reset=False)
        else:
            pre, ext = os.path.splitext(filename)
            self.read_sdds_file(pre+'.flr')
            self.read_sdds_file(pre+'.sig')
            self.read_sdds_file(pre+'.twi')
            z = self.elegant['Z']
            self['z'] = np.concatenate([self['z'], z])
            cp = self.elegant['pCentral0'] * self.E0
            self['cp'] = np.concatenate([self['cp'], cp])
            self['cp_eV'] = np.concatenate([self['cp_eV'], cp / constants.elementary_charge])
            ke = np.array((np.sqrt(self.E0**2 + cp**2) - self.E0**2) / constants.elementary_charge)
            self['kinetic_energy'] = np.concatenate([self['kinetic_energy'], ke])
            gamma = 1 + ke / self.E0_eV
            self['gamma'] = np.concatenate([self['gamma'], gamma])
            self['p'] = np.concatenate([self['p'], cp * self.q_over_c])
            self['enx'] = np.concatenate([self['enx'], self.elegant['enx']])
            self['ex'] = np.concatenate([self['ex'], self.elegant['ex']])
            self['eny'] = np.concatenate([self['eny'], self.elegant['eny']])
            self['ey'] = np.concatenate([self['ey'], self.elegant['ey']])
            self['enz'] = np.concatenate([self['enz'], np.zeros(len(self.elegant['Z']))])
            self['ez'] = np.concatenate([self['ez'], np.zeros(len(self.elegant['Z']))])
            self['beta_x'] = np.concatenate([self['beta_x'], self.elegant['betax']])
            self['alpha_x'] = np.concatenate([self['alpha_x'], self.elegant['alphax']])
            self['gamma_x'] = np.concatenate([self['gamma_x'], (1 + self.elegant['alphax']**2) / self.elegant['betax']])
            self['beta_y'] = np.concatenate([self['beta_y'], self.elegant['betay']])
            self['alpha_y'] = np.concatenate([self['alpha_y'], self.elegant['alphay']])
            self['gamma_y'] = np.concatenate([self['gamma_y'], (1 + self.elegant['alphay']**2) / self.elegant['betay']])
            self['beta_z'] = np.concatenate([self['beta_z'], np.zeros(len(self.elegant['Z']))])
            self['gamma_z'] = np.concatenate([self['gamma_z'], np.zeros(len(self.elegant['Z']))])
            self['alpha_z'] = np.concatenate([self['alpha_z'], np.zeros(len(self.elegant['Z']))])
            self['sigma_x'] = np.concatenate([self['sigma_x'], self.elegant['Sx']])
            self['sigma_y'] = np.concatenate([self['sigma_y'], self.elegant['Sy']])
            self['sigma_t'] = np.concatenate([self['sigma_t'], self.elegant['St']])
            beta = np.sqrt(1-(gamma**-2))
            # print 'len(z) = ', len(z), '  len(beta) = ', len(beta)
            self['t'] = np.concatenate([self['t'], z / (beta * constants.speed_of_light)])
            self['sigma_z'] =  np.concatenate([self['sigma_z'], self.elegant['St'] * (beta * constants.speed_of_light)])
            self['sigma_cp'] = np.concatenate([self['sigma_cp'], self.elegant['Sdelta'] * cp ])
            self['sigma_cp_eV'] = np.concatenate([self['sigma_cp_eV'], self.elegant['Sdelta'] * cp / constants.elementary_charge])
            # print('elegant = ', (self.elegant['Sdelta'] * cp / constants.elementary_charge)[-1])
            self['sigma_p'] = np.concatenate([self['sigma_p'], self.elegant['Sdelta'] ])
            self['mux'] = np.concatenate([self['mux'], self.elegant['psix'] / (2*constants.pi)])
            self['muy'] = np.concatenate([self['muy'], self.elegant['psiy'] / (2*constants.pi)])
            self['eta_x'] = np.concatenate([self['eta_x'], self.elegant['etax']])
            self['eta_xp'] = np.concatenate([self['eta_xp'], self.elegant['etaxp']])
            self['element_name'] = np.concatenate([self['element_name'], self.elegant['ElementName']])

    def read_astra_emit_files(self, filename, reset=True):
        if reset:
            self.reset_dicts()
        if isinstance(filename, (list, tuple)):
            for f in filename:
                self.read_astra_emit_files(f, reset=False)
        else:
            if 'xemit' not in filename.lower():
                filename = filename.replace('Yemit','Xemit').replace('Zemit','Xemit')
            xemit = np.loadtxt(filename, unpack=False)
            if 'yemit' not in filename.lower():
                filename = filename.replace('Xemit','Yemit').replace('Zemit','Yemit')
            yemit = np.loadtxt(filename, unpack=False)
            if 'zemit' not in filename.lower():
                filename = filename.replace('Xemit','Zemit').replace('Yemit','Zemit')
            zemit = np.loadtxt(filename, unpack=False)
            self.interpret_astra_data(xemit, yemit, zemit)

    def read_hdf_summary(self, filename, reset=True):
        if reset:
            self.reset_dicts()
        f = h5py.File(filename, "r")
        xemit = f.get('Xemit')
        yemit = f.get('Yemit')
        zemit = f.get('Zemit')
        for item, params in sorted(xemit.items()):
            self.interpret_astra_data(np.array(xemit.get(item)), np.array(yemit.get(item)), np.array(zemit.get(item)))

    def interpret_astra_data(self, xemit, yemit, zemit):
            z, t, mean_x, rms_x, rms_xp, exn, mean_xxp = np.transpose(xemit)
            z, t, mean_y, rms_y, rms_yp, eyn, mean_yyp = np.transpose(yemit)
            z, t, e_kin, rms_z, rms_e, ezn, mean_zep = np.transpose(zemit)
            e_kin = 1e6 * e_kin
            t = 1e-9 * t
            exn = 1e-6*exn
            eyn = 1e-6*eyn
            rms_x, rms_xp, rms_y, rms_yp, rms_z, rms_e = 1e-3*np.array([rms_x, rms_xp, rms_y, rms_yp, rms_z, rms_e])
            rms_e = 1e6 * rms_e
            self['z'] = np.concatenate([self['z'],z])
            self['t'] = np.concatenate([self['t'],t])
            self['kinetic_energy'] = np.concatenate([self['kinetic_energy'], e_kin])
            gamma = 1 + (e_kin / self.E0_eV)
            self['gamma'] = np.concatenate([self['gamma'], gamma])
            cp = np.sqrt(e_kin * (2 * self.E0_eV + e_kin)) * constants.elementary_charge
            self['cp'] = np.concatenate([self['cp'], cp])
            self['cp_eV'] = np.concatenate([self['cp_eV'], cp / constants.elementary_charge])
            p = cp * self.q_over_c
            self['p'] = np.concatenate([self['p'], p])
            self['enx'] = np.concatenate([self['enx'], exn])
            ex = exn / gamma
            self['ex'] = np.concatenate([self['ex'], ex])
            self['eny'] = np.concatenate([self['eny'], eyn])
            ey = eyn / gamma
            self['ey'] = np.concatenate([self['ey'], ey])
            self['enz'] = np.concatenate([self['enz'], ezn])
            ez = ezn / gamma
            self['ez'] = np.concatenate([self['ez'], ez])
            self['beta_x'] = np.concatenate([self['beta_x'], rms_x**2 / ex])
            self['gamma_x'] = np.concatenate([self['gamma_x'], rms_xp**2 / ex])
            self['alpha_x'] = np.concatenate([self['alpha_x'], (-1 * np.sign(mean_xxp) * rms_x * rms_xp) / ex])
            # self['alpha_x'] = np.concatenate([self['alpha_x'], (-1 * mean_xxp * rms_x) / ex])
            self['beta_y'] = np.concatenate([self['beta_y'], rms_y**2 / ey])
            self['gamma_y'] = np.concatenate([self['gamma_y'], rms_yp**2 / ey])
            self['alpha_y'] = np.concatenate([self['alpha_y'], (-1 * np.sign(mean_yyp) * rms_y * rms_yp) / ey])
            self['beta_z'] = np.concatenate([self['beta_z'], rms_z**2 / ez])
            self['gamma_z'] = np.concatenate([self['gamma_z'], rms_e**2 / ez])
            self['alpha_z'] = np.concatenate([self['alpha_z'], (-1 * np.sign(mean_zep) * rms_z * rms_e) / ez])
            self['sigma_x'] = np.concatenate([self['sigma_x'], rms_x])
            self['sigma_y'] = np.concatenate([self['sigma_y'], rms_y])
            self['sigma_z'] = np.concatenate([self['sigma_z'], rms_z])
            beta = np.sqrt(1-(gamma**-2))
            self['sigma_t'] = np.concatenate([self['sigma_t'], rms_z / (beta * constants.speed_of_light)])
            self['sigma_p'] = np.concatenate([self['sigma_p'], (rms_e / e_kin)])
            self['sigma_cp'] = np.concatenate([self['sigma_cp'], (rms_e / e_kin) * p])
            self['sigma_cp_eV'] = np.concatenate([self['sigma_cp_eV'], (rms_e)])
            # print('astra = ', (rms_e)[-1])
            self['mux'] = np.concatenate([self['mux'], integrate.cumtrapz(x=self['z'], y=1/self['beta_x'], initial=0)])
            self['muy'] = np.concatenate([self['muy'], integrate.cumtrapz(x=self['z'], y=1/self['beta_y'], initial=0)])

    def interpolate(self, z=None, value='z', index='z'):
        f = interpolate.interp1d(self[index], self[value], kind='linear')
        if z is None:
            return f
        else:
            if z > max(self[index]):
                return 10**6
            else:
                return float(f(z))

    def extract_values(self, array, start, end):
        startidx = self.find_nearest(self['z'], start)
        endidx = self.find_nearest(self['z'], end) + 1
        return self[array][startidx:endidx]

    def covariance(self, u, up):
        u2 = u - np.mean(u)
        up2 = up - np.mean(up)
        return np.mean(u2*up2) - np.mean(u2)*np.mean(up2)

    def write_HDF5_beam_file(self, filename, sourcefilename=None):
        with h5py.File(filename, "w") as f:
            inputgrp = f.create_group("Parameters")
            if sourcefilename is not None:
                inputgrp['Source'] = sourcefilename
            inputgrp['code'] = self.beam['code']
            twissgrp = f.create_group("twiss")
            array = np.array([self.s, self.t, self.sigma_x, self.sigma_y, self.sigma_z, self.sigma_p, self.sigma_t, self.beta_x, self.alpha_x, self.gamma_x,
                      self.beta_y, self.alpha_y, self.gamma_y, self.beta_z, self.alpha_z, self.gamma_z, self.mux, self.muy,
                      self.ex, self.enx, self.ey, self.eny]).transpose()
            beamgrp['columns'] = ("s","t","Sx","Sy","Sz","Sp","St","betax","alphax","gammax","betay","alphay","gammay","betaz","alphaz","gammaz","mux","muy")
            beamgrp['units'] = ("m","s","m","m","m","eV/c","s","m","","","m","","","m","","","","")
            beamgrp.create_dataset("twiss", data=array)

    def read_HDF5_beam_file(self, filename, local=False):
        self.reset_dicts()
        with h5py.File(filename, "r") as h5file:
            if h5file.get('beam/reference_particle') is not None:
                self.beam['reference_particle'] = np.array(h5file.get('beam/reference_particle'))
            if h5file.get('beam/longitudinal_reference') is not None:
                self.beam['longitudinal_reference'] = np.array(h5file.get('beam/longitudinal_reference'))
            else:
                self.beam['longitudinal_reference'] = 't'
            if h5file.get('beam/status') is not None:
                self.beam['status'] = np.array(h5file.get('beam/status'))
            x, y, z, cpx, cpy, cpz, t, charge = np.array(h5file.get('beam/beam')).transpose()
            cp = np.sqrt(cpx**2 + cpy**2 + cpz**2)
            self.beam['x'] = x
            self.beam['y'] = y
            self.beam['z'] = z
            # self.beam['cpx'] = cpx
            # self.beam['cpy'] = cpy
            # self.beam['cpz'] = cpz
            self.beam['px'] = cpx * self.q_over_c
            self.beam['py'] = cpy * self.q_over_c
            self.beam['pz'] = cpz * self.q_over_c
            # self.beam['cp'] = cp
            # self.beam['p'] = cp * self.q_over_c
            # self.beam['xp'] = np.arctan(self.px/self.pz)
            # self.beam['yp'] = np.arctan(self.py/self.pz)
            self.beam['clock'] = np.full(len(self.x), 0)
            # self.beam['gamma'] = np.sqrt(1+(self.cp/self.E0_eV)**2)
            # velocity_conversion = 1 / (constants.m_e * self.gamma)
            # self.beam['vx'] = velocity_conversion * self.px
            # self.beam['vy'] = velocity_conversion * self.py
            # self.beam['vz'] = velocity_conversion * self.pz
            # self.beam['Bx'] = self.vx / constants.speed_of_light
            # self.beam['By'] = self.vy / constants.speed_of_light
            # self.beam['Bz'] = self.vz / constants.speed_of_light
            self.beam['t'] = t
            self.beam['charge'] = charge
            self.beam['total_charge'] = np.sum(self.beam['charge'])
            startposition = np.array(h5file.get('/Parameters/Starting_Position'))
            startposition = startposition if startposition is not None else [0,0,0]
            self.beam['starting_position'] = startposition
            theta =  np.array(h5file.get('/Parameters/Rotation'))
            theta = theta if theta is not None else 0
            self.beam['rotation'] = theta
            if local == True:
                self.rotate_beamXZ(self.beam['rotation'], preOffset=self.beam['starting_position'])
