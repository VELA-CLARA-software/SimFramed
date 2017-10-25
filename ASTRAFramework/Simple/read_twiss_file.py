import os
import numpy as np
import scipy.constants as constants
if os.name == 'nt':
    import sdds
import read_gdf_file as rgf

class twiss(dict):

    E0 = constants.m_e * constants.speed_of_light**2
    E0_eV = E0 / constants.elementary_charge
    q_over_c = (constants.e / constants.speed_of_light)

    def reset_dicts(self):
        self['z'] = []
        self['t'] = []
        self['gamma'] = []
        self['cp'] = []
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
        self['sigma_p'] = []
        self['sigma_cp'] = []

    def read_sdds_file(self, fileName, charge=None):
        self.reset_dicts()
        self.sdds = sdds.SDDS(0)
        self.sdds.load(fileName)
        for col in range(len(self.sdds.columnName)):
            if len(self.sdds.columnData[col]) == 1:
                self[self.sdds.columnName[col]] = np.array(self.sdds.columnData[col][0])
            else:
                self[self.sdds.columnName[col]] = np.array(self.sdds.columnData[col])
        self.SDDSparameterNames = list()
        for param in self.sdds.columnName:
            if isinstance(self.beam[param][0], (float, long)):
                self.SDDSparameterNames.append(param)

    def read_astra_emit_files(self, filename, reset=True):
        if reset:
            self.reset_dicts()
        if isinstance(filename, (list, tuple)):
            for f in filename:
                self.read_astra_emit_files(f, reset=False)
        else:
            if 'xemit' not in filename.lower():
                filename = filename.replace('Yemit','Xemit').replace('Zemit','Xemit')
            z, t, mean_x, rms_x, rms_xp, exn, mean_xxp = np.loadtxt(filename, unpack=True)
            if 'yemit' not in filename.lower():
                filename = filename.replace('Xemit','Yemit').replace('Zemit','Yemit')
            z, t, mean_y, rms_y, rms_yp, eyn, mean_yyp = np.loadtxt(filename, unpack=True)
            if 'zemit' not in filename.lower():
                filename = filename.replace('Xemit','Zemit').replace('Yemit','Zemit')
            z, t, e_kin, rms_z, rms_e, ezn, mean_zep = np.loadtxt(filename, unpack=True)
            e_kin = 1e6 * e_kin
            t = 1e-9 * t
            exn = 1e-6*exn
            eyn = 1e-6*eyn
            rms_x, rms_xp, rms_y, rms_yp, rms_z, rms_e = 1e-3*np.array([rms_x, rms_xp, rms_y, rms_yp, rms_z, rms_e])
            self['z'] = np.concatenate([self['z'],z])
            self['t'] = np.concatenate([self['t'],t])
            gamma = 1 + (e_kin / self.E0_eV)
            self['gamma'] = np.concatenate([self['gamma'], gamma])
            cp = np.sqrt(e_kin * (2 * self.E0_eV + e_kin))
            self['cp'] = np.concatenate([self['cp'], cp])
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
            self['sigma_p'] = np.concatenate([self['sigma_p'], (rms_e / e_kin) * p])
            self['sigma_cp'] = np.concatenate([self['sigma_cp'], (rms_e / e_kin) * p])

    def covariance(self, u, up):
        u2 = u - np.mean(u)
        up2 = up - np.mean(up)
        return np.mean(u2*up2) - np.mean(u2)*np.mean(up2)
