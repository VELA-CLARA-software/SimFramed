import numpy as np
import scipy.constants as constants
import sdds
import read_gdf_file as rgf

class twiss(object):

    E0 = constants.m_e * constants.speed_of_light**2
    E0_eV = E0 / constants.elementary_charge
    q_over_c = (constants.e / constants.speed_of_light)

    def __init__(self):
        self.twiss = {}


    def reset_dicts(self):
        self.twiss = {}

    def read_sdds_file(self, fileName, charge=None):
        self.reset_dicts()
        self.sdds = sdds.SDDS(0)
        self.sdds.load(fileName)
        for col in range(len(self.sdds.columnName)):
            if len(self.sdds.columnData[col]) == 1:
                self.twiss[self.sdds.columnName[col]] = np.array(self.sdds.columnData[col][0])
            else:
                self.twiss[self.sdds.columnName[col]] = np.array(self.sdds.columnData[col])
        self.SDDSparameterNames = list()
        for param in self.sdds.columnName:
            if isinstance(self.beam[param][0], (float, long)):
                self.SDDSparameterNames.append(param)

    def read_astra_emit_files(self, filename):
        self.reset_dicts()
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
        exn = 1e-6*eyn
        rms_x, rms_xp, rms_y, rms_yp, rms_z, rms_e = 1e-3*np.array([rms_x, rms_xp, rms_y, rms_yp, rms_z, rms_e])
        self.twiss['z'] = z
        self.twiss['t'] = t
        self.twiss['gamma'] = 1 + (e_kin / self.E0_eV)
        self.twiss['cp'] = np.sqrt(e_kin * (2 * self.E0_eV + e_kin))
        self.twiss['p'] = self.twiss['cp'] * self.q_over_c
        self.twiss['exn'] = exn
        self.twiss['ex'] = exn / self.twiss['gamma']
        self.twiss['eyn'] = eyn
        self.twiss['ey'] = eyn / self.twiss['gamma']
        self.twiss['ezn'] = ezn
        self.twiss['ez'] = ezn / self.twiss['gamma']
        self.twiss['beta_x'] = rms_x**2 / self.twiss['ex']
        self.twiss['gamma_x'] = rms_xp**2 / self.twiss['ex']
        self.twiss['alpha_x'] = (-1 * np.sign(mean_xxp) * rms_x * rms_xp) / self.twiss['ex']
        self.twiss['beta_y'] = rms_y**2 / self.twiss['ey']
        self.twiss['gamma_y'] = rms_yp**2 / self.twiss['ey']
        self.twiss['alpha_y'] = (-1 * np.sign(mean_yyp) * rms_y * rms_yp) / self.twiss['ey']
        self.twiss['beta_z'] = rms_z**2 / self.twiss['ez']
        self.twiss['gamma_z'] = rms_e**2 / self.twiss['ez']
        self.twiss['alpha_z'] = (-1 * np.sign(mean_zep) * rms_z * rms_e) / self.twiss['ez']
        self.twiss['sigma_x'] = rms_x
        self.twiss['sigma_y'] = rms_y
        self.twiss['sigma_z'] = rms_z
        self.twiss['sigma_p'] = (rms_e / e_kin) * self.twiss['p']
        self.twiss['sigma_cp'] = (rms_e / e_kin) * self.twiss['cp']

    def read_gdf_emit_file(self, file, position=None, time=None, block=None, charge=None):
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
