import os
import subprocess
from collections import OrderedDict
from SimulationFramework.FrameworkHelperFunctions import *
from SimulationFramework.Modules.merge_two_dicts import merge_two_dicts

astra_generator_keywords = {
    'keywords':[
        'fname','add','ipart','species','probe','noise_reduc','high_res','cathode','lprompt', 'q_total','ref_zpos','ref_clock','dist_z','ref_ekin','lt','rt','sig_clock','sig_z','lz','rz',
        'dist_pz','le','dist_x','sig_x','dist_y','sig_y','dist_px','nemit', 'C_sig_x', 'C_sig_y', 'x_off', 'y_off',
    ],
    'defaults': {
        'clara_400_3ps':{
            'add': False,'species': 'electrons', 'probe': True,'noise_reduc': False, 'high_res': True, 'cathode': True, 'lprompt': False, 'ref_zpos': 0, 'ref_clock': 0, 'dist_z': 'p',
            'ref_ekin': 0, 'lt': 3e-3, 'rt': 0.2e-3, 'dist_pz': 'i', 'le': 0.62e-3, 'dist_x': 'radial', 'sig_x': 0.25, 'dist_y': 'r', 'sig_y': 0.25,
            'x_off': 0, 'y_off': 0,
        },
        'clara_400_1ps':{
            'add': False,'species': 'electrons', 'probe': True,'noise_reduc': False, 'high_res': True, 'cathode': True, 'lprompt': False, 'ref_zpos': 0, 'ref_clock': 0, 'dist_z': 'p',
            'ref_ekin': 0, 'lt': 1e-3, 'rt': 0.2e-3, 'dist_pz': 'i', 'le': 0.62e-3, 'dist_x': 'radial', 'sig_x': 0.25, 'dist_y': 'r', 'sig_y': 0.25,
            'x_off': 0, 'y_off': 0,
        },
        'clara_400_2ps_Gaussian':{
            'add': False,'species': 'electrons', 'probe': True,'noise_reduc': False, 'high_res': True, 'cathode': True, 'lprompt': False, 'ref_zpos': 0, 'ref_clock': 0, 'dist_z': 'g',
            'sig_clock': 0.85e-3,
            'ref_ekin': 0, 'dist_pz': 'i', 'le': 0.62e-3, 'dist_x': '2DGaussian', 'sig_x': 0.25, 'dist_y': '2DGaussian', 'sig_y': 0.25, 'C_sig_x': 3, 'C_sig_y': 3,
            'x_off': 0, 'y_off': 0,
        },
    },
    'framework_keywords': [
        'number_of_particles', 'charge', 'filename',
    ]
}

elegant_generator_keywords = {
    'keywords':[
        'bunch','n_particles_per_bunch', 'time_start', 'matched_to_cell', 'emit_x', 'emit_nx', 'beta_x', 'alpha_x', 'eta_x', 'etap_x', 'emit_y',
        'emit_ny', 'beta_y', 'alpha_y', 'eta_y', 'etap_y', 'use_twiss_command_values', 'use_moments_output_values', 'Po', 'sigma_dp','sigma_s',
        'dp_s_coupling', 'emit_z', 'beta_z', 'alpha_z', 'momentum_chirp', 'one_random_bunch', 'symmetrize', 'optimized_halton', 'limit_invariants',
        'limit_in_4d', 'first_is_fiducial', 'save_initial_coordinates', 'halton_sequence', 'halton_radix', 'randomize_order', 'enforce_rms_values',
        'distribution_cutoff', 'distribution_type', 'centroid'
    ],
    'defaults': {
    },
    'framework_keywords': [
        'number_of_particles', 'charge', 'filename',
    ]
}

class frameworkGenerator(object):
    def __init__(self, executables, global_parameters, **kwargs):
        super(frameworkGenerator, self).__init__()
        self.global_parameters = global_parameters
        self.executables = executables
        self.objectdefaults = {}
        self.objectproperties = {}

    def run(self):
        command = self.executables['generator'] + [self.objectname+'.in']
        with open(os.devnull, "w") as f:
            subprocess.call(command, stdout=f, cwd=self.global_parameters['master_subdir'])

    def load_defaults(self, defaults):
        if isinstance(defaults, str) and defaults in astra_generator_keywords['defaults']:
            self.__init__(self.executables, **astra_generator_keywords['defaults'][defaults])
        elif isinstance(defaults, dict):
            self.__init__(self.executables, **defaults)

    @property
    def particles(self):
        return self.number_of_particles if self.number_of_particles is not None else 512

    @particles.setter
    def particles(self, npart):
        self.add_property('number_of_particles', npart)

    @property
    def charge(self):
        return float(self.objectproperties['charge']) if 'charge' in self.objectproperties and self.objectproperties['charge'] is not None else 250e-12
    @charge.setter
    def charge(self, q):
        self.objectproperties['charge'] = q

    @property
    def objectname(self):
        return self.objectproperties['name'] if 'name' in self.objectproperties and self.objectproperties['name'] is not None else 'laser'

    def write(self):
        pass

    @property
    def parameters(self):
        return self.objectproperties

    def __getattr__(self, a):
        return None

    def add_property(self, key, value):
        if key.lower() in self.allowedKeyWords:
            self.objectproperties[key.lower()] = value
            self.__setattr__(key.lower(), value)

class ASTRAGenerator(frameworkGenerator):
    def __init__(self, executables, global_parameters, **kwargs):
        super(ASTRAGenerator, self).__init__(executables, global_parameters, **kwargs)
        self.allowedKeyWords = astra_generator_keywords['keywords'] + astra_generator_keywords['framework_keywords']
        self.allowedKeyWords = [x.lower() for x in self.allowedKeyWords]
        for key, value in list(kwargs.items()):
            key = key.lower()
            if key in self.allowedKeyWords:
                try:
                    # print 'key = ', key
                    self.objectproperties[key] = value
                    setattr(self, key, value)
                except:
                    pass
                    # print 'WARNING: Unknown keyword: ', key, value
                    # exit()

    def run(self):
        command = self.executables['ASTRAgenerator'] + [self.objectname+'.in']
        with open(os.devnull, "w") as f:
            subprocess.call(command, stdout=f, cwd=self.global_parameters['master_subdir'])

    def _write_ASTRA(self, d):
        output = ''
        for k, v in list(d.items()):
            val = v['value'] if v['value'] is not None else v['default'] if 'default' in v else None
            if isinstance(val,str):
                param_string = k+' = \''+str(val)+'\',\n'
            else:
                param_string = k+' = '+str(val)+',\n'
            if len((output + param_string).splitlines()[-1]) > 70:
                output += '\n'
            output += param_string
        return output[:-2]

    def write(self):
        output = '&INPUT\n'
        try:
            npart = eval(self.number_of_particles)
        except:
            npart = self.number_of_particles
        if self.filename is None:
            self.filename = 'laser.generator'
        framework_dict = OrderedDict([
            ['FName', {'value': self.filename, 'default': 'laser.generator'}],
            ['q_total', {'value': self.charge*1e9, 'default': 0.25}],
            ['Ipart', {'value': npart, 'default': 2**(3*3)}],
        ])
        keyword_dict = OrderedDict()
        for k in astra_generator_keywords['keywords']:
            k = k.lower()
            if getattr(self, k) is not None:
                try:
                    val = eval(getattr(self, k))
                except:
                    val = getattr(self, k)
                keyword_dict[k] = {'value': val}
        output += self._write_ASTRA(merge_two_dicts(framework_dict, keyword_dict))
        output += '\n/\n'
        saveFile(self.global_parameters['master_subdir']+'/'+self.objectname+'.in', output)

    def astra_to_hdf5(self):
        astrabeamfilename = self.filename
        self.global_parameters['beam'].read_astra_beam_file(self.global_parameters['master_subdir'] + '/' + astrabeamfilename, normaliseZ=False)
        HDF5filename = self.filename.replace('.generator','.hdf5')
        self.global_parameters['beam'].write_HDF5_beam_file(self.global_parameters['master_subdir'] + '/' + HDF5filename, centered=False, sourcefilename=astrabeamfilename)
