import os
import subprocess
from collections import OrderedDict
from SimulationFramework.FrameworkHelperFunctions import *
from SimulationFramework.Modules.merge_two_dicts import merge_two_dicts
from munch import Munch

astra_generator_keywords = {
    'keywords':{
        'filename': 'FName',
        'combine_distributions': 'add',
        'number_of_particles': 'Ipart',
        'probe_particle': 'probe',
        'noise_reduction': 'noise_reduc',
        'high_resolution': 'high_res',
        'charge': ['q_total', 1e9],
        'reference_position': 'ref_zpos',
        'reference_time': 'ref_clock',
        'distribution_type_z': 'dist_z',
        'inital_energy': 'ref_ekin',
        'plateau_bunch_length': ['lt', 1e9],
        'plateau_rise_time': ['rt', 1e9],
        'sigma_t': ['sig_clock', 1e9],
        'sigma_z': ['sig_z', 1e3],
        'bunch_length': ['lz', 1e3],
        'plateau_rise_distance': ['rz', 1e3],
        'distribution_type_pz': 'dist_pz',
        'energy_width': ['le', 1e-3],
        'distribution_type_x': 'dist_x',
        'sigma_x': ['sig_x', 1e3],
        'distribution_type_y': 'dist_y',
        'sigma_y': ['sig_y', 1e3],
        'distribution_type_px': 'dist_px',
        'distribution_type_py': 'dist_py',
        'normalized_horizontal_emittance': ['Nemit_x', 1e6],
        'normalized_vertical_emittance': ['Nemit_y', 1e6],
        'guassian_cutoff_x': 'C_sig_x',
        'guassian_cutoff_y': 'C_sig_y',
        'offset_x': ['x_off', 1e3],
        'offset_y': ['y_off', 1e3],
    },
}
gpt_generator_keywords = {
    'keywords':{
    },
}
generator_keywords = {
    'defaults': {
        'clara_400_3ps':{
            'combine_distributions': False,
            'species': 'electrons',
            'probe_particle': True,
            'noise_reduction': False,
            'high_resolution': True,
            'cathode': True,
            'reference_position': 0,
            'reference_time': 0,
            'distribution_type_z': 'p',
            'inital_energy': 0,
            'plateau_bunch_length': 3e-12,
            'plateau_rise_time': 0.2e-12,
            'distribution_type_pz': 'i',
            'energy_width': 0.62,
            'distribution_type_x': 'radial',
            'sigma_x': 0.25e-3,
            'distribution_type_y': 'r',
            'sigma_y': 0.25e-3,
            'offset_x': 0,
            'offset_y': 0,
        },
        'clara_400_1ps':{
            'combine_distributions': False,'species': 'electrons', 'probe_particle': True,'noise_reduction': False, 'high_resolution': True, 'cathode': True,
            'reference_position': 0, 'reference_time': 0, 'distribution_type_z': 'p',
            'inital_energy': 0, 'plateau_bunch_length': 1e-12, 'plateau_rise_time': 0.2e-12, 'distribution_type_pz': 'i', 'energy_width': 0.62,
            'distribution_type_x': 'radial', 'sigma_x': 0.25e-3, 'distribution_type_y': 'r', 'sigma_y': 0.25e-3,
            'offset_x': 0, 'offset_y': 0,
        },
        'clara_400_2ps_Gaussian':{
            'combine_distributions': False,'species': 'electrons', 'probe_particle': True,'noise_reduction': False, 'high_resolution': True, 'cathode': True,
            'reference_position': 0, 'reference_time': 0, 'distribution_type_z': 'g',
            'sigma_t': 0.85e-12,
            'inital_energy': 0, 'distribution_type_pz': 'i', 'energy_width': 0.62,
            'distribution_type_x': '2DGaussian', 'sigma_x': 0.25e-3, 'distribution_type_y': '2DGaussian', 'sigma_y': 0.25e-3,
            'guassian_cutoff_x': 3, 'guassian_cutoff_y': 3,
            'offset_x': 0, 'offset_y': 0,
        },
    },
    'keywords': [
        'number_of_particles', 'filename',
        'probe_particle', 'noise_reduction', 'high_resolution', 'combine_distributions',
        'cathode', 'cathode_radius',
        'charge', 'species',
        'emission_time', 'energy_width', 'bunch_length', 'inital_energy', 'energy_width',
        'sigma_x', 'sigma_y', 'sigma_z', 'sigma_t',
        'distribution_type_z', 'distribution_type_pz', 'distribution_type_x', 'distribution_type_px', 'distribution_type_y', 'distribution_type_py',
        'guassian_cutoff_x', 'guassian_cutoff_y',
        'plateau_bunch_length', 'plateau_rise_time', 'plateau_rise_distance',
        'offset_x', 'offset_y',
        'reference_position', 'reference_time',
        'normalized_horizontal_emittance', 'normalized_vertical_emittance',
        'image_file'
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
        if isinstance(defaults, str) and defaults in generator_keywords['defaults']:
            self.__init__(self.executables, **generator_keywords['defaults'][defaults])
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
        astra_keywords = list(astra_generator_keywords['keywords'].values())
        keywords = generator_keywords['keywords']
        self.allowedKeyWords = [*astra_keywords, *keywords]
        self.allowedKeyWords = [x.lower() if not isinstance(x, list) else x[0].lower() for x in self.allowedKeyWords]
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
            ['q_total', {'value': self.charge*1e9, 'default': 0.25}],
            ['Lprompt', {'value': False}],
        ])
        keyword_dict = OrderedDict()
        for k in self.allowedKeyWords:
            klower = k.lower()
            if klower not in [fk.lower() for fk in framework_dict.keys()]:
                if getattr(self, klower) is not None:
                    try:
                        val = eval(getattr(self, klower))
                    except:
                        val = getattr(self, klower)
                    if klower in astra_generator_keywords['keywords'].keys():
                        k = astra_generator_keywords['keywords'][klower]
                        if isinstance(k, list):
                            k, m = k
                            val = m * val
                    keyword_dict[k] = {'value': val}
                    print(k, val)
        output += self._write_ASTRA(merge_two_dicts(framework_dict, keyword_dict))
        output += '\n/\n'
        saveFile(self.global_parameters['master_subdir']+'/'+self.objectname+'.in', output)

    def astra_to_hdf5(self):
        astrabeamfilename = self.filename
        self.global_parameters['beam'].read_astra_beam_file(self.global_parameters['master_subdir'] + '/' + astrabeamfilename, normaliseZ=False)
        HDF5filename = self.filename.replace('.generator','.hdf5')
        self.global_parameters['beam'].write_HDF5_beam_file(self.global_parameters['master_subdir'] + '/' + HDF5filename, centered=False, sourcefilename=astrabeamfilename)

class GPTGenerator(frameworkGenerator):
    def __init__(self, executables, global_parameters, **kwargs):
        super(GPTGenerator, self).__init__(executables, global_parameters, **kwargs)
        gpt_keywords = list(gpt_generator_keywords['keywords'].values())
        keywords = generator_keywords['keywords']
        self.allowedKeyWords = [*gpt_keywords, *keywords]
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
        """Run the code with input 'filename'"""
        command = self.executables[self.code] + ['-o', self.objectname+'_out.gdf'] + ['GPTLICENSE='+self.global_parameters['GPTLICENSE']] + [self.objectname+'.in']
        my_env = os.environ.copy()
        my_env["LD_LIBRARY_PATH"] = my_env["LD_LIBRARY_PATH"] + ":/opt/GPT3.3.6/lib/" if "LD_LIBRARY_PATH" in my_env else "/opt/GPT3.3.6/lib/"
        my_env["OMP_WAIT_POLICY"] = "PASSIVE"
        # post_command_t = [self.executables[self.code][0].replace('gpt.exe','gdfa.exe')] + ['-o', self.objectname+'_emit.gdf'] + [self.objectname+'_out.gdf'] + ['time','avgx','avgy','stdx','stdBx','stdy','stdBy','stdz','stdt','nemixrms','nemiyrms','nemizrms','numpar','nemirrms','avgG','avgp','stdG','avgt','avgBx','avgBy','avgBz','CSalphax','CSalphay','CSbetax','CSbetay']
        post_command = [self.executables[self.code][0].replace('gpt','gdfa')] + ['-o', self.objectname+'_emit.gdf'] + [self.objectname+'_out.gdf'] + ['position','avgx','avgy','stdx','stdBx','stdy','stdBy','stdz','stdt','nemixrms','nemiyrms','nemizrms','numpar','nemirrms','avgG','avgp','stdG','avgt','avgBx','avgBy','avgBz','CSalphax','CSalphay','CSbetax','CSbetay']
        post_command_t = [self.executables[self.code][0].replace('gpt','gdfa')] + ['-o', self.objectname+'_emitt.gdf'] + [self.objectname+'_out.gdf'] + ['time','avgx','avgy','stdx','stdBx','stdy','stdBy','stdz','nemixrms','nemiyrms','nemizrms','numpar','nemirrms','avgG','avgp','stdG','avgBx','avgBy','avgBz','CSalphax','CSalphay','CSbetax','CSbetay']
        post_command_traj = [self.executables[self.code][0].replace('gpt','gdfa')] + ['-o', self.objectname+'traj.gdf'] + [self.objectname+'_out.gdf'] + ['time','avgx','avgy','avgz']
        with open(os.path.relpath(self.global_parameters['master_subdir']+'/'+self.objectname+'.log', '.'), "w") as f:
            # print('gpt command = ', command)
            subprocess.call(command, stdout=f, cwd=self.global_parameters['master_subdir'], env=my_env)
            subprocess.call(post_command, stdout=f, cwd=self.global_parameters['master_subdir'])
            subprocess.call(post_command_t, stdout=f, cwd=self.global_parameters['master_subdir'])
            subprocess.call(post_command_traj, stdout=f, cwd=self.global_parameters['master_subdir'])

    def generate_particles(self):
        return """#--Basic beam parameters--
                E0 = """ + self.energy_width + """;
                G = 1-qe*E0/(me*c*c);
                GB = sqrt(G^2-1);
                Qtot = """ + str(-1e12*self.charge) + """e-12;
                npart = """ + str(self.number_of_particles) + """;
                setparticles( "beam", nps, me, qe, Qtot ) ;
                """

    def check_xy_parameters(self, x: str, y: str, default: str):
        if getattr(self, x) is None and getattr(self, y) is not None:
            setattr(self, x, getattr(self, y))
        elif getattr(self, x) is not None and getattr(self, y) is None:
            setattr(self, y, getattr(self, x))
        elif getattr(self, x) is None and getattr(self, y) is None:
            setattr(self, x, default)
            setattr(self, y, default)

    def _uniform_distribution(self, distname: str):
        cutoff = self.guassian_cutoff_x if self.guassian_cutoff_x is not None else 3
        return distname + '( "beam", "g", 0, radius, 0, '+ str(cutoff) +') ;'

    def _gaussian_distribution(self, distname: str):
        cutoff = self.guassian_cutoff_x if self.guassian_cutoff_x is not None else 3
        return distname + '( "beam", "u", 0, radius, 0, '+ str(cutoff) +') ;'

    def _distribution(self, distname):
        if self.distribution_type_x.lower() in ["g","guassian"]:
            pass

    def generate_radial_distribution(self):
        self.check_xy_parameters("sigma_x", "sigma_y", 1)
        self.check_xy_parameters("distribution_type_x", "distribution_type_y", "g")
        if (self.sigma_x == self.sigma_y) and (self.distribution_type_x == self.distribution_type_y):
            output =  "radius = " + str(self.sigma_x) + ";"
            output += self._gaussian_distribution('setrxydist')
            setphidist("beam","u",0,2*pi) ;
            """
        else:
            return """radius_x = """ + str(self.sigma_x) + """;
            setxdist( "beam", " """+ str(self.distribution_type_x) +""" ", radius_x/2, radius_x ) ;
            setphidist("beam","u",0,2*pi) ;
            """

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
