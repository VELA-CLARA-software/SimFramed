from SimulationFramework.Codes.ASTRA.ASTRARules import *
from SimulationFramework.Framework_objects import *
from SimulationFramework.Framework_elements import *

section_header_text_ASTRA = {'cavities': {'header': 'CAVITY', 'bool': 'LEField'},
                             'wakefields': {'header': 'WAKE', 'bool': 'LWAKE'},
                             'solenoids': {'header': 'SOLENOID', 'bool': 'LBField'},
                             'quadrupoles': {'header': 'QUADRUPOLE', 'bool': 'LQuad'},
                             'dipoles': {'header': 'DIPOLE', 'bool': 'LDipole'},
                             'astra_newrun': {'header': 'NEWRUN'},
                             'astra_output': {'header': 'OUTPUT'},
                             'astra_charge': {'header': 'CHARGE'},
                             'global_error': {'header': 'ERROR'},
                             'apertures': {'header': 'APERTURE', 'bool': 'LApert'},
                            }


class astraLattice(frameworkLattice):
    def __init__(self, *args, **kwargs):
        super(astraLattice, self).__init__(*args, **kwargs)
        self.code = 'astra'
        self.bunch_charge = None
        self.headers = OrderedDict()
        self.starting_offset = eval(expand_substitution(self, self.file_block['starting_offset'])) if 'starting_offset' in self.file_block else [0,0,0]

        # This calculated the starting rotation based on the input file and the number of dipoles
        self.starting_rotation = -1*self.allElementObjects[self.start].global_rotation[2] if self.allElementObjects[self.start].global_rotation is not None else 0
        # print 'self.starting_rotation = ', self.starting_rotation
        # Calculate the correct starting offset by adding up the dipole angles
        for d in self.dipoles:
            self.starting_rotation -= d.angle
        # print 'self.starting_rotation after subtraction = ', self.starting_rotation
        self.starting_rotation = eval(expand_substitution(self, str(self.file_block['starting_rotation']))) if 'starting_rotation' in self.file_block else self.starting_rotation
        # print 'self.starting_rotation at end = ', self.starting_rotation

        # Create a "newrun" block
        if 'input' not in self.file_block:
            self.file_block['input'] = {}
        self.headers['newrun'] = astra_newrun(self.starting_offset, self.starting_rotation, global_parameters=self.global_parameters, **merge_two_dicts(self.file_block['input'],self.globalSettings['ASTRAsettings']))
        # If the initial distribution is derived from a generator file, we should use that
        if self.headers['newrun']['particle_definition'] == 'initial_distribution':
            self.headers['newrun']['particle_definition'] = 'laser.astra'
        else:
            self.headers['newrun']['particle_definition'] = self.allElementObjects[self.start].objectname+'.astra'

        # Create an "output" block
        if 'output' not in self.file_block:
            self.file_block['output'] = {}
        self.headers['output'] = astra_output(self.screens_and_bpms, self.starting_offset, self.starting_rotation, global_parameters=self.global_parameters, **merge_two_dicts(self.file_block['output'],self.globalSettings['ASTRAsettings']))

        # Create a "charge" block
        if 'charge' not in self.file_block:
            self.file_block['charge'] = {}
        self.headers['charge'] = astra_charge(global_parameters=self.global_parameters, **merge_two_dicts(self.file_block['charge'],self.globalSettings['ASTRAsettings']))

        # Create an "error" block
        if 'global_errors' not in self.file_block:
            self.file_block['global_errors'] = {}
        if 'global_errors' not in self.globalSettings:
            self.globalSettings['global_errors'] = {}
        if 'global_errors' in self.file_block or 'global_errors' in self.globalSettings:
            self.global_error = global_error(name=self.objectname+'_global_error', global_parameters=self.global_parameters)
            self.headers['global_errors'] = astra_errors(element=self.global_error, global_parameters=self.global_parameters, **merge_two_dicts(self.file_block['global_errors'], self.globalSettings['global_errors']))
        # print 'errors = ', self.file_block, self.headers['global_errors']

    @property
    def sample_interval(self):
        return self._sample_interval
    @sample_interval.setter
    def sample_interval(self, interval):
        # print('Setting new ASTRA sample_interval = ', interval)
        self._sample_interval = interval
        self.headers['newrun'].sample_interval = interval
        self.headers['charge'].sample_interval = interval

    def writeElements(self):
        fulltext = ''
        # Create objects for the newrun, output and charge blocks
        self.headers['output'].start_element = self.allElementObjects[self.start]
        self.headers['output'].end_element = self.allElementObjects[self.end]
        self.headers['output'].screens = self.screens_and_bpms
        # write the headers and their elements
        for header in self.headers:
            fulltext += self.headers[header].write_ASTRA()+'\n'
        # Initialise a counter object
        counter = frameworkCounter(sub={'kicker': 'dipole', 'collimator': 'aperture'})
        for t in [['cavities'], ['wakefields'], ['solenoids'], ['quadrupoles'], ['dipoles', 'dipoles_and_kickers'], ['apertures']]:
            fulltext += '&' + section_header_text_ASTRA[t[0]]['header']+'\n'
            elements = getattr(self, t[-1])
            fulltext += section_header_text_ASTRA[t[0]]['bool']+' = '+str(len(elements) > 0)+'\n'
            for element in elements:
                element.starting_offset = self.starting_offset
                element.starting_rotation = self.starting_rotation
                elemstr = element.write_ASTRA(counter.counter(element.objecttype))
                if elemstr is not None and not elemstr == '':
                    fulltext += elemstr+'\n'
                    if element.objecttype == 'kicker':
                        counter.add(element.objecttype,2)
                    elif element.objecttype == 'longitudinal_wakefield':
                        counter.add(element.objecttype, element.cells)
                    elif element.objecttype == 'aperture' or element.objecttype == 'collimator':
                        counter.add(element.objecttype, element.number_of_elements)
                    else:
                        counter.add(element.objecttype)
            fulltext += '\n/\n'
        return fulltext

    def write(self):
        self.code_file = self.global_parameters['master_subdir']+'/'+self.objectname+'.in'
        saveFile(self.code_file, self.writeElements())

    def preProcess(self):
        prefix = self.file_block['input']['prefix'] if 'input' in self.file_block and 'prefix' in self.file_block['input'] else ''
        self.headers['newrun'].hdf5_to_astra(prefix)
        self.headers['charge'].npart = len(self.global_parameters['beam'].x)

    def postProcess(self):
        for e in self.screens_and_bpms:
            if not self.starting_offset == [0,0,0]:
                e.zstart = self.allElementObjects[self.start].start
            else:
                e.zstart = [0,0,0]
            e.astra_to_hdf5(self.objectname)
        self.astra_to_hdf5()

    def astra_to_hdf5(self):
        master_run_no = self.global_parameters['run_no'] if 'run_no' in self.global_parameters else 1
        if not self.starting_offset == [0,0,0]:
            zstart = self.allElementObjects[self.start].start
        else:
            zstart = [0,0,0]
        startpos = np.array(self.allElementObjects[self.start].start)-np.array(zstart)
        endpos = np.array(self.allElementObjects[self.end].end)-np.array(zstart)
        astrabeamfilename = self.objectname + '.' + str(int(round(endpos[2]*100))).zfill(4) + '.' + str(master_run_no).zfill(3)
        if not os.path.isfile(self.global_parameters['master_subdir'] + '/' + astrabeamfilename):
            # print('Can\'t find ASTRA beam file: ', astrabeamfilename)
            astrabeamfilename = self.objectname + '.' + str(int(round((endpos[2]-startpos[2])*1000))).zfill(4) + '.' + str(master_run_no).zfill(3)
            # print('Trying relative naming convention: ', astrabeamfilename)
        self.global_parameters['beam'].read_astra_beam_file(self.global_parameters['master_subdir'] + '/' + astrabeamfilename, normaliseZ=False)
        self.global_parameters['beam'].rotate_beamXZ(-1*self.starting_rotation, preOffset=[0,0,0], postOffset=-1*np.array(self.starting_offset))
        HDF5filename = self.allElementObjects[self.end].objectname+'.hdf5'
        self.global_parameters['beam'].write_HDF5_beam_file(self.global_parameters['master_subdir'] + '/' + HDF5filename, centered=False, sourcefilename=astrabeamfilename, pos=self.allElementObjects[self.end].middle)

class astra_header(frameworkElement):

    def __init__(self, name=None, type=None, **kwargs):
        super(astra_header, self).__init__(name, type, **kwargs)

    def framework_dict(self):
        return OrderedDict()

    def write_ASTRA(self):
        keyword_dict = OrderedDict()
        for k in elementkeywords[self.objecttype]['keywords']:
            if getattr(self, k.lower()) is not None:
                keyword_dict[k.lower()] = {'value': getattr(self,k.lower())}
        output = '&' + section_header_text_ASTRA[self.objecttype]['header']+'\n'
        output += self._write_ASTRA(merge_two_dicts(self.framework_dict(), keyword_dict), None) + '\n/\n'
        return output

class astra_newrun(astra_header):
    def __init__(self, offset, rotation, **kwargs):
        super(astra_header, self).__init__('newrun', 'astra_newrun', **kwargs)
        self.starting_offset = offset
        self.starting_rotation = rotation
        self.sample_interval = 1
        if 'run' not in kwargs:
            self.run = 1
        if 'head' not in kwargs:
            self.head = 'trial'
        if 'lprompt' not in kwargs:
            self.add_property('lprompt', False)

    def framework_dict(self):
        return OrderedDict([
            ['Distribution', {'value': '\''+self.particle_definition+'\''}],
            ['high_res', {'value': self.high_res, 'default': True}],
            ['n_red', {'value': self.sample_interval, 'default': 1}],
            ['auto_phase', {'value': self.auto_phase, 'default': True}]
        ])

    def hdf5_to_astra(self, prefix=''):
        HDF5filename = prefix+self.particle_definition.replace('.astra','')+'.hdf5'
        self.global_parameters['beam'].read_HDF5_beam_file(self.global_parameters['master_subdir'] + '/' + HDF5filename)
        self.global_parameters['beam'].rotate_beamXZ(self.starting_rotation, preOffset=self.starting_offset)
        astrabeamfilename = self.particle_definition
        self.global_parameters['beam'].write_astra_beam_file(self.global_parameters['master_subdir'] + '/' + astrabeamfilename, normaliseZ=False)

class astra_output(astra_header):
    def __init__(self, screens, offset, rotation, **kwargs):
        super(astra_header, self).__init__('output', 'astra_output', **kwargs)
        self.screens = screens
        self.starting_offset = offset
        self.starting_rotation = rotation

    def framework_dict(self):
        self.start_element.starting_offset = self.starting_offset
        self.end_element.starting_offset = self.starting_offset
        self.start_element.starting_rotation = self.starting_rotation
        self.end_element.starting_rotation = self.starting_rotation
        # print self.end_element.objectname, self.end_element.end, self.start_element.objectname, self.start_element.end
        keyworddict = OrderedDict([
            ['zstart', {'value': self.start_element.start[2]}],
            ['zstop', {'value': self.end_element.end[2]}],
        ])
        for i, element in enumerate(self.screens,1):
            element.starting_offset = self.starting_offset
            element.starting_rotation = self.starting_rotation
            keyworddict['Screen('+str(i)+')'] = {'value': element.middle[2]}
            # if abs(element.theta) > 0:
                # keyworddict['Scr_xrot('+str(i)+')'] = {'value': element.theta}
        return keyworddict

class astra_charge(astra_header):
    def __init__(self, **kwargs):
        super(astra_header, self).__init__('charge', 'astra_charge', **kwargs)
        self.npart = 2**(3*5)
        self.sample_interval = 1
        self.grids = getGrids()

    @property
    def space_charge(self):
        return False if self.space_charge_mode == 'False' or self.space_charge_mode == False else True

    @property
    def space_charge_2D(self):
        return True if self.space_charge and self.space_charge_mode != '3D' else False

    @property
    def space_charge_3D(self):
        return True if self.space_charge and not self.space_charge_2D else False

    @property
    def grid_size(self):
        # print 'asking for grid sizes n = ', self.npart, ' is ', self.grids.getGridSizes(self.npart)
        return self.grids.getGridSizes((self.npart/self.sample_interval))

    def framework_dict(self):
        sc_dict = OrderedDict([
            ['Lmirror', {'value': self.cathode, 'default': False}],
            ['LSPCH', {'value': self.space_charge, 'default': True}],
            ['LSPCH3D', {'value': self.space_charge_3D, 'default': True}]
        ])
        if self.space_charge_2D:
            sc_n_dict = OrderedDict([
                ['nrad', {'value': self.grid_size, 'default': 32}],
                ['nlong_in', {'value': self.grid_size, 'default': 32}],
            ])
        elif self.space_charge_3D:
            sc_n_dict = OrderedDict([
                ['nxf', {'value': self.grid_size, 'default': 8}],
                ['nyf', {'value': self.grid_size, 'default': 8}],
                ['nzf', {'value': 2*self.grid_size, 'default': 8}],
            ])
        else:
            sc_n_dict = OrderedDict([])
        return merge_two_dicts(sc_dict, sc_n_dict)

class astra_errors(astra_header):
    def __init__(self, element=None, **kwargs):
        super(astra_errors, self).__init__('astra_error', 'global_error', **kwargs)
        self._element = element
        self.add_property('global_errors', True)
        self.add_property('Log_Error', True)
        self.add_property('generate_output', True)
        self.add_property('Suppress_output', False)

    def write_ASTRA(self):
        keyword_dict = {}
        conversion = dict([a, b] for a, b in zip(elementkeywords[self.objecttype]['keywords'], elementkeywords[self.objecttype]['astra_keywords']))
        for k in elementkeywords[self.objecttype]['keywords']:
            # print k, conversion[k]
            if getattr(self, k.lower()) is not None:
                keyword_dict[conversion[k].lower()] = {'value': getattr(self,k.lower())}
        joineddict = merge_two_dicts(self.framework_dict(), keyword_dict)
        if len(joineddict) > 0:
            output = '&' + section_header_text_ASTRA[self.objecttype]['header']+'\n'
            output += self._write_ASTRA(merge_two_dicts(self.framework_dict(), keyword_dict), None) + '\n/\n'
        else:
            output = ''
        return output
