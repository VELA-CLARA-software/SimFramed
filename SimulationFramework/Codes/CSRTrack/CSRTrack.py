import os
from ruamel import yaml
from SimulationFramework.Framework_objects import *
from SimulationFramework.Framework_elements import *
from SimulationFramework.FrameworkHelperFunctions import _rotation_matrix

with open(os.path.dirname( os.path.abspath(__file__))+'/csrtrack_defaults.yaml', 'r') as infile:
    csrtrack_defaults = yaml.load(infile, Loader=yaml.UnsafeLoader)

class csrtrackLattice(frameworkLattice):
    def __init__(self, *args, **kwargs):
        super(csrtrackLattice, self).__init__(*args, **kwargs)
        self.code = 'csrtrack'
        self.particle_definition = ''
        self.CSRTrackelementObjects = OrderedDict()
        self.set_particles_filename()

    def endScreen(self, **kwargs):
        return screen(name='end', type='screen', position_start=self.endObject.position_start, position_end=self.endObject.position_start, global_rotation=self.endObject.global_rotation, global_parameters=self.global_parameters, **kwargs)

    def set_particles_filename(self):
        self.CSRTrackelementObjects['particles'] = csrtrack_particles(particle_definition=self.particle_definition, global_parameters=self.global_parameters)
        self.CSRTrackelementObjects['particles'].format = 'astra'
        if self.particle_definition == 'initial_distribution':
            self.CSRTrackelementObjects['particles'].particle_definition = 'laser.astra'
            self.CSRTrackelementObjects['particles'].add_default('array', '#file{name=laser.astra}')
        else:
            self.CSRTrackelementObjects['particles'].particle_definition = self.allElementObjects[self.start].objectname
            self.CSRTrackelementObjects['particles'].add_default('array', '#file{name='+self.allElementObjects[self.start].objectname+'.astra'+'}')

    @property
    def dipoles_screens_and_bpms(self):
        return sorted(self.getElementType('dipole') + self.getElementType('screen') + self.getElementType('beam_position_monitor'), key=lambda x: x.position_end[2])

    def setCSRMode(self):
        if 'csr' in self.file_block and 'csr_mode' in self.file_block['csr']:
            if self.file_block['csr']['csr_mode'] == '3D':
                self.CSRTrackelementObjects['forces'] = csrtrack_forces(type='csr_g_to_p')
            elif self.file_block['csr']['csr_mode'] == '1D':
                self.CSRTrackelementObjects['forces'] = csrtrack_forces(type='projected')
        else:
            self.CSRTrackelementObjects['forces'] = csrtrack_forces()

    def writeElements(self):
        fulltext = 'io_path{logfile = log.txt}\nlattice{\n'
        counter = frameworkCounter(sub={'beam_position_monitor': 'screen'})
        for e in self.dipoles_screens_and_bpms:
            # if not e.type == 'dipole':
                # self.CSRTrackelementObjects[e.name] = csrtrack_online_monitor(filename=e.name+'.fmt2', monitor_type='phase', marker='screen'+str(counter.counter(e.type)), particle='all')
            fulltext += e.write_CSRTrack(counter.counter(e.objecttype))
            counter.add(e.objecttype)
        fulltext += self.endScreen().write_CSRTrack(counter.counter(self.endScreen().objecttype))
        fulltext += '}\n'
        self.set_particles_filename()
        self.setCSRMode()
        self.CSRTrackelementObjects['track_step'] = csrtrack_track_step()
        self.CSRTrackelementObjects['tracker'] = csrtrack_tracker(end_time_marker='screen'+str(counter.counter(self.endScreen().objecttype))+'a')
        self.CSRTrackelementObjects['monitor'] = csrtrack_monitor(name=self.end+'.fmt2', global_parameters=self.global_parameters)
        for c in self.CSRTrackelementObjects:
            fulltext += self.CSRTrackelementObjects[c].write_CSRTrack()
        return fulltext

    def write(self):
        self.code_file = self.global_parameters['master_subdir']+'/csrtrk.in'
        saveFile(self.code_file, self.writeElements())

    def preProcess(self):
        prefix = self.file_block['input']['prefix'] if 'input' in self.file_block and 'prefix' in self.file_block['input'] else ''
        self.CSRTrackelementObjects['particles'].hdf5_to_astra(prefix)

    def postProcess(self):
        self.CSRTrackelementObjects['monitor'].csrtrack_to_hdf5()

class csrtrack_element(frameworkElement):

    def __init__(self, elementName=None, elementType=None, **kwargs):
        super(csrtrack_element, self).__init__(elementName, elementType, **kwargs)
        self.header = ''
        if elementName in csrtrack_defaults:
            for k, v in list(csrtrack_defaults[elementName].items()):
                self.add_default(k, v)

    def CSRTrack_str(self, s):
        if s is True:
            return 'yes'
        elif s is False:
            return 'no'
        else:
            return str(s)

    def write_CSRTrack(self):
        output = str(self.header) + '{\n'
        for k in elementkeywords[self.objecttype]['keywords']:
            k = k.lower()
            if getattr(self,k) is not None:
                output += k+'='+self.CSRTrack_str(getattr(self, k))+'\n'
            elif k in self.objectdefaults:
                output += k+'='+self.CSRTrack_str(self.objectdefaults[k])+'\n'
        output+='}\n'
        return output

class csrtrack_online_monitor(csrtrack_element):

    def __init__(self, marker='', **kwargs):
        super(csrtrack_online_monitor, self).__init__('online_monitor', 'csrtrack_online_monitor', **kwargs)
        self.header = 'online_monitor'
        self.end_time_marker = marker+'b'

class csrtrack_forces(csrtrack_element):

    def __init__(self, **kwargs):
        super(csrtrack_forces, self).__init__('forces', 'csrtrack_forces', **kwargs)
        self.header = 'forces'

class csrtrack_track_step(csrtrack_element):

    def __init__(self, **kwargs):
        super(csrtrack_track_step, self).__init__('track_step', 'csrtrack_track_step', **kwargs)
        self.header = 'track_step'

class csrtrack_tracker(csrtrack_element):

    def __init__(self, end_time_marker='', **kwargs):
        super(csrtrack_tracker, self).__init__('tracker', 'csrtrack_tracker', **kwargs)
        self.header = 'tracker'
        self.end_time_marker = end_time_marker

class csrtrack_monitor(csrtrack_element):

    def __init__(self, **kwargs):
        super(csrtrack_monitor, self).__init__(elementName='monitor', elementType='csrtrack_monitor', **kwargs)
        self.header = 'monitor'

    def csrtrack_to_hdf5(self):
        csrtrackbeamfilename = self.name
        astrabeamfilename = csrtrackbeamfilename.replace('.fmt2','.astra')
        self.global_parameters['beam'].convert_csrtrackfile_to_astrafile(self.global_parameters['master_subdir'] + '/' + csrtrackbeamfilename, self.global_parameters['master_subdir'] + '/' + astrabeamfilename)
        self.global_parameters['beam'].read_astra_beam_file(self.global_parameters['master_subdir'] + '/' + astrabeamfilename, normaliseZ=False)
        HDF5filename = self.name.replace('.fmt2','.hdf5')
        self.global_parameters['beam'].write_HDF5_beam_file(self.global_parameters['master_subdir'] + '/' + HDF5filename, sourcefilename=csrtrackbeamfilename)

class csrtrack_particles(csrtrack_element):

    def __init__(self, **kwargs):
        super(csrtrack_particles, self).__init__('particles', 'csrtrack_particles', **kwargs)
        self.header = 'particles'

    def hdf5_to_astra(self, prefix=''):
        HDF5filename = prefix+self.particle_definition.replace('.astra','')+'.hdf5'
        self.global_parameters['beam'].read_HDF5_beam_file(self.global_parameters['master_subdir'] + '/' + HDF5filename)
        astrabeamfilename = self.particle_definition+'.astra'
        self.global_parameters['beam'].write_astra_beam_file(self.global_parameters['master_subdir'] + '/' + astrabeamfilename, normaliseZ=False)
