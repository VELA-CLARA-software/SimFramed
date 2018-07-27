import time, os, subprocess, re
import yaml
import traceback
import itertools
import copy
from collections import OrderedDict
from ASTRARules import *
from SimulationFramework.FrameworkHelperFunctions import *
import SimulationFramework.Modules.read_beam_file as rbf
beam = rbf.beam()
from operator import add
import progressbar

_mapping_tag = yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG

def dict_representer(dumper, data):
    return dumper.represent_dict(data.iteritems())

def dict_constructor(loader, node):
    return OrderedDict(loader.construct_pairs(node))

yaml.add_representer(OrderedDict, dict_representer)
yaml.add_constructor(_mapping_tag, dict_constructor)


astra_generator_keywords = {
    'keywords':[
        'fname','add','ipart','species','probe','noise_reduc','high_res','cathode','lprompt', 'q_total','ref_zpos','ref_clock','dist_z','ref_ekin','lt','rt','sig_clock','sig_z','lz','rz',
        'dist_pz','le','dist_x','sig_x','dist_y','sig_y','dist_px','nemit',
    ],
    'defaults': {
        'clara_400_3ps':{
            'add': False,'species': 'electrons', 'probe': True,'noise_reduc': False, 'high_res': True, 'cathode': True, 'lprompt': False, 'ref_zpos': 0, 'ref_clock': 0, 'dist_z': 'p',
            'ref_ekin': 0, 'lt': 3e-3, 'rt': 0.2e-3, 'dist_pz': 'i', 'le': 0.62e-3, 'dist_x': 'radial', 'sig_x': 0.25, 'dist_y': 'r', 'sig_y': 0.25,
        },
        'clara_400_1ps':{
            'add': False,'species': 'electrons', 'probe': True,'noise_reduc': False, 'high_res': True, 'cathode': True, 'lprompt': False, 'ref_zpos': 0, 'ref_clock': 0, 'dist_z': 'p',
            'ref_ekin': 0, 'lt': 1e-3, 'rt': 0.2e-3, 'dist_pz': 'i', 'le': 0.62e-3, 'dist_x': 'radial', 'sig_x': 0.25, 'dist_y': 'r', 'sig_y': 0.25,
        }
    },
    'framework_keywords': [
        'number_of_particles', 'charge', 'filename',
    ]
}

with open(os.path.dirname( os.path.abspath(__file__))+'/commands_Elegant.yaml', 'r') as infile:
    commandkeywords = yaml.load(infile)

with open(os.path.dirname( os.path.abspath(__file__))+'/csrtrack_defaults.yaml', 'r') as infile:
    csrtrack_defaults = yaml.load(infile)

with open(os.path.dirname( os.path.abspath(__file__))+'/elementkeywords.yaml', 'r') as infile:
    elementkeywords = yaml.load(infile)

with open(os.path.dirname( os.path.abspath(__file__))+'/elements_Elegant.yaml', 'r') as infile:
    elements_Elegant = yaml.load(infile)

type_conversion_rules_Elegant = {'dipole': 'csrcsbend', 'quadrupole': 'kquad', 'beam_position_monitor': 'moni', 'screen': 'watch', 'aperture': 'rcol',
                         'collimator': 'ecol', 'monitor': 'moni', 'solenoid': 'mapsolenoid', 'wall_current_monitor': 'charge', 'cavity': 'rfcw',
                         'rf_deflecting_cavity': 'rfdf'}

keyword_conversion_rules_Elegant = {'length': 'l','entrance_edge_angle': 'e1', 'exit_edge_angle': 'e2', 'edge_field_integral': 'fint', 'horizontal_size': 'x_max', 'vertical_size': 'y_max',
                            'field_amplitude': 'volt', 'frequency': 'freq', 'output_filename': 'filename', 'csr_bins': 'bins'}

section_header_text_ASTRA = {'cavities': {'header': 'CAVITY', 'bool': 'LEField'},
                             'solenoids': {'header': 'SOLENOID', 'bool': 'LBField'},
                             'quadrupoles': {'header': 'QUADRUPOLE', 'bool': 'LQuad'},
                             'dipoles': {'header': 'DIPOLE', 'bool': 'LDipole'},
                             'astra_newrun': {'header': 'NEWRUN'},
                             'astra_output': {'header': 'OUTPUT'},
                             'astra_charge': {'header': 'CHARGE'},
                            }

def isevaluable(self, s):
    try:
        eval(s)
        return True
    except:
        return False

def expand_substitution(self, param, subs={}):
    if isinstance(param,(str)):
        regex = re.compile('\$(.*)\$')
        s = re.search(regex, param)
        if s:
            if isevaluable(self, s.group(1)) is True:
                replaced_str = str(eval(re.sub(regex, str(eval(s.group(1))), param)))
            else:
                replaced_str = re.sub(regex, s.group(1), param)
            for key in subs:
                replaced_str = replaced_str.replace(key, subs[key])
            if os.path.exists(replaced_str):
                replaced_str = os.path.relpath(replaced_str, master_subdir).replace('\\','/')
            return replaced_str
        else:
            return param
    else:
        return param

def checkValue(self, d, default=None):
    if isinstance(d,dict):
        if 'type' in d and d['type'] == 'list':
            if 'default' in d:
                return [a if a is not None else b for a,b in zip(d['value'],d['default'])]
            else:
                if isinstance(d['value'], list):
                    return [val if val is not None else default for val in d['value']]
                else:
                    return None
        else:
            d['value'] = expand_substitution(self, d['value'])
            return d['value'] if d['value'] is not None else d['default'] if 'default' in d else default
    elif isinstance(d, str):
        return getattr(self, d) if hasattr(self, d) and getattr(self, d) is not None else default

def merge_two_dicts(x, y):
    if x is None and y is None:
        return OrderedDict()
    elif x is None:
        return y
    elif y is None:
        return x
    else:
        z = x.copy()   # start with x's keys and values
        z.update(y)    # modifies z with y's keys and values & returns None
        return OrderedDict(z)

def clean_directory(folder):
    for the_file in os.listdir(folder):
        file_path = os.path.join(folder, the_file)
        try:
            if os.path.isfile(file_path):
                os.unlink(file_path)
            #elif os.path.isdir(file_path): shutil.rmtree(file_path)
        except Exception as e:
            print(e)

class Framework(object):

    def __init__(self, directory='test', master_lattice=None, overwrite=None, runname='CLARA_240', clean=False, verbose=True):
        super(Framework, self).__init__()
        global master_lattice_location, master_subdir
        self.verbose = verbose
        self.subdir = directory
        self.elementObjects = OrderedDict()
        self.latticeObjects = OrderedDict()
        self.commandObjects = OrderedDict()
        self.groupObjects = OrderedDict()

        self.basedirectory = os.getcwd()
        self.filedirectory = os.path.dirname(os.path.abspath(__file__))
        self.overwrite = overwrite
        self.runname = runname
        self.subdirectory = self.basedirectory+'/'+self.subdir
        master_subdir = self.subdirectory
        if not os.path.exists(self.subdirectory):
            os.makedirs(self.subdirectory)
        else:
            if clean == True:
                clean_directory(self.subdirectory)
        if self.overwrite == None:
            self.overwrite = True

        if master_lattice is None:
            master_lattice_location = (os.path.relpath(os.path.dirname(os.path.abspath(__file__)) + '/../MasterLattice/')+'/').replace('\\','/')
        else:
            master_lattice_location = master_lattice

        self.executables = {'generator': [master_lattice_location+'Codes/generator'], 'astra': [master_lattice_location+'Codes/astra'],
                            'elegant': [master_lattice_location+'Codes/elegant'], 'csrtrack': [master_lattice_location+'Codes/csrtrack'],
                            'gpt': [r'C:/Program Files/General Particle Tracer/bin/gpt.exe']}


    def defineASTRACommand(self, command=['astra']):
        """Modify the defined ASTRA command variable"""
        self.executables['astra'] = command

    def defineElegantCommand(self, command=['astra']):
        """Modify the defined Elegant command variable"""
        self.executables['elegant'] = command

    def defineGeneratorCommand(self,command=['generator']):
        self.executables['generator'] = command

    def defineCSRTrackCommand(self,command=['csrtrack']):
        self.executables['csrtrack'] = command

    def load_Elements_File(self, input):
        if isinstance(input,(list,tuple)):
            filename = input
        else:
            filename = [input]
        for f in filename:
            with file(master_lattice_location + f, 'r') as stream:
                elements = yaml.load(stream)['elements']
            for name, elem in elements.iteritems():
                self.read_Element(name, elem)

    def loadSettings(self, filename='short_240.settings'):
        """Load Lattice Settings from file"""
        global master_run_no
        stream = file(master_lattice_location+filename, 'r')
        self.settings = yaml.load(stream)
        self.globalSettings = self.settings['global']
        master_run_no = self.globalSettings['run_no'] if 'run_no' in self.globalSettings else 1
        if 'generator' in self.settings:
            self.generatorSettings = self.settings['generator']
            self.add_Generator(**self.generatorSettings)
        self.fileSettings = self.settings['files']
        elements = self.settings['elements']
        self.groups = self.settings['groups'] if 'groups' in self.settings and self.settings['groups'] is not None else {}
        stream.close()

        for name, elem in self.groups.iteritems():
            group = globals()[elem['type']](name, self.elementObjects, **elem)
            self.groupObjects[name] = group

        for name, elem in elements.iteritems():
            self.read_Element(name, elem)

        for name, lattice in self.fileSettings.iteritems():
            self.read_Lattice(name, lattice)

    def read_Lattice(self, name, lattice):
        code = lattice['code'] if 'code' in lattice else 'astra'
        self.latticeObjects[name] = globals()[lattice['code'].lower()+'Lattice'](name, lattice, self.elementObjects, self.groupObjects, self.settings, self.executables)

    def read_Element(self, name, element, subelement=False):
        if name == 'filename':
            self.load_Elements_File(element)
        else:
            if subelement:
                self.add_Element(name, subelement=True, **element)
            else:
                self.add_Element(name, **element)
            if 'sub_elements' in element:
                for name, elem in element['sub_elements'].iteritems():
                    self.read_Element(name, elem, True)

    def add_Element(self, name=None, type=None, **kwargs):
        if name == None:
            if not 'name' in kwargs:
                raise NameError('Element does not have a name')
            else:
                name = kwargs['name']
        try:
            element = globals()[type](name, type, **kwargs)
            self.elementObjects[name] = element
            return element
        except:
            raise NameError('Element \'%s\' does not exist' % type)

    def getElement(self, element, param=None):
        if element in self.elementObjects:
            if param is not None:
                return self.elementObjects[element][param]
            else:
                return self.elementObjects[element]
        else:
            print 'WARNING: Element ', element,' does not exist'
            return {}

    def getElementType(self, type, setting=None):
        return [self.elementObjects[element] if setting is None else self.elementObjects[element][setting] for element in self.elementObjects.keys() if self.elementObjects[element].objectType.lower() == type.lower()]

    def setElementType(self, type, setting, values):
        elems = self.getElementType(type)
        if len(elems) == len(values):
            for e, v  in zip(elems, values):
                e[setting] = v
        else:
            print len(elems), len(values)
            raise ValueError

    def add_Generator(self, default=None, **kwargs):
        if default in astra_generator_keywords['defaults']:
            self.generator = frameworkGenerator(self.executables, **merge_two_dicts(kwargs,astra_generator_keywords['defaults'][default]))
        else:
            self.generator = frameworkGenerator(self.executables, **kwargs)

    def __getitem__(self,key):
        if key in self.elementObjects:
            return self.elementObjects[key]
        elif key in self.latticeObjects:
            return self.latticeObjects[key]
        elif key in self.groupObjects:
            return self.groupObjects[key]
        elif hasattr(self, key):
            return getattr(self,key.lower())

    @property
    def elements(self):
        return self.elementObjects.keys()

    @property
    def lines(self):
        return self.latticeObjects.keys()

    @property
    def commands(self):
        return self.commandObjects.keys()

    def modifyElement(self, elementName, parameter, value):
        setattr(self.elementObjects[elementName], parameter, value)

    def track(self, files=None, startfile=None, preprocess=True, write=True, track=True, postprocess=True):
        if files is None:
            files = ['generator'] + self.lines if hasattr(self, 'generator') else self.lines
        if startfile is not None:
            index = files.index(startfile)
            files = files[index:]
        if self.verbose:
            format_custom_text = progressbar.FormatCustomText(
                'File: %(running)s', {'running': ''}
            )
            bar = progressbar.ProgressBar(widgets=[format_custom_text, progressbar.Percentage(), progressbar.Bar(), progressbar.Percentage(),], max_value=len(files))
            format_custom_text.update_mapping(running=files[0]+'  ')
            for i in bar(range(len(files))):
                l = files[i]
                if l == 'generator' and hasattr(self, 'generator'):
                    format_custom_text.update_mapping(running='Generator  ')
                    if write:
                        self.generator.write()
                    if track:
                        self.generator.run()
                    if postprocess:
                        self.generator.astra_to_hdf5()
                else:
                    if i == (len(files) - 1):
                        format_custom_text.update_mapping(running='Finished')
                    else:
                        format_custom_text.update_mapping(running=files[i+1]+'  ')
                    if preprocess:
                        self.latticeObjects[l].preProcess()
                    if write:
                        self.latticeObjects[l].write()
                    if track:
                        self.latticeObjects[l].run()
                    if postprocess:
                        self.latticeObjects[l].postProcess()
        else:
            for i in range(len(files)):
                l = files[i]
                if l == 'generator' and hasattr(self, 'generator'):
                    if write:
                        self.generator.write()
                    if track:
                        self.generator.run()
                    if postprocess:
                        self.generator.astra_to_hdf5()
                else:
                    if preprocess:
                        self.latticeObjects[l].preProcess()
                    if write:
                        self.latticeObjects[l].write()
                    if track:
                        self.latticeObjects[l].run()
                    if postprocess:
                        self.latticeObjects[l].postProcess()

class frameworkLattice(object):
    def __init__(self, name, file_block, elementObjects, groupObjects, settings, executables):
        super(frameworkLattice, self).__init__()
        self.objectName = name
        self.allElementObjects = elementObjects
        self.groupObjects = groupObjects
        self.allElements = self.allElementObjects.keys()
        self.file_block = file_block
        self.settings = settings
        self.globalSettings = settings['global']
        self.groupSettings = file_block['groups'] if 'groups' in file_block else {}
        self.update_groups()
        self.executables = executables

    def update_groups(self):
        for g in self.groupSettings.keys():
            if g in self.groupObjects:
                self.groupObjects[g].update(**self.groupSettings[g])

    def getElement(self, element, param=None):
        if element in self.allElements:
            if param is not None:
                return self.allElementObjects[element][param]
            else:
                return self.allElements[element]
        else:
            print 'WARNING: Element ', element,' does not exist'
            return {}

    def getElementType(self, type):
        return [self.elements[element] for element in self.elements.keys() if self.elements[element].objectType.lower() == type.lower()]

    def setElementType(self, type, setting, values):
        elems = self.getElementType(type)
        if len(elems) == len(values):
            for e, v  in zip(elems, values):
                e[setting] = v
        else:
            raise ValueError

    @property
    def quadrupoles(self):
        return self.getElementType('quadrupole')

    @property
    def cavities(self):
        return self.getElementType('cavity')

    @property
    def solenoids(self):
        return self.getElementType('solenoid')

    @property
    def dipoles(self):
        return self.getElementType('dipole')

    @property
    def kickers(self):
        return self.getElementType('kicker')

    @property
    def dipoles_and_kickers(self):
        return sorted(self.getElementType('dipole') + self.getElementType('kicker'), key=lambda x: x.position_end[2])

    @property
    def screens(self):
        return self.getElementType('screen')

    @property
    def screens_and_bpms(self):
        return sorted(self.getElementType('screen') + self.getElementType('beam_position_monitor'), key=lambda x: x.position_end[2])

    @property
    def lines(self):
        return self.lineObjects.keys()

    @property
    def start(self):
        if 'start_element' in self.file_block['output']:
            return self.file_block['output']['start_element']
        elif 'zstart' in self.file_block['output']:
            for e in self.allElementObjects.keys():
                if self.allElementObjects[e].position_start[2] == self.file_block['output']['zstart']:
                    return e
        else:
            return self.allElementObjects[0]
    @property
    def startObject(self):
        return self.allElementObjects[self.start]

    @property
    def end(self):
        if 'end_element' in self.file_block['output']:
            return self.file_block['output']['end_element']
        elif 'zstop' in self.file_block['output']:
            endelems = []
            for e in self.allElementObjects.keys():
                if self.allElementObjects[e]['position_end'] == self.file_block['output']['zstop']:
                    endelems.append(e)
                elif self.allElementObjects[e]['position_end'] > self.file_block['output']['zstop'] and len(endelems) == 0:
                    endelems.append(e)
            return endelems[-1]
        else:
            return self.allElementObjects[0]
    @property
    def endObject(self):
        return self.allElementObjects[self.end]

    # @property
    def endScreen(self, **kwargs):
        return screen(name='end', type='screen', position_start=self.endObject.position_start, position_end=self.endObject.position_start, global_rotation=self.endObject.global_rotation, **kwargs)

    @property
    def elements(self):
        index_start = self.allElements.index(self.start)
        index_end = self.allElements.index(self.end)
        f = OrderedDict([[e,self.allElementObjects[e]] for e in self.allElements[index_start:index_end+1]])
        return f

    def createDrifts(self):
        """Insert drifts into a sequence of 'elements'"""
        positions = []
        elementno = 0
        newelements = OrderedDict()
        for name in self.elements.keys():
            pos = np.array(self.allElementObjects[name].position_start)
            positions.append(pos)
            positions.append(self.allElementObjects[name].position_end)
        positions = positions[1:]
        positions.append(positions[-1])
        driftdata = zip(self.elements.iteritems(), list(chunks(positions, 2)))
        for e, d in driftdata:
            newelements[e[0]] = e[1]
            if len(d) > 1:
                x1, y1, z1 = d[0]
                x2, y2, z2 = d[1]
                length = np.sqrt((x2-x1)**2 + (z2-z1)**2)
                if length > 0:
                    elementno += 1
                    name = 'drift'+str(elementno)
                    newdrift = drift(name, type='drift', **{'length': length,
                     'position_start': list(d[0]),
                     'position_end': list(d[1])
                    })
                    newelements[name] = newdrift
        return newelements

    def write(self):
        pass

    def run(self):
        """Run the code with input 'filename'"""
        command = self.executables[self.code] + [self.objectName]
        with open(os.path.relpath(master_subdir+'/'+self.objectName+'.log', '.'), "w") as f:
            subprocess.call(command, stdout=f, cwd=master_subdir)

    def postProcess(self):
        pass

    def preProcess(self):
        pass

class frameworkObject(object):

    objectDefaults = {'length': 0}

    def __init__(self, objectName=None, objectType=None, **kwargs):
        super(frameworkObject, self).__init__()
        if objectName == None:
            raise NameError('Command does not have a name')
        if objectType == None:
            raise NameError('Command does not have a type')
        super(frameworkObject, self).__setattr__('objectProperties', dict())
        super(frameworkObject, self).__setattr__('objectName', objectName)
        super(frameworkObject, self).__setattr__('objectType', objectType)
        if self.objectType in commandkeywords:
            self.allowedKeyWords = commandkeywords[self.objectType]
        elif self.objectType in elementkeywords:
            self.allowedKeyWords = elementkeywords[self.objectType]['keywords'] + elementkeywords['common']['keywords']
            if 'framework_keywords' in  elementkeywords[self.objectType]:
                 self.allowedKeyWords += elementkeywords[self.objectType]['framework_keywords']
        else:
            print 'Unknown type = ', objectType
            raise NameError
        self.allowedKeyWords = map(lambda x: x.lower(), self.allowedKeyWords)
        for key, value in kwargs.iteritems():
            self.add_property(key, value)

    def add_property(self, key, value):
        key = key.lower()
        if key in self.allowedKeyWords:
            try:
                self.objectProperties[key] = value
                # self.__setattr__(key.lower(), value)
            except:
                print key

    def add_default(self, key, value):
        self.objectDefaults[key] = value

    @property
    def parameters(self):
        return self.objectProperties

    def __getattr__(self, key):
        key = key.lower()
        if key in self.objectProperties:
            return self.objectProperties[key]
        elif key in self.objectDefaults:
            return self.objectDefaults[key]
        else:
            return None

    def __getitem__(self, key):
        key = key.lower()
        if key in self.objectProperties:
            return self.objectProperties[key]
        elif key in self.objectDefaults:
            return self.objectDefaults[key]
        else:
            return None

    def __setattr__(self, key, value):
        self.objectProperties.__setitem__(key.lower(),value)

    def __setitem__(self, key, value):
        self.objectProperties.__setitem__(key.lower(),value)

    def long(self, value):
        return int(value)

    def double(self, value):
        return float(value)

    def string(self, value):
        return str(value)

class elegantLattice(frameworkLattice):
    def __init__(self, *args, **kwargs):
        super(elegantLattice, self).__init__(*args, **kwargs)
        self.code = 'elegant'
        # if 'input' in self.file_block and 'particle_definition' in self.file_block['input'] and self.file_block['input'].particle_definition == 'initial_distribution':
        #         self.particle_definition = 'laser'
        # else:
        self.particle_definition = self.allElementObjects[self.start].objectName

    def writeElements(self):
        self.q = charge(name='START', type='charge',**{'total': 250e-12})
        self.w = self.endScreen(output_filename=self.end+'.SDDS')
        elements = self.createDrifts()
        fulltext = ''
        fulltext += self.q.write_Elegant()
        for element in elements.values():
            fulltext += element.write_Elegant()
        fulltext += self.w.write_Elegant()
        fulltext += self.objectName+': Line=(START, '
        for e in elements.keys():
            if len((fulltext + e).splitlines()[-1]) > 60:
                fulltext += '&\n'
            fulltext += e+', '
        return fulltext[:-2] + ', END )\n'

    def write(self):
        self.lattice_file = master_subdir+'/'+self.objectName+'.lte'
        saveFile(self.lattice_file, self.writeElements())
        commandFile = elegantCommandFile(lattice=self.objectName, elegantbeamfilename=self.particle_definition+'.sdds')
        self.command_file = master_subdir+'/'+self.objectName+'.ele'
        saveFile(self.command_file, commandFile.write())

    def preProcess(self):
        prefix = self.file_block['input']['prefix'] if 'input' in self.file_block and 'prefix' in self.file_block['input'] else ''
        # print 'prefix = ', prefix
        self.hdf5_to_sdds(prefix)

    def postProcess(self):
        for s in self.screens:
            s.sdds_to_hdf5()
        self.w.sdds_to_hdf5()

    def hdf5_to_sdds(self, prefix=''):
        HDF5filename = prefix+self.particle_definition+'.hdf5'
        beam.read_HDF5_beam_file(master_subdir + '/' + HDF5filename)
        sddsbeamfilename = self.particle_definition+'.sdds'
        beam.write_SDDS_file(master_subdir + '/' + sddsbeamfilename)

    def run(self):
        """Run the code with input 'filename'"""
        command = self.executables[self.code] + [self.objectName+'.ele']
        with open(os.path.relpath(master_subdir+'/'+self.objectName+'.log', '.'), "w") as f:
            subprocess.call(command, stdout=f, cwd=master_subdir)

class elegantCommandFile(object):
    def __init__(self, lattice='', elegantbeamfilename='', *args, **kwargs):
        super(elegantCommandFile, self).__init__(*args, **kwargs)
        self.commandObjects = OrderedDict()
        self.lattice_filename = lattice+'.lte'
        self.elegantbeamfilename = elegantbeamfilename
        self.addCommand(type='run_setup',lattice=self.lattice_filename, \
            use_beamline=lattice,p_central=np.mean(beam.BetaGamma), \
            centroid='%s.cen',always_change_p0 = 1, \
            sigma='%s.sig')
        self.addCommand(type='run_control',n_steps=1, n_passes=1)
        self.addCommand(type='twiss_output',matched = 0,output_at_each_step=0,radiation_integrals=1,statistics=1,filename="%s.twi",
        beta_x  = beam.beta_x,
        alpha_x = beam.alpha_x,
        beta_y  = beam.beta_y,
        alpha_y = beam.alpha_y)
        self.addCommand(type='sdds_beam', input=self.elegantbeamfilename)
        self.addCommand(type='track')

    def addCommand(self, name=None, **kwargs):
        if name == None:
            if not 'name' in kwargs:
                if not 'type' in kwargs:
                    raise NameError('Command does not have a name')
                else:
                    name = kwargs['type']
            else:
                name = kwargs['name']
        command = frameworkCommand(name, **kwargs)
        self.commandObjects[name] = command
        return command

    def write(self):
        output = ''
        for c in self.commandObjects.values():
            output += c.write()
        return output

class frameworkCommand(frameworkObject):

    def __init__(self, name=None, type=None, **kwargs):
        super(frameworkCommand, self).__init__(name, type, **kwargs)
        if not type in commandkeywords:
            raise NameError('Command \'%s\' does not exist' % commandname)

    def write(self):
        wholestring=''
        string = '&'+self.objectType+'\n'
        for key, value in self.objectProperties.iteritems():
            if not key =='name' and not key == 'type' and key in commandkeywords[self.objectType]:
                string+='\t'+key+' = '+str(value)+'\n'
        string+='&end\n'
        return string

class astraLattice(frameworkLattice):
    def __init__(self, *args, **kwargs):
        super(astraLattice, self).__init__(*args, **kwargs)
        self.code = 'astra'
        self.headers = OrderedDict()
        self.starting_offset = eval(expand_substitution(self, self.file_block['starting_offset'])) if 'starting_offset' in self.file_block else [0,0,0]
        # self.starting_rotation = eval(expand_substitution(self, str(self.file_block['starting_rotation']))) if 'starting_rotation' in self.file_block else 0
        self.starting_rotation = -1*self.allElementObjects[self.start].global_rotation[2] if self.allElementObjects[self.start].global_rotation is not None else 0
        for d in self.dipoles:
            self.starting_rotation -= d.angle
        if 'input' in self.file_block:
            self.headers['newrun'] = astra_newrun(self.starting_offset, self.starting_rotation, **merge_two_dicts(self.file_block['input'],self.globalSettings['ASTRAsettings']))
        else:
            self.headers['newrun'] = astra_newrun(self.starting_offset, self.starting_rotation, **self.globalSettings['ASTRAsettings'])
        self.headers['output'] = astra_output(self.screens_and_bpms, self.starting_offset, self.starting_rotation, **merge_two_dicts(self.file_block['output'],self.globalSettings['ASTRAsettings']))
        self.headers['charge'] = astra_charge(**merge_two_dicts(self.file_block['charge'],self.globalSettings['ASTRAsettings']))
        if self.headers['newrun'].particle_definition == 'initial_distribution':
            self.headers['newrun'].particle_definition = 'laser.astra'
        else:
            self.headers['newrun'].particle_definition = self.allElementObjects[self.start].objectName+'.astra'

    def writeElements(self):
        fulltext = ''
        ''' Create objects for the newrun, output and charge blocks'''
        self.headers['output'].start_element = self.allElementObjects[self.start]
        self.headers['output'].end_element = self.allElementObjects[self.end]
        ''' write the headers and their elements '''
        for header in self.headers:
            fulltext += self.headers[header].write_ASTRA()+'\n'
        counter = frameworkCounter(sub={'kicker': 'dipole'})
        for t in [['cavities'], ['solenoids'], ['quadrupoles'], ['dipoles', 'dipoles_and_kickers']]:
            fulltext += '&' + section_header_text_ASTRA[t[0]]['header']+'\n'
            elements = getattr(self, t[-1])
            fulltext += section_header_text_ASTRA[t[0]]['bool']+' = '+str(len(elements) > 0)+'\n'
            for element in elements:
                element.starting_offset = self.starting_offset
                element.starting_rotation = self.starting_rotation
                elemstr = element.write_ASTRA(counter.counter(element.objectType))
                if elemstr is not None and not elemstr == '':
                    fulltext += elemstr+'\n'
                    if element.objectType == 'kicker':
                        counter.add(element.objectType,2)
                    else:
                        counter.add(element.objectType)
            fulltext += '\n/\n'
        return fulltext

    def write(self):
        self.code_file = master_subdir+'/'+self.objectName+'.in'
        saveFile(self.code_file, self.writeElements())

    def preProcess(self):
        prefix = self.file_block['input']['prefix'] if 'input' in self.file_block and 'prefix' in self.file_block['input'] else ''
        self.headers['newrun'].hdf5_to_astra(prefix)

    def postProcess(self):
        for e in self.screens_and_bpms:
            if not self.starting_offset == [0,0,0]:
                e.zstart = self.allElementObjects[self.start].start
            else:
                e.zstart = [0,0,0]
            e.astra_to_hdf5(self.objectName)
        self.astra_to_hdf5()

    def astra_to_hdf5(self):
        if not self.starting_offset == [0,0,0]:
            zstart = self.allElementObjects[self.start].start
        else:
            zstart = [0,0,0]
        endpos = np.array(self.allElementObjects[self.end].end)-np.array(zstart)
        astrabeamfilename = self.objectName + '.' + str(int(round(endpos[2]*100))).zfill(4) + '.' + str(master_run_no).zfill(3)
        beam.read_astra_beam_file(master_subdir + '/' + astrabeamfilename, normaliseZ=False)
        beam.rotate_beamXZ(-1*self.starting_rotation, preOffset=[0,0,0], postOffset=-1*np.array(self.starting_offset))
        HDF5filename = self.allElementObjects[self.end].objectName+'.hdf5'
        beam.write_HDF5_beam_file(master_subdir + '/' + HDF5filename, centered=False, sourcefilename=astrabeamfilename, pos=self.allElementObjects[self.end].middle)

class gptLattice(frameworkLattice):
    def __init__(self, *args, **kwargs):
        super(gptLattice, self).__init__(*args, **kwargs)
        self.code = 'gpt'
        self.headers = OrderedDict()
        # self.headers['newrun'] = astra_newrun(**merge_two_dicts(self.file_block['input'],self.globalSettings['ASTRAsettings']))
        # self.headers['output'] = astra_output(self.screens_and_bpms, **merge_two_dicts(self.file_block['output'],self.globalSettings['ASTRAsettings']))
        # self.headers['charge'] = astra_charge(**merge_two_dicts(self.file_block['charge'],self.globalSettings['ASTRAsettings']))
        # if self.headers['newrun'].particle_definition == 'initial_distribution':
        #             self.headers['newrun'].particle_definition = 'laser.astra'
        # else:
        #     self.headers['newrun'].particle_definition = self.allElementObjects[self.start].objectName+'.astra'
    def writeElements(self):
        fulltext = ''
        for element in self.elements.values():
            fulltext += element.write_GPT(1)
        fulltext += self.endScreen().write_GPT(1)
        return fulltext

    def write(self):
        self.code_file = master_subdir+'/'+self.objectName+'.in'
        # saveFile(self.code_file, self.writeElements())
        return self.writeElements()

    def preProcess(self):
        self.headers['newrun'].hdf5_to_astra()

    def postProcess(self):
        for e in self.screens_and_bpms:
            e.astra_to_hdf5(self.objectName)
        self.astra_to_hdf5()

    def astra_to_hdf5(self):
        astrabeamfilename = self.objectName + '.' + str(int(round((self.allElementObjects[self.end].position_end[2])*100))).zfill(4) + '.' + str(master_run_no).zfill(3)
        beam.read_astra_beam_file(master_subdir + '/' + astrabeamfilename, normaliseZ=False)
        HDF5filename = self.allElementObjects[self.end].objectName+'.hdf5'
        beam.write_HDF5_beam_file(master_subdir + '/' + HDF5filename, centered=False, sourcefilename=astrabeamfilename)

class csrtrackLattice(frameworkLattice):
    def __init__(self, *args, **kwargs):
        super(csrtrackLattice, self).__init__(*args, **kwargs)
        self.code = 'csrtrack'
        self.particle_definition = ''
        self.CSRTrackelementObjects = OrderedDict()
        self.set_particles_filename()

    def set_particles_filename(self):
        self.CSRTrackelementObjects['particles'] = csrtrack_particles(particle_definition=self.particle_definition)
        self.CSRTrackelementObjects['particles'].format = 'astra'
        if self.particle_definition == 'initial_distribution':
            self.CSRTrackelementObjects['particles'].particle_definition = 'laser.astra'
            self.CSRTrackelementObjects['particles'].add_default('array', '#file{name=laser.astra}')
        else:
            self.CSRTrackelementObjects['particles'].particle_definition = self.allElementObjects[self.start].objectName
            self.CSRTrackelementObjects['particles'].add_default('array', '#file{name='+self.allElementObjects[self.start].objectName+'.astra'+'}')

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
            fulltext += e.write_CSRTrack(counter.counter(e.objectType))
            counter.add(e.objectType)
        fulltext += self.endScreen().write_CSRTrack(counter.counter(self.endScreen().objectType))
        fulltext += '}\n'
        self.set_particles_filename()
        self.setCSRMode()
        self.CSRTrackelementObjects['track_step'] = csrtrack_track_step()
        self.CSRTrackelementObjects['tracker'] = csrtrack_tracker(end_time_marker='screen'+str(counter.counter(self.endScreen().objectType))+'a')
        self.CSRTrackelementObjects['monitor'] = csrtrack_monitor(name=self.end+'.fmt2')
        for c in self.CSRTrackelementObjects:
            fulltext += self.CSRTrackelementObjects[c].write_CSRTrack()
        return fulltext

    def write(self):
        self.code_file = master_subdir+'/csrtrk.in'
        saveFile(self.code_file, self.writeElements())

    def preProcess(self):
        self.CSRTrackelementObjects['particles'].hdf5_to_astra()

    def postProcess(self):
        self.CSRTrackelementObjects['monitor'].csrtrack_to_hdf5()

class frameworkGroup(object):
    def __init__(self, name, elementObjects, type, elements, **kwargs):
        super(frameworkGroup, self).__init__()
        self.objectName = name
        self.type = type
        self.elements = elements
        self.allElementObjects = elementObjects

    def change_Parameter(self, p ,v):
        for e in self.elements:
            setattr(self.allElementObjects[e], p, v)

    def __repr__(self):
        return [self.allElementObjects[e].theta for e in self.elements]

    def __str__(self):
        return str([self.allElementObjects[e].objectName for e in self.elements])

class chicane(frameworkGroup):
    def __init__(self, name, elementObjects, type, elements, **kwargs):
        super(chicane, self).__init__(name, elementObjects, type, elements, **kwargs)

    def xform(self, theta, tilt, length, x, r):
        """Calculate the change on local coordinates through an element"""
        theta = theta if abs(theta) > 1e-9 else 1e-9
        tiltMatrix = np.matrix([
            [np.cos(tilt), -np.sin(tilt), 0],
            [np.sin(tilt), np.cos(tilt), 0],
            [0, 0, 1]
        ])
        angleMatrix = np.matrix([
            [length/theta*(np.cos(theta)-1)],
            [0],
            [length/theta*np.sin(theta)]
        ])
        dx = np.dot(r, angleMatrix)
        rt = np.transpose(r)
        n = rt[1]*np.cos(tilt)-rt[0]*np.sin(tilt)
        crossMatrix = np.matrix([
            np.cross(rt[0], n),
            np.cross(rt[1], n),
            np.cross(rt[2], n)
        ])*np.sin(theta)
        rp = np.outer(np.dot(rt,n), n)*(1-np.cos(theta))+rt*np.cos(theta)+crossMatrix
        return [np.array(x + dx), np.array(np.transpose(rp))]

    def update(self, **kwargs):
        if 'dipoleangle' in kwargs:
            self.set_angle(kwargs['dipoleangle'])
        if 'width' in kwargs:
            self.change_Parameter('width', kwargs['width'])
        if 'gap' in kwargs:
            self.change_Parameter('gap', kwargs['gap'])

    @property
    def angle(self):
        obj = [self.allElementObjects[e] for e in self.elements]
        return obj[0].angle

    def set_angle(self, a):
        obj = [self.allElementObjects[e] for e in self.elements]
        starting_angle = obj[0].theta
        for i in [0,1,2,3]:
            x1 = np.transpose([self.allElementObjects[self.elements[i]].position_start])
            angle = obj[i].angle
            obj[i].global_rotation[2] = starting_angle
            obj[i].angle = -1*np.sign(angle)*a
            localXYZ = self.xform(starting_angle, 0, 0, x1, np.identity(3))[1]
            x1, localXYZ = self.xform(obj[i].angle, 0, obj[i].length, x1, localXYZ)
            xstart, ystart, zstart = x1
            obj[i].position_end[0] = chop(float(xstart))
            if i < 3:
                xend, yend, zend = obj[i+1].position_start
                angle = starting_angle + obj[i].angle
                length = float((-zstart + zend) * (1.0/np.cos(angle)))
                endx = chop(float(xstart + np.tan(-1*angle)*length))
                obj[i+1].position_start[0] =  endx
                starting_angle += obj[i].angle

    def __str__(self):
        return str([[self.allElementObjects[e].objectName, self.allElementObjects[e].angle, self.allElementObjects[e].position_start[0], self.allElementObjects[e].position_end[0]] for e in self.elements])

class frameworkCounter(dict):
    def __init__(self, sub={}):
        super(frameworkCounter, self).__init__()
        self.sub = sub

    def counter(self, type):
        if type not in self:
            return 1
        return self[type] + 1

    def value(self, type):
        if type not in self:
            return 1
        return self[type]

    def add(self, type, n=1):
        type = self.sub[type] if type in self.sub else type
        if type not in self:
            self[type] = n
        else:
            self[type] += n
        return self[type]

    def subtract(self, type):
        type = self.sub[type] if type in self.sub else type
        if type not in self:
            self[type] = 0
        else:
            self[type] = self[type] - 1 if self[type] > 0 else 0
        return self[type]

class frameworkGenerator(object):
    def __init__(self, executables, **kwargs):
        super(frameworkGenerator, self).__init__()
        self.executables = executables
        self.objectDefaults = {}
        self.objectProperties = {}
        self.allowedKeyWords = astra_generator_keywords['keywords'] + astra_generator_keywords['framework_keywords']
        self.allowedKeyWords = map(lambda x: x.lower(), self.allowedKeyWords)
        for key, value in kwargs.iteritems():
            key = key.lower()
            if key in self.allowedKeyWords:
                try:
                    self.objectProperties[key] = value
                    self.__setattr__(key, value)
                except:
                    print 'WARNING: Unknown keyword: ', key

    def run(self):
        command = self.executables['generator'] + [self.objectName+'.in']
        with open(os.devnull, "w") as f:
            subprocess.call(command, stdout=f, cwd=master_subdir)

    def load_defaults(self, defaults):
        if isinstance(defaults, str) and defaults in astra_generator_keywords['defaults']:
            self.__init__(**astra_generator_keywords['defaults'][defaults])
        elif isinstance(defaults, dict):
            self.__init__(**defaults)

    @property
    def particles(self):
        return self.number_of_particles if self.number_of_particles is not None else 512

    @particles.setter
    def particles(self, npart):
        self.add_property('number_of_particles', npart)

    def _write_ASTRA(self, d):
        output = ''
        for k, v in d.iteritems():
            val = v['value'] if v['value'] is not None else v['default'] if 'default' in v else None
            if isinstance(val,str):
                param_string = k+' = \''+str(val)+'\',\n'
            else:
                param_string = k+' = '+str(val)+',\n'
            if len((output + param_string).splitlines()[-1]) > 70:
                output += '\n'
            output += param_string
        return output[:-2]

    @property
    def charge(self):
        return self.objectProperties['charge'] if 'charge' in self.objectProperties and self.objectProperties['charge'] is not None else 250e-12

    @property
    def objectName(self):
        return self.objectProperties['name'] if 'name' in self.objectProperties and self.objectProperties['name'] is not None else 'laser'

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
        saveFile(master_subdir+'/'+self.objectName+'.in', output)

    @property
    def parameters(self):
        return self.objectProperties

    def add_property(self, key, value):
        if key.lower() in self.allowedKeyWords:
            self.objectProperties[key.lower()] = value
            self.__setattr__(key.lower(), value)

    def __getattr__(self, a):
        return None

    def astra_to_hdf5(self):
        astrabeamfilename = self.filename
        beam.read_astra_beam_file(master_subdir + '/' + astrabeamfilename, normaliseZ=False)
        HDF5filename = self.filename.replace('.generator','.hdf5')
        beam.write_HDF5_beam_file(master_subdir + '/' + HDF5filename, centered=False, sourcefilename=astrabeamfilename)

class frameworkElement(frameworkObject):

    def __init__(self, elementName=None, elementType=None, **kwargs):
        super(frameworkElement, self).__init__(elementName, elementType, **kwargs)

    def __mul__(self, other):
        return [self.objectProperties for x in range(other)]

    def __rmul__(self, other):
        return [self.objectProperties for x in range(other)]

    def __neg__(self):
        return self

    @property
    def theta(self):
        if hasattr(self, 'global_rotation') and self.global_rotation is not None:
            rotation =  self.global_rotation[2] if len(self.global_rotation) is 3 else self.global_rotation
        else:
            rotation = 0
        if hasattr(self, 'starting_rotation') and self.starting_rotation is not None:
            rotation +=  self.starting_rotation
        return rotation

    def _rotation_matrix(self, theta):
        return np.array([[np.cos(theta), 0, np.sin(theta)], [0, 1, 0], [-1*np.sin(theta), 0, np.cos(theta)]])

    @property
    def rotation_matrix(self):
        return self._rotation_matrix(self.theta)

    def rotated_position(self, pos=[0,0,0], offset=None, theta=None):
        if offset is None:
            if not hasattr(self, 'starting_offset') or self.starting_offset is None:
                offset = [0,0,0]
            else:
                offset = self.starting_offset
        if theta is None:
            return chop(np.dot(np.array(pos) - np.array(offset), self.rotation_matrix), 1e-6)
        else:
            return chop(np.dot(np.array(pos) - np.array(offset), self._rotation_matrix(theta)), 1e-6)

    @property
    def start(self):
        start = self.position_start
        length_vector = self.rotated_position([0,0, 0], offset=[0,0,0], theta=self.theta)
        starting_middle = np.array(self.position_start) + length_vector
        return self.rotated_position(starting_middle, offset=self.starting_offset, theta=self.starting_rotation)

    @property
    def middle(self):
        start = self.position_start
        length_vector = self.rotated_position([0,0, self.length / 2.0], offset=[0,0,0], theta=self.theta)
        starting_middle = np.array(self.position_start) + length_vector
        return self.rotated_position(starting_middle, offset=self.starting_offset, theta=self.starting_rotation)

    @property
    def end(self):
        start = self.position_start
        length_vector = self.rotated_position([0,0, self.length], offset=[0,0,0], theta=self.theta)
        starting_middle = np.array(self.position_start) + length_vector
        return self.rotated_position(starting_middle, offset=self.starting_offset, theta=self.starting_rotation)

    def _write_ASTRA(self, d, n=1):
        output = ''
        for k, v in d.iteritems():
            if checkValue(self, v) is not None:
                if 'type' in v and v['type'] == 'list':
                    for i, l in enumerate(checkValue(self, v)):
                        if n is not None:
                            param_string = k+'('+str(i+1)+','+str(n)+') = '+str(l)+', '
                        else:
                            param_string = k+' = '+str(l)+'\n'
                        if len((output + param_string).splitlines()[-1]) > 70:
                            output += '\n'
                        output += param_string
                elif 'type' in v and v['type'] == 'array':
                    if n is not None:
                        param_string = k+'('+str(n)+') = ('
                    else:
                        param_string = k+' = ('
                    for i, l in enumerate(checkValue(self, v)):
                        param_string += str(l)+', '
                        if len((output + param_string).splitlines()[-1]) > 70:
                            output += '\n'
                    output += param_string[:-2] + '),\n'
                else:
                    if n is not None:
                        param_string = k+'('+str(n)+') = '+str(checkValue(self, v))+', '
                    else:
                        param_string = k+' = '+str(checkValue(self, v))+',\n'
                    if len((output + param_string).splitlines()[-1]) > 70:
                        output += '\n'
                    output += param_string
        return output[:-2]

    def write_ASTRA(self, n):
        return ""

    def _write_Elegant(self):
        wholestring=''
        etype = self._convertType_Elegant(self.objectType)
        string = self.objectName+': '+ etype
        for key, value in merge_two_dicts(self.objectDefaults, self.objectProperties).iteritems():
            key = self._convertKeword_Elegant(key)
            if not key is 'name' and not key is 'type' and not key is 'commandtype' and key in elements_Elegant[etype]:
                tmpstring = ', '+key+' = '+str(value)
                if len(string+tmpstring) > 76:
                    wholestring+=string+',&\n'
                    string=''
                    string+=tmpstring[2::]
                else:
                    string+= tmpstring
        wholestring+=string+';\n'
        return wholestring

    def write_Elegant(self):
        if not self.subelement:
            return self._write_Elegant()

    def _convertType_Elegant(self, etype):
        return type_conversion_rules_Elegant[etype] if etype in type_conversion_rules_Elegant else etype

    def _convertKeword_Elegant(self, keyword):
        return keyword_conversion_rules_Elegant[keyword] if keyword in keyword_conversion_rules_Elegant else keyword

    def write_CSRTrack(self, n=0):
        return ""

    def write_GPT(self, Brho):
        return ""

    def gpt_coordinates(self):
        x,y,z = self.middle
        psi, phi, theta = self.global_rotation
        output =''
        for c in [x, y, z]:
            output += str(c)+', '
        output += 'cos('+str(theta)+'), 0, -sin('+str(theta)+'), 0, 1 ,0'
        return output

class dipole(frameworkElement):

    def __init__(self, name=None, type=None, **kwargs):
        super(dipole, self).__init__(name, type, **kwargs)
        self.add_default('bins',20)
        self.add_default('csr', 1)
        self.add_default('isr', 1)
        self.add_default('n_kicks', 10)
        self.add_default('integration_order', 4)

    @property
    def width(self):
        if 'width' in self.objectProperties:
            return self.objectProperties['width']
        else:
            return 0.2

    def __neg__(self):
        newself = copy.deepcopy(self)
        if 'exit_edge_angle' in newself.objectProperties and 'entrance_edge_angle' in newself.objectProperties:
            e1 = newself['entrance_edge_angle']
            e2 = newself['exit_edge_angle']
            newself.objectProperties['entrance_edge_angle'] = e2
            newself.objectProperties['exit_edge_angle'] = e1
        elif 'entrance_edge_angle' in newself.objectProperties:
            newself.objectProperties['exit_edge_angle'] = newself.objectProperties['entrance_edge_angle']
            del newself.objectProperties['entrance_edge_angle']
        elif 'exit_edge_angle' in newself.objectProperties:
            newself.objectProperties['entrance_edge_angle'] = newself.objectProperties['exit_edge_angle']
            del newself.objectProperties['exit_edge_angle']
        newself.objectName = '-'+newself.objectName
        return newself

    def _edge_angles(self, estr):
        if estr in self.objectProperties:
            if isinstance(self.objectProperties[estr], str):
                return checkValue(self, self.objectProperties[estr],0)
            else:
                return self.objectProperties[estr]
        else:
            return 0

    @property
    def rho(self):
        return self.length/self.angle if self.length is not None and abs(self.angle) > 1e-9 else 0

    @property
    def e1(self):
        return self._edge_angles('entrance_edge_angle')
    @property
    def e2(self):
        return self._edge_angles('exit_edge_angle')

    @property
    def corners(self):
        corners = [0,0,0,0]
        if hasattr(self, 'global_rotation') and self.global_rotation is not None:
            rotation =  self.global_rotation[2] if len(self.global_rotation) is 3 else self.global_rotation
        else:
            rotation = 0
        theta = self.e1+rotation
        corners[0] = np.array(map(add, np.transpose(self.position_start), np.dot([-self.width*self.length,0,0], self._rotation_matrix(theta))))
        corners[3] = np.array(map(add, np.transpose(self.position_start), np.dot([self.width*self.length,0,0], self._rotation_matrix(theta))))
        theta = self.angle-self.e2+rotation
        corners[1] = np.array(map(add, np.transpose(self.position_end), np.dot([-self.width*self.length,0,0], self._rotation_matrix(theta))))
        corners[2] = np.array(map(add, np.transpose(self.position_end), np.dot([self.width*self.length,0,0], self._rotation_matrix(theta))))
        corners = map(lambda x: self.rotated_position(x, offset=self.starting_offset, theta=self.starting_rotation), corners)
        return corners

    def write_CSRTrack(self, n):
        z1 = self.position_start[2]
        z2 = self.position_end[2]
        return """dipole{\nposition{rho="""+str(z1)+""", psi="""+str(chop(self.theta+self.e1))+""", marker=d"""+str(n)+"""a}\nproperties{r="""+str(self.rho)+"""}\nposition{rho="""+str(z2)+""", psi="""+str(chop(self.theta+self.angle-self.e2))+""", marker=d"""+str(n)+"""b}\n}\n"""

    def write_ASTRA(self, n):
        if abs(checkValue(self, 'strength', default=0)) > 0 or abs(self.rho) > 0:
            corners = self.corners
            if self.plane is None:
                self.plane = 'horizontal'
            params = OrderedDict([
                ['D_Type', {'value': '\''+self.plane+'\'', 'default': '\'horizontal\''}],
                ['D_Gap', {'type': 'list', 'value': [self.gap, self.gap], 'default': [0.02, 0.02]}],
                ['D1', {'type': 'array', 'value': [corners[3][0],corners[3][2]] }],
                ['D3', {'type': 'array', 'value': [corners[2][0],corners[2][2]] }],
                ['D4', {'type': 'array', 'value': [corners[1][0],corners[1][2]] }],
                ['D2', {'type': 'array', 'value': [corners[0][0],corners[0][2]] }],
                ])
            if abs(checkValue(self, 'strength', default=0)) > 0 or not abs(self.rho) > 0:
                params['D_strength'] = {'value': checkValue(self, 'strength',0), 'default': 1e6}
            else:
                params['D_radius'] =  {'value': self.rho, 'default': 1e6}
            return self._write_ASTRA(params, n)
        else:
            return None

class kicker(dipole):

    def __init__(self, name=None, type=None, **kwargs):
        super(kicker, self).__init__(name, type, **kwargs)

    @property
    def angle(self):
        return 0

    def write_ASTRA(self, n):
        output = ''
        self.plane = 'horizontal'
        hkick = super(kicker, self).write_ASTRA(n)
        if hkick is not None:
            output += hkick
        self.plane = 'vertical'
        vkick = super(kicker, self).write_ASTRA(n+1)
        if vkick is not None:
            output += vkick
        return output

class quadrupole(frameworkElement):

    def __init__(self, name=None, type=None, **kwargs):
        super(quadrupole, self).__init__(name, type, **kwargs)
        self.add_default('k1', 0)
        self.add_default('n_kicks', 20)

    def write_ASTRA(self, n):
        return self._write_ASTRA(OrderedDict([
            ['Q_pos', {'value': self.middle[2], 'default': 0}],
            ['Q_k', {'value': self.k1, 'default': 0}],
            ['Q_length', {'value': self.length, 'default': 0}],
            ['Q_smooth', {'value': self.smooth, 'default': 10}],
            ['Q_bore', {'value': self.bore, 'default': 0.016}],
            ['Q_noscale', {'value': self.scale_field}],
            ['Q_mult_a', {'type': 'list', 'value': self.multipoles}],
        ]), n)

    def write_GPT(self, Brho):
        output = str(self.objectType) + '( "wcs", '+self.gpt_coordinates()+
        ', '+str(self.length)+', '+str(Brho*self.k1)+');\n'
        return output

class cavity(frameworkElement):

    def __init__(self, name=None, type=None, **kwargs):
        super(cavity, self).__init__(name, type, **kwargs)

    def write_ASTRA(self, n):
        if (self.n_cells is 0 or self.n_cells is None) and self.cell_length > 0:
                cells = round((self.length-self.cell_length)/self.cell_length)
                cells = int(cells - (cells % 3))
        elif self.n_cells > 0 and self.cell_length > 0:
            cells = int(self.n_cells - (self.n_cells % 3))
        else:
            cells = None
        return self._write_ASTRA(OrderedDict([
            ['C_pos', {'value': self.start[2], 'default': 0}],
            ['FILE_EFieLD', {'value': '\''+expand_substitution(self, '\''+self.field_definition+'\'')+'\'', 'default': 0}],
            ['C_numb', {'value': cells}],
            ['Nue', {'value': self.frequency / 1e9, 'default': 2998.5}],
            ['MaxE', {'value': self.field_amplitude / 1e6, 'default': 0}],
            ['Phi', {'value': self.phase, 'default': 0.0}],
            ['C_smooth', {'value': self.smooth, 'default': 10}],
        ]), n)

class rf_deflecting_cavity(cavity):

    def __init__(self, name=None, type=None, **kwargs):
        super(rf_deflecting_cavity, self).__init__(name, type, **kwargs)

class solenoid(frameworkElement):

    def __init__(self, name=None, type=None, **kwargs):
        super(solenoid, self).__init__(name, type, **kwargs)

    def write_ASTRA(self, n):
        return self._write_ASTRA(OrderedDict([
            ['S_pos', {'value': self.start[2], 'default': 0}],
            ['FILE_BFieLD', {'value': '\''+expand_substitution(self, '\''+self.field_definition+'\'')+'\''}],
            ['MaxB', {'value': self.field_amplitude, 'default': 0}],
            ['S_smooth', {'value': self.smooth, 'default': 10}],
        ]), n)

class aperture(frameworkElement):

    def __init__(self, name=None, type=None, **kwargs):
        super(aperture, self).__init__(name, type, **kwargs)

    def write_GPT(self, Brho):
        if self.shape == 'elliptical':
            output = 'rmax'
        else:
            output = 'xymax'
        output += '( "wcs", '+self.gpt_coordinates()+', '+str(self.horizontal_size)+', '+str(self.length)+');\n'
        return output

class wall_current_monitor(frameworkElement):

    def __init__(self, name=None, type=None, **kwargs):
        super(wall_current_monitor, self).__init__(name, type, **kwargs)

class screen(frameworkElement):

    def __init__(self, name=None, type=None, **kwargs):
        super(screen, self).__init__(name, type, **kwargs)
        if 'output_filename' not in kwargs:
            self.objectProperties['output_filename'] = self.objectName+'.sdds'

    def write_ASTRA(self, n):
        return self._write_ASTRA(OrderedDict([
            ['Screen', {'value': self.middle[2], 'default': 0}],
        ]), n)

    def write_CSRTrack(self, n):
        z = self.start[2]
        return """quadrupole{\nposition{rho="""+str(z)+""", psi=0.0, marker=screen"""+str(n)+"""a}\nproperties{strength=0.0, alpha=0, horizontal_offset=0,vertical_offset=0}\nposition{rho="""+str(z+1e-6)+""", psi=0.0, marker=screen"""+str(n)+"""b}\n}\n"""

    def write_GPT(self, Brho):
        output = 'screen( "wcs", '+self.gpt_coordinates()+');\n'
        return output

    def astra_to_hdf5(self, lattice):
        for i in [0, -0.001, 0.001]:
            tempfilename = lattice + '.' + str(int(round((self.middle[2]+i-self.zstart[2])*100))).zfill(4) + '.' + str(master_run_no).zfill(3)
            if os.path.isfile(master_subdir + '/' + tempfilename):
                astrabeamfilename = tempfilename
        beam.read_astra_beam_file(master_subdir + '/' + astrabeamfilename, normaliseZ=False)
        beam.rotate_beamXZ(-1*self.starting_rotation, preOffset=[0,0,0], postOffset=-1*np.array(self.starting_offset))

        HDF5filename = self.objectName+'.hdf5'
        beam.write_HDF5_beam_file(master_subdir + '/' + HDF5filename, centered=False, sourcefilename=astrabeamfilename, pos=self.middle)

    def sdds_to_hdf5(self):
        elegantbeamfilename = self.output_filename.replace('.sdds','.SDDS')
        beam.read_SDDS_beam_file(master_subdir + '/' + elegantbeamfilename)
        HDF5filename = self.output_filename.replace('.sdds','.hdf5').replace('.SDDS','.hdf5')
        beam.write_HDF5_beam_file(master_subdir + '/' + HDF5filename, centered=False, sourcefilename=elegantbeamfilename, pos=self.middle)

class monitor(screen):

    def __init__(self, name=None, type=None, **kwargs):
        super(monitor, self).__init__(name, type, **kwargs)

class watch_point(screen):

    def __init__(self, name=None, type=None, **kwargs):
        super(watch_point, self).__init__(name, 'screen', **kwargs)

class beam_position_monitor(screen):

    def __init__(self, name=None, type=None, **kwargs):
        super(beam_position_monitor, self).__init__(name, type, **kwargs)

    def write_ASTRA(self, n):
        return self._write_ASTRA(OrderedDict([
            ['Screen', {'value': self.middle[2], 'default': 0}],
        ]), n)

class collimator(frameworkElement):

    def __init__(self, name=None, type=None, **kwargs):
        super(collimator, self).__init__(name, type, **kwargs)

class marker(screen):

    def __init__(self, name=None, type=None, **kwargs):
        super(marker, self).__init__(name, 'screen', **kwargs)

    def write_CSRTrack(self, n):
        return ''

class drift(frameworkElement):

    def __init__(self, name=None, type=None, **kwargs):
        super(drift, self).__init__(name, type, **kwargs)

class charge(frameworkElement):

    def __init__(self, name=None, type=None, **kwargs):
        super(charge, self).__init__(name, type, **kwargs)

class astra_header(frameworkElement):

    def __init__(self, name=None, type=None, **kwargs):
        super(astra_header, self).__init__(name, type, **kwargs)

    def framework_dict(self):
        return OrderedDict()

    def write_ASTRA(self):
        keyword_dict = {}
        for k in elementkeywords[self.objectType]['keywords']:
            if getattr(self, k.lower()) is not None:
                keyword_dict[k.lower()] = {'value': getattr(self,k.lower())}
        output = '&' + section_header_text_ASTRA[self.objectType]['header']+'\n'
        output += self._write_ASTRA(merge_two_dicts(self.framework_dict(), keyword_dict), None) + '\n/\n'
        return output

class astra_newrun(astra_header):
    def __init__(self, offset, rotation, **kwargs):
        super(astra_header, self).__init__('newrun', 'astra_newrun', **kwargs)
        self.starting_offset = offset
        self.starting_rotation = rotation
        if 'run' not in kwargs:
            self.objectProperties['run'] = 1
        if 'head' not in kwargs:
            self.objectProperties['head'] = 'trial'
        if 'lprompt' not in kwargs:
            self.add_property('lprompt', False)

    def framework_dict(self):
        return OrderedDict([
            ['Distribution', {'value': '\''+self.particle_definition+'\''}],
            ['high_res', {'value': self.high_res, 'default': True}],
            ['auto_phase', {'value': self.auto_phase, 'default': True}]
        ])

    def hdf5_to_astra(self, prefix=''):
        # print 'self.particle_definition = ', self.particle_definition
        HDF5filename = prefix+self.particle_definition.replace('.astra','')+'.hdf5'
        beam.read_HDF5_beam_file(master_subdir + '/' + HDF5filename)
        beam.rotate_beamXZ(self.theta, preOffset=self.starting_offset)
        astrabeamfilename = self.particle_definition
        beam.write_astra_beam_file(master_subdir + '/' + astrabeamfilename, normaliseZ=False)

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
        keyworddict = OrderedDict([
            ['zstart', {'value': self.start_element.start[2]}],
            ['zstop', {'value': self.end_element.end[2]}],
        ])
        for i, element in enumerate(self.screens,1):
            element.starting_offset = self.starting_offset
            element.starting_rotation = self.starting_rotation
            keyworddict['Screen('+str(i)+')'] = {'value': element.middle[2]}
        return keyworddict

class astra_charge(astra_header):
    def __init__(self, **kwargs):
        super(astra_header, self).__init__('charge', 'astra_charge', **kwargs)

    @property
    def space_charge(self):
        return False if self.space_charge_mode == 'False' or self.space_charge_mode == False else True

    @property
    def space_charge_2D(self):
        return True if self.space_charge and self.space_charge_mode != '3D' else False

    @property
    def space_charge_3D(self):
        return True if self.space_charge and not self.space_charge_2D else False

    def framework_dict(self):
        sc_dict = OrderedDict([
            ['Lmirror', {'value': self.cathode, 'default': False}],
            ['LSPCH', {'value': self.space_charge, 'default': True}],
            ['LSPCH3D', {'value': self.space_charge_3D, 'default': True}]
        ])
        if self.space_charge_2D:
            sc_n_dict = OrderedDict([
                ['nrad', {'value': self.sc_2d_nrad, 'default': 32}],
                ['nlong_in', {'value': self.sc_2d_nlong, 'default': 32}],
            ])
        elif self.space_charge_3D:
            sc_n_dict = OrderedDict([
                ['nxf', {'value': self.sc_3d_nxf, 'default': 6}],
                ['nyf', {'value': self.sc_3d_nyf, 'default': 6}],
                ['nzf', {'value': self.sc_3d_nzf, 'default': 6}],
            ])
        else:
            sc_n_dict = OrderedDict([])
        return merge_two_dicts(sc_dict, sc_n_dict)

class csrtrack_element(frameworkElement):

    def __init__(self, elementName=None, elementType=None, **kwargs):
        super(csrtrack_element, self).__init__(elementName, elementType, **kwargs)
        self.header = ''
        if elementName in csrtrack_defaults:
            for k, v in csrtrack_defaults[elementName].iteritems():
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
        for k in elementkeywords[self.objectType]['keywords']:
            k = k.lower()
            if getattr(self,k) is not None:
                output += k+'='+self.CSRTrack_str(getattr(self, k))+'\n'
            elif k in self.objectDefaults:
                output += k+'='+self.CSRTrack_str(self.objectDefaults[k])+'\n'
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

    def __init__(self, end_time_marker=''):
        super(csrtrack_tracker, self).__init__('tracker', 'csrtrack_tracker')
        self.header = 'tracker'
        self.end_time_marker = end_time_marker

class csrtrack_monitor(csrtrack_element):

    def __init__(self, **kwargs):
        super(csrtrack_monitor, self).__init__(elementName='monitor', elementType='csrtrack_monitor', **kwargs)
        self.header = 'monitor'

    def csrtrack_to_hdf5(self):
        csrtrackbeamfilename = self.name
        astrabeamfilename = csrtrackbeamfilename.replace('.fmt2','.astra')
        beam.convert_csrtrackfile_to_astrafile(master_subdir + '/' + csrtrackbeamfilename, master_subdir + '/' + astrabeamfilename)
        beam.read_astra_beam_file(master_subdir + '/' + astrabeamfilename, normaliseZ=False)
        HDF5filename = self.name.replace('.fmt2','.hdf5')
        beam.write_HDF5_beam_file(master_subdir + '/' + HDF5filename, sourcefilename=csrtrackbeamfilename)

class csrtrack_particles(csrtrack_element):

    def __init__(self, **kwargs):
        super(csrtrack_particles, self).__init__('particles', 'csrtrack_particles', **kwargs)
        self.header = 'particles'

    def hdf5_to_astra(self):
        HDF5filename = self.particle_definition+'.hdf5'
        beam.read_HDF5_beam_file(master_subdir + '/' + HDF5filename)
        astrabeamfilename = self.particle_definition+'.astra'
        beam.write_astra_beam_file(master_subdir + '/' + astrabeamfilename, normaliseZ=False)

class gpt_element(frameworkElement):

    def __init__(self, elementName=None, elementType=None, **kwargs):
        super(gpt_element, self).__init__(elementName, elementType, **kwargs)
        if elementName in csrtrack_defaults:
            for k, v in csrtrack_defaults[elementName].iteritems():
                self.add_default(k, v)

    def write_GPT(self):
        output = str(self.elementName) + '('
        for k in elementkeywords[self.objectType]['keywords']:
            k = k.lower()
            if getattr(self,k) is not None:
                output += getattr(self, k)+', '
            elif k in self.objectDefaults :
                output += self.objectDefaults[k]+', '
        output[:-2]+=')\n'
        return output

class gpt_setfile(gpt_element):

    def __init__(self, **kwargs):
        super(gpt_setfile, self).__init__('setfile', 'gpt_setfile', **kwargs)
