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

_mapping_tag = yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG

def dict_representer(dumper, data):
    return dumper.represent_dict(data.iteritems())

def dict_constructor(loader, node):
    return OrderedDict(loader.construct_pairs(node))

yaml.add_representer(OrderedDict, dict_representer)
yaml.add_constructor(_mapping_tag, dict_constructor)


commandkeywords = {}

astra_generator_keywords = {
    'keywords':[
        'fname','add','ipart','species','probe','noise_reduc','high_res','cathode','lprompt', 'q_total','ref_zpos','ref_clock','dist_z','ref_ekin','lt','rt','sig_clock','sig_z','lz','rz',
        'dist_pz','le','dist_x','sig_x','dist_y','sig_y','dist_px','nemit',
    ],
    'defaults': {
        'clara_400_3ps':{
            'add': False,' species': 'electrons', 'probe': True,'noise_reduc': False, 'high_res': True, 'cathode': True, 'lprompt': False, 'ref_zpos': 0, 'ref_clock': 0, 'dist_z': 'p',
            'ref_ekin': 0, 'lt': 3e-3, 'rt': 0.2e-3, 'dist_pz': 'i', 'le': 0.62e-3, 'dist_x': 'radial', 'sig_x': 0.25, 'dist_y': 'r', 'sig_y': 0.25,
        }
    },
    'framework_keywords': [
        'number_of_particles', 'charge', 'filename',
    ]
}


with open(os.path.dirname( os.path.abspath(__file__))+'/elementkeywords.yaml', 'r') as infile:
    elementkeywords = yaml.load(infile)

with open(os.path.dirname( os.path.abspath(__file__))+'/elements_Elegant.yaml', 'r') as infile:
    elements_Elegant = yaml.load(infile)

type_conversion_rules_Elegant = {'dipole': 'csrcsbend', 'quadrupole': 'kquad', 'beam_position_monitor': 'moni', 'screen': 'watch', 'aperture': 'rcol',
                         'collimator': 'ecol', 'monitor': 'moni', 'solenoid': 'mapsolenoid', 'wall_current_monitor': 'charge', 'cavity': 'rfcw',
                         'rf_deflecting_cavity': 'rfdf'}
keyword_conversion_rules_Elegant = {'length': 'l','entrance_edge_angle': 'e1', 'exit_edge_angle': 'e2', 'edge_field_integral': 'fint', 'horizontal_size': 'x_max', 'vertical_size': 'y_max',
                            'field_amplitude': 'volt', 'frequency': 'freq', 'output_filename': 'filename'}

section_header_text_ASTRA = {'cavities': {'header': 'CAVITY', 'bool': 'LEField'},
                             'solenoids': {'header': 'SOLENOID', 'bool': 'LBField'},
                             'quadrupoles': {'header': 'QUADRUPOLE', 'bool': 'LQuad'},
                             'dipoles': {'header': 'DIPOLE', 'bool': 'LDipole'},
                             'newrun': {'header': 'NEWRUN'},
                             'output': {'header': 'OUTPUT'},
                             'charge': {'header': 'CHARGE'},
                            }

def merge_two_dicts(x, y):
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

class framework(object):

    def __init__(self, directory='test', master_lattice=None, overwrite=None, runname='CLARA_240', clean=False):
        super(framework, self).__init__()
        global master_lattice_location, master_subdir
        self.subdir = directory
        self.elementObjects = OrderedDict()
        self.latticeObjects = OrderedDict()
        self.commandObjects = OrderedDict()
        self.executables = {'generator': ['generator'], 'astra': ['astra'], 'elegant': ['elegant']}

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
            if not os.path.exists(self.subdirectory):
                os.makedirs(self.subdirectory)
                self.overwrite = True

        if master_lattice is None:
            master_lattice_location = (os.path.relpath(os.path.dirname(os.path.abspath(__file__)) + '/../MasterLattice/')+'/').replace('\\','/')
        else:
            master_lattice_location = master_lattice

    def defineASTRACommand(self, command=['astra']):
        """Modify the defined ASTRA command variable"""
        self.executables['astra'] = command

    def defineElegantCommand(self, command=['astra']):
        """Modify the defined Elegant command variable"""
        self.executables['elegant'] = command

    def defineGeneratorCommand(self,command=['generator']):
        self.executables['generator'] = command


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
        settings = yaml.load(stream)
        self.globalSettings = settings['global']
        master_run_no = self.globalSettings['run_no'] if 'run_no' in self.globalSettings else 1
        if 'generator' in settings:
            self.generatorSettings = settings['generator']
            self.add_Generator(**self.generatorSettings)
        self.fileSettings = settings['files']
        elements = settings['elements']
        self.groups = settings['groups']
        stream.close()

        for name, elem in elements.iteritems():
            self.read_Element(name, elem)
        for name, elem in self.fileSettings.iteritems():
            self.read_Lattice(name, elem)

    def read_Lattice(self, name, elem):
        self.latticeObjects[name] = frameworkLattice(name, self.elementObjects, elem, self.globalSettings, self.executables)

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

    def track(self, files=None):
        if files is None:
            files = ['generator'] + self.lines
        for l in files:
            if l == 'generator' and hasattr(self, 'generator'):
                print 'Running Generator...'
                self.generator.write()
                self.generator.run()
                self.generator.astra_to_hdf5()
            else:
                print 'Running ',l,'...'
                self.latticeObjects[l].write()
                self.latticeObjects[l].run()
                self.latticeObjects[l].postProcess()

class frameworkLattice(object):
    def __init__(self, name, elementObjects, file_block, globalSettings, executables):
        super(frameworkLattice, self).__init__()
        self.name = name
        self.allElementObjects = elementObjects
        self.allElements = self.allElementObjects.keys()
        self.file_block = file_block
        self.globalSettings = globalSettings
        self.executables = executables

    def getElement(self, element):
        if element in self.elements:
            return self.elements[element]
        else:
            print 'WARNING: Element ', element,' does not exist'
            return {}

    def getElementType(self, type):
        return [self.elements[element] for element in self.elements.keys() if self.elements[element].type.lower() == type.lower()]

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
            pos = np.array(self.getElement(name).position_start)
            positions.append(pos)
            positions.append(self.getElement(name).position_end)
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
        if self.file_block['code'].lower() == 'astra':
            self.code_file = master_subdir+'/'+self.name+'.in'
            saveFile(self.code_file, self.writeElements_ASTRA())
        elif self.file_block['code'].lower() == 'elegant':
            self.code_file = master_subdir+'/'+self.name+'.lte'
            saveFile(self.code_file, self.writeElements_Elegant())

    def writeElements_Elegant(self):
        q = charge(name='START', type='charge',**{'total': 250e-12})
        w = screen(name='END', type='screen',**{'output_filename': self.end+'.sdds'})
        elements = self.createDrifts()
        fulltext = ''
        for element in elements.values():
            fulltext += element.write_Elegant()
        fulltext += q.write_Elegant()
        fulltext += w.write_Elegant()
        fulltext += self.name+' = Line=(START,'
        for e in elements.keys():
            if len((fulltext + e).splitlines()[-1]) > 60:
                fulltext += '&\n'
            fulltext += e+', '
        return fulltext[:-2] + ', END);\n'

    def writeElements_ASTRA(self):
        fulltext = ''
        ''' Create objects for the newrun, output and charge blocks'''
        for t, d in [['newrun','input'], ['output', 'output'], ['charge','charge']]:
            fulltext += '&' + section_header_text_ASTRA[t]['header']+'\n'
            ''' instantiate an object based on the type of header block'''
            header = globals()['astra_'+t](t, 'astra_'+t, **merge_two_dicts(self.file_block[d],self.globalSettings['ASTRAsettings']))
            ''' some special sauce for newrun blocks '''
            if t == 'newrun':
                if header.particle_definition == 'initial_distribution':
                    header.particle_definition = 'laser.astra'
                else:
                    header.particle_definition = self.allElementObjects[self.start].name+'.astra'
                header.hdf5_to_astra()
            ''' some special sauce for output blocks '''
            if t == 'output':
                header.start_element = self.allElementObjects[self.start].position_start[2]
                header.end_element = self.allElementObjects[self.end].position_end[2]
            ''' write the headers and their elements '''
            fulltext += header.write_ASTRA()+'\n'
            if t == 'output':
                counter = frameworkCounter({'beam_position_monitor': 'screen'})
                elements = self.screens_and_bpms
                for element in elements:
                    fulltext += element.write_ASTRA(counter.add(element.type))+'\n'
            fulltext += '\n/\n'
        counter = frameworkCounter(sub={'kicker': 'dipole'})
        for t in [['cavities'], ['solenoids'], ['quadrupoles'], ['dipoles', 'dipoles_and_kickers']]:
            fulltext += '&' + section_header_text_ASTRA[t[0]]['header']+'\n'
            elements = getattr(self, t[-1])
            fulltext += section_header_text_ASTRA[t[0]]['bool']+' = '+str(len(elements) > 0)+'\n'
            for element in elements:
                elemstr = element.write_ASTRA(counter.add(element.type))
                if elemstr is not None and not elemstr == '':
                    fulltext += elemstr+'\n'
                    if element.type == 'kicker':
                        counter.add(element.type)
                else:
                    counter.subtract(element.type)
            fulltext += '\n/\n'
        return fulltext

    def run(self):
        self.runCode(self.file_block['code'].lower(), self.code_file)

    def runCode(self, code, filename=''):
        """Run the code with input 'filename'"""
        command = self.executables[code] + [os.path.relpath(filename, master_subdir)]
        with open(os.path.relpath(filename+'.log', '.'), "w") as f:
            subprocess.call(command, stdout=f, cwd=master_subdir)

    def postProcess(self):
        if self.file_block['code'].lower() == 'astra':
            for e in self.screens_and_bpms:
                e.astra_to_hdf5(self.name)
            self.astra_to_hdf5()

    def astra_to_hdf5(self):
        astrabeamfilename = self.name + '.' + str(int(round((self.allElementObjects[self.end].position_end[2])*100))).zfill(4) + '.' + str(master_run_no).zfill(3)
        beam.read_astra_beam_file(master_subdir + '/' + astrabeamfilename, normaliseZ=False)
        HDF5filename = self.allElementObjects[self.end].name+'.hdf5'
        beam.write_HDF5_beam_file(master_subdir + '/' + HDF5filename, centered=False, sourcefilename=astrabeamfilename)


class frameworkCounter(dict):
    def __init__(self, sub={}):
        super(frameworkCounter, self).__init__()
        self.sub = sub

    def add(self, type):
        type = self.sub[type] if type in self.sub else type
        if type not in self:
            self[type] = 1
        else:
            self[type] += 1
        return self[type]

    def subtract(self, type):
        type = self.sub[type] if type in self.sub else type
        if type not in self:
            self[type] = 0
        else:
            self[type] = self[type] - 1 if self[type] > 0 else 0
        return self[type]

class frameworkObject(object):

    defaults = {'length': 0}

    def __init__(self, name=None, type=None, **kwargs):
        super(frameworkObject, self).__init__()
        if name == None:
            raise NameError('Command does not have a name')
        if type == None:
            raise NameError('Command does not have a type')
        self.name = name
        self.type = type
        self.commandtype = 'command'
        self.properties = {}
        self.properties['name'] = self.name
        self.properties['type'] = self.type
        if type in commandkeywords:
            self.allowedKeyWords = commandkeywords[self.type]['keywords'] + commandkeywords['common']['keywords']
        elif type in elementkeywords:
            self.allowedKeyWords = elementkeywords[self.type]['keywords'] + elementkeywords['common']['keywords']
            if 'framework_keywords' in  elementkeywords[self.type]:
                 self.allowedKeyWords += elementkeywords[self.type]['framework_keywords']
        else:
            print 'Unknown type = ', type
            raise NameError
        self.allowedKeyWords = map(lambda x: x.lower(), self.allowedKeyWords)
        for key, value in kwargs.iteritems():
            self.add_property(key, value)

    def add_property(self, key, value):
        if key.lower() in self.allowedKeyWords:
            try:
                self.properties[key.lower()] = value
                self.__setattr__(key.lower(), value)
            except:
                print key.lower()

    def add_default(self, key, value):
        self.defaults[key.lower()] = value

    @property
    def parameters(self):
        return self.properties

    def __getattr__(self, key):
        if key.lower() in self.defaults:
            return self.defaults[key.lower()]
        else:
            return None

    def __getitem__(self, key):
        return None

    def __setitem__(self,key,value):
        self.properties.__setitem__(key.lower(),value)

    def long(self, value):
        return int(value)

    def double(self, value):
        return float(value)

    def string(self, value):
        return ''+value+''

class frameworkGenerator(object):
    def __init__(self, executables, **kwargs):
        super(frameworkGenerator, self).__init__()
        self.executables = executables
        self.defaults = {}
        self.properties = {}
        self.allowedKeyWords = astra_generator_keywords['keywords'] + astra_generator_keywords['framework_keywords']
        self.allowedKeyWords = map(lambda x: x.lower(), self.allowedKeyWords)
        for key, value in kwargs.iteritems():
            if key.lower() in self.allowedKeyWords:
                try:
                    self.properties[key.lower()] = value
                    self.__setattr__(key.lower(), value)
                except:
                    print 'WARNING: Unknown keyword: ', key.lower()

    def run(self):
        command = self.executables['generator'] + [self.name+'.in']
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
            param_string = k+' = '+str(val)+',\n'
            if len((output + param_string).splitlines()[-1]) > 70:
                output += '\n'
            output += param_string
        return output[:-2]

    @property
    def charge(self):
        return self.properties['charge'] if 'charge' in self.properties and self.properties['charge'] is not None else 250e-12

    @property
    def name(self):
        return self.properties['name'] if 'name' in self.properties and self.properties['name'] is not None else 'laser'

    def write(self):
        output = '&INPUT\n'
        try:
            npart = eval(self.number_of_particles)
        except:
            npart = self.number_of_particles
        if self.filename is None:
            self.filename = 'laser.astra'
        framework_dict = OrderedDict([
            ['FName', {'value': self.filename, 'default': 'laser.astra'}],
            ['q_total', {'value': self.charge*1e9, 'default': 0.25}],
            ['Ipart', {'value': npart, 'default': 2**(3*3)}],
        ])
        keyword_dict = OrderedDict()
        for k in astra_generator_keywords['keywords']:
            if getattr(self, k.lower()) is not None:
                try:
                    val = eval(getattr(self,k.lower()))
                except:
                    val = getattr(self,k.lower())
                keyword_dict[k.lower()] = {'value': val}
        output += self._write_ASTRA(merge_two_dicts(framework_dict, keyword_dict))
        output += '\n/\n'
        saveFile(master_subdir+'/'+self.name+'.in', output)

    @property
    def parameters(self):
        return self.properties

    def add_property(self, key, value):
        if key.lower() in self.allowedKeyWords:
            self.properties[key.lower()] = value
            self.__setattr__(key.lower(), value)
        print self.properties

    def __getattr__(self, a):
        return None

    def astra_to_hdf5(self):
        astrabeamfilename = self.filename
        beam.read_astra_beam_file(master_subdir + '/' + astrabeamfilename, normaliseZ=False)
        HDF5filename = self.filename.replace('.astra','.hdf5')
        beam.write_HDF5_beam_file(master_subdir + '/' + HDF5filename, centered=False, sourcefilename=astrabeamfilename)

class frameworkCommand(frameworkObject):

    def __init__(self, name=None, type=None, **kwargs):
        super(frameworkCommand, self).__init__(name, type, **kwargs)
        if not type in commandkeywords:
            raise NameError('Command \'%s\' does not exist' % commandname)

    def write(self):
        wholestring=''
        string = '&'+self.type+'\n'
        for key, value in self.properties.iteritems():
            if not key =='name' and not key == 'type':
                string+='\t'+key+' = '+str(value)+'\n'
        string+='&end\n'
        return string

class frameworkElement(frameworkObject):

    def __init__(self, name=None, type=None, **kwargs):
        super(frameworkElement, self).__init__(name, type, **kwargs)

    def __mul__(self, other):
        return [self.properties for x in range(other)]

    def __rmul__(self, other):
        return [self.properties for x in range(other)]

    def __neg__(self):
        return self

    @property
    def theta(self):
        if hasattr(self, 'global_rotation') and self.global_rotation is not None:
            return self.global_rotation[2] if len(self.global_rotation) is 3 else self.global_rotation
        else:
            return 0

    def _rotation_matrix(self, theta):
        return np.array([[np.cos(theta), 0, np.sin(theta)], [0, 1, 0], [-1*np.sin(theta), 0, np.cos(theta)]])

    @property
    def rotation_matrix(self):
        return self._rotation_matrix(self.theta)

    def rotated_position(self, pos=[0,0,0], offset=[0,0,0]):
        return chop(np.dot(np.array(pos) - np.array(offset), self.rotation_matrix), 1e-6)

    @property
    def start(self):
        return self.position_start

    @property
    def middle(self):
        return np.array(self.position_start) + self.rotated_position([0,0,self.length / 2.0])

    @property
    def end(self):
        return self.position_end

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
                if self.isevaluable(s.group(1)) is True:
                    replaced_str = eval(re.sub(regex, eval(s.group(1)), param))
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
            d['value'] = self.expand_substitution(d['value'])
            return d['value'] if d['value'] is not None else d['default'] if 'default' in d else default
        elif isinstance(d, str):
            return getattr(self, d) if hasattr(self, d) and getattr(self, d) is not None else default

    def _write_ASTRA(self, d, n=1):
        output = ''
        for k, v in d.iteritems():
            if self.checkValue(v) is not None:
                if 'type' in v and v['type'] == 'list':
                    for i, l in enumerate(self.checkValue(v)):
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
                    for i, l in enumerate(self.checkValue(v)):
                        param_string += str(l)+', '
                        if len((output + param_string).splitlines()[-1]) > 70:
                            output += '\n'
                    output += param_string[:-2] + '),\n'
                else:
                    if n is not None:
                        param_string = k+'('+str(n)+') = '+str(self.checkValue(v))+', '
                    else:
                        param_string = k+' = '+str(self.checkValue(v))+',\n'
                    if len((output + param_string).splitlines()[-1]) > 70:
                        output += '\n'
                    output += param_string
        return output[:-2]

    def write_ASTRA(self, n):
        return ""

    def _write_Elegant(self):
        wholestring=''
        etype = self._convertType_Elegant(self.type)
        string = self.name+': '+ etype
        for key, value in self.properties.iteritems():
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

class dipole(frameworkElement):

    def __init__(self, name=None, type=None, **kwargs):
        super(dipole, self).__init__(name, type, **kwargs)

    @property
    def width(self):
        if 'width' in self.properties:
            return self.properties['width']
        else:
            return 0.2

    def __neg__(self):
        newself = copy.deepcopy(self)
        if 'exit_edge_angle' in newself.properties and 'entrance_edge_angle' in newself.properties:
            e1 = newself['entrance_edge_angle']
            e2 = newself['exit_edge_angle']
            newself.properties['entrance_edge_angle'] = e2
            newself.properties['exit_edge_angle'] = e1
        elif 'entrance_edge_angle' in newself.properties:
            newself.properties['exit_edge_angle'] = newself.properties['entrance_edge_angle']
            del newself.properties['entrance_edge_angle']
        elif 'exit_edge_angle' in newself.properties:
            newself.properties['entrance_edge_angle'] = newself.properties['exit_edge_angle']
            del newself.properties['exit_edge_angle']
        newself.name = '-'+newself.name
        return newself

    def _edge_angles(self, estr):
        if estr in self.properties:
            if isinstance(self.properties[estr], str):
                return self.checkValue(self.properties[estr],0)
            else:
                return self.properties[estr]
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
        theta = -self.theta-self.e1
        corners[0] = np.array(map(add,np.transpose(self.start), np.dot([-self.width*self.length,0,0], self._rotation_matrix(theta))))
        corners[3] = np.array(map(add,np.transpose(self.start), np.dot([self.width*self.length,0,0], self._rotation_matrix(theta))))
        theta = -self.theta-self.angle+self.e2
        corners[1] = np.array(map(add,np.transpose(self.end), np.dot([-self.width*self.length,0,0], self._rotation_matrix(theta))))
        corners[2] = np.array(map(add,np.transpose(self.end), np.dot([self.width*self.length,0,0], self._rotation_matrix(theta))))
        return corners

    def write_ASTRA(self, n):
        if abs(self.checkValue('strength', default=0)) > 0 or abs(self.rho) > 0:
            corners = self.corners
            if self.plane is None:
                self.plane = 'horizontal'
            params = OrderedDict([
                ['D_Type', {'value': '\''+self.plane+'\'', 'default': '\'horizontal\''}],
                ['D_Gap', {'type': 'list', 'value': self.gap, 'default': [0.02, 0.02]}],
                ['D1', {'type': 'array', 'value': [corners[3][0],corners[3][2]] }],
                ['D3', {'type': 'array', 'value': [corners[2][0],corners[2][2]] }],
                ['D4', {'type': 'array', 'value': [corners[1][0],corners[1][2]] }],
                ['D2', {'type': 'array', 'value': [corners[0][0],corners[0][2]] }],
                ])
            if abs(self.checkValue('strength', default=0)) > 0 or not abs(self.rho) > 0:
                params['D_strength'] = {'value': self.checkValue('strength',0), 'default': 1e6}
            else:
                params['D_radius'] =  {'value': -1*self.rho, 'default': 1e6}
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

class cavity(frameworkElement):

    def __init__(self, name=None, type=None, **kwargs):
        super(cavity, self).__init__(name, type, **kwargs)

    def write_ASTRA(self, n):
        if (self.n_cells is 0 or self.n_cells is None) and self.cell_length > 0:
                cells = round((self.length-self.cell_length)/self.cell_length)
                cells = cells - (cells % 3)
        elif self.n_cells > 0 and self.cell_length > 0:
            cells = self.n_cells - (self.n_cells % 3)
        else:
            cells = None
        return self._write_ASTRA(OrderedDict([
            ['C_pos', {'value': self.start[2], 'default': 0}],
            ['FILE_EFieLD', {'value': '\''+self.expand_substitution('\''+self.field_definition+'\'')+'\'', 'default': 0}],
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
            ['FILE_BFieLD', {'value': '\''+self.expand_substitution('\''+self.field_definition+'\'')+'\''}],
            ['MaxB', {'value': self.field_amplitude, 'default': 0}],
            ['S_smooth', {'value': self.smooth, 'default': 10}],
        ]), n)

class aperture(frameworkElement):

    def __init__(self, name=None, type=None, **kwargs):
        super(aperture, self).__init__(name, type, **kwargs)

class wall_current_monitor(frameworkElement):

    def __init__(self, name=None, type=None, **kwargs):
        super(wall_current_monitor, self).__init__(name, type, **kwargs)

class screen(frameworkElement):

    def __init__(self, name=None, type=None, **kwargs):
        super(screen, self).__init__(name, type, **kwargs)
        if 'output_filename' not in kwargs:
            self.properties['output_filename'] = self.name

    def write_ASTRA(self, n):
        return self._write_ASTRA(OrderedDict([
            ['Screen', {'value': self.middle[2], 'default': 0}],
        ]), n)

    def astra_to_hdf5(self, lattice):
        astrabeamfilename = lattice + '.' + str(int(round((self.middle[2])*100))).zfill(4) + '.' + str(master_run_no).zfill(3)
        beam.read_astra_beam_file(master_subdir + '/' + astrabeamfilename, normaliseZ=False)
        HDF5filename = self.name+'.hdf5'
        beam.write_HDF5_beam_file(master_subdir + '/' + HDF5filename, centered=False, sourcefilename=astrabeamfilename)

    def sdds_to_hdf5(self):
        elegantbeamfilename = self.output_filename+'.sdds'
        beam.read_SDDS_beam_file(master_subdir + '/' + elegantbeamfilename)
        print 'total charge = ', beam.beam['total_charge']
        HDF5filename = self.output_filename+'.hdf5'
        beam.write_HDF5_beam_file(master_subdir + '/' + HDF5filename, centered=False, sourcefilename=elegantbeamfilename)

class monitor(screen):

    def __init__(self, name=None, type=None, **kwargs):
        super(monitor, self).__init__(name, type, **kwargs)

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

class marker(frameworkElement):

    def __init__(self, name=None, type=None, **kwargs):
        super(marker, self).__init__(name, type, **kwargs)

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
        for k in elementkeywords[self.type]['keywords']:
            if getattr(self, k.lower()) is not None:
                keyword_dict[k.lower()] = {'value': getattr(self,k.lower())}
        return self._write_ASTRA(merge_two_dicts(self.framework_dict(), keyword_dict), None)

class astra_newrun(astra_header):
    def __init__(self, name=None, type=None, **kwargs):
        super(astra_header, self).__init__(name, type, **kwargs)
        if 'run' not in kwargs:
            self.properties['run'] = 1
        if 'head' not in kwargs:
            self.properties['head'] = 'trial'
        if 'lprompt' not in kwargs:
            self.add_property('lprompt', False)

    def framework_dict(self):
        return OrderedDict([
            ['Distribution', {'value': self.particle_definition}],
        ])

    def hdf5_to_astra(self):
        HDF5filename = self.particle_definition.replace('.astra','')+'.hdf5'
        beam.read_HDF5_beam_file(master_subdir + '/' + HDF5filename)
        astrabeamfilename = self.particle_definition
        beam.write_astra_beam_file(master_subdir + '/' + astrabeamfilename, normaliseZ=False)

class astra_output(astra_header):
    def __init__(self, name=None, type=None, **kwargs):
        super(astra_header, self).__init__(name, type, **kwargs)

    def framework_dict(self):
        return OrderedDict([
            ['zstart', {'value': self.start_element}],
            ['zstop', {'value': self.end_element}],
        ])

class astra_charge(astra_header):
    def __init__(self, name=None, type=None, **kwargs):
        super(astra_header, self).__init__(name, type, **kwargs)

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
