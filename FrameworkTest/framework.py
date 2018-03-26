import time, os, subprocess, re
import yaml
import traceback
import itertools
import copy
from collections import OrderedDict
from ASTRARules import *
from SimulationFramework.FrameworkHelperFunctions import *
from operator import add

_mapping_tag = yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG

def dict_representer(dumper, data):
    return dumper.represent_dict(data.iteritems())

def dict_constructor(loader, node):
    return OrderedDict(loader.construct_pairs(node))

yaml.add_representer(OrderedDict, dict_representer)
yaml.add_constructor(_mapping_tag, dict_constructor)


commandkeywords = {}

with open(os.path.dirname( os.path.abspath(__file__))+'/elementkeywords.yaml', 'r') as infile:
    elementkeywords = yaml.load(infile)

with open(os.path.dirname( os.path.abspath(__file__))+'/elements_Elegant.yaml', 'r') as infile:
    elements_Elegant = yaml.load(infile)

type_conversion_rules_Elegant = {'dipole': 'csrcsbend', 'quadrupole': 'kquad', 'beam_position_monitor': 'moni', 'screen': 'watch', 'aperture': 'rcol',
                         'collimator': 'ecol', 'monitor': 'moni', 'solenoid': 'mapsolenoid', 'wall_current_monitor': 'charge', 'cavity': 'rfcw',
                         'rf_deflecting_cavity': 'rfdf'}
keyword_conversion_rules_Elegant = {'length': 'l','entrance_edge_angle': 'e1', 'exit_edge_angle': 'e2', 'edge_field_integral': 'fint', 'horizontal_size': 'x_max', 'vertical_size': 'y_max',
                            'field_amplitude': 'volt', 'frequency': 'freq'}

section_header_text_ASTRA = {'cavities': {'header': 'CAVITY', 'bool': 'LEField'},
                             'solenoids': {'header': 'SOLENOID', 'bool': 'LBField'},
                             'quadrupoles': {'header': 'QUADRUPOLE', 'bool': 'LQuad'},
                            }

def merge_two_dicts(x, y):
    z = x.copy()   # start with x's keys and values
    z.update(y)    # modifies z with y's keys and values & returns None
    return OrderedDict(z)

class framework(object):

    def __init__(self, directory, master_lattice=None):
        super(framework, self).__init__()
        global master_lattice_location
        self.directory = directory
        self.elementObjects = OrderedDict()
        self.latticeObjects = {}
        self.commandObjects = OrderedDict()

        if master_lattice is None:
            master_lattice_location = (os.path.relpath(os.path.dirname(os.path.abspath(__file__)) + '/../MasterLattice/')+'/').replace('\\','/')
        else:
            master_lattice_location = master_lattice

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
        stream = file(master_lattice_location+filename, 'r')
        settings = yaml.load(stream)
        self.globalSettings = settings['global']
        self.generatorFile = master_lattice_location + self.globalSettings['generatorFile'] if 'generatorFile' in self.globalSettings else None
        self.fileSettings = settings['files']
        elements = settings['elements']
        self.groups = settings['groups']
        stream.close()
        for name, elem in elements.iteritems():
            self.read_Element(name, elem)
        for name, elem in self.fileSettings.iteritems():
            self.read_Lattice(name, elem)

    def read_Lattice(self, name, elem):
        self.latticeObjects[name] = frameworkLattice(name, self.elementObjects, elem, self.globalSettings)

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
        return self.lineObjects.keys()

    @property
    def commands(self):
        return self.commandObjects.keys()

class frameworkLattice(object):
    def __init__(self, name, elementObjects, file_block, globalSettings):
        super(frameworkLattice, self).__init__()
        self.name = name
        self.allElementObjects = elementObjects
        self.allElements = self.allElementObjects.keys()
        self.file_block = file_block
        self.globalSettings = globalSettings

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
    def lines(self):
        return self.lineObjects.keys()

    @property
    def start(self):
        if 'start_element' in self.output_block:
            return self.output_block['start_element']
        elif 'zstart' in self.output_block:
            for e in self.allElementObjects.keys():
                if self.allElementObjects[e].position_start[2] == self.output_block['zstart']:
                    return e
        else:
            return self.allElementObjects[0]

    @property
    def end(self):
        if 'end_element' in self.output_block:
            return self.output_block['end_element']
        elif 'zstop' in self.output_block:
            endelems = []
            for e in self.allElementObjects.keys():
                if self.allElementObjects[e]['position_end'] == self.output_block['zstop']:
                    endelems.append(e)
                elif self.allElementObjects[e]['position_end'] > self.output_block['zstop'] and len(endelems) == 0:
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

    def writeElements_ASTRA(self, file=None):
        counter = frameworkCounter()
        newrun = astra_newrun('newrun', **merge_two_dicts(self.file_block['input'],self.globalSettings['ASTRAsettings']))
        output = astra_output('output', **merge_two_dicts(self.file_block['output'],self.globalSettings['ASTRAsettings']))
        print newrun.parameters
        print output.parameters
        exit()
        fulltext = ''
        for t in ['cavities', 'solenoids', 'quadrupoles']:
            fulltext = '&' + section_header_text_ASTRA[t]['header']+'\n'
            elements = getattr(self, t)
            fulltext += section_header_text_ASTRA[t]['bool']+' = '+str(len(elements) > 0)+'\n'
            for element in elements:
                fulltext += element.write_ASTRA(counter.number(element.type))+'\n'
            fulltext += '\n/\n'
            print fulltext

class frameworkCounter(dict):
    def __init__(self):
        super(frameworkCounter, self).__init__()

    def number(self, type):
        if type not in self:
            self[type] = 1
        else:
            self[type] += 1
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
        else:
            print 'Unknown type = ', type
            raise NameError
        self.allowedKeyWords = map(lambda x: x.lower(), self.allowedKeyWords)
        for key, value in kwargs.iteritems():
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
        if hasattr(self, 'global_rotation'):
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
                    replaced_str = re.sub(regex, eval(s.group(1)), param)
                else:
                    replaced_str = re.sub(regex, s.group(1), param)
                for key in subs:
                    replaced_str = replaced_str.replace(key, subs[key])
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

    def _write_ASTRA(self, d, n):
        output = ''
        for k, v in d.iteritems():
            if self.checkValue(v) is not None:
                if 'type' in v and v['type'] == 'list':
                    for i, l in enumerate(self.checkValue(v)):
                        param_string = k+'('+str(n)+','+str(i+1)+') = '+str(l)+', '
                        if len((output + param_string).splitlines()[-1]) > 60:
                            output += '\n'
                        output += param_string
                elif 'type' in v and v['type'] == 'array':
                    param_string = k+'('+str(n)+') = ('
                    for i, l in enumerate(self.checkValue(v)):
                        param_string += str(l)+', '
                        if len((output + param_string).splitlines()[-1]) > 60:
                            output += '\n'
                    output += param_string[:-2] + '), '
                else:
                    param_string = k+'('+str(n)+') = '+str(self.checkValue(v))+', '
                    if len((output + param_string).splitlines()[-1]) > 60:
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

    def write_ASTRA(self, n):
        corners = [0,0,0,0]
        theta = -self.theta-self.e1
        corners[0] = np.array(map(add,np.transpose(self.start), np.dot([-self.width*self.length,0,0], self._rotation_matrix(theta))))
        corners[3] = np.array(map(add,np.transpose(self.start), np.dot([self.width*self.length,0,0], self._rotation_matrix(theta))))
        theta = -self.theta-self.angle+self.e2
        corners[1] = np.array(map(add,np.transpose(self.end), np.dot([-self.width*self.length,0,0], self._rotation_matrix(theta))))
        corners[2] = np.array(map(add,np.transpose(self.end), np.dot([self.width*self.length,0,0], self._rotation_matrix(theta))))
        params = OrderedDict([
            ['D_Type', {'value': self.plane, 'default': 'horizontal'}],
            ['D_Gap', {'type': 'list', 'value': self.gap, 'default': [0.02, 0.02]}],
            ['D1', {'type': 'array', 'value': [corners[3][0],corners[3][2]] }],
            ['D3', {'type': 'array', 'value': [corners[2][0],corners[2][2]] }],
            ['D4', {'type': 'array', 'value': [corners[1][0],corners[1][2]] }],
            ['D2', {'type': 'array', 'value': [corners[0][0],corners[0][2]] }],
            ])
        if abs(self.checkValue('strength', default=0)) > 0:
            params['D_strength'] = {'value': self.checkValue('strength',0), 'default': 1e6}
        else:
            params['D_radius'] =  {'value': self.rho, 'default': 1e6}
        return self._write_ASTRA(params, n)

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
        return self._write_ASTRA(OrderedDict([
            ['C_pos', {'value': self.start[2], 'default': 0}],
            ['FILE_EFieLD', {'value': self.field_definition, 'default': 0}],
            ['Nue', {'value': self.frequency / 1e6, 'default': 2.9985}],
            ['MaxE', {'value': self.field_amplitude / 1e6, 'default': 0}],
            ['Phi', {'value': self.phase, 'default': 0.0}],
            ['C_smooth', {'value': self.smooth, 'default': 10}],
        ]), n)

class rf_deflecting_cavity(frameworkElement):

    def __init__(self, name=None, type=None, **kwargs):
        super(rf_deflecting_cavity, self).__init__(name, type, **kwargs)

class solenoid(frameworkElement):

    def __init__(self, name=None, type=None, **kwargs):
        super(solenoid, self).__init__(name, type, **kwargs)

    def write_ASTRA(self, n):
        return self._write_ASTRA(OrderedDict([
            ['S_pos', {'value': self.start[2], 'default': 0}],
            ['FILE_BFieLD', {'value': self.field_definition}],
            ['MaxB', {'value': self.field_amplitude / 1e6, 'default': 0}],
            ['S_smooth', {'value': self.smooth, 'default': 10}],
        ]), n)

class aperture(frameworkElement):

    def __init__(self, name=None, type=None, **kwargs):
        super(aperture, self).__init__(name, type, **kwargs)

class beam_position_monitor(frameworkElement):

    def __init__(self, name=None, type=None, **kwargs):
        super(beam_position_monitor, self).__init__(name, type, **kwargs)

class kicker(frameworkElement):

    def __init__(self, name=None, type=None, **kwargs):
        super(kicker, self).__init__(name, type, **kwargs)

class wall_current_monitor(frameworkElement):

    def __init__(self, name=None, type=None, **kwargs):
        super(wall_current_monitor, self).__init__(name, type, **kwargs)

class screen(frameworkElement):

    def __init__(self, name=None, type=None, **kwargs):
        super(screen, self).__init__(name, type, **kwargs)

class monitor(frameworkElement):

    def __init__(self, name=None, type=None, **kwargs):
        super(monitor, self).__init__(name, type, **kwargs)

class collimator(frameworkElement):

    def __init__(self, name=None, type=None, **kwargs):
        super(collimator, self).__init__(name, type, **kwargs)

class marker(frameworkElement):

    def __init__(self, name=None, type=None, **kwargs):
        super(marker, self).__init__(name, type, **kwargs)

class astra_newrun(frameworkObject):

    def __init__(self, name=None, **kwargs):
        super(astra_newrun, self).__init__(name, 'astra_newrun', **kwargs)

class astra_output(frameworkObject):

    def __init__(self, name=None, **kwargs):
        super(astra_output, self).__init__(name, 'astra_output', **kwargs)
