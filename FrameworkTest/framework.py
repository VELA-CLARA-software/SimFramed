import time, os, subprocess
import yaml
import traceback
import itertools
import copy
from collections import OrderedDict
from ASTRARules import *
from SimulationFramework.FrameworkHelperFunctions import *

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

def merge_two_dicts(x, y):
    z = x.copy()   # start with x's keys and values
    z.update(y)    # modifies z with y's keys and values & returns None
    return OrderedDict(z)

class framework(object):

    def __init__(self, directory, master_lattice_location=None):
        super(framework, self).__init__()
        self.directory = directory
        self.elementObjects = OrderedDict()
        self.lineObjects = {}
        self.commandObjects = OrderedDict()

        if master_lattice_location is None:
            self.master_lattice_location = (os.path.relpath(os.path.dirname(os.path.abspath(__file__)) + '/../MasterLattice/')+'/').replace('\\','/')
        else:
            self.master_lattice_location = master_lattice_location

    def loadElementsFile(self, input):
        if isinstance(input,(list,tuple)):
            filename = input
        else:
            filename = [input]
        for f in filename:
            with file(self.master_lattice_location + f, 'r') as stream:
                elements = yaml.load(stream)['elements']
            for name, elem in elements.iteritems():
                self.readElement(name, elem)

    def loadSettings(self, filename='short_240.settings'):
        """Load Lattice Settings from file"""
        stream = file(self.master_lattice_location+filename, 'r')
        settings = yaml.load(stream)
        self.globalSettings = settings['global']
        self.generatorFile = self.master_lattice_location + self.globalSettings['generatorFile'] if 'generatorFile' in self.globalSettings else None
        self.fileSettings = settings['files']
        elements = settings['elements']
        self.groups = settings['groups']
        stream.close()
        for name, elem in elements.iteritems():
            self.readElement(name, elem)

    def readElement(self, name, element, subelement=False):
        if name == 'filename':
            self.loadElementsFile(element)
        else:
            if subelement:
                self.addElement(name, subelement=True, **element)
            else:
                self.addElement(name, **element)
            if 'sub_elements' in element:
                for name, elem in element['sub_elements'].iteritems():
                    self.readElement(name, elem, True)

    def __getitem__(self,key):
        if key in self.elementObjects:
            return self.elementObjects[key].properties
        elif key in self.lineObjects:
            return self.lineObjects[key].line
        elif hasattr(self, key):
            return getattr(self,key.lower())

    @property
    def elements(self):
        return self.elementObjects.keys()

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
    def commands(self):
        return self.commandObjects.keys()

    def addElement(self, name=None, type=None, **kwargs):
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

    def getElement(self, element):
        if element in self.elementObjects:
            return self.elementObjects[element]
        else:
            print 'WARNING: Element ', element,' does not exist'
            return {}

    def getElementType(self, type):
        return [self.elementObjects[element] for element in self.elementObjects if self.elementObjects[element].type.lower() == type.lower()]

    def writeElements_ASTRA(self, file=None):
        counter = frameworkCounter()
        for t in ['cavities', 'solenoids', 'quadrupoles']:
            for element in getattr(self, t):
                print element.write_ASTRA(counter.number(element.type))

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

    @property
    def parameters(self):
        return self.properties

    def __getattr__(self, key):
        return None

    def __getitem__(self, key):
        # if key in self.properties:
        #     return self.properties.__getitem__(key.lower())
        # else:
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

    def rotated_position(self, pos=[0,0,0], offset=[0,0,0]):
        rotation_matrix = np.array([[np.cos(self.theta), 0, np.sin(self.theta)], [0, 1, 0], [-1*np.sin(self.theta), 0, np.cos(self.theta)]])
        return chop(np.dot(np.array(pos) - np.array(offset), rotation_matrix), 1e-6)

    @property
    def start(self):
        return self.position_start

    @property
    def middle(self):
        return np.array(self.position_start) + self.rotated_position([0,0,self.length / 2.0])

    def checkValue(self, d):
        return d['value'] if d['value'] is not None else d['default'] if 'default' in d else None

    def _write_ASTRA(self, d, n):
        output = ''
        for k, v in d.iteritems():
            if self.checkValue(v) is not None:
                if 'type' in v and v['type'] == 'list' or isinstance(self.checkValue(v), (list, tuple)):
                    print v['value']
                    for i, l in enumerate(self.checkValue(v)):
                            param_string = k+'('+str(n)+','+str(i+1)+') = '+str(l)+', '
                            if len((output + param_string).splitlines()[-1]) > 60:
                                output += '\n'
                            output += param_string
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

    def __neg__(self):
        newself = copy.deepcopy(self)
        # print newself.properties
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

class quadrupole(frameworkElement):

    def __init__(self, name=None, type=None, **kwargs):
        super(quadrupole, self).__init__(name, type, **kwargs)

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
