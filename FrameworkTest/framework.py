import time, os, subprocess
import yaml
import traceback
import itertools
import copy
from collections import OrderedDict
from ASTRARules import *
from SimulationFramework.FrameworkHelperFunctions import *

commandkeywords = {}

elementkeywords = {'common': {'keywords':['subelement', 'length', 'Online_Model_Name', 'Controller_Name', 'PV', 'global_rotation', 'position_end', 'position_start', 'buffer_start', 'buffer_end', 'buffer_start_length','buffer_end_length']},
                   'quadrupole': {'keywords':['k1','sr_enable','isr_enable']},
                   'aperture': {'keywords':['shape','horizontal_size', 'vertical_size']},
                   'cavity': {'keywords':['n_cells', 'Structure_Type', 'frequency', 'phase', 'field_definition', 'lsc_cutoff_high']},
                   'solenoid': {'keywords':['field_amplitude', 'field_definition']},
                   'kicker': {'keywords':['sr_enable', 'isr_enable', 'Horizontal_PV', 'Vertical_PV']},
                   'wall_current_monitor': {'keywords':[]},
                   'beam_position_monitor': {'keywords':[]},
                   'monitor': {'keywords':[]},
                   'screen': {'keywords':['output_filename']},
                   'collimator': {'keywords':['shape', 'horizontal_size', 'vertical_size']},
                   'marker': {'keywords':['fit_point']},
                   'dipole': {'keywords':['angle', 'entrance_edge_angle', 'exit_edge_angle', 'half_gap', 'edge_field_integral', 'csr_bins', 'sr_enable', 'isr_enable', 'csr_output_filename']},
                   'rf_deflecting_cavity': {'keywords':['phase']},
                   }

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
            stream = file(self.master_lattice_location + f, 'r')
            elements = yaml.load(stream)['elements']
            stream.close()
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
                subelements = [name]
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
        # if not type in elementkeywords:
        #     raise NameError('Element \'%s\' does not exist' % commandname)
        try:
            element = globals()[type](name, type, **kwargs)
            self.elementObjects[name] = element
            return element
        except:
            raise NameError('Element \'%s\' does not exist' % type)
        # element = frameworkElement(name, **kwargs)

    def getElement(self, element):
        if element in self.elementObjects:
            return self.elementObjects[element]
        else:
            print 'WARNING: Element ', element,' does not exist'
            return {}

    def addLine(self, name=None, line=[]):
        if name == None:
            raise NameError('Line does not have a name')
        line = elegantLine(name, line)
        # setattr(self,name,line)
        # setattr(self.parent,name,line)
        self.lineObjects[name] = line
        return line

    def writeElements(self, file, lattice):
        self._doLineExpansion(lattice)
        self.elementDefinitions = reduce(lambda l, x: l if x in l else l+[x], self.elementDefinitions, [])
        for element in self.elementDefinitions:
            file.write(self.getElement(element).write())

    def screensToWatch(self):
        self._doLineExpansion('cla-ebt')
        self.elementDefinitions = reduce(lambda l, x: l if x in l else l+[x], self.elementDefinitions, [])
        for element in self.elementDefinitions:
            if 'scr' in element.lower():
                print getattr(self,element).properties

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
                    setattr(self,key.lower(), value)
                except:
                    pass

    @property
    def parameters(self):
        return self.properties

    def __getitem__(self,key):
        return self.properties.__getitem__(key.lower())

    def __setitem__(self,key,value):
        self.properties.__setitem__(key.lower(),value)
        setattr(self,key.lower(),value)

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

    def write(self):
        wholestring=''
        string = self.name+': '+self.type
        for key, value in self.properties.iteritems():
            if not key is 'name' and not key is 'type' and not key is 'commandtype':
                tmpstring = ', '+key+' = '+str(value)
                if len(string+tmpstring) > 76:
                    wholestring+=string+',&\n'
                    string=''
                    string+=tmpstring[2::]
                else:
                    string+= tmpstring
        wholestring+=string+';\n'
        return wholestring

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
        if 'smooth' not in kwargs or smooth is None:
            self.smooth = 10
        if 'bore' not in kwargs or bore is None:
            self.bore = 0.016

    def write_ASTRA(self, n):
        quadtext = 'Q_K('+str(n)+')='+str(self.k1)+', Q_length('+str(n)+')='+str(self.length)+',\n'+\
        'Q_pos('+str(n)+')='+str(self.position_start[2]+(self.length/2.))+', Q_smooth('+str(n)+')='+str(self.smooth)+', Q_bore('+str(n)+')='+str(self.bore)+'\n'

        for var in ASTRARules['QUADRUPOLE']:
            soltext += createOptionalString(sol, var, n)

        return quadtext

class cavity(frameworkElement):

    def __init__(self, name=None, type=None, **kwargs):
        super(cavity, self).__init__(name, type, **kwargs)

class rf_deflecting_cavity(frameworkElement):

    def __init__(self, name=None, type=None, **kwargs):
        super(rf_deflecting_cavity, self).__init__(name, type, **kwargs)

class solenoid(frameworkElement):

    def __init__(self, name=None, type=None, **kwargs):
        super(solenoid, self).__init__(name, type, **kwargs)

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
