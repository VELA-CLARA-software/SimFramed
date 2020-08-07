import os
import subprocess
from ruamel import yaml
from munch import Munch, unmunchify
from collections import OrderedDict
from SimulationFramework.Modules.merge_two_dicts import merge_two_dicts
from SimulationFramework.FrameworkHelperFunctions import *
from SimulationFramework.FrameworkHelperFunctions import _rotation_matrix
import numpy as np

with open(os.path.dirname( os.path.abspath(__file__))+'/Codes/type_conversion_rules.yaml', 'r') as infile:
    type_conversion_rules = yaml.load(infile, Loader=yaml.UnsafeLoader)
    type_conversion_rules_Elegant = type_conversion_rules['elegant']
    type_conversion_rules_Names = type_conversion_rules['name']

with open(os.path.dirname( os.path.abspath(__file__))+'/Codes/Elegant/commands_Elegant.yaml', 'r') as infile:
    commandkeywords = yaml.load(infile, Loader=yaml.UnsafeLoader)

with open(os.path.dirname( os.path.abspath(__file__))+'/Codes/Elegant/elementkeywords.yaml', 'r') as infile:
    elementkeywords = yaml.load(infile, Loader=yaml.UnsafeLoader)

with open(os.path.dirname( os.path.abspath(__file__))+'/Codes/Elegant/keyword_conversion_rules_elegant.yaml', 'r') as infile:
    keyword_conversion_rules_elegant = yaml.load(infile, Loader=yaml.UnsafeLoader)

with open(os.path.dirname( os.path.abspath(__file__))+'/Codes/Elegant/elements_Elegant.yaml', 'r') as infile:
    elements_Elegant = yaml.load(infile, Loader=yaml.UnsafeLoader)

class frameworkLattice(Munch):
    def __init__(self, name, file_block, elementObjects, groupObjects, settings, executables, global_parameters):
        super(frameworkLattice, self).__init__()
        self.global_parameters = global_parameters
        self.objectname = name
        for key, value in list(elementObjects.items()):
            setattr(self, key, value)
        self.allElementObjects = elementObjects
        self.groupObjects = groupObjects
        self.allElements = list(self.allElementObjects.keys())
        self.file_block = file_block
        self.settings = settings
        self.globalSettings = settings['global']
        self.groupSettings = file_block['groups'] if 'groups' in file_block else {}
        self.update_groups()
        self.executables = executables
        self._csr_enable = True
        self.csrDrifts = True
        self.lscDrifts = True
        self.lsc_bins = 20
        self.lsc_high_frequency_cutoff_start = -1
        self.lsc_high_frequency_cutoff_end = -1
        self.lsc_low_frequency_cutoff_start = -1
        self.lsc_low_frequency_cutoff_end = -1
        self._sample_interval = self.file_block['input']['sample_interval'] if 'input' in self.file_block and 'sample_interval' in self.file_block['input'] else 1

    def insert_element(self, index, element):
        for i, _ in enumerate(range(len(self.elements))):
            k, v = self.elements.popitem(False)
            self.elements[element.objectname if i == index else k] = element

    @property
    def csr_enable(self):
        return self._csr_enable
    @csr_enable.setter
    def csr_enable(self, csr):
        self.csrDrifts = csr
        self._csr_enable = csr

    @property
    def sample_interval(self):
        return self._sample_interval
    @sample_interval.setter
    def sample_interval(self, interval):
        # print('Setting new sample_interval = ', interval)
        self._sample_interval = interval

    @property
    def prefix(self):
        if 'input' not in self.file_block:
            self.file_block['input'] = {}
        if 'prefix' not in self.file_block['input']:
            self.file_block['input']['prefix'] = ''
        return self.file_block['input']['prefix']
    @prefix.setter
    def prefix(self, prefix):
        self.file_block['input']['prefix'] = prefix

    def update_groups(self):
        for g in list(self.groupSettings.keys()):
            if g in self.groupObjects:
                self.groupObjects[g].update(**self.groupSettings[g])

    def getElement(self, element, param=None):
        if element in self.allElements:
            if param is not None:
                return getattr(self.allElementObjects[element], param.lower())
            else:
                return self.allElements[element]
        elif element in list(self.groupObjects.keys()):
            if param is not None:
                return getattr(self.groupObjects[element], param.lower())
            else:
                return self.groupObjects[element]
        else:
            print(( 'WARNING: Element ', element,' does not exist'))
            return {}

    def getElementType(self, type, setting=None):
        return [self.elements[element] if setting is None else self.elements[element][setting] for element in list(self.elements.keys()) if self.elements[element].objecttype.lower() == type.lower()]

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
    def wakefields(self):
        return self.getElementType('longitudinal_wakefield')

    @property
    def screens(self):
        return self.getElementType('screen')

    @property
    def screens_and_bpms(self):
        return sorted(self.getElementType('screen') + self.getElementType('beam_position_monitor'), key=lambda x: x.position_start[2])

    @property
    def apertures(self):
        return sorted(self.getElementType('aperture') + self.getElementType('collimator'), key=lambda x: x.position_start[2])

    @property
    def lines(self):
        return list(self.lineObjects.keys())

    @property
    def start(self):
        if 'start_element' in self.file_block['output']:
            return self.file_block['output']['start_element']
        elif 'zstart' in self.file_block['output']:
            for e in list(self.allElementObjects.keys()):
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
            for e in list(self.allElementObjects.keys()):
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

    @property
    def elements(self):
        index_start = self.allElements.index(self.start)
        index_end = self.allElements.index(self.end)
        f = OrderedDict([[e,self.allElementObjects[e]] for e in self.allElements[index_start:index_end+1]])
        return f

    def write(self):
        pass

    def run(self):
        """Run the code with input 'filename'"""
        command = self.executables[self.code] + [self.objectname]
        with open(os.path.relpath(self.global_parameters['master_subdir']+'/'+self.objectname+'.log', '.'), "w") as f:
            subprocess.call(command, stdout=f, cwd=self.global_parameters['master_subdir'])

    def postProcess(self):
        pass

    def preProcess(self):
        pass

    def __repr__(self):
        return self.elements

    def __str__(self):
        str = self.objectname + ' = ('
        for e in self.elements:
            if len((str + e).splitlines()[-1]) > 60:
                str += '&\n'
            str += e+', '
        return str + ')'

class frameworkObject(Munch):

    def __init__(self, objectname=None, objecttype=None, **kwargs):
        super(frameworkObject, self).__init__()
        if 'global_parameters' in kwargs:
            self.global_parameters = kwargs['global_parameters']
        if objectname == None:
            raise NameError('Command does not have a name')
        if objecttype == None:
            raise NameError('Command does not have a type')
        setattr(self, 'objectdefaults', OrderedDict())
        setattr(self, 'objectname', objectname)
        setattr(self, 'objecttype', objecttype)
        if self.objecttype in commandkeywords:
            self.allowedkeywords = commandkeywords[self.objecttype]
        elif self.objecttype in elementkeywords:
            self.allowedkeywords = merge_two_dicts(elementkeywords[self.objecttype]['keywords'], elementkeywords['common']['keywords'])
            if 'framework_keywords' in  elementkeywords[self.objecttype]:
                 self.allowedkeywords = merge_two_dicts(self.allowedkeywords, elementkeywords[self.objecttype]['framework_keywords'])
        else:
            print(( 'Unknown type = ', objecttype))
            raise NameError
        self.allowedkeywords = [x.lower() for x in self.allowedkeywords]
        for key, value in list(kwargs.items()):
            self.add_property(key, value)

    def add_property(self, key, value):
        key = key.lower()
        if key in self.allowedkeywords:
            try:
                setattr(self, key, value)
            except Exception as e:
                print((self.objecttype,'[', key, ']: ', e))

    def add_default(self, key, value):
        self.objectdefaults[key] = value

    @property
    def parameters(self):
        return list(self.keys())

    @property
    def objectproperties(self):
        return self

    def __getitem__(self, key):
        lkey = key.lower()
        defaults = super(frameworkObject, self).__getitem__('objectdefaults')
        if lkey in defaults:
            try:
                return super(frameworkObject, self).__getitem__(lkey)
            except:
                return defaults[lkey]
        else:
            try:
                return super(frameworkObject, self).__getitem__(lkey)
            except:
                try:
                    return super(frameworkObject, self).__getattribute__(key)
                except:
                    return None

    def __repr__(self):
        string = ''
        for k,v in list(self.items()):
            string+="{} ({})".format(k, v)+'\n'
        return string

class frameworkCommand(frameworkObject):

    def __init__(self, name=None, type=None, **kwargs):
        super(frameworkCommand, self).__init__(name, type, **kwargs)
        if not type in commandkeywords:
            raise NameError('Command \'%s\' does not exist' % commandname)

    def write_Elegant(self):
        wholestring=''
        string = '&'+self.objecttype+'\n'
        # print(self.objecttype, self.objectproperties)
        for key in commandkeywords[self.objecttype]:
            if key.lower() in self.objectproperties and not key =='name' and not key == 'type' and not self.objectproperties[key.lower()] is None:
                string+='\t'+key+' = '+str(self.objectproperties[key.lower()])+'\n'
        string+='&end\n'
        return string

    def write_MAD8(self):
        wholestring=''
        string = self.objecttype
        # print(self.objecttype, self.objectproperties)
        for key in commandkeywords[self.objecttype]:
            if key.lower() in self.objectproperties and not key =='name' and not key == 'type' and not self.objectproperties[key.lower()] is None:
                e = ',' + key +'=' + str(self.objectproperties[key.lower()])
                if len((string + e).splitlines()[-1]) > 79:
                    string += ',&\n'
                string += e
        string+=';\n'
        return string

class frameworkGroup(object):
    def __init__(self, name, elementObjects, type, elements, **kwargs):
        super(frameworkGroup, self).__init__()
        self.objectname = name
        self.type = type
        self.elements = elements
        self.allElementObjects = elementObjects

    def get_Parameter(self, p):
        try:
            isinstance(type(self).p, p)
            return getattr(self, p)
        except:
            return self.allElementObjects[self.elements[0]][p]

    def change_Parameter(self, p, v):
        # print 'p = ', getattr(self, p)
        try:
            getattr(self, p)
            setattr(self, p, v)
            if p == 'angle':
                self.set_angle(v)
            # print ('Changing group ', self.objectname, ' ', p, ' = ', v, '  result = ', self.get_Parameter(p))
        except:
            for e in self.elements:
                setattr(self.allElementObjects[e], p, v)
                # print ('Changing group elements ', self.objectname, ' ', p, ' = ', v, '  result = ', self.allElementObjects[self.elements[0]].objectname, self.get_Parameter(p))


    # def __getattr__(self, p):
    #     return self.get_Parameter(p)

    def __repr__(self):
        return [self.allElementObjects[e].objectname for e in self.elements]

    def __str__(self):
        return str([self.allElementObjects[e].objectname for e in self.elements])

    def __getitem__(self, key):
        return self.get_Parameter(key)

class element_group(frameworkGroup):
    def __init__(self, name, elementObjects, type, elements, **kwargs):
        super(element_group, self).__init__(name, elementObjects, type, elements, **kwargs)

class chicane(frameworkGroup):
    def __init__(self, name, elementObjects, type, elements, **kwargs):
        super(chicane, self).__init__(name, elementObjects, type, elements, **kwargs)
        self.ratios = (1,-1,-1,1)

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
        return float(obj[0].angle)
    @angle.setter
    def angle(self, theta):
        'using setter! angle = ', theta
        self.set_angle(theta)

    def set_angle(self, a):
        indices = list(sorted([list(self.allElementObjects).index(e) for e in self.elements]))
        dipole_objs = [self.allElementObjects[e] for e in self.elements]
        obj = [self.allElementObjects[list(self.allElementObjects)[e]] for e in range(indices[0],indices[-1]+1)]
        starting_angle = obj[0].theta
        dipole_number = 0
        for i in range(len(obj)):
            x1 = np.transpose([obj[i].position_start])
            obj[i].global_rotation[2] = starting_angle
            if obj[i] in dipole_objs:
                obj[i].angle = a*self.ratios[dipole_number]
                dipole_number += 1
                elem_angle = obj[i].angle
            else:
                elem_angle = obj[i].angle if obj[i].angle is not None else 0
            obj[i].position_end = list(obj[i].end)
            xstart, ystart, zstart = obj[i].position_end
            if i < len(obj)-1:
                xend, yend, zend = obj[i+1].position_start
                angle = starting_angle + elem_angle
                length = float((-zstart + zend))
                endx = chop(float(xstart + np.tan(angle)*length))
                obj[i+1].position_start[0] =  endx
                obj[i+1].global_rotation[2] =  angle
                starting_angle += elem_angle

    def __str__(self):
        return str([[self.allElementObjects[e].objectname, self.allElementObjects[e].angle, self.allElementObjects[e].global_rotation[2], self.allElementObjects[e].position_start[0], self.allElementObjects[e].position_end[0]] for e in self.elements])

class s_chicane(chicane):
    def __init__(self, name, elementObjects, type, elements, **kwargs):
        super(s_chicane, self).__init__(name, elementObjects, type, elements, **kwargs)
        self.ratios = (-1,2,-2,1)

class frameworkCounter(dict):
    def __init__(self, sub={}):
        super(frameworkCounter, self).__init__()
        self.sub = sub

    def counter(self, type):
        type = self.sub[type] if type in self.sub else type
        if type not in self:
            return 1
        return self[type] + 1

    def value(self, type):
        type = self.sub[type] if type in self.sub else type
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

class frameworkElement(frameworkObject):

    def __init__(self, elementName=None, elementType=None, **kwargs):
        super(frameworkElement, self).__init__(elementName, elementType, **kwargs)
        self.add_default('length', 0)
        self.add_property('position_errors', [0,0,0])
        self.add_property('rotation_errors', [0,0,0])
        self.add_default('global_rotation', [0,0,0])
        self.add_default('starting_rotation', 0)
        self.keyword_conversion_rules_elegant = keyword_conversion_rules_elegant['general']
        if elementType in keyword_conversion_rules_elegant:
            self.keyword_conversion_rules_elegant = merge_two_dicts(self.keyword_conversion_rules_elegant, keyword_conversion_rules_elegant[elementType])

    def __mul__(self, other):
        return [self.objectproperties for x in range(other)]

    def __rmul__(self, other):
        return [self.objectproperties for x in range(other)]

    def __neg__(self):
        return self

    @property
    def x(self):
        return self.position_start[0]
    @x.setter
    def x(self, x):
        self.position_start[0] = x
        self.position_end[0] = x
    @property
    def y(self):
        return self.position_start[1]
    @y.setter
    def y(self, y):
        self.position_start[1] = y
        self.position_end[1] = y
    @property
    def z(self):
        return self.position_start[2]
    @z.setter
    def z(self, z):
        self.position_start[2] = z
        self.position_end[2] = z

    @property
    def dx(self):
        return self.position_errors[0]
    @dx.setter
    def dx(self, x):
        self.position_errors[0] = x
    @property
    def dy(self):
        return self.position_errors[1]
    @dy.setter
    def dy(self, y):
        self.position_errors[1] = y
    @property
    def dz(self):
        return self.position_errors[2]
    @dz.setter
    def dz(self, z):
        self.position_errors[2] = z

    @property
    def x_rot(self):
        return self.global_rotation[1]
    @property
    def y_rot(self):
        return self.global_rotation[2] + self.starting_rotation
    @property
    def z_rot(self):
        return self.global_rotation[0]

    @property
    def dx_rot(self):
        return self.rotation_errors[1]
    @dx_rot.setter
    def dx_rot(self, x):
        self.rotation_errors[1] = x
    @property
    def dy_rot(self):
        return self.rotation_errors[2]
    @dy_rot.setter
    def dy_rot(self, y):
        self.rotation_errors[2] = y
    @property
    def dz_rot(self):
        return self.rotation_errors[0]
    @dz_rot.setter
    def dz_rot(self, z):
        self.rotation_errors[0] = z
    @property
    def tilt(self):
        return self.dz_rot

    @property
    def PV(self):
        if hasattr(self, 'PV_root'):
            return self.PV_root
        else:
            return self.objectName

    @property
    def theta(self):
        if hasattr(self, 'global_rotation') and self.global_rotation is not None:
            rotation =  self.global_rotation[2] if len(self.global_rotation) is 3 else self.global_rotation
        else:
            rotation = 0
        # if hasattr(self, 'starting_rotation') and self.starting_rotation is not None:
        #     rotation +=  self.starting_rotation
        return rotation

    @property
    def rotation_matrix(self):
        return _rotation_matrix(self.theta)

    def rotated_position(self, pos=[0,0,0], offset=None, theta=None):
        if offset is None:
            if not hasattr(self, 'starting_offset') or self.starting_offset is None:
                offset = [0,0,0]
            else:
                offset = self.starting_offset
        if theta is None:
            return chop(np.dot(np.array(pos) - np.array(offset), self.rotation_matrix), 1e-6)
        else:
            return chop(np.dot(np.array(pos) - np.array(offset), _rotation_matrix(theta)), 1e-6)

    @property
    def start(self):
        start = np.array(self.position_start)
        return start

    @property
    def middle(self):
        start = np.array(self.position_start)
        return start + self.rotated_position(np.array([0,0, self.length / 2.0]), offset=self.starting_offset, theta=self.y_rot)

    @property
    def end(self):
        start = np.array(self.position_start)
        # print(self.objectName, start, start + self.rotated_position(np.array([0,0, self.length]), offset=self.starting_offset, theta=-self.y_rot))
        return start + self.rotated_position(np.array([0,0, self.length]), offset=self.starting_offset, theta=self.y_rot)

    def relative_position_from_centre(self, vec=[0,0,0]):
        start = np.array(self.start)
        return start + self.rotated_position(np.array([0,0, self.length / 2.0]) + np.array(vec), offset=self.starting_offset, theta=self.y_rot)
    def relative_position_from_start(self, vec=[0,0,0]):
        start = np.array(self.start)
        return start + self.rotated_position(np.array(vec), offset=self.starting_offset, theta=self.y_rot)

    def _write_ASTRA(self, d, n=1):
        output = ''
        for k, v in list(d.items()):
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
        etype = self._convertType_Elegant(self.objecttype)
        string = self.objectname+': '+ etype
        k1 = self.k1 if self.k1 is not None else 0
        for key, value in list(merge_two_dicts({'k1': k1}, merge_two_dicts(self.objectproperties, self.objectdefaults)).items()):
            if not key is 'name' and not key is 'type' and not key is 'commandtype' and self._convertKeword_Elegant(key) in elements_Elegant[etype]:
                value = getattr(self, key) if hasattr(self, key) and getattr(self, key) is not None else value
                key = self._convertKeword_Elegant(key)
                value = 1 if value is True else value
                value = 0 if value is False else value
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
        return self.keyword_conversion_rules_elegant[keyword] if keyword in self.keyword_conversion_rules_elegant else keyword

    def write_CSRTrack(self, n=0):
        return ""

    def write_GPT(self, Brho, ccs="wcs", *args, **kwargs):
        return ''

    def gpt_coordinates(self, position, rotation):
        x,y,z = chop(position, 1e-6)
        psi, phi, theta = rotation
        output =''
        for c in [-x, y, z]:
            output += str(c)+', '
        output += 'cos('+str(theta)+'), 0, -sin('+str(theta)+'), 0, 1 ,0'
        return output

    def gpt_ccs(self, ccs):
        return ccs


class getGrids(object):

    def __init__(self):
        self.powersof8 = np.asarray([ 2**(j) for j in range(1,20) ])

    def getGridSizes(self, x):
        self.x = abs(x)
        self.cuberoot = int(round(self.x ** (1. / 3)))
        return max([4,self.find_nearest(self.powersof8, self.cuberoot)])

    def find_nearest(self, array, value):
        self.array = array
        self.value = value
        self.idx = (np.abs(self.array - self.value)).argmin()
        return self.array[self.idx]
