import time, os, subprocess, re
from ruamel import yaml
import copy
from collections import OrderedDict
from SimulationFramework.Modules.merge_two_dicts import merge_two_dicts
import SimulationFramework.Modules.read_beam_file as rbf
import SimulationFramework.Codes.Executables as exes
from SimulationFramework.Codes.ASTRA.ASTRA import *
from SimulationFramework.Codes.CSRTrack.CSRTrack import *
from SimulationFramework.Codes.Elegant.Elegant import *
from SimulationFramework.Codes.Generators.Generators import *
from SimulationFramework.Codes.GPT.GPT import *
from SimulationFramework.Codes.MAD8.MAD8 import *
import progressbar
from munch import Munch, unmunchify

_mapping_tag = yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG

def dict_representer(dumper, data):
    return dumper.represent_dict(iter(list(data.items())))

def dict_constructor(loader, node):
    return OrderedDict(loader.construct_pairs(node))

yaml.add_representer(OrderedDict, dict_representer)
yaml.add_constructor(_mapping_tag, dict_constructor)

class Framework(Munch):

    def __init__(self, directory='test', master_lattice=None, overwrite=None, runname='CLARA_240', clean=False, verbose=True, sddsindex=0):
        super(Framework, self).__init__()
        # global master_lattice_location
        self.global_parameters = {'beam': rbf.beam(sddsindex=sddsindex)}
        self.verbose = verbose
        self.subdir = directory
        self.clean = clean
        self.elementObjects = OrderedDict()
        self.latticeObjects = OrderedDict()
        self.commandObjects = OrderedDict()
        self.groupObjects = OrderedDict()
        self.progress = 0
        self.basedirectory = os.getcwd()
        self.filedirectory = os.path.dirname(os.path.abspath(__file__))
        self.overwrite = overwrite
        self.runname = runname
        if self.subdir is not None:
            self.setSubDirectory(self.subdir)
        self.setMasterLatticeLocation(master_lattice)

        self.executables = exes.Executables(self.global_parameters['master_lattice_location'])
        self.defineASTRACommand = self.executables.define_astra_command
        self.defineElegantCommand = self.executables.define_elegant_command
        self.defineASTRAGeneratorCommand = self.executables.define_ASTRAgenerator_command
        self.defineCSRTrackCommand = self.executables.define_csrtrack_command
        self.define_gpt_command = self.executables.define_gpt_command

    def setSubDirectory(self, dir):
        # global self.global_parameters['master_subdir'], self.global_parameters['master_lattice_location']
        self.subdirectory = os.path.abspath(self.basedirectory+'/'+dir)
        self.global_parameters['master_subdir'] = self.subdirectory
        if not os.path.exists(self.subdirectory):
            os.makedirs(self.subdirectory, exist_ok=True)
        else:
            if self.clean == True:
                clean_directory(self.subdirectory)
        if self.overwrite == None:
            self.overwrite = True
        self.setMasterLatticeLocation()

    def setMasterLatticeLocation(self, master_lattice=None):
        # global master_lattice_location
        if master_lattice is None:
            self.global_parameters['master_lattice_location'] = (os.path.abspath(os.path.dirname(os.path.abspath(__file__)) + '/../MasterLattice/')+'/').replace('\\','/')
        else:
            self.global_parameters['master_lattice_location'] = master_lattice
        self.master_lattice_location = self.global_parameters['master_lattice_location']

    def load_Elements_File(self, input):
        if isinstance(input,(list,tuple)):
            filename = input
        else:
            filename = [input]
        for f in filename:
            if os.path.isfile(f):
                with open(f, 'r') as stream:
                    elements = yaml.load(stream, Loader=yaml.UnsafeLoader)['elements']
            else:
                with open(self.global_parameters['master_lattice_location'] + f, 'r') as stream:
                    elements = yaml.load(stream, Loader=yaml.UnsafeLoader)['elements']
            for name, elem in list(elements.items()):
                self.read_Element(name, elem)

    def loadSettings(self, filename='short_240.settings'):
        """Load Lattice Settings from file"""
        self.settingsFilename = filename
        # print 'self.settingsFilename = ', self.settingsFilename
        if os.path.exists(filename):
            stream = open(filename, 'r')
        else:
            stream = open(self.global_parameters['master_lattice_location']+filename, 'r')
        self.settings = yaml.load(stream, Loader=yaml.UnsafeLoader)
        self.globalSettings = self.settings['global']
        master_run_no = self.globalSettings['run_no'] if 'run_no' in self.globalSettings else 1
        if 'generator' in self.settings:
            self.generatorSettings = self.settings['generator']
            self.add_Generator(**self.generatorSettings)
        self.fileSettings = self.settings['files']
        elements = self.settings['elements']
        self.groups = self.settings['groups'] if 'groups' in self.settings and self.settings['groups'] is not None else {}
        changes = self.settings['changes'] if 'changes' in self.settings and self.settings['changes'] is not None else {}
        stream.close()

        for name, elem in list(self.groups.items()):
            group = globals()[elem['type']](name, self.elementObjects, global_parameters=self.global_parameters, **elem)
            self.groupObjects[name] = group

        for name, elem in list(elements.items()):
            self.read_Element(name, elem)

        for name, lattice in list(self.fileSettings.items()):
            self.read_Lattice(name, lattice)

        self.apply_changes(changes)

        self.original_elementObjects = {}
        for e in self.elementObjects:
            self.original_elementObjects[e] = unmunchify(self.elementObjects[e])


    def read_Lattice(self, name, lattice):
        code = lattice['code'] if 'code' in lattice else 'astra'
        self.latticeObjects[name] = globals()[lattice['code'].lower()+'Lattice'](name, lattice, self.elementObjects, self.groupObjects, self.settings, self.executables, self.global_parameters)

    def convert_numpy_types(self, v):
        if isinstance(v, (np.ndarray, list, tuple)):
            return [self.convert_numpy_types(l) for l in v]
        elif isinstance(v, (np.float64, np.float32, np.float16, np.float_ )):
            return float(v)
        elif isinstance(v, (np.int_, np.intc, np.intp, np.int8, np.int16, np.int32, np.int64, np.uint8, np.uint16, np.uint32, np.uint64)):
            return int(v)
        else:
            return v

    def detect_changes(self, elementtype=None, elements=None, function=None):
        start = time.time()
        changedict = {}
        if elementtype is not None:
            changeelements = self.getElementType(elementtype, 'objectname')
        elif elements is not None:
            changeelements = elements
        else:
            changeelements = list(self.elementObjects.keys())
        # print('changeelements = ', changeelements)
        # print(changeelements[0])
        if len(changeelements) > 0 and isinstance(changeelements[0], (list, tuple, dict)) and len(changeelements[0]) > 1:
                for ek in changeelements:
                    new = None
                    e, k = ek[:2]
                    if e in self.elementObjects:
                        new = unmunchify(self.elementObjects[e])
                    elif e in self.groupObjects:
                        new = self.groupObjects[e]
                        # print 'detecting group = ', e, new, new[k]
                    if new is not None:
                        # print (new)
                        if e not in changedict:
                            changedict[e] = {}
                        changedict[e][k] = self.convert_numpy_types(new[k])
        else:
            for e in changeelements:
                # print 'saving element: ', e
                if not self.original_elementObjects[e] == unmunchify(self.elementObjects[e]):
                    orig = self.original_elementObjects[e]
                    new = unmunchify(self.elementObjects[e])
                    try:
                        changedict[e] = {k: self.convert_numpy_types(new[k]) for k in new if k in orig and not new[k] == orig[k]}
                        changedict[e].update({k: self.convert_numpy_types(new[k]) for k in new if k not in orig})
                    except:
                        print ('##### ERROR IN CHANGE ELEMS: ', e, new)
                        pass
        return changedict

    def save_changes_file(self, filename=None, type=None, elements=None, function=None):
        if filename is None:
            pre, ext = os.path.splitext(os.path.basename(self.settingsFilename))
            filename =  pre     + '_changes.yaml'
        changedict = self.detect_changes(elementtype=type, elements=elements, function=function)
        with open(filename,"w") as yaml_file:
            yaml.dump(changedict, yaml_file, default_flow_style=False)

    def save_lattice(self, lattice=None, filename=None, directory='.'):
        if filename is None:
            pre, ext = os.path.splitext(os.path.basename(self.settingsFilename))
        dict = OrderedDict({'elements': OrderedDict()})
        latticedict = dict['elements']
        if lattice is None:
            elements = list(self.elementObjects.keys())
            filename =  pre     + '_lattice.yaml'
        else:
            if self.latticeObjects[lattice].elements is None:
                return
            elements = list(self.latticeObjects[lattice].elements.keys())
            filename =  pre + '_' + lattice + '_lattice.yaml'
        disallowed = ['allowedkeywords', 'keyword_conversion_rules_elegant', 'objectdefaults','global_parameters']
        for e in elements:
            new = unmunchify(self.elementObjects[e])
            if ('subelement' in new and not new['subelement']) or not 'subelement' in new:
                try:
                    latticedict[e] = {k.replace('object',''): self.convert_numpy_types(new[k]) for k in new if not k in disallowed}
                    # latticedict[e].update({k: self.convert_numpy_types(new[k]) for k in new})
                except:
                    print ('##### ERROR IN CHANGE ELEMS: ', e, new)
                    pass
        print('#### Saving Lattice - ', filename)
        with open(directory + '/' + filename,"w") as yaml_file:
            yaml.dump(dict, yaml_file, default_flow_style=False)

    def load_changes_file(self, filename=None, apply=True):
        if isinstance(filename, (tuple, list)):
            for c in filename:
                self.load_changes_file(c)
        else:
            if filename is None:
                pre, ext = os.path.splitext(os.path.basename(self.settingsFilename))
                filename =  pre     + '_changes.yaml'
            with open(filename,"r") as infile:
                changes = dict(yaml.load(infile, Loader=yaml.UnsafeLoader))
            if apply:
                self.apply_changes(changes)
            else:
                return changes

    def apply_changes(self, changes):
        for e, d in list(changes.items()):
            # print 'found change element = ', e
            if e in self.elementObjects:
                # print 'change element exists!'
                for k, v in list(d.items()):
                    self.modifyElement(e, k, v)
                    # print ('modifying ',e,'[',k,']', ' = ', v)
            if e in self.groupObjects:
                # print ('change group exists!')
                for k, v in list(d.items()):
                    self.groupObjects[e].change_Parameter(k, v)
                    # print ('modifying ',e,'[',k,']', ' = ', v)

    def check_lattice(self, decimals=4):
        for elem in self.elementObjects.values():
            start = elem.position_start
            end = elem.position_end
            length = elem.length
            theta = elem.global_rotation[2]
            if elem.objecttype == 'dipole':
                angle = elem.angle
                rho = length / angle
                clength = np.array([rho * (np.cos(angle) - 1), 0, rho * np.sin(angle)])
            else:
                clength = np.array([0, 0, length])
            cend = start + np.dot(clength, _rotation_matrix(theta))
            if not np.round(cend - end, decimals=decimals).any() == 0:
                print (elem.objectname, cend - end)

    def change_Lattice_Code(self, name, code, exclude=None):
        if name == 'All':
            [self.change_Lattice_Code(l, code, exclude) for l in self.latticeObjects]
        elif isinstance(name, (tuple, list)):
            [self.change_Lattice_Code(l, code, exclude) for l in name]
        else:
            if not name == 'generator' and not (name == exclude or (isinstance(exclude, (list, tuple)) and name in exclude)):
                # print('Changing lattice ', name, ' to ', code.lower())
                currentLattice = self.latticeObjects[name]
                self.latticeObjects[name] = globals()[code.lower()+'Lattice'](currentLattice.objectname, currentLattice.file_block, self.elementObjects, self.groupObjects, self.settings, self.executables, self.global_parameters)

    def read_Element(self, name, element, subelement=False):
        if name == 'filename':
            self.load_Elements_File(element)
        else:
            if subelement:
                self.add_Element(name, subelement=True, **element)
            else:
                self.add_Element(name, **element)
            if 'sub_elements' in element:
                for name, elem in list(element['sub_elements'].items()):
                    self.read_Element(name, elem, subelement=True)

    def add_Element(self, name=None, type=None, **kwargs):
        if name == None:
            if not 'name' in kwargs:
                raise NameError('Element does not have a name')
            else:
                name = kwargs['name']
        # try:
        element = globals()[type](name, type, global_parameters=self.global_parameters, **kwargs)
        # print element
        self.elementObjects[name] = element
        return element
        # except Exception as e:
        #     raise NameError('Element \'%s\' does not exist' % type)

    def getElement(self, element, param=None):
        if self.__getitem__(element) is not None:
            if param is not None:
                param = param.lower()
                return getattr(self.__getitem__(element), param)
            else:
                return self.__getitem__(element)
        else:
            print(( 'WARNING: Element ', element,' does not exist'))
            return {}

    def getElementType(self, type, setting=None):
        if isinstance(type, (list, tuple)):
            all_elements = [self.getElementType(t) for t in type]
            return [item for sublist in all_elements for item in sublist]
        return [self.elementObjects[element] if setting is None else self.elementObjects[element][setting] for element in list(self.elementObjects.keys()) if self.elementObjects[element].objecttype.lower() == type.lower()]

    def setElementType(self, type, setting, values):
        elems = self.getElementType(type)
        if len(elems) == len(values):
            for e, v  in zip(elems, values):
                e[setting] = v
        else:
            # print(( len(elems), len(values)))
            raise ValueError

    def modifyElement(self, elementName, parameter, value):
        if elementName in self.groupObjects:
            self.groupObjects[elementName].change_Parameter(parameter,value)
        elif elementName in self.elementObjects:
            setattr(self.elementObjects[elementName], parameter, value)
        # set(getattr(self.elementObjects[elementName], parameter), value)

    def add_Generator(self, default=None, **kwargs):
        if default in astra_generator_keywords['defaults']:
            self.generator = ASTRAGenerator(self.executables, self.global_parameters, **merge_two_dicts(kwargs,astra_generator_keywords['defaults'][default]))
        else:
            self.generator = ASTRAGenerator(self.executables, self.global_parameters, **kwargs)
        self.latticeObjects['generator'] = self.generator

    def loadParametersFile(self, file):
        pass

    def saveParametersFile(self, file, parameters):
        output = {}
        if isinstance(parameters, dict):
            for k,v in list(parameters.items()):
                output[k] = {}
                if isinstance(v, (list, tuple)):
                    for p in v:
                        output[k][p] = getattr(self[k],p)
                else:
                    output[k][v] = getattr(self[k],v)
        elif isinstance(parameters, (list,tuple)):
            for k, v in parameters:
                output[k] = {}
                if isinstance(v, (list, tuple)):
                    for p in v:
                        output[k][p] = getattr(self[k],p)
                else:
                    output[k][v] = getattr(self[k],v)
        with open(file,"w") as yaml_file:
            yaml.dump(output, yaml_file, default_flow_style=False)
        # try:
        # elem = self.getelement(k, v)
        # outputfile.write(k+' '+v+' ')

    def __getitem__(self,key):
        if key in super(Framework, self).__getitem__('elementObjects'):
            return self.elementObjects.get(key)
        elif key in super(Framework, self).__getitem__('latticeObjects'):
            return self.latticeObjects.get(key)
        elif key in super(Framework, self).__getitem__('groupObjects'):
            return self.groupObjects.get(key)
        else:
            try:
                return super(Framework, self).__getitem__(key)
            except:
                return None

    @property
    def elements(self):
        return list(self.elementObjects.keys())

    @property
    def lines(self):
        return list(self.latticeObjects.keys())

    @property
    def commands(self):
        return list(self.commandObjects.keys())

    def track(self, files=None, startfile=None, endfile=None, preprocess=True, write=True, track=True, postprocess=True):
        self.progress = 0
        if files is None:
            files = ['generator'] + self.lines if hasattr(self, 'generator') else self.lines
        if startfile is not None and startfile in files:
            index = files.index(startfile)
            files = files[index:]
        if endfile is not None and endfile in files:
            index = files.index(endfile)
            files = files[:index+1]
        if self.verbose:
            format_custom_text = progressbar.FormatCustomText(
                'File: %(running)s', {'running': ''}
            )
            bar = progressbar.ProgressBar(widgets=[format_custom_text, progressbar.Percentage(), progressbar.Bar(), progressbar.Percentage(),], max_value=len(files))
            format_custom_text.update_mapping(running=files[0]+'  ')
            for i in bar(list(range(len(files)))):
                l = files[i]
                self.progress = 100. * (i+1)/len(files)
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
                self.progress = 100. * (i)/len(files)
                if l == 'generator' and hasattr(self, 'generator'):
                    if write:
                        self.generator.write()
                    self.progress = 100. * (i+0.33)/len(files)
                    if track:
                        self.generator.run()
                    self.progress = 100. * (i+0.66)/len(files)
                    if postprocess:
                        self.generator.astra_to_hdf5()
                else:
                    if preprocess:
                        self.latticeObjects[l].preProcess()
                    self.progress = 100. * (i+0.25)/len(files)
                    if write:
                        self.latticeObjects[l].write()
                    self.progress = 100. * (i+0.5)/len(files)
                    if track:
                        self.latticeObjects[l].run()
                    self.progress = 100. * (i+0.75)/len(files)
                    if postprocess:
                        self.latticeObjects[l].postProcess()
            self.progress = 100
