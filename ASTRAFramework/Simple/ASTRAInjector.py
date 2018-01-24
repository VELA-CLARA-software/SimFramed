import yaml, collections, subprocess, os, re, time
import numpy
import h5py
import glob
from copy import deepcopy

_mapping_tag = yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG

def dict_representer(dumper, data):
    return dumper.represent_dict(data.iteritems())

def dict_constructor(loader, node):
    return collections.OrderedDict(loader.construct_pairs(node))

yaml.add_representer(collections.OrderedDict, dict_representer)
yaml.add_constructor(_mapping_tag, dict_constructor)

class ASTRAInjector(object):

    def __init__(self, subdir='test', overwrite=None, runname='CLARA_240'):
        super(ASTRAInjector, self).__init__()
        self.lineIterator = 0
        self.astraCommand = ['astra']
        self.basedirectory = os.getcwd()
        self.subdir = subdir
        self.overwrite = overwrite
        self.runname = runname
        self.subdirectory = self.basedirectory+'/'+subdir
        if not os.path.exists(self.subdirectory):
            os.makedirs(self.subdirectory)
        if self.overwrite == None:
            if not os.path.exists(self.subdirectory):
                os.makedirs(self.subdirectory)
                self.overwrite = True
            else:
                response = raw_input('Overwrite existing directory ? [Y/n]')
                self.overwrite = True if response in {'','y','Y','yes','Yes','YES'} else False
        self.astraFiles = []

    def loadSettings(self, filename='short_240.settings'):
        self.settings = {}
        if isinstance(filename, dict):
            settings = filename
        else:
            self.settingsFile = filename
            stream = file(filename, 'r')
            settings = yaml.load(stream)
            stream.close()
        self.settings = settings
        self.globalSettings = settings['global']
        self.fileSettings = settings['files']

    def saveSettings(self, filename='example.settings'):
        self.settings = {}
        if isinstance(filename, dict):
            settings = filename
        else:
            self.settingsFile = filename
            stream = file(filename, 'r')
            settings = yaml.load(stream)
            stream.close()
        self.globalSettings = settings['global']
        self.fileSettings = settings['files']

    def findKey(self, key, var):
        if hasattr(var,'iteritems'):
            for k, v in var.iteritems():
                if k == key:
                    yield k, v
                if isinstance(v, dict):
                    for result in self.findKey(key, v):
                        yield result
                elif isinstance(v, list):
                    for d in v:
                        for result in self.findKey(key, d):
                            yield result

    def replaceKeyValue(self, container, key, value):
        copy = container
        try:
            if key in container:
                container[key] = value
                return container
        except:
            return container
        else:
            if isinstance(container, list):
                for i, v in enumerate(container):
                    container[i] = self.replaceKeyValue(v, key, value)
            elif isinstance(container, tuple):
                container = list(container)
                for i, v in enumerate(container):
                    container[i] = self.replaceKeyValue(v, key, value)
                container = tuple(container)
            elif isinstance(container, dict):
                for k in container.keys():
                    container[k] = self.replaceKeyValue(container[k], key, value)
        return container

    def modifySetting(self, key, value):
        self.replaceKeyValue(self.fileSettings, key, value)
        self.replaceKeyValue(self.globalSettings, key, value)

    def readFile(self, fname=None):
        with open(fname) as f:
            content = f.readlines()
        return content

    def lineReplaceFunction(self, line, findString, replaceString, i=None):
        if findString in line:
            self.lineIterator += 1
            if not i == None:
                return line.replace('$'+findString+'$', str(replaceString[i]))
            else:
                return line.replace('$'+findString+'$', str(replaceString))
        else:
            return line

    def replaceString(self, lines=[], findString=None, replaceString=None):
        if isinstance(replaceString, (list, tuple, set)):
            self.lineIterator = 0
            return [self.lineReplaceFunction(line, findString, replaceString, self.lineIterator) for line in lines]
        else:
            return [self.lineReplaceFunction(line, findString, replaceString) for line in lines]

    def saveFile(self, lines=[], filename='runastratemp.in'):
        # print 'filename = ', self.subdir+'/'+filename
        stream = file(self.subdir+'/'+filename, 'w')
        for line in lines:
            stream.write(line)
        stream.close()

    def applySettings(self, filelines=None):
        if self.overwrite:
            for f, s in self.fileSettings.iteritems():
                newfilename = f
                # print 'newfilename = ', newfilename
                if not newfilename[-3:-1] == '.in':
                    newfilename = newfilename+'.in'
                if filelines is None:
                    lines = self.readFile(f)
                else:
                    lines = filelines[f]
                for var, val in s.iteritems():
                    lines = self.replaceString(lines, var, val)
                for var, val in self.globalSettings.iteritems():
                    lines = self.replaceString(lines, var, val)
                self.saveFile(lines, newfilename)
                self.astraFiles.append(newfilename)

    def runASTRAFiles(self, files=[]):
        if self.overwrite:
            if files == []:
                for f in self.astraFiles:
                    self.runASTRA(f)
            else:
                for f in files:
                    self.runASTRA(f)

    def runASTRA(self, filename=''):
        command = self.astraCommand + [filename]
        with open(os.devnull, "w") as f:
            subprocess.call(command, stdout=f, cwd=self.subdir)

    def defineASTRACommand(self,command=['astra']):
        self.astraCommand = command

    def setInitialDistribution(self, filename='../1k-250pC-76fsrms-1mm_TE09fixN12.ini'):
        self.globalSettings['initial_distribution'] = filename

    def createInitialDistribution(self, npart=1000, charge=250, generatorCommand=None):
        self.globalSettings['npart'] = npart
        self.globalSettings['charge'] = charge/1000.0
        if self.overwrite:
            astragen = ASTRAGenerator(self.subdir, charge, npart)
            if not generatorCommand is None:
                astragen.defineGeneratorCommand(generatorCommand)
            elif os.name == 'nt':
                astragen.defineGeneratorCommand(['generator_7June2007'])
            else:
                astragen.defineGeneratorCommand(['/opt/ASTRA/generator.sh'])
            inputfile = astragen.generateBeam()
            self.setInitialDistribution(inputfile)
            scgrid = getGrids(npart)
            self.globalSettings['SC_2D_Nrad'] = max([scgrid.gridSizes,4])
            self.globalSettings['SC_2D_Nlong'] = max([scgrid.gridSizes,4])
            for scvar in ['SC_3D_Nxf','SC_3D_Nyf','SC_3D_Nzf']:
                self.globalSettings[scvar] = scgrid.gridSizes

    def getScreenFiles(self):
        self.screenpositions = {}
        for f, s in self.fileSettings.iteritems():
            filenumber = re.search('\d\d\d',f).group(0)
            files = glob.glob(self.subdir+'/'+'test.in.'+filenumber+'.????.'+filenumber)
            screenpositions = [re.search('\d\d\d\d',s).group(0) for s in files]
            self.screenpositions[filenumber] = screenpositions
        return self.screenpositions

    def createHDF5Summary(self, screens=None):
        savescreens = [str(s) for s in screens]
        screenpositions = self.getScreenFiles()
        filename = '_'.join(map(str,[self.runname,self.globalSettings['charge'],self.globalSettings['npart']])) + '.hdf5'
        # print filename
        f = h5py.File(filename, "w")
        inputgrp = f.create_group("Input")
        inputgrp['charge'] = self.globalSettings['charge']
        inputgrp['npart'] = self.globalSettings['npart']
        inputgrp['subdirectory'] = self.subdirectory
        xemitgrp = f.create_group("Xemit")
        yemitgrp = f.create_group("Yemit")
        zemitgrp = f.create_group("Zemit")
        screengrp = f.create_group("screens")
        inputgrp.create_dataset('initial_distribution',data=numpy.loadtxt(self.subdir+'/'+self.globalSettings['initial_distribution']))
        for n, screens in screenpositions.iteritems():
            inputfile = file(self.subdir+'/'+'test.in.'+n+'.in','r')
            inputfilecontents = inputfile.read()
            inputgrp.create_dataset(n,data=inputfilecontents)
            for emit, grp in {'X': xemitgrp,'Y': yemitgrp,'Z': zemitgrp}.iteritems():
                emitfile = self.subdir+'/'+'test.in.'+n+'.'+emit+'emit.'+n
                grp.create_dataset(n,data=numpy.loadtxt(emitfile))
            for s in screens:
                screenfile = self.subdir+'/'+'test.in.'+n+'.'+s+'.'+n
                if screens == None or s in savescreens:
                    # print 'screen = ', s
                    screengrp.create_dataset(s,data=numpy.loadtxt(screenfile))

class feltools(object):

    def __init__(self, subdir='.'):
        super(feltools, self).__init__()
        self.subdirectory = subdir
        self.basedirectory = os.getcwd()
        self.astra2elegantCommand = ['astra2elegant']
        self.sddsprocessCommand = ['sddsprocess']
        self.sddsMatchTwissCommand = ['sddsmatchtwiss']

    def convertToSDDS(self, inputfile='test.in.128.4929.128', outputfile=None):
        outputfile = outputfile if not outputfile == None else inputfile+'.sdds'
        command = self.astra2elegantCommand + [inputfile, outputfile]
        comm = subprocess.call(command, cwd=self.subdirectory)
        return outputfile

    def filterInput(self, inputstr, parameters, **kwargs):
        # print kwargs
        # print parameters
        kwargs = {k.lower():v for k,v in kwargs.items()}
        if any(p.lower() in kwargs for p in parameters):
            string = inputstr+'='
            stringsep = ''
            for p, s in parameters.iteritems():
                string = string+stringsep+s+'='+str(kwargs[p.lower()]) if p.lower() in kwargs else string
                if p in kwargs:
                    stringsep = ','
            return string
        else:
            return ''

    def sddsMatchTwiss(self, bunch, outputbunch, **kwargs):
        xstr = self.filterInput('-xPlane',{'betax': 'beta','alphax':'alpha','etax':'etaValue','etaxslope':'etaSlope'},**kwargs)
        ystr = self.filterInput('-yPlane',{'betay': 'beta','alphay':'alpha','etay':'etaValue','etayslope':'etaSlope'},**kwargs)
        zstr = self.filterInput('-zPlane',{'deltaStDev': 'deltaStDev','tStDev':'tStDev','correlation':'correlation','chirp':'chirp','betaGamma':'betaGamma'},**kwargs)
        outputcommand = self.sddsMatchTwissCommand + [bunch, outputbunch]
        outputcommand = outputcommand + [xstr] if not xstr == '' else outputcommand
        outputcommand = outputcommand + [ystr] if not ystr == '' else outputcommand
        outputcommand = outputcommand + [zstr] if not zstr == '' else outputcommand
        outputcommand += ['-saveMatrices='+kwargs['saveMatrices']] if 'saveMatrices' in kwargs else []
        outputcommand += ['-loadMatrices='+kwargs['loadMatrices']] if 'loadMatrices' in kwargs else []
        outputcommand += ['-nowarnings ='+kwargs['nowarnings']] if 'nowarnings' in kwargs else []
        comm = subprocess.call(outputcommand, cwd=self.subdirectory)

    def compressOutputBunch(self, bunch='test.in.128.4929.128.sdds', dt=5e-13, outputbunch='compressed.sdds'):
        self.sddsMatchTwiss(bunch,outputbunch,tStDev=dt)

class ASTRAGenerator(object):

    def __init__(self, subdir='test', charge=250, npart=1000):
        super(ASTRAGenerator, self).__init__()
        self.lineIterator = 0
        self.generatorBaseFile = 'generator.in'
        self.charge = charge
        self.npart = npart
        self.generatorCommand = ['generator']
        self.basedirectory = os.getcwd()
        self.subdir = subdir
        self.subdirectory = self.basedirectory+'/'+self.subdir
        if not os.path.exists(self.subdirectory):
            os.makedirs(self.subdirectory)

    def generateBeam(self):
        self.createSettings()
        self.createGeneratorInput()
        return self.settings['filename']

    def readFile(self, fname=None):
        with open(fname) as f:
            content = f.readlines()
        return content

    def lineReplaceFunction(self, line, findString, replaceString):
        if findString in line:
            return line.replace('$'+findString+'$', str(replaceString))
        else:
            return line

    def replaceString(self, lines=[], findString=None, replaceString=None):
        return [self.lineReplaceFunction(line, findString, replaceString) for line in lines]

    def saveFile(self, lines=[], filename='generatortemp.in'):
        stream = file(filename, 'w')
        for line in lines:
            stream.write(line)
        stream.close()

    def particleSuffix(self):
        suffix = str(int(round(self.npart/1e9))) + 'G'
        if self.npart < 1e9:
            suffix = str(int(round(self.npart/1e6))) + 'M'
        if self.npart < 1e6:
            suffix = str(int(round(self.npart/1e3))) + 'k'
        if self.npart < 1e3:
            suffix = str(int(round(self.npart)))
        return suffix

    def createSettings(self):
        self.settings = {}
        self.settings['charge'] = self.charge/1000.0
        self.settings['number_particles'] = self.npart
        self.settings['filename'] = self.particleSuffix() + '-' + str(self.charge) + 'pC-76fsrms-1mm_TE09fixN12.ini'

    def createGeneratorInput(self):
        lines = self.readFile(self.generatorBaseFile)
        for var, val in self.settings.iteritems():
            lines = self.replaceString(lines, var, val)
        self.saveFile(lines, self.subdir+'/'+'generator.in')
        self.runGenerator('generator.in')

    def runGenerator(self, filename=''):
        command = self.generatorCommand + [filename]
        with open(os.devnull, "w") as f:
            subprocess.call(command, stdout=f, cwd=self.subdir)

    def defineGeneratorCommand(self,command=['generator']):
        self.generatorCommand = command

class getGrids(object):

    def __init__(self, npart=1000):
        self.powersof8 = numpy.asarray([ 2**(j) for j in range(1,20) ])
        self.n = npart
        self.gridSizes = self.getGridSizes(int(self.n))

    def getGridSizes(self):
        return self.gridSizes

    def getGridSizes(self, x):
        self.x = abs(x)
        self.cuberoot = int(round(self.x ** (1. / 3)))
        return max([4,self.find_nearest(self.powersof8, self.cuberoot)])

    def find_nearest(self, array, value):
        self.array = array
        self.value = value
        self.idx = (numpy.abs(self.array - self.value)).argmin()
        return self.array[self.idx]
