import yaml, collections, subprocess, os, math, re, sys
from shutil import copyfile
import numpy as np
from operator import add
from FrameworkHelperFunctions import *
import ASTRAGenerator as GenPart
from ASTRARules import *
from getGrids import *
sys.path.append('..')
import read_beam_file as rbf
from collections import defaultdict

_mapping_tag = yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG

def dict_representer(dumper, data):
    return dumper.represent_dict(data.iteritems())

def dict_constructor(loader, node):
    return collections.OrderedDict(loader.construct_pairs(node))

yaml.add_representer(collections.OrderedDict, dict_representer)
yaml.add_constructor(_mapping_tag, dict_constructor)

def merge_two_dicts(x, y):
    z = x.copy()   # start with x's keys and values
    z.update(y)    # modifies z with y's keys and values & returns None
    return dict(z)

class Framework(object):

    def __init__(self, subdir='test', overwrite=None, runname='CLARA_240'):
        super(Framework, self).__init__()
        self.lineIterator = 0
        self.basedirectory = os.getcwd()
        self.subdir = subdir
        self.overwrite = overwrite
        self.runname = runname
        self.subdirectory = self.basedirectory+'/'+subdir
        self.globalSettings = dict()
        self.fileSettings = dict()
        self.elements = dict()
        self.groups = dict()
        self.astra = ASTRA(parent=self, directory=self.subdir)
        self.CSRTrack = CSRTrack(parent=self, directory=self.subdir)
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

    def loadElementsFile(self, input):
        if isinstance(input,(list,tuple)):
            filename = input
        else:
            filename = [input]
        for f in filename:
            stream = file(f, 'r')
            elements = yaml.load(stream)['elements']
            stream.close()
            if 'filename' in elements:
                self.loadElementsFile(elements['filename'])
            self.elements = merge_two_dicts(self.elements, elements)

    def loadSettings(self, filename='short_240.settings'):
        """Load Lattice Settings from file"""
        stream = file(filename, 'r')
        settings = yaml.load(stream)
        self.globalSettings = settings['global']
        self.generatorFile = self.globalSettings['generatorFile'] if 'generatorFile' in self.globalSettings else None
        self.fileSettings = settings['files']
        self.elements = settings['elements']
        self.groups = settings['groups']
        stream.close()
        if 'filename' in self.elements:
            self.loadElementsFile(self.elements['filename'])

    def getFileSettings(self, file, block):
        """Return the correct settings 'block' from 'file' dict if exists else return empty dict"""
        if file in self.fileSettings and block in self.fileSettings[file]:
            return self.fileSettings[file][block]
        else:
            return {}

    def getSettingsBlock(self, dict, block):
        """Return the correct settings 'block' from dict if exists else return empty dict"""
        if block in dict:
            return dict[block]
        else:
            return {}

    def modifyFile(self, filename, setting, value):
        if filename in self.fileSettings:
            if isinstance(setting, (list,tuple)):
                dic = self.fileSettings[filename]
                for key in setting[:-1]:
                    dic = dic.setdefault(key, {})
                dic[setting[-1]] = value
            elif not setting in self.fileSettings:
                self.fileSettings[setting] = {}
                self.fileSettings[filename][setting] = value

    def getElement(self, element='', setting=None):
        """return 'element' from the main elements dict"""
        if element in self.elements:
            if setting is not None:
                return self.elements[element][setting]
            else:
                return self.elements[element]
        else:
            return []

    def modifyElement(self, element='', setting='', value=''):
        """return 'element' from the main elements dict"""
        element = self.getElement(element)
        if setting in element:
            element[setting] = value

    def getElementType(self, type='', setting=None):
        """return 'element' from the main elements dict"""
        elems = []
        for name, element in self.elements.viewitems():
            if 'type' in element and element['type'] == type:
                elems.append(name)
                # print element
        elems = sorted(elems, key=lambda x: self.elements[x]['position_start'][2])
        if setting is not None:
            return [self.elements[x][setting] for x in elems]
        else:
            return elems

    def setElementType(self, type='', setting=None, values=[]):
        """return 'element' from the main elements dict"""
        elems = self.getElementType(type)
        if len(elems) == len(values):
            for e, v  in zip(elems, values):
                self.elements[e][setting] = v
        else:
            raise ValueError

    def getElementsBetweenS(self, elementtype, output={}, zstart=None, zstop=None):
        # zstart = zstart if zstart is not None else getParameter(output,'zstart',default=0)
        if zstart is None:
            zstart = getParameter(output,'zstart',default=None)
            if zstart is None:
                startelem = getParameter(output,'start_element',default=None)
                if startelem is None or startelem not in self.elements:
                    zstart = 0
                else:
                    zstart = self.elements[startelem]['position_start'][2]
        # zstop = zstop if zstop is not None else getParameter(output,'zstop',default=0)
        if zstop is None:
            zstop = getParameter(output,'zstop',default=None)
            if zstop is None:
                endelem = getParameter(output,'end_element',default=None)
                if endelem is None or endelem not in self.elements:
                    zstop = 0
                else:
                    zstop = self.elements[endelem]['position_end'][2]

        elements = findSetting('type',elementtype,dictionary=self.elements)
        elements = sorted([[s[1]['position_start'][2],s[0]] for s in elements if s[1]['position_start'][2] >= zstart and s[1]['position_start'][2] <= zstop])

        return [e[1] for e in elements]

    def getGroup(self, group=''):
        """return all elements in a group from the main elements dict"""
        elements = []
        if group in self.groups:
            groupelements = self.groups[group]['elements']
            for e in groupelements:
                elements.append([e,self.getElement(e)])
        return elements

    def xform(self, theta, tilt, length, x, r):
        """Calculate the change on local coordinates through an element"""
        theta = theta if abs(theta) > 1e-6 else 1e-6
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

    def elementPositions(self, elements, startpos=None):
        """Calculate element positions for the given 'elements'"""
        anglesum = [0]
        localXYZ = np.identity(3)
        if startpos == None:
            startpos = elements[0][1]['position_start']
            if len(startpos) == 1:
                startpos = [0,0,startpos]
        x1 = np.matrix(np.transpose([startpos]))
        x = [np.array(x1)]
        for name, d in elements:
            angle = getParameter(d,'angle',default=1e-9)
            anglesum.append(anglesum[-1]+angle)
            x1, localXYZ = self.xform(angle, 0, getParameter(d,'length'), x1, localXYZ)
            x.append(x1)
        return zip(x, anglesum[:-1], elements), localXYZ

    def createDrifts(self, elements, startpos=None):
        """Insert drifts into a sequence of 'elements'"""
        positions = []
        elementno = 0
        for name, e in elements:
            pos = np.array(e['position_start'])
            positions.append(pos)
            length = np.array([0,0,e['length']])
            positions.append(pos+length)
        if not startpos == None:
            positions.prepend(startpos)
        else:
            positions = positions[1:]
            positions.append(positions[-1]+[0,0,0.1])
        driftdata = list(chunks(positions, 2))
        for d in driftdata:
            if len(d) > 1:
                elementno += 1
                length = d[1][2] - d[0][2]
                elements.append(['drift'+str(elementno),
                                {'length': length, 'type': 'drift',
                                 'position_start': list(d[0])
                                }])
        return sorted(elements, key=sortByPositionFunction)

    def setDipoleAngle(self, dipole, angle=0):
        """Set the dipole angle for a given 'dipole'"""
        name, d = dipole
        if getParameter(d,'entrance_edge_angle') == 'angle':
            d['entrance_edge_angle'] = np.sign(d['angle'])*angle
        if getParameter(d,'exit_edge_angle') == 'angle':
            d['exit_edge_angle'] = np.sign(d['angle'])*angle
        d['angle'] = np.sign(d['angle'])*angle
        return [name, d]

    def createInputFiles(self):
        for f in self.fileSettings.keys():
            filename = self.subdirectory+'/'+f+'.in'
            if 'code' in self.fileSettings[f]:
                code = self.fileSettings[f]['code']
                if code.upper() == 'ASTRA':
                    saveFile(filename, lines=self.astra.createASTRAFileText(f))
                if code.upper() == 'CSRTRACK':
                    saveFile(filename, lines=self.CSRTrack.createCSRTrackFileText(f))
            else:
                saveFile(filename, lines=self.astra.createASTRAFileText(f))

    def runInputFiles(self, files=None):
        if not isinstance(files, (list, tuple)):
            files = self.fileSettings.keys()
        for f in files:
            if f in self.fileSettings.keys():
                filename = f+'.in'
                # print 'Running file: ', f
                if 'code' in self.fileSettings[f]:
                    code = self.fileSettings[f]['code']
                    if code.upper() == 'ASTRA':
                        self.astra.runASTRA(filename)
                    if code.upper() == 'CSRTRACK':
                        self.CSRTrack.runCSRTrack(filename)
                        self.CSRTrack.convertCSRTrackOutput(f)
            else:
                print 'File does not exist! - ', f

class ASTRA(object):

    def __init__(self, parent=None, directory='test'):
        super(ASTRA, self).__init__()
        self.subdir = directory
        self.parent = parent
        self.astraCommand = ['astra']

    def runASTRA(self, filename=''):
        """Run the ASTRA program with input 'filename'"""
        command = self.astraCommand + [filename]
        with open(os.devnull, "w") as f:
            subprocess.call(command, stdout=f, cwd=self.subdir)

    def defineASTRACommand(self,command=['astra']):
        """Modify the defined ASTRA command variable"""
        self.astraCommand = command

    def setInitialDistribution(self, filename='../1k-250pC-76fsrms-1mm_TE09fixN12.ini'):
        """Modify the 'initial_distribution' global setting"""
        self.parent.globalSettings['initial_distribution'] = filename

    def createInitialDistribution(self, npart=1000, charge=250, generatorCommand=None, generatorFile=None):
        """Create an initiail dostribution of 'npart' particles of 'charge' pC"""
        self.parent.globalSettings['npart'] = npart
        self.parent.globalSettings['charge'] = charge/1000.0
        if generatorFile is None:
            if self.parent.generatorFile is not None:
                generatorFile = self.parent.generatorFile
            else:
                generatorFile = 'generator.in'
        astragen = GenPart.ASTRAGenerator(self.subdir, charge, npart, overwrite=self.parent.overwrite, generatorFile=generatorFile)
        if not generatorCommand is None:
            astragen.defineGeneratorCommand(generatorCommand)
        elif os.name == 'nt':
            astragen.defineGeneratorCommand(['generator_7June2007'])
        else:
            astragen.defineGeneratorCommand(['/opt/ASTRA/generator.sh'])
        inputfile = astragen.generateBeam()
        self.setInitialDistribution(inputfile)
        scgrid = getGrids(npart)
        self.parent.globalSettings['SC_2D_Nrad'] = max([scgrid.gridSizes,4])
        self.parent.globalSettings['SC_2D_Nlong'] = max([scgrid.gridSizes,4])
        for scvar in ['SC_3D_Nxf','SC_3D_Nyf','SC_3D_Nzf']:
            self.parent.globalSettings[scvar] = scgrid.gridSizes

    def createASTRAChicane(self, group, dipoleangle=None, width=0.2, gap=0.02):
        """Create a 4 dipole chicane in ASTRA with the correct edge points"""
        chicanetext = ''
        dipoles = self.parent.getGroup(group)
        if not dipoleangle is None:
            dipoleangle = float(dipoleangle)
            dipoles = [self.parent.setDipoleAngle(d, dipoleangle) for d in dipoles]
        dipoles = self.parent.createDrifts(dipoles)
        dipolepos, localXYZ = self.parent.elementPositions(dipoles)
        dipolepos = list(chunks(dipolepos,2))
        corners = [0,0,0,0]
        dipoleno = 0
        for elems in dipolepos:
            dipoleno += 1
            p1, psi1, nameelem1 = elems[0]
            p2, psi2, nameelem2 = elems[1]
            # print 'psi = ', psi1, psi2
            name, d = nameelem1
            angle1 = getParameter(nameelem1[1],'angle')
            angle2 = getParameter(nameelem2[1],'angle')
            length = getParameter(d,'length')
            e1 = getParameter(d,'entrance_edge_angle')
            e2 = getParameter(d,'exit_edge_angle')
            width = getParameter(d,'width',default=width)
            gap = getParameter(d,'gap',default=gap)
            rbend = 1 if d['type'] == 'rdipole' else 0
            rho = d['length']/angle1 if 'length' in d and abs(angle1) > 1e-9 else 0
            theta = -1*psi1-e1-rbend*np.sign(rho)*angle1/2.0
            corners[0] = np.array(map(add,np.transpose(p1),np.dot([-width*length,0,0], rotationMatrix(theta))))[0,0]
            corners[3] = np.array(map(add,np.transpose(p1),np.dot([width*length,0,0], rotationMatrix(theta))))[0,0]
            theta = -1*psi2+e2-rbend*np.sign(rho)*angle1/2.0
            corners[1] = np.array(map(add,np.transpose(p2),np.dot([-width*length,0,0], rotationMatrix(theta))))[0,0]
            corners[2] = np.array(map(add,np.transpose(p2),np.dot([width*length,0,0], rotationMatrix(theta))))[0,0]
            dipolenostr = str(dipoleno)
            chicanetext += "D_Type("+dipolenostr+")='horizontal',\n"+\
            "D_Gap(1,"+dipolenostr+")="+str(gap)+",\n"+\
            "D_Gap(2,"+dipolenostr+")="+str(gap)+",\n"+\
            "D1("+dipolenostr+")=("+str(corners[3][0])+","+str(corners[3][2])+"),\n"+\
            "D3("+dipolenostr+")=("+str(corners[2][0])+","+str(corners[2][2])+"),\n"+\
            "D4("+dipolenostr+")=("+str(corners[1][0])+","+str(corners[1][2])+"),\n"+\
            "D2("+dipolenostr+")=("+str(corners[0][0])+","+str(corners[0][2])+"),\n"+\
            "D_radius("+dipolenostr+")="+str(rho)+"\n"
        return chicanetext

    def createASTRAQuad(self, quadname, n=1):
        """Create an ASTRA quadrupole string"""
        quad        = self.parent.getElement(quadname)
        k1          = str(getParameter(quad,'k1'))
        length      = str(getParameter(quad,'length'))
        x,y,s       =     getParameter(quad,'position_start')
        bore        = str(getParameter(quad,'bore_size', default=0.01))
        smooth      = str(getParameter(quad,'smooth', default=3))

        quadtext = 'Q_K('+str(n)+')='+k1+', Q_length('+str(n)+')='+length+',\n'+\
        'Q_pos('+str(n)+')='+str(s)+', Q_smooth('+str(n)+')='+smooth+', Q_bore('+str(n)+')='+bore+'\n'
        return quadtext

    def createASTRASolenoid(self, solname, n=1):
        """Create an ASTRA solenoid string"""
        sol         = self.parent.getElement(solname)
        definition  = str(getParameter(sol,'field_definition'))
        length      = str(getParameter(sol,'length'))
        x,y,s       =     getParameter(sol,'position_start')
        amplitude   = str(getParameter(sol,'field_amplitude'))
        smooth      = str(getParameter(sol,'smooth', default=10))

        soltext = 'FILE_BFieLD('+str(n)+')=\''+definition+'\',\nMaxB('+str(n)+')='+amplitude+',\n'+\
        'S_pos('+str(n)+')='+str(s)+', S_xoff('+str(n)+')='+str(x)+', S_yoff('+str(n)+')='+str(y)+', S_smooth('+str(n)+')='+smooth+'\n'
        for var in ASTRARules['SOLENOID']:
            soltext += createOptionalString(sol, var, n)
        return soltext

    def createASTRACavity(self, cavname, n=1):
        """Create an ASTRA cavity string"""
        cav         = self.parent.getElement(cavname)
        definition  = str(getParameter(cav,'field_definition'))
        length      =     getParameter(cav,'length')
        x,y,s       =     getParameter(cav,'position_start')
        amplitude   = str(float(getParameter(cav,'field_amplitude'))/1e6)
        frequency   = str(float(getParameter(cav,'frequency'))/1e9)
        phase       = str(getParameter(cav,'phase'))
        cells       =     getParameter(cav,'number_of_cells')
        celllength  =     getParameter(cav,'cell_length')
        if cells is 0 and celllength > 0:
                cells = round((length-celllength)/celllength)
                cells = cells - (cells % 3)
        smooth      = str(getParameter(cav,'smooth', default=10))

        cavtext = 'FILE_EFieLD('+str(n)+')=\''+definition+'\', Nue('+str(n)+')='+str(frequency)+',\n'+\
        'MaxE('+str(n)+')='+amplitude+', Phi('+str(n)+')='+phase+', \n'+\
        'C_pos('+str(n)+')='+str(s)+', C_xoff('+str(n)+')='+str(x)+', C_yoff('+str(n)+')='+str(y)+', C_smooth('+str(n)+')='+smooth
        if cells > 0:
            cavtext += ', C_numb('+str(n)+')='+str(int(cells))+'\n'
        cavtext += '\n'
        for var in ASTRARules['CAVITY']:
            cavtext += createOptionalString(cav, var, n)
        return cavtext

    def createASTRAScreen(self, screenname, n=1):
        """Create an ASTRA screen string"""
        screen         = self.parent.getElement(screenname)
        x,y,s          =     getParameter(screen,'position_start')

        screentext = 'Screen('+str(n)+')='+str(s)+'\n'
        return screentext

    def formatASTRAStartElement(self, name):
        return str(int(round(self.parent.elements[name]['position_end'][2]*100))).zfill(4)

    def createASTRANewRunBlock(self, settings={}, input={}, output={}):
        """Create an ASTRA NEWRUN Block string"""
        title           = str(getParameter(settings,'title',default='trial'))
        runno           = str(getParameter(settings,'run_no',default=1))
        loop            = str(getParameter(settings,'Loop',default=False))
        lprompt         = str(getParameter(settings,'Lprompt',default=False))
        distribution    = str(getParameter(input,'particle_definition',default=''))
        if distribution == 'initial_distribution':
            distribution = self.parent.globalSettings['initial_distribution']
        else:
            regex = re.compile('\$(.*)\$')
            s = re.search(regex, distribution)
            if s:
                distribution = re.sub(regex, self.formatASTRAStartElement(eval(s.group(1))), distribution)
        # print 'qbunch = ', getParameter([self.globalSettings,settings],'total_charge',default=250)
        Qbunch          = str(getParameter([self.parent.globalSettings,settings],'total_charge',default=250))
        zstart          =     getParameter(settings,'zstart',default=0)
        zstop           =     getParameter(settings,'zstop',default=0)
        accuracy        = str(getParameter([self.parent.globalSettings,settings],'accuracy',default=4))
        highres = True if accuracy > 4 else False

        newruntext = '&NEWRUN\n' +\
        ' Loop='+str(loop)+'\n' + \
        ' Lprompt='+str(lprompt)+'\n' + \
        ' Head=\''+str(title)+'\'\n' + \
        ' Run='+str(runno)+'\n' + \
        ' Distribution=\''+str(distribution)+'\'\n' + \
        ' high_res='+str(highres)+'\n' + \
        ' Qbunch='+str(Qbunch)+'\n'
        for var in ASTRARules['NEWRUN']:
            newruntext += createOptionalString([self.parent.globalSettings['ASTRAsettings'],settings], var)
        newruntext += '/\n'

        return newruntext

    def createASTRAOutputBlock(self, output={}, settings={}):
        """Create an ASTRA OUTPUT Block string"""

        screens = self.parent.getElementsBetweenS('screen', output=output)
        # print 'screens = ', screens

        zstart = getParameter(output,'zstart',default=None)
        if zstart is None:
            startelem = getParameter(output,'start_element',default=None)
            if startelem is None or startelem not in self.parent.elements:
                zstart = 0
            else:
                # print self.parent.elements[startelem]
                zstart = self.parent.elements[startelem]['position_start'][2]
                output['zstart'] = zstart
        zstop = getParameter(output,'zstop',default=None)
        if zstop is None:
            endelem = getParameter(output,'end_element',default=None)
            if endelem is None or endelem not in self.parent.elements:
                zstop = 0
            else:
                zstop = self.parent.elements[endelem]['position_end'][2]
                output['zstop'] = zstop
        outputtext = '&OUTPUT\n'
        for var in ASTRARules['OUTPUT']:
            outputtext += createOptionalString([self.parent.globalSettings['ASTRAsettings'],settings, output], var)
        for i,s in enumerate(screens):
            outputtext += ' '+self.createASTRAScreen(s,i+1)
        outputtext += '/\n'

        return outputtext

    def createASTRAChargeBlock(self, charge={}, settings={}):
        """Create an ASTRA CHARGE Block string"""
        loop        = str(getParameter(charge,'Loop',default=False))
        mode        = str(getParameter(charge,'space_charge_mode',default='2D'))
        mirror_charge = str(getParameter(charge,'mirror_charge',default=False))

        lspch       = False if mode == False else True
        lspch2d     = True if lspch and mode != '3D' else False
        lspch3d     = True if lspch and not lspch2d else False
        if lspch2d:
            nrad    = str(getParameter([charge,self.parent.globalSettings],'SC_2D_Nrad',default=6))
            nlong   = str(getParameter([charge,self.parent.globalSettings],'SC_2D_Nlong',default=6))
        else:
            nxf     = str(getParameter([charge,self.parent.globalSettings],'SC_3D_Nxf',default=6))
            nyf     = str(getParameter([charge,self.parent.globalSettings],'SC_3D_Nyf',default=6))
            nzf     = str(getParameter([charge,self.parent.globalSettings],'SC_3D_Nzf',default=6))

        chargetext = '&CHARGE\n' +\
        ' Loop='+str(loop)+'\n' + \
        ' LSPCH='+str(lspch)+'\n' + \
        ' LSPCH3D='+str(lspch3d)+'\n' + \
        ' Lmirror='+str(mirror_charge)+'\n'
        if lspch and lspch2d:
            chargetext += ' Nrad='+nrad+', Nlong_in='+nlong+'\n'
        elif lspch and lspch3d:
            chargetext += ' Nxf='+nxf+', Nyf='+nyf+', Nzf='+nzf+'\n'

        for var in ASTRARules['CHARGE']:
            chargetext += createOptionalString([self.parent.globalSettings['ASTRAsettings'], settings, charge], var)
        chargetext += '/\n'

        return chargetext

    def createASTRAScanBlock(self, scan={}, settings={}):
        """Create an ASTRA SCAN Block string"""
        loop        = str(getParameter(scan,'Loop',default=False))
        lscan       = str(getParameter(scan,'LScan',default=False))

        scantext = '&SCAN\n' +\
        ' Loop='+str(loop)+'\n' +\
        ' LScan='+str(lscan)+'\n'
        for var in ASTRARules['SCAN']:
            scantext += createOptionalString([self.parent.globalSettings['ASTRAsettings'], settings, scan], var)
        scantext += '/\n'

        return scantext

    def createASTRAApertureBlock(self, aperture={}, settings={}):
        """Create an ASTRA APERTURE Block string"""
        loop        = str(getParameter(aperture,'Loop',default=False))
        lapert      = str(getParameter(aperture,'LApert',default=False))

        aperturetext = '&APERTURE\n' +\
        ' Loop='+str(loop)+'\n' +\
        ' LApert='+str(lapert)+'\n'
        for var in ASTRARules['APERTURE']:
            aperturetext += createOptionalString([self.parent.globalSettings['ASTRAsettings'], settings, aperture], var)
        aperturetext += '/\n'

        return aperturetext

    def createASTRACavityBlock(self, cavity={}, output={}):
        """Create an ASTRA APERTURE Block string"""
        loop        = str(getParameter(cavity,'Loop',default=False))
        lefield        = str(getParameter(cavity,'LEField',default=True))

        cavitytext = '&CAVITY\n' +\
        ' Loop='+str(loop)+'\n' +\
        ' LEField='+str(lefield)+'\n'

        cavities = self.parent.getElementsBetweenS('cavity', output=output)

        for i,s in enumerate(cavities):
            cavitytext += ' '+self.createASTRACavity(s,i+1)
        cavitytext += '/\n'

        return cavitytext

    def createASTRASolenoidBlock(self, solenoid={}, output={}):
        """Create an ASTRA SOLENOID Block string"""
        loop        = str(getParameter(solenoid,'Loop',default=False))
        lbfield        = str(getParameter(solenoid,'LBField',default=True))


        solenoidtext = '&SOLENOID\n' +\
        ' Loop='+str(loop)+'\n' +\
        ' LBField='+str(lbfield)+'\n'

        solenoids = self.parent.getElementsBetweenS('solenoid', output=output)

        for i,s in enumerate(solenoids):
            solenoidtext += ' '+self.createASTRASolenoid(s,i+1)
        solenoidtext += '/\n'

        return solenoidtext

    def createASTRAQuadrupoleBlock(self, quad={}, output={}):
        """Create an ASTRA QUADRUPOLE Block string"""
        loop        = str(getParameter(quad,'Loop',default=False))
        lquad        = str(getParameter(quad,'LQuad',default=True))

        quadrupoletext = '&QUADRUPOLE\n' +\
        ' Loop='+str(loop)+'\n' +\
        ' LQuad='+str(lquad)+'\n'

        quadrupoles = self.parent.getElementsBetweenS('quadrupole', output=output)

        for i,s in enumerate(quadrupoles):
            quadrupoletext += ' '+self.createASTRAQuad(s,i+1)
        quadrupoletext += '/\n'

        return quadrupoletext

    def createASTRAChicaneBlock(self, groups, dipoles):
        """Create an ASTRA DIPOLE Block string"""
        loop = False
        ldipole = False

        # print 'groups = ', groups

        for g in groups:
            if g in self.parent.groups:
                # print 'group!'
                if self.parent.groups[g]['type'] == 'chicane':
                    if all([i for i in self.parent.groups[g] if i in dipoles]):
                        ldipole = True

        dipoletext = '&DIPOLE\n' +\
        ' Loop='+str(loop)+'\n' +\
        ' LDipole='+str(ldipole)+'\n'

        for g in groups:
            if g in self.parent.groups:
                # print 'group!'
                if self.parent.groups[g]['type'] == 'chicane':
                    if all([i for i in self.parent.groups[g] if i in dipoles]):
                        dipoletext += self.createASTRAChicane(g, **groups[g])

        dipoletext += '/\n'

        return dipoletext

    def createASTRAFileText(self, file):
        settings = self.parent.getFileSettings(file,'ASTRA_Options')
        input = self.parent.getFileSettings(file,'input')
        output = self.parent.getFileSettings(file,'output')
        charge = self.parent.getFileSettings(file,'charge')
        scan = self.parent.getFileSettings(file,'scan')
        aperture = self.parent.getFileSettings(file,'aperture')
        cavity = self.parent.getFileSettings(file,'cavity')
        solenoid = self.parent.getFileSettings(file,'solenoid')
        quadrupole = self.parent.getFileSettings(file,'quadrupole')

        dipoles = self.parent.getElementsBetweenS('dipole', output)
        groups = self.parent.getFileSettings(file,'groups')

        astrafiletext = ''
        astrafiletext += self.createASTRANewRunBlock(settings, input, output)
        astrafiletext += self.createASTRAOutputBlock(output, settings)
        astrafiletext += self.createASTRAChargeBlock(charge, settings)
        astrafiletext += self.createASTRAScanBlock(scan, settings)
        astrafiletext += self.createASTRAApertureBlock(aperture, settings)
        astrafiletext += self.createASTRACavityBlock(cavity, output)
        astrafiletext += self.createASTRASolenoidBlock(solenoid, output)
        astrafiletext += self.createASTRAQuadrupoleBlock(quadrupole, output)
        astrafiletext += self.createASTRAChicaneBlock(groups, dipoles)

        return astrafiletext

    def createASTRAFiles(self):
        for f in self.fileSettings.keys():
            filename = self.subdirectory+'/'+f+'.in'
            # print filename
            saveFile(filename, lines=self.createASTRAFileText(f))

    def runASTRAFiles(self, files=None):
        if isinstance(files, (list, tuple)):
            for f in files:
                if f in self.fileSettings.keys():
                    filename = f+'.in'
                    # print 'Running file: ', f
                    self.runASTRA(filename)
                else:
                    print 'File does not exist! - ', f
        else:
            for f in self.fileSettings.keys():
                filename = f+'.in'
                # print 'Running file: ', f
                self.runASTRA(filename)

class CSRTrack(object):

    def __init__(self, parent=None, directory='test'):
        super(CSRTrack, self).__init__()
        self.subdir = directory
        self.parent = parent
        self.beam = rbf.beam()
        self.CSRTrackCommand = ['csrtrack']

    def runCSRTrack(self, filename=''):
        """Run the CSRTrack program with input 'filename'"""
        copyfile(self.subdir+'/'+filename, self.subdir+'/'+'csrtrk.in')
        command = self.CSRTrackCommand + [filename]
        with open(os.devnull, "w") as f:
            subprocess.call(command, stdout=f, cwd=self.subdir)

    def convertCSRTrackOutput(self, f):
        output = self.parent.getFileSettings(f,'output')
        outputdistribution    = f+'.$output[\'end_element\']$.001'
        regex = re.compile('\$(.*)\$')
        s = re.search(regex, outputdistribution)
        if s:
            sub = self.parent.astra.formatASTRAStartElement(eval(s.group(1)))
            outputdistribution = re.sub(regex, sub, outputdistribution)

        options = self.parent.getFileSettings(f,'CSRTrack_Options')
        monitor = self.parent.getSettingsBlock(options,'monitor')
        # print monitor
        inputdistribution    = str(getParameter(monitor,'name',default=''))
        regex = re.compile('\$(.*)\$')
        s = re.search(regex, inputdistribution)
        if s:
            sub = self.parent.astra.formatASTRAStartElement(eval(s.group(1)))
            inputdistribution = re.sub(regex, sub, inputdistribution)
        # print 'conversion = ', self.subdir+'/'+inputdistribution, self.subdir+'/'+outputdistribution
        self.beam.convert_csrtrackfile_to_astrafile(self.subdir+'/'+inputdistribution, self.subdir+'/'+outputdistribution)

    def defineCSRTrackCommand(self, command=['csrtrack']):
        """Modify the defined ASTRA command variable"""
        self.CSRTrackCommand = command

    def setInitialDistribution(self, filename='../1k-250pC-76fsrms-1mm_TE09fixN12.ini'):
        """Modify the 'initial_distribution' global setting"""
        self.parent.globalSettings['initial_distribution'] = filename

    def createCSRTrackChicane(self, group, dipoleangle=None, width=0.2, gap=0.02):
        """Create a 4 dipole chicane in CSRTrack with the correct edge points"""
        chicanetext = ''
        dipoles = self.parent.getGroup(group)
        if not dipoleangle is None:
            dipoleangle = float(dipoleangle)
            dipoles = [self.parent.setDipoleAngle(d, dipoleangle) for d in dipoles]
        dipoles = self.parent.createDrifts(dipoles)
        dipolepos, localXYZ = self.parent.elementPositions(dipoles)
        dipolepos = list(chunks(dipolepos,2))
        corners = [0,0,0,0]
        dipoleno = 0
        lastangle = 0
        for elems in dipolepos:
            dipoleno += 1
            p1, psi1, nameelem1 = elems[0]
            p2, psi2, nameelem2 = elems[1]
            name, d = nameelem1
            angle1 = getParameter(nameelem1[1],'angle')
            angle2 = getParameter(nameelem2[1],'angle')
            length = getParameter(d,'length')
            rbend = 1 if d['type'] == 'rdipole' else 0
            e1 = getParameter(d,'entrance_edge_angle')
            e1 = e1 if rbend is 0 else e1 + angle1/2.0
            e2 = getParameter(d,'exit_edge_angle')
            e2 = e2 if rbend is 0 else e2 + angle2/2.0
            width = getParameter(d,'width',default=width)
            gap = getParameter(d,'gap',default=gap)
            rho = d['length']/angle1 if 'length' in d and abs(angle1) > 1e-9 else 0
            np.transpose(p1)
            chicanetext += """    dipole{
    position{rho="""+str(p1[2,0])+""", psi="""+str(chop(e1+psi1))+""", marker=d"""+str(dipoleno)+"""a}
    properties{r="""+str(rho)+"""}
    position{rho="""+str(p2[2,0])+""", psi="""+str(chop(e2-psi2))+""", marker=d"""+str(dipoleno)+"""b}
    }
"""
        return chicanetext

    def createCSRTrackParticlesBlock(self, particles={}, output={}):
        """Create an CSRTrack Particles Block string"""

        distribution    = str(getParameter(particles,'array',default=''))
        regex = re.compile('\$(.*)\$')
        s = re.search(regex, distribution)
        if s:
            sub = self.parent.astra.formatASTRAStartElement(eval(s.group(1)))
            distribution = re.sub(regex, sub, distribution)

        del particles['array']
        particlestext = ''

        particlestext += 'particles{\n'
        for param, val in particles.iteritems():
            particlestext+= str(param)+' = '+str(val)+'\n'
        particlestext+= str('array')+' = '+str(distribution)+'\n'
        particlestext+= '}\n'

        return particlestext

    def createCSRTrackForcesBlock(self, forces={}):
        """Create an CSRTrack forces Block string"""

        forcestext = ''

        forcestext += 'forces{\n'

        for param, val in forces.iteritems():
            forcestext+= str(param)+' = '+str(val)+'\n'

        forcestext+= '}\n'

        return forcestext

    def createCSRTrackTrackStepBlock(self, trackstep={}):
        """Create an CSRTrack trackstep Block string"""

        tracksteptext = ''

        tracksteptext += 'track_step{\n'

        for param, val in trackstep.iteritems():
            tracksteptext+= str(param)+' = '+str(val)+'\n'

        tracksteptext+= '}\n'

        return tracksteptext

    def createCSRTrackTrackerBlock(self, tracker={}):
        """Create an CSRTrack tracker Block string"""

        trackertext = ''

        trackertext += 'tracker{\n'

        for param, val in tracker.iteritems():
            trackertext+= str(param)+' = '+str(val)+'\n'

        trackertext+= '}\n'

        return trackertext

    def createCSRTrackMonitorBlock(self, monitor={}, output={}):
        """Create an CSRTrack monitor Block string"""

        monitortext = ''

        monitortext += 'monitor{\n'

        distribution    = str(getParameter(monitor,'name',default=''))
        regex = re.compile('\$(.*)\$')
        s = re.search(regex, distribution)
        if s:
            sub = self.parent.astra.formatASTRAStartElement(eval(s.group(1)))
            distribution = re.sub(regex, sub, distribution)

        monitor = {i:monitor[i] for i in monitor if i!='name'}
        self.monitorName = distribution
        monitortext+= str('name')+' = '+str(distribution)+'\n'
        for param, val in monitor.iteritems():
            monitortext+= str(param)+' = '+str(val)+'\n'

        monitortext+= '}\n'

        return monitortext

    def createCSRTrackOnlineMonitorBlock(self, omon={}):
        """Create an CSRTrack OnlineMonitor Block string"""

        omontext = ''

        for name, params in omon.iteritems():

            omontext += 'online_monitor{\n'+'name = '+name+',\n'

            for param, val in params.iteritems():
                omontext+= str(param)+' = '+str(val)+'\n'
            omontext+= '}\n'

        return omontext

    def createCSRTrackChicaneBlock(self, groups, dipoles):
        """Create an CSRTrack DIPOLE Block string"""
        loop = False
        ldipole = False

        for g in groups:
            if g in self.parent.groups:
                # print 'group!'
                if self.parent.groups[g]['type'] == 'chicane':
                    if all([i for i in self.parent.groups[g] if i in dipoles]):
                        ldipole = True

        dipoletext = "io_path{logfile = log.txt}\n\n    lattice{\n"

        for g in groups:
            if g in self.parent.groups:
                # print 'group!'
                if self.parent.groups[g]['type'] == 'chicane':
                    if all([i for i in self.parent.groups[g] if i in dipoles]):
                        dipoletext += self.createCSRTrackChicane(g, **groups[g])

        dipoletext += "    }\n"

        return dipoletext

    def createCSRTrackFileText(self, file):
        options = self.parent.getFileSettings(file,'CSRTrack_Options')
        # input = self.parent.getFileSettings(file,'input')
        output = self.parent.getFileSettings(file,'output')
        online_monitors = self.parent.getSettingsBlock(options,'online_monitors')
        particles = self.parent.getSettingsBlock(options,'particles')
        forces = self.parent.getSettingsBlock(options,'forces')
        trackstep = self.parent.getSettingsBlock(options,'track_step')
        tracker = self.parent.getSettingsBlock(options,'tracker')
        monitor = self.parent.getSettingsBlock(options,'monitor')

        dipoles = self.parent.getElementsBetweenS('dipole', output)
        groups = self.parent.getFileSettings(file,'groups')

        CSRTrackfiletext = ''
        CSRTrackfiletext += self.createCSRTrackChicaneBlock(groups, dipoles)
        CSRTrackfiletext += self.createCSRTrackOnlineMonitorBlock(online_monitors)
        CSRTrackfiletext += self.createCSRTrackParticlesBlock(particles, output)
        CSRTrackfiletext += self.createCSRTrackForcesBlock(forces)
        CSRTrackfiletext += self.createCSRTrackTrackStepBlock(trackstep)
        CSRTrackfiletext += self.createCSRTrackTrackerBlock(tracker)
        CSRTrackfiletext += self.createCSRTrackMonitorBlock(monitor, output)
        CSRTrackfiletext += 'exit\n'
        # print 'CSRTrackfiletext = ', CSRTrackfiletext
        return CSRTrackfiletext

    def createCSRTrackFiles(self):
        for f in self.fileSettings.keys():
            filename = self.subdirectory+'/'+f+'.in'
            # print filename
            saveFile(filename, lines=self.createCSRTrackFileText(f))

    def runCSRTrackFiles(self, files=None):
        if isinstance(files, (list, tuple)):
            for f in files:
                if f in self.fileSettings.keys():
                    filename = f+'.in'
                    # print 'Running file: ', f
                    self.runCSRTrack(filename)
                else:
                    print 'File does not exist! - ', f
        else:
            for f in self.fileSettings.keys():
                filename = f+'.in'
                # print 'Running file: ', f
                self.runCSRTrack(filename)
