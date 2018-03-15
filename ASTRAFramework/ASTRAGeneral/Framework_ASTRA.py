import yaml, collections, subprocess, os, math, re, sys, copy
import numpy as np
import ASTRAGenerator as GenPart
from ASTRARules import *
from getGrids import *
from FrameworkHelperFunctions import *
from operator import add

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

    def rotateAndOffset(self, start_pos, offset, theta):
        rotation_matrix = np.array([[np.cos(theta), 0, np.sin(theta)], [0, 1, 0], [-1*np.sin(theta), 0, np.cos(theta)]])
        return chop(np.dot(np.array(start_pos)-np.array(offset), rotation_matrix))

    def createASTRADipole(self, dipolename, n=1, width=0.2, gap=0.02, plane='horizontal'):
        """Create an ASTRA dipole string"""
        dipole        = self.parent.getElement(dipolename)
        length      = getParameter(dipole,'length')
        e1          = getParameter(dipole,'entrance_edge_angle')
        e2          = getParameter(dipole,'exit_edge_angle')
        width       = getParameter(dipole,'width',default=width)
        gap         = getParameter(dipole,'gap',default=gap)
        plane       = getParameter(dipole,'plane',default=plane)
        angle       = getParameter(dipole,'angle')
        x,y,z       = getParameter(dipole,'position_start')

        corners = [0,0,0,0]
        dipoles = self.parent.createDrifts([[dipolename, dipole]], zerolengthdrifts=True)
        dipolepos, localXYZ = self.parent.elementPositions(dipoles)#, startpos=[x,y,z])
        dipolepos = list(chunks(dipolepos,2))[0]
        p1, psi1, nameelem1 = dipolepos[0]
        p2, psi2, nameelem2 = dipolepos[1]
        rbend = 1 if getParameter(dipole,'type') == 'rdipole' else 0
        rho = getParameter(dipole,'length')/angle if getParameter(dipole,'length') is not None and abs(angle) > 1e-9 else 0
        theta = -1*psi1-e1-rbend*np.sign(rho)*angle/2.0
        corners[0] = np.array(map(add,np.transpose(p1),np.dot([-width*length,0,0], rotationMatrix(theta))))[0,0]
        corners[0] = self.rotateAndOffset(corners[0], self.global_offset, self.global_rotation)
        corners[3] = np.array(map(add,np.transpose(p1),np.dot([width*length,0,0], rotationMatrix(theta))))[0,0]
        corners[3] = self.rotateAndOffset(corners[3], self.global_offset, self.global_rotation)
        theta = -1*psi2+e2-rbend*np.sign(rho)*angle/2.0
        corners[1] = np.array(map(add,np.transpose(p2),np.dot([-width*length,0,0], rotationMatrix(theta))))[0,0]
        corners[1] = self.rotateAndOffset(corners[1], self.global_offset, self.global_rotation)
        corners[2] = np.array(map(add,np.transpose(p2),np.dot([width*length,0,0], rotationMatrix(theta))))[0,0]
        corners[2] = self.rotateAndOffset(corners[2], self.global_offset, self.global_rotation)

        # print 'corners = ', corners

        dipoletext = "D_Type("+str(n)+")='"+str(plane)+"',\n"+\
        "D_Gap(1,"+str(n)+")="+str(gap)+",\n"+\
        "D_Gap(2,"+str(n)+")="+str(gap)+",\n"+\
        "D1("+str(n)+")=("+str(corners[3][0])+","+str(corners[3][2])+"),\n"+\
        "D3("+str(n)+")=("+str(corners[2][0])+","+str(corners[2][2])+"),\n"+\
        "D4("+str(n)+")=("+str(corners[1][0])+","+str(corners[1][2])+"),\n"+\
        "D2("+str(n)+")=("+str(corners[0][0])+","+str(corners[0][2])+"),\n"+\
        "D_radius("+str(n)+")="+str(rho)+"\n"

        return dipoletext

    def createASTRAQuad(self, quadname, n=1):
        """Create an ASTRA quadrupole string"""
        quad        = self.parent.getElement(quadname)
        k1          = str(getParameter(quad,'k1'))
        length      = str(getParameter(quad,'length'))
        x,y,z       =     getParameter(quad,'position_start')
        x,y,z =  self.rotateAndOffset([x,y,z], self.global_offset, self.global_rotation)
        bore        = str(getParameter(quad,'bore_size', default=0.01))
        smooth      = str(getParameter(quad,'smooth', default=3))

        quadtext = 'Q_K('+str(n)+')='+k1+', Q_length('+str(n)+')='+length+',\n'+\
        'Q_pos('+str(n)+')='+str(z)+', Q_smooth('+str(n)+')='+smooth+', Q_bore('+str(n)+')='+bore+'\n'
        return quadtext

    def createASTRASolenoid(self, solname, n=1):
        """Create an ASTRA solenoid string"""
        sol         = self.parent.getElement(solname)
        definition  = str(getParameter(sol,'field_definition'))
        length      = str(getParameter(sol,'length'))
        x,y,z       =     getParameter(sol,'position_start')
        x,y,z =  self.rotateAndOffset([x,y,z], self.global_offset, self.global_rotation)
        amplitude   = str(getParameter(sol,'field_amplitude'))
        smooth      = str(getParameter(sol,'smooth', default=10))

        soltext = 'FILE_BFieLD('+str(n)+')=\''+definition+'\',\nMaxB('+str(n)+')='+amplitude+',\n'+\
        'S_pos('+str(n)+')='+str(z)+', S_xoff('+str(n)+')='+str(x)+', S_yoff('+str(n)+')='+str(y)+', S_smooth('+str(n)+')='+smooth+'\n'
        for var in ASTRARules['SOLENOID']:
            soltext += createOptionalString(sol, var, n)
        return soltext

    def createASTRACavity(self, cavname, n=1):
        """Create an ASTRA cavity string"""
        cav         = self.parent.getElement(cavname)
        definition  = str(getParameter(cav,'field_definition'))
        length      =     getParameter(cav,'length')
        x,y,z       =     getParameter(cav,'position_start')
        x,y,z =  self.rotateAndOffset([x,y,z], self.global_offset, self.global_rotation)
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
        'C_pos('+str(n)+')='+str(z)+', C_xoff('+str(n)+')='+str(x)+', C_yoff('+str(n)+')='+str(y)+', C_smooth('+str(n)+')='+smooth
        if cells > 0:
            cavtext += ', C_numb('+str(n)+')='+str(int(cells))+'\n'
        cavtext += '\n'
        for var in ASTRARules['CAVITY']:
            cavtext += createOptionalString(cav, var, n)
        return cavtext

    def createASTRAScreen(self, screenname, n=1):
        """Create an ASTRA screen string"""
        screen         = self.parent.getElement(screenname)
        x,y,z          =     getParameter(screen,'position_start')
        x,y,z =  self.rotateAndOffset([x,y,z], self.global_offset, self.global_rotation)

        screentext = 'Screen('+str(n)+')='+str(z)+'\n'
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
        # zstart          =     getParameter(settings,'zstart',default=0)
        # zstop           =     getParameter(settings,'zstop',default=0)
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

    def createASTRAOutputBlock(self, originaloutput={}, settings={}):
        """Create an ASTRA OUTPUT Block string"""

        output = copy.deepcopy(originaloutput)
        screens = self.parent.getElementsBetweenS('screen', output=output)
        # print 'screens = ', screens

        zstart = getParameter(output,'zstart',default=None)
        if zstart is None:
            startelem = getParameter(output,'start_element',default=None)
            if startelem is None or startelem not in self.parent.elements:
                zstart = [0,0,0]
            else:
                # print self.parent.elements[startelem]
                zstart = self.parent.elements[startelem]['position_start']
                originaloutput['zstart'] = zstart[2]
        elif not isinstance(zstart, (list, tuple)):
            zstart = [0,0, zstart]
        zstop = getParameter(output,'zstop',default=None)
        if zstop is None:
            endelem = getParameter(output,'end_element',default=None)
            if endelem is None or endelem not in self.parent.elements:
                zstop = [0,0,0]
            else:
                zstop = self.parent.elements[endelem]['position_end']
                originaloutput['zstop'] = zstop[2]
        elif not isinstance(zstop, (list, tuple)):
            zstop = [0,0,zstop]
        # print 'zstart = ', zstart
        # print 'zstop = ', zstop
        zstart = self.rotateAndOffset(zstart, self.global_offset, self.global_rotation)
        output['zstart'] = zstart[2]
        # print 'zstop = ', self.parent.elements[endelem]['position_end']
        zstop = self.rotateAndOffset(zstop, self.global_offset, self.global_rotation)
        # print 'zstop after = ', zstop
        output['zstop'] = zstop[2]
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
        lscan = True if len(scan) > 0 else False
        lscan       = str(getParameter(scan,'LScan',default=lscan))

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
        lapert = True if len(aperture) > 0 else False
        lapert      = str(getParameter(aperture,'LApert',default=lapert))

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
        lefield = True if len(cavity) > 0 else False
        lefield        = str(getParameter(cavity,'LEField',default=lefield))

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
        lbfield = True if len(solenoid) > 0 else False
        lbfield        = str(getParameter(solenoid,'LBField',default=lbfield))


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
        lquad = True if len(quad) > 0 else False
        lquad        = str(getParameter(quad,'LQuad',default=lquad))

        quadrupoletext = '&QUADRUPOLE\n' +\
        ' Loop='+str(loop)+'\n' +\
        ' LQuad='+str(lquad)+'\n'

        quadrupoles = self.parent.getElementsBetweenS('quadrupole', output=output)

        for i,s in enumerate(quadrupoles):
            quadrupoletext += ' '+self.createASTRAQuad(s,i+1)
        quadrupoletext += '/\n'

        return quadrupoletext

    def createASTRADipoleBlock(self, dipole={}, output={}, groups={}):
        """Create an ASTRA DIPOLE Block string"""
        loop        = str(getParameter(dipole,'Loop',default=False))
        ldipole = True if len(dipole) > 0 else False
        ldipole     = str(getParameter(dipole,'LDipole', default=ldipole))

        dipoletext = '&DIPOLE\n' +\
        ' Loop='+str(loop)+'\n' +\
        ' LDipole='+str(ldipole)+'\n'

        dipoles = self.parent.getElementsBetweenS('dipole', output=output)

        for g in groups:
            if g in self.parent.groups:
                if self.parent.groups[g]['type'] == 'chicane':
                    if all([i for i in self.parent.groups[g] if i in dipoles]):
                        dipoles = [i for i in dipoles if i not in self.parent.groups[g]]
                        dipoletext += self.createASTRAChicane(g, **groups[g])

        for i,s in enumerate(dipoles):
            dipoletext += ' '+self.createASTRADipole(s,i+1)

        dipoletext += '/\n'

        return dipoletext

    # def createASTRAChicaneBlock(self, groups, dipoles):
    #     """Create an ASTRA DIPOLE Block string for a chicane"""
    #     for g in groups:
    #         if g in self.parent.groups:
    #             # print 'group!'
    #             if self.parent.groups[g]['type'] == 'chicane':
    #                 if all([i for i in self.parent.groups[g] if i in dipoles]):
    #                     ldipole = True
    #
    #     for g in groups:
    #         if g in self.parent.groups:
    #             # print 'group!'
    #             if self.parent.groups[g]['type'] == 'chicane':
    #                 if all([i for i in self.parent.groups[g] if i in dipoles]):
    #                     dipoletext += self.createASTRAChicane(g, **groups[g])
    #
    #     dipoletext += '/\n'
    #
    #     return dipoletext

    def createASTRAFileText(self, file):
        settings = self.parent.getFileSettings(file,'ASTRA_Options')
        self.global_offset = self.parent.getFileSettings(file,'global_offset', [0,0,0])
        if isinstance(self.global_offset,(str)):
            regex = re.compile('\$(.*)\$')
            s = re.search(regex, self.global_offset)
            if s:
                self.global_offset = eval(s.group(1))
        self.global_rotation = self.parent.getFileSettings(file,'global_rotation', 0)
        if isinstance(self.global_rotation,(str)):
            regex = re.compile('\$(.*)\$')
            s = re.search(regex, self.global_rotation)
            if s:
                self.global_rotation = eval(s.group(1))
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
        astrafiletext += self.createASTRADipoleBlock(dipoles, output, groups)

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
