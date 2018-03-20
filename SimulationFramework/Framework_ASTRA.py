import yaml, collections, subprocess, os, math, re, sys, copy, ast
import numpy as np
import ASTRAGenerator as GenPart
from ASTRARules import *
from getGrids import *
from FrameworkHelperFunctions import *
import SimulationFramework.Modules.read_beam_file as rbf
from operator import add

class ASTRA(object):

    beam = rbf.beam()

    def __init__(self, framework=None, directory='test'):
        super(ASTRA, self).__init__()
        self.subdir = directory
        self.framework = framework
        self.astraCommand = ['astra']

    def runASTRA(self, filename=''):
        """Run the ASTRA program with input 'filename'"""
        command = self.astraCommand + [filename]
        with open(os.devnull, "w") as f:
            subprocess.call(command, stdout=f, cwd=self.subdir)

    def postProcesssASTRA(self):
        if hasattr(self, 'screens'):
            for s in self.screens:
                self.convert_astra_beam_to_HDF5_beam(self.subdir, self.filename, s, self.runno)
        self.convert_astra_beam_to_HDF5_beam(self.subdir, self.filename, self.zstop, self.runno)

    def convert_astra_beam_to_HDF5_beam(self, subdir, filename, pos, runno=1):
        astrabeamfilename = filename + '.' + str(int(round(pos[2]*100))).zfill(4) + '.' + str(runno.zfill(3))
        self.beam.read_astra_beam_file(subdir + '/' + astrabeamfilename)
        HDF5filename = filename + '.' + str(int(round(pos[2]*100))).zfill(4) + '.hdf5'
        self.beam.write_HDF5_beam_file(subdir + '/' + HDF5filename, centered=False, sourcefilename=astrabeamfilename, pos=pos)

    def defineASTRACommand(self,command=['astra']):
        """Modify the defined ASTRA command variable"""
        self.astraCommand = command

    def setInitialDistribution(self, filename='../1k-250pC-76fsrms-1mm_TE09fixN12.ini'):
        """Modify the 'initial_distribution' global setting"""
        self.framework.globalSettings['initial_distribution'] = filename

    def createInitialDistribution(self, npart=1000, charge=250, generatorCommand=None, generatorFile=None):
        """Create an initiail dostribution of 'npart' particles of 'charge' pC"""
        self.framework.globalSettings['npart'] = npart
        self.framework.globalSettings['charge'] = charge/1000.0
        if generatorFile is None:
            if self.framework.generatorFile is not None:
                generatorFile = self.framework.generatorFile
            else:
                generatorFile = 'generator.in'
        astragen = GenPart.ASTRAGenerator(self.subdir, charge, npart, overwrite=self.framework.overwrite, generatorFile=generatorFile)
        if not generatorCommand is None:
            astragen.defineGeneratorCommand(generatorCommand)
        elif os.name == 'nt':
            astragen.defineGeneratorCommand(['../../ASTRA/generator_7June2007'])
        else:
            astragen.defineGeneratorCommand(['/opt/ASTRA/generator.sh'])
        inputfile = astragen.generateBeam()
        self.setInitialDistribution(inputfile)
        scgrid = getGrids(npart)
        self.framework.globalSettings['SC_2D_Nrad'] = max([scgrid.gridSizes,4])
        self.framework.globalSettings['SC_2D_Nlong'] = max([scgrid.gridSizes,4])
        for scvar in ['SC_3D_Nxf','SC_3D_Nyf','SC_3D_Nzf']:
            self.framework.globalSettings[scvar] = scgrid.gridSizes

    def createASTRAChicane(self, group, dipoleangle=None, width=0.2, gap=0.02):
        """Create a 4 dipole chicane in ASTRA with the correct edge points"""
        chicanetext = ''
        dipoles = self.framework.getGroup(group)
        if not dipoleangle is None:
            dipoleangle = float(dipoleangle)
            dipoles = [self.framework.setDipoleAngle(d, dipoleangle) for d in dipoles]
        dipoles = self.framework.createDrifts(dipoles)
        dipolepos, localXYZ = self.framework.elementPositions(dipoles)
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
        return chop(np.dot(np.array(start_pos)-np.array(offset), rotation_matrix), 1e-6)

    def createASTRACorrector(self, kickername, n=1, width=0.2, gap=0.02):
        """Create an ASTRA dipole string"""
        kicker        = self.framework.getElement(kickername)
        length      = getParameter(kicker,'length')
        e1          = getParameter(kicker,'entrance_edge_angle')
        e2          = getParameter(kicker,'exit_edge_angle')
        width       = getParameter(kicker,'width', default=width)
        gap         = getParameter(kicker,'gap', default=gap)
        plane       = getParameter(kicker,'plane', default='combined')
        angle       = getParameter(kicker,'angle', default=0)
        x,y,z       = getParameter(kicker,'position_start')
        strength_H  = getParameter(kicker,'strength_H')# strength for horizontal kicker
        strength_V  = getParameter(kicker,'strength_V')# strength for vertical kicker

        dipoletext = ""

        corners = [0,0,0,0]
        kickers = self.framework.createDrifts([[kickername, kicker]], zerolengthdrifts=True)
        kickerpos, localXYZ = self.framework.elementPositions(kickers, startangle=self.starting_rotation)
        kickerpos = list(chunks(kickerpos,2))[0]
        p1, psi1, nameelem1 = kickerpos[0]
        p2, psi2, nameelem2 = kickerpos[1]

        rbend = 0 #if getParameter(kicker,'type') == 'rdipole' else 0
        rho = getParameter(kicker,'length')/angle if getParameter(kicker,'length') is not None and abs(angle) > 1e-9 else 0
        theta = -1*psi1 - e1 - rbend*np.sign(rho)*angle/2.0
        corners[0] = np.array(map(add,np.transpose(p1),np.dot([-width,0,0], rotationMatrix(theta))))[0,0]
        corners[0] = self.rotateAndOffset(corners[0], self.global_offset, self.global_rotation)
        corners[3] = np.array(map(add,np.transpose(p1),np.dot([width,0,0], rotationMatrix(theta))))[0,0]
        corners[3] = self.rotateAndOffset(corners[3], self.global_offset, self.global_rotation)
        theta = -1*psi2+e2-rbend*np.sign(rho)*angle/2.0
        corners[1] = np.array(map(add,np.transpose(p2),np.dot([-width,0,0], rotationMatrix(theta))))[0,0]
        corners[1] = self.rotateAndOffset(corners[1], self.global_offset, self.global_rotation)
        corners[2] = np.array(map(add,np.transpose(p2),np.dot([width,0,0], rotationMatrix(theta))))[0,0]
        corners[2] = self.rotateAndOffset(corners[2], self.global_offset, self.global_rotation)

        if plane is 'horizontal' or plane is 'combined':
            dipoletext += "D_Type("+str(n)+")='horizontal',\n"+\
            "D_Gap(1,"+str(n)+")="+str(gap)+",\n"+\
            "D_Gap(2,"+str(n)+")="+str(gap)+",\n"+\
            "D1("+str(n)+")=("+str(corners[3][0])+","+str(corners[3][2])+"),\n"+\
            "D3("+str(n)+")=("+str(corners[2][0])+","+str(corners[2][2])+"),\n"+\
            "D4("+str(n)+")=("+str(corners[1][0])+","+str(corners[1][2])+"),\n"+\
            "D2("+str(n)+")=("+str(corners[0][0])+","+str(corners[0][2])+"),\n"+\
            "D_strength("+str(n)+")="+str(strength_H)+"\n"
        if plane is 'vertical' or plane is 'combined':
            dipoletext += "D_Type("+str(n+1)+")='vertical',\n"+\
            "D_Gap(1,"+str(n+1)+")="+str(gap)+",\n"+\
            "D_Gap(2,"+str(n+1)+")="+str(gap)+",\n"+\
            "D1("+str(n+1)+")=("+str(corners[3][0])+","+str(corners[3][2])+"),\n"+\
            "D3("+str(n+1)+")=("+str(corners[2][0])+","+str(corners[2][2])+"),\n"+\
            "D4("+str(n+1)+")=("+str(corners[1][0])+","+str(corners[1][2])+"),\n"+\
            "D2("+str(n+1)+")=("+str(corners[0][0])+","+str(corners[0][2])+"),\n"+\
            "D_strength("+str(n+1)+")="+str(strength_V)+"\n"

        return dipoletext

    def createASTRADipole(self, dipolename, n=1, width=0.2, gap=0.02, plane='horizontal'):
        """Create an ASTRA dipole string"""
        dipole        = self.framework.getElement(dipolename)
        length      = getParameter(dipole,'length')
        e1          = getParameter(dipole,'entrance_edge_angle')
        e2          = getParameter(dipole,'exit_edge_angle')
        width       = getParameter(dipole,'width',default=width)
        gap         = getParameter(dipole,'gap',default=gap)
        plane       = getParameter(dipole,'plane',default=plane)
        angle       = self.framework.getElement(dipolename,'angle', default=0)
        x,y,z       = getParameter(dipole,'position_start')

        corners = [0,0,0,0]
        dipoles = self.framework.createDrifts([[dipolename, dipole]], zerolengthdrifts=True)
        dipolepos, localXYZ = self.framework.elementPositions(dipoles, startangle=self.starting_rotation)
        dipolepos = list(chunks(dipolepos,2))[0]
        p1, psi1, nameelem1 = dipolepos[0]
        p2, psi2, nameelem2 = dipolepos[1]
        self.starting_rotation += angle
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
        "D2("+str(n)+")=("+str(corners[0][0])+","+str(corners[0][2])+"),\n"
        if abs(getParameter(dipole,'strength')) > 0:
            dipoletext += "D_strength("+str(n)+")="+str(getParameter(dipole,'strength'))+"\n"
        else:
            dipoletext += "D_radius("+str(n)+")="+str(rho)+"\n"

        return dipoletext

    def createASTRAQuad(self, quadname, n=1):
        """Create an ASTRA quadrupole string"""
        quad        = self.framework.getElement(quadname)
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
        sol         = self.framework.getElement(solname)
        definition  = str(getParameter(sol,'field_definition'))
        definition = self.framework.expand_substitution(definition,{'master_lattice_location': self.framework.master_lattice_location})
        definition = os.path.relpath(definition, './'+self.subdir).replace('\\','/').replace('//','/')

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
        cav         = self.framework.getElement(cavname)
        definition  = str(getParameter(cav,'field_definition'))
        definition = self.framework.expand_substitution(definition,{'master_lattice_location': self.framework.master_lattice_location})
        definition = os.path.relpath(definition, './'+self.subdir).replace('\\','/').replace('//','/')

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
        screen         = self.framework.getElement(screenname)
        x,y,z          =     getParameter(screen,'position_start')
        self.screens.append([x,y,z])
        x,y,z =  self.rotateAndOffset([x,y,z], self.global_offset, self.global_rotation)

        screentext = 'Screen('+str(n)+')='+str(z)+'\n'

        return screentext

    def formatASTRAZValue(self, z, fill=4):
        return str(int(round(z))).zfill(fill)

    def formatASTRAStartElement(self, name):
        return self.formatASTRAZValue(self.framework.elements[name]['position_end'][2]*100)

    def createASTRANewRunBlock(self, settings={}, input={}, output={}):
        """Create an ASTRA NEWRUN Block string"""
        title           = str(getParameter(settings,'title',default='trial'))
        runno           = str(getParameter(settings,'run_no',default=1))
        self.runno = runno
        loop            = str(getParameter(settings,'Loop',default=False))
        lprompt         = str(getParameter(settings,'Lprompt',default=False))
        distribution    = str(getParameter(input,'particle_definition',default=''))
        if distribution == 'initial_distribution':
            distribution = self.framework.globalSettings['initial_distribution']
        else:
            regex = re.compile('\$(.*)\$')
            s = re.search(regex, distribution)
            if s:
                distribution = re.sub(regex, self.formatASTRAStartElement(eval(s.group(1))), distribution)
        # print 'qbunch = ', getParameter([self.globalSettings,settings],'total_charge',default=250)
        Qbunch          = str(getParameter([self.framework.globalSettings,settings],'total_charge',default=250))
        # zstart          =     getParameter(settings,'zstart',default=0)
        # zstop           =     getParameter(settings,'zstop',default=0)
        accuracy        = str(getParameter([self.framework.globalSettings,settings],'accuracy',default=4))
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
            newruntext += createOptionalString([self.framework.globalSettings['ASTRAsettings'],settings], var)
        newruntext += '/\n'

        return newruntext

    def createASTRAOutputBlock(self, originaloutput={}, settings={}):
        """Create an ASTRA OUTPUT Block string"""

        output = copy.deepcopy(originaloutput)
        screens = self.framework.getElementsBetweenS('screen', output=output)
        self.screens = []# print 'screens = ', screens

        zstart = getParameter(output,'zstart',default=None)
        if zstart is None:
            startelem = getParameter(output,'start_element',default=None)
            if startelem is None or startelem not in self.framework.elements:
                zstart = [0,0,0]
            else:
                # print self.framework.elements[startelem]
                zstart = self.framework.elements[startelem]['position_start']
                originaloutput['zstart'] = zstart[2]
        elif not isinstance(zstart, (list, tuple)):
            zstart = [0,0, zstart]
        zstop = getParameter(output,'zstop',default=None)
        if zstop is None:
            endelem = getParameter(output,'end_element',default=None)
            if endelem is None or endelem not in self.framework.elements:
                zstop = [0,0,0]
            else:
                zstop = self.framework.elements[endelem]['position_end']
                originaloutput['zstop'] = zstop[2]
        elif not isinstance(zstop, (list, tuple)):
            zstop = [0,0,zstop]
        # print 'zstart = ', zstart
        # print 'zstop = ', zstop
        zstart = self.rotateAndOffset(zstart, self.global_offset, self.global_rotation)
        output['zstart'] = zstart[2]
        # print 'zstop = ', self.framework.elements[endelem]['position_end']
        zstop = self.rotateAndOffset(zstop, self.global_offset, self.global_rotation)
        # print 'zstop after = ', zstop
        output['zstop'] = zstop[2]

        self.zstop = zstop
        outputtext = '&OUTPUT\n'
        for var in ASTRARules['OUTPUT']:
            outputtext += createOptionalString([self.framework.globalSettings['ASTRAsettings'],settings, output], var)
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
            nrad    = str(getParameter([charge,self.framework.globalSettings],'SC_2D_Nrad',default=6))
            nlong   = str(getParameter([charge,self.framework.globalSettings],'SC_2D_Nlong',default=6))
        else:
            nxf     = str(getParameter([charge,self.framework.globalSettings],'SC_3D_Nxf',default=6))
            nyf     = str(getParameter([charge,self.framework.globalSettings],'SC_3D_Nyf',default=6))
            nzf     = str(getParameter([charge,self.framework.globalSettings],'SC_3D_Nzf',default=6))

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
            chargetext += createOptionalString([self.framework.globalSettings['ASTRAsettings'], settings, charge], var)
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
            scantext += createOptionalString([self.framework.globalSettings['ASTRAsettings'], settings, scan], var)
        scantext += '/\n'

        return scantext

    def createASTRAApertureBlock(self, aperture={}, settings={}):
        """Create an ASTRA APERTURE Block string"""

        loop        = str(getParameter(aperture,'Loop',default=False))
        lapert = True if len(aperture.keys()) > 0 else False
        lapert      = str(getParameter(aperture,'LApert',default=lapert))

        aperturetext = '&APERTURE\n' +\
        ' Loop='+str(loop)+'\n' +\
        ' LApert='+str(lapert)+'\n'
        for var in ASTRARules['APERTURE']:
            aperturetext += createOptionalString([self.framework.globalSettings['ASTRAsettings'], settings, aperture], var)
        aperturetext += '/\n'

        return aperturetext

    def createASTRACavityBlock(self, cavity={}, output={}):
        """Create an ASTRA APERTURE Block string"""
        cavities = self.framework.getElementsBetweenS('cavity', output=output)

        loop        = str(getParameter(cavity,'Loop',default=False))
        lefield = True if len(cavities) > 0 else False
        lefield        = str(getParameter(cavity,'LEField',default=lefield))

        cavitytext = '&CAVITY\n' +\
        ' Loop='+str(loop)+'\n' +\
        ' LEField='+str(lefield)+'\n'


        for i,s in enumerate(cavities):
            cavitytext += ' '+self.createASTRACavity(s,i+1)
        cavitytext += '/\n'

        return cavitytext

    def createASTRASolenoidBlock(self, solenoid={}, output={}):
        """Create an ASTRA SOLENOID Block string"""
        solenoids = self.framework.getElementsBetweenS('solenoid', output=output)

        loop        = str(getParameter(solenoid,'Loop',default=False))
        lbfield = True if len(solenoids) > 0 else False
        lbfield        = str(getParameter(solenoid,'LBField',default=lbfield))


        solenoidtext = '&SOLENOID\n' +\
        ' Loop='+str(loop)+'\n' +\
        ' LBField='+str(lbfield)+'\n'


        for i,s in enumerate(solenoids):
            solenoidtext += ' '+self.createASTRASolenoid(s,i+1)
        solenoidtext += '/\n'

        return solenoidtext

    def createASTRAQuadrupoleBlock(self, quad={}, output={}):
        """Create an ASTRA QUADRUPOLE Block string"""
        quadrupoles = self.framework.getElementsBetweenS('quadrupole', output=output)

        loop        = str(getParameter(quad,'Loop',default=False))
        lquad = True if len(quadrupoles) > 0 else False
        lquad        = str(getParameter(quad,'LQuad',default=lquad))

        quadrupoletext = '&QUADRUPOLE\n' +\
        ' Loop='+str(loop)+'\n' +\
        ' LQuad='+str(lquad)+'\n'


        for i,s in enumerate(quadrupoles):
            quadrupoletext += ' '+self.createASTRAQuad(s,i+1)
        quadrupoletext += '/\n'

        return quadrupoletext

    def createASTRADipoleBlock(self, dipole={}, output={}, groups={}):
        """Create an ASTRA DIPOLE Block string"""

        dipoles = self.framework.getElementsBetweenS('dipole', output=output)
        kickers = self.framework.getElementsBetweenS('kicker', output=output)
        zerokickers = []
        for k in kickers:
            if not abs(self.framework.getElement(k,'strength_H', 0)) > 0 and not abs(self.framework.getElement(k,'strength_V',0)) > 0:
                zerokickers.append(k)
        for k in zerokickers:
            kickers.remove(k)

        loop        = str(getParameter(dipole,'Loop',default=False))
        ldipole = True if len(dipoles) > 0 or len(kickers) > 0 else False
        ldipole     = str(getParameter(dipole,'LDipole', default=ldipole))

        dipoletext = '&DIPOLE\n' +\
        ' Loop='+str(loop)+'\n' +\
        ' LDipole='+str(ldipole)+'\n'

        for g in groups:
            if g in self.framework.groups:
                if self.framework.groups[g]['type'] == 'chicane':
                    if all([i for i in self.framework.groups[g] if i in dipoles]):
                        dipoles = [i for i in dipoles if i not in self.framework.groups[g]]
                        dipoletext += self.createASTRAChicane(g, **groups[g])
        counter = 1
        for i,s in enumerate(dipoles):
            dipoletext += ' '+self.createASTRADipole(s,counter)
            counter += 1

        # Add in correctors
        for i,s in enumerate(kickers):
            dipoletext += ' '+ self.createASTRACorrector(s, counter)
            counter += 2 # two dipole one horizontal one vertical

        dipoletext += '/\n'

        return dipoletext

    def createASTRAFileText(self, file):
        self.filename = file
        settings = self.framework.getFileSettings(file,'ASTRA_Options')
        output = self.framework.getFileSettings(file,'output')

        self.global_offset = self.framework.getFileSettings(file,'global_offset', [0,0,0])
        self.global_offset = self.framework.expand_substitution(self.global_offset)

        self.starting_rotation = self.framework.getFileSettings(file,'starting_rotation', 0)
        self.starting_rotation = self.framework.expand_substitution(self.starting_rotation)


        dipoles = self.framework.getElementsBetweenS('dipole', output)
        self.global_rotation = -self.starting_rotation
        for d in dipoles:
            self.global_rotation -= self.framework.getElement(d,'angle')

        input = self.framework.getFileSettings(file,'input')

        charge = self.framework.getFileSettings(file,'charge')
        scan = self.framework.getFileSettings(file,'scan')
        aperture = self.framework.getFileSettings(file,'aperture')
        cavity = self.framework.getFileSettings(file,'cavity')
        solenoid = self.framework.getFileSettings(file,'solenoid')
        quadrupole = self.framework.getFileSettings(file,'quadrupole')


        groups = self.framework.getFileSettings(file,'groups')

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
