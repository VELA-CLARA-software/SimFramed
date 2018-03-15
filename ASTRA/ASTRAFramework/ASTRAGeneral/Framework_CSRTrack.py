import yaml, collections, subprocess, os, math, re, sys
import numpy as np
import read_beam_file as rbf
from FrameworkHelperFunctions import *

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
