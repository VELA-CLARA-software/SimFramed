from PyQt4.QtCore import QThread
import createBeam as cb
import createBeamline as cbl
import yaml
import os
from epics import caget, caput
# from sourceCode.SAMPLcore.SAMPLlab import PhysicalUnits
'''Assumes  you have the a OnlineModel in you path directory'''
from SimulationFramework import Framework


class Setup(QThread):
    # CONSTRUCTOR
    def __init__(self, V_MAG_Ctrl=None, C_S01_MAG_Ctrl=None,
                 C_S02_MAG_Ctrl=None, C2V_MAG_Ctrl=None, V_RF_Ctrl=None,
                 C_RF_Ctrl=None, L01_RF_Ctrl=None, messages=False):
        QThread.__init__(self)
        self.showMessages = messages
        self.V_MAG_Ctrl = V_MAG_Ctrl
        self.C_S01_MAG_Ctrl = C_S01_MAG_Ctrl
        self.C_S02_MAG_Ctrl = C_S02_MAG_Ctrl
        self.C2V_MAG_Ctrl = C2V_MAG_Ctrl
        self.V_RF_Ctrl = V_RF_Ctrl
        self.C_RF_Ctrl = C_RF_Ctrl
        self.L01_RF_Ctrl = L01_RF_Ctrl
        self.startElement = None
        self.stopElement = None
        self.initDistribFile = 'Null'
        self.initDistrib = None
        self.initCharge = 0.0
        self.pathway = Framework.Framework(subdir='.', overwrite='overwrite')

    # DESTRUCTOR
    def __del__(self):
        self.wait()

    def loadPathway(self):
        ans = False
        stream = file(str(os.path.abspath(__file__)).split('sampl')[0] +
                      "\\..\\..\\MasterLattice\\Lattices\\allPathways.yaml", 'r')
        settings = yaml.load(stream)
        for path in settings['pathways']:
            # print path
            self.pathway.loadSettings(filename='\\Lattices\\' + path)
            hasStart = False
            hasStop = False
            # print self.pathway.elements.keys()
            for element in self.pathway.elements.keys():
                #print element
                if element == self.startElement:
                    #print('found start')
                    hasStart = True
                if element == self.stopElement:
                    #print ('found stop')
                    hasStop = True
            if (hasStart and hasStop):
                ans = True
                print '    Loaded pathway: ', path
                return ans
            else:
                print '    Could not find a pathway...'
                print '    Exiting ...'
                return ans

    def go(self, startElement, stopElement, initDistrib, charge=0.25):
            self.startElement = startElement
            self.stopElement = stopElement
            self.initDistrib = initDistrib
            self.initCharge = charge  # in nC
            # Run in Thread
            self.start()

    def run(self):
        print('------------------------------')
        print('-------------SAMPL------------')
        print('------NEW SIMULTAION RUN------')
        print('------------------------------')

        print('1. Create a beam ...')
        createBeam = cb.createBeam()
        if 'ini' in self.initDistribFile:
            self.initDistrib = createBeam.useASTRAFile(self.initDistribFile)
        else:
            print('Creating default beam.')
            xOffset = caget('VM-EBT-INJ-DIA-DUMMY-01:DDOUBLE8')
            yOffset = caget('VM-EBT-INJ-DIA-DUMMY-01:DDOUBLE9')
            self.initDistrib = createBeam.guassian(x=xOffset, y=yOffset)

        print('2. Create a beamline ...')
        self.loadPathway()
        lineCreator = cbl.createBeamline(V_MAG_Ctrl=self.V_MAG_Ctrl,
                                         C_S01_MAG_Ctrl=self.C_S01_MAG_Ctrl,
                                         C_S02_MAG_Ctrl=self.C_S02_MAG_Ctrl,
                                         C2V_MAG_Ctrl=self.C2V_MAG_Ctrl,
                                         V_RF_Ctrl=self.V_RF_Ctrl,
                                         C_RF_Ctrl=self.C_RF_Ctrl,
                                         L01_RF_Ctrl=self.L01_RF_Ctrl)
        beamLine = lineCreator.create(self.pathway, self.pathway.elements)
        # Run simulation
        for element in self.pathway.elements.keys():
            if element == self.startElement:
                startName = element
            if element == self.stopElement:
                stopName = element

        startIndex = [i for i, x in enumerate(beamLine.componentlist) if x.name == startName]
        stopIndex = [i for i, x in enumerate(beamLine.componentlist) if x.name == stopName]
        print('3. Running SAMPL simulation from ' +
              startName + ' to ' + stopName + '.')
        beamLine.TrackMatlab([startIndex[0], stopIndex[0]], self.initDistrib)

        # SAMPL to EPICS (look how short it is it is!!!!!!)
        # this take the most time to complete
        print('4. Writing data to EPICS ...')
        """CHANGED FOR NEW CAMER PVs"""
        for i in beamLine.componentlist:
            if beamLine.componentlist.index(i) >= startIndex[0] and beamLine.componentlist.index(i) <= stopIndex[0]:
                if 'SCR' in i.name or 'YAG' in i.name:
                    if 'CLu' in i.name:
                        caput('VM-' + self.elements[i.name]['camPV'] +
                              ':ANA:X_RBV', i.x)
                        caput('VM-' + self.elements[i.name]['camPV'] +
                              ':ANA:Y_RBV', i.y)
                        caput('VM-' + self.elements[i.name]['camPV'] +
                              ':ANA:SigmaX_RBV', i.xSigma)
                        caput('VM-' + self.elements[i.name]['camPV'] +
                              ':ANA:SigmaY_RBV', i.ySigma)
                    else:
                        caput('VM-' + self.elements[i.name]['camPV'] +
                              ':X', i.x)
                        caput('VM-' + self.elements[i.name]['camPV'] +
                              ':Y', i.y)
                        caput('VM-' + self.elements[i.name]['camPV'] +
                              ':SigmaX', i.xSigma)
                        caput('VM-' + self.elements[i.name]['camPV'] +
                              ':SigmaY', i.ySigma)
                    print '    Written data for ', self.elements[i.name]['Online_Model_Name']
                if 'BPM'in i.name:
                    caput('VM-' + self.elements[i.name]['pv'] + ':X', i.x)
                    caput('VM-' + self.elements[i.name]['pv'] + ':Y', i.y)
                    print '    Written data for ', self.elements[i.name]['omName']
                    print 'x Value:', str(i.x)
