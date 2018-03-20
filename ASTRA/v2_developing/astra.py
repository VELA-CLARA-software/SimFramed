'''Assumes  you have the a OnlineModel in you path directory'''
from PyQt4.QtCore import QThread
import os
import yaml
import modifyASTRABeamline as mabl
from SimulationFramework import Framework


class Setup(QThread):
    # CONSTRUCTOR
    def __init__(self, V_MAG_Ctrl=None, C_S01_MAG_Ctrl=None,
                 C_S02_MAG_Ctrl=None, C2V_MAG_Ctrl=None, LRRG_RF_Ctrl=None,
                 HRRG_RF_Ctrl=None, L01_RF_Ctrl=None,
                 messages=False, subdir='.'):
        QThread.__init__(self)
        self.showMessages = messages
        self.V_MAG_Ctrl = V_MAG_Ctrl
        self.C_S01_MAG_Ctrl = C_S01_MAG_Ctrl
        self.C_S02_MAG_Ctrl = C_S02_MAG_Ctrl
        self.C2V_MAG_Ctrl = C2V_MAG_Ctrl
        self.LRRG_RF_Ctrl = LRRG_RF_Ctrl
        self.HRRG_RF_Ctrl = HRRG_RF_Ctrl
        self.L01_RF_Ctrl = L01_RF_Ctrl
        self.startElement = 'Null'
        self.stopElement = 'Null'
        self.initDistribFile = 'Null'
        self.initCharge = 0.0
        self.pathway = Framework.Framework(subdir=subdir, overwrite='overwrite')


    # DESTRUCTOR
    def __del__(self):
        self.wait()

    def loadPathway(self):
        ans = False
        stream = file(str(os.path.abspath(__file__)).split('astra')[0] +
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

    def go(self, startElement, stopElement, initDistribFile, charge=0.25):
            self.startElement = startElement
            self.stopElement = stopElement
            self.initDistribFile = initDistribFile
            # charge in nC
            self.initCharge = charge
            # Run in Thread
            self.start()
            # Don't run in thread
            # self.run()

    def run(self):
        print('------------------------------')
        print('-------------ASTRA------------')
        print('------NEW SIMULTAION RUN------')
        print('------------------------------')

        print('1. Define Beam...')
        print('    Using file: ' + self.initDistribFile)

        print('2. Create a Beamline...')
        print('    Find aproriate pathway...')
        sucessful = self.loadPathway()
        self.pathway.astra.setInitialDistribution(filename=self.initDistribFile)
        if sucessful is False:
            return
        print('    Modify pathway using Virtual EPICS...')
        modifier = mabl.beamline(V_MAG_Ctrl=self.V_MAG_Ctrl,
                                 C_S01_MAG_Ctrl=self.C_S01_MAG_Ctrl,
                                 C_S02_MAG_Ctrl=self.C_S02_MAG_Ctrl,
                                 C2V_MAG_Ctrl=self.C2V_MAG_Ctrl,
                                 LRRG_RF_Ctrl=self.LRRG_RF_Ctrl,
                                 HRRG_RF_Ctrl=self.HRRG_RF_Ctrl,
                                 L01_RF_Ctrl=self.L01_RF_Ctrl)
        modifier.modfiy(self.pathway, self.startElement, self.stopElement)
        print('    Crop pathway for ASTRA run...')
        startIndex = self.pathway.elementIndex(self.startElement)
        stopIndex = self.pathway.elementIndex(self.stopElement)

        for section in self.pathway.fileSettings.keys():
            startOfSection = self.pathway.fileSettings[section]['output']['start_element']
            stopOfSection = self.pathway.fileSettings[section]['output']['end_element']
            sectionStartIndex = self.pathway.elementIndex(startOfSection)
            sectionStopIndex = self.pathway.elementIndex(stopOfSection)
            if startIndex > sectionStartIndex and startIndex > sectionStopIndex:
                print '      Deleting section', section
                del self.pathway.fileSettings[section]
            if startIndex > sectionStartIndex and startIndex < sectionStopIndex:
                print '      Make start element the start of section', section
                self.pathway.fileSettings[section]['output']['start_element'] = self.startElement
            if stopIndex > sectionStartIndex and stopIndex < sectionStopIndex:
                print '      Make stop element the stop of section', section
                self.pathway.fileSettings[section]['output']['end_element'] = self.stopElement
            if stopIndex < sectionStartIndex and stopIndex < sectionStopIndex:
                print '      Deleting section', section
                del self.pathway.fileSettings[section]
        print self.pathway.fileSettings.keys()
        print('    Creating .in files...')

        # Write .in files
        self.pathway.createInputFiles()

        print('3. Running ASTRA simulation from ' +
              self.startElement + ' to ' + self.stopElement)

        #self.pathway.runInputFiles()
