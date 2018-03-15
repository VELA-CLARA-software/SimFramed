'''Assumes  you have the a OnlineModel in you path directory'''
from PyQt4.QtCore import QThread
import os
import yaml
import modifyASTRABeamline as mabl
from SimulationFramework import Framework


class Setup(QThread):
    # CONSTRUCTOR
    def __init__(self, V_MAG_Ctrl=None, C_S01_MAG_Ctrl=None,
                 C_S02_MAG_Ctrl=None, C2V_MAG_Ctrl=None, V_RF_Ctrl=None,
                 C_RF_Ctrl=None, L01_RF_Ctrl=None, messages=False):
        QThread.__init__(self)
        self.showMessages = messages
        self.VELA_MAG_Controller = V_MAG_Ctrl
        self.CLA_MAG_S01_Controller = C_S01_MAG_Ctrl
        self.CLA_MAG_S02_Controller = C_S02_MAG_Ctrl
        self.C2V_MAG_Controller = C2V_MAG_Ctrl
        self.VELA_LLRF_Controller = V_RF_Ctrl
        self.CLA_LLRF_Controller = C_RF_Ctrl
        self.CLA_L01_LLRF_Controller = L01_RF_Ctrl
        self.startElement = 'Null'
        self.stopElement = 'Null'
        self.initDistribFile = 'Null'
        self.initCharge = 0.0
        self.pathway = Framework.Framework(subdir='.', overwrite='overwrite')

    # DESTRUCTOR
    def __del__(self):
        self.wait()

    def loadPathway(self):
        stream = file(str(os.path.abspath(__file__)).split('astra')[0] +
                      "\\..\\..\\MasterLattice\\YAML\\allPathways.yaml", 'r')
        settings = yaml.load(stream)
        for path in settings['pathways']:
            hasStart = any(self.startElement in path.elements[s]['Online_Model_Name']
                           for s, value in path.elements.iteritems())
            hasStop = any(self.stopElement in path.elements[s]['Online_Model_Name']
                          for s, value in path.elements.iteritems())
            if (hasStart and hasStop):
                self.pathway.loadSettings(filename=path)
                print '    Loading pathway: ', path

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
        self.loadPathway()
        print('    Modify pathway using Virtual EPICS...')
        modifier = mabl.beamline(V_MAG_Ctrl=self.V_MAG_Ctrl,
                                 C_S01_MAG_Ctrl=self.C_S01_MAG_Ctrl,
                                 C_S02_MAG_Ctrl=self.C_S02_MAG_Ctrl,
                                 C2V_MAG_Ctrl=self.C2V_MAG_Ctrl,
                                 V_RF_Ctrl=self.V_RF_Ctrl,
                                 C_RF_Ctrl=self.C_RF_Ctrl,
                                 L01_RF_Ctrl=self.L01_RF_Ctrl)
        modifier.modfiy(self.pathway)
        print('    Creating .in files...')
        # Write .in files
        self.pathway.createInputFiles()
        # Once written copy the files to virtual Machine
        for section in self.pathway.fileSettings.keys():
            os.system('VBoxManage --nologo guestcontrol "VE-11g" copyto ' +
                      '--username "vmsim" --password "password" ' +
                      '--target-directory "/home/vmsim/Desktop/V2/ASTRA/" "' +
                      os.getcwd() + '\\temp-' + section + '"')
            os.system('VBoxManage --nologo guestcontrol "VE-11g" copyto ' +
                      '--username "vmsim" --password "password" ' +
                      '--target-directory "/home/vmsim/Desktop/V2/ASTRA/" "' +
                      os.getcwd() + '\\' + self.initDistribFile + '"')

        print('3. Running ASTRA simulation from ' +
              self.startElement + ' to ' + self.stopElement)
        # Now run Python script in Virtual Machine to run ASTRA
        if self.showMessages is True:
            os.system('VBoxManage --nologo guestcontrol "VE-11g" run ' +
                      '"usr/bin/python" --username "vmsim" --password ' +
                      '"password" -- /home/vmsim/Desktop/V2/ASTRA/' +
                      'runASTRA.py %s' % (','.join(self.pathway.fileSettings.keys())))
        else:
            os.system('VBoxManage --nologo guestcontrol "VE-11g" run ' +
                      '"usr/bin/python" --username "vmsim" --password ' +
                      '"password" -- /home/vmsim/Desktop/V2/ASTRA/runASTRA.py' +
                      ' >> OverallSimMessages.txt %s' % (','.join(inFiles)))
