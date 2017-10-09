'''
This Class is the one User's are to initialize in there code.
You need to pass in your controllers, specifying which element of the
machine it is for
'''

from PyQt4.QtCore import QThread
import createBeam
import createBeamline as cbl
import yaml
import sys
import os
import numpy as np
sys.path.append(str(os.path.dirname(os.path.abspath(__file__)))+'\\..\\SAMM3.0\\Python\\')
from SAMMcore.SAMMlab import Beam
from SAMMcore.Particles import Electron
from SAMMcore.SAMMlab import PhysicalUnits
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
        self.startElement = None
        self.stopElement = None
        self.initDistrib = None
        self.initCharge = 0.0

        stream = file("VELA.yaml", 'r')
        settings = yaml.load(stream)
        self.elements = settings['elements']
        self.groups = settings['groups']

    # DESTRUCTOR
    def __del__(self):
        self.wait()

    # Shell function to run AStra simulations in a thread is need.
    # Using this 'shell' function alows me to pass in agurments
    def go(self, startElement, stopElement, initDistrib, charge=0.25):
            self.startElement = startElement
            self.stopElement = stopElement
            self.initDistrib = initDistrib
            self.initCharge = charge  # in nC
            # Run in Thread
            self.start()

    # Main functions (has to be called run if I want to use in a thread)
    def run(self):
        # create a beam
        print('Create a beam ...')
        # initDistrib = createBeam.guassian(x=0.0, y=0.0,
        #                                  sigmaX=0.001, sigmaY=0.001,
        #                                  particle='Electron', number=1000)
        self.initDistrib = Beam.Beam(species=Electron.Electron, energy=4.5 * PhysicalUnits.MeV)
        ptcle1 = [0.001, 0, 0, 0, 0, 0]
        ptcle2 = [0, 0, 0.001, 0, 0, 0]
        ptcle3 = [0, 0, 0, 0, 0, 0]
        self.initDistrib.particles = np.array([ptcle1, ptcle2, ptcle3])
        selectedGroup = None

        # Find which line to use
        print('Create a beamline ...')
        for line in self.groups:
            if (any(self.startElement in s for s in self.groups[line]) and
               any(self.stopElement in s for s in self.groups[line])):
                selectedGroup = self.groups[line]
                print 'Using group: ', line
        lineCreator = cbl.createBeamline(V_MAG_Ctrl=self.VELA_MAG_Controller,
                                         C_S01_MAG_Ctrl=self.CLA_MAG_S01_Controller,
                                         C_S02_MAG_Ctrl=self.CLA_MAG_S02_Controller,
                                         C2V_MAG_Ctrl=self.C2V_MAG_Controller,
                                         V_RF_Ctrl=self.VELA_LLRF_Controller,
                                         C_RF_Ctrl=self.CLA_LLRF_Controller,
                                         L01_RF_Ctrl=self.CLA_L01_LLRF_Controller)
        beamLine = lineCreator.create(selectedGroup, self.elements)

        # Run simulation
        print('------------------------------')
        print('------NEW SIMULTAION RUN------')
        print('------------------------------')
        startIndex = selectedGroup.index(self.startElement)
        stopIndex = selectedGroup.index(self.stopElement)
        finalDistrib = beamLine.TrackMatlab([startIndex, stopIndex], self.initDistrib)
        print(finalDistrib.particles)

        # #SAMPL to EPICS
