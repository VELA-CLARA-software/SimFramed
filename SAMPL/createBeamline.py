import sys
import os
sys.path.append(str(os.path.dirname(os.path.abspath(__file__)))+'\\..\\SAMM3.0\\Python\\')
from SAMMcore.Components import Drift as d
from SAMMcore.Components import Dipole as D
from SAMMcore.Components import Quadrupole as Q
from SAMMcore.Components import Screen as S
from SAMMcore.Components import OrbitCorrector as C
from SAMMcore.Components import BeamPositionMonitor as BPM
from SAMMcore.SAMMlab import Beamline
from SAMMcore.SAMMlab import PhysicalUnits
import numpy as np


class createBeamline():

    def __init__(self, V_MAG_Ctrl=None, C_S01_MAG_Ctrl=None,
                 C_S02_MAG_Ctrl=None, C2V_MAG_Ctrl=None, V_RF_Ctrl=None,
                 C_RF_Ctrl=None, L01_RF_Ctrl=None):
        self.V_MAG_Ctrl = V_MAG_Ctrl
        self.C_S01_MAG_Ctrl = C_S01_MAG_Ctrl
        self.C_S02_MAG_Ctrl = C_S02_MAG_Ctrl
        self.C2V_MAG_Ctrl = C2V_MAG_Ctrl
        self.V_RF_Ctrl = V_RF_Ctrl
        self.C_RF_Ctrl = C_RF_Ctrl
        self.L01_RF_Ctrl = L01_RF_Ctrl
    def getObject(self, element,component):
        print component
        if 'EBT' in component and 'MAG' in component:
            if 'COR' in component:
                return self.V_MAG_Ctrl.getMagObjConstRef('V'+element['name']), self.V_MAG_Ctrl.getMagObjConstRef('H'+element['name'])
            else:
                return self.V_MAG_Ctrl.getMagObjConstRef(element['name'])
        elif 'S01' in component and 'MAG' in  component:
            if 'COR' in component:
                return self.C_S01_MAG_Ctrl.getMagObjConstRef('V'+element['name']), self.C_S01_MAG_Ctrl.getMagObjConstRef('H'+element['name'])
            else:
                return self.C_S01_MAG_Ctrl.getMagObjConstRef(element['name'])
        elif 'S02' in component and 'MAG' in  component:
            if 'COR' in component:
                return self.C_S02_MAG_Ctrl.getMagObjConstRef('V'+element['name']), self.C_S02_MAG_Ctrl.getMagObjConstRef('H'+element['name'])
            else:
                return self.C_S02_MAG_Ctrl.getMagObjConstRef(element['name'])
        elif 'C2V' in component and 'MAG' in component:
            if 'COR' in component:
                return self.C2V_MAG_Ctrl.getMagObjConstRef('V'+element['name']), self.C2V_MAG_Ctrl.getMagObjConstRef('H'+element['name'])
            else:
                return self.C2V_MAG_Ctrl.getMagObjConstRef(element['name'])
        elif 'L01' in  component:
            return self.L01_RF_Ctrl.getLLRFObjConstRef(element['name'])
        elif 'GUN' in  component:
            return self.V_RF_Ctrl.getLLRFObjConstRef(element['name'])
        else:
            print ("Trying to get unrecognised object.")

    def create(self, selectedGroup, elements):
        line = Beamline.Beamline(componentlist=[])
        driftCounter = 0
        for index, name in enumerate(selectedGroup):
            element = elements[name]
            component = None
            # Check element type and add accordingly
            if element['type'] == 'dipole':
                dip = self.getObject(element, name)
                angle = element['angle']
                print 'What I read from epics: ', dip.siWithPol
                field = 0.0
                if dip.siWithPol != 0.0:
                    field = np.copysign(np.polyval(dip.fieldIntegralCoefficients,
                                        abs(dip.siWithPol))/dip.magneticLength,dip.siWithPol)
                print field
                component = D.Dipole(name=name, length=element['length'],
                                    theta=angle * (np.pi / 180), field=field)
            elif element['type'] == 'quadrupole':
                quad = self.getObject(element, name)
                grad = 0.0
                if quad.siWithPol != 0.0:
                    grad = 1000*np.copysign(np.polyval(quad.fieldIntegralCoefficients,
                                            abs(quad.siWithPol))/quad.magneticLength,quad.siWithPol)
                component = Q.Quadrupole(name=name,
                                        length=element['length'],
                                        gradient=grad)
            elif element['type'] == 'kicker':
                vObj, hObj = self.getObject(element, name)
                vField = 0.0
                hField = 0.0
                if vObj.siWithPol != 0.0:
                    vField = 1000*np.copysign(np.polyval(vObj.fieldIntegralCoefficients,
                                              abs(vObj.siWithPol))/vObj.magneticLength,vObj.siWithPol)
                if hObj.siWithPol != 0.0:
                    hField = 1000*np.copysign(np.polyval(hObj.fieldIntegralCoefficients,
                                              abs(hObj.siWithPol))/hObj.magneticLength,hObj.siWithPol)
                component = C.OrbitCorrector(name=name,
                                            field=[hField, vField],
                                            length=element['length'])
            elif element['type'] == 'bpm':
                component = BPM.BeamPositionMonitor(name=name, length=element['length'])
            elif element['type'] == 'screen':
                component = S.Screen(name=name)
            elif element['type'] == 'wcm':
                component = d.Drift(name=name, length=element['length'])
            elif element['type'] == 'tdc':
                component = d.Drift(name=name, length=element['length'])
            else:
                print ('ERROR: This reader doesn\'t recognise element type of ', name)

            if index != 0:
                lastElement = elements[selectedGroup[index - 1]]
                backOfLast = lastElement['global_position'][-1]
                frontOfCurrent = element['global_position'][-1]
                if element['type'] == 'dipole':
                    frontOfCurrent = element['global_front'][-1]
                else:
                    frontOfCurrent = (frontOfCurrent -
                                      element['length'] * np.cos(element['global_rotation'][-1]*np.pi / 180))

                if frontOfCurrent < backOfLast:
                    print ('Elements ', selectedGroup[index - 1],
                           ' and ', name, ' overlap!!')
                elif frontOfCurrent > backOfLast:
                    # Add a drift before adding component
                    b = frontOfCurrent
                    a = backOfLast
                    # print 'drift: ', (b - a)/np.cos(element['global_rotation'][-1]*np.pi / 180)
                    driftCounter = driftCounter + 1
                    driftComponent = d.Drift(name='drift' + str(driftCounter),
                                             length=((b - a)
                                             /np.cos(element['global_rotation'][-1]
                                             * np.pi / 180)))
                    line.componentlist.append(driftComponent)
                else:
                    print 'No drift required', index

           # Append component
            line.componentlist.append(component)
        return line
