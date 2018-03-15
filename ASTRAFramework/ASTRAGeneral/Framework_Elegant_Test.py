import sys
from Framework import *
from FrameworkHelperFunctions import *
sys.path.append('../../ElegantWriter')
from elegantWriter_objects import *

type_conversion_rules = {'dipole': 'csrcsbend', 'quadrupole': 'kquad', 'beam_position_monitor': 'moni', 'screen': 'watch', 'aperture': 'rcol',
                         'collimator': 'ecol', 'monitor': 'moni', 'solenoid': 'mapsolenoid', 'wall_current_monitor': 'charge', 'cavity': 'rfcw',
                         'rf_deflecting_cavity': 'rfdf'}
keyword_conversion_rules = {'length': 'l','entrance_edge_angle': 'e1', 'exit_edge_angle': 'e2', 'edge_field_integral': 'fint', 'horizontal_size': 'x_max', 'vertical_size': 'y_max'}

def addElement(name):
    elem = framework.getElement(name)
    if getParameter(elem,'entrance_edge_angle') == 'angle':
        framework.modifyElement(name,'entrance_edge_angle', getParameter(elem,'angle') )
    if getParameter(elem,'exit_edge_angle') == 'angle':
        framework.modifyElement(name,'exit_edge_angle', getParameter(elem,'angle') )
    element_dict = dict(framework.getElement(name))
    if 'type' in element_dict:
        if element_dict['type'] in type_conversion_rules:
            element_dict['type'] = type_conversion_rules[element_dict['type']]
        element_dict = dict((keyword_conversion_rules[k] if k in keyword_conversion_rules else k, v) for k, v in element_dict.iteritems())
        try:
            testlattice.addElement(name=name, **element_dict)
        # print getattr(testlattice, name).properties
        except:
            print name
    else:
        print name

testlattice = elegantLattice()

globalsettings = testlattice.addCommand(name='global',type='global_settings',log_file="elegant.log",error_log_file="elegant.err")
runsetup = testlattice.addCommand(name='runsetup',type='run_setup',lattice="doublering.lte",use_beamline="doublering",p_central_mev=700,centroid='%s.cen')

framework = Framework('longitudinal_best', overwrite=False)
framework.loadSettings('clara400_v12.def')
# print framework.elements
lattice = framework.createDrifts(list(framework.elements.iteritems()))
map(addElement, lattice.keys())
testlattice.addLine(name='CLA', line=lattice.keys())
print testlattice['CLA-L01-APER']
print testlattice['CLA']
testlattice.writeLatticeFile('CLA.lte', 'CLA')
exit()
