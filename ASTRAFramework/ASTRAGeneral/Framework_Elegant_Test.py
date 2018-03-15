import sys
from Framework import *
from FrameworkHelperFunctions import *
sys.path.append('../../ElegantWriter')
from elegantWriter_objects import *

type_conversion_rules = {'dipole': 'csrcsbend', 'quadrupole': 'kquad', 'beam_position_monitor': 'moni', 'screen': 'watch', 'aperture': 'ecol',
                         'collimator': 'ecol', 'monitor': 'moni', 'solenoid': 'mapsolenoid', 'wall_current_monitor': 'charge', 'cavity': 'rfcw',
                         'rf_deflecting_cavity': 'rfdf'}
keyword_conversion_rules = {'length': 'l','entrance_edge_angle': 'e1', 'exit_edge_angle': 'e2', 'edge_field_integral': 'fint'}

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

testlattice = elegantLattice()

globalsettings = testlattice.addCommand(name='global',type='global_settings',log_file="elegant.log",error_log_file="elegant.err")
runsetup = testlattice.addCommand(name='runsetup',type='run_setup',lattice="doublering.lte",use_beamline="doublering",p_central_mev=700,centroid='%s.cen')

framework = Framework('longitudinal_best', overwrite=False)
framework.loadSettings('clara400_v12.def')
print np.array(framework.createDrifts(list(framework.elements.iteritems())))[:40,0]
map(addElement, framework.elements.keys())
# addElement('CLA-VBC-MAG-DIP-01')
exit()

testlattice.addElement(name='Q1',type='kquad',l=0.5,k1=2.3)
testlattice.addElement(name='D1',type='drift',l=0.5)
testlattice.addElement(name='BEND1',type='ksbend',l=0.5,angle=0.36,E1=0.01,E2=0.022)
#
# # print D1*2
latt1 = testlattice.addLine(name='latt1',line=[2*testlattice.D1,testlattice.Q1*2])
latt3 = testlattice.addLine(name='latt3',line=['2*latt1','-BEND1','latt1','D1'])
latt2 = testlattice.addLine(name='latt2',line=['2*latt1','latt3','D1'])

# print testlattice.Q1
# print testlattice['elements']['Q1']
print testlattice.latt1.elements
# testlattice.writeCommandFile('test.ele')
# testlattice.writeLatticeFile('test.lte')
