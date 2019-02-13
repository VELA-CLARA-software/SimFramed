import os, sys
import numpy as np
sys.path.append(os.path.abspath(__file__+'/../../../../'))
import SimulationFramework.Framework as Fw

Framework = Fw.Framework('phw_test_elegant_VBC')
if not os.name == 'nt':
    framework.defineASTRACommand(['mpiexec','-np',str(3*scaling),'/opt/ASTRA/astra_MPICH2.sh'])
    framework.defineGeneratorCommand(['/opt/ASTRA/generator.sh'])
    framework.defineCSRTrackCommand(['/opt/OpenMPI-1.4.3/bin/mpiexec','-n',str(3*scaling),'/opt/CSRTrack/csrtrack_openmpi.sh'])

Framework.loadSettings('Lattices/clara400_v12_v3_elegantVBC.def')
for scaling in [3,4,5,6]:
    Framework['L02'].file_block['input']['prefix'] = '../basefiles_'+str(scaling)+'/'
    oldl02grad = Framework.getElement('CLA-L02-CAV', 'field_amplitude')
    oldl02phase = Framework.getElement('CLA-L02-CAV', 'phase')
    oldl03grad = Framework.getElement('CLA-L03-CAV', 'field_amplitude')
    oldl03phase = Framework.getElement('CLA-L03-CAV', 'phase')
    oldl04grad = Framework.getElement('CLA-L04-CAV', 'field_amplitude')
    oldl04phase = Framework.getElement('CLA-L04-CAV', 'phase')
    oldl0xgrad = Framework.getElement('CLA-L4H-CAV', 'field_amplitude')
    oldl0xphase = Framework.getElement('CLA-L4H-CAV', 'phase')
    print 'oldl02grad = ', oldl02grad
    print 'oldl02phase = ', oldl02phase
    print 'oldl03grad = ', oldl03grad
    print 'oldl03phase = ', oldl03phase
    print 'oldl04grad = ', oldl04grad
    print 'oldl04phase = ', oldl04phase
    print 'oldl0xgrad = ', oldl0xgrad
    print 'oldl0xphase = ', oldl0xphase
    Framework.modifyElement('CLA-L02-CAV', 'field_amplitude', np.sqrt(2)*11.25e6 )
    Framework.modifyElement('CLA-L03-CAV', 'field_amplitude', np.sqrt(2)*11.25e6 )
    Framework.modifyElement('CLA-L04-CAV', 'field_amplitude', np.sqrt(2)*12e6 )
    Framework.modifyElement('CLA-L4H-CAV', 'field_amplitude', np.sqrt(2)*15.75e6 )
    Framework.modifyElement('CLA-L02-CAV', 'phase', oldl02phase - 5.0)
    newl02grad = Framework.getElement('CLA-L02-CAV', 'field_amplitude')
    newl02phase = Framework.getElement('CLA-L02-CAV', 'phase')
    newl03grad = Framework.getElement('CLA-L03-CAV', 'field_amplitude')
    newl03phase = Framework.getElement('CLA-L03-CAV', 'phase')
    newl04grad = Framework.getElement('CLA-L04-CAV', 'field_amplitude')
    newl04phase = Framework.getElement('CLA-L04-CAV', 'phase')
    newl0xgrad = Framework.getElement('CLA-L4H-CAV', 'field_amplitude')
    newl0xphase = Framework.getElement('CLA-L4H-CAV', 'phase')
    print 'newl02grad = ', newl02grad
    print 'newl02phase = ', newl02phase
    print 'newl03grad = ', newl03grad
    print 'newl03phase = ', newl03phase
    print 'newl04grad = ', newl04grad
    print 'newl04phase = ', newl04phase
    print 'newl0xgrad = ', newl0xgrad
    print 'newl0xphase = ', newl0xphase
    Framework.track(startfile='L02', endfile='S07')#startfile='FMS')

