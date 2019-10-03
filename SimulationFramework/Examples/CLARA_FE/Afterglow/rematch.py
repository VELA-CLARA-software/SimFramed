import os, sys
sys.path.append('../../../../')
import SimulationFramework.Framework as fw
import SimulationFramework.Modules.read_beam_file as rbf

framework = fw.Framework(None)
framework.loadSettings('./afterglow.def')
middle = framework['EBT-BA1-COFFIN-FOC'].middle
end = framework['EBT-BA1-COFFIN-FOC'].end
print('middle = ', middle)
print('end = ', end)
beam = rbf.beam()

#######  Rematch -10 file from TOMP_SC #########
# beam.read_SDDS_beam_file('EBT-BA1-COFFIN-FOC.SDDS')
# beam.rematchXPlane(beta=1, alpha=0)
# beam.rematchYPlane(beta=1, alpha=0)
# beam.write_HDF5_beam_file('EBT-BA1-COFFIN-FOC.hdf5', sourcefilename='EBT-BA1-COFFIN-FOC.SDDS', pos=middle, zoffset=end)

#######  Convert VSIM H5 file to my hdf5 beam #########

beam.read_vsim_h5_beam_file('CLARA_microlens_BeamElectrons_50.h5', charge=70e-12, interval=1)
beam.write_HDF5_beam_file('EBT-BA1-COFFIN-FOC.hdf5', sourcefilename='CLARA_microlens_BeamElectrons_50.h5', pos=middle, zoffset=end)
