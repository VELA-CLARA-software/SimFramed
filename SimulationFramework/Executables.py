import socket
import os

class Executables(object):

    def __init__(self, master_lattice):
        super(Executables, self).__init__()
        self.osname = os.name
        self.hostname = socket.gethostname()
        if master_lattice is None:
            master_lattice_location = (os.path.relpath(os.path.dirname(os.path.abspath(__file__)) + '/../MasterLattice/')+'/').replace('\\','/')
        else:
            master_lattice_location = master_lattice
        self.define_generator_command()
        self.define_astra_command()
        self.define_elegant_command()
        self.define_csrtrack_command()
        self.define_gpt_command()

    def __getitem__(self, item):
        return getattr(self, item)

    def getNCPU(self, ncpu, scaling):
        if scaling is not None and ncpu == 1:
            return 3*scaling
        else:
            return ncpu

    def define_generator_command(self):
        if not self.osname == 'nt':
            if 'apclara1' in self.hostname:
                self.generator =  ['/opt/ASTRA/generator.sh']
            elif 'apclara2' in self.hostname:
                self.generator =  ['/opt/ASTRA/generator.sh']
            elif 'apclara3' in self.hostname:
                self.generator =  ['/opt/ASTRA/generator.sh']
        else:
            self.generator =  [master_lattice_location+'Codes/generator']

    def define_astra_command(self, ncpu=1, scaling=None):
        ncpu = self.getNCPU(ncpu, scaling)
        if not self.osname == 'nt':
            if 'apclara1' in self.hostname:
                self.astra = ['mpiexec','-np',str(ncpu),'/opt/ASTRA/astra_MPICH2.sh']
            elif 'apclara2' in self.hostname:
                self.astra =  ['/usr/lib64/compat-openmpi16/bin/mpiexec','-np',str(ncpu),'/opt/ASTRA/astra_r62_Linux_x86_64_OpenMPI_sld6']
            elif 'apclara3' in self.hostname:
                self.astra =  ['/usr/lib64/compat-openmpi16/bin/mpiexec','-np',str(ncpu),'/opt/ASTRA/astra_r62_Linux_x86_64_OpenMPI_sld6']
        else:
            self.astra =  [master_lattice_location+'Codes/astra']

    def define_elegant_command(self, ncpu=1, scaling=None):
        ncpu = self.getNCPU(ncpu, scaling)
        if ncpu > 1:
            if not self.osname == 'nt':
                if 'apclara1' in self.hostname:
                    self.elegant = ['/opt/MPICH2-3.2/bin/mpiexec','-np',str(ncpu),'Pelegant']
                elif 'apclara2' in self.hostname:
                    self.elegant = ['/usr/lib64/openmpi/bin/mpiexec','-np',str(ncpu),'Pelegant']
                elif 'apclara3' in self.hostname:
                    self.elegant = ['/usr/lib64/openmpi/bin/mpiexec','-np',str(ncpu),'Pelegant']
            else:
                self.elegant = [master_lattice_location+'Codes/elegant']
        else:
            if not self.osname == 'nt':
                if 'apclara1' in self.hostname:
                    self.elegant = ['elegant']
                elif 'apclara2' in self.hostname:
                    self.elegant = ['elegant']
                elif 'apclara3' in self.hostname:
                    self.elegant = ['elegant']
            else:
                self.elegant = [master_lattice_location+'Codes/elegant']

    def define_csrtrack_command(self, ncpu=1, scaling=None):
        ncpu = self.getNCPU(ncpu, scaling)
        if not self.osname == 'nt':
            if 'apclara1' in self.hostname:
                self.csrtrack = ['/opt/OpenMPI-1.4.3/bin/mpiexec','-n',str(3*ncpu),'/opt/CSRTrack/csrtrack_openmpi.sh']
            elif 'apclara2' in self.hostname:
                self.csrtrack = ['/usr/lib64/compat-openmpi16/bin/mpiexec','-n',str(3*ncpu),'/opt/CSRTrack/csrtrack_1.204_Linux_x86_64_OpenMPI_1.5.4_sld62']
            elif 'apclara3' in self.hostname:
                self.csrtrack = ['/usr/lib64/compat-openmpi16/bin/mpiexec','-n',str(3*ncpu),'/opt/CSRTrack/csrtrack_1.204_Linux_x86_64_OpenMPI_1.5.4_sld62']
        else:
            self.csrtrack = [master_lattice_location+'Codes/csrtrack']

    def define_gpt_command(self, ncpu=1, scaling=None):
        self.gpt = [r'C:/Program Files/General Particle Tracer/bin/gpt.exe']
