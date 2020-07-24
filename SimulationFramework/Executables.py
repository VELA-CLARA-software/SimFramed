import socket
import os

class Executables(object):

    def __init__(self, master_lattice):
        super(Executables, self).__init__()
        self.osname = os.name
        self.hostname = socket.gethostname()
        if master_lattice is None:
            self.master_lattice_location = (os.path.relpath(os.path.dirname(os.path.abspath(__file__)) + '/../MasterLattice/')+'/').replace('\\','/')
        else:
            self.master_lattice_location = master_lattice
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

    def define_generator_command(self, location=None):
        if location is not None:
            if isinstance(location,str):
                self.generator = [location]
            elif isinstance(location,list):
                self.generator = location
        elif not self.osname == 'nt':
            if 'apclara1' in self.hostname:
                self.generator =  ['/opt/ASTRA/generator.sh']
            elif 'apclara2' in self.hostname:
                self.generator =  ['/opt/ASTRA/generator.sh']
            elif 'apclara3' in self.hostname:
                self.generator =  ['/opt/ASTRA/generator.sh']
        else:
            self.generator =  [self.master_lattice_location+'Codes/generator']

    def define_astra_command(self, location=None, ncpu=1, scaling=None):
        ncpu = self.getNCPU(ncpu, scaling)
        if location is not None:
            if isinstance(location,str):
                self.astra = [location]
            elif isinstance(location,list):
                self.astra = location
        elif not self.osname == 'nt':
            if 'apclara1' in self.hostname:
                self.astra = ['mpiexec','-np',str(ncpu),'/opt/ASTRA/astra_MPICH2.sh']
            elif 'apclara2' in self.hostname:
                self.astra =  ['salloc','-n',str(ncpu),'/usr/lib64/compat-openmpi16/bin/mpirun','/opt/ASTRA/astra_r62_Linux_x86_64_OpenMPI_sld6']
            elif 'apclara3' in self.hostname:
                self.astra =  ['salloc','-n',str(ncpu),'/usr/lib64/compat-openmpi16/bin/mpirun','/opt/ASTRA/astra_r62_Linux_x86_64_OpenMPI_sld6']
        else:
            self.astra =  [self.master_lattice_location+'Codes/astra']

    def define_elegant_command(self, location=None, ncpu=1, scaling=None):
        ncpu = self.getNCPU(ncpu, scaling)
        if location is not None:
            if isinstance(location,str):
                self.elegant = [location]
            elif isinstance(location,list):
                self.elegant = location
        elif ncpu > 1:
            if not self.osname == 'nt':
                if 'apclara1' in self.hostname:
                    self.elegant = ['/opt/MPICH2-3.2/bin/mpiexec','-np',str(ncpu),'Pelegant']
                elif 'apclara2' in self.hostname:
                    self.elegant = ['srun','--mpi=pmi2','-n',str(ncpu),'Pelegant']
                elif 'apclara3' in self.hostname:
                    self.elegant = ['srun','--mpi=pmi2','-n',str(ncpu),'Pelegant']
            else:
                self.elegant = ['mpiexec','-np',str(ncpu),'Pelegant']
        else:
            if not self.osname == 'nt':
                if 'apclara1' in self.hostname:
                    self.elegant = ['elegant']
                elif 'apclara2' in self.hostname:
                    self.elegant = ['srun','elegant']
                elif 'apclara3' in self.hostname:
                    self.elegant = ['srun','elegant']
            else:
                self.elegant = [self.master_lattice_location+'Codes/elegant']

    def define_csrtrack_command(self, location=None, ncpu=1, scaling=None):
        ncpu = self.getNCPU(ncpu, scaling)
        if location is not None:
            if isinstance(location,str):
                self.csrtrack = [location]
            elif isinstance(location,list):
                self.csrtrack = location
        elif not self.osname == 'nt':
            if 'apclara1' in self.hostname:
                self.csrtrack = ['/opt/OpenMPI-1.4.3/bin/mpiexec','-n',str(ncpu),'/opt/CSRTrack/csrtrack_openmpi.sh']
            elif 'apclara2' in self.hostname:
                self.csrtrack = ['/opt/OpenMPI-1.4.5/bin/mpiexec','-n',str(ncpu),'/opt/CSRTrack/csrtrack_1.204_Linux_x86_64_OpenMPI_1.4.3']
            elif 'apclara3' in self.hostname:
                self.csrtrack = ['/opt/OpenMPI-1.4.5/bin/mpiexec','-n',str(ncpu),'/opt/CSRTrack/csrtrack_1.204_Linux_x86_64_OpenMPI_1.4.3']
        else:
            self.csrtrack = [self.master_lattice_location+'Codes/csrtrack']

    def define_gpt_command(self, location=None, ncpu=1, scaling=None):
        ncpu = self.getNCPU(ncpu, scaling)
        if location is not None:
            if isinstance(location,str):
                self.gpt = [location]
            elif isinstance(location,list):
                self.gpt = location
        elif not self.osname == 'nt':
            # print('gpt on apclara3')
            self.gpt = ['/opt/GPT3.3.6/bin/gpt', '-j',str(ncpu)]
            # print('gpt on apclara3', self.gpt)
        else:
            self.gpt = ['C:/Program Files/General Particle Tracer/bin/gpt.exe']
