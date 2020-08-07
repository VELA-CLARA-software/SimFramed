from SimulationFramework.Framework_objects import *
from SimulationFramework.Framework_elements import *

class mad8Lattice(frameworkLattice):
    def __init__(self, *args, **kwargs):
        super(mad8Lattice, self).__init__(*args, **kwargs)
        self.code = 'mad'
        self.particle_definition = self.allElementObjects[self.start].objectname
        self.bunch_charge = None
        self.trackBeam = True
        self.betax = None
        self.betay = None
        self.alphax = None
        self.alphay = None

    def writeElements(self):
        elements = self.createDrifts()
        fulltext = ''
        for element in list(elements.values()):
            if not element.subelement:
                fulltext += element.write_MAD8()
        fulltext += self.w.write_Elegant()
        fulltext += self.objectname+': Line=('
        for e, element in list(elements.items()):
            if not element.subelement:
                if len((fulltext + e).splitlines()[-1]) > 79:
                    fulltext += '&\n'
                fulltext += e+', '
        return fulltext[:-2] + ')\n'

    def write(self):
        self.lattice_file = self.global_parameters['master_subdir']+'/'+self.objectname+'.mad'
        saveFile(self.lattice_file, self.writeElements())
        self.command_file = self.global_parameters['master_subdir']+'/'+self.objectname+'.mff'
        saveFile(self.command_file, self.commandFile.write())

    def preProcess(self):
        prefix = self.file_block['input']['prefix'] if 'input' in self.file_block and 'prefix' in self.file_block['input'] else ''
        if self.trackBeam:
            self.hdf5_to_tfs(prefix)
        self.commandFile = mad8TrackFile(lattice=self, trackBeam=self.trackBeam, mad8beamfilename=self.objectname+'.tfs', sample_interval=self.sample_interval,
        betax=self.betax,
        betay=self.betay,
        alphax=self.alphax,
        alphay=self.alphay, global_parameters=self.global_parameters)

    def postProcess(self):
        if self.trackBeam:
            for s in self.screens:
                s.tfs_to_hdf5()

    def hdf5_to_tfs(self, prefix=''):
        HDF5filename = prefix+self.particle_definition+'.hdf5'
        self.global_parameters['beam'].read_HDF5_beam_file(self.global_parameters['master_subdir'] + '/' + HDF5filename)
        tfsbeamfilename = self.objectname+'.tfs'
        self.global_parameters['beam'].write_TFS_file(self.global_parameters['master_subdir'] + '/' + tfsbeamfilename, xyzoffset=self.startObject.position_start)

    def run(self):
        """Run the code with input 'filename'"""
        command = self.executables[self.code] + [self.objectname+'.mff']
        # print ('run command = ', command)
        with open(os.path.relpath(self.global_parameters['master_subdir']+'/'+self.objectname+'.log', '.'), "w") as f:
            subprocess.call(command, stdout=f, cwd=self.global_parameters['master_subdir'])

class mad8CommandFile(frameworkObject):
    def __init__(self, lattice='', *args, **kwargs):
        super(mad8CommandFile, self).__init__()
        self.commandObjects = OrderedDict()
        self.lattice_filename = lattice.objectname+'.mff'

    def addCommand(self, name=None, **kwargs):
        if name == None:
            if not 'name' in kwargs:
                if not 'type' in kwargs:
                    raise NameError('Command does not have a name')
                else:
                    name = kwargs['type']
            else:
                name = kwargs['name']
        command = frameworkCommand(name, global_parameters=self.global_parameters, **kwargs)
        self.commandObjects[name] = command
        return command

    def write(self):
        output = ''
        for c in list(self.commandObjects.values()):
            output += c.write_mad8()
        return output

class mad8TrackFile(mad8CommandFile):
    def __init__(self, lattice='', trackBeam=True, mad8beamfilename='', betax=None, betay=None, alphax=None, alphay=None, *args, **kwargs):
        super(mad8TrackFile, self).__init__(lattice, *args, **kwargs)
        self.mad8beamfilename = mad8beamfilename
        self.sample_interval = kwargs['sample_interval'] if 'sample_interval' in kwargs else 1
        self.trackBeam = trackBeam
        self.betax = betax if betax is not None else self.global_parameters['beam'].beta_x
        self.betay = betay if betay is not None else self.global_parameters['beam'].beta_y
        self.alphax = alphax if alphax is not None else self.global_parameters['beam'].alpha_x
        self.alphay = alphay if alphay is not None else self.global_parameters['beam'].alpha_y
        self.addCommand(type='global_settings')#, mpi_io_read_buffer_size=262144, mpi_io_write_buffer_size=262144)
        self.addCommand(type='run_setup',lattice=self.lattice_filename, \
            use_beamline=lattice.objectname,p_central=np.mean(self.global_parameters['beam'].BetaGamma), \
            centroid='%s.cen',always_change_p0 = 1, \
            sigma='%s.sig', default_order=3)
        self.addCommand(type='run_control',n_steps=1, n_passes=1)
        self.addCommand(type='twiss_output',matched = 0,output_at_each_step=0,radiation_integrals=1,statistics=1,filename="%s.twi",
        beta_x  = self.betax,
        alpha_x = self.alphax,
        beta_y  = self.betay,
        alpha_y = self.alphay)
        flr = self.addCommand(type='floor_coordinates', filename="%s.flr",
        X0  = lattice.startObject['position_start'][0],
        Z0 = lattice.startObject['position_start'][2],
        theta0 = 0)
        mat = self.addCommand(type='matrix_output', SDDS_output="%s.mat",
        full_matrix_only=1, SDDS_output_order=2)
        if self.trackBeam:
            self.addCommand(type='sdds_beam', input=self.elegantbeamfilename, sample_interval=self.sample_interval)
            self.addCommand(type='track')

class mad8Optimisation(mad8CommandFile):

    def __init__(self, lattice='', variables={}, constraints={}, terms={}, settings={}, *args, **kwargs):
        super(mad8Optimisation, self).__init__(lattice, *args, **kwargs)
        for k, v in list(variables.items()):
            self.add_optimisation_variable(k, **v)

    def add_optimisation_variable(self, name, item=None, lower=None, upper=None, step=None, restrict_range=None):
        self.addCommand(name=name, type='optimization_variable', item=item, lower_limit=lower, upper_limit=upper, step_size=step, force_inside=restrict_range)

    def add_optimisation_constraint(self, name, item=None, lower=None, upper=None):
        self.addCommand(name=name, type='optimization_constraint', quantity=item, lower=lower, upper=upper)

    def add_optimisation_term(self, name, item=None, **kwargs):
        self.addCommand(name=name, type='optimization_term', term=item, **kwargs)
