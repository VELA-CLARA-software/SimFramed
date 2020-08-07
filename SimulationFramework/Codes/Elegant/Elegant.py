import os
from ruamel import yaml
from SimulationFramework.Framework_objects import *
from SimulationFramework.Framework_elements import *

class elegantLattice(frameworkLattice):
    def __init__(self, *args, **kwargs):
        super(elegantLattice, self).__init__(*args, **kwargs)
        self.code = 'elegant'
        self.particle_definition = self.allElementObjects[self.start].objectname
        self.bunch_charge = None
        self.q = charge(name='START', type='charge', global_parameters=self.global_parameters,**{'total': 250e-12})
        self.trackBeam = True
        self.betax = None
        self.betay = None
        self.alphax = None
        self.alphay = None

    def endScreen(self, **kwargs):
        return screen(name='end', type='screen', position_start=self.endObject.position_start, position_end=self.endObject.position_start, global_rotation=self.endObject.global_rotation, global_parameters=self.global_parameters, **kwargs)

    def createDrifts(self):
        """Insert drifts into a sequence of 'elements'"""
        positions = []
        originalelements = OrderedDict()
        elementno = 0
        newelements = OrderedDict()
        for name in list(self.elements.keys()):
            if not self.elements[name].subelement:
                originalelements[name] = self.elements[name]
                pos = np.array(self.allElementObjects[name].position_start)
                positions.append(pos)
                positions.append(self.allElementObjects[name].position_end)
        positions = positions[1:]
        positions.append(positions[-1])
        driftdata = list(zip(iter(list(originalelements.items())), list(chunks(positions, 2))))

        lscbins = self.lsc_bins if self.lscDrifts is True else 0
        csr = 1 if self.csrDrifts is True else 0
        lsc = 1 if self.lscDrifts is True else 0
        drifttype = csrdrift if self.csrDrifts else lscdrift

        for e, d in driftdata:
            if e[1]['objecttype'] == 'screen' and e[1]['length'] > 0:
                name = e[0]+'-drift-01'
                newdrift = drifttype(name, global_parameters=self.global_parameters, **{'length': e[1]['length']/2,
                 'csr_enable': csr,
                 'lsc_enable': lsc,
                 'use_stupakov': 1,
                 'csrdz': 0.01,
                 'lsc_bins': lscbins,
                 'lsc_high_frequency_cutoff_start': self.lsc_high_frequency_cutoff_start,
                 'lsc_high_frequency_cutoff_end': self.lsc_high_frequency_cutoff_end,
                 'lsc_low_frequency_cutoff_start': self.lsc_low_frequency_cutoff_start,
                 'lsc_low_frequency_cutoff_end': self.lsc_low_frequency_cutoff_end,
                })
                newelements[name] = newdrift
                newelements[e[0]] = e[1]
                name = e[0]+'-drift-02'
                newdrift = drifttype(name, global_parameters=self.global_parameters, **{'length': e[1]['length']/2,
                 'csr_enable': csr,
                 'lsc_enable': lsc,
                 'use_stupakov': 1,
                 'csrdz': 0.01,
                 'lsc_bins': lscbins,
                 'lsc_high_frequency_cutoff_start': self.lsc_high_frequency_cutoff_start,
                 'lsc_high_frequency_cutoff_end': self.lsc_high_frequency_cutoff_end,
                 'lsc_low_frequency_cutoff_start': self.lsc_low_frequency_cutoff_start,
                 'lsc_low_frequency_cutoff_end': self.lsc_low_frequency_cutoff_end,
                })
                newelements[name] = newdrift
            else:
                newelements[e[0]] = e[1]
            if len(d) > 1:
                x1, y1, z1 = d[0]
                x2, y2, z2 = d[1]
                try:
                    length = np.sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)
                except Exception as exc:
                    print('Element with error = ', e[0])
                    print(d)
                    raise exc
                if length > 0:
                    elementno += 1
                    name = 'drift'+str(elementno)
                    newdrift = drifttype(name, global_parameters=self.global_parameters, **{'length': length,
                     'position_start': list(d[0]),
                     'position_end': list(d[1]),
                     'csr_enable': csr,
                     'lsc_enable': lsc,
                     'use_stupakov': 1,
                     'csrdz': 0.01,
                     'lsc_bins': lscbins,
                     'lsc_high_frequency_cutoff_start': self.lsc_high_frequency_cutoff_start,
                     'lsc_high_frequency_cutoff_end': self.lsc_high_frequency_cutoff_end,
                     'lsc_low_frequency_cutoff_start': self.lsc_low_frequency_cutoff_start,
                     'lsc_low_frequency_cutoff_end': self.lsc_low_frequency_cutoff_end,
                    })
                    newelements[name] = newdrift
                elif length < 0:
                    raise Exception('Lattice has negative drifts!', name, length)
                    exit()
        return newelements

    def getSValues(self):
        elems = self.createDrifts()
        s = [0]
        for e in list(elems.values()):
            s.append(s[-1]+e.length)
        return s[1:]

    def getNames(self):
        elems = self.createDrifts()
        return [e.objectname for e in list(elems.values())]

    def getSNames(self):
        s = self.getSValues()
        names = self.getNames()
        return list(zip(names, s))

    def findS(self, elem):
        if elem in self.allElements:
            sNames = self.getSNames()
            return [a for a in sNames if a[0] == elem]

    def writeElements(self):
        self.w = self.endScreen(output_filename=self.end+'.SDDS')
        elements = self.createDrifts()
        fulltext = ''
        fulltext += self.q.write_Elegant()
        for element in list(elements.values()):
            # print(element.write_Elegant())
            if not element.subelement:
                fulltext += element.write_Elegant()
        fulltext += self.w.write_Elegant()
        fulltext += self.objectname+': Line=(START, '
        for e, element in list(elements.items()):
            if not element.subelement:
                if len((fulltext + e).splitlines()[-1]) > 60:
                    fulltext += '&\n'
                fulltext += e+', '
        return fulltext[:-2] + ', END )\n'

    def write(self):
        self.lattice_file = self.global_parameters['master_subdir']+'/'+self.objectname+'.lte'
        saveFile(self.lattice_file, self.writeElements())
        self.command_file = self.global_parameters['master_subdir']+'/'+self.objectname+'.ele'
        saveFile(self.command_file, self.commandFile.write())

    def preProcess(self):
        prefix = self.file_block['input']['prefix'] if 'input' in self.file_block and 'prefix' in self.file_block['input'] else ''
        if self.trackBeam:
            self.hdf5_to_sdds(prefix)
        self.commandFile = elegantTrackFile(lattice=self, trackBeam=self.trackBeam, elegantbeamfilename=self.objectname+'.sdds', sample_interval=self.sample_interval,
        betax=self.betax,
        betay=self.betay,
        alphax=self.alphax,
        alphay=self.alphay,
        global_parameters=self.global_parameters)

    def postProcess(self):
        if self.trackBeam:
            for s in self.screens:
                s.sdds_to_hdf5()
            self.w.sdds_to_hdf5()

    def hdf5_to_sdds(self, prefix=''):
        HDF5filename = prefix+self.particle_definition+'.hdf5'
        self.global_parameters['beam'].read_HDF5_beam_file(self.global_parameters['master_subdir'] + '/' + HDF5filename)
        if self.bunch_charge is not None:
            self.q = charge(name='START', type='charge', global_parameters=self.global_parameters,**{'total': abs(self.bunch_charge)})
        else:
            self.q = charge(name='START', type='charge', global_parameters=self.global_parameters,**{'total': abs(self.global_parameters['beam'].charge)})
        # print('mean cpz = ', np.mean(self.global_parameters['beam'].cpz), ' prefix = ', prefix)
        sddsbeamfilename = self.objectname+'.sdds'
        self.global_parameters['beam'].write_SDDS_file(self.global_parameters['master_subdir'] + '/' + sddsbeamfilename, xyzoffset=self.startObject.position_start)

    def run(self):
        """Run the code with input 'filename'"""
        if not os.name == 'nt':
            command = self.executables[self.code] + ['-rpnDefns='+os.path.relpath(self.global_parameters['master_lattice_location'],self.global_parameters['master_subdir'])+'/Codes/defns.rpn'] + [self.objectname+'.ele']
        else:
            command = self.executables[self.code] + ['-rpnDefns='+os.path.relpath(self.global_parameters['master_lattice_location'],self.global_parameters['master_subdir'])+'/Codes/defns.rpn'] + [self.objectname+'.ele']
            command = [c.replace('/','\\') for c in command]
        # print ('run command = ', command)
        with open(os.path.relpath(self.global_parameters['master_subdir']+'/'+self.objectname+'.log', '.'), "w") as f:
            subprocess.call(command, stdout=f, cwd=self.global_parameters['master_subdir'])

class elegantCommandFile(object):
    def __init__(self, lattice='', *args, **kwargs):
        super(elegantCommandFile, self).__init__()
        self.global_parameters = kwargs['global_parameters']
        self.commandObjects = OrderedDict()
        self.lattice_filename = lattice.objectname+'.lte'

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
            output += c.write_Elegant()
        return output

class elegantTrackFile(elegantCommandFile):
    def __init__(self, lattice='', trackBeam=True, elegantbeamfilename='', betax=None, betay=None, alphax=None, alphay=None, etax=None, etaxp=None, *args, **kwargs):
        super(elegantTrackFile, self).__init__(lattice, *args, **kwargs)
        self.elegantbeamfilename = elegantbeamfilename
        self.sample_interval = kwargs['sample_interval'] if 'sample_interval' in kwargs else 1
        self.trackBeam = trackBeam
        self.betax = betax if betax is not None else self.global_parameters['beam'].beta_x
        self.betay = betay if betay is not None else self.global_parameters['beam'].beta_y
        self.alphax = alphax if alphax is not None else self.global_parameters['beam'].alpha_x
        self.alphay = alphay if alphay is not None else self.global_parameters['beam'].alpha_y
        self.etax = etax if etax is not None else self.global_parameters['beam'].eta_x
        self.etaxp = etaxp if etaxp is not None else self.global_parameters['beam'].eta_xp
        self.addCommand(type='global_settings', mpi_io_read_buffer_size=16777216, mpi_io_write_buffer_size=16777216, inhibit_fsync=1)
        self.addCommand(type='run_setup',lattice=self.lattice_filename, \
            use_beamline=lattice.objectname,p_central=np.mean(self.global_parameters['beam'].BetaGamma), \
            centroid='%s.cen',always_change_p0 = 1, \
            sigma='%s.sig', default_order=3)
        self.addCommand(type='run_control',n_steps=1, n_passes=1)
        self.addCommand(type='twiss_output',matched = 0,output_at_each_step=0,radiation_integrals=1,statistics=1,filename="%s.twi",
        beta_x  = self.betax,
        alpha_x = self.alphax,
        beta_y  = self.betay,
        alpha_y = self.alphay,
        eta_x = self.etax,
        etap_x = self.etaxp)
        flr = self.addCommand(type='floor_coordinates', filename="%s.flr",
        X0  = lattice.startObject['position_start'][0],
        Z0 = lattice.startObject['position_start'][2],
        theta0 = 0)
        mat = self.addCommand(type='matrix_output', SDDS_output="%s.mat",
        full_matrix_only=1, SDDS_output_order=2)
        if self.trackBeam:
            self.addCommand(type='sdds_beam', input=self.elegantbeamfilename, sample_interval=self.sample_interval)
            self.addCommand(type='track')

class elegantOptimisation(elegantCommandFile):

    def __init__(self, lattice='', variables={}, constraints={}, terms={}, settings={}, *args, **kwargs):
        super(elegantOptimisation, self).__init__(lattice, *args, **kwargs)
        for k, v in list(variables.items()):
            self.add_optimisation_variable(k, **v)

    def add_optimisation_variable(self, name, item=None, lower=None, upper=None, step=None, restrict_range=None):
        self.addCommand(name=name, type='optimization_variable', item=item, lower_limit=lower, upper_limit=upper, step_size=step, force_inside=restrict_range)

    def add_optimisation_constraint(self, name, item=None, lower=None, upper=None):
        self.addCommand(name=name, type='optimization_constraint', quantity=item, lower=lower, upper=upper)

    def add_optimisation_term(self, name, item=None, **kwargs):
        self.addCommand(name=name, type='optimization_term', term=item, **kwargs)
