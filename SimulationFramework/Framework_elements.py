from operator import add
from SimulationFramework.Framework_objects import *
from SimulationFramework.FrameworkHelperFunctions import *
from SimulationFramework.FrameworkHelperFunctions import _rotation_matrix

class dipole(frameworkElement):

    def __init__(self, name=None, type='dipole', **kwargs):
        super(dipole, self).__init__(name, type, **kwargs)
        self.add_default('csr_bins', 100)
        self.add_default('deltaL', 0)
        self.add_default('csr_enable', 1)
        self.add_default('isr_enable', True)
        self.add_default('n_kicks', 10)
        self.add_default('sr_enable', True)
        self.add_default('integration_order', 4)
        self.add_default('nonlinear', 1)
        self.add_default('smoothing_half_width', 1)
        self.add_default('edge_order', 2)
        self.add_default('edge1_effects', 1)
        self.add_default('edge2_effects', 1)

    # @property
    # def middle(self):
    #     start = self.position_start
    #     length_vector = self.rotated_position([0,0, self.length / 2.0], offset=[0,0,0], theta=self.theta)
    #     starting_middle = length_vector
    #     # print(self.objectname, self.theta, self.starting_rotation, self.rotated_position(starting_middle, offset=self.starting_offset, theta=self.starting_rotation)[0])
    #     return np.array(start) + self.rotated_position(starting_middle, offset=self.starting_offset, theta=self.starting_rotation)

    @property
    def middle(self):
        sx, sy, sz = self.position_start
        angle = -self.angle
        l = self.length
        if abs(angle) > 0:
            cx = 0
            cy = 0
            cz = (l * np.tan(angle/2.0)) / angle
            vec = [cx, cy, cz]
        else:
            vec = [0,0,l/2.0]
        # print (vec)
        return np.array(self.position_start) + self.rotated_position(np.array(vec), offset=self.starting_offset, theta=self.y_rot)

    @property
    def arc_middle(self):
        sx, sy, sz = self.position_start
        angle = -self.angle
        l = self.length
        r = l / angle
        if abs(angle) > 0:
            cx = r * (np.cos(angle/2.0) - 1)
            cy = 0
            cz = r * np.sin(angle/2.0)
            vec = [cx, cy, cz]
        else:
            vec = [0,0,l/2.0]
        # print (vec)
        return np.array(self.position_start) + self.rotated_position(np.array(vec), offset=self.starting_offset, theta=self.y_rot)

    @property
    def line_middle(self):
        sx, sy, sz = self.position_start
        angle = -self.angle
        l = self.length
        r = l / angle
        if abs(angle) > 0:
            cx = 0.5 * r * (np.cos(angle) - 1)
            cy = 0
            cz = 0.5 * r * np.sin(angle)
            vec = [cx, cy, cz]
        else:
            vec = [0,0,l/2.0]
        # print (vec)
        return np.array(self.position_start) + self.rotated_position(np.array(vec), offset=self.starting_offset, theta=self.y_rot)

    @property
    def TD_middle(self):
        sx, sy, sz = self.position_start
        angle = -self.angle
        l = self.length
        r = l / angle
        if abs(angle) > 0:
            cx = 0.25 * r * (2.0 * np.cos(angle/2.0) + np.cos(angle) - 3)
            cy = 0
            cz = 0.25 * r * (2 * np.sin(angle/2.0) + np.sin(angle))
            vec = [cx, cy, cz]
        else:
            vec = [0,0,l/2.0]
        # print (vec)
        return np.array(self.position_start) + self.rotated_position(np.array(vec), offset=self.starting_offset, theta=self.y_rot)

    @property
    def end(self):
        start = self.position_start
        if abs(self.angle) > 1e-9:
            ex = -1. * (self.length * (np.cos(self.angle) - 1)) / self.angle
            ey = 0
            ez = (self.length * (np.sin(self.angle))) / self.angle
            return np.array(self.position_start) + self.rotated_position(np.array([ex, ey, ez]), offset=self.starting_offset, theta=-1*self.y_rot)
        else:
            return np.array(self.position_start) + self.rotated_position(np.array([0,0,self.length]), offset=self.starting_offset, theta=-1*self.y_rot)

    @property
    def width(self):
        if 'width' in self.objectproperties:
            return self.objectproperties['width']
        else:
            return 0.2
    @width.setter
    def width(self, w):
        self.objectproperties['width'] = w

    def __neg__(self):
        newself = copy.deepcopy(self)
        if 'exit_edge_angle' in newself.objectproperties and 'entrance_edge_angle' in newself.objectproperties:
            e1 = newself['entrance_edge_angle']
            e2 = newself['exit_edge_angle']
            newself.objectproperties['entrance_edge_angle'] = e2
            newself.objectproperties['exit_edge_angle'] = e1
        elif 'entrance_edge_angle' in newself.objectproperties:
            newself.objectproperties['exit_edge_angle'] = newself.objectproperties['entrance_edge_angle']
            del newself.objectproperties['entrance_edge_angle']
        elif 'exit_edge_angle' in newself.objectproperties:
            newself.objectproperties['entrance_edge_angle'] = newself.objectproperties['exit_edge_angle']
            del newself.objectproperties['exit_edge_angle']
        newself.objectname = '-'+newself.objectname
        return newself

    def check_value(self, estr, default=0):
        if estr in self.objectproperties:
            if isinstance(self.objectproperties[estr], str):
                return checkValue(self, self.objectproperties[estr],default)
            else:
                return self.objectproperties[estr]
        else:
            return default

    @property
    def intersect(self):
        return self.rho * np.tan(self.angle / 2.0)
    @property
    def rho(self):
        return -1*self.length/self.angle if self.length is not None and abs(self.angle) > 1e-9 else 0

    @property
    def e1(self):
        return self.check_value('entrance_edge_angle')
    @property
    def e2(self):
        return self.check_value('exit_edge_angle')

    def _write_Elegant(self):
        wholestring=''
        etype = self._convertType_Elegant(self.objecttype)
        string = self.objectname+': '+ etype
        k1 = self.k1 if self.k1 is not None else 0
        for key, value in list(merge_two_dicts({'k1': k1}, merge_two_dicts(self.objectproperties, self.objectdefaults)).items()):
            if not key is 'name' and not key is 'type' and not key is 'commandtype' and self._convertKeword_Elegant(key) in elements_Elegant[etype]:
                # if 'bins' in key or 'bins' in self._convertKeword_Elegant(key):
                    # print('BINS KEY ', key, '  ', self._convertKeword_Elegant(key))
                if 'edge_angle' in key:
                    key = self._convertKeword_Elegant(key)
                    value = getattr(self, key) if hasattr(self, key) and getattr(self, key) is not None else value
                else:
                    value = getattr(self, key) if hasattr(self, key) and getattr(self, key) is not None else value
                    key = self._convertKeword_Elegant(key)
                value = 1 if value is True else value
                value = 0 if value is False else value
                tmpstring = ', '+key+' = '+str(value)
                if len(string+tmpstring) > 76:
                    wholestring+=string+',&\n'
                    string=''
                    string+=tmpstring[2::]
                else:
                    string+= tmpstring
        wholestring+=string+';\n'
        return wholestring

    @property
    def corners(self):
        corners = [0,0,0,0]
        if hasattr(self, 'global_rotation') and self.global_rotation is not None:
            rotation =  self.global_rotation[2] if len(self.global_rotation) is 3 else self.global_rotation
        else:
            rotation = 0
        theta = self.e1+rotation
        corners[0] = np.array(list(map(add, np.transpose(self.position_start), np.dot([-self.width*self.length,0,0], _rotation_matrix(theta)))))
        corners[3] = np.array(list(map(add, np.transpose(self.position_start), np.dot([self.width*self.length,0,0], _rotation_matrix(theta)))))
        theta = self.angle-self.e2+rotation
        corners[1] = np.array(list(map(add, np.transpose(self.end), np.dot([-self.width*self.length,0,0], _rotation_matrix(theta)))))
        corners[2] = np.array(list(map(add, np.transpose(self.end), np.dot([self.width*self.length,0,0], _rotation_matrix(theta)))))
        corners = [self.rotated_position(x, offset=self.starting_offset, theta=self.starting_rotation) for x in corners]
        return corners

    def write_CSRTrack(self, n):
        z1 = self.position_start[2]
        z2 = self.position_end[2]
        return """dipole{\nposition{rho="""+str(z1)+""", psi="""+str(chop(self.theta+self.e1))+""", marker=d"""+str(n)+"""a}\nproperties{r="""+str(self.rho)+"""}\nposition{rho="""+str(z2)+""", psi="""+str(chop(self.theta+self.angle-self.e2))+""", marker=d"""+str(n)+"""b}\n}\n"""

    def write_ASTRA(self, n):
        if abs(checkValue(self, 'strength', default=0)) > 0 or abs(self.rho) > 0:
            corners = self.corners
            if self.plane is None:
                self.plane = 'horizontal'
            params = OrderedDict([
                ['D_Type', {'value': '\''+self.plane+'\'', 'default': '\'horizontal\''}],
                ['D_Gap', {'type': 'list', 'value': [self.gap, self.gap], 'default': [0.0001, 0.0001]}],
                ['D1', {'type': 'array', 'value': [corners[3][0],corners[3][2]] }],
                ['D3', {'type': 'array', 'value': [corners[2][0],corners[2][2]] }],
                ['D4', {'type': 'array', 'value': [corners[1][0],corners[1][2]] }],
                ['D2', {'type': 'array', 'value': [corners[0][0],corners[0][2]] }],
                # ['D_xoff', {'value': self.start[0] + self.dx, 'default': 0}],
                # ['D_yoff', {'value': self.start[1] + self.dy, 'default': 0}],
                # ['D_zoff', {'value': self.dz, 'default': 0}],
                # ['D_xrot', {'value': self.y_rot + self.dy_rot, 'default': 0}],
                # ['D_yrot', {'value': self.x_rot + self.dx_rot, 'default': 0}],
                ['D_zrot', {'value': self.z_rot + self.dz_rot, 'default': 0}],
                ])
            if abs(checkValue(self, 'strength', default=0)) > 0 or not abs(self.rho) > 0:
                params['D_strength'] = {'value': checkValue(self, 'strength', 0), 'default': 1e6}
            else:
                params['D_radius'] =  {'value': self.rho, 'default': 1e6}
            return self._write_ASTRA(params, n)
        else:
            return None

    def gpt_coordinates(self, position, rotation):
        x,y,z = chop(position, 1e-6)
        psi, phi, theta = rotation
        output =''
        for c in [0, 0, z]:
            output += str(c)+', '
        output += 'cos('+str(self.angle)+'), 0, -sin('+str(self.angle)+'), 0, 1 ,0'
        return output

    def write_GPT(self, Brho, ccs="wcs", *args, **kwargs):
        # field = Brho/self.rho if abs(self.rho) > 0 else 0
        field = self.angle * Brho / self.length
        if abs(field) > 0 and abs(self.rho) < 100:
            relpos, relrot = ccs.relative_position(self.position_start, self.global_rotation)
            relpos = relpos + [0, 0, self.intersect]
            coord = self.gpt_coordinates(relpos, relrot)
            new_ccs = self.gpt_ccs(ccs).name
            b1 = 1.0 / (2 * self.check_value(self.half_gap, default=0.02) * self.check_value(self.edge_field_integral, default=0.4))
            dl = 0 if self.deltaL is None else self.deltaL
            # print(self.objectname, ' - deltaL = ', dl)
            # b1 = 0.
            '''
            ccs( "wcs", 0, 0, startofdipole +  intersect1, Cos(theta), 0, -Sin(theta), 0, 1, 0, "bend1" ) ;
            sectormagnet( "wcs", "bend1", rho, field, e1, e2, 0., 100., 0 ) ;
            '''
            output = 'ccs( ' + ccs.name + ', '+ coord + ', ' + new_ccs + ');\n'
            output += 'sectormagnet( ' + ccs.name + ', '+ new_ccs +', '+str(abs(self.rho))+', '+str(abs(field))+', ' + str(abs(self.e1)) + ', ' + str(abs(self.e2)) + ', ' + str(abs(dl)) + ', ' + str(b1) + ', 0);\n'
        else:
            output = ''
        return output

    def gpt_ccs(self, ccs):
        if abs(self.angle) > 0 and abs(self.rho) < 100:
            # print('Creating new CCS')
            number = str(int(ccs._name.split('_')[1])+1) if ccs._name is not "wcs" else "1"
            name = 'ccs_' + number if ccs._name is not "wcs" else "ccs_1"
            # print('middle position = ', self.end)
            return gpt_ccs(name, self.end, self.global_rotation + np.array([0, 0, self.angle]), self.intersect)
        else:
            return ccs

class kicker(dipole):

    def __init__(self, name=None, type='kicker', **kwargs):
        super(kicker, self).__init__(name, type, **kwargs)

    @property
    def angle(self):
        hkick = self.horizontal_kick if self.horizontal_kick is not None else 0
        vkick = self.vertical_kick if self.vertical_kick is not None else 0
        return np.sqrt(hkick**2 + vkick**2)

    @property
    def z_rot(self):
        hkick = self.horizontal_kick if self.horizontal_kick is not None else 0
        vkick = self.vertical_kick if self.vertical_kick is not None else 0
        return self.global_rotation[0] + np.arctan2(vkick, hkick)

    def write_ASTRA(self, n):
        output = ''
        output = super(kicker, self).write_ASTRA(n)
        return output

    def write_GPT(self, Brho, ccs="wcs", *args, **kwargs):
        return ''

    def gpt_ccs(self, ccs):
        return ccs

class quadrupole(frameworkElement):

    def __init__(self, name=None, type='quadrupole', **kwargs):
        super(quadrupole, self).__init__(name, type, **kwargs)
        self.add_default('k1l', 0)
        self.add_default('n_kicks', 20)
        self.strength_errors = [0]


    @property
    def k1(self):
        return self.k1l / self.length
    @k1.setter
    def k1(self, k1):
        self.k1l = self.length * k1

    @property
    def dk1(self):
        return self.strength_errors[0]
    @dk1.setter
    def dk1(self, dk1):
        self.strength_errors[0] = dk1

    def write_ASTRA(self, n):
        if abs(self.k1 + self.dk1) > 0:
            return self._write_ASTRA(OrderedDict([
                ['Q_pos', {'value': self.middle[2] + self.dz, 'default': 0}],
                ['Q_xoff', {'value': self.middle[0], 'default': 0}],
                ['Q_yoff', {'value': self.middle[1] + self.dy, 'default': 0}],
                ['Q_xrot', {'value': -1*self.y_rot + self.dy_rot, 'default': 0}],
                ['Q_yrot', {'value': -1*self.x_rot + self.dx_rot, 'default': 0}],
                ['Q_zrot', {'value': -1*self.z_rot + self.dz_rot, 'default': 0}],
                ['Q_k', {'value': self.k1 + self.dk1, 'default': 0}],
                ['Q_length', {'value': self.length, 'default': 0}],
                ['Q_smooth', {'value': self.smooth, 'default': 10}],
                ['Q_bore', {'value': self.bore, 'default': 0.016}],
                ['Q_noscale', {'value': self.scale_field}],
                ['Q_mult_a', {'type': 'list', 'value': self.multipoles}],
            ]), n)
        else:
            return None

    def write_GPT(self, Brho, ccs="wcs", *args, **kwargs):
        # print(self.objectname)
        # print('self.start = ', self.position_start)
        relpos, relrot = ccs.relative_position(self.position_start, self.global_rotation)
        relpos = relpos + [0, 0, self.length/2.]
        coord = self.gpt_coordinates(relpos, relrot)
        output = str(self.objecttype) + '( ' + ccs.name + ', "z", '+ str(relpos[2]) +', '+str(self.length)+', '+str(-Brho*self.k1)+');\n'
        # coord = self.gpt_coordinates(self.middle, self.global_rotation)
        # output = str(self.objecttype) + '( "wcs", ' + coord + ', '+str(self.length)+', '+str(-Brho*self.k1)+');\n'
        return output

class cavity(frameworkElement):

    def __init__(self, name=None, type='cavity', **kwargs):
        super(cavity, self).__init__(name, type, **kwargs)
        self.add_default('tcolumn', '"t"')
        self.add_default('wzcolumn', '"W"')
        self.add_default('wxcolumn', '"W"')
        self.add_default('wycolumn', '"W"')
        self.add_default('wcolumn', '"Ez"')
        self.add_default('change_p0', 1)
        self.add_default('n_kicks', self.n_cells)
        # self.add_default('method', '"non-adaptive runge-kutta"')
        self.add_default('end1_focus', 1)
        self.add_default('end2_focus', 1)
        self.add_default('body_focus_model', "SRS")
        self.add_default('lsc_bins', 100)
        self.add_default('current_bins', 0)
        self.add_default('interpolate_current_bins', 1)
        self.add_default('smooth_current_bins', 1)

    def update_field_definition(self):
        if hasattr(self, 'field_definition') and self.field_definition is not None:
            self.field_definition = '"' + expand_substitution(self, '\''+self.field_definition+'\'').strip('\'"')+'"'
        if hasattr(self, 'field_definition_sdds') and self.field_definition_sdds is not None:
            self.field_definition_sdds = '"' + expand_substitution(self, '\''+self.field_definition_sdds+'\'').strip('\'"')+'"'
        if hasattr(self, 'field_definition_gdf') and self.field_definition_gdf is not None:
            self.field_definition_gdf = '"' + expand_substitution(self, '\''+self.field_definition_gdf+'\'').strip('\'"')+'"'
        if hasattr(self, 'longitudinal_wakefield_sdds') and self.longitudinal_wakefield_sdds is not None:
            self.longitudinal_wakefield_sdds = '"' + expand_substitution(self, '\''+self.longitudinal_wakefield_sdds+'\'').strip('\'"')+'"'
        if hasattr(self, 'transverse_wakefield_sdds') and self.transverse_wakefield_sdds is not None:
            self.transverse_wakefield_sdds = '"' + expand_substitution(self, '\''+self.transverse_wakefield_sdds+'\'').strip('\'"')+'"'

    @property
    def cells(self):
        if (self.n_cells is 0 or self.n_cells is None) and self.cell_length > 0:
                cells = round((self.length-self.cell_length)/self.cell_length)
                cells = int(cells - (cells % 3))
        elif self.n_cells > 0 and (self.cell_length is not None and self.cell_length) > 0:
            if self.cell_length == self.length:
                cells = 1
            else:
                cells = int(self.n_cells - (self.n_cells % 3))
        else:
            cells = None
        return cells

    def write_ASTRA(self, n):
        return self._write_ASTRA(OrderedDict([
            ['C_pos', {'value': self.start[2] + self.dz, 'default': 0}],
            ['FILE_EFieLD', {'value': ('\''+expand_substitution(self, '\''+self.field_definition+'\'').strip('\'"')+'\'').replace('\\','/'), 'default': 0}],
            ['C_numb', {'value': self.cells}],
            ['Nue', {'value': self.frequency / 1e9, 'default': 2998.5}],
            ['MaxE', {'value': self.field_amplitude / 1e6, 'default': 0}],
            ['Phi', {'value': -self.phase, 'default': 0.0}],
            ['C_smooth', {'value': self.smooth, 'default': 10}],
            ['C_xoff', {'value': self.start[0] + self.dx, 'default': 0}],
            ['C_yoff', {'value': self.start[1] + self.dy, 'default': 0}],
            ['C_xrot', {'value': self.y_rot + self.dy_rot, 'default': 0}],
            ['C_yrot', {'value': self.x_rot + self.dx_rot, 'default': 0}],
            ['C_zrot', {'value': self.z_rot + self.dz_rot, 'default': 0}],
        ]), n)

    def _write_Elegant(self):
        self.update_field_definition()
        wholestring=''
        etype = self._convertType_Elegant(self.objecttype)
        if (not hasattr(self, 'longitudinal_wakefield_sdds') or self.longitudinal_wakefield_sdds == None) and (not hasattr(self, 'transverse_wakefield_sdds') or self.transverse_wakefield_sdds == None):
            # print('cavity ', self.objectname, ' is an RFCA!')
            etype = 'rfca'
        string = self.objectname+': '+ etype
        for key, value in list(merge_two_dicts(self.objectproperties, self.objectdefaults).items()):
            if not key is 'name' and not key is 'type' and not key is 'commandtype' and self._convertKeword_Elegant(key) in elements_Elegant[etype]:
                value = getattr(self, key) if hasattr(self, key) and getattr(self, key) is not None else value
                key = self._convertKeword_Elegant(key)
                if self.objecttype == 'cavity':
                    # In ELEGANT all phases are +90degrees!!
                    value = 90 - value if key.lower() == 'phase' else value
                    # If using rftmez0 or similar
                    # value = ((90+value)/360.0)*(2*3.14159) if key.lower() == 'phase' else value
                    # In ELEGANT the voltages  need to be compensated
                    value = (self.cells+4.7) * self.cell_length * (1 / np.sqrt(2)) * value if key.lower() == 'volt' else value
                    # If using rftmez0 or similar
                    value = 1/(2**0.5) * value if key.lower() == 'ez' else value
                    # In CAVITY NKICK = n_cells
                    value = 3*self.cells if key.lower() == 'n_kicks' else value
                    if key.lower() == 'n_bins' and value > 0:
                        print('WARNING: Cavity n_bins is not zero - check log file to ensure correct behaviour!')
                    value = 1 if value is True else value
                    value = 0 if value is False else value
                tmpstring = ', '+key+' = '+str(value)
                if len(string+tmpstring) > 76:
                    wholestring+=string+',&\n'
                    string=''
                    string+=tmpstring[2::]
                else:
                    string+= tmpstring
        wholestring+=string+';\n'
        return wholestring

    def write_GPT(self, Brho, ccs="wcs", *args, **kwargs):
        self.update_field_definition()
        relpos, relrot = ccs.relative_position(self.start, self.global_rotation)
        relpos = relpos + [0, 0, 0]
        coord = self.gpt_coordinates(relpos, relrot)
        '''
        map1D_TM("wcs","z",linacposition,"mockup2m.gdf","Z","Ez",ffacl,phil,w);
        '''
        if self.gpt_phase_offset is None:
            self.gpt_phase_offset = 0
        output = 'f = ' + str(self.frequency) +';\n' + \
        'w = 2*pi*f;\n' + \
        'phi = ' + str(self.gpt_phase_offset + self.phase) + '/deg;\n' + \
        'ffac = ' + str(self.field_amplitude)+';\n' + \
        'map1D_TM' + '( ' + ccs.name + ', "z", '+ str(relpos[2]) +', \"'+str(expand_substitution(self,self.field_definition_gdf).strip('\'"')).replace('\\','/')+'", "Z","Ez", ffac, phi, w);\n'
        return output

class rf_deflecting_cavity(cavity):

    def __init__(self, name=None, type='rf_deflecting_cavity', **kwargs):
        super(rf_deflecting_cavity, self).__init__(name, type, **kwargs)
        self.add_default('n_kicks', 10)

class solenoid(frameworkElement):

    def __init__(self, name=None, type='solenoid', **kwargs):
        super(solenoid, self).__init__(name, type, **kwargs)

    def write_ASTRA(self, n):
        return self._write_ASTRA(OrderedDict([
            ['S_pos', {'value': self.start[2] + self.dz, 'default': 0}],
            ['FILE_BFieLD', {'value': (''+expand_substitution(self, '\''+self.field_definition+'\'')+'').replace('\\','/')}],
            ['MaxB', {'value': self.field_amplitude, 'default': 0}],
            ['S_smooth', {'value': self.smooth, 'default': 10}],
            ['S_xoff', {'value': self.start[0] + self.dx, 'default': 0}],
            ['S_yoff', {'value': self.start[1] + self.dy, 'default': 0}],
            ['S_xrot', {'value': self.y_rot + self.dy_rot, 'default': 0}],
            ['S_yrot', {'value': self.x_rot + self.dx_rot, 'default': 0}],
        ]), n)

    def write_GPT(self, Brho, ccs="wcs", *args, **kwargs):
        relpos, relrot = ccs.relative_position(self.start, self.global_rotation)
        relpos = relpos + [0, 0, 0]
        coord = self.gpt_coordinates(relpos, relrot)
        '''
        map1D_B("wcs",xOffset,0,zOffset+0.,cos(angle),0,-sin(angle),0,1,0,"bas_sol_norm.gdf","Z","Bz",gunSolField);
        '''
        if self.gpt_phase_offset is None:
            self.gpt_phase_offset = 0
        output = 'map1D_B' + '( ' + ccs.name + ', "z", '+ str(relpos[2]) + \
        ', \"'+str(expand_substitution(self, self.field_definition_gdf).strip('\'"')).replace('\\','/') + \
        '", "Z","Bz", '+ str(self.field_amplitude) + ');\n'
        return output

class aperture(frameworkElement):

    def __init__(self, name=None, type='aperture', **kwargs):
        super(aperture, self).__init__(name, type, **kwargs)
        self.number_of_elements = 1

    def write_GPT(self, Brho, ccs="wcs", *args, **kwargs):
        return ''
        # if self.shape == 'elliptical':
        #     output = 'rmax'
        # else:
        #     output = 'xymax'
        # output += '( "wcs", '+self.gpt_coordinates()+', '+str(self.horizontal_size)+', '+str(self.length)+');\n'
        # return output

    def write_ASTRA_Common(self, dic):
        dic['Ap_Z1'] = {'value': self.start[2] + self.dz, 'default': 0}
        end = self.end[2] + self.dz if self.end[2] > self.start[2] else self.start[2] + self.dz + 1e-3
        dic['Ap_Z2'] = {'value': end, 'default': 0}
        dic['A_xrot'] = {'value': self.y_rot + self.dy_rot, 'default': 0}
        dic['A_yrot'] = {'value': self.x_rot + self.dx_rot, 'default': 0}
        dic['A_zrot'] = {'value': self.z_rot + self.dz_rot, 'default': 0}
        return dic

    def write_ASTRA_Circular(self, n):
        dic = OrderedDict()
        dic['File_Aperture'] = {'value': 'RAD'}
        if self.radius is not None:
            radius = self.radius
        elif self.horizontal_size > 0 and self.vertical_size > 0:
            radius = min([self.horizontal_size, self.vertical_size])
        elif self.horizontal_size > 0:
            radius = self.horizontal_size
        elif self.vertical_size > 0:
            radius = self.vertical_size
        else:
            radius = 1
        dic['Ap_R'] = {'value': 1e3*radius}
        return self.write_ASTRA_Common(dic)

    def write_ASTRA_Planar(self, n, plane, width):
        dic = OrderedDict()
        dic['File_Aperture'] = {'value': plane}
        dic['Ap_R'] = {'value': width}
        return self.write_ASTRA_Common(dic)

    def write_ASTRA(self, n):
        self.number_of_elements = 1
        if self.shape == 'elliptical' or self.shape == 'circular':
            dic = self.write_ASTRA_Circular(n)
            return self._write_ASTRA(dic, n)
        elif self.shape == 'planar' or self.shape == 'rectangular':
            text = ''
            if self.horizontal_size > 0:
                dic = self.write_ASTRA_Planar(n, 'Col_X', 1e3*self.horizontal_size)
                text += self._write_ASTRA(dic, n)
                n = n + 1
                self.number_of_elements = self.number_of_elements + 1
            if self.vertical_size > 0:
                dic = self.write_ASTRA_Planar(n, 'Col_Y', 1e3*self.vertical_size)
                if self.number_of_elements > 1:
                    text += '\n'
                text += self._write_ASTRA(dic, n)
            return text

class scatter(frameworkElement):

    def __init__(self, name=None, type='scatter', **kwargs):
        super(scatter, self).__init__(name, type, **kwargs)
        # print('Scatter object ', self.objectname,' - DP = ', self.objectproperties)

    def _write_Elegant(self):
        wholestring=''
        etype = 'scatter'
        string = self.objectname+': '+ etype
        k1 = self.k1 if self.k1 is not None else 0
        for key, value in list(merge_two_dicts({'k1': k1}, merge_two_dicts(self.objectproperties, self.objectdefaults)).items()):
            if not key is 'name' and not key is 'type' and not key is 'commandtype' and self._convertKeword_Elegant(key) in elements_Elegant[etype]:
                value = getattr(self, key) if hasattr(self, key) and getattr(self, key) is not None else value
                key = self._convertKeword_Elegant(key)
                tmpstring = ', '+key+' = '+str(value)
                if len(string+tmpstring) > 76:
                    wholestring+=string+',&\n'
                    string=''
                    string+=tmpstring[2::]
                else:
                    string+= tmpstring
        wholestring+=string+';\n'
        return wholestring

class cleaner(frameworkElement):

    def __init__(self, name=None, type='scatter', **kwargs):
        super(cleaner, self).__init__(name, type, **kwargs)
        # print('Scatter object ', self.objectname,' - DP = ', self.objectproperties)

    def _write_Elegant(self):
        wholestring=''
        etype = 'clean'
        string = self.objectname+': '+ etype
        for key, value in merge_two_dicts(self.objectproperties, self.objectdefaults).items():
            if not key is 'name' and not key is 'type' and not key is 'commandtype' and self._convertKeword_Elegant(key) in elements_Elegant[etype]:
                value = getattr(self, key) if hasattr(self, key) and getattr(self, key) is not None else value
                key = self._convertKeword_Elegant(key)
                tmpstring = ', '+key+' = '+str(value)
                if len(string+tmpstring) > 76:
                    wholestring+=string+',&\n'
                    string=''
                    string+=tmpstring[2::]
                else:
                    string+= tmpstring
        wholestring+=string+';\n'
        return wholestring

class wall_current_monitor(frameworkElement):

    def __init__(self, name=None, type='wall_current_monitor', **kwargs):
        super(wall_current_monitor, self).__init__(name, type, **kwargs)

class integrated_current_transformer(wall_current_monitor):

    def __init__(self, name=None, type='integrated_current_transformer', **kwargs):
        super(integrated_current_transformer, self).__init__(name, type, **kwargs)

class screen(frameworkElement):

    def __init__(self, name=None, type='screen', **kwargs):
        super(screen, self).__init__(name, type, **kwargs)
        if 'output_filename' not in kwargs:
            self.output_filename = str(self.objectname)+'.sdds'

    def write_ASTRA(self, n):
        return self._write_ASTRA(OrderedDict([
            ['Screen', {'value': self.middle[2], 'default': 0}],
            ['Scr_xrot', {'value': self.y_rot + self.dy_rot, 'default': 0}],
            ['Scr_yrot', {'value': self.x_rot + self.dx_rot, 'default': 0}],
        ]), n)

    def _write_Elegant(self):
        wholestring=''
        etype = self._convertType_Elegant(self.objecttype)
        string = self.objectname+': '+ etype
        # if self.length > 0:
        #     d = drift(self.objectname+'-drift-01', type='drift', **{'length': self.length/2})
        #     wholestring+=d._write_Elegant()
        for key, value in list(merge_two_dicts(self.objectproperties, self.objectdefaults).items()):
            if not key is 'name' and not key is 'type' and not key is 'commandtype' and self._convertKeword_Elegant(key) in elements_Elegant[etype]:
                value = getattr(self, key) if hasattr(self, key) and getattr(self, key) is not None else value
                key = self._convertKeword_Elegant(key)
                tmpstring = ', '+key+' = '+str(value)
                if len(string+tmpstring) > 76:
                    wholestring+=string+',&\n'
                    string=''
                    string+=tmpstring[2::]
                else:
                    string+= tmpstring
        wholestring+=string+';\n'
        # if self.length > 0:
        #     d = drift(self.objectname+'-drift-02', type='drift', **{'length': self.length/2})
        #     wholestring+=d._write_Elegant()
        return wholestring

    def write_CSRTrack(self, n):
        z = self.middle[2]
        return """quadrupole{\nposition{rho="""+str(z)+""", psi=0.0, marker=screen"""+str(n)+"""a}\nproperties{strength=0.0, alpha=0, horizontal_offset=0,vertical_offset=0}\nposition{rho="""+str(z+1e-6)+""", psi=0.0, marker=screen"""+str(n)+"""b}\n}\n"""

    def write_GPT(self, Brho, ccs="wcs", *args, **kwargs):
        relpos, relrot = ccs.relative_position(self.position_start, self.global_rotation)
        relpos = relpos + [0, 0, self.length/2.]
        coord = self.gpt_coordinates(relpos, relrot)
        self.gpt_screen_position = relpos[2]
        output = 'screen( ' + ccs.name + ', "I", '+ str(relpos[2]) +');\n'
        return output

    def astra_to_hdf5(self, lattice):
        master_run_no = self.global_parameters['run_no'] if 'run_no' in self.global_parameters else 1
        astrabeamfilename = None
        for i in [0, -0.001, 0.001]:
            tempfilename = lattice + '.' + str(int(round((self.middle[2]+i-self.zstart[2])*100))).zfill(4) + '.' + str(master_run_no).zfill(3)
            if os.path.isfile(self.global_parameters['master_subdir'] + '/' + tempfilename):
                astrabeamfilename = tempfilename
        if astrabeamfilename is None:
            print(( 'Screen Error: ', lattice, self.middle[2], self.zstart[2]))
        else:
            self.global_parameters['beam'].read_astra_beam_file((self.global_parameters['master_subdir'] + '/' + astrabeamfilename).strip('\"'), normaliseZ=False)
            self.global_parameters['beam'].rotate_beamXZ(-1*self.starting_rotation, preOffset=[0,0,0], postOffset=-1*np.array(self.starting_offset))
            HDF5filename = (self.objectname+'.hdf5').strip('\"')
            self.global_parameters['beam'].write_HDF5_beam_file(self.global_parameters['master_subdir'] + '/' + HDF5filename, centered=False, sourcefilename=astrabeamfilename, pos=self.middle)

    def sdds_to_hdf5(self):
        elegantbeamfilename = self.output_filename.replace('.sdds','.SDDS').strip('\"')
        self.global_parameters['beam'].read_SDDS_beam_file(self.global_parameters['master_subdir'] + '/' + elegantbeamfilename)
        HDF5filename = self.output_filename.replace('.sdds','.hdf5').replace('.SDDS','.hdf5').strip('\"')
        self.global_parameters['beam'].write_HDF5_beam_file(self.global_parameters['master_subdir'] + '/' + HDF5filename, centered=False, sourcefilename=elegantbeamfilename, pos=self.middle, zoffset=self.end)

    def gdf_to_hdf5(self, gptbeamfilename):
        # gptbeamfilename = self.objectname + '.' + str(int(round((self.allElementObjects[self.end].position_end[2])*100))).zfill(4) + '.' + str(master_run_no).zfill(3)
        try:
            # print('Converting screen', self.objectname,'at', self.gpt_screen_position)
            self.global_parameters['beam'].read_gdf_beam_file(self.global_parameters['master_subdir'] + '/' + gptbeamfilename, position=self.gpt_screen_position)
            HDF5filename = self.objectname+'.hdf5'
            self.global_parameters['beam'].write_HDF5_beam_file(self.global_parameters['master_subdir'] + '/' + HDF5filename, centered=False, sourcefilename=gptbeamfilename)
        except:
            print('Error with screen', self.objectname,'at', self.gpt_screen_position)

class monitor(screen):

    def __init__(self, name=None, type='monitor', **kwargs):
        super(monitor, self).__init__(name, type, **kwargs)

class watch_point(screen):

    def __init__(self, name=None, type='watch_point', **kwargs):
        super(watch_point, self).__init__(name, 'screen', **kwargs)

class beam_position_monitor(screen):

    def __init__(self, name=None, type='beam_position_monitor', **kwargs):
        super(beam_position_monitor, self).__init__(name, type, **kwargs)

    def write_ASTRA(self, n):
        return self._write_ASTRA(OrderedDict([
            ['Screen', {'value': self.middle[2], 'default': 0}],
            ['Scr_xrot', {'value': self.y_rot + self.dy_rot, 'default': 0}],
            ['Scr_yrot', {'value': self.x_rot + self.dx_rot, 'default': 0}],
        ]), n)

class beam_arrival_monitor(screen):

    def __init__(self, name=None, type='beam_arrival_monitor', **kwargs):
        super(beam_arrival_monitor, self).__init__(name, type, **kwargs)

    def write_ASTRA(self, n):
        return ''

class collimator(aperture):

    def __init__(self, name=None, type='collimator', **kwargs):
        super(collimator, self).__init__(name, type, **kwargs)

class marker(screen):

    def __init__(self, name=None, type='marker', **kwargs):
        super(marker, self).__init__(name, 'screen', **kwargs)

    def write_CSRTrack(self, n):
        return ''

class drift(frameworkElement):

    def __init__(self, name=None, type='drift', **kwargs):
        super(drift, self).__init__(name, type, **kwargs)

    # def _write_Elegant(self):
    #     wholestring=''
    #     etype = self._convertType_Elegant(self.objecttype)
    #     string = self.objectname+': '+ etype
    #     for key, value in list(merge_two_dicts(self.objectproperties, self.objectdefaults).items()):
    #         if not key is 'name' and not key is 'type' and not key is 'commandtype' and self._convertKeword_Elegant(key) in elements_Elegant[etype]:
    #             value = getattr(self, key) if hasattr(self, key) and getattr(self, key) is not None else value
    #             key = self._convertKeword_Elegant(key)
    #             value = 1 if value is True else value
    #             value = 0 if value is False else value
    #             tmpstring = ', '+key+' = '+str(value)
    #             if len(string+tmpstring) > 76:
    #                 wholestring+=string+',&\n'
    #                 string=''
    #                 string+=tmpstring[2::]
    #             else:
    #                 string+= tmpstring
    #     wholestring+=string+';\n'
    #     return wholestring

class csrdrift(frameworkElement):

    def __init__(self, name=None, type='csrdrift', **kwargs):
        super(csrdrift, self).__init__(name, type, **kwargs)
        self.add_default('lsc_interpolate', 1)

    def _write_Elegant(self):
        wholestring=''
        etype = self._convertType_Elegant(self.objecttype)
        string = self.objectname+': '+ etype
        for key, value in list(merge_two_dicts(self.objectproperties, self.objectdefaults).items()):
            if not key is 'name' and not key is 'type' and not key is 'commandtype' and self._convertKeword_Elegant(key) in elements_Elegant[etype]:
                value = getattr(self, key) if hasattr(self, key) and getattr(self, key) is not None else value
                key = self._convertKeword_Elegant(key)
                value = 1 if value is True else value
                value = 0 if value is False else value
                tmpstring = ', '+key+' = '+str(value)
                if len(string+tmpstring) > 76:
                    wholestring+=string+',&\n'
                    string=''
                    string+=tmpstring[2::]
                else:
                    string+= tmpstring
        wholestring+=string+';\n'
        return wholestring

class lscdrift(frameworkElement):

    def __init__(self, name=None, type='lscdrift', **kwargs):
        super(lscdrift, self).__init__(name, type, **kwargs)

    def _write_Elegant(self):
        wholestring=''
        etype = self._convertType_Elegant(self.objecttype)
        string = self.objectname+': '+ etype
        for key, value in list(merge_two_dicts(self.objectproperties, self.objectdefaults).items()):
            if not key is 'name' and not key is 'type' and not key is 'commandtype' and self._convertKeword_Elegant(key) in elements_Elegant[etype]:
                value = getattr(self, key) if hasattr(self, key) and getattr(self, key) is not None else value
                key = self._convertKeword_Elegant(key)
                value = 1 if value is True else value
                value = 0 if value is False else value
                tmpstring = ', '+key+' = '+str(value)
                if len(string+tmpstring) > 76:
                    wholestring+=string+',&\n'
                    string=''
                    string+=tmpstring[2::]
                else:
                    string+= tmpstring
        wholestring+=string+';\n'
        return wholestring

class shutter(csrdrift):

    def __init__(self, name=None, type='shutter', **kwargs):
        super(shutter, self).__init__(name, type, **kwargs)

class valve(csrdrift):

    def __init__(self, name=None, type='valve', **kwargs):
        super(valve, self).__init__(name, type, **kwargs)

class bellows(csrdrift):

    def __init__(self, name=None, type='bellows', **kwargs):
        super(bellows, self).__init__(name, type, **kwargs)

class fel_modulator(frameworkElement):

    def __init__(self, name=None, type='modulator', **kwargs):
        super(fel_modulator, self).__init__(name, type, **kwargs)
        self.add_default('k1l', 0)
        self.add_default('n_steps', 1*self.periods)

    def write_ASTRA(self, n):
        return self._write_ASTRA(OrderedDict([
            ['Q_pos', {'value': self.middle[2] + self.dz, 'default': 0}],
        ]), n)

    def _write_Elegant(self):
        wholestring=''
        etype = self._convertType_Elegant(self.objecttype)
        string = self.objectname+': '+ etype
        for key, value in list(merge_two_dicts(self.objectproperties, self.objectdefaults).items()):
            if not key is 'name' and not key is 'type' and not key is 'commandtype' and self._convertKeword_Elegant(key) in elements_Elegant[etype]:
                value = getattr(self, key) if hasattr(self, key) and getattr(self, key) is not None else value
                key = self._convertKeword_Elegant(key)
                tmpstring = ', '+key+' = '+str(value)
                if len(string+tmpstring) > 76:
                    wholestring+=string+',&\n'
                    string=''
                    string+=tmpstring[2::]
                else:
                    string+= tmpstring
        wholestring+=string+';\n'
        return wholestring

class wiggler(frameworkElement):

    def __init__(self, name=None, type='wiggler', **kwargs):
        super(wiggler, self).__init__(name, type, **kwargs)
        # self.add_default('k1l', 0)
        # self.add_default('n_steps', 1*self.periods)

    def write_ASTRA(self, n):
        return self._write_ASTRA(OrderedDict([
            ['Q_pos', {'value': self.middle[2] + self.dz, 'default': 0}],
        ]), n)

    def _write_Elegant(self):
        wholestring=''
        if ('k' in self and abs(self.k) > 0) or ('peak_field' in self and abs(self.peak_field) > 0) or ('radius' in self and abs(self.radius) > 0):
            etype = self._convertType_Elegant(self.objecttype)
        else:
            etype = 'drift'
        string = self.objectname+': '+ etype
        for key, value in list(merge_two_dicts(self.objectproperties, self.objectdefaults).items()):
            if not key is 'name' and not key is 'type' and not key is 'commandtype' and self._convertKeword_Elegant(key) in elements_Elegant[etype]:
                value = getattr(self, key) if hasattr(self, key) and getattr(self, key) is not None else value
                key = self._convertKeword_Elegant(key)
                tmpstring = ', '+key+' = '+str(value)
                if len(string+tmpstring) > 76:
                    wholestring+=string+',&\n'
                    string=''
                    string+=tmpstring[2::]
                else:
                    string+= tmpstring
        wholestring+=string+';\n'
        return wholestring

class charge(frameworkElement):
    def __init__(self, name=None, type='charge', **kwargs):
        super(charge, self).__init__(name, 'charge', **kwargs)

class global_error(frameworkElement):

    def __init__(self, name=None, type='global_error', **kwargs):
        super(global_error, self).__init__(name, 'global_error', **kwargs)
        # self._errordict = {}

    def add_Error(self, type, sigma):
        if type in global_Error_Types:
            self.add_property(type, sigma)

    def write_ASTRA(self):
        return self._write_ASTRA(OrderedDict([[key, {'value': value}] for key, value in self._errordict]))

    def write_GPT(self, Brho, ccs="wcs", *args, **kwargs):
        relpos, relrot = ccs.relative_position(self.middle, [0,0,0])
        coord = self.gpt_coordinates(relpos, relrot)
        output = str(self.objecttype) + '( '+ ccs.name +', '+ coord +', '+str(self.length)+', '+str(Brho*self.k1)+');\n'
        return output

class longitudinal_wakefield(cavity):

    def __init__(self, name=None, type='longitudinal_wakefield', **kwargs):
        super(longitudinal_wakefield, self).__init__(name, type, **kwargs)
        self.add_default('coupling_cell_length', 0)

    def write_ASTRA(self, startn):
        self.update_field_definition()
        current_bins = self.current_bins if self.current_bins > 0 else 11
        output = ''
        if self.scale_kick > 0:
            for n in range(startn, startn+self.cells):
                output += self._write_ASTRA(OrderedDict([
                    ['Wk_Type', {'value': self.waketype, 'default': '\'Taylor_Method_F\''}],
                    ['Wk_filename', {'value': ('\''+expand_substitution(self, '\''+self.field_definition+'\'').strip('\'"')+'\'').replace('\\','/'), 'default': 0}],
                    ['Wk_x', {'value': self.x_offset, 'default': 0}],
                    ['Wk_y', {'value': self.y_offset, 'default': 0}],
                    ['Wk_z', {'value': self.start[2] + self.coupling_cell_length + (n-1)*self.cell_length}],
                    ['Wk_ex', {'value': self.scale_field_ex, 'default': 0}],
                    ['Wk_ey', {'value': self.scale_field_ey, 'default': 0}],
                    ['Wk_ez', {'value': self.scale_field_ez, 'default': 1}],
                    ['Wk_hx', {'value': self.scale_field_hx, 'default': 1}],
                    ['Wk_hy', {'value': self.scale_field_hy, 'default': 0}],
                    ['Wk_hz', {'value': self.scale_field_hz, 'default': 0}],
                    ['Wk_equi_grid', {'value': self.equal_grid, 'default': 0}],
                    ['Wk_N_bin', {'value': current_bins, 'default': 11}],
                    ['Wk_ip_method', {'value': self.interpolation_method, 'default': 2}],
                    ['Wk_smooth', {'value': self.smooth, 'default': 0.5}],
                    ['Wk_sub', {'value': self.subbins, 'default': 4}],
                    ['Wk_scaling', {'value': self.scale_kick, 'default': 1}],
                ]), n)
                output += '\n'
        return output

    def _write_Elegant(self):
        self.update_field_definition()
        wholestring=''
        etype = self._convertType_Elegant(self.objecttype)
        string = self.objectname+': '+ etype
        if self.length > 0:
            d = drift(self.objectname+'-drift', type='drift', **{'length': self.length})
            wholestring+=d._write_Elegant()
        for key, value in list(merge_two_dicts(self.objectproperties, self.objectdefaults).items()):
            if not key is 'name' and not key is 'type' and not key is 'commandtype' and self._convertKeword_Elegant(key) in elements_Elegant[etype]:
                value = getattr(self, key) if hasattr(self, key) and getattr(self, key) is not None else value
                key = self._convertKeword_Elegant(key)
                tmpstring = ', '+key+' = '+str(value)
                if len(string+tmpstring) > 76:
                    wholestring+=string+',&\n'
                    string=''
                    string+=tmpstring[2::]
                else:
                    string+= tmpstring
        wholestring+=string+';\n'
        return wholestring

class gpt_ccs(Munch):

    def __init__(self, name, position, rotation, intersect=0):
        super(gpt_ccs, self).__init__()
        self._name = name
        self.intersect = intersect
        self.x, self.y, self.z = position
        self.psi, self.phi, self.theta = rotation

    def relative_position(self, position, rotation):
        x, y, z = position
        psi, phi, theta = rotation
        # print(self.name, [x - self.x, y - self.y, z - self.z])
        # print(self.name, [psi - self.psi, phi - self.phi, theta - self.theta])
        newpos = [x - self.x, y - self.y, z - self.z]
        # print('newpos = ', self.name,  x, self.x, y, self.y, z, self.z)
        finalrot = [psi - self.psi, phi - self.phi, theta - self.theta]
        finalpos = np.array([0,0,self.intersect]) + np.dot(np.array(newpos), _rotation_matrix(-self.theta))
        return finalpos, finalrot

    @property
    def name(self):
        return '"' + self._name + '"'
    @property
    def position(self):
        return self.x, self.y, self.z
    @property
    def rotation(self):
        return self.psi, self.phi, self.theta
