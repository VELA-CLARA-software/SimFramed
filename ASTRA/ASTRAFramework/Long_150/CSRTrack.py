import os, subprocess
import numpy as np

class CSRTrack(object):

    def __init__(self, subdir='test', overwrite=None, runname='CLARA_240'):
        super(CSRTrack, self).__init__()
        self.defineCSRTrackCommand(['/opt/CSRTrack/csrtrack_1.204_Linux_x86_64_serial csrtrk.in'])
        self.basedirectory = os.getcwd()
        self.subdir = subdir
        self.overwrite = overwrite
        self.subdirectory = self.basedirectory+'/'+subdir
        if not os.path.exists(self.subdirectory):
            os.makedirs(self.subdirectory)
        if self.overwrite == None:
            if not os.path.exists(self.subdirectory):
                os.makedirs(self.subdirectory)
                self.overwrite = True
            else:
                response = raw_input('Overwrite existing directory ? [Y/n]')
                self.overwrite = True if response in {'','y','Y','yes','Yes','YES'} else False

    def defineCSRTrackCommand(self,command=['csrtrack']):
        self.CSRTrackCommand = command

    def runCSRTrackFile(self, filename):
        if self.overwrite:
            command = self.CSRTrackCommand + [filename]
            with open(os.devnull, "w") as f:
                subprocess.call(command, stdout=f, cwd=self.subdir)

    def csrtrackinputtext(self, angle, forces='projected', inputfile='./long_150.4.2573.001'):
        ANGLE = angle
        return """io_path{logfile = log.txt}

    lattice{
    dipole{
    position{rho="""+str(25.969833)+""", psi=0.0, marker=d1a}
    properties{r="""+str(-0.2009814751607801/ANGLE)+"""}
    position{rho="""+str(25.969833 + (0.2009814751607801*np.sin(1.*ANGLE))/ANGLE)+""", psi=0.0, marker=d1b}
    }
    dipole{
    position{rho="""+str(25.969833 + 1.5067942970813246*np.cos(ANGLE) + (0.2009814751607801*np.sin(ANGLE))/ANGLE)+""", psi=0.0, marker=d2a}
    properties{r="""+str(0.2009814751607801/ANGLE)+"""}
    position{rho="""+str((25.969833*ANGLE + 1.5067942970813246*ANGLE*np.cos(ANGLE) + 0.2009814751607801*np.sin(ANGLE) + 0.2009814751607801*np.sin(1.*ANGLE))/ANGLE)+""", psi=0.0, marker=d2b}
    }
    dipole{
    position{rho="""+str(1.3 + 25.969833 + 1.5067942970813246*np.cos(ANGLE) + (0.4019629503215602*np.sin(ANGLE))/ANGLE)+""", psi=0.0, marker=d3a}
    properties{r="""+str(0.2009814751607801/ANGLE)+"""}
    position{rho="""+str((27.26978776235*ANGLE + 1.5067942970813246*ANGLE*np.cos(ANGLE) + 0.4019629503215602*np.sin(ANGLE) + 0.2009814751607801*np.sin(1.*ANGLE))/ANGLE)+""", psi=0.0, marker=d3b}
    }
    dipole{
    position{rho="""+str(((1.3 + 25.969833)*ANGLE + 1.5067942970813246*ANGLE*np.cos(ANGLE) + 1.5067942970813246*ANGLE*np.cos(1.*ANGLE) + 0.4019629503215602*np.sin(ANGLE) + 0.2009814751607801*np.sin(1.*ANGLE))/ANGLE)+""", psi=0.0, marker=d4a}
    properties{r="""+str(-0.2009814751607801/ANGLE)+"""}
    position{rho="""+str(((1.3 + 25.969833)*ANGLE + 1.5067942970813246*ANGLE*np.cos(ANGLE) + 1.5067942970813246*ANGLE*np.cos(1.*ANGLE) + 0.4019629503215602*np.sin(ANGLE) + 0.4019629503215602*np.sin(1.*ANGLE))/ANGLE)+""", psi=0.0, marker=d4b}
    }

    }
    particles{
    reference_momentum = reference_particle
    reference_point_x   = 0.0
    reference_point_y   = 0.0
    reference_point_phi = 0.0
    format = astra, array = #file{name="""+ inputfile + """}}
    online_monitor{name = sub_bunch.dat, type = subbunch
    start_time_c0 = now
    end_time_marker = d4b, end_time_shift_c0 = 0.116726
    time_step_c0    = all
    }
    online_monitor{name = steps.dat, type = steps
    start_time_c0 = now
    end_time_marker = d4b, end_time_shift_c0 = 0.116726
    time_step_c0 = all
    }
    online_monitor{name = p1.fmt2, type = phase, format = fmt2, particle = 1
    start_time_c0 = now
    end_time_marker = d4b, end_time_shift_c0 = 0.116726
    time_step_c0 = all
    }
    forces{type = """ + forces + """
    shape = ellipsoid
    sigma_long = relative, relative_long = 0.1}
    track_step{ precondition=yes
    iterative=2
    error_per_ct=0.001
    error_weight_momentum=0.0
    ct_step_min=0.002
    ct_step_max=0.010
    ct_step_first=0.010
    increase_factor=1.5
    arc_factor=0.3
    duty_steps=yes}

    tracker{end_time_marker = d4b, end_time_shift_c0 = 0.116726}
    monitor{format = fmt2, name = end.fmt2}
    exit
    """

    def writeCSRTrackFile(self, filename, angle=0.105, forces='projected', inputfile='./long_150.4.2573.001'):
        with file(self.subdir+'/'+filename, 'w') as f:
            f.write(self.csrtrackinputtext(angle, forces=forces, inputfile=inputfile))
