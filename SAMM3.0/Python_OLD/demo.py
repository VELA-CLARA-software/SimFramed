from Components import *
from SAMMlab  import Beamline
from Particles  import *


drift1 = Drift.Drift(length = 1.5)
quad1 = Quadrupole.Quadrupole(length=0.3, gradient=0.8)

beamline1 = Beamline.Beamline( componentlist=[drift1,quad1])

beam1 = Beam.Beam(Positron)