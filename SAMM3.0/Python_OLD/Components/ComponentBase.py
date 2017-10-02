# SAM to Python Conversion
# DJS August 2017
# Version 0.1
#

class ComponentBase(object):
    def __init__(self, length=0, name="", aperture=[]):
        #super(ComponentBase, self).__init__(**kwargs)
        # device length, in meters
        self.length  = length

        # device name, string
        self.name = name

        # 1x2 array of elliptical aperture half-axes, in metres
        self.aperture = aperture

