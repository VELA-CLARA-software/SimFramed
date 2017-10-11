# SAM to Python Conversion
# DJS August 2017
# Version 0.1
#
from ComponentBase import ComponentBase
import numpy

import Drift


class Quadrupole(ComponentBase):
    def __init__(self, length=0, name="", aperture=[], gradient=0):
        ComponentBase.__init__(self, length, name, aperture)

        # quadrupole gradient, in tesla/metre
        self.gradient = gradient
        self.__bx = None
        self.__by = None
        self.__bz = None
        self.k1 = None
        self.absk1 = None
        self.ds = self.length

    def Track(self, beam):
        # normalised gradient
        self.k1 = numpy.divide(self.gradient, beam.rigidity)
        self.absk1 = abs(self.k1)
        # print('Quad = ',self.gradient,beam.rigidity)
        if self.k1 > 0:
            self.TrackFQuad(beam)
        elif self.k1 < 0:
            self.TrackDQuad(beam)
        else:
            # print('Quad TrackDrift')
            self.TrackDrift(beam)

        # save
        self.lastTrackedBeam = beam

    def GetBField(self, beam):
        # [bx, by, bz] = Quadrupole.GetBField(beam)
        # Returns the magnetic field (in tesla) at the locations of
        # the particles in the beam.
        self.__bx = beam.y * self.gradient
        self.__by = beam.x * self.gradient
        self.__bz = [0.0] * len(beam.y)

    def TrackFQuad(self, beam):
        # print 'QUAD_TRACK'
        d1 = numpy.sqrt(1 + numpy.divide(2 * beam.dp, beam.beta) +
                        beam.dp * beam.dp)
        w = numpy.sqrt(numpy.divide(self.k1, d1))

        xs = numpy.sin(w * self.ds)
        xc = numpy.cos(w * self.ds)
        ys = numpy.sinh(w * self.ds)
        yc = numpy.cosh(w * self.ds)
        xs2 = numpy.sin(2 * w * self.ds)
        ys2 = numpy.sinh(2 * w * self.ds)

        d0 = 1 / beam.beta + beam.dp
        d2 = -d0 / d1 / d1 / d1 / 2

        c0 = (1 / beam.beta - d0 / d1) * self.ds
        c11 = self.k1 * self.k1 * d2 * (xs2 / w - 2 * self.ds) / w / w / 4
        c12 = -self.k1 * d2 * xs * xs / w / w
        c22 = d2 * (xs2 / w + 2 * self.ds) / 4
        c33 = self.k1 * self.k1 * d2 * (ys2 / w - 2 * self.ds) / w / w / 4
        c34 = self.k1 * d2 * ys * ys / w / w
        c44 = d2 * (ys2 / w + 2 * self.ds) / 4
        beam.ct = (beam.ct + c0
                   + c11 * beam.x * beam.x
                   + c12 * beam.x * beam.px
                   + c22 * beam.px * beam.px
                   + c33 * beam.y * beam.y
                   + c34 * beam.y * beam.py
                   + c44 * beam.py * beam.py)

        x0 = beam.x
        beam.x = beam.x * xc + beam.px * (xs * w / self.absk1)
        beam.px = -x0 * self.k1 * xs / w + beam.px * xc
        y0 = beam.y
        beam.y = beam.y * yc + beam.py * (ys * w / self.absk1)
        beam.py = y0 * self.k1 * ys / w + beam.py * yc

        # i'm sure the following could be tidied up ++ CHECK!

        # d0  = 1/beam.beta + beam.dp
        # d2  = -d0/d1/d1/d1/2
        #
        # c0  = (1/beam.beta - d0/d1 )*self.ds
        # c11 = self.k1*self.k1*d2*(xs2/w - 2*self.ds)/w/w/4
        # c12 =-self.k1*d2*xs*xs/w/w
        # c22 = d2*(xs2/w + 2*self.ds)/4
        # c33 = self.k1*self.k1*d2*(ys2/w - 2*self.ds)/w/w/4
        # c34 = self.k1*d2*ys*ys/w/w
        # c44 = d2*(ys2/w + 2*self.ds)/4
        # beam.ct = beam.ct + c0 \
        #           + c11*beam.x*beam.x   \
        #           + c12*beam.x*beam.px  \
        #           + c22*beam.px*beam.px \
        #           + c33*beam.y*beam.y   \
        #           + c34*beam.y*beam.py  \
        #           + c44*beam.py*beam.py
        # # update particles
        # beam.x  = x1
        # beam.y  = y1
        # beam.px = px1
        # beam.py = py1

    def TrackDQuad(self, beam):
        # print 'QUAD_TRACK'
        # [x0, px0, y0, py0, ct0, dp0] = beam.GetParticles();
        d1 = numpy.sqrt(1 + 2 * numpy.divide(beam.dp, beam.beta) + numpy.multiply(beam.dp, beam.dp))
        w = numpy.sqrt(self.absk1 / d1)

        xs = numpy.sinh(w * self.ds)
        xc = numpy.cosh(w * self.ds)
        ys = numpy.sin(w * self.ds)
        yc = numpy.cos(w * self.ds)
        xs2 = numpy.sinh(2 * w * self.ds)
        ys2 = numpy.sin(2 * w * self.ds)

        d0 = 1 / beam.beta + beam.dp
        d2 = -d0 / d1 / d1 / d1 / 2
        c0 = (1 / beam.beta - d0 / d1) * self.ds
        c11 = self.k1 * self.k1 * d2 * (xs2 / w - 2 * self.ds) / w / w / 4
        c12 = -self.k1 * d2 * xs * xs / w / w
        c22 = d2 * (xs2 / w + 2 * self.ds) / 4
        c33 = self.k1 * self.k1 * d2 * (ys2 / w - 2 * self.ds) / w / w / 4
        c34 = self.k1 * d2 * ys * ys / w / w
        c44 = d2 * (ys2 / w + 2 * self.ds) / 4
        beam.ct = (beam.ct + c0
                   + c11 * beam.x * beam.x
                   + c12 * beam.x * beam.px
                   + c22 * beam.px * beam.px
                   + c33 * beam.y * beam.y
                   + c34 * beam.y * beam.py
                   + c44 * beam.py * beam.py)
        x0 = beam.x
        beam.x = beam.x * xc + beam.px * (xs * w / self.absk1)
        beam.px = -x0 * self.k1 * xs / w + beam.px * xc
        y0 = beam.y
        beam.y = beam.y * yc + beam.py * (ys * w / self.absk1)
        beam.py = y0 * self.k1 * ys / w + beam.py * yc

    def TrackDrift(self, beam):
        # print 'QUAD_TRACK'
        drift = Drift.Drift(length=self.length)
        drift.Track(beam)


# classdef Quadrupole < handle
#     % Quadrupole class
#     %
#     % Properties:
#     %   name
#     %   length
#     %   gradient
#     %   aperture
#     %
#     % Methods:
#     %   Track
#     %   TrackSpin
#     %   GetBField
#
#     properties
#         name     = ''; % string
#         length   = 0;  % in metres
#         gradient = 0;  % in tesla/metre
#         aperture = []; % 1x2 array of elliptical aperture half-axes, in metres
#     end % properties
#
#     methods
#
#         function beam = Track(quadrupole,beam)
#             % beam2 = Quadrupole.Track(beam1)
#             % Applies the transfer map for a quadrupole to the particles
#             % in beam1.
#
#             k1 = quadrupole.gradient / beam.rigidity; % normalised gradient
#
#             if k1>0
#                 beam = TrackFQuad(quadrupole,beam);
#             elseif k1<0
#                 beam = TrackDQuad(quadrupole,beam);
#             else
#                 beam = TrackDrift(quadrupole,beam);
#             end
#
#         end % function Track
#
#
#         function TrackLibrary(quadrupole,trackingMethod,libroutine)
#
#             quadrupole1.length   = quadrupole.length;
#             quadrupole1.gradient = quadrupole.gradient;
#             if(~isempty(quadrupole.aperture))
#                quadrupole1.apertureX = quadrupole.aperture(1);
#                quadrupole1.apertureY = quadrupole.aperture(2);
#             end
#
#             calllib(trackingMethod,[libroutine 'Quadrupole'],quadrupole1);
#
#         end % function TrackLibrary
#
#
#         function [bx, by, bz] = GetBField(quadrupole,beam)
#             % [bx, by, bz] = Quadrupole.GetBField(beam)
#             % Returns the magnetic field (in tesla) at the locations of
#             % the particles in the beam.
#
#             [x, ~, y] = beam.GetParticles();
#
#             bx = y*quadrupole.gradient;
#             by = x*quadrupole.gradient;
#             bz = 0;
#
#         end % function GetBField
#
#
#         function beam = TrackSpin(quadrupole,beam)
#             % beam2 = Quadrupole.Trackspin(beam1)
#             % Tracks particle spins through a Quadrupole.
#
#             [bx, by, bz] = quadrupole.GetBField(beam);
#
#             Utilities.SpinRotation(beam,bx,by,bz,quadrupole.length);
#
#         end % function TrackSpin
#
#
#     end % methods
#
#
#     methods (Access='private')
#
#
#         function beam = TrackFQuad(quadrupole,beam)
#
#             [x0, px0, y0, py0, ct0, dp0] = beam.GetParticles();
#
#             beta0 = beam.beta;
#
#             ds = quadrupole.length;
#             k1 = quadrupole.gradient / beam.rigidity; % normalised gradient
#
#             d1 = sqrt(1 + 2*dp0/beta0 + dp0.*dp0);
#             w  = sqrt(k1 ./ d1);
#
#             xs  = sin(w*ds);
#             xc  = cos(w*ds);
#             ys  = sinh(w*ds);
#             yc  = cosh(w*ds);
#             xs2 = sin(2*w*ds);
# %           xc2 = cos(2*w*ds);
#             ys2 = sinh(2*w*ds);
# %           yc2 = cosh(2*w*ds);
#
#             x1  =  x0.*xc       + px0.*xs.*w/k1;
#             px1 = -k1*x0.*xs./w + px0.*xc;
#             y1  =  y0.*yc       + py0.*ys.*w/k1;
#             py1 =  k1*y0.*ys./w + py0.*yc;
#
#             d0  = 1/beta0 + dp0;
#             d2  =-d0./d1./d1./d1/2;
#
#             c0  = (1/beta0 - d0./d1)*ds;
#             c11 = k1*k1*d2.*(xs2./w - 2*ds)./w./w/4;
#             c12 =-k1*d2.*xs.*xs./w./w;
#             c22 = d2.*(xs2./w + 2*ds)/4;
#             c33 = k1*k1*d2.*(ys2./w - 2*ds)./w./w/4;
#             c34 = k1*d2.*ys.*ys./w./w;
#             c44 = d2.*(ys2./w + 2*ds)/4;
#
#             ct1 = ct0 + c0 ...
#                       + c11.* x0.* x0 ...
#                       + c12.* x0.*px0 ...
#                       + c22.*px0.*px0 ...
#                       + c33.* y0.* y0 ...
#                       + c34.* y0.*py0 ...
#                       + c44.*py0.*py0;
#
#             beam.SetParticles(x1, px1, y1, py1, ct1, dp0);
#
#         end
#
#
#         function beam = TrackDQuad(quadrupole,beam)
#
#             [x0, px0, y0, py0, ct0, dp0] = beam.GetParticles();
#
#             beta0 = beam.beta;
#
#             ds = quadrupole.length;
#             k1 = quadrupole.gradient / beam.rigidity; % normalised gradient
#
#             d1 = sqrt(1 + 2*dp0/beta0 + dp0.*dp0);
#             w  = sqrt(abs(k1) ./ d1);
#
#             xs  = sinh(w*ds);
#             xc  = cosh(w*ds);
#             ys  = sin(w*ds);
#             yc  = cos(w*ds);
#             xs2 = sinh(2*w*ds);
# %           xc2 = cosh(2*w*ds);
#             ys2 = sin(2*w*ds);
# %           yc2 = cos(2*w*ds);
#
#             x1  =  x0.*xc       + px0.*xs.*w/abs(k1);
#             px1 = -k1*x0.*xs./w + px0.*xc;
#             y1  =  y0.*yc       + py0.*ys.*w/abs(k1);
#             py1 =  k1*y0.*ys./w + py0.*yc;
#
#             d0  = 1/beta0 + dp0;
#             d2  =-d0./d1./d1./d1/2;
#
#             c0  = (1/beta0 - d0./d1)*ds;
#             c11 = k1*k1*d2.*(xs2./w - 2*ds)./w./w/4;
#             c12 =-k1*d2.*xs.*xs./w./w;
#             c22 = d2.*(xs2./w + 2*ds)/4;
#             c33 = k1*k1*d2.*(ys2./w - 2*ds)./w./w/4;
#             c34 = k1*d2.*ys.*ys./w./w;
#             c44 = d2.*(ys2./w + 2*ds)/4;
#
#             ct1 = ct0 + c0 ...
#                       + c11.* x0.* x0 ...
#                       + c12.* x0.*px0 ...
#                       + c22.*px0.*px0 ...
#                       + c33.* y0.* y0 ...
#                       + c34.* y0.*py0 ...
#                       + c44.*py0.*py0;
#
#             beam.SetParticles(x1, px1, y1, py1, ct1, dp0);
#
#         end
#
#
#         function beam = TrackDrift(quadrupole,beam)
#
#             [x0, px0, y0, py0, ct0, dp0] = beam.GetParticles();
#
#             beta0  = beam.beta;
#
#             ds  = quadrupole.length;
#
#             d1  = sqrt(1 - px0.*px0 - py0.*py0 + 2*dp0/beta0 + dp0.*dp0);
#
#             x1  = x0  + ds*px0./d1;
#             y1  = y0  + ds*py0./d1;
#             ct1 = ct0 + ds*(1 - (1 + beta0*dp0)./d1)/beta0;
#
#             beam.SetParticles(x1, px0, y1, py0, ct1, dp0);
#
#         end
#
#
#     end % private methods
#
# end % classdef Quadrupole
