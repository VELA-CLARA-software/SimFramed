function [closedorbit, m] = ComputeClosedOrbit(beamline,beam)

    m   = zeros(6,6);
    xco = zeros(6,1);
    
    lf  = zeros(6,6);
    lf(6,5) = 1e-12;
    
    itern = 1;
    resdl = 1;
    precn = beamline.precision;
    
    while (itern<10) && (resdl>10*precn*precn)
        
        p  = precn*eye(6,7);
        for n = 1:7
            p(:,n) = p(:,n) + xco;
        end

        beam.particles = p;
        beam.globaltime = 0;
        beam = beamline.Track([1 length(beamline.componentlist)],beam);

        p1 = beam.particles;
        for n = 1:6
            m(:,n) = (p1(:,n) - p1(:,7))/precn;
        end

        if m(6,5)==0 % add some longitudinal focusing
            m = (eye(6) - sign(m(5,6))*lf)*m;
        end
      
        dt = beam.globaltime * MasterOscillator.GetFrequency();
        p1(5,7) = p1(5,7) - ...
            (dt - round(dt))*beam.beta*PhysicalConstants.SpeedOfLight/MasterOscillator.GetFrequency();

        d     = p1(:,7) - xco;
        dco   = (eye(6) - m) \ d;
        resdl = dco'*dco;
        xco   = xco + dco;
        
        itern = itern + 1;
        
    end
    
    beam.particles = xco;
    beam.globaltime = 0;
    
    closedorbit = zeros(6,length(beamline.componentlist)+1);
    closedorbit(:,1) = xco;
    for n = 1:length(beamline.componentlist)
        beam = beamline.Track([n n],beam);
        closedorbit(:,n+1) = beam.particles;
    end
    
return