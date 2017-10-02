function [I1 I2 I3 I4 I5] = ComputeRadiationIntegrals(beamline,beam)

    I1 = 0;
    I2 = 0;
    I3 = 0;
    I4 = 0;
    I5 = 0;
    
    beta = ComputeMatchedTwiss(beamline,beam);
    
    for n = 1:length(beamline.componentlist)
       
        component = beamline.componentlist{n};
        
        if(isa(component,'Dipole'))
            
            ds = component.length;
            h  = component.curvature;
            k1 = component.gradient / beam.rigidity;
            e1 = component.e1;
            
            theta  = h*ds;
            kx     = sqrt(h*h + k1);
            
            fringe = [1         0; ...
                      h*tan(e1) 1];

            hvec0  = [beta(1,6,3,n)/beta(6,6,3,n); ...
                      beta(2,6,3,n)/beta(6,6,3,n)];
                  
            hvec0  = fringe*hvec0;
            eta0   = hvec0(1);
            etap0  = hvec0(2);

            a0     = beta(1:2,1:2,1,n); 
            
            a0     = fringe*a0*fringe';
            alpha0 =-a0(1,2);
            beta0  = a0(1,1);
            gamma0 = a0(2,2);
            
            inthds = ds - sin(ds*h)/h ...
                  + eta0*h*sin(kx*ds)/kx ...
                  + etap0*h*(1 - cos(kx*ds))/kx/kx;
            
            I1 = I1 + inthds;
                
            I2 = I2 + ds * h * h;
            
            I3 = I3 + ds * abs(h)^3;
            
            I4 = I4 + inthds * (h*h + 2*k1);
            
            intHds = (gamma0*eta0*eta0 + 2*alpha0*eta0*etap0 + beta0*etap0*etap0)*theta/h...
                    + (alpha0*eta0 + beta0*etap0)*theta^2/h ...
                    + (beta0 - alpha0*etap0/h - gamma0*eta0/h)*theta^3/3/h ...
                    - (alpha0/h - alpha0*eta0 - beta0*etap0)*theta^4/12/h;
            
            I5 = I5 + intHds * abs(h)^3;
            
        end % if
        
    end % for
    
return
    