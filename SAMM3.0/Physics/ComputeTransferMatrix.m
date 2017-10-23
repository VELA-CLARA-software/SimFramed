function [m, eref] = ComputeTransferMatrix(beamline,range,beam,trajectory)

    eref     = beam.energy * ones(1,range(2)-range(1)+2);
    m        = zeros(6,6,range(2)-range(1)+2);
    m(:,:,1) = eye(6);

    precn = beamline.precision;
    p  = precn*eye(6,7);
    for n = 1:7
        p(:,n) = p(:,n) + trajectory;
    end

    beam.particles = p;
    p2 = p(:,7)

    for n = range(1):range(2)
        beam = beamline.Track([n n],beam);
        p1 = beam.particles
        for v = 1:6
            m(:,v,n-range(1)+2) = (p1(:,v) - p1(:,7)) / precn;
        end
        eref(n-range(1)+2) = beam.energy;
    end

return
