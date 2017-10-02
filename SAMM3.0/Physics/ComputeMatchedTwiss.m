function [beta, tune, closedorbit] = ComputeMatchedTwiss(beamline,beam)

% [beta tune closedorbit] = FindMatchedTwiss(beamline,rigidity)
%   Find the lattice functions around the ring, the betatron and
%   synchrotron tunes, and the closed orbit.
% 
% beta(n,i,j,k) is a four dimensional array containing the lattice functions:
%   - n is an element index specifying the location in the lattice;
%   - i,j are indices of a 6x6 matrix (see below);
%   - k (=1,2,3) is an index specifying a degree of freedom.
% 
% tune(k) is a list of the three tunes, with k specifying a degree of freedom.
%
% closedorbit(6,n) contains the closed orbit at each point around the ring.
%
% In an uncoupled lattice:
%   - beta(:,1,1,1) is beta_x around the ring
%   - beta(:,1,2,1) is -alpha_x
%   - beta(:,2,2,1) is gamma_x
%   - beta(:,3,3,2) is beta_y
%   - beta(:,3,4,2) is -alpha_y
%   - beta(:,4,4,2) is gamma_y

    S = [ 0  1  0  0  0  0;...
         -1  0  0  0  0  0;...
          0  0  0  1  0  0;...
          0  0 -1  0  0  0;...
          0  0  0  0  0  1;...
          0  0  0  0 -1  0];

    T1= [ 0  1  0  0  0  0;...
          1  0  0  0  0  0;...
          0  0  0  0  0  0;...
          0  0  0  0  0  0;...
          0  0  0  0  0  0;...
          0  0  0  0  0  0];
      
    T2= [ 0  0  0  0  0  0;...
          0  0  0  0  0  0;...
          0  0  0  1  0  0;...
          0  0  1  0  0  0;...
          0  0  0  0  0  0;...
          0  0  0  0  0  0];
      
    T3= [ 0  0  0  0  0  0;...
          0  0  0  0  0  0;...
          0  0  0  0  0  0;...
          0  0  0  0  0  0;...
          0  0  0  0  0  1;...
          0  0  0  0  1  0];
      
    
    closedorbit = ComputeClosedOrbit(beamline,beam);
    m   = ComputeTransferMatrix(beamline,[1 length(beamline.componentlist)],beam,closedorbit(:,1));
    
    m1  = m(:,:,end);
    
    if m1(6,5)==0 % add some longitudinal focusing
        lf  = zeros(6,6);
        lf(6,5) = 1e-12;
        m1  = (eye(6) - sign(m1(5,6))*lf)*m1;
    end

    % Find the eigensystem of the transfer matrix
    [V,D] = eig(m1(:,:,end));
    
    % Sort the eigenvalues and eigenvectors into conjugate pairs
    [~,indx] = sort(abs(angle(diag(D))),'descend');
    V = V(:,indx);
    
    % Normalise the eigenvectors
    w = diag(V'*S*V);
    V = V*diag(1./sqrt(w));
    
    % Attempt to sort the eigenvectors into horizontal, vertical, and
    % longitudinal
    vsort = zeros(3,3);
    for mi = 2:2:6
        normv = V(:,mi)'*V(:,mi);
        for ni = 2:2:6
            Vsub = V((ni-1):(ni),mi);
            vsort(ni/2,mi/2) = (Vsub'*Vsub) / normv; 
        end
    end
    [~,ix] = max(vsort,[],2);
    ix1 = zeros(1,6);
    ix1([1 3 5]) = 2*ix - 1;
    ix1([2 4 6]) = 2*ix;
    V = V(:,ix1);
    
    % Calculate the Twiss parameters at the start of the beamline...
    beta = zeros(6,6,3,length(beamline.componentlist)+1);
    beta(:,:,1,1) = V*T1*V.';
    beta(:,:,2,1) = V*T2*V.';
    beta(:,:,3,1) = V*T3*V.';
    
    % Check to see if the "transverse" beta matrices are in the expected order.
    % If not, swap them round.
    if (beta(3,3,1,1)>beta(1,1,1,1)) && (beta(1,1,2,1)>beta(3,3,2,1))
        beta = beta(:,:,[2 1 3],:);
        V    = V(:,[3 4 1 2 5 6]);
    end
    
    mu = zeros(length(beamline.componentlist)+1,3);
    
    nmat = sqrt(2)*[...
        real(V(:,1)), imag(V(:,1)), ...
        real(V(:,3)), imag(V(:,3)), ...
        real(V(:,5)), imag(V(:,5)) ];
    
    nmat = standardphase(nmat);
    
    % ...and propagate along the beamline.
    for n = 2:length(beamline.componentlist)+1
        m1 = m(:,:,n);
        beta(:,:,1,n) = m1*beta(:,:,1,1)*m1';
        beta(:,:,2,n) = m1*beta(:,:,2,1)*m1';
        beta(:,:,3,n) = m1*beta(:,:,3,1)*m1';
        
        nmat1   = standardphase(m1*nmat);
        R       = nmat1\m1*nmat;
        mu(n,:) = atan2([R(1,2),R(3,4),R(5,6)],[R(1,1),R(3,3),R(5,5)]);
    end;
    
    tune = abs(unwrap(mu))/2/pi;
    
return

function mrot = standardphase(m)

    psix = angle(m(1,2) + 1i*m(1,1)) - pi/2;
    psiy = angle(m(3,4) + 1i*m(3,3)) - pi/2;
    psiz = angle(m(5,6) + 1i*m(5,5)) - pi/2;

    rot  = [ cos(psix) sin(psix)   0         0         0         0        ;...
            -sin(psix) cos(psix)   0         0         0         0        ;...
               0         0       cos(psiy) sin(psiy)   0         0        ;...
               0         0      -sin(psiy) cos(psiy)   0         0        ;...
               0         0         0         0       cos(psiz) sin(psiz)  ;...
               0         0         0         0      -sin(psiz) cos(psiz)  ];

    mrot = m * rot;

return
