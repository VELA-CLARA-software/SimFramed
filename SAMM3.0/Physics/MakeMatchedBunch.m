function particles = MakeMatchedBunch(orbit,beta,emittances,nparticles)

    S = [ 0  1  0  0  0  0;...
         -1  0  0  0  0  0;...
          0  0  0  1  0  0;...
          0  0 -1  0  0  0;...
          0  0  0  0  0  1;...
          0  0  0  0 -1  0];
      
      
    % Construct the correlation (sigma) matrix for the given beta functions and emittances  
    sigma = beta(:,:,1)*emittances(1) + ...
            beta(:,:,2)*emittances(2) + ...
            beta(:,:,3)*emittances(3);
      
    % Find the eigensystem of sigma*S
    [V,D] = eig(sigma*S);
    % Sort the eigenvalues and eigenvectors into conjugate pairs
    [D,indx] = sort(abs(diag(D)),'descend');

    % Normalise the eigenvectors
    w = diag(V'*S*V);
    V = V*diag(1./sqrt(w));
    
    nmat = sqrt(2)*[...
        real(V(:,1)), imag(V(:,1)), ...
        real(V(:,3)), imag(V(:,3)), ...
        real(V(:,5)), imag(V(:,5)) ];
    
    xn = diag(sqrt(D))*randn(6,nparticles);
    
    particles = orbit*ones(1,nparticles) + nmat*xn;

return
