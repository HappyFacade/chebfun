function [x, disc] = mldivide(disc, M, b)
% ULTRAS.MLDIVIDE
% Solver for ultraspherical spectral method.

if ( length( disc.dimension ) > 1 )

    x = M \ b;

else
    % We use the Woodbury formula, i.e., when M = A + UV
    % we have
    %    M\b  = A\b - A\(U*((I+V*(A\U))\(V*(A\b))
    % We split M into banded + low-rank to get an efficient solver.
    
    % How many boundary rows?
    k = sum( disc.projOrder );
    [m, n] = size( M );
    
    % Extract banded part:
    A = M;
    A(1:k,:) = speye(k, n);
    
    % Extract low-rank part:
    U = speye(m, k);
    V = M(1:k, :) - speye(k, n);
    I = speye(k);
    
    % Woodbury formula applied to M\b = (A+U*V)\b:
    Ab = A \ b;
    cfs = (I+V*(A\U))\(V*Ab);
    x = Ab - A\(U*cfs);
    
end

end