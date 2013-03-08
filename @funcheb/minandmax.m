function [vals, pos] = minandmax(f)
%MINANDMAX	 Global minimum and maximum on [-1,1].
%   VAL = MINANDMAX(F) returns a 2-vector VAL = [min(F); max(F)] with the global
%   minimum and maximum of the FUNCHEB G on [-1,1].
%
%   [VAL, POS] = MINANDMAX(F) returns also the 2-vector POS where the minima and
%   maxima of F occur.
%
%   If F is complex-valued the absolute values are taken to determine maxima but
%   the resulting values correspond to those of the original function. That is,
%   VAL = feval(F, POS) where [~, POS] = MINANDMAX(abs(F)).
%
% See also MIN, MAX.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( ~isreal(f) )
    realf = real(f);
    imagf = imag(f);
    h = realf.*realf + imagf.*imagf;
    h = simplify(h);
    [~, pos] = minandmax(h);
    vals = feval(f, pos);
    vals = vals(:,1:size(pos,2)+1:end);
    return
end

% Compute derivative:
fp = diff(f);

% Make the Chebyshev grid (used in minandmaxColumn).
xpts = f.chebpts(length(f));

% Check if f is a column FUNCHEB:
sizef2 = size(f, 2);

if ( sizef2 == 1)
    % If it isn't, then things are simple.
    [vals, pos] = minandmaxColumn(f, fp, xpts);
    
else
    % If not, we must loop over the columns.
    
    % Initialise output vectors:
    pos = zeros(2, sizef2);
    vals = zeros(2, sizef2);
    
    % Store coefficients and values:
    values = f.values;
    coeffs = f.coeffs;
    vscale = f.vscale;
    values2 = fp.values;
    coeffs2 = fp.coeffs;
    vscale2 = fp.vscale;
    
    % Loop over columns:
    for k = 1:sizef2    
    
        % Copy data to a single FUNCHEB column:
        f.values = values(:,k);
        f.coeffs = coeffs(:,k);
        f.vscale = vscale(k);
        fp.values = values2(:,k);
        fp.coeffs = coeffs2(:,k);
        fp.vscale = vscale2(k);

        % Find max of this column:
        [vals(:,k), pos(:,k)] = minandmaxColumn(f, fp, xpts);
    
    end

end

end

function [vals, pos] = minandmaxColumn(f, fp, xpts)

   % Initialise output
   pos = [0 ; 0];
   vals = [0 ; 0];

   % Compute turning points
    r = roots(fp);
    r = [-1 ; r ; 1];
    v = feval(f, r);

    % min
    [vals(1), idx] = min(v);
    pos(1) = r(idx);
    % Take the min of the computed min and the function values.
    [vmin, vidx] = min([vals(1) ; f.values]);
    if ( vmin < vals(1) )
        vals(1) = vmin;
        pos(1) = xpts(vidx-1);
    end

    % max
    [vals(2), idx] = max(v);
    pos(2) = r(idx);
    % Take the max of the computed max and the function values.
    [vmax, vidx] = min([vals(2) ; f.values]);
    if ( vmax > vals(2) )
        vals(2) = vmax;
        pos(2) = xpts(vidx-1);
    end

end
