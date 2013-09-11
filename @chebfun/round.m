function f = round(f)
%ROUND    Round a CHEBFUN pointwise to nearest integer.
%   G = ROUND(F) returns the CHEBFUN G such that G(x) = ROUND(F(x)) for each x
%   in F.domain.
%
%   If F is complex, then the G = ROUND(REAL(F)) + 1i*ROUND(IMAG(F)).
%
% See also FIX, FLOOR, CEIL.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Deal with the trivial empty case:
if ( isempty(f) )
    return
end

% Deal with unbounded functions:
if ( ~isfinite(f) )
    error('CHEBFUN:floor:inf', ...
        'Floor is not defined for functions which diverge to infinity.');
end

% Deal with complex-valued functions:
if ( ~isreal(f) )
    if ( isreal(1i*f) )
        f = 1i*round(imag(f));
    else
        f = round(real(f)) + 1i*round(imag(f));
    end
    return
end

% Find all the integer crossings for f:
[minf, maxf] = minandmax(f); % [TODO]: Only need a good bound?
range = floor([minf, maxf]);
for k = (range(1)+1):range(2)
    f = addBreaksAtRoots(f - k + .5) + k - .5;
end

% Loop over the FUNs:
for k = 1:numel(f.funs)
    f.funs{k} = round(f.funs{k});
end

% Floor the impulses:
f.impulses = round(f.impulses(:,:,1));

% SImplify the result:
f = merge(f);

end

