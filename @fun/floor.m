function f = floor(f)
%FLOOR   Pointwise floor function of a FUN.
%   G = FLOOR(F) returns the FUN G such that G(X) = FLOOR(F(x)) for each x in
%   F.domain. 
%
%   If F is complex, then the G = FLOOR(REAL(F)) + 1i*FLOOR(IMAG(F)).
%
%   Note that FLOOR() assumes the output G(X) is a constant. If it is not, then
%   garbage is returned with no warning.
%
% See also CEIL, ROUND, FIX.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% FLOOR() the ONEFUN:
f.onefun = floor(f.onefun);

end