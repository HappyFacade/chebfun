function r = roots( varargin )
%ROOTS      Zero contours of a SPHEREFUN
%   R = ROOTS(F), returns the zero contours of F as a quasimatrix of 
%   array-valued chebfuns. Each column of R is one zero contour. This 
%   command only finds contours when there is a change of sign and it may 
%   group intersecting contours in a non-optimal way. Contours are computed
%   to, roughly, four digits of precision. In particular, this command 
%   cannot reliably compute isolated real roots of F. 
%
%   In the special case when F is of length 1 then the zero contours are 
%   found to full precision.
%  
% See also CHEBFUN2/ROOTS.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Check for empty:
if ( isempty( varargin{ 1 }  ) )
    r = []; 
    return
end

% rts = roots@separableApprox(chebfun2(varargin{1},[-pi-0.1 pi+0.1 0 pi]));
rts = roots@separableApprox(chebfun2(varargin{1},1.5*[-pi pi -pi pi]));

% Now make into a collection of array-valued chebfuns ready for plotting on
% the sphere.  Make sure we have enough points to sample the contours.
x = chebpts(max(length(rts),17) + 1);

vals = feval(rts, x);
r = cell(size(vals,2), 1);

% Go through each component and make it an array-valued chebfun: 
cnt = 1;
for k = 1:size(vals, 2)
    comp = feval(rts(:, k), x); 
    lam = real(comp);
    th = imag(comp);
    idlam = lam <= pi & lam >= -pi;
    idth = th <= pi & th >= 0;
    id = idlam & idth;
    if any(id)
        AX = cos(lam).*sin(th);
        AY = sin(lam).*sin(th);
        AZ = cos(th);
        r{cnt} = chebfun([AX, AY, AZ]);
        cnt = cnt + 1;
    end
end
r = r{1:cnt-1};

end