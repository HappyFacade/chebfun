function out = deflationFunH1(Nu, u, r, p, alp)
% DEFLATIONFUN    Wrapper for CHEBMATRIX/DEFLATIONFUN
%
% See also chebmatrix.deflationFun

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Call the CHEBMATRIX method, which is what CHEBOP uses
out = deflationFunH1(Nu, u, chebmatrix(r), p, alp);

end