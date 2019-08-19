function varargout = chebcoeffs2(f, m, n)
%CHEBCOEFFS2    Bivariate Chebyshev coefficients
%   X = CHEBCOEFFS2(F) returns the matrix of bivariate coefficients such that
%       F = sum_{i=0}^{m-1} ( sum_{j=0}^{n-1} X(i+1,j+1) T_i(y) T_j(x) ),
%   where m and n are the degree of F. 
%
%   [A, D, B] = CHEBCOEFFS2( f ) returns the same coefficients keeping them in
%   low rank form, i.e., X = A * D * B'.
%
%   CHEBCOEFFS2(F, M, N) is the same as CHEBCOEFFS2(F) except the matrix of
%   coefficients is returned as an MxN matrix. CHEBCOEFFS2(F, M) is the
%   same as CHEBCOEFFS2(F, M, M).
%
% See also PLOTCOEFFS2, CHEBCOEFFS.

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) )
    varargout = { [ ] }; 
    return
end

if ( iszero(f) ) 
    varargout = { 0 } ; 
    return
end

% Get the low rank representation for f:
[cols, d, rows] = cdr(f);

% Get the coeffs of the rows and the columns:
if ( nargin == 1 )
    cols_coeffs = chebcoeffs(cols);
    rows_coeffs = chebcoeffs(rows);
elseif ( nargin == 2 )
    cols_coeffs = chebcoeffs(cols, m);
    rows_coeffs = chebcoeffs(rows, m);
elseif ( nargin == 3 )
    cols_coeffs = chebcoeffs(cols, m);
    rows_coeffs = chebcoeffs(rows, n);
else 
    error('CHEBFUN2:CHEBCOEFFS2:ARGIN', 'Too many input arguments');
end

if ( nargout <= 1 )
    % Return the matrix of coefficients
    varargout = { cols_coeffs * d * rows_coeffs.' }; 
elseif ( nargout <= 3 )
    varargout = { cols_coeffs, d, rows_coeffs };
else
    % Two output variables are not allowed.
    error('CHEBFUN:CHEBFUN2:chebcoeffs2:outputs', ...
        'Incorrect number of outputs.');
end

end