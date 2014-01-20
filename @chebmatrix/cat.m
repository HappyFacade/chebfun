function C = cat(dim, varargin)
%CAT    Concatenation of chebmatrices.
%   C = CAT(dim,A,B,...) concatenates chebmatrices along the indicated
%   dimension.
%
%   See also CHEBMATRIX.HORZCAT, CHEBMATRIX.VERTCAT.

%   Copyright 2013 by The University of Oxford and The Chebfun Developers.
%   See http://www.chebfun.org for Chebfun information.

% Any singleton operator blocks must be encased in cells to match up with
% chebmatrix blocks.

if ( dim == 1 )
    C = vertcat(varargin{:});
elseif ( dim == 2 )
    C = horzcat(varargin{:});
else
    error('Must concatenate along dimension 1 or 2.')
end  

end
