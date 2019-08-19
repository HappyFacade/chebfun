function pass = test_chebcoeffs2( pref ) 
% Test chebpoly2 

if ( nargin == 0) 
    pref = chebfunpref; 
end

tol = 1000*pref.cheb2Prefs.chebfun2eps;
j = 1; 

% Rank-2 function
n = 10;
m = 8;
Tn = chebpoly(n);
Tm = chebpoly(m);
f = Tn * Tn' + Tm * Tn'; 
X = chebcoeffs2( f );
Exact = zeros(n+1); Exact(n+1,n+1) = 1; Exact(m+1,n+1) = 1; 
pass(1) = norm( X - Exact ) < tol; 

X = chebcoeffs2( f, 20, 30);
Exact = zeros(20, 30); 
Exact(n+1,n+1) = 1; Exact(m+1,n+1) = 1; 
pass(2) = norm( X - Exact ) < tol; 

%f = Tm * Tn';
% check inverses
Exact = Exact(1:n+1,1:n+1); 
Z = chebfun2.vals2coeffs( chebpolyval2( f ) );
pass(3) = norm( Z - Exact ) < tol; 

f = chebfun2( @(x,y) x + y );
pass(4) = isequal( size(coeffs2(f,2,1)), [2 1] ); 

end