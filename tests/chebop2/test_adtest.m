function pass = test_adtest( prefs ) 
% Test that the AD machinery is work correctly. 
% Alex Townsend, August 2013. 

if ( nargin < 1 ) 
    prefs = chebfunpref(); 
end 
tol = prefs.cheb2Prefs.eps; 

% simple case. 
N = chebop2(@(u) diffx(u,2) + diffy(u,2));
pass(1) = ( norm( cell2mat(N.coeffs) - [0 0 1; 0 0 0; 1 0 0]) < tol );

% Helmoltz 
N = chebop2(@(u) diffx(u,2) + diffy(u,2) + pi*u);
pass(2) = ( norm( cell2mat(N.coeffs) - [pi 0 1; 0 0 0; 1 0 0]) < tol );

% Higher order 
N = chebop2(@(u) diffx(u,3) + diffx(diffy(u,1),2) + diffy(u,2) + pi*u);
pass(3) = ( norm( cell2mat(N.coeffs) - [pi 0 0 1; 0 0 1 0; 1 0 0 0]) < tol );

% Variable coefficients on zero derivatives.
x = chebfun2(@(x,y) x); 
y = chebfun2(@(x,y) y);
N = chebop2(@(x,y,u) diff(u,1,2) + diff(u,2,1)  + x.*u);
CC = N.coeffs;
pass(4) = ( norm( CC{1,1} - x) < 2*tol ); 
pass(5) = ( norm( CC{1,2} - 0*x - 1 ) < 10*tol ); 
pass(6) = ( norm( CC{3,1} - 0*x - 1 ) < 10*tol ); 

% % Variable coefficients on derivatives
x = chebfun2(@(x,y) x); 
y = chebfun2(@(x,y) y);
N = chebop2(@(x,y,u) diff(u,1,2) + x.*diff(u,2,1));
CC = N.coeffs;
pass(7) = ( norm( CC{3,1} - x) < 2*tol ); 
pass(8) = ( norm( CC{1,2} - 0*x - 1 ) < 10*tol ); 
pass(9) = ( norm( CC{1,1} - 0*x ) < 2*tol ); 

end