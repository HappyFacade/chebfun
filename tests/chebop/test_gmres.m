function pass = test_gmres() 
% Test preconditioned conjugate gradient method. 

% Tolerance: 
tol = 5e-13; 
n = 50;

%% Example 1: (-u_xx=f, bc=0, sum(f)==0)
a = @(x) 1;
c = @(x) 1;
L = @(x,u) -diff(a(x).*diff(u)) + c(x).*u;

% Chebop solve:
N = chebop( L );
N.bc = 0; 
f = chebfun(@(x) 1-3*x.^2);

u = N\f;
% Chebfun-ready gmres:

v = gmres2(L, chebfun(f));

% Errors:
pass(1) = ( norm( u - v ) < tol );

%% Example 2: (-(a(x)*u_x)_x=f, bc=0, sum(f)==0)
a = @(x) 2+cos(pi*x);
c = @(x) 0*x;
L = @(x,u) -diff(a(x).*diff(u)) + c(x).*u;

% Chebop solve:
N = chebop( L );
N.bc = 0; 
f = chebfun(@(x) 1-3*x.^2);
u = N \ f; 

% Chebfun-ready gmres:
v = gmres2(L, chebfun(f, n));

% Errors:
pass(2) = ( norm( u - v ) < tol );

%% Example 3: (-(a(x)*u_x)_x + c(x)*u=f, bc=0, sum(f)==0)

a = @(x) 2+cos(pi*x);
c = @(x) 1+10*x.^2;
L = @(x,u) -diff(a(x).*diff(u)) + c(x).*u;

% Chebop solve:
N = chebop( L );
N.bc = 0; 
f = chebfun(@(x) 1-3*x.^2);
u = N \ f; 
% Chebfun-ready gmres:
v = gmres2(L, chebfun(f, n));

% Errors:
pass(3) = ( norm( u - v ) < tol );

%% Example 4: (-(a(x)*u_x)_x + c(x)*u=f, bc=0, sum(f)~=0)
a = @(x) 2+cos(pi*abs((x-0.5)));
c = @(x) 1+10*x.^2;
L = @(x,u) -diff(a(x).*diff(u)) + c(x).*u;

% Chebop solve:
N = chebop( L );
N.bc = 0; 
%f = chebfun(@(x) 1-2*x.^2);
f = chebfun(@(x) abs(x) - 0.5 ,'splitting','on');
u = N \ f; 
% Chebfun-ready gmres:
v = gmres2(L, f);
pass(4) = ( norm(u - v) < tol );

%% Example 5: (-(a(x)*u_x)_x + c(x)*u=f, bc~=0, sum(f)~=0)
a = @(x) 2+cos(pi*x);
c = @(x) 1+10*x.^2;
L = @(x,u) -diff(a(x).*diff(u)) + c(x).*u;

% Chebop solve:
N = chebop( L );
N.lbc = 1; 
N.rbc = -1; 
f = chebfun(@(x) 1-2*x.^2);
u = N \ f; 

% Chebfun-ready gmres:
v = gmres2(N, chebfun(f));

pass(5) = ( norm(u - v) < tol );

%% Example 6: (-(a(x)*u_x)_x + c(x)*u=f, bc~=0, sum(f)~=0)
a = @(x) 2+cos(pi*x);
b = @(x) 1;
c = @(x) -10*x.^2;
L = @(x,u) -diff(a(x).*diff(u)) + b(x).*diff(u,1) + c(x).*u;

% Chebop solve:
N = chebop( L );
N.lbc = 1; 
N.rbc = -1; 
f = chebfun(@(x) 1-2*x.^2);
u = N \ f; 

% Chebfun-ready gmres:
v = gmres2(N, chebfun(f));

pass(6) = ( norm(u - v) < tol );


if ( all(pass) == 1 )
    pass = 1;
end

end