function [A, B, C, D, Xtrue] = conv_diff_3d_sep(n,r,exact_sol)
% CONV_DIFF_3D_SEP creates matrices A and B of dimension n^3 and random data C, D of
% rank r for setting up the Sylvetser equation A X + X B + C D' = 0.  Data
% are saved in an eponymous .mat file with n^3 and r specified. By default,
% n^3 =  32768, r = 2.
% Example 5-6 from Palitta and Simoncini (2016) (We use zero Dirichlet bc
% here)

%%
if nargin == 0
    n = 32;
    r = 2;
    exact_sol = true;
elseif nargin == 1
    r = 2;
    exact_sol = true;
elseif nargin == 2
    exact_sol = true;
end

% epsilon: viscosity parameter. For epsilon >> 0 the problem becomes harder
% as we have dominant convection
epsilon = 5e-3;
h = 1 / (n - 1);

% Negative Laplacian
T = -epsilon / h^2 *...
    spdiags([-ones(n,1), 2 * ones(n,1), -ones(n,1)], -1:1, n, n);
% negative first derivative
N = -1 / (2 * h) *...
    spdiags([-ones(n, 1), zeros(n, 1), ones(n, 1)], -1:1, n, n);

x = linspace(0, 1, n);

phi1 = x.*sin(x); 
psi2 = x.*cos(x);
upsilon3=exp(x.^2-1);  
%psi2 = 0;

PHI1 = spdiags(phi1',0,n,n);    
%PHI2 = spdiags(phi2',0,n,n);
%PSI1 = spdiags(psi1',0,n,n);    
PSI2 = spdiags(psi2',0,n,n);                             
UPSILON3=spdiags(upsilon3',0,n,n);

B1 = PHI1 * N;
B2= PSI2*N;
B3= UPSILON3* N';

I = speye(n);

A=kron(kron(I,T+B2'),I)+kron(T+B3,kron(I,I))+kron(kron(I,I),T+B1);


phi1=1-x.^2;
psi1=x;
upsilon1=x;
%phi2=x;
%psi2=-2*(1-x.^2);
upsilon3=exp(x);

PHI1=spdiags(phi1',0,n,n);    
PSI1=spdiags(psi1',0,n,n);    
UPSILON1=spdiags(upsilon1',0,n,n);
UPSILON3=spdiags(upsilon3',0,n,n);

I=speye(n);

B1=PHI1*N;
B3= UPSILON3* N';

B =kron(kron(I,T),I)+kron(T+B3,kron(I,I))+kron(kron(I,I),T)+kron(UPSILON1,kron(PSI1,B1));

%A=T+B1;
%B=T;

% Random RHS
rng('default')
C = randn(n^3,r);
D = randn(n^3,r);
normCD=sqrt(trace((C'*C)*(D'*D)));
C = C/sqrt(normCD);
D = D/sqrt(normCD);

if exact_sol
    Xtrue = lyap(A,B,-C*D');
end

end
