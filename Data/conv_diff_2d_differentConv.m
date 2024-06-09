function [A, B, C, D, Xtrue] = conv_diff_2d_differentConv(n,r,epsilon,exact_sol)
% CONV_DIFF_2D_constconv creates matrices A and B of dimension n^2 and random data C, D of
% rank r for setting up the Sylvetser equation A X + X B + C D' = 0.  Data
% are saved in an eponymous .mat file with n^3 and r specified. By default,
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

% epsilon: visosity parameter. For epsilon >> 0 the problem becomes harder
% as we have dominant convction

h = 1 / (n - 1);

% Negative Laplacian
T = -epsilon / h^2 *...
    spdiags([-ones(n,1), 2 * ones(n,1), -ones(n,1)], -1:1, n, n);
% negative first derivative
N = -1 / (2 * h) *...
    spdiags([-ones(n, 1), zeros(n, 1), ones(n, 1)], -1:1, n, n);

x = linspace(0, 1, n);
% 
% phi1 = x.*sin(x); 
% psi2 = x.*cos(x);
% upsilon3=exp(x.^2-1);  
% %psi2 = 0;
% 
% PHI1 = spdiags(phi1',0,n,n);    
% %PHI2 = spdiags(phi2',0,n,n);
% %PSI1 = spdiags(psi1',0,n,n);    
% PSI2 = spdiags(psi2',0,n,n);                             
% UPSILON3=spdiags(upsilon3',0,n,n);
%B1 = PHI1 * N;
%B2= PSI2*N;
%B3= UPSILON3* N';

B1 = N;
B2= N;

I = speye(n);

A=kron(I,T+B2')+kron(T+B1,I);

%epsilon = 1e-
%epsilon = 1e-2

T = -epsilon / h^2 *...
    spdiags([-ones(n,1), 2 * ones(n,1), -ones(n,1)], -1:1, n, n);

phi1=3*(1-x.^2);
psi1=x;
%upsilon1=x;
phi2=x;
psi2=-2*(1-x.^2);
%upsilon3=exp(x);

PHI1=spdiags(phi1',0,n,n);    
PSI1=spdiags(psi1',0,n,n);    
PHI2=spdiags(phi2',0,n,n);    
PSI2=spdiags(psi2',0,n,n);   
I=speye(n);

B1=PHI1*N;
B2= PSI2* N;

B =kron(I,T)+kron(T,I)+kron(B2,PHI2)+kron(PSI1,B1');



%A=T+B1;
%B=T;

% Random RHS
rng('default')
C = randn(n^2,r);
%C=A*(A*(A*C));
D = randn(n^2,r);
%D=B*(B*(B*D));
normCD=sqrt(trace((C'*C)*(D'*D)));
C = C/sqrt(normCD);
D = D/sqrt(normCD);

if exact_sol
    Xtrue = lyap(A,B,-C*D');
end

end
