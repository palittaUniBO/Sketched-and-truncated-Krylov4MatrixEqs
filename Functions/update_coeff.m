function Hnew=update_coeff(Hold,R,H,h)
% function Hnew=update_coeff(Hold,R,H,h)
% function that cheaply computes the update of (part of) the projected
% coefficient matrices for the sketched-and-truncated Krylov subspace
% method presented in Section 5 in
%  
% Sketched and Truncated Polynomial Krylov Subspace Methods: Matrix Equations 
% Davide Palitta, Marcel Schweitzer, Valeria Simoncini 
% ArXiv: 2311.16019
%
% Please, acknowledge our work by citing our manuscript whenever you use 
% the software provided here in your research.

s=size(h,1);
ds=size(Hold,1);
Hnew=zeros(ds+s);

tau_d=R(ds-s+1:ds,ds-s+1:ds);
tau_d1=R(ds+1:ds+s,ds+1:ds+s);
T_H=R(1:ds,ds+1:ds+s);
hh=H(ds+1:ds+s,:);
ed=zeros(ds,s);
ed(ds-s+s:ds,:)=eye(s);

Hnew(1:ds,1:ds)=Hold+T_H*(h/tau_d)*ed';
Hnew(1:ds,ds+1:ds+s)=-Hold(1:ds,1:ds)*(T_H/tau_d1)+R(1:ds,1:ds)*(H(1:ds,:)/tau_d1)+T_H*(hh/tau_d1)-T_H*(h/tau_d*(ed'*T_H)/tau_d1);
Hnew(ds+1:ds+s,ds-s+1:ds)=tau_d1*h/tau_d;
Hnew(ds+1:ds+s,ds+1:ds+s)=-(tau_d1*h/tau_d)*T_H(ds-s+1:ds,:)/tau_d1+tau_d1*hh/tau_d1;


