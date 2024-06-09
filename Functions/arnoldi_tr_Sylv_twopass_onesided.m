function [out]=arnoldi_tr_Sylv_twopass_onesided(eqn,opts)
% function [out]=arnoldi_tr_Sylv_whitening_twopass_onesided(eqn,opts)
% Function that approximately solve the Sylvester matrix equation
%
%    A1 X B2 + A2 X B1 = C1 C2^T
%
% by the one-sided (only left projection is performed), possibly truncated Krylov subspace method presented in
%
% Sketched and Truncated Polynomial Krylov Subspace Methods: Matrix Equations 
% Davide Palitta, Marcel Schweitzer, Valeria Simoncini 
% ArXiv: 2311.16019
%
% Please, acknowledge our work by citing our manuscript whenever you use 
% the software provided here in your research.
%
% INPUT:
% eqn: structure that contains the coefficients defining the equation to be
%      solved.
%      eqn.A1: A1
% function [out]=arnoldi_tr_Sylv_whitening_twopass(eqn,opts)
% Function that approximately solve the Sylvester matrix equation
%
%    A X + X B = C1 C2^T
%
% by the sketched-and-truncated Krylov subspace method presented in
%
% Sketched and Truncated Polynomial Krylov Subspace Methods: Matrix Equations 
% Davide Palitta, Marcel Schweitzer, Valeria Simoncini 
% ArXiv: 2311.16019
%
% Please, acknowledge our work by citing our manuscript whenever you use 
% the software provided here in your research.
%
% INPUT:
% eqn: structure that contains the coefficients defining the equation to be
%      solved.
%      eqn.A1: A1
%      eqn.A2: A2
%      eqn.B1: B1
%      eqn.B2: B2
%      eqn.C1: C1
%      eqn.C2: C2
% opts: structure that contains all the parameters needed by the solver
%       opts.m: maximum number of iterations to be performed
%       opts.ktrunc: truncation parameter for the truncated
%                    orthogonalization step
%       opts.tol: tolerance on the relative residual (sketched) norm
%       opts.check_res: frequency used to solved the projected problem and
%                       check the residual norm; these operations are
%                       carried our every opts.check_res iterations.
%                       Default value: opts.check_res=20.
%
% OUTPUT:
% out: structure containing different quantities
%      out.Z1, out.Z2: low-rank factors of the approximate solution
%                      X \approx Z1 Z2^T
%      out.res_vec: array containing the residual norm history
%      out.it: number of iterations that have been performed
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
% WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR 
% IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

A1=eqn.A1;
A2=eqn.A2;
B1=eqn.B1;
B2=eqn.B2;
C1=eqn.C1;
C2=eqn.C2;
m=opts.m;
check_res=opts.check_res;
ktrunc=opts.ktrunc;
tol=opts.tol;
tol_cg=opts.tol_cg;

[n]=size(A1,1); 

B=B1/B2;
C2=B2'\C2;

s=size(C1,2);
H=zeros((m+1)*s+1,(m+1)*s);
H1=zeros((m+1)*s+1,(m+1)*s);
H2=zeros((m+1)*s+1,(m+1)*s);

D=spdiags(diag(A2),0,n,n);
[wrk,~]=pcg2(A2,C1,0*C1,100,tol_cg,D);

norm_rel=sqrt(trace((wrk'*wrk)*(C2'*C2)));
if s==1
    V(:,1)=wrk/norm(wrk);
    beta=norm(wrk);
else
    [V(:,1:s),beta]=qr(wrk,0);
end


e1=eye(m*s,s);

i=0;

res_vec=zeros(m+1,1);
res=1;
l=1;

ktruncA=ktrunc;

while i<m && res>=tol  %*normC1*normC2 

    i=i+1; 
     
     
    it = 0;

    % we first apply A1 and then solve by A2
    UpA = (A1*V(:,end-s+1:end));
    [UpA,~]=pcg2(A2,UpA,0*UpA,100,tol_cg,D);

    % full orthogonalization for ktruncA=m
    % otherwise we perform a truncated orthogonalization step
    while ( it < 2 )
       it = it+1;
       tt=1;
       for j=max(1,i-ktruncA+1):i
           k1=(j-1)*s+1; k2=j*s;
     
           t = V(1:n,(tt-1)*s+1:tt*s)'*UpA;
           if it==1
                H1(k1:k2,(i-1)*s+1:i*s)=t;
            else
                H2(k1:k2,(i-1)*s+1:i*s)=t;
            end
            UpA=UpA-V(1:n,(tt-1)*s+1:tt*s)*t;
            tt=tt+1;
       end
    end
    
    % normalization
    if s==1
       H2(i+1,i)=norm(UpA);
       UU=UpA/H2(i+1,i);
    else
       [UU, H2(i*s+1:(i+1)*s,(i-1)*s+1:i*s)]=qr(UpA,0);
    end
    H(1:(i+1)*s,1:i*s)=H2(1:(i+1)*s,1:i*s)+H1(1:(i+1)*s,1:i*s);
     
    % move the computed basis blocks if necessary
    if i-ktruncA+1<1
        V(:,(tt-1)*s+1:tt*s)=UU;
    else
        V(:,1:s*(tt-2))=V(:,s+1:s*(tt-1));
        V(:,(tt-2)*s+1:(tt-1)*s)=UU;
    end
        
    if mod(i,check_res)==0
        
        coeff=H(1:i*s,1:i*s);
          
        Y= lyap(coeff,B,-(e1(1:i*s,:)*beta)*C2');

        if ktrunc==m
             resabs=norm(H(i*s+1:(i+1)*s,1:i*s)*Y,'fro');
        else
            resabs=sqrt(i*s)*norm(H(i*s+1:(i+1)*s,1:i*s)*Y,'fro');
        end
        res=resabs/norm_rel;

        res_vec(l)=res;
        l=l+1;
     end


end

out.it=i;
out.res_vec=res_vec(1:l-1);

% Truncation strategy
[uY,sY,vY]=svd(full(Y)); 
sY=diag(sY);
is=sum(abs(sY/sY(1))>1e-12);
Y01 = (uY(:,1:is)*diag(sqrt(sY(1:is))));
Y02 = (vY(:,1:is)*diag(sqrt(sY(1:is))));

% don't perform the two-pass in case of a full method
if ktruncA>=m 
    out.Z1=V(:,1:out.it*s)*Y01;
    out.Z2=Y02;
    return
end

% Two-pass
out.Z1=two_pass_onesided(A1,A2,H1(1:(i+1)*s,1:i*s),H2(1:(i+1)*s,1:i*s),Y01,C1,ktruncA);
out.Z2=Y02;
