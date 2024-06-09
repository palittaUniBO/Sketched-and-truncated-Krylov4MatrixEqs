function [out]=arnoldi_tr_Sylv_twopass(eqn,opts)
%function [Z1,Z2,res_vec,it,time,res_actualres,VV,WW]=arnoldi_tr_Sylv_twopass(A,B,C1,C2,m,ktrunc,tol,check_res,check_actualres,store)
% function [out]=arnoldi_tr_Sylv_whitening_twopass(eqn,opts)
% Function that approximately solve the Sylvester matrix equation
%
%    A X + X B = C1 C2^T
%
% by the (possibly truncated) polynomial Krylov subspace method. 
% The truncated method is presented in
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
%      eqn.A: A
%      eqn.B: B
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
%       opts.store: boolean value ruling whether the whole bases are stored
%                   or not. Default value: opts.store=0.
%
% OUTPUT:
% out: structure containing different quantities
%      out.Z1, out.Z2: low-rank factors of the approximate solution
%                      X \approx Z1 Z2^T
%      out.res_vec: array containing the residual norm history
%      out.it: number of iterations that have been performed
%      out.condA, out.condB: history of the condition numbers of the 
%                            constructed with A (out.condA) and B^T (out.condB). 
%                            These arrays are defined iff opts.store==1. 
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
% WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR 
% IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

A=eqn.A;
B=eqn.B;
C1=eqn.C1;
C2=eqn.C2;
m=opts.m;
ktrunc=opts.ktrunc;
tol=opts.tol;

[n]=size(A,1); 
nB=size(B,1);

s=size(C1,2);

% if check_res is not set, we set it equal to 20
if isfield(opts,'check_res')
    check_res=opts.check_res;
else
    check_res=20;
end



H=zeros((m+1)*s,m*s);
H1=zeros((m+1)*s,m*s);
H2=zeros((m+1)*s,m*s);

T=zeros((m+1)*s,m*s);
T1=zeros((m+1)*s,m*s);
T2=zeros((m+1)*s,m*s);

if s==1
    normC1 = norm(C1);
    V(:,1) = C1/normC1;
    beta1=normC1;
    
    W(:,1)=C2/norm(C2);
    normC2 = norm(C2);
    beta2=normC2;
else
    [V(:,1:s),beta1]=qr(C1,0);
    [W(:,1:s),beta2]=qr(C2,0);
    
end


e1=eye(m*s,s);

i=0;

res_vec=zeros(m+1,1);
res=1;
l=1;

if length(ktrunc)==1
    ktruncA=ktrunc;
    ktruncB=ktrunc;
else
    ktruncA=ktrunc(1);
    ktruncB=ktrunc(2);
end

% if store is not set, we set it equal to 0
if isfield(opts,'store')
    store=opts.store;
    out.condA=zeros(m,1);
    out.condA=zeros(m,1);  
    wholeV=zeros(n,(m+1)*s);
    wholeW=wholeV;
    wholeV(:,1:s)=V(:,1:s);
    wholeW(:,1:s)=W(:,1:s);

else
    store=0;
    out.condA=0;
    out.condB=0;
end



while i<m && res>=tol*norm(beta1*beta2','fro') 

     i=i+1; 

    
    it = 0; 

    %Sequence in A
    UpA = A*V(:,end-s+1:end);


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
     


    %Sequence in B
    UpB = B'*W(:,end-s+1:end);

    it = 0; 

    % full orthogonalization for ktruncB=m
    % otherwise we perform a truncated orthogonalization step

     while ( it < 2 )
        it = it+1;
        tt=1;
        for j=max(1,i-ktruncB+1):i
            k1=(j-1)*s+1; k2=j*s;
     
            t = W(1:nB,(tt-1)*s+1:tt*s)'*UpB;
            if it==1
                T1(k1:k2,(i-1)*s+1:i*s)=t;
            else
                T2(k1:k2,(i-1)*s+1:i*s)=t;
            end
            UpB=UpB-W(1:nB,(tt-1)*s+1:tt*s)*t;
            tt=tt+1;
        end
     end

     % normalization
     if s==1
        T2(i+1,i)=norm(UpB);
        WW=UpB/T2(i+1,i);
     else
        [WW, T2(i*s+1:(i+1)*s,(i-1)*s+1:i*s)]=qr(UpB,0);
     end
     T(1:(i+1)*s,1:i*s)=T2(1:(i+1)*s,1:i*s)+T1(1:(i+1)*s,1:i*s);


     % move the computed basis blocks if necessary
     if i-ktruncB+1<1
          W(:,(tt-1)*s+1:tt*s)=WW;
     else
          W(:,1:s*(tt-2))=W(:,s+1:s*(tt-1));
          W(:,(tt-2)*s+1:(tt-1)*s)=WW;
     end
    
     
     if mod(i,check_res)==0
  
    
          %solve the projected problem
          Y = lyap(H(1:i*s,1:i*s),T(1:i*s,1:i*s)', -(e1(1:i*s,:)*beta1)*(e1(1:i*s,:)*beta2)');
          

          % residual norm
          if ktruncA < m || ktruncB < m
              % for the truncated method this is an estimate from above of
              % the actual residual norm
              res=sqrt(i*s)*sqrt(norm(H(i*s+1:(1+i)*s,1:i*s)*Y,'fro')^2+norm(Y*T(i*s+1:(1+i)*s,1:i*s)','fro')^2);
          else
              res=sqrt(norm(H(i*s+1:(1+i)*s,1:i*s)*Y,'fro')^2+norm(Y*T(i*s+1:(1+i)*s,1:i*s)','fro')^2);
          end
          res_vec(l)=res;
          
              

          l=l+1;
     end

     if store
         wholeV(:,i*s+1:(i+1)*s)=UU;
         out.condA(i)=cond(wholeV(:,1:i*s));
         
         wholeW(:,i*s+1:(i+1)*s)=WW;
             out.condB(i)=cond(wholeW(:,1:i*s));
     end
end

out.it=i;
out.res_vec=res_vec(1:l-1);

% Truncation strategy
[uY,sY,vY]=svd(full(Y)); 
sY=diag(sY);
is=sum(abs(sY/sY(1))>1e-12);
Y01 = uY(:,1:is)*diag(sqrt(sY(1:is)));
Y02= vY(:,1:is)*diag(sqrt(sY(1:is)));

% don't perform the two-pass in case of a full method
if (ktruncA==m && ktruncB==m) 
    out.Z1=V(:,1:out.it*s)*Y01;
    out.Z2=W(:,1:out.it*s)*Y02;
elseif store
    out.Z1=wholeV(:,1:out.it*s)*Y01;
    out.Z2=wholeW(:,1:out.it*s)*Y02;
    return
end

% Two-pass

[out.Z1,out.Z2]=two_pass(A,B,H1(1:(i+1)*s,1:i*s),H2(1:(i+1)*s,1:i*s),T1(1:(i+1)*s,1:i*s),...
    T2(1:(i+1)*s,1:i*s),Y01,Y02,C1,C2,ktruncA,ktruncB);

% %sequence in A
% V=[];
% Z1=0;
% 
% % sequence in B
% W=[];
% Z2=0;
% 
% if s==1
%     V(1:n,1:s)=C1/normC1;
%     W(1:nB,1:s)=C2/normC2;
% else
%     [V(1:n,1:s),~]=qr(C1,0);
%     [W(1:nB,1:s),~]=qr(C2,0);
% end
% 
% 
% 
% for j=1:it-1
% 
%     % sequence in A and B  
%     if j<itA
%         kminA=max(1,j-ktruncA+1);kmax=j;
%         tmin=1;tmaxA=kmax-kminA+1;
%         UpA = A*V(1:n,end-s+1:end);
%         UpA=UpA-V(1:n,(tmin-1)*s+1:tmaxA*s)*H((kminA-1)*s+1:kmax*s,(j-1)*s+1:j*s);
% 
%         ttA=tmaxA+1;
% 
%         % normalize
%         UpA=UpA/H(j*s+1:(j+1)*s,(j-1)*s+1:j*s);
% 
% %         itt=0;
% % 
% %         while ( itt < 2 )
% %         itt = itt+1;
% %             tt=1;
% %             for jj=max(1,j-ktruncA+1):j
% %                 k1=(jj-1)*s+1; k2=jj*s;
% %                 t = V(1:n,(tt-1)*s+1:tt*s)'*UpA;
% %                 H_nosket(k1:k2,(j-1)*s+1:j*s)=H_nosketched(k1:k2,(j-1)*s+1:j*s)+t;
% %                 UpA=UpA-V(1:n,(tt-1)*s+1:tt*s)*t;
% %                 tt=tt+1;
% %             end
% %         end
% %         if s==1
% %             H_nosketched(j+1,j)=norm(UpA);
% %             UpA=UpA/H_nosketched(j+1,j);
% %         else
% %             [UpA, H_nosketched(j*s+1:(j+1)*s,(j-1)*s+1:j*s)]=qr(UpA,0);
% %         end
% 
%         
%         % move the computed basis blocks if necessary
%         if j-ktruncA+1<1
%             V(:,(ttA-1)*s+1:ttA*s)=UpA;
%         else
%             V(:,1:s*(ttA-2))=V(:,s+1:s*(ttA-1));
%             V(:,(ttA-2)*s+1:(ttA-1)*s)=UpA;
%         end
% 
%         % compute ktrunc-rank components of the solution
%         if j-ktruncA+1==0
%             Z1=Z1+V*Y01(1:(j+1)*s,:);
%         elseif mod(j+1,ktruncA)==0
%             Z1=Z1+V*Y01((j-ktruncA+1)*s+1:(j+1)*s,:);
%             %Z1=Z1+V*Y01(j-ktruncA+2:j+1,:);
%         end
%   
%     end
% 
%     
%     if j<itB
%         kminB=max(1,j-ktruncB+1);kmax=j;
%         tmin=1;tmaxB=kmax-kminB+1;
%         UpB = B'*W(1:nB,end-s+1:end);
%         UpB=UpB-W(1:nB,(tmin-1)*s+1:tmaxB*s)*T((kminB-1)*s+1:kmax*s,(j-1)*s+1:j*s);
%         ttB=tmaxB+1;
% 
%         UpB=UpB/T(j*s+1:(j+1)*s,(j-1)*s+1:j*s);
% %         itt = 0; 
% % 
% %         while ( itt < 2 )
% %             itt = itt+1;
% %             tt=1;
% %             for jj=max(1,j-ktruncB+1):j
% %                 k1=(jj-1)*s+1; k2=jj*s;
% %      
% %                 t = W(1:nB,(tt-1)*s+1:tt*s)'*UpB;
% %                 T_nosketched(k1:k2,(j-1)*s+1:j*s)=T_nosketched(k1:k2,(j-1)*s+1:j*s)+t;
% %                 UpB=UpB-W(1:nB,(tt-1)*s+1:tt*s)*t;
% %                 tt=tt+1;
% %             end
% %         end
% %         if s==1
% %             T_nosketched(i+1,i)=norm(UpB);
% %             UpB=UpB/T_nosketched(i+1,i);
% %         else
% %             [UpB, T_nosketched(i*s+1:(i+1)*s,(i-1)*s+1:i*s)]=qr(UpB,0);
% %         end
% 
% 
%         % move the computed basis blocks if necessary
%         if j-ktruncB+1<1
%             W(:,(ttB-1)*s+1:ttB*s)=UpB;
%         else
%             W(:,1:s*(ttB-2))=W(:,s+1:s*(ttB-1));
%             W(:,(ttB-2)*s+1:(ttB-1)*s)=UpB;
%         end
%     
%         % compute ktrunc-rank components of the solution
%         if j-ktruncB+1==0
%             Z2=Z2+W*Y02(1:(j+1)*s,:);
%         elseif mod(j+1,ktruncB)==0
%             Z2=Z2+W*Y02((j-ktruncB+1)*s+1:(j+1)*s,:);
%             %Z2=Z2+W*Y02(j-ktruncB+2:j+1,:);
%         end
%     end
% 
% end
% 
% k=mod(itA,ktruncA);
% if k~=0
%     Z1=Z1+V(:,(end-k)*s+1:end)*Y01((end-k)*s+1:end,:); 
%     %Z1=Z1+V(:,end-k+1:end)*Y01(end-k+1:end,:);
% end    
% k=mod(itB,ktruncB);
% if k~=0
%     Z2=Z2+W(:,(end-k)*s+1:end)*Y02((end-k)*s+1:end,:);
%     %Z2=Z2+W(:,end-k+1:end)*Y02(end-k+1:end,:);
% end    
% % further truncation
% [Q1,R1]=qr(Z1,0);
% [Q2,R2]=qr(Z2,0);
% 
% [uY,sY,vY]=svd(full(R1*R2')); 
% sY=diag(sY);
% is=sum(abs(sY/sY(1))>eps)
% 
% Z1 = Q1*(uY(:,1:is)*diag(sqrt(sY(1:is))));
% Z2 = Q2*(vY(:,1:is)*diag(sqrt(sY(1:is))));
% 
% time.tot=toc(tt_tot);
% time.basis=time_basis;
% time.proj=time_proj;
% time.it=time_it;
% 
