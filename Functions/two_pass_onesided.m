function [Z1]=two_pass_onesided(A1,A2,H,H1,Y01,C1,ktruncA)
% function [Z1]=two_pass_onesided(A1,A2,H,H1,Y01,C1,ktruncA)
% function that implements the two-pass strategy to recover the solution
% factor Z1.
%
% INPUT: 
% A1, A2, C1: coefficient matrices definining the Sylvester equation to
%               be solved
% Y01: low-rank factor of the solution to the last projected equation
% H, H1: Banded upper Hessenberg matrices containing the orthogonalization 
%        coefficients computed during the first Arnoldi pass in A
% ktruncA: truncation parameter adopted in the first Arnoldi pass
%
% OUTPUT:
% Z1: left low-rank factor of the approximate solution X \approx Z1 Z2^T
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
% WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR 
% IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

[n,s]=size(C1);
d=size(H,2)/s;
tol=1e-8;

D=spdiags(diag(A2),0,n,n);
[wrk,~]=pcg2(A2,C1,0*C1,100,tol,D);

if s==1
    V(1:n,1:s)=wrk/norm(wrk,'fro');
else
    [V(1:n,1:s),~]=qr(wrk,0);
end

l=size(Y01,2);
Z1=zeros(n,l);

for j=1:d-1
    
    % sequence in A: first apply A1 then solve by A2
    UpA = A1*V(:,end-s+1:end);
    [UpA,~]=pcg2(A2,UpA,0*UpA,100,tol,D);
       


    % we compute the same exact operations performed in the first Arnoldi
    % loop. We notice that this leads to a more stable procedure,
    % especially in the block case (size(C1,2)>1)
    tt=1;
    for jj=max(1,j-ktruncA+1):j
        k1=(jj-1)*s+1; k2=jj*s;
        t = H(k1:k2,(j-1)*s+1:j*s);
        UpA=UpA-V(1:n,(tt-1)*s+1:tt*s)*t;
        tt=tt+1;
    end

    tt=1;
    for jj=max(1,j-ktruncA+1):j
        k1=(jj-1)*s+1; k2=jj*s;
        t = H1(k1:k2,(j-1)*s+1:j*s);
        UpA=UpA-V(1:n,(tt-1)*s+1:tt*s)*t;
        tt=tt+1;
    end

    % normalize 
    if s==1
       UpA=UpA/norm(UpA);
    else
       [UpA, ~]=qr(UpA,0);
    end
       
    % move the computed basis blocks if necessary
    if j-ktruncA+1<1
        V(:,(tt-1)*s+1:tt*s)=UpA;
    else
        V(:,1:s*(tt-2))=V(:,s+1:s*(tt-1));
        V(:,(tt-2)*s+1:(tt-1)*s)=UpA;
    end

    % compute ktrunc-rank components of the solution
    if j-ktruncA+1==0
        Z1=Z1+V*Y01(1:(j+1)*s,:);
    elseif mod(j+1,ktruncA)==0
        Z1=Z1+V*Y01((j-ktruncA+1)*s+1:(j+1)*s,:);
    end


end

k=mod(d,ktruncA);
if k~=0
    Z1=Z1+V(:,end-k*s+1:end)*Y01(end-k*s+1:end,:);
end    



  
    
