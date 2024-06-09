function [Z1,Z2]=two_pass(A,B,H,H1,T,T1,Y01,Y02,C1,C2,ktruncA,ktruncB)
% function [Z1,Z2]=two_pass(A,B,H,H1,T,T1,Y01,Y02,C1,C2,ktruncA,ktruncB)
% function that implements the two-pass strategy to recover the solution
% factors Z1 and Z2.
%
% INPUT: 
% A, B, C1, C2: coefficient matrices definining the Sylvester equation to
%               be solved
% Y01,Y02: low-rank factors of the solution to the last projected equation
% H, H1: Banded upper Hessenberg matrices containing the orthogonalization 
%        coefficients computed during the first Arnoldi pass in A
% T, T1: Banded upper Hessenberg matrices containing the orthogonalization 
%        coefficients computed during the first Arnoldi pass in B^T
% ktruncA, ktruncB: truncation parameter adopted in the first Arnoldi pass
%                   for A and B^T, respectively
%
% OUTPUT:
% Z1,Z2: low-rank factors of the approximate solution X \approx Z1 Z2^T
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
% WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR 
% IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

[n,s]=size(C1);
nB=size(C2,1);
d=size(H,2)/s;

if s==1
    V(1:n,1:s)=C1/norm(C1,'fro');
    W(1:nB,1:s)=C2/norm(C2,'fro');
else
    [V(1:n,1:s),~]=qr(C1,0);
    [W(1:nB,1:s),~]=qr(C2,0);
end

l=size(Y01,2);
Z1=zeros(n,l);
Z2=zeros(n,l);

for j=1:d-1
    
    % sequence in A
    UpA = A*V(:,end-s+1:end);


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
        Z1=Z1+V*Y01((j-ktruncA+1)*s+ 1:(j+1)*s,:);
    end

    %Sequence in B
    UpB = B'*W(:,end-s+1:end);

    % we compute the same exact operations performed in the first Arnoldi
    % loop. We notice that this leads to a more stable procedure,
    % especially in the block case (size(C2,2)>1)
    tt=1;
    for jj=max(1,j-ktruncB+1):j
        k1=(jj-1)*s+1; k2=jj*s;
        t=T(k1:k2,(j-1)*s+1:j*s);
        UpB=UpB-W(1:nB,(tt-1)*s+1:tt*s)*t;
        tt=tt+1;
    end
        
    tt=1;
    for jj=max(1,j-ktruncB+1):j
        k1=(jj-1)*s+1; k2=jj*s;
        t=T1(k1:k2,(j-1)*s+1:j*s);
        UpB=UpB-W(1:nB,(tt-1)*s+1:tt*s)*t;
        tt=tt+1;
    end

    if s==1
       UpB=UpB/norm(UpB);
    else
       [UpB, ~]=qr(UpB,0);
    end
        
    % move the computed basis blocks if necessary
    if j-ktruncB+1<1
        W(:,(tt-1)*s+1:tt*s)=UpB;
    else
        W(:,1:s*(tt-2))=W(:,s+1:s*(tt-1));
        W(:,(tt-2)*s+1:(tt-1)*s)=UpB;
    end

    % compute ktrunc-rank components of the solution
    if j-ktruncB+1==0
        Z2=Z2+W*Y02(1:(j+1)*s,:);
    elseif mod(j+1,ktruncB)==0
        Z2=Z2+W*Y02((j-ktruncB+1)*s+1:(j+1)*s,:);
    end

end

k=mod(d,ktruncA);
if k~=0
    Z1=Z1+V(:,end-k*s+1:end)*Y01(end-k*s+1:end,:);
end    
k=mod(d,ktruncB);
if k~=0
    Z2=Z2+W(:,end-k*s+1:end)*Y02(end-k*s+1:end,:);
end    
% further truncation
[Q1,R1]=qr(Z1,0);
[Q2,R2]=qr(Z2,0);

[uY,sY,vY]=svd(full(R1*R2')); 
sY=diag(sY);
is=sum(abs(sY/sY(1))>eps);

Z1 = Q1*(uY(:,1:is)*diag(sqrt(sY(1:is))));
Z2 = Q2*(vY(:,1:is)*diag(sqrt(sY(1:is))));

  
    
