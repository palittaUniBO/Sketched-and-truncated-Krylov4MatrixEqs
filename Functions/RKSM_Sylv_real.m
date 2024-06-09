function [Z1,Z2,nrmrestotnew]=RKSM_Sylv_real(A,B,C1,C2,params)
% function [Z1,Z2,nrmrestotnew]=RKSM_Sylv_real(A,B,C1,C2,params    
% Approximately Solve  
%                A X + X B + C1 C2' = 0
%     X \approx Z Z'
%
% by the Rational Krylov subspace method 
%
% Input:

% A               coeff. matrix.  n x n
% B               rhs factor   n x p

% params.m        max space dimension allowed
% params.tol      Algebraic stopping tolerance 
% params.smin,params.smax estimates for real spectral interval
%                 associated with field of values of A'
%                 e.g., smax=norm(A,1); smin=smax/condest(A);
% params.ch       ch=1  complex poles  ch=0 real poles
% params.period   how often check convergence (period=1 means at each iteration)
% params.iter     iterative inner solver
%

% Hints:
% 2) Provide "comfortable" (loose bounds) estimates s1, emax

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If this code is used, please cite:
% For the all-real version of the code:
%
% Kirsten, G., & Simoncini, V. (2019). Order reduction methods for solving 
% large-scale differential matrix Riccati equations. 
% arXiv preprint arXiv:1905.12119.
%
% For the overall rational Krylov Lyapunov eqn solver:
%
% V. Druskin & V.Simoncini 
% Adaptive rational Krylov subspaces for large-scale dynamical systems 
% Systems & Control Letters, 60 (2011), pp. 546-560. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m=params.m;
tol=params.tol;
s1=params.smin;
emax=params.smax;
s1B=params.sminB;
emaxB=params.smaxB;
ch=params.ch;
period=params.period;
   
dim1 = [];

n=size(A,1);
C1=full(C1);
C2=full(C2);
p=size(C1,2);
In=speye(n);

csize = size(C1,2);

[V,rr1]=qr(C1,0);
[W,rr2]=qr(C2,0);

nrmb=norm(rr1*rr2','fro'); 
if params.iter
    setup = struct('type','nofill'); %1e-1);
    [M1A,M2A]=ilu(A,setup); [M1B,M2B]=ilu(B',setup);
 
    opt.Tol= tol*1e-2;
    opt.MaxIt= 100;
    opt.Disp= 0;
end



VV=V;
WW=W;

H=zeros(p*(m+2),p*(m+1));
T=zeros(p*(m+2),p*(m+1));

nrmrestotnew=[];

if (norm(A-A',1)<1e-14)
    symm=1;
else
    symm=0;
end

if (norm(B-B',1)<1e-14)
    symmB=1;
else
    symmB=0;
end


% iteration in A
newAv=A*V;
K=full(V'*newAv); 
s=s1(1);
eH=eig(K);
eHpoints = sort([s1(:).',emax]);
snew=newpolei(eHpoints,eH,s*ones(p,1));
       
 if real(snew)<0
     snew=-real(snew)+sqrt(-1)*imag(snew);
 end
s=[s,snew];

% iteration in B
newBv=B*W;
KK=full(W'*newBv); 
sB=s1B(1);
eHB=eig(KK);
eHpointsB = sort([s1B(:).',emaxB]);
snewB=newpolei(eHpointsB,eHB,s1B(1)*ones(p,1));

 if real(snewB)<0
     snewB=-real(snewB)+sqrt(-1)*imag(snewB);
 end
sB=[sB,snewB];
 
 
% additional steps
cmplxflag=0;
cmplxflagB=0;


i=0;
iB=0;

while i < m


  % iteration in A
  i=i+1;

  
  paired=0;
  itp=1;
  cflag = 0;
  Vwrk = V;

  while (paired==0)

    i1=i+1;
    w=Vwrk;
 
    if params.iter
     [wrk1,~,~,ITER(i)] = bicgstabl(A-snew*In,w,opt.Tol,opt.MaxIt,M1A,M2A);
    else
     wrk1 = (A-snew*In)\w; 
    end
     
     %%%% All real basis implementation for RKSM %%%%%
     
     if norm(imag(wrk1)) ~= 0 && cflag == 0
         wrk = real(wrk1);
         cflag = 1;
     elseif norm(imag(wrk1)) ~= 0 && cflag == 1
         wrk = imag(wrk1);
     else
         wrk = wrk1;
     end
     
     
    % Gram-Schmidt step
    jms=(i-1)*p+1;j1s=(i+1)*p;js=i*p;js1=js+1;
    C1_m(jms:js,1:csize)= V'*C1;
 
    for it=1:2
      for kk=1:i
        k1=(kk-1)*p+1; k2=kk*p; 
        ww=VV(1:n,k1:k2);
        gamma=ww'*wrk;
        H(k1:k2,jms:js) = H(k1:k2,jms:js)+ gamma;
        wrk = wrk - ww*gamma;
      end
    end
    
    [V,hinv]=qr(wrk,0); H(js1:j1s,jms:js)=hinv; %hinv = inv(hinv);
    
    if (cmplxflag)
        snew=conj(snew); 
        s=[s,snew];
        cmplxflag=0;
        newAv=A*V;
        
        g = VV'*newAv;
        g1 = g; 
        wrkVA=V'*A;
        g2 = wrkVA*VV; g3 = wrkVA*V;
        
        K = [K g1; g2, g3];
        VV=[VV,V];
        i=i+1; itp=itp+1;
    else 
        paired=1; 
    end
  end

  ih1=i1; ih=i;
  newAv=A*V;

  g = VV'*newAv;
    
    
  if (symm), K=(K+K')/2; end

  
  % iteration in B
  iB=iB+1;

  pairedB=0;
  itpB=1;
  cflagB = 0;
  Wwrk = W;

  while (pairedB==0)

    i1B=iB+1;
    w=Wwrk;
 
    if params.iter
        [wrk1,~,~,ITERB(i)] = bicgstabl(B'-snew*In,w,opt.Tol,opt.MaxIt,M1B,M2B);
    else
        wrk1 = (B'-snew*In)\w; 
    end
     
     %%%% All real basis implementation for RKSM %%%%%
     
     if norm(imag(wrk1)) ~= 0 && cflagB == 0
         wrk = real(wrk1);
         cflagB = 1;
     elseif norm(imag(wrk1)) ~= 0 && cflagB == 1
         wrk = imag(wrk1);
     else
         wrk = wrk1;
     end
     
     
    % Gram-Schmidt step
    jms=(iB-1)*p+1;j1s=(iB+1)*p;js=iB*p;js1=js+1;
    C2_m(jms:js,1:csize)= W'*C2;
 
    for it=1:2
      for kk=1:iB
        k1=(kk-1)*p+1; k2=kk*p; 
        ww=WW(1:n,k1:k2);
        gamma=ww'*wrk;
        T(k1:k2,jms:js) = T(k1:k2,jms:js)+ gamma;
        wrk = wrk - ww*gamma;
      end
    end
    
    [W,tinv]=qr(wrk,0); T(js1:j1s,jms:js)=tinv; %hinv = inv(hinv);
    
    if (cmplxflagB)
        snewB=conj(snewB); 
        sB=[sB,snewB];
        cmplxflagB=0;
        newBv=B*W;
        
        gB = WW'*newBv;
        g1B = gB; 
        g2B = W'*B*WW; g3B = W'*B*W;
        
        KK = [KK g1B; g2B, g3B];
        WW=[WW,W];
        iB=iB+1; itpB=itpB+1;
    else 
        pairedB=1; 
    end
  end

  ih1B=i1B; ihB=iB;
  newBv=B*W;

  gB = WW'*newBv;
    
    
  if (symmB), KK=(KK+KK')/2; end

  
  if (rem(i,period)==0)

      % Solve the projected problem
      Y=lyap(K,KK,-C1_m*C2_m');

      % computed residual   (exact, in exact arithmetic) cheaper computation possible
      uA=newAv-VV*g;
      uB=newBv-WW*gB;
     
      UA=(snew*V-uA)*(H(p*ih+1:p*ih1,1:p*ih)/H(1:ih*p,1:ih*p));
      UB =(snew*W-uB)*(T(p*ihB+1:p*ih1B,1:p*ihB)/T(1:ihB*p,1:ihB*p));
      rrA=qr(UA,0); rrA=triu(rrA(1:size(rrA,2),:));
      rrB=qr(UB,0); rrB=triu(rrB(1:size(rrB,2),:));   

      nrmres = sqrt(norm(rrA*Y,'fro')^2+norm(Y*rrB','fro')^2);



     % relative residual norm
     nrmresnew = (nrmres)/nrmb;
 
     nrmrestotnew = [nrmrestotnew, nrmresnew];

     dim = size(VV,2);
     dim1 = [dim1,dim];

     
     if (nrmresnew<tol)
        break
     end
  end 


  % iteration in A
  % New poles and zeros
  eH=sort(eig(K));
  eHorig=eH;
  
  if (ch)                     % Complex poles. Compute set for next complex pole of r_m

    if (any(imag(eH)) ~=0 && max(abs(imag(eH)))>1e-5 && length(eH)>2) % Roots lambdas come from convex hull too
        eH=[eH;-emax];
        ij=convhull(real(eH),imag(eH)); eH=eH(ij);
        ieH=length(eH); missing=ih*p-ieH;
        while missing>0                         % include enough points from the border
            neweH=(eH(1:ieH-1)+eH(2:ieH))/2;missing=ih*p-length(eH);
            eH=[eH;neweH];
        end
    
        eHpoints=-eH;
        eH=eHorig;
        else                                  % if all real eigs, no convex hull possible
            eHpoints = sort([s1; emax.';-real(eH)]);
        end


    else   % Real poles s from real set. Compute complex roots of r_m via Ritz convex hull
     
        if (any(imag(eH)) ~=0 && length(eH)>2)    % Roots lambdas come from convex hull too
            eH=[eH;-s1;-emax.'];
            ij=convhull(real(eH),imag(eH)); eH=eH(ij);
            ieH=length(eH); missing=ih*p-ieH;
            while missing>0 % include enough points from the border
                neweH=(eH(1:ieH-1)+eH(2:ieH))/2;
                eH=[eH;neweH];
                missing=ih*p-length(eH);
            end
            eH=eH(1:ih*p);
        end
        eHpoints = sort([s1; emax.';-real(eH)]);
        eH=eHorig;
    end


    gs=kron(s(2:end),ones(1,p))';

    snew = newpolei(eHpoints,eH,gs);
    if real(snew)<0, snew=-real(snew)+sqrt(-1)*imag(snew);end


    % If pole is complex, include its conjugate

    if (imag(snew) ~=0), cmplxflag=1;end
    s=[s,snew];

    g1 = g; 
    wrkVA=V'*A;
    g2 = wrkVA*VV; g3 = wrkVA*V;

    K = [K g1; g2, g3];
    VV=[VV,V];
    

    % iteration in B
  
    % New poles and zeros
    eHB=sort(eig(KK));
    eHorigB=eHB;
  
    if (ch)                     % Complex poles. Compute set for next complex pole of r_m

        if (any(imag(eHB)) ~=0 && max(abs(imag(eHB)))>1e-5 && length(eHB)>2) % Roots lambdas come from convex hull too
            eHB=[eHB;-emaxB];
            ij=convhull(real(eHB),imag(eHB)); eHB=eHB(ij);
            ieHB=length(eHB); missing=ihB*p-ieHB;
            while missing>0                         % include enough points from the border
                neweHB=(eHB(1:ieHB-1)+eHB(2:ieHB))/2;missing=ihB*-length(eHB);
                eHB=[eHB;neweHB];
            end
    
            eHpointsB=-eHB;
            eHB=eHorigB;
        else                                  % if all real eigs, no convex hull possible
            eHpointsB = sort([s1B; emaxB.';-real(eHB)]);
        end


    else   % Real poles s from real set. Compute complex roots of r_m via Ritz convex hull
     
        if (any(imag(eHB)) ~=0 && length(eHB)>2)    % Roots lambdas come from convex hull too
            eHB=[eHB;-s1B;-emaxB.'];
            ij=convhull(real(eHB),imag(eHB)); eHB=eHB(ij);
            ieHB=length(eHB); missing=ihB*p-ieHB;
            while missing>0 % include enough points from the border
                neweHB=(eHB(1:ieHB-1)+eHB(2:ieHB))/2;
                eHB=[eHB;neweHB];
                missing=ihB*p-length(eHB);
            end
            eHB=eHB(1:ihB*p);
        end
            eHpointsB = sort([s1B; emaxB.';-real(eHB)]);
            eHB=eHorigB;
    end


    gsB=kron(sB(2:end),ones(1,p))';

    snewB = newpolei(eHpointsB,eHB,gsB);
    if real(snewB)<0, snewB=-real(snewB)+sqrt(-1)*imag(snewB);end


    % If pole is complex, include its conjugate

    if (imag(snewB) ~=0), cmplxflagB=1;end
    sB=[sB,snewB];

    g1B = gB; 
    wrkVA=W'*B;
    g2B = wrkVA*WW; g3B = wrkVA*W;

    KK = [KK g1B; g2B, g3B];
    WW=[WW,W];

end


% factored solution 
% Truncation strategy
[uY,sY,vY]=svd(full(Y),0); 
sY=diag(sY);
is=sum(abs(sY/sY(1))>1e-12);
Z1 = VV(:,1:size(Y,1))*(uY(:,1:is)*diag(sqrt(sY(1:is))));
Z2= WW(:,1:size(Y,2))*(vY(:,1:is)*diag(sqrt(sY(1:is))));

fprintf(' \n')
fprintf('BicgSTAB average iteration count in A: %f and B %f',sum(ITER)/length(ITER),sum(ITERB)/length(ITERB))
fprintf(' \n')

structC=whos('C1');
structM1A=whos('M1A');
structM2A=whos('M2A');

structM1B=whos('M1B');
structM2B=whos('M2B');


fprintf('Storing the preconditioner for A is equivalent to storing about %d vectors of length n\n',round((structM1A.bytes+structM2A.bytes)/structC.bytes))
fprintf('Storing the preconditioner for B is equivalent to storing about %d vectors of length n\n',round((structM1B.bytes+structM2B.bytes)/structC.bytes))

 

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r=ratfun(x,eH,s)
 
for j=1:length(x)
r(j)=abs(prod( (x(j)-s)./(x(j)-eH) ));
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r=ratfuni(x,eH,s)

for j=1:length(x)
r(j)=1./abs(prod( (x(j)-s)./(x(j)-eH) ));
end

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function snew=newpolei(eHpoints,eH,s)

for j=1:length(eHpoints)-1
%   snew(j) = fminbnd( @(x) ratfuni(x,eH,s), eHpoints(j),eHpoints(j+1),optimset('TolX',1e-3));
    sval=linspace(eHpoints(j),eHpoints(j+1),20);
   %sval=linspace(eHpoints(j),eHpoints(j+1),100);
    [sf,jx] = max (abs(ratfun(sval,eH,s)));
    snew(j)=sval(jx);
end
[sn,jx]=max(abs(ratfun(snew,eH,s)));
snew=snew(jx);
return
