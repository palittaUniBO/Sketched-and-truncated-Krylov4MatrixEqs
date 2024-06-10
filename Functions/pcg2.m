function [x,k]=pcg2(a,b,x,maxit,tol,D)
%function [x,k]=pcg2(a,b,x,maxit,tol,D)

%tic
n=length(b);
beta=0; 
r=b-a*x; ap=0*r;
gamma=r'*r;
res0=sqrt(trace(gamma));
res=res0;
k=0;
while (res/res0 > tol & k<maxit)

  z= D\r;
  k=k+1;
  gamma=r'*z;
  if (k==1), p=z;else, beta=gamma/gamma0;p=z+p*beta;end
  ap=a*p;
  delta=p.'*ap;
  alfa = gamma/delta;
  x = x + p*alfa;
  r = r - ap*alfa;
  gamma0=gamma;
  res=norm(r,'fro');

end	
