% Code to reproduce the results collected in Table 5 in
%
% Sketched and Truncated Polynomial Krylov Subspace Methods: Matrix Equations 
% Davide Palitta, Marcel Schweitzer, Valeria Simoncini 
% ArXiv: 2311.16019
%
% Please, acknowledge our work by citing our manuscript whenever you use 
% the software provided here in your research.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
% WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR 
% IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
clear all

addpath('../Functions')
addpath(genpath('../Data'))

% set up problem
r = 1;

for N = 50:10:100


    
    % 3D example with nonconstant convection vector and B\neq A'. 
    [eqn.A, eqn.B, eqn.C1, eqn.C2] = conv_diff_3d_sep(N,r,false);

    n = size(eqn.A,1);

    fprintf('****************** \n n= %d\n ******************\n',n);



    % Maximum Krylov dimension
    opts.m = 250;
    M =opts.m;
    % solve the projected problem and compute the residual norm every check_res
    % iterations
    opts.check_res=20;

    % tolerance
    opts.tol = 1e-6;

    opts.ktrunc=3;
    
    s = 2*r*M;

    % Construct subsampled DCT
    opts.hS = setup_sketching_handle(n,s);



    %% SKETCHING

    tt=tic;
    [out_sketched]=arnoldi_tr_Sylv_whitening_twopass(eqn,opts);
    time_whitening=toc(tt);

    [~,R1]=qr(opts.hS([eqn.A*out_sketched.Z1, out_sketched.Z1,-eqn.C1]),0);
    [~,R2]=qr(opts.hS([out_sketched.Z2, eqn.B'*out_sketched.Z2,eqn.C2]),0);
    actualres=norm(R1*R2','fro')/(norm(opts.hS(eqn.C1),'fro')*norm(opts.hS(eqn.C2),'fro'));
    fprintf('sketched-and-truncated Arnoldi \n its: %d, CPU Time: %e, Computed res norm (sketched norm): %e, \n Actual res norm (sketched norm): %e, ktrunc (A): %d, ktrunc (B): %d\n Rank of the computed solution: %d, Storage demand: %d\n',...
            out_sketched.it,time_whitening,out_sketched.res_vec(end), actualres,opts.ktrunc,opts.ktrunc,size(out_sketched.Z1,2),max(2*r*sum(opts.ktrunc),2*size(out_sketched.Z1,2)))


    %% RKSM
    params.m=opts.m;
    params.tol=opts.tol;
    params.ch=0;
    params.period=1;
    opts.tol=1e-1;
    if ~isfield(params,'smax') && ~isfield(params,'smin')
        s1=real(eigs(-eqn.A,1,'largestreal',opts));
        if isnan(s1)
            s1=eigs(-.5*(eqn.A+eqn.A'),1,'largestreal',opts);
        end
        s2=real(eigs(-eqn.A,1,'smallestreal',opts));
        if isnan(s2)
            s2=eigs(-.5*(eqn.A+eqn.A'),1,'smallestreal',opts);
        end
        params.smax=s1;
        params.smin=s2;
    end
    if ~isfield(params,'smaxB') && ~isfield(params,'sminB')
        s1=real(eigs(-eqn.B,1,'largestreal',opts));
        if isnan(s1)
            s1=eigs(-.5*(eqn.B+eqn.B'),1,'largestreal',opts);
        end        
        s2=real(eigs(-eqn.B,1,'smallestreal',opts));
        if isnan(s2)
            s2=eigs(-.5*(eqn.B+eqn.B'),1,'smallestreal',opts);
        end
        params.smaxB=s1;
        params.sminB=s2;
    end

    params.iter=1;
    params.period=1;
    tt=tic;
    [Z1_rksm,Z2_rksm,nrmrestotnew]=RKSM_Sylv_real(eqn.A,eqn.B,eqn.C1,eqn.C2,params);
    time_RKSM=toc(tt);

    [~,R1]=qr([eqn.A*Z1_rksm, Z1_rksm,-eqn.C1],0);
    [~,R2]=qr([Z2_rksm, eqn.B'*Z2_rksm,eqn.C2],0);
    actualres=norm(R1*R2','fro')/norm(R1*R1');

    fprintf('RKSM \n its: %d, CPU Time: %e, Computed res norm : %e, \n Actual res norm (sketched norm): %e\n Rank of the computed solution: %d, Storage demand: %d\n',...
            length(nrmrestotnew),time_RKSM,nrmrestotnew(end), actualres,size(Z1_rksm,2),max(2*(length(nrmrestotnew)+1),2*size(Z1_rksm,2)))

end


