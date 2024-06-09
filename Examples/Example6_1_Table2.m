% Code to reproduce the results collected in Table 2 in
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
N = 300;

% nu indicates the viscosity parameter
for r=[1 3]
    for nu=[1e-1 1e-2 1e-3]
        for check_res=[1 20]


            fprintf('****************** \n r= %d, nu=%e, p=%d\n ******************\n',r,nu,check_res);

            [eqn.A, eqn.B, eqn.C1, eqn.C2] = conv_diff_2d_differentConv(N,r,nu,false);

            n = size(eqn.A,1);

            % Maximum Krylov dimension
            opts.m = 800;
            M =opts.m;
            s = 2*r*M;


            % Construct subsampled DCT
            opts.hS = setup_sketching_handle(n,s);


            opts.check_res=check_res;


            % tolerance
            opts.tol = 1e-6;


            %% FULL ARNOLDI
            opts.ktrunc=[opts.m opts.m];

            tt=tic;
            [out_fullArnoldi]=arnoldi_tr_Sylv_twopass(eqn,opts);
            time_full=toc(tt);
                        
            [~,R1]=qr([eqn.A*out_fullArnoldi.Z1, out_fullArnoldi.Z1,-eqn.C1],0);
            [~,R2]=qr([out_fullArnoldi.Z2, eqn.B'*out_fullArnoldi.Z2,eqn.C2],0);
            actualres=norm(R1*R2','fro')/(norm(eqn.C1,'fro')*norm(eqn.C2,'fro'));
            fprintf('full Arnoldi \n its: %d, CPU Time: %e, Computed res norm: %e, \n Actual res norm: %e, ktrunc (A): %d, ktrunc (B): %d\n Rank of the computed solution: %d, Storage demand: %d\n',...
                out_fullArnoldi.it,time_full,out_fullArnoldi.res_vec(end), actualres,opts.ktrunc(1),opts.ktrunc(2),size(out_fullArnoldi.Z1,2),max(2*r*(out_fullArnoldi.it+1),2*size(out_fullArnoldi.Z1,2)))


            %% SKETCHING

            opts.ktrunc=10;
            

            tt=tic;
            [out_sketched]=arnoldi_tr_Sylv_whitening_twopass(eqn,opts);
            time_whitening=toc(tt);

            [~,R1]=qr(opts.hS([eqn.A*out_sketched.Z1, out_sketched.Z1,-eqn.C1]),0);
            [~,R2]=qr(opts.hS([out_sketched.Z2, eqn.B'*out_sketched.Z2,eqn.C2]),0);
            actualres=norm(R1*R2','fro')/(norm(opts.hS(eqn.C1),'fro')*norm(opts.hS(eqn.C2),'fro'));
            fprintf('sketched-and-truncated Arnoldi \n its: %d, CPU Time: %e, Computed res norm (sketched norm): %e, \n Actual res norm (sketched norm): %e, ktrunc (A): %d, ktrunc (B): %d\n Rank of the computed solution: %d, Storage demand: %d\n',...
                out_sketched.it,time_whitening,out_sketched.res_vec(end), actualres,opts.ktrunc,opts.ktrunc,size(out_sketched.Z1,2),max(r*sum(opts.ktrunc),2*size(out_sketched.Z1,2)))


        end
    end
end