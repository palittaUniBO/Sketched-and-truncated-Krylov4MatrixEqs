% Code to reproduce the results collected in Table 4 in
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


% Specify the right hand side
% f = f_space{1}*t_time{1} ( + f_space{2}*f_time{2} + ...)
u_exact = @(x,y,z,t) (x-1).*(x+1).*(y-1).*(y+1).*(z-1).*(z+1).*t.^2;

f_space = {}; f_time = {};
f_space{1} = @(location, state) (location.x-1).*(location.x+1) ...
    .*(location.y-1).*(location.y+1).*(location.z-1).*(location.z+1);
f_space{2} = @(location, state) ...
    (location.y-1).*(location.y+1).*(location.z-1).*(location.z+1) ...
    + (location.x-1).*(location.x+1).*(location.z-1).*(location.z+1) ...
    + (location.x-1).*(location.x+1).*(location.y-1).*(location.y+1);
f_time{1} = @(t) 2*t.*(t>0);
f_time{2} = @(t) -2*t.^2.*(t>0);

% f should always be zero for t <= 0!
T = 1;
%geometricOrderSpace = 'quadratic'; % 'linear' or 'quadratic',
geometricOrderSpace = 'linear'; % 'linear' or 'quadratic',
geometry = 'cube.stl'; % = [-1, 1]^3

% Specify the problem resolution (should be higher for actual
% benchmarking)
Hmax=[ .03 0.028 0.026 0.024]; % Maximum element size in space
nvect=[393968, 482768, 685214, 880370];
K=800; % Number of time steps

for i=1:length(Hmax)
    for j=1:length(K)
        fprintf('***************** n=%d ************\n',nvect(i))
        % Space time matrices
        fprintf('create space-time data\n')
        [N_time, M_time, A_space, M_space, F_ST, ~, ~] = ...
            calculate_space_time_matrices(f_space, f_time, T, K(j), ...
            geometry, Hmax(i),  geometricOrderSpace);
        fprintf('done \n')

   

        eqn.A1=A_space;
        eqn.A2=M_space;
        eqn.B1=N_time';
        eqn.B2=M_time';

        [m1,~]=size(M_space); [m2,~]=size(M_time);
       
        FF=reshape(F_ST,m1,m2);
        [uu,ss,vv]=svds(FF,1);
        eqn.C1=uu*ss; eqn.C2=vv;
        
        opts.tol=1e-6;
        opts.m=600;
        opts.check_res=20;
        opts.tol_cg=1e-8;
        
        sdim=2*opts.m;
        opts.ktrunc=3;
        opts.hS = setup_sketching_handle(size(A_space,1),sdim);

 
        %% SKETCHED-AND-TRUNCATED
        
        time=tic;
        [out_sketched] =arnoldi_tr_Sylv_whitening_twopass_onesided(eqn,opts);
        time_sketched=toc(time);

        [~,R1]=qr([A_space*out_sketched.Z1, M_space*out_sketched.Z1,-eqn.C1],0);
        [~,R2]=qr([M_time*out_sketched.Z2, N_time*out_sketched.Z2,eqn.C2],0);
        actualres=norm(R1*R2','fro')/norm(eqn.C1*eqn.C2','fro');
        fprintf('sketched-and-truncated \n its: %d, CPU Time: %e, Computed res norm: %e, \n Actual res norm: %e, ktrunc: %d\n Rank of the computed solution: %d, Storage demand: %d\n',...
              out_sketched.it,time_sketched,out_sketched.res_vec(end), actualres,opts.ktrunc(1),size(out_sketched.Z1,2),max(opts.ktrunc,size(out_sketched.Z1,2)))

         %% FULL ARNOLDI
         opts.ktrunc=opts.m;
 
         tt=tic;
         [out_fullArnoldi]=arnoldi_tr_Sylv_twopass_onesided(eqn,opts);
         time_full=toc(tt);
                         
         [~,R1]=qr([A_space*out_fullArnoldi.Z1, M_space*out_fullArnoldi.Z1,-eqn.C1],0);
         [~,R2]=qr([M_time*out_fullArnoldi.Z2, N_time*out_fullArnoldi.Z2,eqn.C2],0);
         actualres=norm(R1*R2','fro')/norm(eqn.C1*eqn.C2','fro');
         fprintf('full Arnoldi \n its: %d, CPU Time: %e, Computed res norm: %e, \n Actual res norm: %e, ktrunc: %d\n Rank of the computed solution: %d, Storage demand: %d\n',...
               out_fullArnoldi.it,time_full,out_fullArnoldi.res_vec(end), actualres,opts.ktrunc(1),size(out_fullArnoldi.Z1,2),max((out_fullArnoldi.it+1),size(out_fullArnoldi.Z1,2)))
       
        %% Truncated ARNOLDI
        opts.ktrunc=3;
        opts.m=1000;

        tt=tic;
        [out_truncated]=arnoldi_tr_Sylv_twopass_onesided(eqn,opts);
        time_truncated=toc(tt);
                        
        [~,R1]=qr([A_space*out_truncated.Z1, M_space*out_truncated.Z1,-eqn.C1],0);
        [~,R2]=qr([M_time*out_truncated.Z2, N_time*out_truncated.Z2,eqn.C2],0);
        actualres=norm(R1*R2','fro')/norm(eqn.C1*eqn.C2','fro');
        fprintf('truncated Arnoldi \n its: %d, CPU Time: %e, Computed res norm: %e, \n Actual res norm: %e, ktrunc: %d\n Rank of the computed solution: %d, Storage demand: %d\n',...
             out_truncated.it,time_truncated,out_truncated.res_vec(end), actualres,opts.ktrunc(1),size(out_truncated.Z1,2),max(opts.ktrunc,size(out_truncated.Z1,2)))
    end
end


return

