% Code to reproduce the results collected in Figure 4 in
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

%mydefaults

% set up problem
N = 300;

% nu indicates the viscosity parameter
nu=1e-3;
ktrunc=10;
r=1;

[eqn.A, eqn.B, eqn.C1, eqn.C2] = conv_diff_2d_differentConv(N,r,nu,false);

n = size(eqn.A,1);

opts.check_res=1;
opts.store=1;
opts-ktrunc=10;
% Maximum Krylov dimension
opts.m = 800;
M =opts.m;
s = 2*r*M;


% Construct subsampled DCT
opts.hS = setup_sketching_handle(n,s);

% tolerance
opts.tol = 1e-6;

[out_sketched]=arnoldi_tr_Sylv_whitening_twopass(eqn,opts);
            
semilogy(out_sketched.condB*eps,'*b-')
hold on
semilogy(out_sketched.condA*eps,'or-')
semilogy(out_sketched.res_vec,'k-')
leg=legend('$\mathtt{u}\cdot\kappa(\mathbf{V}_d)$','$\mathtt{u}\cdot\kappa(\mathbf{U}_d)$','Rel. Res.');
set(leg,'interpreter','latex')

