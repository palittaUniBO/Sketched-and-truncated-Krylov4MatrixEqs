% Calculates all matrices for the linear equations system
% (N_time (x) M_space + M_time (x) A_space) x = F, as well as
% B and a mesh. The actual solution of the heat equation is defined by
% U = x * B on the returned mesh.
function [N_time, M_time, A_space, M_space, F, B, mesh] = ...
    calculate_space_time_matrices(f_space, f_time, T, K, geometry, Hmax,...
    geometricOrder)
model = createpde(1);
importGeometry(model,geometry);
delta_t = T / (K - 1); 

N_time = eye(K) + diag(-ones(K-1,1),1);
N_time = sparse(N_time)';

M_time = eye(K) + diag(ones(K-1,1),1);
M_time = sparse(1/2 * delta_t * M_time)';

applyBoundaryCondition(model,'dirichlet','Face',...
    1:model.Geometry.NumFaces,'u',0);
specifyCoefficients(model,'m',0,'d',1,'c',1,'a',0,'f',f_space{1});
mesh = generateMesh(model,'GeometricOrder',geometricOrder, ...
    'Hmax', Hmax);

M_space = [];
A_space = [];
B = [];

% Calculate the right hand side
Fs = {size(f_space,2),1};
for l=1:size(f_space,2)
    specifyCoefficients(model,'m',0,'d',1,'c',1,'a',0,'f',f_space{l});
    FEM = assembleFEMatrices(model,'nullspace');
    Fs{l} = FEM.Fc;
    if l == 1
        M_space = FEM.M;
        A_space = FEM.Kc;
        B = FEM.B;
    end
end

F = zeros(size(Fs{1},1), K);
for l=1:size(f_space,2)
    f_time_current = f_time{l};
    F_current = zeros(size(Fs{1},1), K);
    for k=1:K
        F_current(:,k) = delta_t .* Fs{l} * 0.5 ...
            * (f_time_current(delta_t * (k - 2)) ...
            + f_time_current(delta_t * (k-1)));
    end
    F = F + F_current;
end
end