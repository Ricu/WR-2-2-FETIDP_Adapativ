clear; clc;
addpath('libs')

%% Create grid
n = 10;
N = 5;
numSD = N^2;
xLim = [0,1];
yLim = [0,1];

[x,tri] = genMeshSquare(N,n);
[x__sd,tri__sd,l2g__sd] = meshPartSquare(N,x,tri);
dirichlet = or(ismember(x(:,1),xLim), ismember(x(:,2),yLim));
f = @(x,y) ones(size(x));

[cu,u_FETIDP_glob] = fetidp(x,x__sd,tri__sd,l2g__sd,f,dirichlet,true);

                
%% compare residuals
[K,~,b] = assemble(tri,x,1,f);
K_II = K(~dirichlet,~dirichlet);
b_I = b(~dirichlet);

u_global = zeros(size(x,1),1);
u_global(~dirichlet) = K_II\b_I;

diff = u_FETIDP_glob-u_global;
fprintf("Norm der Differenz: %e\n", norm(diff))

% %% Global system PCG
% tol = 10^(-8);
% hK = @(x) K_II * x;
% [lambda,resid,iter,kappa_est] = preCG(hK,speye(size(K_II)),b_I,zeros(length(K_II),1),tol);
% fprintf("#### Global assembliertes System ####\n")
% fprintf("Anzahl Iterationen: %i\n",iter)
% fprintf("Schaetzung Konditionszahl: %e\n",kappa_est)


% %% Teil b)
% ploth = @(lambda,iter) plotiter(lambda,iter,cB_B,cK_BB,cK_PiB,cb_B,cPrimalMap, ...
%                                 l2g__sd,cPrimal,cIDual,S_PiPi,f_PiTilde,f_B,tri__sd,x__sd);
% % preCG(hF,speye(n_LM),d,zeros(n_LM,1),tol,ploth);





