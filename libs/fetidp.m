function [cu,u_FETIDP_glob] = fetidp(x,x__sd,tri__sd,l2g__sd,f,dirichlet,plot)

numSD = length(l2g__sd);
%% Create logical vectors
% dirichlet = or(ismember(x(:,1),xLim), ismember(x(:,2),yLim));
cDirichlet = cell(numSD,1); % Dirichlet Knoten
cInner = cell(numSD,1); % Innere Knoten
cGamma = cell(numSD,1); % Interface Knoten
cDual = cell(numSD,1);
cPrimal = cell(numSD,1);
cIDual = cell(numSD,1);

multiplicity = zeros(size(x,1),1);
for i = 1:numSD
    multiplicity(l2g__sd{i}) = multiplicity(l2g__sd{i}) + 1; % Count occurences
end
gamma = (multiplicity > 1) & ~dirichlet; % Extract interface vertices
primal = (multiplicity == 4) & ~dirichlet; % Extract primal vertices
dual = gamma & ~primal; % Extract dual vertices

%% Assign local logical vectors
for i = 1:numSD
    cDirichlet{i} = dirichlet(l2g__sd{i});
    cGamma{i} = gamma(l2g__sd{i});
    cInner{i} = ~(cGamma{i} | cDirichlet{i});
    cPrimal{i} = primal(l2g__sd{i});
    cDual{i} = dual(l2g__sd{i});
    cIDual{i} = cInner{i} | cDual{i};
end

%% FETI DP start

%% Lokal interface -> globale Nummer
cGammaMap = cell(numSD,1);
mapGamma = zeros(length(gamma),1);
mapGamma(gamma)=1:nnz(gamma);

cPrimalMap = cell(numSD,1);
mapPi = zeros(length(primal),1);
mapPi(primal)=1:nnz(primal);

cDualMap = cell(numSD,1);
mapDual = zeros(length(dual),1);
mapDual(dual)=1:nnz(dual);

for i = 1:numSD
    cGammaMap{i} = mapGamma((l2g__sd{i}(cGamma{i})));
    cPrimalMap{i} = mapPi((l2g__sd{i}(cPrimal{i})));
    cDualMap{i} = mapDual((l2g__sd{i}(cDual{i})));
end

%% Lagrange multiplikatioren info
cLM = cell(sum(dual),1);
for i = 1:numSD
    i_ind = find(cDual{i});
    for j = 1:length(i_ind)
        s=cDualMap{i}(j);
        cLM{s}=[cLM{s},[i;i_ind(j)]];
    end
end

%% B^(i) initialisieren
n_LM = sum(multiplicity(dual) - 1);
cB = cell(numSD,1);
for i = 1:numSD
    cB{i} = sparse(n_LM,length(x__sd{i}));
end

%% Sprungoperator aufstellen
row_ind_LM = 1;
for i = 1:length(cLM)
    for j = 2:size(cLM{i},2)
        cB{cLM{i}(1,1)}(row_ind_LM,cLM{i}(2,1)) = 1;
        cB{cLM{i}(1,j)}(row_ind_LM,cLM{i}(2,j)) = -1;
        row_ind_LM = row_ind_LM + 1;
    end
end
%% Assemble stiffness matrices and load vector

cK = cell(numSD,1); % Steifigkeitsmatrizen
cb = cell(numSD,1); % Rechte Seiten

for i = 1:numSD
    [cK{i},~,cb{i}] = assemble(tri__sd{i}, x__sd{i},1,f);
end

%% Assemble global K
K_PiPiTilde = sparse(sum(primal),sum(primal));
f_PiTilde = sparse(sum(primal),1);

for i = 1:numSD
    piInd = mapPi(l2g__sd{i}(cPrimal{i}));
    K_PiPiTilde(piInd,piInd) = K_PiPiTilde(piInd,piInd) + cK{i}(cPrimal{i},cPrimal{i});
    f_PiTilde(piInd) = f_PiTilde(piInd) + cb{i}(cPrimal{i});
end

%% Extract matrices
cB_B = cell(numSD,1);
cK_BB = cell(numSD,1);
cb_B = cell(numSD,1);
cK_PiB = cell(numSD,1);

for i = 1:numSD
    cB_B{i} = cB{i}(:,cIDual{i});
    cK_BB{i} = cK{i}(cIDual{i},cIDual{i});
    cb_B{i} = cb{i}(cIDual{i});
    cK_PiB{i} = sparse(nnz(primal),nnz(cIDual{i}));
    piInd = mapPi(l2g__sd{i}(cPrimal{i}));
    cK_PiB{i}(piInd,:) = cK{i}(cPrimal{i},cIDual{i}); 
end

%% Compute schurcomplement
S_PiPi = K_PiPiTilde;
for i = 1:numSD
    S_PiPi = S_PiPi - cK_PiB{i} *(cK_BB{i}\cK_PiB{i}');
end

%% create handle
hF = @(lambda) F(cB_B,cK_BB,cK_PiB,S_PiPi,lambda);

%% compute d
f_B = cell2mat(cb_B);
cb_B_trans = cellfun(@transpose,cb_B,'UniformOutput', false);
d = apply_1(cB_B,cK_BB,cb_B_trans,1);
temp = f_PiTilde - apply_2(cb_B_trans,cK_BB,cK_PiB,1);
temp = apply_1(cB_B,cK_BB,cK_PiB,S_PiPi \ temp);
d = d - temp;


%% PCG
tol = 10^(-8);
[lambda,~,iter,kappa_est] = preCG(hF,speye(n_LM),d,zeros(n_LM,1),tol);
fprintf("#### FETI-DP ####\n")
fprintf("Anzahl Iterationen: %i\n",iter)
fprintf("Schaetzung Konditionszahl: %e\n",kappa_est)

%% Extract solution u_i & plot

if plot
    [cu,u_FETIDP_glob] = plotiter(lambda,iter,cB_B,cK_BB,cK_PiB,cb_B,cPrimalMap,...
                    l2g__sd,cPrimal,cIDual,S_PiPi,f_PiTilde,f_B,tri__sd,x__sd);
else
    [cu,u_FETIDP_glob] = extract_u(lambda,cB_B,cK_BB,cK_PiB,cb_B,cPrimalMap,...
                    l2g__sd,cPrimal,cIDual,S_PiPi,f_PiTilde,f_B);
end


end

function [cu,u_glob] = extract_u(lambda,cB_B,cK_BB,cK_PiB,cb_B,cPrimalMap,...
                    l2g__sd,cPrimal,cIDual,S_PiPi,f_PiTilde,f_B)
    numSD = length(cB_B);
    cb_B_trans = cellfun(@transpose,cb_B,'UniformOutput', false);
    % u_pi_tilde
    temp1 = apply_2(cb_B_trans,cK_BB,cK_PiB,1);
    temp2 = apply_2(cB_B,cK_BB,cK_PiB,lambda);
    u_pi_tilde = S_PiPi \ (f_PiTilde - temp1 + temp2);

    % u_B
    temp1 = cell2mat(cK_PiB')' * u_pi_tilde;    %K_BPi_tilde * u_Pi_tilde
    temp2 = cell2mat(cB_B')' * lambda;          %B_B^T * lambda
    u_B = blkdiag(cK_BB{:}) \ (f_B - temp1 - temp2);

    pointer = 1;
    cu = cell(numSD,1);
    for i = 1:numSD
        cu{i} = zeros(length(l2g__sd{i}),1); % Setze Dirichletwerte
        skip = nnz(cIDual{i});
        cu{i}(cIDual{i}) = u_B(pointer : pointer + skip - 1);
        cu{i}(cPrimal{i}) = u_pi_tilde(cPrimalMap{i});
        pointer = pointer + skip;
    end

    u_glob = zeros(length(l2g__sd{1}),1);
    for i = 1:numSD
        u_glob(l2g__sd{i}) = cu{i};
    end
end

function [cu,u_FETIDP_glob] = plotiter(x,iter,cB_B,cK_BB,cK_PiB,cb_B,cPrimalMap,...
                    l2g__sd,cPrimal,cIDual,S_PiPi,f_PiTilde,f_B,tri__sd,x__sd)
    [cu,u_FETIDP_glob] = extract_u(x,cB_B,cK_BB,cK_PiB,cb_B,...
                       cPrimalMap, l2g__sd,cPrimal,cIDual,S_PiPi,f_PiTilde,f_B);
    figure()
    hold on
    for i = 1:length(tri__sd)
        trisurf(tri__sd{i},x__sd{i}(:,1),x__sd{i}(:,2),cu{i});
    end
    xlabel("x"); ylabel("y"); zlabel("z");
    title(sprintf("Plot der FETI-DP Loesung - Iteration: %i",iter));
    view(3)
    hold off
end


%% Define functions 
function y = apply_1(cB_B,cK_BB,cK_PiB,x)
    temp = (cB_B{1} * (cK_BB{1}\cK_PiB{1}')) * x;
    for i = 2:length(cB_B)
        temp = temp + (cB_B{i} * (cK_BB{i}\cK_PiB{i}')) * x;
    end
    y = temp;
end

function y = apply_2(cB_B,cK_BB,cK_PiB,x)
    temp = (cK_PiB{1} * (cK_BB{1}\cB_B{1}')) * x;
    for i = 2:length(cB_B)
        temp = temp + (cK_PiB{i} * (cK_BB{i}\cB_B{i}')) * x;
    end
    y = temp;
end

function y = F(cB_B,cK_BB,cK_PiB,S_PiPi,x)
    temp1 = S_PiPi \ apply_2(cB_B,cK_BB,cK_PiB,x);  %S_PiPi^-1*K_PiB*K_BB^-1*B_B^T * x
	temp1 = apply_1(cB_B,cK_BB,cK_PiB,temp1);       %B_B*K_BB^-1*K_BPi * temp1
    
    temp2 = apply_1(cB_B,cK_BB,cB_B,x);             %B_B*K_BB^-1*B_B^T * x
    y = temp1 + temp2;
end
