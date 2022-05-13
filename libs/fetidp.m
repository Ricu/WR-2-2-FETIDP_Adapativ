function [cu,u_FETIDP_glob] = fetidp(numSD,vert,numVert,vert__sd,tri__sd,edges,numEdges,l2g__sd,f,dirichlet,VK,maxRhoSD,maxRhoVert,plot)

%% Create logical vectors
cDirichlet = cell(numSD,1); % Dirichlet Knoten
cInner = cell(numSD,1);     % Innere Knoten Lokal
cGamma = cell(numSD,1);     % Interface Knoten Lokal
cDual = cell(numSD,1);      % Duale Knoten Lokal
cPrimal = cell(numSD,1);    % Primale Knoten Lokal
cIDual = cell(numSD,1);     % Innere und Duale Knoten Lokal

multiplicity = zeros(size(vert,1),1);
for i = 1:numSD
    multiplicity(l2g__sd{i}) = multiplicity(l2g__sd{i}) + 1; % zaehle Vorkommnisse
end
gamma = (multiplicity > 1) & ~dirichlet; % extrahiere Interfaceknoten
primal = (multiplicity == 4) & ~dirichlet; % extrahiere primale Knoten
dual = gamma & ~primal; % extrahiere duale Knoten

%% Ordne lokale logische Vektoren zu
for i = 1:numSD
    cDirichlet{i} = dirichlet(l2g__sd{i}); % Dirichletknoten teilgebietsweise
    cGamma{i} = gamma(l2g__sd{i}); % Interfaceknoten teilgebietsweise
    cInner{i} = ~(cGamma{i} | cDirichlet{i}); % innere Knoten teilgebietsweise
    cPrimal{i} = primal(l2g__sd{i}); % primale Knoten teilgebietsweise
    cDual{i} = dual(l2g__sd{i}); % duale Knoten teilgebietsweise
    cIDual{i} = cInner{i} | cDual{i}; % innere/duale Knoten teilgebietsweise
end

%% FETI DP start

%% Lokales Interface -> globale Nummerierung
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
    cGammaMap{i} = mapGamma((l2g__sd{i}(cGamma{i}))); % Interfaceknoten: lokal zu global
    cPrimalMap{i} = mapPi((l2g__sd{i}(cPrimal{i})));
    cDualMap{i} = mapDual((l2g__sd{i}(cDual{i})));
end

%% Bestimme Indexabbildungen global -> dual und global -> primal
mapDual=zeros(numVert,1);
mapDual(dual) = 1:nnz(dual);
mapPi=zeros(numVert,1);
mapPi(primal) = 1:nnz(primal);

%% Lagrange Multiplikatoren info
cLM = cell(sum(dual),1);
for i = 1:numSD
    i_ind = find(cDual{i});
    for j = 1:length(i_ind)
        s=cDualMap{i}(j);
        cLM{s}=[cLM{s},[i;i_ind(j)]]; % Enthaelt TG-Nummer und lokale Knotennummer
    end
end

%% B^(i) initialisieren: mit und ohne Skalierung
n_LM = sum(multiplicity(dual) - 1);
cB = cell(1,numSD);
cBskal=cell(1,numSD);
for i = 1:numSD
    cB{i}=sparse(n_LM,length(vert__sd{i}));
    cBskal{i}=sparse(n_LM,length(vert__sd{i}));
end

%% Sprungoperator aufstellen: mit und ohne Skalierung
row_ind_LM = 1;
Nx=cell(n_LM,1);   % Liste der Teilgebiete, die LM verbindet
sumRho=zeros(length(cLM),1);
for i = 1:length(cLM)
    cB{cLM{i}(1,1)}(row_ind_LM,cLM{i}(2,1)) = 1;
    Nx{i}=cLM{i}(1,1);
    for j = 2:size(cLM{i},2)
        cB{cLM{i}(1,j)}(row_ind_LM,cLM{i}(2,j)) = -1;
        row_ind_LM = row_ind_LM + 1;
        Nx{i}=[Nx{i},cLM{i}(1,j)];
    end
    sumRho(i) = sum(maxRhoSD(Nx{i}));
end


row_ind_LM = 1;
for i = 1:length(cLM)
    for j = 2:size(cLM{i},2)
        cBskal{cLM{i}(1,1)}(i,:)=cB{cLM{i}(1,1)}(i,:)*(maxRhoSD(cLM{j}(1,j))/sumRho(i)); % jeweils rho des anderen Teilgebiets verwenden
        cBskal{cLM{i}(1,j)}(i,:)=cB{cLM{i}(1,j)}(i,:)*(maxRhoSD(cLM{i}(1,1))/sumRho(i));
        row_ind_LM = row_ind_LM + 1;
    end
end

%% Assembliere die Steifigkeitsmatrix und den Lastvektor

cK = cell(numSD,1); % Steifigkeitsmatrizen
cb = cell(numSD,1); % Rechte Seiten

for i = 1:numSD
    [cK{i},~,cb{i}] = assemble(tri__sd{i}, vert__sd{i},1,f);
end

%% Assembliere globales K in primalen Variablen
K_PiPiTilde = sparse(sum(primal),sum(primal));
f_PiTilde = sparse(sum(primal),1);

for i = 1:numSD
    piInd = mapPi(l2g__sd{i}(cPrimal{i}));
    K_PiPiTilde(piInd,piInd) = K_PiPiTilde(piInd,piInd) + cK{i}(cPrimal{i},cPrimal{i});
    f_PiTilde(piInd) = f_PiTilde(piInd) + cb{i}(cPrimal{i});
end

%% Extrahiere Matrizen
cBskal_Delta=cell(numSD,1);
cB_B = cell(numSD,1);
cK_BB = cell(numSD,1);
cK_DeltaDelta=cell(numSD,1);
cK_II=cell(numSD,1);
cK_DeltaI=cell(numSD,1);
cb_B = cell(numSD,1);
cK_PiB = cell(numSD,1);

for i = 1:numSD
    cBskal_Delta{i}=cBskal{i}(:,cDual{i});
    cB_B{i} = cB{i}(:,cIDual{i});
    cK_BB{i} = cK{i}(cIDual{i},cIDual{i});
    cK_DeltaDelta{i}=cK{i}(cDual{i},cDual{i});
    cK_II{i}=cK{i}(cInner{i},cInner{i});
    cK_DeltaI{i}=cK{i}(cDual{i},cInner{i});
    cb_B{i} = cb{i}(cIDual{i});
    cK_PiB{i} = sparse(nnz(primal),nnz(cIDual{i}));
    piInd = mapPi(l2g__sd{i}(cPrimal{i}));
    cK_PiB{i}(piInd,:) = cK{i}(cPrimal{i},cIDual{i});
end

%% Berechne das Schurkomplement
S_PiPi = K_PiPiTilde;
for i = 1:numSD
    S_PiPi = S_PiPi - cK_PiB{i} *(cK_BB{i}\cK_PiB{i}');
end

%% Erstelle function handle
hF = @(lambda) F(cB_B,cK_BB,cK_PiB,S_PiPi,lambda);

%% Definiere Matrix F

% temp1 = (cK_PiB{1} * (cK_BB{1}\cB_B{1}'));
% for i = 2:length(cB_B)
%     temp1 = temp1 + (cK_PiB{i} * (cK_BB{i}\cB_B{i}'));
% end
% temp1 = S_PiPi \ temp1;  %S_PiPi^-1*K_PiB*K_BB^-1*B_B^T
% temp2 = (cB_B{1} * (cK_BB{1}\cK_PiB{1}'));
% for i = 2:length(cB_B)
%     temp2 = temp2 + (cB_B{i} * (cK_BB{i}\cK_PiB{i}'));
% end
% F1 = temp2 * temp1;     %B_B*K_BB^-1*K_BPi * temp1
% temp3 = (cB_B{1} * (cK_BB{1}\cB_B{1}'));
% for i = 2:length(cB_B)
%     temp3 = temp3 + (cB_B{i} * (cK_BB{i}\cB_B{i}'));
% end
% F2 = temp3;     %B_B*K_BB^-1*B_B^T
% F = F1+F2;

% K_BB = blkdiag(cK_BB{:});
% K_PiB = cell2mat(cK_PiB');
% K_Tilde = [K_BB,K_PiB';
%            K_PiB,K_PiPiTilde];
% B=cell2mat(cB);
% F = B*(K_Tilde\ B');

%% Berechne d
f_B = cell2mat(cb_B);
cb_B_trans = cellfun(@transpose,cb_B,'UniformOutput', false);
d = apply_1(cB_B,cK_BB,cb_B_trans,1);
temp = f_PiTilde - apply_2(cb_B_trans,cK_BB,cK_PiB,1);
temp = apply_1(cB_B,cK_BB,cK_PiB,S_PiPi \ temp);
d = d - temp;


%% Definiere Matrix U
U=zeros(n_LM,numEdges);
dualVertNum = find(dual); % Globale Knotennnummern der dualen Knoten
for i = 1:numEdges % Iteriere ueber Kanten
    for j = 1:n_LM % Iteriere ueber duale Knoten
        xH = dualVertNum(j);
        if ismember(xH,edges(i,:)) % Dualer Knoten gehoert zur Kante
            U(i,j) = maxRhoVert(xH);
        end
    end
end


%% Definiere Projektion P
UFU = @(x) U'*F(cB_B,cK_BB,cK_PiB,S_PiPi,U*x); % TODO: UFU explizit
% aufstellen?! % F explizit aufstellen?!
%UFU = U'*F*U;
invUFU = @(x) UFU(x)\eye(size(UFU(x)));
P = @(x) U*invUFU(U'*F(cB_B,cK_BB,cK_PiB,S_PiPi,x));

%% Definiere Vorkonditionierer
% Vorkonditionierer
idVK= @(x) eye(size(x,1),size(x,1))*x;
dirichletVK = @(x) dirVKfunction(numSD,cBskal_Delta,cK_DeltaDelta,cK_II,cK_DeltaI,x);
deflationVK = @(x) (dirichletVK(idVK(x)'-P(x)'))-P(dirichletVK(idVK(x)'-P(x)'));  % invM = dirichletVK
balancingVK = @(x) deflationVK(x)+U*((U'*F*U)\eye(size(U'*F*U)))*U'*x;


if strcmp('Deflation',VK)  % Deflation-VK M^-1_PP
    invM  = @(x) deflationVK(x);
elseif strcmp('Balancing',VK) % Balancing-VK M^-1_BP
    invM  = @(x) balancingVK(x);
elseif strcmp('Dirichlet',VK)  % Dirichlet-VK
    invM  = @(x) dirichletVK(x);
elseif strcmp('Identitaet',VK)    % Identitaet
    invM  = @(x) idVK(x);
end

%% PCG
tol = 10^(-8);
[lambda,~,iter,kappa_est] = preCG(hF,invM,d,zeros(n_LM,1),tol);
fprintf("#### FETI-DP ####\n")
fprintf("Anzahl Iterationen: %i\n",iter)
fprintf("Schaetzung Konditionszahl: %e\n",kappa_est)

%% Extrahiere Loesung u_i & plot

if plot
    [cu,u_FETIDP_glob] = plotiter(lambda,iter,cB_B,cK_BB,cK_PiB,cb_B,cPrimalMap,...
        l2g__sd,cPrimal,cIDual,S_PiPi,f_PiTilde,f_B,tri__sd,vert__sd);
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

function [cu,u_FETIDP_glob] = plotiter(vert,iter,cB_B,cK_BB,cK_PiB,cb_B,cPrimalMap,...
    l2g__sd,cPrimal,cIDual,S_PiPi,f_PiTilde,f_B,tri__sd,vert__sd)
[cu,u_FETIDP_glob] = extract_u(vert,cB_B,cK_BB,cK_PiB,cb_B,...
    cPrimalMap, l2g__sd,cPrimal,cIDual,S_PiPi,f_PiTilde,f_B);
figure()
hold on
for i = 1:length(tri__sd)
    trisurf(tri__sd{i},vert__sd{i}(:,1),vert__sd{i}(:,2),cu{i});
end
xlabel("x"); ylabel("y"); zlabel("z");
title(sprintf("Plot der FETI-DP Loesung - Iteration: %i",iter));
view(3)
hold off
end


%% Definiere Hilfsfunktionen
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

%Schurkomplement
function [ergebnis]=S_DeltaDeltaiFct(cK_DeltaDelta,cK_II,cK_DeltaI,x)
invcK_II=cK_II\eye(size(cK_II));
ergebnis=(cK_DeltaDelta-cK_DeltaI*invcK_II*cK_DeltaI')*x;
end

function [ergebnis] = dirVKfunction(numSD,cBskal_Delta,cK_DeltaDelta,cK_II,cK_DeltaI,x)
ergebnis=zeros(size(cBskal_Delta{1},1),1);
for i=1:numSD
    Vec1=cBskal_Delta{i}'*x;
    Vec2=S_DeltaDeltaiFct(cK_DeltaDelta{i},cK_II{i},cK_DeltaI{i},Vec1);
    ergebnis=ergebnis+cBskal_Delta{i}*Vec2;
end
end