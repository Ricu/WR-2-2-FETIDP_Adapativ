function [cu,u_FETIDP_glob,lambda,iter,kappa_est] = fetidp(vert__sd,tri__sd,l2g__sd,f,dirichlet,VK,rhoTriSD,maxRhoVert,maxRhoVertSD,tol,x0,resid)
numSD = length(vert__sd);
numVert = length(dirichlet);

%% Create logical vectors
cDirichlet = cell(numSD,1); % Dirichlet Knoten
cInner = cell(numSD,1);     % Innere Knoten Lokal
cGamma = cell(numSD,1);     % Interface Knoten Lokal
cDual = cell(numSD,1);      % Duale Knoten Lokal
cPrimal = cell(numSD,1);    % Primale Knoten Lokal
cIDual = cell(numSD,1);     % Innere und Duale Knoten Lokal

multiplicity = zeros(numVert,1);
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

%% Sprungoperator aufstellen: ohne Skalierung
row_ind_LM = 1;
for i = 1:length(cLM)
    cB{cLM{i}(1,1)}(row_ind_LM,cLM{i}(2,1)) = 1;
    for j = 2:size(cLM{i},2)
        cB{cLM{i}(1,j)}(row_ind_LM,cLM{i}(2,j)) = -1;
        row_ind_LM = row_ind_LM + 1;
    end
end

%% Sprungoperator aufstellen: mit Skalierung
row_ind_LM = 1;
for i = 1:length(cLM) % Iteriere ueber duale Knoten
    globNum = l2g__sd{cLM{i}(1,1)}(cLM{i}(2,1));    % Globale Knotennummer des dualen Knotens
    sumMaxRhoDual = sum(maxRhoVertSD{globNum}); % Summer ueber maximale Koeffizienten der TG zu diesem Knoten
    % Verwende zur Skalierung jeweils den Koeffizienten des anderen TG
    cBskal{cLM{i}(1,1)}(row_ind_LM,cLM{i}(2,1)) = (maxRhoVertSD{globNum}(2)/sumMaxRhoDual)*cB{cLM{i}(1,1)}(row_ind_LM,cLM{i}(2,1));
    cBskal{cLM{i}(1,2)}(row_ind_LM,cLM{i}(2,2)) = (maxRhoVertSD{globNum}(1)/sumMaxRhoDual)*cB{cLM{i}(1,2)}(row_ind_LM,cLM{i}(2,2));
    row_ind_LM = row_ind_LM + 1;
end


%% Assembliere die Steifigkeitsmatrix und den Lastvektor
cK = cell(numSD,1); % Steifigkeitsmatrizen
cb = cell(numSD,1); % Rechte Seiten

for i = 1:numSD
    [cK{i},~,cb{i}] = assemble(tri__sd{i}, vert__sd{i},1,f,rhoTriSD{i});
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

%% Erstelle function handle auf Systemmatrix F
hF = @(lambda) F(cB_B,cK_BB,cK_PiB,S_PiPi,lambda);

%% Berechne d
f_B = cell2mat(cb_B);
cb_B_trans = cellfun(@transpose,cb_B,'UniformOutput', false);
d = apply_1(cB_B,cK_BB,cb_B_trans,1);
temp = f_PiTilde - apply_2(cb_B_trans,cK_BB,cK_PiB,1);
temp = apply_1(cB_B,cK_BB,cK_PiB,S_PiPi \ temp);
d = d - temp;

%% Erstelle Kantenlisten
cEdgesSD = cell(1,1);
for i = 1:length(cLM)
    cEdgesSD{1} = [cEdgesSD{1};cLM{i}(1,:)];
end
edgesSD = unique(cEdgesSD{1},'rows'); % Enthaelt die beiden angrenzenden Teilgebietsnummern pro TG-Kante
numEdges = size(edgesSD,1);

edgesDualGlobalAll = cell(size(edgesSD));
edgesDual = cell(numEdges,1);
edgesDualGlobal = cell(numEdges,1);
for i = 1:numEdges % Iteriere ueber TG-Kanten
    for j = 1:size(edgesSD,2) % Iteriere ueber angrenzende TG
        SD = edgesSD(i,j);
        edgesDualGlobalAll{i,j} = l2g__sd{SD}(cDual{SD});   % Enthaelt fuer jedes angrenzende TG die dualen Knoten (GLOBALE Knotennummern)
    end
    edgesDualGlobal{i} = intersect(edgesDualGlobalAll{i,1},edgesDualGlobalAll{i,2}); % Enthaelt fuer die TG-Kanten die dualen Knoten (GLOBALE Knotennummern)
    edgesDual{i} = mapDual(edgesDualGlobal{i}); % Enthaelt fuer die TG-Kanten die dualen Knoten (DUALE Knotennummern)
end

%% Definiere Matrix U
U = zeros(n_LM,numEdges);
for i = 1:numEdges % Iteriere ueber Kanten
    cnt=1;
    for j = edgesDual{i}'    % Iteriere ueber duale Knoten der Kante (DUALE Knotennummern)
        U(j,i) = maxRhoVert(edgesDualGlobal{i}(cnt));
        cnt = cnt+1;
    end
end

%% Definiere Projektion P
UFU = U'*hF(U);
invUFU = UFU\eye(size(UFU));
P = @(x) U*invUFU*U'*F(cB_B,cK_BB,cK_PiB,S_PiPi,x);
P_transpose = @(x) F(cB_B,cK_BB,cK_PiB,S_PiPi,U*invUFU*U'*x);
IminusP = @(x) x-P(x);
IminusP_transpose = @(x) x-P_transpose(x);


%% PCG mit Vorkonditionierer
if strcmp('Identitaet',VK)
    % Identitaet VK
    invM  = @(x) x;
else
    % Dirichlet VK
    invM = @(x) dirVKfunction(numSD,cBskal_Delta,cK_DeltaDelta,cK_II,cK_DeltaI,x);
    if strcmp('Deflation',VK) || strcmp('Balancing',VK)
        % Deflation VK
        invM = @(x) IminusP(invM(IminusP_transpose(x)));
        if strcmp('Balancing',VK)
            % Balancing VK
            invM = @(x) invM(x)+U*invUFU*U'*x;
        end
    end
end

%% PCG
ploth = @(lambda,iter,VK) plotiter(lambda,iter,VK,cB_B,cK_BB,cK_PiB,cb_B,cPrimalMap, ...
    l2g__sd,cPrimal,cIDual,S_PiPi,f_PiTilde,f_B,tri__sd,vert__sd);
x0Vec = x0(n_LM);
[lambda,~,iter,kappa_est] = preCG(hF,invM,d,x0Vec,tol,resid,VK,ploth,U,invUFU,d);

% Korrektur bei Deflation-VK notwendig
if strcmp('Deflation',VK) 
    lambdaBar = U*invUFU*U'*d;
    lambda = lambdaBar+lambda;
end

%% Extrahiere Loesung u_i
[cu,u_FETIDP_glob] = extract_u(lambda,cB_B,cK_BB,cK_PiB,cb_B,cPrimalMap,...
    l2g__sd,cPrimal,cIDual,S_PiPi,f_PiTilde,f_B);


%% PCG + Vorkonditionierer alt
% %% Definiere Vorkonditionierer
% dirichletVK = @(x) dirVKfunction(numSD,cBskal_Delta,cK_DeltaDelta,cK_II,cK_DeltaI,x);
% deflationVK = @(x) IminusP(dirichletVK(IminusP_transpose(x)));
% balancingVK = @(x) deflationVK(x)+U*invUFU*U'*x;
%
% %% PCG
% tol = 10^(-8);
% x0 = zeros(n_LM,1);
% % Funktion zum Plotten der Loesungen waehrend der Iteration von PCG
% ploth = @(lambda,iter,VK) plotiter(lambda,iter,VK,cB_B,cK_BB,cK_PiB,cb_B,cPrimalMap, ...
%     l2g__sd,cPrimal,cIDual,S_PiPi,f_PiTilde,f_B,tri__sd,vert__sd);
%
% lambda = cell(length(VK),1);
% iter = cell(length(VK),1);
% kappa_est = cell(length(VK),1);
% cu = cell(length(VK),1);
% u_FETIDP_glob = cell(length(VK),1);
% for i=1:length(VK)
%     if strcmp('Deflation',VK{i})  % Deflation-VK M^-1_PP
%         invM  = @(x) deflationVK(x);
%     elseif strcmp('Balancing',VK{i}) % Balancing-VK M^-1_BP
%         invM  = @(x) balancingVK(x);
%     elseif strcmp('Dirichlet',VK{i})  % Dirichlet-VK
%         invM  = @(x) dirichletVK(x);
%     elseif strcmp('Identitaet',VK{i})    % Identitaet
%
%     end
%
%     %% Eigenwerte berechnen (50 groessten)
%     invMF = invM(hF(eye(size(U,1))));
%     eigenwerte = eig(invMF);
%     eigenwerte = sort(eigenwerte,'descend');
%     EW50 = eigenwerte(1:50);
%
%     [lambda{i},~,iter{i},kappa_est{i}] = preCG(hF,invM,d,x0,tol,VK{i},ploth,U,invUFU,d);
%
%     % Korrektur bei Deflation-VK notwendig
%     if strcmp('Deflation',VK{i})  % Deflation-VK M^-1_PP
%         lambdaBar = U*invUFU*U'*d;
%         lambda{i} = lambdaBar+lambda{i};
%     end
%
%     %% Extrahiere Loesung u_i
%     [cu{i},u_FETIDP_glob{i}] = extract_u(lambda{i},cB_B,cK_BB,cK_PiB,cb_B,cPrimalMap,...
%         l2g__sd,cPrimal,cIDual,S_PiPi,f_PiTilde,f_B);
%
% end

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

function [cu,u_FETIDP_glob] = plotiter(lambda,iter,VK,cB_B,cK_BB,cK_PiB,cb_B,cPrimalMap,...
    l2g__sd,cPrimal,cIDual,S_PiPi,f_PiTilde,f_B,tri__sd,vert__sd)
[cu,u_FETIDP_glob] = extract_u(lambda,cB_B,cK_BB,cK_PiB,cb_B,...
    cPrimalMap, l2g__sd,cPrimal,cIDual,S_PiPi,f_PiTilde,f_B);
subplot(2,2,iter+1)
hold on
for i = 1:length(tri__sd)
    trisurf(tri__sd{i},vert__sd{i}(:,1),vert__sd{i}(:,2),cu{i});
end
xlabel("x"); ylabel("y"); zlabel("z");
title(sprintf("Plot der Loesung: %s-VK",VK))
subtitle(sprintf("Iteration %g",iter))
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