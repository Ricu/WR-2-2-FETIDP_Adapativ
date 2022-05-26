function [cu,u_FETIDP_glob,lambda,iter,kappa_est] = fetidp_eig(numSD,vert,numVert,vert__sd,tri__sd,l2g__sd,f,dirichlet,VK,rhoTri,rhoTriSD,maxRhoVert,vertTris,logicalTri__sd,plot)

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
% Definiere maximalen Koeffizienten pro dualen Knoten (teilgebietsweise)
maxRhoDual = zeros(length(cLM),2);
dualTris = cell(length(cLM),1);
row_ind_LM = 1;
for i = 1:length(cLM) % Iteriere ueber duale Knoten
    globNum = l2g__sd{cLM{i}(1,1)}(cLM{i}(2,1));    % Globale Knotennummer des dualen Knotens
    dualTris{i} = vertTris{globNum};   % Enthaelt fuer den dualen Knoten die Dreiecke in denen er liegt
    cnt = 1;
    for j = cLM{i}(1,:)    % Iteriere ueber TG in denen dualer Knoten liegt
        dualSDTris = intersect(dualTris{i},find(logicalTri__sd{j})); % Dreiecke des dualen Knotens UND im entsprechenden TG
        maxRhoDual(i,cnt)=max(rhoTri(dualSDTris)); % Maximaler Koeffizient
        cnt = cnt+1;
    end
    sumMaxRhoDual = sum(maxRhoDual(i,:)); % Summer ueber maximale Koeffizienten zu diesem Knoten
    % Verwende zur Skalierung jeweils den Koeffizienten des anderen TG
    cBskal{cLM{i}(1,1)}(row_ind_LM,cLM{i}(2,1)) = (maxRhoDual(i,2)/sumMaxRhoDual)*cB{cLM{i}(1,1)}(row_ind_LM,cLM{i}(2,1));
    cBskal{cLM{i}(1,2)}(row_ind_LM,cLM{i}(2,2)) = (maxRhoDual(i,1)/sumMaxRhoDual)*cB{cLM{i}(1,2)}(row_ind_LM,cLM{i}(2,2));
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

edgesGammaGlobalAll = cell(size(edgesSD));
edgesGammaGlobal = cell(numEdges,1);
edgesGamma = cell(numEdges,1);

edgesPrimalGlobalAll = cell(size(edgesSD));
edgesPrimalGlobal = cell(numEdges,1);
edgesPrimal = cell(numEdges,1);

for i = 1:numEdges % Iteriere ueber TG-Kanten
    for j = 1:size(edgesSD,2) % Iteriere ueber angrenzende TG
        SD = edgesSD(i,j);
        edgesDualGlobalAll{i,j} = l2g__sd{SD}(cDual{SD});   % Enthaelt fuer jedes angrenzende TG die dualen Knoten (GLOBALE Knotennummern)
        edgesGammaGlobalAll{i,j} = l2g__sd{SD}(cGamma{SD});
        edgesPrimalGlobalAll{i,j} =  l2g__sd{SD}(cPrimal{SD});
    end
    edgesGammaGlobal{i} = intersect(edgesGammaGlobalAll{i,1},edgesGammaGlobalAll{i,2});
    edgesGamma{i} = mapGamma(edgesGammaGlobal{i});
    edgesDualGlobal{i} = intersect(edgesDualGlobalAll{i,1},edgesDualGlobalAll{i,2}); % Enthaelt fuer die TG-Kanten die dualen Knoten (GLOBALE Knotennummern)
    edgesDual{i} = mapDual(edgesDualGlobal{i}); % Enthaelt fuer die TG-Kanten die dualen Knoten (DUALE Knotennummern)
    edgesPrimalGlobal{i} = intersect(edgesPrimalGlobalAll{i,1},edgesPrimalGlobalAll{i,2});
    edgesPrimal{i} = mapPi(edgesPrimalGlobal{i});
end
% Umkehrabbildung
g2l__sd = cell(numSD,1);
for sd = 1:numSD
    g2l__sd{sd} = zeros(numVert,1);
    ind = l2g__sd{sd};
    g2l__sd{sd}(ind) = 1:length(l2g__sd{sd});
end

cLocalPrimal = cell(size(edgesSD));
cLocalDual = cell(size(edgesSD));
cLocalGamma = cell(size(edgesSD));
for i = 1:numEdges
    for j = 1:2
        sd = edgesSD(i,j);
        cLocalPrimal{i,j} = g2l__sd{sd}(edgesPrimalGlobal{i});
        cLocalDual{i,j} = g2l__sd{sd}(edgesDualGlobal{i});
        cLocalGamma{i,j} = g2l__sd{sd}(edgesGammaGlobal{i});
    end
end



%% Neu Schurkomplemente aufstellen (eventuell nicht noetig)
cS = cell(numSD,1);
for sd = 1:numSD
    gamma_selected = cGamma{sd};
    inner_local = cInner{sd};
    cS{sd} = cK{sd}(gamma_selected,gamma_selected);
    temp = cK{sd}(inner_local,inner_local) \ cK{sd}(inner_local,gamma_selected);
    cS{sd} = cS{sd} - cK{sd}(gamma_selected,inner_local) * temp;
end

%% NEU Eigenwertproblem (Spaeter Definition von U)
cU=cell(1,numEdges);
% Brauchen pi_bar_ij, pi_ij, P_D_ij, S_ij, I
for edgeID = 1:numEdges
%     edgeID = 4;
    ind = edgesDual{edgeID};
    eInd_gamma = edgesGammaGlobal{edgeID};
    eInd_dual = edgesDualGlobal{edgeID};
    eInd_primal = edgesPrimalGlobal{edgeID};
    
    %% R
    nPrimal = length(eInd_primal);
    nGamma = [nnz(cGamma{edgesSD(edgeID,1)}),nnz(cGamma{edgesSD(edgeID,2)})];
    nGammaUnass = sum(nGamma);
    nDual = nGamma-nPrimal;
    P_e = zeros(nGammaUnass,nPrimal);
    R_1 = zeros(nGammaUnass,nDual(1));
    R_2 = zeros(nGammaUnass,nDual(2));
%     R = [P_e, R_1, R_2];
    % Knotensortierung in Spalten ist:
    % primal, dual(1), dual(2)
    % Knotensortierung in Zeilen ist: 
    % primal(1), primal(2), dual(1), dual(2)
    % Alternative Knotensortierung in Zeilen wäre:
    % gamma(1), gamma(2) (=primal(1),dual(1),primal(2),dual(2))
    P_e(1:nPrimal,1:nPrimal) = eye(nPrimal);
    R_1(nPrimal+1 : nGamma(1),:) = eye(nDual(1));

    P_e(nGamma(1) + 1:nGamma(1)+nPrimal,1:nPrimal) = eye(nPrimal);
    R_2(nGamma(1)+nPrimal+1 : nGammaUnass,:) = eye(nDual(2));
    R = [P_e, R_1, R_2];
    pi = R*((R'*R)\R');
    % pi hat Dimension nGammaUnass x nGammaUnass






    %% P_D und B


    % P_D wird Dimension dim(B_D_e,2) x dim(B_e,2) haben
    % Also muss dim(B_D_e,2) = nGammaUnass entsprechen
    % Also muss dim(B_e,2) = nGammaUnass entsprechen
    B_D_e = cell(1,2);
    B_e = cell(1,2);
    
    for k = 1:2
        sd = edgesSD(edgeID,k);
        B_e{k} = zeros(n_LM,nGamma(k));
        B_D_e{k} = zeros(n_LM,nGamma(k));

        % Der folgende Index listet alle lokalen Knotenindizes vom
        % aktuellen Teilgebiet, welche zum Interface gehoeren aber kein
        % primaler Knoten auf der betrachteten Kante ist. Die primalen
        % Knoten auf der Kante sind lediglich Nullspalten welche vorne
        % angefügt werden.
        relevantGamma = setdiff(find(cGamma{sd}),cLocalPrimal{edgeID,k});
        B_e{k}(:,nPrimal+1:end) = cB{sd}(:,relevantGamma);
        B_D_e{k}(:,nPrimal+1:end) = cBskal{sd}(:,relevantGamma);        
    end
    B_e = cell2mat(B_e);
    B_D_e = cell2mat(B_D_e);

    % Loesche nun alle Zeilen welche nicht zu einem LM auf der
    % betrachteten Kante gehoeren.
    subset = (full(sum(abs(B_e),2)) == 2);
    B_e = B_e(subset,:);
    B_D_e = B_D_e(subset,:);
    P_D_e = B_D_e'*B_e;



    %% S
    % Option 1: erst die Schurkomplemente normal Teilgebietsweise
    % berechenen (siehe oben) und anschliessend umsortieren
%     S_e = cell(2,1);
%     for k = 1:2
%         sd = edgesSD(edgeID,k);
%         S
%     end
%     S_alternativ = blkdiag(S_e{:})


    % Option 2: die Schurkomplemente von Anfang an umsortiert berechenen
    S_e = cell(2,1);
    for k = 1:2
        sd = edgesSD(edgeID,k);
        inner_local = cInner{sd};
        relevantPrimal = cLocalPrimal{edgeID,k};
        relevantGamma = setdiff(find(cGamma{sd}),relevantPrimal);


        S = cell(2,2);

        % Nur in primalen
        gamma_1 = relevantPrimal;
        gamma_2 = relevantPrimal;
        
        S{1,1} = cK{sd}(gamma_1,gamma_2);
        temp = cK{sd}(inner_local,inner_local) \ cK{sd}(inner_local,gamma_2);
        S{1,1} = S{1,1} - cK{sd}(gamma_1,inner_local) * temp;

        % In primalen und restlichen
        gamma_1 = relevantPrimal;
        gamma_2 = relevantGamma;
        
        S{1,2} = cK{sd}(gamma_1,gamma_2);
        temp = cK{sd}(inner_local,inner_local) \ cK{sd}(inner_local,gamma_2);
        S{1,2} = S{1,2} - cK{sd}(gamma_1,inner_local) * temp;
        S{2,1} = S{1,2}';

        % Nur in restlichen
        gamma_1 = relevantGamma;
        gamma_2 = relevantGamma;
        
        S{2,2} = cK{sd}(gamma_1,gamma_2);
        temp = cK{sd}(inner_local,inner_local) \ cK{sd}(inner_local,gamma_2);
        S{2,2} = S{2,2} - cK{sd}(gamma_1,inner_local) * temp;

        S = cell2mat(S);
        S_e{k} = S;
    end
    S = blkdiag(S_e{:});

    sigma = max(diag(S));


    %% Pibar
    c = ones(nGammaUnass,1) / norm(ones(nGammaUnass,1));
%     cct = ones(nGammaUnass) / nGammaUnass;
    
    [eigenvalues, eigenvectors] = adaptiveEigenvalues(c, pi, P_D_e, S,sigma);
    TOL = 100;
    eigenvectors = eigenvectors(:,diag(eigenvalues) > TOL);

    % Extrahiere u
    n_EV = size(eigenvectors,2);
    U_temp = zeros(n_LM,n_EV);
    for k = 1:n_EV
        U_temp(subset,:) = B_D_e * S * P_D_e * eigenvectors;
    end
    cU{edgeID} = U_temp;





















%     %% Alt
%     pi_bar_e = getPiBar(eInd_gamma);
%     pi_e = getPi(eInd_gamma,eInd_primal,eInd_dual);
% 
% 
% 
%     %% P_D
%     B_D_e = cell(1,2);
%     B_e = cell(1,2);
%     nPrimal = length(eInd_primal);
%     nDual = length(eInd_dual);
%     nDual = length(eInd_gamma);
%     for k = 1:2
%         sd = edgesSD(edgeID,k);
%         B_e{k} = zeros(nDual,nDual);
%         B_D_e{k} = zeros(nDual,nDual);
% 
%         i_dual = mapSD{k}(eInd_dual);
%         j_primal = mapSD{k}(eInd_primal);
%         j_dual = mapSD{k}(eInd_dual);
%         
%         B_e{k}(:,1:nPrimal) = cB{sd}(i_dual,j_primal);
%         B_e{k}(:,nPrimal+1:end) = cB{sd}(i_dual,j_dual);
% 
%         B_D_e{k}(:,1:nPrimal) = cBskal{sd}(i_dual,j_primal);
%         B_D_e{k}(:,nPrimal+1:end) = cBskal{sd}(i_dual,j_dual);
% % 
% %         B_D_e{k} = zeros(nDual,nGamma);
% %         i_dual = mapSD{k}(edge_indices_dual);
% %         j_primal = mapSD{k}(edge_indices_primal);
% %         j_dual = mapSD{k}(edge_indices_dual);
% %         B_D_e{k}(:,1:nPrimal) = cBskal{sd}(i_dual,j_primal);
% %         B_D_e{k}(:,nPrimal+1:end) = cBskal{sd}(i_dual,j_dual);
% % 
% % 
% % 
% %         B_D_e{k} = cBskal{sd}(mapSD{k}(edge_indices_dual),j_primal);
% %         
%     end
%     P_D_e = cell2mat(B_D_e)'*cell2mat(B_e);
% 
% %     P_D_e = getPD;
% 
%     %% S_ij
%     S_e_primal = cell(2,1);
%     S_e_dual
%     for k = 1:2
%         sd = edgesSD(edgeID,k);
% %         S_GG = zeros(nGamma,nGamma);
%         gamma_loc = cGamma{sd};
%         inner_loc = cInner{sd};
% 
%         S_GG = cK{sd}(gamma_loc,gamma_loc);
%         S_GG = S_GG - cK{sd}(gamma_loc,inner_loc)*(cK{sd}(inner_loc,inner_loc)\cK{sd}(inner_loc,gamma_loc));
% %         S_k(1:nPrimal,1:nPrimal) = S_GG();
% 
%         temp = zeros(length(vert));
%         temp(l2g__sd{sd}(gamma_loc),l2g__sd{sd}(gamma_loc)) = S_GG;
%         S_e{k} = zeros(nDual);
%         S_e{k}(1:nPrimal,1:nPrimal) = temp(eInd_primal,eInd_primal);
%         S_e{k}(nPrimal+1:end,nPrimal+1:end) = temp(eInd_dual,eInd_dual);
%     end
%     S_e = blkdiag(S_e{:});
% 
%     S_e = getS();
% 
%     sigma = nnz(diag(S_e));

    
end
U = cell2mat(cU);

% %% Definiere Matrix U
% U = zeros(n_LM,numEdges);
% for edgeID = 1:numEdges % Iteriere ueber Kanten
%     cnt=1;
%     for j = edgesDual{edgeID}'    % Iteriere ueber duale Knoten der Kante (DUALE Knotennummern)
%         U(j,edgeID) = maxRhoVert(edgesDualGlobal{edgeID}(cnt));
%         cnt = cnt+1;
%     end
% end


%% Definiere Projektion P
UFU = U'*hF(U);
invUFU = UFU\eye(size(UFU));
P = @(x) U*invUFU*U'*F(cB_B,cK_BB,cK_PiB,S_PiPi,x);
P_transpose = @(x) F(cB_B,cK_BB,cK_PiB,S_PiPi,U*invUFU*U'*x);
IminusP = @(x) x-P(x);
IminusP_transpose = @(x) x-P_transpose(x);

%% Definiere Vorkonditionierer
% Vorkonditionierer
idVK= @(x) x;
dirichletVK = @(x) dirVKfunction(numSD,cBskal_Delta,cK_DeltaDelta,cK_II,cK_DeltaI,x);
deflationVK = @(x) IminusP(dirichletVK(IminusP_transpose(x)));
balancingVK = @(x) deflationVK(x)+U*invUFU*U'*x;

%% PCG
tol = 10^(-8);
x0 = zeros(n_LM,1);
% Funktion zum Plotten der Loesungen waehrend der Iteration von PCG
ploth = @(lambda,iter,VK) plotiter(lambda,iter,VK,cB_B,cK_BB,cK_PiB,cb_B,cPrimalMap, ...
    l2g__sd,cPrimal,cIDual,S_PiPi,f_PiTilde,f_B,tri__sd,vert__sd);

lambda = cell(length(VK),1);
iter = cell(length(VK),1);
kappa_est = cell(length(VK),1);
cu = cell(length(VK),1);
u_FETIDP_glob = cell(length(VK),1);
for i=1:length(VK)
    if strcmp('Deflation',VK{i})  % Deflation-VK M^-1_PP
        invM  = @(x) deflationVK(x);
    elseif strcmp('Balancing',VK{i}) % Balancing-VK M^-1_BP
        invM  = @(x) balancingVK(x);
    elseif strcmp('Dirichlet',VK{i})  % Dirichlet-VK
        invM  = @(x) dirichletVK(x);
    elseif strcmp('Identitaet',VK{i})    % Identitaet
        invM  = @(x) idVK(x);
    end

    [lambda{i},~,iter{i},kappa_est{i}] = preCG(hF,invM,d,x0,tol,VK{i},ploth,U,invUFU,d);

    % Korrektur bei Deflation-VK notwendig
    if strcmp('Deflation',VK{i})  % Deflation-VK M^-1_PP
        lambdaBar = U*invUFU*U'*d;
        lambda{i} = lambdaBar+lambda{i};
    end

    %% Extrahiere Loesung u_i
    [cu{i},u_FETIDP_glob{i}] = extract_u(lambda{i},cB_B,cK_BB,cK_PiB,cb_B,cPrimalMap,...
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

%% NEUE
function pi_bar = getPiBar(edge_indices)

n_edge_vert = 2*length(edge_indices);
% Falls problematisch: n_edge_vert = length(edge_indices);
cct = ones(n_edge_vert) / n_edge_vert;
pi_bar = eye(size(cct)) - cct;

end

function pi = getPi(edge_indices_gamma,edge_indices_primal,edge_indices_dual)
nPrimal = length(edge_indices_primal);
nGamma = length(edge_indices_gamma);
nDual = length(edge_indices_dual);
R_primal = cell(2,1);
R_dual = cell(2,1);
[~,primal_loc] = ismember(edge_indices_primal,edge_indices_gamma);
[~,dual_loc] = ismember(edge_indices_dual,edge_indices_gamma);
for k = 1:2
    R_primal{k} = zeros(nGamma,nPrimal);
    % Nummerierung funktioniert entlang der 1. Achse nur da Knoten in
    % der Kante aufsteigend sortiert sind (genauso wie die lokale
    % Knotennummerierung.
    primal_ind = sub2ind(size(R_primal{k}),primal_loc',1:nPrimal);
    R_primal{k}(primal_ind) = 1;

    R_dual{k} = zeros(nGamma,nDual);
    dual_ind = sub2ind(size(R_dual{k}),dual_loc',1:nDual);
    R_dual{k}(dual_ind) = 1;
end
R = [cell2mat(R_primal),blkdiag(R_dual{:})];
pi = R*((R'*R)\R');

% edge_indices_primal = setdiff(edge_indices_gamma,edge_indices_dual);
% i_ind = []
% j_ind = [edge_indices_primal,edge_indices_primal,edge_indices_dual];
% values = 1:length(i_ind);
%
%
%
% R_1 = zeros(length(x),length(x__sd{edgeSDs(1)})+length(x__sd{edgeSDs(1)}));
%
%
%
%
%
%
%
% R_e = cell(2,1);
% edge_indices_primal = setdiff(edge_indices_gamma,edge_indices_dual);
%
% R_primal = zeros(2*length(edge_indices_gamma),2*length(edge_indices_primal));
% R_primal([mapPi;mapPi])
% R_dual_1 = zeros(2*length(edge_indices_gamma),length(edge_indices_dual));
% R_dual_2 = zeros(2*length(edge_indices_gamma),length(edge_indices_dual));
%
% R = [R_primal, R_dual_1, R_dual_2];
%
% % R_k^T Spalten: doppelte Knoten
% edge_indices
% edge_indices_dual
% for i = 1:2
%     sd = edgeSDs(i);
%     R_i = zeros(length(cPrimal{edgeSDs(i)}));
%     %         R_i(~,) =
%     R_e{i} =  R_i;
% end

end

function [eigenvalues, eigenvectors] = adaptiveEigenvalues(c, pi, P_D, S,sigma)
% A = Pi*P_D^T*S*P_D*Pi (Korrektur notwendig?)
% Fall 1: c ist NICHT im Kern von S:
% Setze B = B_tilde = Pi * S * Pi + sigma * (I-Pi)
% Fall 2: c ist im Kern von S:
% Setze B = Pibar * B_tilde * Pibar + sigma * (I-Pibar)

    

% LHS = pi_bar * pi * P_D' * S * P_D * pi * pi_bar;
A =  pi * P_D' * S * P_D * pi;

Btilde = pi * S * pi + sigma*(eye(size(pi))-pi);
if norm(S*c) < 10^(-8)
    pibar = eye(nGammaUnass) - c*c';
    B = pibar * Btilde * pibar + sigma*(eye(size(pibar))-pibar); 
else
    B = Btilde;
end


% Korrigiere um evenuelle Asymmetrien
A = 1/2*(A + A');
B = 1/2*(B + B');

% Loese das Eigenwertproblem
[eigenvectors,eigenvalues] = eig(A,B);
end