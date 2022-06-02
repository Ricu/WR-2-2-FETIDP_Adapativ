function [cu,u_FETIDP_glob,lambda,iter,kappa_est,residual,preconditioned_system] = fetidp_constraint(grid_struct,f,pc_param,rho_struct,pcg_param,plot_iteration)
% Input: grid_struct: Structure mit allen Gitterkomponenten:
%        Komponenten: vert__sd,tri__sd,l2g__sd,dirichlet
% Input: f: Function handle fuer rechte Seite der DGL
% Input: pc_param: Structure mit allen Vorkonditionierer
%        Komponenten: VK, constraint_type, adaptiveTOL 
%        constraint types: 'none','non-adaptive','adaptive'
% Input: rho_struct: Structure mit allen Koeffizientenkomponenten
%        Komponenten: rhoTriSD, maxRhoVert, maxRhoVertSD
% Input: pcg_param: Structure mit allen PCG-Parametern
% Input: plot_iteration: Boolean, ob Loesungen in den Iterationen von PCG geplottet werden

% Output: cu: Cell-Array mit Loesungen auf den Teilgebieten
% Output: u_FETIDP_glob: Globaler Loesungsvektor
% Output: lambda: Loesungsvektor auf den LM
% Output: iter: Anzahl Iterationen aus PCG
% Output: kappa_est: Konditionszahlschaetzung aus PCG
% Output: residual: Residuum aus PCG
% Output: preconditioned_system: Explizit aufgestellte Matrix M^(-1)F

%% Structures entpacken
rhoTriSD = rho_struct.rhoTriSD;
maxRhoVert = rho_struct.maxRhoVert;
maxRhoVertSD = rho_struct.maxRhoVertSD;

vert__sd = grid_struct.vert__sd;
tri__sd  = grid_struct.tri__sd;
l2g__sd  = grid_struct.l2g__sd;
dirichlet = grid_struct.dirichlet;

VK = pc_param.VK;
if ~(strcmp('Deflation',VK) || strcmp('Balancing',VK))
    constraint_type = 'none';
else
    constraint_type = pc_param.constraint_type;
    if isfield(pc_param,'adaptiveTol')
        adaptiveTOL = pc_param.adaptiveTol;
    elseif strcmp('adaptive',constraint_type)
       error('Error: Adaptive Nebenbedigungen gewaehlt aber keine Toleranz angegeben') 
    end
end

numSD = length(vert__sd);       % Anzahl Teilgebiete
numVert = length(dirichlet);    % Anzahl Knoten Global

%% Partition der Knoten in primale, duale und Interfaceknoten
% Zaehle Anzahl Teilgebiete, in denen Knoten enthalten ist 
multiplicity = zeros(numVert,1);
for i = 1:numSD
    multiplicity(l2g__sd{i}) = multiplicity(l2g__sd{i}) + 1;
end
gamma = (multiplicity > 1) & ~dirichlet;    % Extrahiere Interfaceknoten
primal = (multiplicity == 4) & ~dirichlet;  % Extrahiere primale Knoten
dual = gamma & ~primal;                     % Extrahiere duale Knoten

%% Partition der Knotengruppen auf die Teilgebiete
cDirichlet = cell(numSD,1); 
cInner = cell(numSD,1);     
cGamma = cell(numSD,1);     
cDual = cell(numSD,1);      
cPrimal = cell(numSD,1);    
cIDual = cell(numSD,1);
for i = 1:numSD
    cDirichlet{i} = dirichlet(l2g__sd{i});      % Dirichletknoten pro TG
    cGamma{i} = gamma(l2g__sd{i});              % Interfaceknoten pro TG
    cInner{i} = ~(cGamma{i} | cDirichlet{i});   % Innere Knoten pro TG
    cPrimal{i} = primal(l2g__sd{i});            % Primale Knoten pro TG
    cDual{i} = dual(l2g__sd{i});                % Duale Knoten pro TG
    cIDual{i} = cInner{i} | cDual{i};           % Innere+Duale Knoten pro TG
end

%% Mappings 
% Mapping: Global -> Gamma
mapGamma = zeros(numVert,1);
mapGamma(gamma)=1:nnz(gamma);

% Mapping: Global -> Primal
mapPrimal = zeros(numVert,1);
mapPrimal(primal)=1:nnz(primal);

% Mapping: Global -> Dual
mapDual = zeros(numVert,1);
mapDual(dual)=1:nnz(dual);

% Mapping: Interface lokal -> Interface global
cGammaMap = cell(numSD,1);
cPrimalMap = cell(numSD,1);
cDualMap = cell(numSD,1);
for i = 1:numSD
    cGammaMap{i} = mapGamma((l2g__sd{i}(cGamma{i})));       
    cPrimalMap{i} = mapPrimal((l2g__sd{i}(cPrimal{i})));    
    cDualMap{i} = mapDual((l2g__sd{i}(cDual{i})));          
end

%% Lagrange Multiplikatoren
cLM = cell(sum(dual),1);
for i = 1:numSD
    i_ind = find(cDual{i}); % Lokale Knotennummern der dualen Knoten des TG
    for j = 1:length(i_ind)
        s=cDualMap{i}(j);   % Lagrangescher-Multipikator
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
    piInd = mapPrimal(l2g__sd{i}(cPrimal{i}));
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
    piInd = mapPrimal(l2g__sd{i}(cPrimal{i}));
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

%% Nebenbedingung Vorarbeit
if strcmp(constraint_type,'adaptive') || strcmp(constraint_type,'non-adaptive')
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
        edgesPrimal{i} = mapPrimal(edgesPrimalGlobal{i});
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

    %% Definiere Matrix U
    if strcmp(constraint_type,'non-adaptive')
        U = zeros(n_LM,numEdges);
        for edgeID = 1:numEdges % Iteriere ueber Kanten
            cnt=1;
            for j = edgesDual{edgeID}'    % Iteriere ueber duale Knoten der Kante (DUALE Knotennummern)
                U(j,edgeID) = maxRhoVert(edgesDualGlobal{edgeID}(cnt));
                cnt = cnt+1;
            end
        end
    end
    if strcmp(constraint_type,'adaptive')
        cU=cell(1,numEdges);
        for edgeID = 1:numEdges
            %% R
            nPrimal = length(edgesPrimalGlobal{edgeID});
            nGamma = [nnz(cGamma{edgesSD(edgeID,1)}),nnz(cGamma{edgesSD(edgeID,2)})];
            nGammaUnass = sum(nGamma);
            nDual = nGamma-nPrimal;
            P_e = zeros(nGammaUnass,nPrimal);
            R_1 = zeros(nGammaUnass,nDual(1));
            R_2 = zeros(nGammaUnass,nDual(2));
            % Knotensortierung in Spalten ist:
            % primal, rest(1), rest(2)
            % Knotensortierung in Zeilen ist:
            % gamma(1), gamma(2) (=primal(1),rest(1),primal(2),rest(2))
            % Alternative Knotensortierung in Zeilen wäre:
            % primal(1), primal(2), rest(1), rest(2)
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


            %% c
            c = ones(nGammaUnass,1) / norm(ones(nGammaUnass,1));

            %% Verallgmeinertes Eigenwertproblem loesen
            [eigenvalues, eigenvectors] = adaptiveEigenvalues(c, pi, P_D_e, S,sigma);
            
            eigenvectors = eigenvectors(:,diag(eigenvalues) > adaptiveTOL);

            % Extrahiere u
            n_EV = size(eigenvectors,2);
            U_temp = zeros(n_LM,n_EV);
            for k = 1:n_EV
                U_temp(subset,:) = B_D_e * S * P_D_e * eigenvectors;
            end
            cU{edgeID} = U_temp;
        end
        U = cell2mat(cU);
    end
    fprintf('Anzahl Nebenbedingungen = %i\n', size(U,2));

    %% Definiere Projektion P
    UFU = U'*hF(U);
    invUFU = UFU\eye(size(UFU));
    P = @(x) U*invUFU*U'*F(cB_B,cK_BB,cK_PiB,S_PiPi,x);
    P_transpose = @(x) F(cB_B,cK_BB,cK_PiB,S_PiPi,U*invUFU*U'*x);
    IminusP = @(x) x-P(x);
    IminusPtranspose = @(x) x-P_transpose(x);
end

%% PCG + Vorkonditionierer
if strcmp('Identitaet',VK)
    % Identitaet VK
    invM  = @(x) x;
else
    % Dirichlet VK
    invM = @(x) dirVKfunction(numSD,cBskal_Delta,cK_DeltaDelta,cK_II,cK_DeltaI,x);
    if strcmp('Deflation',VK) || strcmp('Balancing',VK)
        % Deflation VK
        invM = @(x) IminusP(invM(IminusPtranspose(x)));
        if strcmp('Balancing',VK)
            % Balancing VK
            invM = @(x) invM(x)+U*invUFU*U'*x;
        end
    end
end

preconditioned_system = invM(hF(eye(n_LM)));
% ew = eig(preconditioned_system);
% ew = sort(ew,'descend');
% topEW = ew(1:min(length(ew),50));
% % Plot der 50 groessten EW
% if strcmp('Deflation',VK) || strcmp('Balancing',VK)
%     figure("Name","Die 50 groessten Eigenwerte von invM*F");
%     scatter(1:50,topEW)
%     set(gca,'Yscale','log')
%     title(sprintf("invM: %s-VK",VK));
% end

%% PCG

ploth = @(lambda,iter,VK) plotiter(lambda,iter,VK,cB_B,cK_BB,cK_PiB,cb_B,cPrimalMap, ...
    l2g__sd,cPrimal,cIDual,S_PiPi,f_PiTilde,f_B,tri__sd,vert__sd);
plot_struct = struct("plot_iteration",plot_iteration,"ploth",ploth);


if strcmp(constraint_type,'adaptive') || strcmp(constraint_type,'non-adaptive')
    constraint_struct = struct('U',U,'invUFU',invUFU,'IminusPtranspose',IminusPtranspose);
else
    constraint_struct = struct('U',[],'invUFU',[],'IminusPtranspose',[]);
end
[lambda,iter,kappa_est,residual] = preCG(hF,invM,d,pcg_param,VK,plot_struct,constraint_struct);

% Korrektur bei Deflation-VK notwendig
if strcmp('Deflation',VK)  % Deflation-VK M^-1_PP
    lambdaBar = U*invUFU*U'*d;
    lambda = lambdaBar+lambda;
end

%% Extrahiere Loesung u
[cu,u_FETIDP_glob] = extract_u(lambda,cB_B,cK_BB,cK_PiB,cb_B,cPrimalMap,...
    l2g__sd,cPrimal,cIDual,S_PiPi,f_PiTilde,f_B);

end


%% Plot und Extraktionsfunktionen
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
    pibar = eye(length(c)) - c*c';
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