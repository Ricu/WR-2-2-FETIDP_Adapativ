function [cu,u_FETIDP_glob,lambda,iter,kappa_est,residual,preconditioned_system] = fetidp(grid_struct,f,pc_param,rho_struct,pcg_param,plot_iteration)
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

%% Partition der Knoten in Interfaceknoten und in primale und duale Knoten
% Zaehle Anzahl Teilgebiete, in denen Knoten enthalten ist
multiplicity = zeros(numVert,1);
for i = 1:numSD
    multiplicity(l2g__sd{i}) = multiplicity(l2g__sd{i}) + 1;
end
gamma = (multiplicity > 1) & ~dirichlet;    % Extrahiere Interfaceknoten
primal = (multiplicity == 4) & ~dirichlet;  % Extrahiere primale Knoten
dual = gamma & ~primal;                     % Extrahiere duale Knoten

%% Partition der Knotengruppen teilgebietsweise
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
% Mapping: Global -> Interface
mapGamma = zeros(numVert,1);
mapGamma(gamma)=1:nnz(gamma);

% Mapping: Global -> Primal
mapPrimal = zeros(numVert,1);
mapPrimal(primal)=1:nnz(primal);

% Mapping: Global -> Dual
mapDual = zeros(numVert,1);
mapDual(dual)=1:nnz(dual);

% Mappings: Interface lokal -> Interface global
%           Primal lokal -> Primal global
%           Dual lokal -> Dual global
cGammaMap = cell(numSD,1);
cPrimalMap = cell(numSD,1);
cDualMap = cell(numSD,1);
for i = 1:numSD
    cGammaMap{i} = mapGamma((l2g__sd{i}(cGamma{i})));
    cPrimalMap{i} = mapPrimal((l2g__sd{i}(cPrimal{i})));
    cDualMap{i} = mapDual((l2g__sd{i}(cDual{i})));
end

%% Lagrange Multiplikatoren Info
cLM = cell(sum(dual),1);
for i = 1:numSD
    i_ind = find(cDual{i}); % Lokale Knotennummern der dualen Knoten des TG
    for j = 1:length(i_ind)
        s=cDualMap{i}(j);   % Lagrangescher-Multipikator (duale globale Knotennummer)
        cLM{s}=[cLM{s},[i;i_ind(j)]]; % Enthaelt TG-Nummer und lokale Knotennummer des LM
    end
end

%% Lokale Sprungoperatoren: mit und ohne Skalierung
% Initialisierung
n_LM = sum(multiplicity(dual)-1);   % Anzahl Lagrangescher Multiplikatoren
cB = cell(1,numSD);
cBskal=cell(1,numSD);
for i = 1:numSD
    cB{i}=sparse(n_LM,length(vert__sd{i}));
    cBskal{i}=sparse(n_LM,length(vert__sd{i}));
end

% Sprungoperator ohne Skalierung
row_ind_LM = 1; % Jede Zeile von cB repraesentiert einen LM
for i = 1:length(cLM)   % Iteriere ueber LM
    % Erstes zugehoeriges TG des LM
    cB{cLM{i}(1,1)}(row_ind_LM,cLM{i}(2,1)) = 1;   % cLM{i}(1,.) TG-Nummer, cLM{i}(2,.) lokale Knotennummer
    for j = 2:size(cLM{i},2)    % Iteriere ueber restliche zugehoerige TG des LM
        cB{cLM{i}(1,j)}(row_ind_LM,cLM{i}(2,j)) = -1;
        row_ind_LM = row_ind_LM + 1;
    end
end

% Sprungoperator mit Skalierung
row_ind_LM = 1;
for i = 1:length(cLM) % Iteriere ueber LM
    globNum = l2g__sd{cLM{i}(1,1)}(cLM{i}(2,1));    % Globale Knotennummer des dualen Knotens
    sumMaxRhoDual = sum(maxRhoVertSD{globNum}); % Summer ueber maximale Koeffizienten der zugehoerigen TG zu diesem Knoten
    % Verwende zur Skalierung jeweils den Koeffizienten des anderen TG
    cBskal{cLM{i}(1,1)}(row_ind_LM,cLM{i}(2,1)) = (maxRhoVertSD{globNum}(2)/sumMaxRhoDual)*cB{cLM{i}(1,1)}(row_ind_LM,cLM{i}(2,1));
    cBskal{cLM{i}(1,2)}(row_ind_LM,cLM{i}(2,2)) = (maxRhoVertSD{globNum}(1)/sumMaxRhoDual)*cB{cLM{i}(1,2)}(row_ind_LM,cLM{i}(2,2));
    row_ind_LM = row_ind_LM + 1;
end

%% Assembliere die lokalen Steifigkeitsmatrizen und die lokalen Lastvektoren
cK = cell(numSD,1); % Steifigkeitsmatrizen
cb = cell(numSD,1); % Lastvektoren
for i = 1:numSD
    [cK{i},~,cb{i}] = assemble(tri__sd{i}, vert__sd{i},1,f,rhoTriSD{i});
end

%% Assembliere globale Steifigkeitsmatrix in primalen Variablen
K_PiPiTilde = sparse(sum(primal),sum(primal));
f_PiTilde = sparse(sum(primal),1);
for i = 1:numSD
    K_PiPiTilde(cPrimalMap{i},cPrimalMap{i}) = K_PiPiTilde(cPrimalMap{i},cPrimalMap{i}) + cK{i}(cPrimal{i},cPrimal{i});
    f_PiTilde(cPrimalMap{i}) = f_PiTilde(cPrimalMap{i}) + cb{i}(cPrimal{i});
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

%% Berechne rechte Seite d
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
    numEdges = size(edgesSD,1); % Anzahl der Kanten

    edgesDualGlobalAll = cell(size(edgesSD)); % Enthaelt fuer jedes angrenzende TG die dualen Knoten (GLOBALE Knotennummern)
    edgesDual = cell(numEdges,1);  % Enthaelt fuer die TG-Kanten die dualen Knoten (DUALE Knotennummern)
    edgesDualGlobal = cell(numEdges,1); % Enthaelt fuer die TG-Kanten die dualen Knoten (GLOBALE Knotennummern)

    edgesPrimalGlobalAll = cell(size(edgesSD)); % Enthaelt fuer jedes angrenzende TG die primalen Knoten (GLOBALE Knotennummern)
    edgesPrimalGlobal = cell(numEdges,1); % Enthaelt fuer die TG-Kanten die primalen Knoten (GLOBALE Knotennummern)

    for i = 1:numEdges % Iteriere ueber TG-Kanten
        for j = 1:size(edgesSD,2) % Iteriere ueber angrenzende TG
            SD = edgesSD(i,j);
            edgesDualGlobalAll{i,j} = l2g__sd{SD}(cDual{SD});   % Enthaelt fuer jedes angrenzende TG die dualen Knoten (GLOBALE Knotennummern)
            edgesPrimalGlobalAll{i,j} =  l2g__sd{SD}(cPrimal{SD}); % Enthaelt fuer jedes angrenzende TG die primalen Knoten (GLOBALE Knotennummern)
        end
        edgesDualGlobal{i} = intersect(edgesDualGlobalAll{i,1},edgesDualGlobalAll{i,2}); % Enthaelt fuer die TG-Kanten die dualen Knoten (GLOBALE Knotennummern)
        edgesDual{i} = mapDual(edgesDualGlobal{i}); % Enthaelt fuer die TG-Kanten die dualen Knoten (DUALE Knotennummern)
        edgesPrimalGlobal{i} = intersect(edgesPrimalGlobalAll{i,1},edgesPrimalGlobalAll{i,2}); % Enthaelt fuer die TG-Kanten die primalen Knoten (GLOBALE Knotennummern)
    end
    % Umkehrabbildung von globaler zu lokaler Nummerierung
    g2l__sd = cell(numSD,1);
    for sd = 1:numSD
        g2l__sd{sd} = zeros(numVert,1);
        ind = l2g__sd{sd};
        g2l__sd{sd}(ind) = 1:length(l2g__sd{sd});
    end

    cLocalPrimal = cell(size(edgesSD)); % Enthaelt die auf einer TG-Kante liegenden primalen Knoten des TG in lokaler Nummerierung
    for i = 1:numEdges
        for j = 1:2
            sd = edgesSD(i,j);
            cLocalPrimal{i,j} = g2l__sd{sd}(edgesPrimalGlobal{i}); % Enthaelt die auf einer TG-Kante liegenden primalen Knoten des TG in lokaler Nummerierung
        end
    end

    %% Definiere Matrix U
    if strcmp(constraint_type,'non-adaptive')
        U = zeros(n_LM,numEdges);
        for edgeID = 1:numEdges % Iteriere ueber Kanten
            cnt = 1;
            for j = edgesDual{edgeID}'    % Iteriere ueber duale Knoten der Kante (j = DUALE Knotennummern)
                U(j,edgeID) = maxRhoVert(edgesDualGlobal{edgeID}(cnt)); % Maximaler Koeffizient des dualen Knotens
                cnt = cnt+1;    % Iteriere ueber Anzahl dualer Knoten der Kante
            end
        end
    end
    if strcmp(constraint_type,'adaptive')
        cU=cell(1,numEdges);
        for edgeID = 1:numEdges
            %% Assemblierungmatrix R_ij
            nPrimal = length(edgesPrimalGlobal{edgeID}); % Anzahl primaler Knoten auf der Kante
            nGamma = [nnz(cGamma{edgesSD(edgeID,1)}),nnz(cGamma{edgesSD(edgeID,2)})]; % Anzahl Interfaceknoten pro TG
            nGammaUnass = sum(nGamma); % Anzahl unassemblierte Interface Knoten beider TG
            nRest = nGamma-nPrimal; % Anzahl unassemblierter dualer Knoten beider TG
            % Teile die Matrix R in drei Teile auf:
            % Knotensortierung in Spalten ist:
            % primal, rest(1), rest(2)
            % Knotensortierung in Zeilen ist:
            % gamma(1), gamma(2) (=primal(1),rest(1),primal(2),rest(2))
            P_e = zeros(nGammaUnass,nPrimal); % primal
            R_1 = zeros(nGammaUnass,nRest(1)); % rest(1)
            R_2 = zeros(nGammaUnass,nRest(2)); % rest(2)

            P_e(1:nPrimal,1:nPrimal) = eye(nPrimal); %primal(1),primal
            R_1(nPrimal+1 : nGamma(1),:) = eye(nRest(1)); %rest(1),rest(1)

            P_e(nGamma(1) + 1:nGamma(1)+nPrimal,1:nPrimal) = eye(nPrimal); %primal(2),primal
            R_2(nGamma(1)+nPrimal+1 : nGammaUnass,:) = eye(nRest(2)); % rest(2),rest(2)
            R = [P_e, R_1, R_2]; % Zusammensetzen
            pi = R*((R'*R)\R'); 
            % pi hat Dimension nGammaUnass x nGammaUnass

            %% Sprungoperator und Projektion B und P_D
            % P_D wird Dimension dim(B_D_e,2) x dim(B_e,2) haben
            % Also muss dim(B_D_e,2) = nGammaUnass entsprechen
            % Also muss dim(B_e,2) = nGammaUnass entsprechen
            B_D = cell(1,2);
            B = cell(1,2);

            for k = 1:2
                sd = edgesSD(edgeID,k);
                B{k} = zeros(n_LM,nGamma(k));
                B_D{k} = zeros(n_LM,nGamma(k));

                % Der folgende Index listet alle lokalen Knotenindizes vom
                % aktuellen Teilgebiet, welche zum Interface gehoeren aber kein
                % primaler Knoten auf der betrachteten Kante ist. Die primalen
                % Knoten auf der Kante sind lediglich Nullspalten welche vorne
                % angefuegt werden.
                relevantGamma = setdiff(find(cGamma{sd}),cLocalPrimal{edgeID,k});
                B{k}(:,nPrimal+1:end) = cB{sd}(:,relevantGamma);
                B_D{k}(:,nPrimal+1:end) = cBskal{sd}(:,relevantGamma);
            end
            B = cell2mat(B);
            B_D = cell2mat(B_D);

            % Loesche nun alle Zeilen welche nicht zu einem LM auf der
            % betrachteten Kante gehoeren.
            subset = (full(sum(abs(B),2)) == 2);
            B = B(subset,:);
            B_D = B_D(subset,:);
            P_D = B_D'*B;

            %% Schurkomplement S
            S_temp = cell(2,1);
            % Erstelle nacheinander die lokalen Schurkomplemente
            for k = 1:2
                sd = edgesSD(edgeID,k);
                inner_local = cInner{sd};
                % Finde die lokalen Indizes die zu primal und rest gehoeren
                relevantPrimal = cLocalPrimal{edgeID,k};
                relevantGamma = setdiff(find(cGamma{sd}),relevantPrimal);
                
                % Schreibe die Kombination aus primal und rest in ein 2x2
                % Cell
                S = cell(2,2);

                % Nur in primalen (primal,primal)
                index_1 = relevantPrimal;
                index_2 = relevantPrimal;

                S{1,1} = cK{sd}(index_1,index_2);
                temp = cK{sd}(inner_local,inner_local) \ cK{sd}(inner_local,index_2);
                S{1,1} = S{1,1} - cK{sd}(index_1,inner_local) * temp;

                % In primalen und restlichen (primal,rest)
                index_1 = relevantPrimal;
                index_2 = relevantGamma;

                S{1,2} = cK{sd}(index_1,index_2);
                temp = cK{sd}(inner_local,inner_local) \ cK{sd}(inner_local,index_2);
                S{1,2} = S{1,2} - cK{sd}(index_1,inner_local) * temp;
                S{2,1} = S{1,2}';

                % Nur in restlichen (rest,rest)
                index_1 = relevantGamma;
                index_2 = relevantGamma;

                S{2,2} = cK{sd}(index_1,index_2);
                temp = cK{sd}(inner_local,inner_local) \ cK{sd}(inner_local,index_2);
                S{2,2} = S{2,2} - cK{sd}(index_1,inner_local) * temp;

                S = cell2mat(S);
                S_temp{k} = S;
            end
            S = blkdiag(S_temp{:});
            % Berechne sigma
            sigma = max(diag(S));


            %% Stelle den Vector c auf
            c = ones(nGammaUnass,1) / norm(ones(nGammaUnass,1));

            %% Verallgmeinertes Eigenwertproblem loesen
            [eigenvalues, eigenvectors] = adaptiveEigenvalues(c, pi, P_D, S,sigma);

            % Eigenwerte entsprechend Toleranz auswählen
            eigenvectors = eigenvectors(:,diag(eigenvalues) > adaptiveTOL);

            % Extrahiere U
            % Anzahl an neuen Nebenbedingungen entspricht Anzahl
            % verbliebener Eigenvektoren/-werte
            U_temp = zeros(n_LM,size(eigenvectors,2));
            U_temp(subset,:) = B_D * S * P_D * eigenvectors;
            cU{edgeID} = U_temp;
        end
        % Stelle U auf indem die lokalen U's entlang der 2. Dimension
        % verkettet werden
        U = cell2mat(cU);
    end
    fprintf('%s-Vorkonditionierer\n',VK);
    fprintf('Anzahl Nebenbedingungen = %i, Spaltenrang = %i\n', size(U,2),rank(U,1e-16));
    fprintf('Anzahl Zeilen von U/Anzahl LM= %i\n',size(U,1) );

    %% Definiere Projektion P (und weitere nuetzliche function handles)
    UFU = U'*hF(U);
    invUFU = UFU\eye(size(UFU));
    P = @(x) U*invUFU*U'*F(cB_B,cK_BB,cK_PiB,S_PiPi,x);
    P_transpose = @(x) F(cB_B,cK_BB,cK_PiB,S_PiPi,U*invUFU*U'*x);
    IminusP = @(x) x-P(x);
    IminusPtranspose = @(x) x-P_transpose(x);
end

%% Aufstellen des Vorkonditionierers
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

% Stelle das vorkonditionierte System explizit auf. Anhand dessen EW laesst
% sich die Kkonditionszahl abschaetzen
preconditioned_system = invM(hF(eye(n_LM)));

%% PCG
% Plotfunktion zum Plotten der Loesungen
ploth = @(lambda,iter,VK) plotiter(lambda,iter,VK,cB_B,cK_BB,cK_PiB,cb_B,cPrimalMap, ...
    l2g__sd,cPrimal,cIDual,S_PiPi,f_PiTilde,f_B,tri__sd,vert__sd);
plot_struct = struct("plot_iteration",plot_iteration,"ploth",ploth);

% Definiere constraint-Komponenten
if strcmp(constraint_type,'adaptive') || strcmp(constraint_type,'non-adaptive')
    constraint_struct = struct('U',U,'invUFU',invUFU,'IminusPtranspose',IminusPtranspose);
else
    constraint_struct = struct('U',[],'invUFU',[],'IminusPtranspose',[]);
end

% Loesen des Sytems mit PCG
[lambda,iter,kappa_est,residual] = preCG(hF,invM,d,pcg_param,VK,plot_struct,constraint_struct);

% Korrektur der Loesung lambda bei Deflation-VK notwendig
if strcmp('Deflation',VK)
    lambda = lambda+U*invUFU*U'*d;
end

%% Extrahiere finale Loesung u
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
[cu,u_FETIDP_glob] = extract_u(lambda,cB_B,cK_BB,cK_PiB,cb_B,cPrimalMap, ...
    l2g__sd,cPrimal,cIDual,S_PiPi,f_PiTilde,f_B);
% subplot(1,4,iter+1)
nexttile
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