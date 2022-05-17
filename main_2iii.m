clear; clc;
addpath('libs')

%% Definiere Vorkonditionierer
VK={'Deflation'};
% VK={'Dirichlet'};

%% Create grid
n = 10; % 2*n^2 Elemente pro Teilgebiet
N = 5;  % Partition in NxN quadratische Teilgebiete
numSD = N^2; % Anzahl Teilgebiete
% Gebiet: Einheitsquadrat
vertLim = [0,1];
yLim = [0,1];

[vert,tri] = genMeshSquare(N,n); % Erstelle Knoten- und Elementliste
numVert=size(vert,1);   numTri=size(tri,1); % Anzahl Knoten und Dreiecke
[vert__sd,tri__sd,l2g__sd,logicalTri__sd] = meshPartSquare(N,vert,tri); % Erstelle Knoten- und Elementlisten pro Teilgebiet
% Dirichletrand fehlt in Aufgabenstellung?!
dirichlet = or(ismember(vert(:,1),vertLim), ismember(vert(:,2),yLim)); % Dirichletknoten, logischer Vektor
[edges,elements_byEdgeIDs,adjacentElements__e] = mesh_edgeList(tri); % Erstelle Kantenliste, etc.
numEdges=size(edges,1); % Anzahl Kanten

%% PDE
f = @(vert,y) ones(size(vert));   % Rechte Seite der DGL
% Definiere Koeffizientenfunktion
rhoMax = 10^6;
rhoMin = 1;
% Setze alle Koeffizienten der Elemente auf 'rhoMin'
pho = rhoMin*ones(numTri,1);
% Setze nun alle 'rhoMax' Koeffizienten.
rho(rand(length(rho),1) < 0.25) = rhoMax;

indElementsrhoMax = (pho == rhoMax); % Logischer Vektor, welche Elemente in rhoMax liegen

%% Definiere maximalen Koeffizienten pro TG 
maxRhoSD = zeros(numSD,1);
for i = 1:numSD
    maxRhoSD(i) = max(pho(logicalTri__sd{i}));
end

%% Definiere maximalen Koeffizienten pro Knoten
maxRhoVert = zeros(numVert,1);
vertTris = cell(numVert,1); % Enthaelt fuer jeden Knoten die Dreiecke in denen er liegt
for i = 1:numVert % Iteriere ueber Knoten
    iVec = i*ones(1,size(tri,2));
    cnt = 1;
    for j = 1:numTri % Iteriere ueber Dreiecke
        testMembership = ismember(iVec,tri(j,:)); 
        if nnz(testMembership) > 1    % Pruefe, ob Knoten im Dreieck liegt
            vertTris{i}(cnt) = j;
            cnt = cnt+1;
        end
    end
    maxRhoVert(i) = max(pho(vertTris{i}));
end


%% Plotten des Gitters mit Kanal
figure()
patch('vertices',vert,'faces',tri,'edgecol','k','facecol',[1,1,1]); hold on; axis equal tight;
patch('vertices',vert,'faces',tri(indElementsrhoMax,:),'edgecol','k','facecol',[.8,.9,1]);

%% Loesen des Systems mit FETI-DP erstmal Identitaet
[cu,u_FETIDP_glob] = fetidp(numSD,vert,numVert,vert__sd,tri__sd,edges,numEdges,l2g__sd,f,dirichlet,VK,maxRhoSD,maxRhoVert,true);
                
%% compare residuals
[K,~,b] = assemble(tri,vert,1,f);
K_II = K(~dirichlet,~dirichlet);
b_I = b(~dirichlet);

u_global = zeros(size(vert,1),1);
u_global(~dirichlet) = K_II\b_I;

diff = u_FETIDP_glob-u_global;
fprintf("Norm der Differenz: %e\n", norm(diff))

% %% Global system PCG
% tol = 10^(-8);
% hK = @(vert) K_II * vert;
% [lambda,resid,iter,kappa_est] = preCG(hK,speye(size(K_II)),b_I,zeros(length(K_II),1),tol);
% fprintf("#### Global assembliertes System ####\n")
% fprintf("Anzahl Iterationen: %i\n",iter)
% fprintf("Schaetzung Konditionszahl: %e\n",kappa_est)


% %% Teil b)
% ploth = @(lambda,iter) plotiter(lambda,iter,cB_B,cK_BB,cK_PiB,cb_B,cPrimalMap, ...
%                                 l2g__sd,cPrimal,cIDual,S_PiPi,f_PiTilde,f_B,tri__sd,vert__sd);
% % preCG(hF,speye(n_LM),d,zeros(n_LM,1),tol,ploth);
