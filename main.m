clear; clc;
addpath('libs')

%% Definiere Vorkonditionierer
VK = {'Deflation','Balancing','Dirichlet','Identitaet'};
% VK={'Deflation'};
% VK={'Balancing'};
% VK = {'Identitaet'};
% VK={'Dirichlet'};

%% Erstelle das Gitter
n = 10; % 2*n^2 Elemente pro Teilgebiet
N = 3;  % Partition in NxN quadratische Teilgebiete
numSD = N^2; % Anzahl Teilgebiete
% Gebiet: Einheitsquadrat
vertLim = [0,1];
yLim = [0,1];

[vert,tri] = genMeshSquare(N,n); % Erstelle Knoten- und Elementliste
numVert=size(vert,1);   numTri=size(tri,1); % Anzahl Knoten und Dreiecke
[vert__sd,tri__sd,l2g__sd,logicalTri__sd] = meshPartSquare(N,vert,tri); % Erstelle Knoten- und Elementlisten pro Teilgebiet
% und logische Liste, welche Dreiecke in welchem TG sind
% Dirichletrand fehlt in Aufgabenstellung?!
dirichlet = or(ismember(vert(:,1),vertLim), ismember(vert(:,2),yLim)); % Dirichletknoten, logischer Vektor
%[edges,elements_byEdgeIDs,adjacentElements__e] = mesh_edgeList(tri); % Erstelle Kantenliste, etc.
%numEdges=size(edges,1); % Anzahl Kanten

%% PDE
f = @(vert,y) ones(size(vert));   % Rechte Seite der DGL

%% Definiere Koeffizientenfunktion (pro Element)
rhoTri = ones(numTri,1);
% Definiere Kanal
xMin=14/30; xMax=16/30;
yMin=3/30;  yMax=27/30;
indVertCanal = (xMin <= vert(:,1) & vert(:,1) <= xMax & yMin <= vert(:,2) & vert(:,2) <= yMax);  % Logischer Vektor, welche Knoten im Kanal liegen
numVertCanal = 1:numVert;
numVertCanal = numVertCanal(indVertCanal); % Knotennummern der Knoten, die im Kanal liegen
for i=1:numTri % Iteriere ueber die Elemente
    if ismember(tri(i,:),numVertCanal) % Alle Knoten des Elements liegen im Kanal
        rhoTri(i)=10^6;    % Im Kanal entspricht die Koeffizientenfunktion 10^6
    end
end
indElementsCanal = rhoTri > 1; % Logischer Vektor, welche Elemente im Kanal liegen

%% Definiere maximalen Koeffizienten pro TG 
rhoTriSD = cell(numSD,1);
maxRhoSD = zeros(numSD,1);
for i = 1:numSD
    rhoTriSD{i} = rhoTri(logicalTri__sd{i});
    maxRhoSD(i) = max(rhoTriSD{i});
end

%% Definiere maximalen Koeffizienten pro Knoten
maxRhoVert = zeros(numVert,1);
vertTris = cell(numVert,1); 
for i = 1:numVert % Iteriere ueber Knoten
    iVec = i*ones(1,size(tri,2));
    cnt = 1;
    for j = 1:numTri % Iteriere ueber Dreiecke
        testMembership = ismember(iVec,tri(j,:)); 
        if nnz(testMembership) > 1    % Pruefe, ob Knoten im Dreieck liegt
            vertTris{i}(cnt) = j;   % Enthaelt fuer jeden Knoten die Dreiecke in denen er liegt
            cnt = cnt+1;
        end
    end
    maxRhoVert(i) = max(rhoTri(vertTris{i}));
end

interface = ((mod(vert(:,1),1/N)==0) | (mod(vert(:,2),1/N)==0))  & (vert(:,1) ~= 0) & (vert(:,2) ~= 0) & (vert(:,1) ~= 1) & (vert(:,2) ~= 1);

%% Plotten des Gitters mit Kanal
figure("Name","Triangulierung des Gebiets mit Koeffizientenfunktion");
patch('vertices',vert,'faces',tri,'edgecol','k','facecol',[1,1,1]); hold on; axis equal tight;
patch('vertices',vert,'faces',tri(indElementsCanal,:),'edgecol','k','facecol',[.8,.9,1]); hold on; axis equal tight;
%plot(vert(interface,1),vert(interface,2),'r.');
for i = 1:N-1
    line([0,1],[i/N,i/N],'LineWidth', 1, 'color', 'g')
    line([i/N,i/N],[0,1],'LineWidth', 1, 'color', 'g')
end
legend('\rho = 1','\rho = 10^6','Interface','','','')
title("Triangulierung mit Koeffizientenfunktion")

%% Loesen des Systems mit FETI-DP erstmal Identitaet
[cu,u_FETIDP_glob,lambda,iter,kappa_est] = fetidp(numSD,vert,numVert,vert__sd,tri__sd,l2g__sd,f,dirichlet,VK,rhoTri,rhoTriSD,maxRhoVert,vertTris,logicalTri__sd,true);

%% Vergleich der Loesung mit Referenzloesung
% Als Referenzloesung dient die Loesung des global assemblierten Sysmtems
% mit Backslash-Operator
[K,~,b] = assemble(tri,vert,1,f,rhoTri);
K_II = K(~dirichlet,~dirichlet);
b_I = b(~dirichlet);

u_global = zeros(size(vert,1),1);
u_global(~dirichlet) = K_II\b_I;

diff = cell(length(VK),1);
for i=1:length(VK)
    diff{i} = norm(u_FETIDP_glob{i}-u_global);
end

%% Ergebnistabelle
T_results = cell2table([iter';kappa_est';diff'],"RowNames",["Anzahl Iterationen","Konditionszahl","Abweichung von Referenzloesung"],"VariableNames",VK)

%% Plot der Loesungen

figure("Name","Loesungen fuer verschiedene Vorkonditionierer")
for j=1:length(VK)
    subplot(2,2,j)
    hold on
    for i = 1:length(tri__sd)
        trisurf(tri__sd{i},vert__sd{i}(:,1),vert__sd{i}(:,2),cu{j}{i});
    end
    xlabel("x"); ylabel("y"); zlabel("z");
    title(sprintf("Finale Loesung: %s-VK",VK{j}));
    view(3)
    hold off
end

% %% Global system PCG
% tol = 10^(-8);
% hK = @(vert) K_II * vert;
% [lambda,resid,iter,kappa_est] = preCG(hK,speye(size(K_II)),b_I,zeros(length(K_II),1),tol);
% fprintf("#### Global assembliertes System ####\n")
% fprintf("Anzahl Iterationen: %i\n",iter)
% fprintf("Schaetzung Konditionszahl: %e\n",kappa_est)
