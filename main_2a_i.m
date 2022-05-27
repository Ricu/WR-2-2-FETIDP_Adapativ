clear; clc;
addpath('libs')

%% Definiere Vorkonditionierer
VK_vec = {...
          'Identitaet',...
          'Dirichlet',...
          'Deflation',...
          'Balancing'...
          };
% VK_vec = {'Deflation'};

TOL_vec = [1,5,10,50,100,500];

%% Parameter fuer PCG
x0 = @(dim) zeros(dim,1); % Startwert
tol = 10^(-7); % Toleranz
% Residuum fuer die Abbruchbedingung
resid = {'vorkonditioniert'}; 
%resid = {'nicht-vorkonditioniert'};


%% Erstelle das Gitter
n = 10; % 2*n^2 Elemente pro Teilgebiet
N = 3;  % Partition in NxN quadratische Teilgebiete
numSD = N^2; % Anzahl Teilgebiete
xyLim = [0,1]; % Gebiet: Einheitsquadrat

[vert,tri] = genMeshSquare(N,n); % Erstelle Knoten- und Elementliste
numVert=size(vert,1);   numTri=size(tri,1); % Anzahl Knoten und Dreiecke
% Erstelle Knoten- und Elementlisten pro Teilgebiet und logische Liste,
% welche Dreiecke in welchem TG sind
[vert__sd,tri__sd,l2g__sd,logicalTri__sd] = meshPartSquare(N,vert,tri); 

% Markiere Dirichletknoten in logischem Vektor
dirichlet = or(ismember(vert(:,1),xyLim), ismember(vert(:,2),xyLim)); 

%% PDE
f = @(vert,y) ones(size(vert));   % Rechte Seite der DGL
% Definiere Koeffizientenfunktion
% rhoMax zum Vergleich der rhoMin/rhoMax zur Konditionszahl
%rhoMax_vec = [1,10,10^2,10^3,10^4,10^5,10^6,10^7,10^8];
rhoMax_vec = 10^6;
rhoMin = 1;
rhoTri = rhoMin*ones(numTri,1);
% Definiere Kanal
xMin=14/30; xMax=16/30;
yMin=3/30;  yMax=27/30;
indVertCanal = (xMin <= vert(:,1) & vert(:,1) <= xMax & yMin <= vert(:,2) & vert(:,2) <= yMax);  % Logischer Vektor, welche Knoten im Kanal liegen
numVertCanal = 1:numVert;
numVertCanal = numVertCanal(indVertCanal); % Knotennummern der Knoten, die im Kanal liegen
for r = 1:length(rhoMax_vec)
    rhoMax = rhoMax_vec(r);
    for i=1:numTri % Iteriere ueber die Elemente
        if ismember(tri(i,:),numVertCanal) % Alle Knoten des Elements liegen im Kanal
        rhoTri(i)=rhoMax;    % Im Kanal entspricht die Koeffizientenfunktion 10^6
        end
    end
end
indElementsCanal = (rhoTri == rhoMax); % Logischer Vektor, welche Elemente im Kanal liegen

rhoTriSD = cell(numSD,1);
for i = 1:numSD
    rhoTriSD{i} = rhoTri(logicalTri__sd{i});
end

%% Definiere maximalen Koeffizienten pro Knoten
maxRhoVert = zeros(numVert,1);
vertTris = cell(numVert,1); 
maxRhoVertSD = cell(numVert,1);
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
    maxRhoVert(i) = max(rhoTri(vertTris{i})); % maximaler Koeffizient pro Knoten
    
    %% Definiere maximalen Koeffizienten pro Knoten teilgebietsweise
    for k = 1:numSD
        vertTrisSD = logicalTri__sd{k}(vertTris{i}); % logischer Vektor welche Dreiecke des Knotens im TG liegen
        maxRhoVertSD{i} = [maxRhoVertSD{i},max(rhoTri(vertTris{i}(vertTrisSD)))]; % maximalen Koeffizienten pro Knoten teilgebietsweise
    end
end


%% Plotten des Gitters mit Kanal
figure("Name","Triangulierung des Gebiets mit Koeffizientenfunktion");
patch('vertices',vert,'faces',tri,'edgecol','k','facecol',[1,1,1]); hold on; axis equal tight;
patch('vertices',vert,'faces',tri(indElementsCanal,:),'edgecol','k','facecol',[.8,.9,1]);
for i = 1:N-1
    line([0,1],[i/N,i/N],'LineWidth', 1, 'color', 'r')
    line([i/N,i/N],[0,1],'LineWidth', 1, 'color', 'r')
end
legend('\rho = 1','\rho = 10^6','Interface','','','')
title("Triangulierung mit Koeffizientenfunktion")

%% Aufstellen der Referenzloesung
% Als Referenzloesung dient die Loesung des global assemblierten Sysmtems
% mit Backslash-Operator
[K,~,b] = assemble(tri,vert,1,f,rhoTri);
K_II = K(~dirichlet,~dirichlet);
b_I = b(~dirichlet);

u_ref = zeros(size(vert,1),1);
u_ref(~dirichlet) = K_II\b_I;

%% Loesen des Systems mit FETI-DP fuer versch. VK
diffs = cell(length(VK_vec),1);
iters = cell(length(VK_vec),1);
kappa_ests = cell(length(VK_vec),1);

fig_VK_comp = figure("Name","Loesungen fuer verschiedene Vorkonditionierer");
tiledlayout('flow')
for vk_ind = 1:length(VK_vec)
    VK = VK_vec{vk_ind};

    [cu,u_FETIDP_glob,~,iters{vk_ind},kappa_ests{vk_ind}] = fetidp_constraint(vert__sd,tri__sd,l2g__sd,f,...
                                                 dirichlet,VK,'adaptive',rhoTriSD,...
                                                 maxRhoVert,maxRhoVertSD,tol,x0,resid);
    diffs{vk_ind} = norm(u_FETIDP_glob-u_ref);

    figure(fig_VK_comp)
    nexttile
    hold on
    for sd = 1:length(tri__sd)
        trisurf(tri__sd{sd},vert__sd{sd}(:,1),vert__sd{sd}(:,2),cu{sd});
    end
    xlabel("x"); ylabel("y"); zlabel("z");
    title(sprintf("Finale Loesung: %s-VK",VK));
    view(3)
    hold off
end

%% Ergebnistabelle
rowNames = ["Anzahl Iterationen","Konditionszahl","Abweichung von Referenzloesung"];
T_results = cell2table([iters';kappa_ests';diffs'],"RowNames",rowNames,"VariableNames",VK_vec)