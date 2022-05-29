clear; clc;
addpath('libs')

%% Definiere Vorkonditionierer
VK_vec = {'Deflation',...
          };

tol_vec = {'TOL=10',...
           'TOL=10^3',...
           'TOL=10^6',...
           'TOL=10^10',...
          };

diffs = cell(length(tol_vec),1);
iters = cell(length(tol_vec),1);
kappa_ests = cell(length(tol_vec),1);

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

%% Definiere Kanal-Koeffizientenfunktion
% Definiere Bereich des Kanals
xMin=14/30; xMax=16/30;
yMin=3/30;  yMax=27/30;

% Definiere rho im Kanal und sonst
% rhoCanal zum Vergleich der rhoMin/rhoMax zur Konditionszahl

rhoCanal = 10^6;
rhoNotCanal = 1;
    
% Definiere Koeffizient auf den Elementen (und teilgebietsweise);
% maximalen Koeffizienten pro Knoten (und teilgebietsweise)
[rhoTri,rhoTriSD,indElementsCanal,maxRhoVert,maxRhoVertSD] = coefficient_1(xMin,xMax,yMin,yMax,rhoCanal,rhoNotCanal,vert,tri,numVert,numTri,numSD,logicalTri__sd);
    
    
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
    
%% Aufstellen der Referenzloesun
% Als Referenzloesung dient die Loesung des global assemblierten Sysmtems
% mit Backslash-Operator
[K,~,b] = assemble(tri,vert,1,f,rhoTri);
K_II = K(~dirichlet,~dirichlet);
b_I = b(~dirichlet);
    
u_ref = zeros(size(vert,1),1);
u_ref(~dirichlet) = K_II\b_I;
    
%% Loesen des Systems mit FETI-DP fuer versch. VK
    
TOL_vec = [10,10^3,10^6,10^10];
for t = 1 : length(TOL_vec)
    TOL = TOL_vec(t);
    fig_VK_comp = figure("Name","Loesungen fuer verschiedene Vorkonditionierer");
    tiledlayout('flow')
    VK = VK_vec{1};
        
    [cu,u_FETIDP_glob,~,iters{t},kappa_ests{t}] = fetidp_constraint(TOL,vert__sd,tri__sd,l2g__sd,f,...
                                                         dirichlet,VK,'adaptive',rhoTriSD,...
                                                         maxRhoVert,maxRhoVertSD,tol,x0,resid);
    diffs{t} = norm(u_FETIDP_glob-u_ref);
        
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
fprintf('Deflation-Vorkonditionierer \n')
fprintf('RhoCanal= %g \n',rhoCanal)
fprintf('Konditionszahl fuer unterschiedliche TOL zur Auswahl der EW: \n')
T_results = cell2table([iters';kappa_ests';diffs'],"RowNames",rowNames,"VariableNames",tol_vec);
disp(T_results)