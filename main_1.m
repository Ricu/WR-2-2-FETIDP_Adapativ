clear; clc;
addpath('libs')

%% Definiere Vorkonditionierer
VK_vec = {'Identitaet',...
          'Dirichlet',...
          'Deflation',...
          'Balancing'...
          };
      
%% Initialisiere Parameter fuer PCG
x0 = @(dim) zeros(dim,1); % Startvektor
tol = 10^(-8); % Toleranz fuer die Abbruchbedingung
% Residuum fuer die Abbruchbedingung
%resid = {'vorkonditioniert'}; 
%resid = {'nicht-vorkonditioniert'};
resid = {'nicht-vorkonditioniert,alternativ'};

%% Erstelle das Gitter
n = 10; % 2*n^2 Elemente pro Teilgebiet
N = 3;  % Partition in NxN quadratische Teilgebiete
numSD = N^2; % Anzahl Teilgebiete
xyLim = [0,1]; % Gebiet: Einheitsquadrat

[vert,tri] = genMeshSquare(N,n); % Erstelle Knoten- und Elementliste
numVert=size(vert,1);   numTri=size(tri,1); % Anzahl Knoten und Dreiecke
% Erstelle Knoten- und Elementlisten pro Teilgebiet und logischer Vektor,
% welche Dreiecke in welchem TG enthalten sind
[vert__sd,tri__sd,l2g__sd,logicalTri__sd] = meshPartSquare(N,vert,tri); 

% Markiere Dirichletknoten in logischem Vektor
dirichlet = or(ismember(vert(:,1),xyLim), ismember(vert(:,2),xyLim)); 

%% PDE
f = @(vert,y) ones(size(vert));   % Rechte Seite der DGL

%% Definiere Kanal-Koeffizientenfunktion
% Definiere Bereich des Kanals
xMin=14/30; xMax=16/30;
yMin=3/30;  yMax=27/30;

% Definiere rho im Kanal und außerhalb des Kanals
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

%% Aufstellen der Referenzloesung
% Als Referenzloesung dient die Loesung des global assemblierten Sysmtems
% mit Backslash-Operator
[K,~,b] = assemble(tri,vert,1,f,rhoTri);
K_II = K(~dirichlet,~dirichlet);
b_I = b(~dirichlet);

u_ref = zeros(size(vert,1),1);
u_ref(~dirichlet) = K_II\b_I;

%% Loesen des Systems mit FETI-DP fuer versch. VK & Plots
diffs = cell(length(VK_vec),1);
iters = cell(length(VK_vec),1);
kappa_ests = cell(length(VK_vec),1);
termCond = cell(length(VK_vec),1);

fig_VK_comp_solution = figure("Name","Loesungen fuer verschiedene Vorkonditionierer");
fig_VK_comp_termCond = figure("Name",sprintf("Verlauf des %s Residuums fuer verschiedene Vorkonditionierer",append(resid{1},"en")));
tiledlayout('flow')
for vk_ind = 1:length(VK_vec) %Iteriere uber VK
    VK = VK_vec{vk_ind};
    % Loesen des Systems mit FETI-DP mit dem entsprechenden VK
    [cu,u_FETIDP_glob,~,iters{vk_ind},kappa_ests{vk_ind},termCond{vk_ind}] = fetidp_defl_bal(vert__sd,tri__sd,l2g__sd,...
                                                                             f,dirichlet,VK,rhoTriSD,maxRhoVert,maxRhoVertSD, ...
                                                                             tol,x0,resid);
                                                 
    diffs{vk_ind} = norm(u_FETIDP_glob-u_ref); % Abweichung der Loesung von der Referenzloesung
    
    % Plotten der finalen Loesung
    figure(fig_VK_comp_solution)
    nexttile
    hold on
    for sd = 1:length(tri__sd)
        trisurf(tri__sd{sd},vert__sd{sd}(:,1),vert__sd{sd}(:,2),cu{sd});
    end
    xlabel("x"); ylabel("y"); zlabel("z");
    title(sprintf("Finale Loesung: %s-VK",VK));
    view(3)
    hold off
    
    % Plotten des Verlaufs des relativen Residuums der Abbruchbedingung von PCG
    figure(fig_VK_comp_termCond)
    nexttile
    hold on
    plot(1:iters{vk_ind},termCond{vk_ind});
    xlabel("Iteration"); ylabel("Relatives Residuum");
    if strcmp('vorkonditioniert',resid) && strcmp('Dirichlet',VK)
        xlim([0,iters{vk_ind}])
        ylim([0 termCond{vk_ind}(1)+0.1])
    end
    title(sprintf("%s Residuum: %s-VK",append(resid{1},"es"),VK));
    view(2)
    hold off   
end

%% Ergebnistabelle
rowNames = ["Anzahl Iterationen","Konditionszahl","Abweichung von Referenzloesung"];
T_results = cell2table([iters';kappa_ests';diffs'],"RowNames",rowNames,"VariableNames",VK_vec)