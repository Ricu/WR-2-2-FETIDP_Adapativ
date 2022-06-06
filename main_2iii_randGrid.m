clear; clc;
addpath('libs')

%% Definiere Vorkonditionierer
VK_vec = {%'Dirichlet',...
          'Deflation',...
          %'Balancing',...
          };

RandomPercentage_vec = {'random_percentage = 0.0',...
                  'random_percentage = 0.1',...
                  'random_percentage = 0.2',...
                  'random_percentage = 0.3',...
                  'random_percentage = 0.4',...
                  'random_percentage = 0.5',...
                  'random_percentage = 0.6',...
                  'random_percentage = 0.7',...
                  'random_percentage = 0.8',...
                  'random_percentage = 0.9',...
                  'random_percentage = 1.0',...
                   };

constraint_type = 'adaptive';

diffs = cell(length(VK_vec),1);
iters = cell(length(VK_vec),1);
kappa_ests = cell(length(VK_vec),1);

%% Parameter fuer PCG
x0 = @(dim) zeros(dim,1); % Startwert
tol = 10^(-8); % Toleranz
% Residuum fuer die Abbruchbedingung
resid_type = {'vorkonditioniert'}; 
% resid_type = {'nicht-vorkonditioniert'};

%% Erstelle das Gitter
n = 10; % 2*n^2 Elemente pro Teilgebiet
N = 5;  % Partition in NxN quadratische Teilgebiete
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

%% Definiere Koeffizientenfunktion
% Definiere maximale unbd minimale Koeffizientenfunktion
% rhoMax zum Vergleich der rhoMin/rhoMax zur Konditionszahl
%rhoMax_vec = [1,10,10^2,10^3,10^4,10^5,10^6,10^7,10^8];
rhoMax_vec = 10^6;
rhoMin = 1;

% Plot-Auswahl
plot_grid = false;
random_state = 0;
random_percentage = 0:0.1:1;
for rand = 1 : length(random_percentage)
    % Definiere Koeffizient auf den Elementen (und teilgebietsweise);
    % maximalen Koeffizienten pro Knoten (und teilgebietsweise)
    [rhoTri,rhoTriSD,maxRhoVert,maxRhoVertSD] = coefficient_2iii(rhoMax_vec,rhoMin,vert,tri,logicalTri__sd,random_percentage(rand),random_state,plot_grid);
                                                                 
    %% Aufstellen der Referenzloesung
    % Als Referenzloesung dient die Loesung des global assemblierten Sysmtems
    % mit Backslash-Operator
    [K,~,b] = assemble(tri,vert,1,f,rhoTri);
    K_II = K(~dirichlet,~dirichlet);
    b_I = b(~dirichlet);
    
    u_ref = zeros(size(vert,1),1);
    u_ref(~dirichlet) = K_II\b_I;
    
    %% Loesen des Systems mit FETI-DP fuer versch. VK
    fig_VK_comp = figure("Name","Loesungen fuer verschiedene Vorkonditionierer");
    tiledlayout('flow')
    TOL = 100;
    for vk_ind = 1:length(VK_vec)
        VK = VK_vec{vk_ind};
    
        rho_struct = struct('rhoTriSD',{rhoTriSD},'maxRhoVert',{maxRhoVert},'maxRhoVertSD',{maxRhoVertSD});
        grid_struct = struct('vert__sd',{vert__sd},'tri__sd',{tri__sd},'l2g__sd',{l2g__sd},'dirichlet',{dirichlet});
        pc_param = struct('VK',VK,'constraint_type',constraint_type,'adaptiveTol',TOL);
        pcg_param = struct('tol', tol, 'x0',x0, 'resid_type',resid_type);
    
        [cu,u_FETIDP_glob,~,iters{vk_ind}{rand},kappa_ests{vk_ind}{rand}] = fetidp_constraint(grid_struct,f,pc_param,rho_struct,pcg_param,false);
        diffs{vk_ind}{rand} = norm(u_FETIDP_glob-u_ref);
    
%         figure(fig_VK_comp)
%         nexttile
%         hold on
%         for sd = 1:length(tri__sd)
%             trisurf(tri__sd{sd},vert__sd{sd}(:,1),vert__sd{sd}(:,2),cu{sd});
%         end
%         xlabel("x"); ylabel("y"); zlabel("z");
%         title(sprintf("Finale Loesung: %s-VK",VK));
%         view(3)
%         hold off
    end
end
%% Ergebnistabelle
rowNames = ["Anzahl Iterationen","Konditionszahl","Abweichung von Referenzloesung"];
T_results = cell2table([iters{1,:};kappa_ests{1,:};diffs{1,:}],"RowNames",rowNames,"VariableNames",RandomPercentage_vec)

figure("Name","Konditionszahl unter verschiedenen random_percentage");
tiledlayout('flow')
    nexttile
    plot(random_percentage,cell2mat(iters{1,:}))
    title("Anzahl Iterationen");
    nexttile
    plot(random_percentage,cell2mat(kappa_ests{1,:}))
    title("Konditionszahl");
    nexttile
    plot(random_percentage,cell2mat(diffs{1,:}))
    title("Abweichung von der Musterlösung");
