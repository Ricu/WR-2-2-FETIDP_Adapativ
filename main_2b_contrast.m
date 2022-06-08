clear; clc;
addpath('libs')

%% Definiere zu testenden Vorkonditionierer
VK = 'Dirichlet';

%% Definiere constraint-Typ
constraint_type = 'adaptive';

TOL= 100;   % Toleranz zur Auswahl der Eigenwerte

% Structure fuer PC-Parameter
pc_param = struct('VK',VK,'constraint_type',constraint_type,'adaptiveTol',TOL);
    
%% Initialisiere Parameter fuer PCG
x0 = @(dim) zeros(dim,1);    % Startvektor
tol = 10^(-8);               % Toleranz fuer die Abbruchbedingung

% Residuum fuer die Abbruchbedingung
resid_type = {'vorkonditioniert'};

% Structure fuer PCG-Parameter
pcg_param = struct('tol', tol, 'x0',x0, 'resid_type',resid_type);

%% PDE
f = @(vert,y) ones(size(vert));   % Rechte Seite der DGL

%% Erstelle das Gitter
n = 10;         % 2*n^2 Elemente pro Teilgebiet
N = 5;          % Partition in NxN quadratische Teilgebiete
numSD = N^2;    % Anzahl Teilgebiete
xyLim = [0,1];  % Gebiet: Einheitsquadrat

[vert,tri] = genMeshSquare(N,n);            % Erstelle Knoten- und Elementliste
numVert=size(vert,1);   numTri=size(tri,1); % Anzahl Knoten und Dreiecke

% Erstelle Knoten- und Elementlisten pro Teilgebiet und logischen Vektor,
% welche Dreiecke in welchem TG enthalten sind
[vert__sd,tri__sd,l2g__sd,logicalTri__sd] = meshPartSquare(N,vert,tri);

% Markiere Dirichletknoten in logischem Vektor
dirichlet = or(ismember(vert(:,1),xyLim), ismember(vert(:,2),xyLim));

% Structure fuer grid-Variablen
grid_struct = struct('vert__sd',{vert__sd},'tri__sd',{tri__sd},'l2g__sd',{l2g__sd},'dirichlet',{dirichlet});


%% Definiere Kanal-Koeffizientenfunktion
% Definiere Bereich des Kanals
xCanalLim = [14/30,16/30];
yCanalLim  = [3/30,27/30];

% Definiere maximale unbd minimale Koeffizientenfunktion
% Variiere rhoMax, um den Kontrast der Koeffizientenfunktion zu veraendern
rhoMax_vec = 10.^(1:8);
rhoMin = 1;

plot_grid = false;   % Auswahl: Plotten der Triangulierung mit Kanal-Koeffizientenfunktion

%% Loesen des Systems mit FETI-DP fuer versch.Kontraste & Ersetllen der Ergebnisplots
diffs = cell(length(rhoMax_vec),1);
iters = cell(length(rhoMax_vec),1);
kappa_ests = cell(length(rhoMax_vec),1);
cond_vec = zeros(length(rhoMax_vec),1);

plot_solution = false;  % Auswahl: Plotten der Loesung
if plot_solution
    fig_solutions = figure("Name","Loesungen");
    tiledlayout('flow')
end

fig_ew = figure("Name", "Darstellung der Top 50 Eigenwerte des vorkonditionierten Systems");
tiledlayout('flow','TileSpacing','tight')

plot_iteration = false; % Auswahl: Plotten der Loesung nach den ersten Iterationen von PCG

for rhoInd = 1: length(rhoMax_vec)
    rhoMax = rhoMax_vec(rhoInd);
    
    % Definiere Koeffizient auf den Elementen (und teilgebietsweise);
    % maximalen Koeffizienten pro Knoten (und teilgebietsweise)
    [rhoTri,rhoTriSD,maxRhoVert,maxRhoVertSD] = coefficient_1(xCanalLim,yCanalLim,rhoMax,rhoMin,vert,tri,logicalTri__sd,plot_grid);
    %[rhoTri,rhoTriSD,maxRhoVert,maxRhoVertSD] = coefficient_2(rhoMax,rhoMin,[2,3,14,17,24,25],vert,tri,logicalTri__sd,plot_grid);
    %[rhoTri,rhoTriSD,maxRhoVert,maxRhoVertSD] = coefficient_3(rhoMax,rhoMin,vert,tri,logicalTri__sd,0.25,0,plot_grid);
    
    % Structure fuer rho-Variablen
    rho_struct = struct('rhoTriSD',{rhoTriSD},'maxRhoVert',{maxRhoVert},'maxRhoVertSD',{maxRhoVertSD});
    
    % Aufstellen der Referenzloesung
    % Als Referenzloesung dient die Loesung des global assemblierten Sysmtems
    % mit Backslash-Operator
    [K,~,b] = assemble(tri,vert,1,f,rhoTri);
    K_II = K(~dirichlet,~dirichlet);
    b_I = b(~dirichlet);
    u_ref = zeros(size(vert,1),1);
    u_ref(~dirichlet) = K_II\b_I;
    
    % Loesen des Systems mit FETI-DP mit entsprechendem VK
    [cu,u_FETIDP_glob,~,iters{rhoInd},kappa_ests{rhoInd},~,preconditioned_system] = fetidp(grid_struct,f,pc_param,rho_struct,pcg_param,plot_iteration);
    
    % Abweichung der Loesung von der Referenzloesung
    diffs{rhoInd} = norm(u_FETIDP_glob-u_ref);
    
    % Plotten der finalen Loesung
    if plot_solution
        figure(fig_solutions)
        nexttile
        hold on
        for sd = 1:length(tri__sd)
            trisurf(tri__sd{sd},vert__sd{sd}(:,1),vert__sd{sd}(:,2),cu{sd});
        end
        xlabel("x"); ylabel("y"); zlabel("z");
        title(sprintf("Finale Loesung: rhoMax = %i",rhoMax));
        view(3)
        hold off
    end

    %% Analysiere EW und Kondition von invmF
    % Bestimme die 50 groessten EW von invmF
    ew = abs(eig(preconditioned_system));
    ew = sort(ew,'descend');
    topEW = ew(1:min(length(ew),50));
    
    % Der groesste EW dient als Approximation der Konditionszahl
    cond_vec(rhoInd) = topEW(1);
    
    % Plot der 50 groessten EW
    figure(fig_ew)
    nexttile
    scatter(1:50,topEW,18,'k','filled')
    xlabel("Rang der Eigenwerte"); ylabel("Groesse der Eigenwerte");
    title(sprintf("\\rho_{max}/\\rho_{min} = %.0e",rhoMax/rhoMin));
end

% Plot der Konditionszahl in Abhaengigkeit vom Kontrast der Koeffizientenfunktion
fig_condition = figure("Name","Kondition in Abhaengigkeit vom Kontrast der Koeffizientenfunktion");
loglog(rhoMax_vec./rhoMin,cond_vec);
xlabel("$\displaystyle\frac{\rho_{max}}{\rho_{min}}$","interpreter","latex"); ylabel("Konditionszahl");
title("Konditionszahl in Abhaengigkeit vom Kontrast der Koeffizientenfunktion");

%% Ergebnistabelle
rowNames = ["Anzahl Iterationen","Konditionszahl","Abweichung von Referenzloesung"];
variableNames = cellfun(@num2str,num2cell(rhoMax_vec),'UniformOutput',false);
T_results = cell2table([iters';kappa_ests';diffs'],"RowNames",rowNames,"VariableNames",variableNames);
disp(T_results)