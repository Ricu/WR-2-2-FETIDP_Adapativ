clear; clc;
addpath('libs')
plot_iteration = 0;
plot_solution = 0;
plot_grid = 0;

%% Parameter fuer PCG
x0 = @(dim) zeros(dim,1);   % Startwert
tol = 10^(-8);              % Toleranz
% Residuum fuer die Abbruchbedingung
resid_type = {'vorkonditioniert'}; 

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

%% Definiere Kanal-Koeffizientenfunktion
% Definiere Bereich des Kanals
xCanalLim = [14/30,16/30];
yCanalLim  = [3/30,27/30];

%% Deflation
VK = 'Deflation';
constraint_type = 'adaptive';
% TOL_vec = 10.^(0);
TOL_vec = [1.1,10.^(1:7)];
rhoMax = 10^6;
rhoMin= 1;

% Definiere Koeffizient auf den Elementen (und teilgebietsweise);
% maximalen Koeffizienten pro Knoten (und teilgebietsweise)
[rhoTri,rhoTriSD,maxRhoVert,maxRhoVertSD] = coefficient_1(xCanalLim,yCanalLim,rhoMax,rhoMin,vert,tri,logicalTri__sd,plot_grid);
% [rhoTri,rhoTriSD,maxRhoVert,maxRhoVertSD] = coefficient_2(rhoMax,rhoMin,[2,3,14,17,24,25],vert,tri,logicalTri__sd,plot_grid);
% [rhoTri,rhoTriSD,maxRhoVert,maxRhoVertSD] = coefficient_3(rhoMax,rhoMin,vert,tri,logicalTri__sd,0.25,0,plot_grid);

diffs = cell(length(TOL_vec),1);
iters = cell(length(TOL_vec),1);
kappa_ests = cell(length(TOL_vec),1);
cond_vec = zeros(length(TOL_vec),1);
if plot_solution
    fig_solutions = figure("Name","Loesungen");
    tiledlayout('flow')
end
fig_ew = figure("Name", "Darstellung der Top 50 Eigenwerte des vorkonditionierten Systems");
tiledlayout('flow','TileSpacing','tight')


% Als Referenzloesung dient die Loesung des global assemblierten Sysmtems
% mit Backslash-Operator
[K,~,b] = assemble(tri,vert,1,f,rhoTri);
K_II = K(~dirichlet,~dirichlet);
b_I = b(~dirichlet);

u_ref = zeros(size(vert,1),1);
u_ref(~dirichlet) = K_II\b_I;

for tolInd = 1 : length(TOL_vec)
    TOL = TOL_vec(tolInd);

    % Uebergabe strukturen erstellen
    rho_struct = struct('rhoTriSD',{rhoTriSD},'maxRhoVert',{maxRhoVert},'maxRhoVertSD',{maxRhoVertSD});
    grid_struct = struct('vert__sd',{vert__sd},'tri__sd',{tri__sd},'l2g__sd',{l2g__sd},'dirichlet',{dirichlet});
    pc_param = struct('VK',VK,'constraint_type',constraint_type,'adaptiveTol',TOL);
    pcg_param = struct('tol', tol, 'x0',x0, 'resid_type',resid_type);

    [cu,u_FETIDP_glob,~,iters{tolInd},kappa_ests{tolInd},~,preconditioned_system] = fetidp(grid_struct,f,pc_param,rho_struct,pcg_param,plot_iteration);
    
    % Abweichung der Loesung von der Referenzloesung
    diffs{tolInd} = norm(u_FETIDP_glob-u_ref);

    if plot_solution
        figure(fig_solutions)
        nexttile
        hold on
        for sd = 1:length(tri__sd)
            trisurf(tri__sd{sd},vert__sd{sd}(:,1),vert__sd{sd}(:,2),cu{sd});
        end
        xlabel("x"); ylabel("y"); zlabel("z");
        title(sprintf("Finale Loesung: TOL = %i",TOL));
        view(3)
        hold off
    end

    % Analysiere EW und Kondition von invmF
    ew = abs(eig(preconditioned_system));
    ew = sort(ew,'descend');
    topEW = ew(1:min(length(ew),50));
    cond_vec(tolInd) = topEW(1); % Resymmetrisiere
    % Plot der 50 groessten EW
    figure(fig_ew)
    nexttile
    scatter(1:50,topEW,18,'k','filled')
    title(sprintf("TOL = %.1f",TOL));

end

fig_condition = figure("Name","Kondition in Abhängigkeit von TOL");
loglog(TOL_vec,cond_vec);
xlabel(sprintf("TOL")); ylabel("Konditionszahl");

% Ergebnistabelle
rowNames = ["Anzahl Iterationen","Konditionszahl","Abweichung von Referenzloesung"];
variableNames = cellfun(@num2str,num2cell(TOL_vec),'UniformOutput',false);
fprintf('RhoCanal: %g \n',rhoMax)
fprintf('TOL zur Auswahl der EW: %g \n',TOL)
T_results = cell2table([iters';kappa_ests';diffs'],"RowNames",rowNames,"VariableNames",variableNames);
disp(T_results)