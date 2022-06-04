clear; clc;
addpath('libs')

%% Definiere zu testende Vorkonditionierer
VK_vec = {'Identitaet',...
          'Dirichlet',...
          'Deflation',...
          'Balancing'...
          };

%% Definiere constraint-Typ
constraint_type = 'non-adaptive';
      
%% Initialisiere Parameter fuer PCG
x0 = @(dim) zeros(dim,1);    % Startvektor
tol = 10^(-8);               % Toleranz fuer die Abbruchbedingung

% Residuum fuer die Abbruchbedingung
resid_type = {'vorkonditioniert'}; 
% resid_type = {'nicht-vorkonditioniert'};
% resid_type = {'nicht-vorkonditioniert,alternativ'}; % Alternative fuer Deflation

pcg_param = struct('tol', tol, 'x0',x0, 'resid_type',resid_type); % Structure fuer PCG-Parameter

%% PDE
f = @(vert,y) ones(size(vert));   % Rechte Seite der DGL

%% Erstelle das Gitter
n = 10;         % 2*n^2 Elemente pro Teilgebiet
N = 3;          % Partition in NxN quadratische Teilgebiete
numSD = N^2;    % Anzahl Teilgebiete
xyLim = [0,1];  % Gebiet: Einheitsquadrat

[vert,tri] = genMeshSquare(N,n);            % Erstelle Knoten- und Elementliste
numVert=size(vert,1);   numTri=size(tri,1); % Anzahl Knoten und Dreiecke

% Erstelle Knoten- und Elementlisten pro Teilgebiet und logischen Vektor,
% welche Dreiecke in welchem TG enthalten sind
[vert__sd,tri__sd,l2g__sd,logicalTri__sd] = meshPartSquare(N,vert,tri); 

% Markiere Dirichletknoten in logischem Vektor
dirichlet = or(ismember(vert(:,1),xyLim), ismember(vert(:,2),xyLim)); 

%% Definiere Kanal-Koeffizientenfunktion
% Definiere Bereich des Kanals
xCanalLim = [14/30,16/30];
yCanalLim  = [3/30,27/30];

% Definiere rho im Kanal und außerhalb des Kanals
rhoMax = 10^6;
rhoMin = 1;

% Definiere Koeffizient auf den Elementen (und teilgebietsweise);
% Maximalen Koeffizienten pro Knoten (und teilgebietsweise)
plot = true;   % Auswahl: Plotten der Triangulierung mit Kanal-Koeffizientenfunktion
[rhoTri,rhoTriSD,maxRhoVert,maxRhoVertSD] = coefficient_1(xCanalLim,yCanalLim,rhoMax, ...
                                                           rhoMin,vert,tri,logicalTri__sd,plot);
                                                      
% Structure fuer rho-Variablen
rho_struct = struct('rhoTriSD',{rhoTriSD},'maxRhoVert',{maxRhoVert},'maxRhoVertSD',{maxRhoVertSD}); 
% Structure fuer grid-Variablen
grid_struct = struct('vert__sd',{vert__sd},'tri__sd',{tri__sd},'l2g__sd',{l2g__sd},'dirichlet',{dirichlet});

%% Aufstellen der Referenzloesung
% Als Referenzloesung dient die Loesung des global assemblierten Sysmtems
% mit Backslash-Operator
[K,~,b] = assemble(tri,vert,1,f,rhoTri);
K_II = K(~dirichlet,~dirichlet);
b_I = b(~dirichlet);

u_ref = zeros(size(vert,1),1);
u_ref(~dirichlet) = K_II\b_I;

%% Loesen des Systems mit FETI-DP fuer versch. VK & Ersetllen der Ergebnisplots
diffs = cell(length(VK_vec),1);
iters = cell(length(VK_vec),1);
kappa_ests = cell(length(VK_vec),1);
residuals = cell(length(VK_vec),1);

fig_VK_comp_solution = figure("Name","Loesungen fuer verschiedene Vorkonditionierer");
fig_VK_comp_residuals = figure("Name",sprintf("Verlauf des %s Residuums fuer verschiedene Vorkonditionierer",append(resid_type{1},"en")));
tiledlayout('flow')

for vk_ind = 1:length(VK_vec) %Iteriere ueber VK
    VK = VK_vec{vk_ind};
    % Loesen des Systems mit FETI-DP mit entsprechendem VK
    pc_param = struct('VK',VK,'constraint_type',constraint_type);
    [cu,u_FETIDP_glob,~,iters{vk_ind},kappa_ests{vk_ind},residuals{vk_ind}] = fetidp(grid_struct,f,pc_param,rho_struct,pcg_param,true);
    
    % Abweichung der Loesung von der Referenzloesung
    diffs{vk_ind} = norm(u_FETIDP_glob-u_ref);
    
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
    
    % Plotten des Verlaufs des relativen Residuums
    figure(fig_VK_comp_residuals)
    nexttile
    
    semilogy(1:iters{vk_ind},residuals{vk_ind});
    hold on
    xlabel("Iteration"); ylabel("Relatives Residuum");
    title(sprintf("%s Residuum: %s-VK",append(resid_type{1},"es"),VK));
    view(2)
    hold off   
end

%% Ergebnistabelle
rowNames = ["Anzahl Iterationen","Konditionszahl","Abweichung von Referenzloesung"];
T_results = cell2table([iters';kappa_ests';diffs'],"RowNames",rowNames,"VariableNames",VK_vec);
disp(T_results)