clear; clc;
addpath('libs')

%% Definiere Vorkonditionierer
VK_vec = {'Deflation',...
          };
Triangulierung_vec = {'Traingulierung i)',...
                      'Traingulierung ii)',...
                      'Traingulierung ii)',...
                      };

cTOL_vec = {'TOL=10',...
            'TOL=10^3',...
            'TOL=10^6',...
            'TOL=10^10',...
              };

ckappa_ests = cell(length(Triangulierung_vec),length(cTOL_vec));


%% Parameter fuer PCG
x0 = @(dim) zeros(dim,1); % Startwert
tol = 10^(-7); % Toleranz
% Residuum fuer die Abbruchbedingung
resid = {'vorkonditioniert'}; 
%resid = {'nicht-vorkonditioniert'};

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
% Definiere Bereich des Kanals für i)
xMin=14/30; xMax=16/30;
yMin=3/30;  yMax=27/30;

% Definiere rhoMax Teilgebiete für ii)
rhoMax_sd = [2,3,14,17,24,25];
rhoMax = 10^6;
rhoMin = 1;
    
% Plot-Auswahl
plot_grid = true;
% Definiere Koeffizient auf den Elementen (und teilgebietsweise);
% maximalen Koeffizienten pro Knoten (und teilgebietsweise)
crhoTri = cell(length(Triangulierung_vec),1);
crhoTriSD = cell(length(Triangulierung_vec),1);
cmaxRhoVert = cell(length(Triangulierung_vec),1);
cmaxRhoVertSD = cell(length(Triangulierung_vec),1);


[crhoTri{1},crhoTriSD{1},cmaxRhoVert{1},cmaxRhoVertSD{1}] = coefficient_1(xMin,xMax,yMin,yMax,rhoMax,rhoMin,vert,tri,numVert,numTri,numSD,logicalTri__sd,N,plot_grid);

[crhoTri{2},crhoTriSD{2},cmaxRhoVert{2},cmaxRhoVertSD{2}] = coefficient_2ii(rhoMax,rhoMin,rhoMax_sd,l2g__sd,tri,vert,numVert,numTri,numSD,logicalTri__sd,N,plot_grid);

[crhoTri{3},crhoTriSD{3},cmaxRhoVert{3},cmaxRhoVertSD{3}] = coefficient_2iii(rhoMax,rhoMin,tri,vert,numVert,numTri,numSD,logicalTri__sd,N,plot_grid);

for i = 1 : length(Triangulierung_vec)
    
  %% Aufstellen der Referenzloesung
  % Als Referenzloesung dient die Loesung des global assemblierten Sysmtems
  % mit Backslash-Operator
  [K,~,b] = assemble(tri,vert,1,f,crhoTri{i});
  K_II = K(~dirichlet,~dirichlet);
  b_I = b(~dirichlet);
        
  u_ref = zeros(size(vert,1),1);
  u_ref(~dirichlet) = K_II\b_I;
        
  %% Loesen des Systems mit FETI-DP fuer versch. VK
    
  TOL_vec = [10,10^3,10^6,10^10]; 
  for t = 1 : length(TOL_vec)
%     fig_VK_comp = figure("Name","Loesungen fuer verschiedene Vorkonditionierer");
%     tiledlayout('flow')
      VK = VK_vec{1};
                
      [cu,u_FETIDP_glob,~,~,ckappa_ests{i}{t}] = fetidp_constraint(TOL_vec(t),vert__sd,tri__sd,l2g__sd,f,...
                                                                 dirichlet,VK,'adaptive',crhoTriSD{i},...
                                                                 cmaxRhoVert{i},cmaxRhoVertSD{i},tol,x0,resid);
%     figure(fig_VK_comp)
%     nexttile
%     hold on
%     for sd = 1:length(tri__sd)
%         trisurf(tri__sd{sd},vert__sd{sd}(:,1),vert__sd{sd}(:,2),cu{sd});
%     end
%     xlabel("x"); ylabel("y"); zlabel("z");
%     title(sprintf("Finale Loesung: %s-VK",VK));
%     view(3)
%     hold off
  end
end
        
%% Ergebnistabelle
rowNames = ["Konditionszahl i)","Konditionszahl ii)","Konditionszahl iii)"];
fprintf('Deflation-Vorkonditionierer \n')
fprintf('RhoMax= %g \n',rhoMax)
fprintf('Konditionszahl fuer unterschiedliche TOL und Triangulierungen: \n')
T_results = cell2table([ckappa_ests{1};ckappa_ests{2};ckappa_ests{3}],"RowNames",rowNames,"VariableNames",cTOL_vec);
disp(T_results)