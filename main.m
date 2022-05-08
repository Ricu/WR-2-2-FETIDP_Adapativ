clear; clc;
addpath('libs')

%% Create grid
n = 10; % 2*n^2 Elemente pro Teilgebiet
N = 3;  % Partition in NxN quadratische Teilgebiete
numSD = N^2; % Anzahl Teilgebiete
% Gebiet: Einheitsquadrat
vertLim = [0,1];
yLim = [0,1];

[vert,tri] = genMeshSquare(N,n); % Erstelle Knoten- und Elementliste
numVert=size(vert,1);   numTri=size(tri,1); % Anzahl Knoten und Dreiecke
[vert__sd,tri__sd,l2g__sd] = meshPartSquare(N,vert,tri); % Erstelle Knoten- und Elementlisten pro Teilgebiet
% Dirichletrand fehlt in Aufgabenstellung?!
dirichlet = or(ismember(vert(:,1),vertLim), ismember(vert(:,2),yLim)); % Dirichletknoten, logischer Vektor

%% PDE
f = @(vert,y) ones(size(vert));   % Rechte Seite der DGL
% Definiere Koeffizientenfunktion
pho = ones(numTri,1);
% Definiere Kanal
xMin=14/30; xMax=16/30;
yMin=3/30;  yMax=27/30;
indVert = (xMin <= vert(:,1) & vert(:,1) <= xMax & yMin <= vert(:,2) & vert(:,2) <= yMax); % Definiere Kanal, logischer Vektor
canalVert = 1:numVert;
canalVert = canalVert(indVert); % Finde Knotennummern der Knoten, die im Kanal liegen
for i=1:numTri % Iteriere ueber die Elemente
    if ismember(tri(i,:),canalVert) % Alle Knoten des Elements liegen im Kanal
        pho(i)=10^6;    % Im Kanal entspricht die Koeffizientenfunktion 10^6
    end
end

canalElements = pho > 1;
figure()

patch('vertices',vert,'faces',tri,'edgecol','k','facecol',[1,1,1]); hold on; axis equal tight;
patch('vertices',vert,'faces',tri(canalElements,:),'edgecol','k','facecol',[.8,.9,1]);

[cu,u_FETIDP_glob] = fetidp(vert,vert__sd,tri__sd,l2g__sd,f,dirichlet,true);

                
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





