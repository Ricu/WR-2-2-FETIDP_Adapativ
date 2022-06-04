clear; clc;
addpath('libs')

%% Definiere Gitterfeinheiten
n_vec = [10,20];%,30,40];    % 2*n^2 Elemente pro Teilgebiet

%% Definiere zu testende Vorkonditionierer
% VK_vec = {'Deflation',...
%     'Balancing'...
%     };
VK_vec = {'Identitaet',...
          'Dirichlet'};

%% Definiere constraint-Typ
constraint_type_vec = {'non-adaptive',...
    'adaptive'...
    };

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
N = 3;                  % Partition in NxN quadratische Teilgebiete
h_vec = 1./(N.*n_vec);  % Gitterfeinheit
numSD = N^2;            % Anzahl Teilgebiete
xyLim = [0,1];          % Gebiet: Einheitsquadrat

%% Definiere Kanal-Koeffizientenfunktion
% Definiere Bereich des Kanals
xCanalLim = [14/30,16/30];
yCanalLim  = [3/30,27/30];

% Definiere rho im Kanal und außerhalb des Kanals
rhoMax = 10^6;
rhoMin = 1;

color = {'blue','red'};

for constraint_type_ind = 1:length(constraint_type_vec)
    constraint_type = constraint_type_vec{constraint_type_ind};
    
    fig = figure("Name",sprintf("Gitterkonvergenz: %s",constraint_type));
    
    diffs = cell(length(VK_vec),length(n_vec));
    iters = cell(length(VK_vec),length(n_vec));
    kappa_ests = cell(length(VK_vec),length(n_vec));
    
    for n_ind = 1:length(n_vec)
        n = n_vec(n_ind);
        [vert,tri] = genMeshSquare(N,n);            % Erstelle Knoten- und Elementliste
        numVert=size(vert,1);   numTri=size(tri,1); % Anzahl Knoten und Dreiecke
        
        % Erstelle Knoten- und Elementlisten pro Teilgebiet und logischen Vektor,
        % welche Dreiecke in welchem TG enthalten sind
        [vert__sd,tri__sd,l2g__sd,logicalTri__sd] = meshPartSquare(N,vert,tri);
        
        % Markiere Dirichletknoten in logischem Vektor
        dirichlet = or(ismember(vert(:,1),xyLim), ismember(vert(:,2),xyLim));
        
        % Definiere Koeffizient auf den Elementen (und teilgebietsweise);
        % Maximalen Koeffizienten pro Knoten (und teilgebietsweise)
        plot_grid = false;   % Auswahl: Plotten der Triangulierung mit Kanal-Koeffizientenfunktion
        [rhoTri,rhoTriSD,maxRhoVert,maxRhoVertSD] = coefficient_1(xCanalLim,yCanalLim,rhoMax,rhoMin,vert,tri,logicalTri__sd,plot_grid);

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
        
        %% Loesen des Systems mit FETI-DP fuer versch. VK
        for vk_ind = 1:length(VK_vec) %Iteriere ueber VK
            VK = VK_vec{vk_ind};
            
            % Loesen des Systems mit FETI-DP mit entsprechendem VK
            if strcmp('adaptive',constraint_type)
                TOL = 100;
                pc_param = struct('VK',VK,'constraint_type',constraint_type,'adaptiveTol',TOL);
            else
                pc_param = struct('VK',VK,'constraint_type',constraint_type);
            end 
            [~,u_FETIDP_glob,~,iters{vk_ind,n_ind},kappa_ests{vk_ind,n_ind},~] = fetidp(grid_struct,f,pc_param,rho_struct,pcg_param,false);
            
            % Abweichung der Loesung von der Referenzloesung
            diffs{vk_ind,n_ind} = norm(u_FETIDP_glob-u_ref);
        end
    end
    
    diffs = cell2mat(diffs);
    iters = cell2mat(iters);
    kappa_ests = cell2mat(kappa_ests);

    %% Plotten der Ergebnisse
    for vk_ind = 1:length(VK_vec) 
        subplot(1,3,1)
        hold on
        plot(h_vec,iters(vk_ind,:),color{vk_ind},'DisplayName',VK_vec{vk_ind})
        hold off
        
        subplot(1,3,2)
        hold on
        plot(h_vec,kappa_ests(vk_ind,:),color{vk_ind},'DisplayName',VK_vec{vk_ind})
        hold off
        
        subplot(1,3,3)
        hold on
        plot(h_vec,diffs(vk_ind,:),color{vk_ind},'DisplayName',VK_vec{vk_ind})
        hold off
    end
    subplot(1,3,1)
    xlabel("Gitterfeinheit",FontWeight='bold'); ylabel("Anzahl Iterationen",FontWeight='bold');
    title("Anzahl Iterationen");
    legend
    subplot(1,3,2)
    xlabel("Gitterfeinheit",FontWeight='bold'); ylabel("Konditionszahl",FontWeight='bold');
    title("Konditionszahl");
    legend
    subplot(1,3,3)
    xlabel("Gitterfeinheit",FontWeight='bold'); ylabel("Abweichung von Referenzloesung",FontWeight='bold');
    title("Abweichung von Referenzloesung");
    legend
end
