function [K,M,b] = assemble(tri,x,order,f,rhoTri)
% Input: tri, x: Element- und Knotenliste
% Input: order: Ordnung der Elementdiskretisierung
% Input: f: rechte Seite der DGL
% Input: rhoTri: Koeffizient pro Element

% Output: K: Elementsteifigkeitsmatrix
% Output: M: Massenmatrix
% Output: b: Lastvektor

%% Definiere nuetzliche Parameter
numElements = size(tri,1);          % Anzahl Elemente
numBaseFun = (order+2)*(order+1)/2; % Anzahl Basisfunktionen
nBF2 = numBaseFun^2;                % Anzahl Elemente der lokalen Steifigkeitsmatrix

%% Initialisiere Matrizen
K_val = zeros(nBF2*numElements,1); 
M_val = K_val;
iIndex = K_val; 
jIndex = K_val;
b = zeros(size(x,1),1);

%% Initialisiere Basisfunktionen und Quadraturformeln
[phi,d_phi] = baseFun(order);
% phi = load(sprintf("p%ibaseFun.mat",order)).baseFun;
% d_phi = load(sprintf("p%ibaseFun.mat",order)).baseDer;
if order == 1
%     quad_low = load("ueb10_programm_daten.mat").quadratur_P5;
    quad_low = load("order1_quad.mat").quadratur_P1;
%     quad_high = load("order1_quad.mat").quadratur_P1;
    quad_high = load("ueb10_programm_daten.mat").quadratur_P5;
elseif order == 2
    quad_low = load("ueb10_programm_daten.mat").quadratur_P2;
    quad_high = load("ueb10_programm_daten.mat").quadratur_P5;
end

%% Assemblierung
for i = 1:size(tri,1)
    [B,d] = aff_map(x,tri(i,:)); % Bestimme affin-lineare Abbildung
    
    % Stelle Elementsteifigekitsmatrix, -massenmatrix und -lastvektor auf
    [K_T,M_T,b_T] = getMatrices(B,d,f,phi,d_phi,quad_low,quad_high);
    K_T = rhoTri(i)*K_T; % Koeffizienten in Elementsteifigkeitsmatrix einbinden
    
    % sparse()-Indizierung
    iIndex((i-1)*nBF2+1:i*nBF2) = reshape(repmat(tri(i,:),numBaseFun,1),nBF2,1);
    jIndex((i-1)*nBF2+1:i*nBF2) = repmat(tri(i,:),1,numBaseFun); 
    K_val((i-1)*nBF2+1:i*nBF2) = reshape(K_T,nBF2,1);
    M_val((i-1)*nBF2+1:i*nBF2) = reshape(M_T,nBF2,1);
    b(tri(i,:)) = b(tri(i,:)) + b_T;
end
n = size(x,1);
K = sparse(iIndex,jIndex,K_val,n,n);
M = sparse(iIndex,jIndex,M_val,n,n);
end

%% Hilfsfunktion
function [K,M,b] = getMatrices(B,d,f,phi,d_phi,quad_low,quad_high)

numBaseFunc = length(phi);
K = zeros(numBaseFunc);
M = zeros(numBaseFunc);
b = zeros(numBaseFunc,1);

detb = abs(det(B));
invb = B\eye(size(B));

for i = 1:numBaseFunc
    for j = 1:numBaseFunc
        %% Matrix K_T
        x = quad_low.knoten(:,1);
        y = quad_low.knoten(:,2);
        temp = dot([d_phi{1,i}(x,y), d_phi{2,i}(x,y)]*invb,...
                   [d_phi{1,j}(x,y), d_phi{2,j}(x,y)]*invb,2);
        K(i,j) = detb* sum(quad_low.gewichte .* temp);
    end  
    
    %% Vektor b_T
    v_ref = quad_high.knoten';
    v = B*v_ref +d;
    temp = f(v(1,:),v(2,:)).*phi{i}(v_ref(1,:),v_ref(2,:));
    b(i) = detb* sum(quad_high.gewichte'.*temp);
end
end

%% Affin-lineare Funktion
function [B,d] = aff_map(coords,x)
a1 = coords(x(1),1:2)';
a2 = coords(x(2),1:2)';
a3 = coords(x(3),1:2)';

d = a1;
B = [a2-a1 , a3-a1];
end
