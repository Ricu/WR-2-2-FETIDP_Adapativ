function [cEdges,adjSd__e,cPrimalNodesOfEdges] = findEdges(cPrimal,primal,cDual,dual,l2g__sd)
% Finde Kanten (alle Eckknoten sind primal).
%
% INPUT:
%    cPrimal: Cell-array: Liste (fuer jedes Teilgebiet) aller lokaler primaler
%             Knoten mit lokaler Nummerierung (cPrimal{i} darf ein logischer 
%             Vektor sein oder eine Nummernliste).
%    primal: Logischer Vektor: Globale Liste aller primaler Knoten.
%    cDual: Cell-array: Liste (fuer jedes Teilgebiet) aller lokaler dualer 
%           Knoten mit lokaler Nummerierung (cDual{i} darf ein logischer 
%           Vektor sein oder eine Nummernliste).
%    dual: Logischer Vektor: Globale Liste aller dualer Knoten.
%    l2g__sd: Cell-array: Local-2-global-map fuer KnotenIDs.
%
% OUTPUT:
%    cEdges: Cell-array: Liste aller Kanten; globale Knotennummerierung.
%    adjSd__e: Liste (n x 2): Zur (offenen) Kante benachbarte Teilgebiete.
%    cPrimalNodesOfEdges: Cell-array (n x 2): 
%       In cPrimalNodesOfEdges{i,j} steht zur i-ten Kante eine Liste aller 
%       inzidenten primalen Knoten des j=1. bzw. j=2. benachbarten 
%       Teilgebiets. Die Nummerierung ist lokal bzgl. des Teilgebiets.
%       
%
    %
    
    % Bestimme Teilgebietszugehoerigkeit dualer Knoten / Kantenknoten.
    dual_sd = findAdjacentSubdomainsOfDualNodes(cDual,l2g__sd,dual);
    
    % Sortiere Teilgebietsnummern bzgl. der Spalten, sodass
    %    dual_sd(:,1) < dual_sd(:,2).
    dual_sd = sort(dual_sd,2);
    
    % Sortiere Zeilen.
    % [ 1,2   <-- Kante 1
    %   1,2   <-- Kante 1
    %   1,2   <-- Kante 1
    %   1,3   <-- Kante 2
    %   1,3   <-- Kante 2
    %   ...
    %   3,4   <-- Kante m
    %   3,4]  <-- Kante m
    [dual_sd,sortmap] = sortrows(dual_sd);
    
    % Ziehe Zeilen voneinander ab (von oben nach unten).
    % [ 1,2   <-- Kante 1  --|
    %   1,2   <-- Kante 1  <-|  --|
    %   1,2   <-- Kante 1       <-|  --|
    %   1,3   <-- Kante 2  --|       <-|
    %   1,3   <-- Kante 2  <-|  --|
    %   ...                     <-|  --|
    %   3,4   <-- Kante m  --|       <-|
    %   3,4]  <-- Kante m  <-|
    % Ergibt:
    % [ 0,0
    %   0,0
    %   0,1
    %   0,0
    %   ...
    %   ?,?
    %   0,0]
    tmp = diff(dual_sd);
    
    % Alle Zeilen, die vollstaendig Null sind, werden auf 'false' 
    % gesetzt (true --> eindeutige Zeile).
    % [ 0,0  --> false
    %   0,0  --> false
    %   0,1  --> true
    %   0,0  --> false
    %   ...
    %   ?,?  --> ?
    %   0,0]  --> false
    tmp = any(tmp,2);
    
    % Die erste Zeile von 'dual_sd' ist per Definition eindeutig.
    % Fuege also 'true' vorne ein.
    ind = [true;
           tmp];
    
    % Fuege am Ende noch ein dummy-'true' ein; eine Art Stoppmarker.
    ind = [ind;
           true];
    
    % Suche nun alle eindeutigen Zeilen.
    % find(indS) = [1;3;...]
    % dual_sd(indS(1:end-1),:) = [1,2
    %                             1,3
    %                             ...
    %                             3,4]
    ind_unique = find(ind);
    nEdges = length(ind_unique)-1;
    
    % Bestimme globale KnotenIDs der dualen Knoten und sortiere diese 
    % genauso wie 'dual_sd'.
    nodesDual = find(dual);
    nodesDual = nodesDual(sortmap);
    
    % Wir koennen nun die Kantenknoten extrahieren.
    cEdges = cell(nEdges,1);
    adjSd__e = zeros(nEdges,2);  % Zur (offenen) Kante benachbarte Teilgebiete
    for i = 1:nEdges
        cEdges{i} = nodesDual(ind_unique(i):ind_unique(i+1)-1);
        adjSd__e(i,:) = dual_sd(ind_unique(i),:);
    end
    
    % Bestimme zu jeder Kante die inzidenten primalen Knoten.
    % (globale Nummerierung)
    % Algorithmus: Finde benachbarte Gebiete der primalen Knoten.
    %              Vergleiche dies mit den benachbarten Gebieten einer Kante.
    %              Sind die Gebiete der Kante eine Teilmenge der Gebiete 
    %                   des primalen Knotens, so ist der Knoten inzident 
    %                   zur Kante.
    primal_sd = findAdjacentSubdomainsOfPrimalNodes(cPrimal,l2g__sd,primal);
    primal_sd = sort(primal_sd,2);
    nodesPrimal = find(primal);
    cPrimalNodesOfEdges = cell(nEdges,1);
    for i = 1:nEdges
        cnt = 0;
        v_primal = zeros(2,1);
        for j = 1:size(primal_sd,1)
            if isSubsetOf(adjSd__e(i,:),primal_sd(j,:))
                cnt = cnt + 1;
                v_primal(cnt) = nodesPrimal(j);
            end
            if cnt == 2
                % Es kann max. 2 primale Endknoten pro Kante geben.
                break
            end
        end
        v_primal = v_primal(1:cnt);
        cPrimalNodesOfEdges{i} = v_primal;
    end
    
    % Mappe globale Nummerierung auf lokale Nummerierung.
    cPrimalNodesOfEdges_loc = cell(nEdges,2);
    for i = 1:nEdges
        for j = 1:2
            sd = adjSd__e(i,j);
            map = zeros(length(dual),1);
            map(l2g__sd{sd}) = 1:length(l2g__sd{sd});
            cPrimalNodesOfEdges_loc{i,j} = map(cPrimalNodesOfEdges{i});
        end
    end
    cPrimalNodesOfEdges = cPrimalNodesOfEdges_loc;
end

function dual_sd = findAdjacentSubdomainsOfDualNodes(cDual,l2g__sd,dual)
    % Bestimme Teilgebietszugehoerigkeit.
    n = length(dual);
    column_cnt = zeros(n,1);
    dual_sd = zeros(n,2);  % zwei Teilgebiete pro dualem Knoten
    for i = 1:length(l2g__sd)
        nodeIDs = l2g__sd{i}(cDual{i});  % globale KnotenIDs dualer Knoten
        column = column_cnt(nodeIDs)+1;
        for j = 1:length(column)
            dual_sd(nodeIDs(j),column(j)) = i;
        end
        column_cnt(nodeIDs) = column_cnt(nodeIDs) + 1;
    end
    dual_sd = dual_sd(dual,:);
end

function primal_sd = findAdjacentSubdomainsOfPrimalNodes(cPrimal,l2g__sd,primal)
    % Bestimme Teilgebietszugehoerigkeit.
    n = length(primal);
    column_cnt = zeros(n,1);
    primal_sd = zeros(n,4);  % vier Teilgebiete pro primalem Knoten
    for i = 1:length(l2g__sd)
        nodeIDs = l2g__sd{i}(cPrimal{i});  % globale KnotenIDs primaler Knoten
        column = column_cnt(nodeIDs)+1;
        for j = 1:length(column)
            primal_sd(nodeIDs(j),column(j)) = i;
        end
        column_cnt(nodeIDs) = column_cnt(nodeIDs) + 1;
    end
    primal_sd = primal_sd(primal,:);
end

function isSubset = isSubsetOf(A,B)
% A and B must be sorted ascendingly and length(A) <= length(B).
    index_B = 1;
    for i = 1:length(A)
        bEqual = false;
        while not(bEqual)
            if A(i) == B(index_B)
                bEqual = true;
            else
                index_B = index_B + 1;
                remainingSlots_B = length(B)-index_B+1;
                remainingSlots_A = length(A)-i+1;
                if remainingSlots_B < remainingSlots_A
                    break
                end
            end
        end
        if not(bEqual), isSubset = false; return; end
    end
    isSubset = true;
end