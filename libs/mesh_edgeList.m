function [edges,elements_byEdgeIDs,adjacentElements__e] = mesh_edgeList(elements)
% Given a triangulation, generate an edge list, a list of adjacent 
% triangles (to an edge), and describe the triangulation by edge IDs.
%
% INPUT:
%    elements: n x 3 matrix where each row describes a triangle by its 
%              node IDs.
%
% OUTPUT:
%    edges: Matrix where each row describes an edge by its node IDs.
%    elements_byEdgeIDs: Matrix where each row describes an element by its
%       edge IDs. The edge IDs correspond to the 'edges' matrix.
%    adjacentElements__e: cell array; each entry {i} contains a list for 
%       the edge 'i' of IDs of all adjacent triangles.
%
    %
    assert(size(elements,2) == 3,'Passed elements matrix must have three columns.')
    
    edges_notUnique = zeros(numel(elements),2);
    elementID__e = zeros(numel(elements),1);
    
    nel = size(elements,1);
    
    % 1) Extract all edges from the element matrix and generate the 
    %    element matrix based on edge IDs.
    % Note: There will be many duplicate edges (all interior edges).
    edges_elementIndices = [1,2; ...  % edge tri(i,[1,2])
                            1,3; ...  % edge tri(i,[1,3])
                            2,3];     % edge tri(i,[2,3])
    ne = size(edges_elementIndices,1);
    elements_byEdgeIDs = zeros(size(elements,1),ne);
    cnt = 0;
    for i = 1:ne
        column1 = edges_elementIndices(i,1);
        column2 = edges_elementIndices(i,2);
        
        preliminaryEdgeIDs = cnt+1:cnt+nel;
        edges_notUnique(preliminaryEdgeIDs,:) = elements(:,[column1,column2]);
        
        % To which element does the edge belong:
        elementID__e(preliminaryEdgeIDs) = 1:nel;
        
        elements_byEdgeIDs(:,i) = preliminaryEdgeIDs;
        cnt = cnt+nel;
    end
    
    % 2) Sort the edges and related instances accordingly.
    % --> Pre-processing step to allow duplicate removal in the next step.
    edges_notUnique = sort(edges_notUnique,2);
    [edges_notUnique,sort_map] = sortrows(edges_notUnique);
    elementID__e = elementID__e(sort_map);
    %
    % Set-up the inverse sorting map, to map the old edge IDs in 
    % 'elements_byEdgeIDs' to the new ones.
    sort_map_inv(sort_map) = 1:length(sort_map);
    elements_byEdgeIDs = sort_map_inv(elements_byEdgeIDs);
    
    % Find and extract unique edges.
    ind_unique_edges = [true; ...
                        any( diff(edges_notUnique) ,2)];
    edges = edges_notUnique(ind_unique_edges,:);
    
    % Generate adjacency list.
    adjacentElements__e = getAdjacentElements(elementID__e,ind_unique_edges);
    
    % Map preliminary (and sorted) edge IDs to unique IDs.
    map_to_unique_ID = cumsum(ind_unique_edges)';
    elements_byEdgeIDs = map_to_unique_ID(elements_byEdgeIDs);
    
end

function adjacentElements__e = getAdjacentElements(elementID__e,ind_unique_edges)
    nUniqueEdges = nnz(ind_unique_edges);
    adjacentElements__e = cell(nUniqueEdges,1);
    
    % Determine the number of elements adjacent to the edge 'i'.
    nAdjacentElements = diff([ find(ind_unique_edges); ...
                               length(ind_unique_edges)+1 ]);
    
    % Generate adjacency list.
    cnt = 0;
    for i = 1:length(nAdjacentElements)
        adjacentElements__e{i} = elementID__e(cnt+1:cnt+nAdjacentElements(i));
        cnt = cnt+nAdjacentElements(i);
    end
end