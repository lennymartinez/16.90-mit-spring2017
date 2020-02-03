Mesh = TestMesh;

% function T = heat(Mesh)
%--------------------------------------------------------------------%
% function T = heat(Mesh);
% 
% This function calculates nodal temperature values for an input Mesh
% using the finite element method.
%
% INPUT    Mesh: structure containing nodes, elements, etc.
% OUTPUT      T: vector of temperature values at Mesh nodes.
%
%
% The input Mesh is a structure with the following fields
%
% nElem:  Number of elements (i.e. triangles) in the mesh
%
% nNode:  Number of nodes (i.e. vertices) in the mesh
%
% Elem2Node(nElem, 3): for each element, a list of the 3 node numbers
%                      which form the element.  Thus, Elem2Node(i,1)
%                      is the 1st node of the i'th element,
%                      Elem2Node(i,2) is the 2nd node of the i'th
%                      element, etc.
%
% Elem2Material(nElem, 1): for each element a flag indicating which 
%                          material the element is in.  0 is the
%                          magnesia paint, 1 is steel.
%
% Coord(nNode,2): list of the x and y locations of each node.  Thus,
%                 Coord(1,i) is the x-location of the i'th node and
%                 Coord(2,i) is the y-location of the i'th node.
%
% nEdge: Number of edges which lie on a boundary of the computational
%        domain.
%
% BC(nEdge,3):  For each boundary edge, i, BC(i,1) and BC(i,2) 
%               are the node numbers for the nodes at the endpoints
%               points of the edge.  BC(i,3) is an
%               integer which identifies which boundary the edge is
%               on.  This third value can be:
%
%               BC(i,3) = 0: edge is on the hot-gas side
%               BC(i,3) = 1: edge is on the cold-gas side
%               BC(i,3) = 2: edge is on the symmetry plane
%               BC(i,3) = 3: edge is on the top insulated wall
% 
%--------------------------------------------------------------------%



% Initialize stiffness matrix and right-hand side
K = sparse(Mesh.nNode, Mesh.nNode); % global stiffness matrix
F = zeros(Mesh.nNode, 1);           % right-hand side


%-----------------------------------------------%
% Construct Stiffness (K) and Mass (M) matrices %
%-----------------------------------------------%

for elem = 1:Mesh.nElem % loop over elements in the mesh
    
  % I(1), I(2), I(3) are the node numbers of the current element
  I = Mesh.Elem2Node(elem,:);
  
  % x(i,:) is the (x,y) coordinate of node i, 1 <= i <= 3
  x = Mesh.Coord(I,:);
  
  
  % Calculate (using the reference element) all of the necessary
  % shape function derivatives, the Jacobian for the element,
  % etc.
  
  % unpack matrix 
  x1 = x(1,1);
  x2 = x(2,1);
  x3 = x(3,1);
  y1 = x(1,2);
  y2 = x(2,2);
  y3 = x(3,2);
  
  % got this from the notes 
  dp1dz1 = -1;
  dp2dz1 = 1;
  dp3dz1 = 0;
  
  dp1dz2 = -1;
  dp2dz2 = 0;
  dp3dz2 = 1;
  
  % Jacobian
  J = (x2-x1)*(y3-y1)-(x3-x1)*(y2-y1);
  
  dzdx = (1/J)*[(y3-y1),-(x3-x1);-(y2-y1),(x2-x1)];
  
  dp1dx = dp1dz1*dzdx(1,1)+dp1dz2*dzdx(2,1);
  dp2dx = dp2dz1*dzdx(1,1)+dp2dz2*dzdx(2,1);
  dp3dx = dp3dz1*dzdx(1,1)+dp3dz2*dzdx(2,1);
  
  dp1dy = dp1dz1*dzdx(1,2)+dp1dz2*dzdx(2,2);
  dp2dy = dp2dz1*dzdx(1,2)+dp2dz2*dzdx(2,2);
  dp3dy = dp3dz1*dzdx(1,2)+dp3dz2*dzdx(2,2);
  
  
  % Increment the stiffness matrix at all locations dependent on
  % current element
  %
  % This will include a bunch of lines of the following form
  %
  K(I(1), I(1)) = K(I(1), I(1)) + dot([dp1dx,dp1dy],[dp1dx,dp1dy]);
  K(I(1), I(2)) = K(I(1), I(2)) + dot([dp1dx,dp1dy],[dp2dx,dp2dy]);
  K(I(1), I(3)) = K(I(1), I(3)) + dot([dp1dx,dp1dy],[dp3dx,dp3dy]);
  
  K(I(2), I(1)) = K(I(2), I(1)) + dot([dp2dx,dp2dy],[dp1dx,dp1dy]);
  K(I(2), I(2)) = K(I(2), I(2)) + dot([dp2dx,dp2dy],[dp2dx,dp2dy]);
  K(I(2), I(3)) = K(I(2), I(3)) + dot([dp2dx,dp2dy],[dp3dx,dp3dy]);
  
  K(I(3), I(1)) = K(I(3), I(1)) + dot([dp3dx,dp3dy],[dp1dx,dp1dy]);
  K(I(3), I(2)) = K(I(3), I(2)) + dot([dp3dx,dp3dy],[dp2dx,dp2dy]);
  K(I(3), I(3)) = K(I(3), I(3)) + dot([dp3dx,dp3dy],[dp3dx,dp3dy]);
  
  %
  %  or, if you're comfortable working with matrices, one line:
  %
  %  K(I,I) = K(I,I) + FILL THIS IN
  %
  %  or, if you're attempting the bonus question, you can construct
  %  the list of indices and values used to construct a sparse matrix
  %  using sparse(i,j,v); 
  
end


%----------------------------------------------------------%
% Add boundary contribution to K and build right-hand side %
%----------------------------------------------------------%
  
for edge = 1:Mesh.nEdge % loop over boundary edges
    
  % I(1) and I(2) are the nodes on the boundary edge
  I = Mesh.BC(edge,1:2);
  
  % (x,y) coordinates of edge nodes
  x = Mesh.Coord(I,:);
  
  % Calculate edge length and any required shape function integrals
  
  % FILL THIS IN
  
  iBC = Mesh.BC(edge,3);
  if (iBC == 0)     % hot gas side
    
     F(I) = 3000;
     K(I,I) = 1;
    
  elseif (iBC == 1) % cold gas side
    
     F(I) = 300;
     K(I,I) = 1;
    
  elseif (iBC == 2) % insulated top wall 
      
     % need to fill in
    
  elseif (iBC == 3) % symmetry
    
    % need to fill in
    
  else
    fprintf(1, 'Error.\n'); % Just in case Mesh is invalid
    return;
  end
  
end


%--------------%
% Solve system %
%--------------%

% Estimate the condition number of K.  If this is "inf" or a very very
% large number (> 10^10), K is likely singular (check your code).
fprintf(1, 'Condition number of K = %e\n', condest(K));

T = K\F; % solve system


%--------------------------------%
% Calculate integrated heat flux %
%--------------------------------%

qhot  = 0;  % This will be the heat flux on the hot gas side
qcold = 0;  % This will be the heat flux on the cold gas side

for edge = 1:Mesh.nEdge % loop over boundary edges
  
  % I(1) and I(2) are the nodes on the boundary edge
  I = Mesh.BC(edge,1:2);
  
  % (x,y) coordinates of edge nodes
  x = Mesh.Coord(I,:);
  
  % Calculate edge length and any other required values
  
  % FILL THIS IN
  
  iBC = Mesh.BC(edge,3);
  if (iBC == 0)     % hot gas side
    
    % FILL THIS IN
    
  elseif (iBC == 1) % cold gas side
    
    % FILL THIS IN

  end
end

fprintf(1, 'qhot = %g, qcold = %g\n', qhot, qcold)
