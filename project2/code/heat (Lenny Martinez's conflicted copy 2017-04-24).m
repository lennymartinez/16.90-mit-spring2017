function T = heat(Mesh)
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

  J = (x(2,1)-x(1,1))*(x(3,2)-x(1,2))-(x(3,1)-x(1,1))*(x(2,2)-x(1,2));
  dksi = (1/J).*[(x(3,2)-x(1,2)), -(x(3,1)-x(1,1)); -(x(2,2)-x(2,1)), x(2,1)-x(1,1)]; 
  dphi= [-1,-1; 1,0;0,1];
  
  gradphi = dphi*dksi;
  
  % Increment the stiffness matrix at all locations dependent on
  % current element
  material = Mesh.Elem2Material(elem,1);
  if material == 1 %steel
      k = 40;
  elseif material == 0 %magnesia paint
      k = 2.5;
  else %material is some third thing that's wrong
      break
      fprintf('Something with the material assignments is wrong.');
  end
  Ae = J/2;
  L = gradphi*transpose(gradphi);
  L = L.*(k*Ae);
  
%   K(I(1), I(1)) = K(I(1), I(1)) + L(1);
%   K(I(1), I(2)) = K(I(1), I(2)) + L(2);
%   K(I(1), I(3)) = K(I(1), I(3)) + L(3);
%   K(I(2), I(1)) = K(I(2), I(1)) + L(4);
%   K(I(2), I(2)) = K(I(2), I(2)) + L(5);
%   K(I(2), I(3)) = K(I(2), I(3)) + L(6);
%   K(I(3), I(1)) = K(I(3), I(1)) + L(7);
%   K(I(3), I(2)) = K(I(3), I(2)) + L(8);
%   K(I(3), I(3)) = K(I(3), I(3)) + L(9);
  
  K(I,I) = K(I,I) + L;
  
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
  
  ds = sqrt((x(2,1)-x(1,1))^2+(x(2,2)-x(1,2))^2);
  
  iBC = Mesh.BC(edge,3);
  if (iBC == 0)     % hot gas side
      h = 10000;
      T = 3000;
  elseif (iBC == 1) % cold gas side
      h = 20000;
      T = 300;
  elseif (iBC == 2) % insulated top wall 
    
    % FILL THIS IN
    
  elseif (iBC == 3) % symmetry
    
    % FILL THIS IN
    
  else
    fprintf(1, 'Error.\n'); % Just in case Mesh is invalid
    return;
  end
  M = [1/3, 1/6; 1/6, 1/3];
  F(I) = F(I) + 0.5*h*T*ds;
  K(I,I) = K(I,I) + h*ds*M;
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
  
  ds = sqrt((x(2,1)-x(1,1))^2+(x(2,2)-x(1,2))^2);
  
  iBC = Mesh.BC(edge,3);
  if (iBC == 0)     % hot gas side
    
    % FILL THIS IN
    
  elseif (iBC == 1) % cold gas side
    
    % FILL THIS IN

  end
end

fprintf(1, 'qhot = %g, qcold = %g\n', qhot, qcold);
