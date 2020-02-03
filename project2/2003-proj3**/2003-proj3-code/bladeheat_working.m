% Clear variables

clear all;


% Load in the grid file

fname = input('Enter gridfile name: ','s');
load(fname);

% NOTE:  after loading a gridfile using the load(fname) command,
%        three important grid variables and data arrays exist.  These are:
%
% Nt: Number of triangles (i.e. elements) in mesh
%
% Nv: Number of nodes (i.e. vertices) in mesh
%
% Nbc: Number of edges which lie on a boundary of the computational
%      domain.
% 
% tri2nod(3,Nt):  list of the 3 node numbers which form the current
%                 triangle.  Thus, tri2nod(1,i) is the 1st node of
%                 the i'th triangle, tri2nod(2,i) is the 2nd node
%                 of the i'th triangle, etc.
%
% xy(2,Nv): list of the x and y locations of each node.  Thus,
%           xy(1,i) is the x-location of the i'th node, xy(2,i)
%           is the y-location of the i'th node, etc.
%
% bedge(3,Nbc): For each boundary edge, bedge(1,i) and bedge(2,i) 
%               are the node numbers for the nodes at the end
%               points of the i'th boundary edge.  bedge(3,i) is an
%               integer which identifies which boundary the edge is
%               on. In this solver, the third value has the
%               following meaning:
%
%               bedge(3,i) = 0: edge is on the airfoil
%               bedge(3,i) = 1: edge is on the first cooling passage
%               bedge(3,i) = 2: edge is on the second cooling passage
%               bedge(3,i) = 3: edge is on the third cooling passage
%               bedge(3,i) = 4: edge is on the fourth cooling passage
% 

% Zero stiffness matrix

K = zeros(Nv, Nv);
b = zeros(Nv, 1);


% Loop over elements and calculate residual and stiffness matrix

for ii = 1:Nt,
  
  % Get node numbers for current element
  kn(1) = tri2nod(1,ii);
  kn(2) = tri2nod(2,ii);
  kn(3) = tri2nod(3,ii);
    
  % Get x location of nodes
  xe(1) = xy(1,kn(1));
  xe(2) = xy(1,kn(2));
  xe(3) = xy(1,kn(3));

  % Get y location of nodes
  ye(1) = xy(2,kn(1));
  ye(2) = xy(2,kn(2));
  ye(3) = xy(2,kn(3));

  
  % Calculate all of the necessary shape function derivatives, the
  % Jacobian of the element, etc.

  % FILL THIS IN
  
  
  % Increment the stiffness matrix at all locations dependent on
  % current element
  %
  % This will include a bunch of lines of the following form
  
  %K(kn(1), kn(1)) = K(kn(1), kn(1)) + FILL THIS IN
  %K(kn(1), kn(2)) = K(kn(1), kn(2)) + FILL THIS IN
  %K(kn(1), kn(3)) = K(kn(1), kn(3)) + FILL THIS IN
  %
  %K(kn(2), kn(1)) = K(kn(2), kn(1)) + FILL THIS IN
  %     etc
  
end


% Loop over boundary edges and account for bc's

for ii = 1:Nbc,
  
  kn(1) = bedge(1,ii);
  kn(2) = bedge(2,ii);
  
  % FILL IN CODE HERE TO SET CONVECTIVE HEAT TRANSFER BC 
  % NOTE: THE EXACT BOUNDARY IS STORED IN bedge(3,ii)

end

% Solve for temperature

Tsol = K\b;


% Plot solution

bladeplot;