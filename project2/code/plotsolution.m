function plotsolution(Mesh, T)
% function plotsolution(Mesh, T);
% Plots a mesh and solution

% change these as necessary
plotsol  = 1;
plotmesh = 0;

figure; hold on;

if (plotsol == 1)
  h = patch('Vertices',1000*Mesh.Coord,'Faces',Mesh.Elem2Node,'FaceVertexCData',T, ...
      'FaceColor','interp','EdgeColor','none');
  axis('equal');
  axis([-.25, 2, -.1, 2.1]);
  axis off;
  caxis([500, 1300]);
  title('Temperature');
  HC = colorbar;
end

if (plotmesh == 1)
  for elem = 1:Mesh.nElem
    xy = Mesh.Coord(Mesh.Elem2Node(elem,:),:);
    xy = [xy; xy(1,:)];
    plot(1000*xy(:,1), 1000*xy(:,2), 'k-'); % mm
  end
  axis('equal');
  axis([-.25, 2, -.1, 2.1]);
  axis off;
  return;
end
