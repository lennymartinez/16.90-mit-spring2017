% Plot T in triangles
figure;
for ii = 1:Nt,
  for nn = 1:3,
    xtri(nn,ii) = xy(1,tri2nod(nn,ii));
    ytri(nn,ii) = xy(2,tri2nod(nn,ii));
    Ttri(nn,ii) = Tsol(tri2nod(nn,ii));
  end
end
HT = patch(xtri,ytri,Ttri);
axis('equal');
set(HT,'LineStyle','none');
title('Temperature');
caxis([900,1300]);
HC = colorbar;
hold on; bladeplotgrid; hold off;