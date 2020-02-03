% This loops of the boundary edges and plots them

hold on;
for ii = 1:Nbc,
  kn1 = bedge(1,ii);
  kn2 = bedge(2,ii);
  plot([xy(1,kn1), xy(1,kn2)],[xy(2,kn1), xy(2,kn2)]);
end
hold off;