km = 2.5; 
ks = 40;
L = 0.0005;
hh = 10000;
hc = 20000;
TL = 3000;
TR = 300;


syms T1 T2 T0;

T0 = (3000*hh + (km/L)*T1)/((km/L)+hh);
T2 = (300*hc + (ks/L)*T1)/((ks/L)+hc);
T1 = (km*T0+ ks*T2)/(km+ks);


T1 = double(22200/29);
T0 = (3000*hh + (km/L)*T1)/((km/L)+hh);
T2 = (300*hc + (ks/L)*T1)/((ks/L)+hc);

T = [T0 T1 T2]


q = [hh*(3000-T0); (km/L)*(T0-T1); (ks/L)*(T1-T2); hc*(-300+T2) ];