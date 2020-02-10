close all;

%% Error for T(1)
T1=heat(NozzleMesh1,1,1,1);
T2=heat(NozzleMesh2,1,1,1);
T3=heat(NozzleMesh3,1,1,1);
T4=heat(NozzleMesh4,1,1,1);

T = [T1(1),T2(1),T3(1)];
eT = abs(T-T4(1));
m = log(eT(end)/eT(1))/log(1504/94)
loglog([94,376,1504],eT,'r*'); grid on;
ylabel('Error');
xlabel('Number of elements in Mesh')


%% Error for q
qtruth = 28258.9;
q1 = 28341.6;
q2 = 28289.4;
q3 = 28267.7;

e = [q1,q2,q3] - qtruth;

loglog([94,376,1504],e,'r*'); grid on;
ylabel('Error');
xlabel('Number of elements in Mesh')
m2 = log(e(end)/e(1))/log(1504/94)