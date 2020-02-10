clear all; close all;

load Meshes
mesh = NozzleMesh3;
vary = linspace(0.75,1.25,5);

%% Question 1

var1 = []; %magnesia
var2 = []; %hhot
var3 = []; %hcold


for i = 1:length(vary)
    k = vary(i);
    T1 = heat(mesh,k,1,1);
    T2 = heat(mesh,1,k,1);
    T3 = heat(mesh,1,1,k);
    
    var1 = [var1, max(T1(117:833))];
    var2 = [var2, max(T2(117:833))];
    var3 = [var3, max(T3(117:833))];
    
end
dx = size(vary);
m1 = var1(end)-var1(1); m2 = var2(end)-var2(1); m3 = var3(end)-var3(1); 
m = [m1,m2,m3]./dx(2);
x1 = vary*2.5;
x2 = vary*10000;
x3 = vary*20000;


figure(1); plot(x1,var1,'*','LineWidth',3); grid on; 
xlabel('k_{magnesia}');ylabel('T_{max} in steel')
figure(2); plot(x2,var2,'*','LineWidth',3); grid on;
xlabel('h_{hot}');ylabel('T_{max} in steel')

figure(3); plot(x3,var3,'*','LineWidth',3); grid on;
xlabel('h_{cold}');ylabel('T_{max} in steel')



%% Question 2

% z = [vary(1), vary(3), vary(end)];
% ind = [];
% 
% for j = 1:3,
%     T1 = heat(mesh,z(j),1,1);
%     T2 = heat(mesh,1,z(j),1);
%     T3 = heat(mesh,1,1,z(j));
%     
%     
%     % Create temeperature plots.
% %     plotsolution(mesh,T1)
% %     plotsolution(mesh,T2)
% %     plotsolution(mesh,T3)
% % Figure out where location of max Temperature in steel is.
%     Tmax1 = max(T1(117:833));
%     Tmax2 = max(T2(117:833));
%     Tmax3 = max(T3(117:833));
%     index1 = find(Tmax1 == T1);
%     index2 = find(Tmax2 == T2);
%     index3 = find(Tmax3 == T3);
%     ind = [ind; index1, index2, index3];
% end


%% Bonus Q

% A = [];
% n = length(vary);
% for a = 1:n
%     for b = 1:n
%         T = heat(mesh,1,a,b);
%         A(a,b) = max(T(117:833));
%     end
% end
% 
% contourf(A);