% 
% Compare results of hessSmith.m, xfoil and openfoam

clear all; clc; close all

folder_hs = './output/' ;
folder_xf = './xfoil/' ;
al = 0;
re = 1e+5;

filen_cp = [ 're+1e', num2str( log10(re), '%1.1d' ) , '_' , ...
             'al+'  , num2str(        al, '%1.1d' ) , '_cp.dat' ] ;
filen_bl = [ 're+1e', num2str( log10(re), '%1.1d' ) , '_' , ...
             'al+'  , num2str(        al, '%1.1d' ) , '_bl.dat' ] ;

% % Re = 1e+6, al = 0°
% filen_of1 = [ './output/', 'komsstlm_naca0012_re+1e6_al+0.dat' ] ;
% filen_of2 = [ './output/',          'naca0012_re+1e6_al+0.dat' ] ;
% of1 = dlmread( filen_of1, ',', 1, 0 );
% of2 = dlmread( filen_of2, ',', 1, 0 );
% % modification to the LM cp solution
% of1(:,4) = of1(:,4)-0.025;
% 
% of1_tauw = sqrt( sum( of1(:, 7: 9).^2, 2 ) );
% of2_tauw = sqrt( sum( of2(:, 4: 6).^2, 2 ) );
% of1_cf   = 2.0 .* of1_tauw ;
% of2_cf   = 2.0 .* of2_tauw ;

%> Pressure coefficient ===
cp_hs = dlmread( [ folder_hs, filen_cp] , '', 0, 0); % x, cp
cp_xf = dlmread( [ folder_xf, filen_cp] , '', 1, 0); % x, cp
%> Boundary layer ===
% hessSmith.m: x, cf, tauw, delta_star, theta, H
bl_hs = dlmread( [ folder_hs, filen_bl] , '', 0, 0);
% xFoil: s, x, y, Ue/Uinf, delta_star, theta, cf, H, H*, p, m 
bl_xf = dlmread( [ folder_xf, filen_bl] , '', 1, 0);

% === pressure coefficient ===
figure; hold on
plot(cp_hs(:, 1),cp_hs(:,2),'LineWidth',2)
plot(cp_xf(:, 1),cp_xf(:,2),'LineWidth',2)
% plot(  of1(:,14),2*of1(:,4),'LineWidth',2)
% plot(  of2(:,11),2*of2(:,3),'LineWidth',2)
legend('lab','xfoil','OF-LM','OF'), xlabel('x/c'), ylabel('c_P'), title('Re=1e+6, cP')
grid on, hold off

% === friction coefficient ===
figure; hold on
plot(bl_hs(:, 1),2.*bl_hs(:,3),'LineWidth',2)
plot(bl_xf(:, 2),   bl_xf(:,7),'LineWidth',2)
% plot(  of1(1:150,14),     of1_cf(1:150)  ,'LineWidth',2)
% plot(  of2(:,11),     of2_cf  ,'LineWidth',2)
legend('lab','xfoil','OF-LM','OF'), xlabel('x/c'), ylabel('c_F')
xlim([-.1 1.1]), title('Re=1e+6, cF')
grid on, hold off

% === delta^* ===
figure; hold on
plot(bl_hs(:,1),bl_hs(:,4),'LineWidth',2)
plot(bl_xf(:,2),bl_xf(:,5),'LineWidth',2)
legend('lab','xfoil'), xlabel('x/c'), ylabel('\delta^*')
xlim([-.1 1.1]), title('Re=1e+5, \delta^*')
grid on, hold off

% === theta ===
figure; hold on
plot(bl_hs(:,1),bl_hs(:,5),'LineWidth',2)
plot(bl_xf(:,2),bl_xf(:,6),'LineWidth',2)
legend('lab','xfoil'), xlabel('x/c'), ylabel('\theta')
xlim([-.1 1.1]), title('Re=1e+6, \theta')
grid on, hold off

% === H ===
figure; hold on
plot(bl_hs(:,1),bl_hs(:,6),'LineWidth',2)
plot(bl_xf(:,2),bl_xf(:,8),'LineWidth',2)
legend('lab','xfoil'), xlabel('x/c'), ylabel('H')
xlim([-.1 1.1]), title('Re=1e+6, H')
grid on, hold off
