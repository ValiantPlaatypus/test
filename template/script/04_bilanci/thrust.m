%
%
%

clear all ; clc ; close all

%> Flight velocity and ambient TD state
u =   0.0 ;  Pa = 101325 ;  Ta = 288.15 ;
u = 250.0 ;  Pa =  18750 ;  Ta = 216.70 ;
u = 250.0 ;  Pa =  18750 ;  Ta = 216.70 ;

% p ratio vector
pratio = [ 1.1:0.2:1.8, 2.0:1.0:100.0 ] ;
nratio = length(pratio) ;
Tm_v   = zeros(size(pratio));
tsfc_v = zeros(size(pratio));
temp04 = [ 1500, 1600, 1700 ] ;
ntemp04 = length(temp04) ;

Pcomp_m = zeros(length(pratio), length(temp04)) ;
hexp_m  = zeros(length(pratio), length(temp04)) ;
ue_m    = zeros(length(pratio), length(temp04)) ;

figure ; hold on

for j = 1 : ntemp04
for i = 1 : nratio

  %> Compressor p ration and turbine inlet T
  % P03_P02 = 20.0 ;   % Compressor pressure ratio
  % T04     = 1500.0 ; % Inlet turbine temperature [K]
  P03_P02 = pratio(i) ;   % Compressor pressure ratio
  T04     = temp04(j) ;
  P04_P03 =  1.0 ;   % Burner pressure ratio
  
  %> Efficiencies and specific heat ratio
  eta_d = 0.97 ;  gam_d = 1.40 ;
  eta_c = 0.85 ;  gam_c = 1.37 ;
  eta_b = 1.00 ;  gam_b = 1.35 ;
  eta_t = 0.90 ;  gam_t = 1.33 ;
  eta_n = 0.98 ;  gam_n = 1.36 ;
  
  %> Fuel heating value
  qr = 45.0e+6 ; % [J/K]
  
  %> a. Free stream conditions 
  R = 287.0   ; % J/kg K
  gamma = 1.4 ;
  cv = R / ( gamma - 1.0 ) ;
  cp = gamma * cv ;
  ca = sqrt( gamma * R * Ta );
  M = u / ca ;
  
  %> a. -> 2. Diffuser
  %> 2. Total TD state at the compressor inlet
  T02 = Ta * ( 1.0 + 0.5*( gamma-1.0 )*M^2 ) ; 
  P02 = Pa * ( 1.0 + eta_d * ( T02/Ta - 1.0 ) ) ^ ( gam_d / ( gam_d - 1.0 ) ) ; 
  
  % 2. -> 3. Compressor
  %> 3.
  P03 = P03_P02 * P02 ;
  T03 = T02 * ( 1.0 + ( P03_P02^((gam_c-1)/gam_c) - 1.0 ) / eta_c ) ;
  
  % 3. -> 4. Combustor
  f = ( T04/T03 - 1.0 ) / ( qr/(cp*T03) - T04/T03 ) ;
  P04 = P04_P03 * P03 ;
  
  % 4. -> 5. Turbine
  T05 = T04 - ( T03 - T02 ) ;
  P05 = P04 * ( 1.0 + ( T05/T04 - 1.0 )/eta_t ) ^ ( gam_t / ( gam_t - 1.0 ) ) ;
  
  % 5. -> 6. Nozzle
  T06 = T05 ;
  ue = sqrt( 2.0 * eta_n * gam_n / ( gam_n - 1.0 ) * R * T06 * ...
             ( 1.0 - ( Pa / P05 )^( (gam_n-1)/gam_n ) ) ) ;
  
  Tm   = ( 1.0 + f ) * ue - u ;
  tsfc = f / Tm ;

  % Save results
  Tm_v(  i) = Tm ;
  tsfc_v(i) = tsfc ;
  
  Pcomp = cp * ( T03 - T02 ) ;
  Pcomp_m(i,j) = Pcomp ; 
  hexp_m( i,j) = cp * T06 ;
  ue_m(   i,j) = ue ;

end

semilogx(pratio,Tm_v/1000 ,'b-o','LineWidth',2)
semilogx(pratio,tsfc_v*10e4,'r-o','LineWidth',2)

end

hold off, grid on


figure('Position',[ 50 50 1200 500])
subplot(1,3,1), plot(pratio, Pcomp_m ), grid on, title('Compressor power') 
subplot(1,3,2), plot(pratio,  hexp_m ), grid on, title('D enthalpy for exp.')
subplot(1,3,3), plot(pratio,    ue_m ), grid on, title('u_e')
