%
%
%

clear all ; clc ; close all

%> Gas constants 
R = 287.0   ; % J/kg K
gamma = 1.4 ;

%> Flight velocity and ambient TD state
u = 250.0 ;  Pa =  26500 ;  Ta = 223.25 ;
rhoa = Pa / ( R * Ta ) ;
cv = R / ( gamma - 1.0 ) ;
cp = gamma * cv ;
ca = sqrt( gamma * R * Ta );
M = u / ca ;
hta = u^2/2 + cp * Ta ;

A2 = 1.0 ; % m^2

%> Fuel heating value
qr = 46.0e+6 ; % [J/K]

% p ratio vector
pratio = [ 40.0 ] ;
nratio = length(pratio) ;
Tm_v   = zeros(size(pratio));
tsfc_v = zeros(size(pratio));
temp04 = [ 1600 ] ;
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
  eta_d = 1.00 ;  gam_d = 1.40 ;
  eta_c = 1.00 ;  gam_c = 1.40 ;
  eta_b = 1.00 ;  gam_b = 1.40 ;
  eta_t = 1.00 ;  gam_t = 1.40 ;
  eta_n = 1.00 ;  gam_n = 1.40 ;

  %> a. -> 2. Diffuser
  P2_Pa = 1.5 ;
  P2 = P2_Pa * Pa ;
  rho2 = ( P2_Pa )^(1.0/gam_d) * rhoa ;
  T2 = P2 / ( R * rho2 ) ;
  ht2 = hta ;
  h2 = cp * T2 ;
  V2 = sqrt( 2*( ht2 - h2 ) ) ;
  mdot2 = rho2 * A2 * V2 ;

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
