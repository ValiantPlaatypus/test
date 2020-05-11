%
% this part of the code is meant only for one airfoil
% to be cleaned...
%

% === Find stagnantion point ===
n_stg =  0  ;
i_stg = [ ] ;
for i = 1 : nelems-1
  if ( vTi(i) * vTi(i+1) <= 0 ) % stg point found
     n_stg = n_stg + 1 ;
     i_stg = [ i_stg , i ] ;

     % Find the coordinate along the airfoil of the stg point
     csi_stg = - vTi(i) / ( vTi(i+1) - vTi(i) ) ;
  end
end

if ( n_stg ~= 1 )
    fprintf(' More than one stg point found on the airfoil. \n')
    fprintf(' Something strange happend. Stop.              \n')
    stop
end

% arc-length from the stagnation point 
s = zeros(nelems,1) ;
s(i_stg+1) =         csi_stg   * ( len(i_stg) + len(i_stg+1) ) * 0.5 ;
s(i_stg  ) = ( 1.0 - csi_stg ) * ( len(i_stg) + len(i_stg+1) ) * 0.5 ;
for i = i_stg + 2 : nelems
   s(i) = s(i-1) + 0.5 * ( len(i) + len(i-1) ) ;
end
for i = i_stg -1 : -1 : 1
   s(i) = s(i+1) + 0.5 * ( len(i) + len(i+1) ) ;
end


% === Initial condition of the momentum thickness ===
th = acos( nvers(:,i_stg)' * nvers(:,i_stg+1) ) ;
r0 = len(i+1)*sin(th) + (len(i)+len(i+1)*cos(th))/tan(th) ;
k = freeStream.v / r0 ;
theta0 = sqrt( 0.075 * freeStream.kin_visc / k ) ;
% theta0 = 0.0 ;


% === Integrate Thwaites' equation in both the directions ===
Ue        = abs(vTi) ;
theta2Ue6 = zeros(nelems,1) ;
theta     = zeros(nelems,1) ;
ReTheta   = zeros(nelems,1) ;
lambda    = zeros(nelems,1) ;
ell       = zeros(nelems,1) ;
H         = zeros(nelems,1) ;
delta     = zeros(nelems,1) ;
cf        = zeros(nelems,1) ;
UeThetaH1 = zeros(nelems,1) ; % for turbulent b.l.
H1        = zeros(nelems,1) ; % for turbulent b.l.
% pressure and its derivative along the surface
Res    = Ue' .* s / freeStream.kin_visc ;
P      = cP .* ( 0.5 * freeStream.rho * freeStream.v^2.0 ) ;
dPdx   = zeros(nelems,1) ;
dUedx  = zeros(nelems,1) ;

% Pressure derivative ---
Ptot = freeStream.P + 0.5 * freeStream.rho * freeStream.v^2.0 ; % Ptot at stg pt 
dPdx(i_stg+1) = ( P(i_stg+1) - Ptot ) / ...  
                ( 0.5 * csi_stg * ( len(i_stg) + len(i_stg+1) ) ) ;
dPdx(i_stg  ) = ( P(i_stg  ) - Ptot ) / ...  
                ( 0.5 * (1.0-csi_stg) * ( len(i_stg) + len(i_stg+1) ) ) ;
dUedx(i_stg+1) = ( Ue(i_stg+1) - 0.0 ) / ...  
                 ( 0.5 * csi_stg * ( len(i_stg) + len(i_stg+1) ) ) ;
dUedx(i_stg  ) = ( Ue(i_stg  ) - 0.0 ) / ...  
                 ( 0.5 * (1.0-csi_stg) * ( len(i_stg) + len(i_stg+1) ) ) ;

% initial conditions ---
theta2Ue6(i_stg+1) = 0.45 * freeStream.kin_visc * ...
                     Ue(i_stg+1)^5.0 * ...
                   ( 0.25 * csi_stg * ( len(i_stg) + len(i_stg+1) ) ) ;
theta2Ue6(i_stg  ) = 0.45 * freeStream.kin_visc * ...
                     Ue(i_stg  )^5.0 * ...
                   ( 0.25 * (1.0-csi_stg) * ( len(i_stg) + len(i_stg+1) ) ) ;
% theta2Ue6(i_stg+1) = 0.0 ;
% theta2Ue6(i_stg  ) = 0.0 ;
theta(i_stg+1)  = ( theta2Ue6(i_stg+1) / Ue(i_stg+1)^6.0 )^0.5 ;
theta(i_stg  )  = ( theta2Ue6(i_stg  ) / Ue(i_stg  )^6.0 )^0.5 ;
ReTheta(i_stg+1) = Ue(i_stg+1) * theta(i_stg+1) / freeStream.kin_visc ;
ReTheta(i_stg  ) = Ue(i_stg  ) * theta(i_stg  ) / freeStream.kin_visc ;
lambda(i_stg+1) = theta(i_stg+1)^2.0 * dUedx(i_stg+1) / freeStream.kin_visc ; 
lambda(i_stg  ) = theta(i_stg  )^2.0 * dUedx(i_stg  ) / freeStream.kin_visc ; 
ell(i_stg+1) = thwaites_ell( lambda(i_stg+1) ) ;
ell(i_stg  ) = thwaites_ell( lambda(i_stg  ) ) ;
H(  i_stg+1) = thwaites_H(   lambda(i_stg+1) ) ;
H(  i_stg  ) = thwaites_H(   lambda(i_stg  ) ) ;
delta(i_stg+1) = theta(i_stg+1) * H(i_stg+1) ;
delta(i_stg  ) = theta(i_stg  ) * H(i_stg  ) ;
cf( i_stg+1) = 2.0 * ell(i_stg+1) / ReTheta(i_stg+1) ;
cf( i_stg  ) = 2.0 * ell(i_stg  ) / ReTheta(i_stg  ) ;

% Integrate equation
regime = 'laminar' ;
transition = 'no'  ; 
for i = i_stg+2 : nelems

   if ( i < nelems ) % centered diff.
      dPdx( i) = (  P(i+1) -  P(i-1) ) / ( 0.5 * ( len(i+1) + len(i-1) ) + len(i) ) ;
      dUedx(i) = ( Ue(i+1) - Ue(i-1) ) / ( 0.5 * ( len(i+1) + len(i-1) ) + len(i) ) ;
   else
      dPdx( i) = (  P(i) -  P(i-1) ) / ( 0.5 * ( len(i) + len(i-1) ) ) ;
      dUedx(i) = ( Ue(i) - Ue(i-1) ) / ( 0.5 * ( len(i) + len(i-1) ) ) ;
   end

   if ( strcmp(regime,'laminar') )
 
      theta2Ue6(i) = theta2Ue6(i-1) + ...
                     0.45 * freeStream.kin_visc * ...
                   ( Ue(i)^5.0 + Ue(i-1)^5.0 ) *  ...
                   ( 0.25 * ( len(i) + len(i-1) ) ) ;
      
      theta(i)  = ( theta2Ue6(i) / Ue(i)^6.0 )^0.5 ;
      ReTheta(i) = Ue(i) .* theta(i) / freeStream.kin_visc ;
      lambda(i) = theta(i)^2.0 * dUedx(i) / freeStream.kin_visc ; 
      ell(i) = thwaites_ell( lambda(i) ) ;
      H(  i) = thwaites_H(   lambda(i) ) ;
      delta(i) = theta(i) * H(i) ;
      cf( i) = 2.0 * ell(i) / ReTheta(i) ;

      %%%% % Granville's method (WRONG) 
      %%%% ReThetaTr = ( 54.2124/(H(i)*(H(i)-2.48)) + 31.6 / H(i) ) * ( H(i) >= 2.591 ) + ...
      %%%%             ( 520.0/H(i) + 2.5e+6/H(i) * ( 1/H(i)-1/2.591 )^1.95 ) * ( H(i) < 2.591 ) ;
      % Michel's criterion for transition ---
      % ReThetaTr = 1.174 * ( 1.0 + 22400.0 / Res(i) ) * Res(i) ^ 0.46 ;
      ReThetaTr = 1.35 * ( 1.0 + 22400.0 / Res(i) ) * Res(i) ^ 0.46 ;

      if ( ReTheta(i) > ReThetaTr )
          regime     = 'turbulent' ;
          transition = 'transition' ;
          H(i) = 1.35 ;
          delta(i) = theta(i) * H(i) ;
          H1(i) = head_HtoH1(H(i)) ;
          UeThetaH1(i) = Ue(i) * theta(i) * H1(i) ;
          i_trans = i ; 
      end 
      

   elseif( strcmp(regime,'turbulent') )

      dx = ( 0.5 * ( len(i) + len(i-1) ) ) ;
      UeThetaH1(i) = UeThetaH1(i-1) + dx * ...
          Ue(i-1) * 0.0306 / ( H1(i-1) - 3.0 )^0.6169 ;
      theta(i) = theta(i-1) + dx * ...
          ( 0.5 * cf(i-1) - dUedx(i-1)/Ue(i-1) * ( 2.0 + H(i-1) ) * theta(i-1) ) ;

%     lambda(i) = theta(i)^2.0 * dUedx(i) / freeStream.kin_visc ; 

      H1(i) = max( UeThetaH1(i) / ( Ue(i) * theta(i) ) , 3.0 + 1e-3 ) ;
      H(i) = head_H1toH(H1(i)) ;
      delta(i) = H(i) * theta(i) ;

      ReTheta(i) = Ue(i) .* theta(i) / freeStream.kin_visc ;
      cf(i) = ludweig_tillman_cf(H(i),ReTheta(i)) ;
      
      transition = 'occurred' ;

   end

end

regime = 'laminar' ;
transition = 'no' ;
for i = i_stg-1 : -1 : 1

   if ( i > 1 ) % centered diff.
      dPdx( i) = (  P(i-1) -  P(i+1) ) / ( 0.5 * ( len(i+1) + len(i-1) ) + len(i) ) ;
      dUedx(i) = ( Ue(i-1) - Ue(i+1) ) / ( 0.5 * ( len(i+1) + len(i-1) ) + len(i) ) ;
   else
      dPdx( i) = (  P(i) -  P(i+1) ) / ( 0.5 * ( len(i) + len(i+1) ) ) ;
      dUedx(i) = ( Ue(i) - Ue(i+1) ) / ( 0.5 * ( len(i) + len(i+1) ) ) ;
   end

   if ( strcmp(regime,'laminar') )
 
      theta2Ue6(i) = theta2Ue6(i+1) + ...
                     0.45 * freeStream.kin_visc * ...
                   ( Ue(i)^5.0 + Ue(i+1)^5.0 ) *  ...
                   ( 0.25 * ( len(i) + len(i+1) ) ) ;
      
      theta(i)  = ( theta2Ue6(i) / Ue(i)^6.0 )^0.5 ;
      ReTheta(i) = Ue(i) .* theta(i) / freeStream.kin_visc ;
      lambda(i) = theta(i)^2.0 * dUedx(i) / freeStream.kin_visc ; 
      ell(i) = thwaites_ell( lambda(i) ) ;
      H(  i) = thwaites_H(   lambda(i) ) ;
      delta(i) = theta(i) * H(i) ;
      cf( i) = 2.0 * ell(i) / ReTheta(i) ;

      %%%% % Granville's method (WRONG) 
      %%%% ReThetaTr = ( 54.2124/(H(i)*(H(i)-2.48)) + 31.6 / H(i) ) * ( H(i) >= 2.591 ) + ...
      %%%%             ( 520.0/H(i) + 2.5e+6/H(i) * ( 1/H(i)-1/2.591 )^1.95 ) * ( H(i) < 2.591 ) ;
      % Michel's criterion for transition ---
      % ReThetaTr = 1.174 * ( 1.0 + 22400.0 / Res(i) ) * Res(i) ^ 0.46 ;
      ReThetaTr = 1.35 * ( 1.0 + 22400.0 / Res(i) ) * Res(i) ^ 0.46 ;

      if ( ReTheta(i) > ReThetaTr )
          regime     = 'turbulent' ;
          transition = 'transition' ;
          H(i) = 1.35 ;
          delta(i) = theta(i) * H(i) ;
          H1(i) = head_HtoH1(H(i)) ;
          UeThetaH1(i) = Ue(i) * theta(i) * H1(i) ;
          i_trans = i ; 
      end 
      

   elseif( strcmp(regime,'turbulent') )

      dx = ( 0.5 * ( len(i) + len(i+1) ) ) ;
      UeThetaH1(i) = UeThetaH1(i+1) + dx * ...
          Ue(i+1) * 0.0306 / ( H1(i+1) - 3.0 )^0.6169 ;
      theta(i) = theta(i+1) + dx * ...
          ( 0.5 * cf(i+1) - dUedx(i+1)/Ue(i+1) * ( 2.0 + H(i+1) ) * theta(i+1) ) ;

%     lambda(i) = theta(i)^2.0 * dUedx(i) / freeStream.kin_visc ; 

      H1(i) = max( UeThetaH1(i) / ( Ue(i) * theta(i) ) , 3.0 + 1e-3 ) ;
      H(i) = head_H1toH(H1(i)) ;
      delta(i) = H(i) * theta(i) ;

      ReTheta(i) = Ue(i) .* theta(i) / freeStream.kin_visc ;
      cf(i) = ludweig_tillman_cf(H(i),ReTheta(i)) ;
      
      transition = 'occurred' ;

   end

end


%  Drag ---
tauW = 0.5 * freeStream.rho * (Ue').^2 .* cf ;

dF_visc = ( ones(2,1) * ( tauW' .* len .* vTi./Ue ) ) .* tvers ;
F_visc = sum(dF_visc') ;

L_visc = - F_visc(1) * sin(freeStream.alpha) + F_visc(2) * cos(freeStream.alpha) ;
D_visc =   F_visc(1) * cos(freeStream.alpha) + F_visc(2) * sin(freeStream.alpha) ;

dF_pres = - ( ( P' .* len' ) * ones(1,2) )' .* nvers ;
F_pres  = sum(dF_pres') ;

L_pres = - F_pres(1) * sin(freeStream.alpha) + F_pres(2) * cos(freeStream.alpha) ;
D_pres =   F_pres(1) * cos(freeStream.alpha) + F_pres(2) * sin(freeStream.alpha) ;

L = L_pres + L_visc ;
D = D_pres + D_visc ;

cL = L / ( 0.5 * freeStream.rho * freeStream.v^2 ) 
cD = D / ( 0.5 * freeStream.rho * freeStream.v^2 ) 

 
figure 
subplot(3,1,1), plot(rrc(1,:),delta,'-','LineWidth',1), ylabel('\delta') , grid on 
subplot(3,1,2), plot(rrc(1,:),theta,'-','LineWidth',1), ylabel('\theta') , grid on
subplot(3,1,3), plot(rrc(1,:),    H,'-','LineWidth',1), ylabel('H')      , grid on , xlabel('x')
  
% figure ; hold on
% plot(rrc(:,1),dPdx ,'LineWidth',2)
% ylabel('dP/dx')
% grid on, hold off
 
figure
subplot(1,2,1), plot(rrc(1,:),lambda,'-','LineWidth',1), ylabel('\lambda') , grid on , xlabel('x')
subplot(1,2,2), plot(rrc(1,:),   ell,'-','LineWidth',1), ylabel('l')       , grid on , xlabel('x')

figure
subplot(1,2,1), plot(rrc(1,:),cf  ,'-','LineWidth',1), ylabel('C_F')   , xlabel('x'), grid on
subplot(1,2,2), plot(rrc(1,:),tauW,'-','LineWidth',1), ylabel('\tau_W'), xlabel('x'), grid on


[ delta(1) , delta(end) ]






