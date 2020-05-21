% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ %
%                                                                              %
% Hess-Smith method                                                            %
%  for 2-dimensional steady incompressible irrotational flows around airfoils  %
%                                                                              %
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ %

close all ; clear all ; clc

% ===============================================================
% Execution controls
% ===============================================================
check_geometry = false; stop_after_check_geometry = false;
check_velocity = false; stop_after_check_velocity = false;

filen_cp    = './output/re+1e6_al+0_cp.dat'; 
filen_bl    = './output/re+1e6_al+0_bl.dat';
filen_polar = './output/re+1e6_al+0_polar.dat';

% ===============================================================
% Input
% ===============================================================
deg2rad = pi/180.0 ;

% Viscous corrections
viscous_corrections = true ; % Thwaites(laminar)+Head(turbulent) b.l. computation
integral_balance    = false;

% === Free-stream conditions ===
freeStream.P        = 0.0   ;          % freestream pressure
freeStream.rho      = 1.0   ;          % freestream density
freeStream.v        = 1.0   ;          % freestream velocity absolute value
freeStream.alpha    = 0.0 * deg2rad ;  % freestream velocity direction
freeStream.vvec     = freeStream.v * [ cos(freeStream.alpha) ; ...
                                       sin(freeStream.alpha) ] ;
%> Reynolds number and viscosity for viscous corrections
Re_c1 = 1.e+6;
freeStream.dyn_visc = freeStream.rho * freeStream.v / Re_c1 ;
freeStream.kin_visc = freeStream.dyn_visc / freeStream.rho ;

% === Airfoil input ===
% Airfoil structure definition
airfoil.id           = 1 ;
airfoil.airfoil_str  = 'NACA0012' ;
airfoil.chord        = 1.0 ;
airfoil.theta        = 0.0 * deg2rad ;
airfoil.xcRefPoint   = 0.25 ;
airfoil.refPoint     = [ 0.0 ; 0.0 ] ;
airfoil.nChordPanels = 30 ;

% Reference Reynolds number
Reynolds = freeStream.rho * freeStream.v * airfoil.chord / freeStream.dyn_visc ;


% ===============================================================
% Geometry
% ===============================================================
%> Read airfoil structure and build geometry data:
% elems is an array of objects, containing the geometrical data of
% the elems(or panels): elems(ie).ver1, .ver2, .cen, .len, .tver, .nver
[ ee , rr , ee_te , elems , nelems , npoints ] = build_geometry( airfoil ) ;

% Fill some auxiliary arrays for plots
rrc  = zeros(2,nelems); len  = zeros(1,nelems); 
tvers= zeros(2,nelems); nvers= zeros(2,nelems);
for ie = 1 : nelems
   rrc(  :,ie) = elems(ie).cen ;
   len(    ie) = elems(ie).len ;
   tvers(:,ie) = elems(ie).tver ;
   nvers(:,ie) = elems(ie).nver ;
end

%> Some checks
if ( check_geometry )
  check_geometry_script
  if ( stop_after_check_geometry ); stop; end
end

if ( check_velocity )
  check_velocity_script
  if ( stop_after_check_velocity ); stop; end
end

% ===============================================================
% Linear aerodynamic problem: Hess-Smith method
% ===============================================================
% === Build the linear system ===
% linear system: A*sing_intensity = b -> sing_intensity = A \ b
% Au, Av: auxiliary matrices for on-body velocity computation
% s.t. onbody velocity: u_onbody = Uinf_x + Au * sing_intensity
%                       v_onbody = Uinf_x + Av * sing_intensity
[ A , b , Au , Av ] = build_linsys( freeStream , elems , ee_te )  ;

% === Solve the linear system ===
% vector vv contains the intensity of the singularities used to model solid bodies
sing_intensity = A \ b;


% ===============================================================
% Retrieve physical quantities
% ===============================================================
%> === On-body velocity ===
%> x and y components
u = Au * sing_intensity + freeStream.v*cos(freeStream.alpha);
v = Av * sing_intensity + freeStream.v*sin(freeStream.alpha);
uvec = [ u , v ] ;

%> Tangent and normal components
vTi =  sum( uvec' .* tvers ) ;
vNi =  sum( uvec' .* nvers ) ;
if ( max(abs(vNi)) > 1.0e-6 ) % check if u.n = 0 at the control points
  print('Error: something could be wrong.')
  print('Boundary conditions maby be not satified, since')
  print('  max(abs(vNi)) = %6.3e \n', max(abs(vNi)) )
end

%> === Pressure coefficient ===
cP = 1.0 - vTi.^2 ./freeStream.v^2;

%> === Integral aerodynamic coefficients ===
%> cFx, cFy
cF = -sum( ( ones(2,1)*( cP .* len ) ) .* nvers, 2) / airfoil.chord ;
%> cL, cD
al = freeStream.alpha;
cA = [ cos(al) , sin(al) ; -sin(al) , cos(al) ] * cF ;

%> === Output ===
%> Write integral coefficients
fprintf('Inviscid forces \n')
fprintf(' cfx, cfy = %12.6f, %12.6f \n', cF(1), cF(2))
fprintf(' cd , cl  = %12.6f, %12.6f \n', cA(1), cA(2))

%> Lift using Kutta-Jukowski theorem
Gam = sing_intensity(end) * sum(len);
lKJ = -freeStream.rho * freeStream.v * Gam;
clKJ = lKJ / (0.5*freeStream.rho*freeStream.v^2*airfoil.chord);
fprintf(' Kutta Jukowski \n')
fprintf(' KJ   cl  =               %12.6f \n', clKJ)

%> Plot cP
figure ; hold on
plot(rr(1,:),rr(2,:),'-','LineWidth',2) 
plot(rrc(1,:),cP,'-o','LineWidth',2)
grid on, axis equal, hold off

% ===============================================================
% Integral balance
% ===============================================================
if ( integral_balance )
  force_integral_balance
end


% ===============================================================
% Viscous correction
% ===============================================================
if ( viscous_corrections )
  thwaites_method
end

% ===
% Save to files
% ===
matcp = [ rrc(1,:)' , cP' ] ;
matbl = [ rrc(1,:)' , cf, tauW, delta, theta, H ] ;
matcf = [ cL, cD, ...
          delta([1,size(delta,1)])' , ...
          theta([1,size(theta,1)])'  ];  

save(filen_cp, 'matcp', '-ascii')
save(filen_bl, 'matbl', '-ascii')
