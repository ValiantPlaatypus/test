% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ %
%                                                                              %
% Hess-Smith method                                                            %
%  for 2-dimensional steady incompressible irrotational flows around airfoils  %
%                                                                              %
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ %

close all ; clear all ; clc

deg2rad = pi/180 ;

% Viscous corrections
viscous_corrections = false ; % Thwaites(laminar)+Head(turbulent) b.l. computation

% === Free-stream conditions ===
freeStream.P        = 0.0   ;          % free-stream pressure
freeStream.rho      = 1.0   ;          % free-stream density
freeStream.v        = 1.0   ;          %             velocity absolute value
freeStream.alpha    = 0.0 * deg2rad ;  %             velocity direction
freeStream.vvec     = freeStream.v * [ cos(freeStream.alpha) ; sin(freeStream.alpha) ] ;
Re_c1 = 1.e+5;
freeStream.dyn_visc = freeStream.rho * freeStream.v / Re_c1 ;         % free-stream dynamic viscosity
% freeStream.dyn_visc = 1.8e-5 ;         % free-stream dynamic viscosity
freeStream.kin_visc = ...              % free-stream kinematic viscosity
                      freeStream.dyn_visc / freeStream.rho ;

% === Airfoil input ===
% Airfoil structure definition
airfoil(1).id           = 1 ;
airfoil(1).airfoil_str  = 'NACA0012' ;
airfoil(1).chord        = 1.0 ;
airfoil(1).theta        = 0.0 * deg2rad ;
airfoil(1).xcRefPoint   = 0.25 ;
airfoil(1).refPoint     = [ 0.0 ; 0.0 ] ;
airfoil(1).nChordPanels = 30 ;

% Reference Reynolds number
Reynolds = freeStream.rho * freeStream.v * airfoil(1).chord / freeStream.dyn_visc ;

% % Add other airfoils ---
% airfoil(2).airfoil_str  = 'NACA0012' ;
% airfoil(2).chord        = 0.5 ;
% airfoil(2).theta        = 7.0 * deg2rad ;
% airfoil(2).xcRefPoint   = 0.25 ;
% airfoil(2).refPoint     = [ 1.00 ; -0.1 ] ;
% airfoil(2).nChordPanels = 40 ;


% === Build Geometry ===
% OBS: all the airfoils have TE -> nAirfoils = nTE = n. Kutta conditions required
% 1. Initialise the needed objects
nelems = 0 ; npoints = 0 ;
rr = [] ; ee = [] ; ee_te = [] ; elems = [] ;
% 2. Update them by looping on the airfoils
nAirfoils = length(airfoil) ;
for ia = 1 : nAirfoils
  [ ee , rr , ee_te , elems , nelems , npoints ] = ...
            build_geometry( ee , rr , ee_te , elems , nelems , npoints , airfoil(ia) ) ;
end

% 3. Fill some matrices with the elems info, because MATLAB is meant to efficiently
% work with matrices
rrc  = zeros(nelems,2) ;
len  = zeros(nelems,1) ; 
tvers= zeros(2,nelems) ;
nvers= zeros(2,nelems) ;
for ie = 1 : nelems
   rrc(  ie,:) = elems(ie).cen ;
   len(  ie  ) = elems(ie).len ;
   tvers(:,ie) = elems(ie).tver ;
   nvers(:,ie) = elems(ie).nver ;
end

% === Build the linear system ===
[ A , b , Au , Av ] = build_linsys( freeStream , elems , ee_te )  ;

% === Solve the linear system ===
% vector vv contains the intensity of the singularities used to model solid bodies
vv = A \ b;


% === Retrieve physical quantities === 
% x, y components of the velocity field on the body ----------------------------
u = Au * vv + freeStream.v*cos(freeStream.alpha);
v = Av * vv + freeStream.v*sin(freeStream.alpha);
uvec = [ u , v ] ;
% modu = sqrt(u.^2 + v.^2);
% angolo = atan2(v,u);

% Tangent and normal components of the velocity field w.r.t. the surface -------
% !!! check that vNi is = 0 (approximately)                                  !!!
vTi =  sum( uvec' .* tvers ) ;
vNi =  sum( uvec' .* nvers ) ;

% Pressure coefficient ---------------------------------------------------------
% !!! check that cP is approximately equal to 1.0 at the stagnation point    !!!
% !!! close to the leading edge                                              !!!
cP = 1.0 - vTi.^2 ./freeStream.v^2;

fprintf('Inviscid forces \n')
cF = -sum( ( ones(2,1)*( cP .* len' ) ) .* nvers, 2) ;
al = freeStream.alpha;
cA = [ cos(al) , sin(al) ; -sin(al) , cos(al) ] * cF ;
fprintf(' cfx, cfy = %12.6f, %12.6f \n', cF(1), cF(2))
fprintf(' cd , cl  = %12.6f, %12.6f \n', cA(1), cA(2))

% Plot cP -----
figure ; hold on
plot(rr(:,1),rr(:,2),'-') 
plot(rrc(:,1),cP,'-o')
grid on, axis equal, hold off

if ( viscous_corrections )
  thwaites_method
end
