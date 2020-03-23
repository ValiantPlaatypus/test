% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ %
%                                                                              %
% Hess-Smith method                                                            %
% for 2-dimensional steady incompressible irrotational flows around an airfoil %
%                                                                              %
%                                                                              %
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ %

close all ; clear all ; clc

deg2rad = pi/180 ; % <- this variable should be global (...)


% === Free-stream conditions ===
freeStream.P        = 0.0   ;          % free-stream pressure
freeStream.rho      = 1.225 ;          % free-stream density
freeStream.v        = 10.0  ;          %             velocity absolute value
freeStream.alpha    = 2.0 * deg2rad ;  %             velocity direction
freeStream.vvec     = freeStream.v * [ cos(freeStream.alpha) ; sin(freeStream.alpha) ] ;

% === Airfoil input ===
% Airfoil structure definition
airfoil.id           = 1 ;
airfoil.airfoil_str  = 'NACA0012' ;
airfoil.chord        = 1.0 ;
airfoil.theta        = 0.0 * deg2rad ;
airfoil.xcRefPoint   = 0.25 ;
airfoil.refPoint     = [ 0.0 ; 0.0 ] ;
airfoil.nChordPanels = 30  ;

% === Build Geometry ===
% OBS: all the airfoils have TE -> nAirfoils = nTE = n. Kutta conditions required
% 1. Build the required arrays and structures
[ rr , ee , ii_te , elems , nelems , npoints ] = build_geometry( airfoil ) ;

% % 2. Fill some matrices with the elems info, because MATLAB is meant to efficiently
% % work with matrices
% rrc  = zeros(2,nelems) ;
% tvers= zeros(2,nelems) ;
% nvers= zeros(2,nelems) ;
% for ie = 1 : nelems
%    rrc(  :,ie) = elems(ie).cen ;
%    tvers(:,ie) = elems(ie).tver ;
%    nvers(:,ie) = elems(ie).nver ;
% end

% % === Build the linear system ===
% [ A , b , Au , Av ] = build_linsys( freeStream , elems , ii_te )  ;
% 
% % === Solve the linear system ===
% % vector vv contains the intensity of the singularities used to model solid bodies
% vv = A \ b;
% 
% % === Retrieve physical quantities === 
% % x, y components of the velocity field on the body ----------------------------
% u = Au * vv + freeStream.v*cos(freeStream.alpha);
% v = Av * vv + freeStream.v*sin(freeStream.alpha);
% uvec = [ u , v ] ;
% % modu = sqrt(u.^2 + v.^2);
% % angolo = atan2(v,u);
% 
% % Tangent and normal components of the velocity field w.r.t. the surface -------
% % !!! check that vNi is = 0 (approximately)                                  !!!
% vTi =  sum( uvec' .* tvers ) ;
% vNi =  sum( uvec' .* nvers ) ;
% 
% % Pressure coefficient ---------------------------------------------------------
% % !!! check that cP is approximately equal to 1.0 at the stagnation point    !!!
% % !!! close to the leading edge                                              !!!
% cP = 1 - vTi.^2 ./freeStream.v^2;
% 
% % Plot cP -----
% figure ; hold on
% plot(rr(1,:),rr(2,:),'-') 
% plot(rrc(1,:),cP,'-o')
% grid on, axis equal, hold off

