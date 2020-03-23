% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ %
%                                                                              %
% Hess-Smith method                                                            %
% for 2-dimensional steady incompressible irrotational flow around an airfoil  %
%                                                                              %
% todo:                                                                        %
%  - home/classroom:                                                           %
%    1. check build_geometry() function                                        %
%  - classroom:                                                                %
%    2. write build_linsys() function and uncomment its respective code lines  %
%    3. uncomment the following lines and fill the lines indicated by <<<<     %
%                                                                              %
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ %

close all ; clear all ; clc

deg2rad = pi/180 ; % <- this variable should be global (...)
tol = 1.0e-9 ;

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
airfoil.refPoint     = [ 0.25 ; 0.0 ] ;
airfoil.nChordPanels = 20  ;

% === Build Geometry ===
% Build the required arrays and structures
[ rr , ee , ii_te , elems , nelems , npoints ] = build_geometry( airfoil ) ;

% Fill some matrices with the elems info for plots
rrc  = zeros(2,nelems) ; tvers= zeros(2,nelems) ; nvers= zeros(2,nelems) ;
for ie = 1 : nelems
   rrc(  :,ie) = elems(ie).cen ;
   tvers(:,ie) = elems(ie).tver ;
   nvers(:,ie) = elems(ie).nver ;
end

% -> check geometry
figure ; hold on
plot(  rr (1,:),rr (2,:),'-x','LineWidth',2) , axis equal , grid on
plot(  rrc(1,:),rrc(2,:), 'o','LineWidth',2)
quiver(rrc(1,:),rrc(2,:),nvers(1,:),nvers(2,:),'LineWidth',2)
% quiver(rrc(1,:),rrc(2,:),tvers(1,:),tvers(2,:),'LineWidth',2)
title('Airfoil geometry'), xlabel('x'), ylabel('y')
hold off


% === Build the linear system ===
[ A , b , Au , Av ] = build_linsys( freeStream , elems , ii_te )  ;

% -> check linsys
figure('Position', [200 300 900 400] )
subplot(1,2,1), spy(A)    , title(' "Sparsity" pattern ')
subplot(1,2,2), imshow(A) , title(' AIC ') , colorbar

% % === Solve the linear system ===
% vector vv contains the intensity of the singularities used to model solid bodies
% vv = A \ b;

% % === Retrieve physical quantities === 
% % x, y components of the velocity field on the body ----------------------------
% u = <<<<<
% v = <<<<<

% % Tangent and normal components of the velocity field w.r.t. the surface -------
% vTi = <<<<< 
% vNi = <<<<< 
% % -> check that vNi is = 0 (approximately)
% if ( max(abs(vNi)) > tol )
%   error('\n ERROR: Normal velocity larger than tolerance, tol: %e\n\n', tol) ;
% end

% % Pressure coefficient ---------------------------------------------------------
% % -> check that cP is approximately equal to 1.0 at the stagnation point
% % -> close to the leading edge
% cP = <<<<<

% % Plot cP -----
% figure ; hold on
% plot(rr(1,:),rr(2,:),'-') 
% plot(rrc(1,:),cP,'-o')
% grid on, axis equal, hold off

% % ...
% % ... further analyses:
% % ... - control volume;
% % ... - propagation of uncertainty: from uncertainty on velocity measurements
% % ...   to uncertainty of forces (by means of RSS);
% % ... - ... whatever you like (having in mind the hypotesis of the mathematical
% % ...       model used to describe the physical process) ...
% % ...


