% Hess-Smith Method for steady flow around a 2d airfoil
%
% Sketch of the Hess-Smith method:
% 1. the body (here, the airfoil) is divided into N panels
% 2. on each panel, a constant-strength line distribution of
%    source and irrotational vortex is assumed, s.t on the
%    i_th panel: q_j is the strength of the source
%                g_j is the strength of the vortex
% -> Two unknowns per each panel have been introduced
%    HESS-SMITH METHOD: the vortex strength g_i is the same
%    on all the panels (g_j = g, for j=1:N).
%    We end up with N+1 unknowns: q_i (j=1:N) and g.
%    We need N+1 equations to get a determined linear system
% 3. N+1 equations can be written:
%    N equations: b.c. assigned in the centre of the panels
%    1 equation : Kutta condition
%
% Since the problem is linear, the flow is obtained as a 
% combination of the free-stream and the sources/vortices
% distributed over the panels                            
% 
%   u(r) = U_infty + sum_{j=1:N} us_j(r) + uv(r)                (1)
%
%  where U_infty is the free-stream velocity
%        us_j(r) is the velocity induced by the j_th source
%                in the point r
%        uv  (r) is the velocity induced by the all the 
%                vortices in the point r
% 
% The intensity of the sources and the vortices can be explicitly
% written in eqn. (1)
%
%   u(r) = U_infty + sum_{j=1:N} us1_j(r) * q_j + uv1(r) * g    (2)
%
% being
%   us_j(r) = us1_j(r) * q_j
%   uv  (r) = uv1  (r) * g  
% where us1_j and uv1 are the velocity induced by a unitary-strength
% source or vortex. 
%
% An abstract description of the linear system is given.
% first N equations: b.c. on the i_th collocation point r_i, i=1:N
%   0 = n_i . u(r_i) = n_i . U_infty + ... 
%                      n_i . sum_{j=1:N} us1_j(r_i) * q_j
%                      n_i . uv1{r_i} * g  
% last equation: Kutta condition. One way to assign Kutta condition is
%   u(r_1) . t_1 = - u(r_N) . t_N
%
% The unknown intensity of the sources and the vortex are listed in 
%  the unknonw vector x = [ q_1 ; q_2 ; ... ; q_N ; g ] and the N+1 
%  equations above are written as the following linear system
%
%
% [                            |             ][ q_1 ]     [           ]
% [                            |             ][ q_2 ]     [           ]
% [                            |             ][ ... ]     [           ]
% [      AIC_bc_s(i,j)         | AIC_bc_v(i) ][     ]     [ RHS_bc(i) ]
% [                            |             ][     ]  =  [           ]
% [                            |             ][     ]     [           ]
% [                            |             ][     ]     [           ]
% [                            |             ][ q_N ]     [           ]
% [----------------------------+-------------][-----]     [-----------]
% [      AIC_kutta_s(j)        | AIC_kutta_v ][  g  ]     [ RHS_kutta ] 
%                                                            
%                                                             in this script 
% where AIC_bc_s(i,j) = n_i . us1_j(r_i)                   ->   A(1:2n,1:2n) 
%       AIC_bc_v(i)   = n_i . uv1  (r_i)                   ->   A(1:2n,2n+1) 
%       RHS_bc(i)     = -n_i . U_infty                     -> rhs(1:2n)
%       AIC_kutta_s(j)= t_1 . us1_j(r_1) + t_N . us1_j(r_N)->   A(2n+1,1:2n)
%       AIC_kutta_v   = t_1 . uv1  (r_1) + t_N . uv1  (r_N)->   A(2n+1,2n+1)
%       RHS_kutta    =-(t_1+t_N) . U_infty                 -> rhs(2n+1)
%
% 

clear all ; clc ; close all ; 

% Free stream data
rho   = 1.225;         % density
Vinf  = 30;            % velocity (abs value)
aldeg = 0.0;           % AOA [deg]
alpha = aldeg*pi/180;  % AOA [rad]

% Build naca 4-digits airfoil. Ex: NACA 2412  MPSS
M  = 2;   % 0 <= M <= 9    % max camber  ()
P  = 4;   % 0 <= P <= 9    % chord coord. of the max camber
SS = 12;  % 0 <= SS<= 99   % thickness
if ( M == 0 ) ; P = 0   ; end  % symm.airfoil -> no camber
% but analytical eqns for naca 4-digit airfoil become singular if P=0
if ( P == 0 ) ; P = 0.1 ; end  %
fprintf('\nNACA %d%d%2d \n\n',M,P,SS)

c  = 1.0 ; % chord
n  = 100 ; % n.panel / 2

[x,y] = naca4digit(M,P,SS,c,n); % x,y (1,1:2*n+1)  
r = [ x ; y ] ;

% Panels: the i-th panel goes from the i_th to the (i+1)_th point
% panel length
len = sqrt((x(2:2*n+1)-x(1:2*n)).^2 + (y(2:2*n+1)-y(1:2*n)).^2);
% normal and tangent unit vector
sintheta = (y(2:2*n+1)-y(1:2*n))./len;
costheta = (x(2:2*n+1)-x(1:2*n))./len;
theta    = atan2(sintheta,costheta);
nvers = [-sintheta; costheta];
tvers = [ costheta; sintheta];
% centre of the panels 
% they are used as collocation points, points where the b.c. are assigned
xc = (x(2:2*n+1)+x(1:2*n))./2;
yc = (y(2:2*n+1)+y(1:2*n))./2;
rc = [ xc ; yc ] ;

% Build the linear system and some auxiliary matrix
% uvi  = zeros(2*n,1);
% vvi  = zeros(2*n,1);

A = zeros(2*n+1);
rhs = zeros(2*n+1,1);
%
Au = zeros(2*n,2*n+1);
Av = zeros(2*n,2*n+1);

for ii = 1 : 2*n   % loop over the collocation points

  % Initialise the contribution of the vortex on the i_th panel equal to zero
  uvi = 0 ;
  vvi = 0 ;
  for jj = 1 : 2*n   % loop over the panels ( source intensity )

    r1 = [ rc(:,ii) - r(:,jj  ) ] ; 
    r2 = [ rc(:,ii) - r(:,jj+1) ] ; 
    
    if ( ii == jj )
      s = 1 ; 
      b = pi ;
    else
      s = norm(r2)/norm(r1) ;
      b = atan2( r1(1)*r2(2) - r2(1)*r1(2) , ...
                 r1(1)*r2(1) + r1(2)*r2(2)       ) ;
    end

    usij_ = [ -log(s);      b  ] / (2*pi) ;
    uvij_ = [     -b ; -log(s) ] / (2*pi) ;
    
    usij =   usij_' * nvers(2:-1:1,jj) ;
    vsij =   usij_' * tvers(2:-1:1,jj) ;

    uvij =   uvij_' * nvers(2:-1:1,jj) ;
    vvij =   uvij_' * tvers(2:-1:1,jj) ;

    % Update the influence AIC_bc_v of the vortex
    uvi = uvi + uvij ;
    vvi = vvi + vvij ;

    % Fill AIC_bc_s
    A (ii,jj) = nvers(:,ii)' * [ usij ; vsij ] ;
    % Fill AIC_kutta_s ( only if ii = 1 or 2*n )
    if ( ( ii == 1 ) || ( ii == 2*n ) )
      A(2*n+1,jj) = A(2*n+1,jj) + tvers(:,ii)' * [ usij ; vsij ] ;
    end 

    % Fill auxiliary matrices to recover the velocity field on the airfoil (sources)
    Au(ii,jj) = usij ;
    Av(ii,jj) = vsij ;
    
  end

  % Fill AIC_bc_v and RHS_bc
  A( ii,2*n+1) = nvers(:,ii)' * [ uvi ; vvi ] ;
  rhs( ii)       = - Vinf * [ cos(alpha) ; sin(alpha)]' * nvers(:,ii) ;
  % Fill AIC_kutta_v and RHS_kutta ( only if ii = 1 or 2*n )
  if ( ( ii == 1 ) || ( ii == 2*n ) )
    A(2*n+1,2*n+1) = A(2*n+1,2*n+1) + tvers(:,ii)' * [ uvij ; vvij ] ;
    rhs(2*n+1) = rhs(2*n+1) - Vinf * [ cos(alpha) ; sin(alpha) ]' * tvers(:,ii) ;
  end
  % Fill auxiliary matrices to recover the velocity field on the airfoil (vortex)
  Au(ii,2*n+1) = uvi ;
  Av(ii,2*n+1) = vvi ;
  
end

% Solve the linear system sol = [ q_1 ; q_2 ; ... ; q_N ; g ] ;
sol = A \ rhs;

% Recover the velocity field
u = Au * sol + Vinf*cos(alpha);
v = Av * sol + Vinf*sin(alpha);

% tangent and normal velocity in the collocation points
vTi =  u.*tvers(1,:)' + v.*tvers(2,:)' ;
vNi =  u.*nvers(1,:)' + v.*nvers(2,:)' ;

fprintf(' check the b.c. u.n=0 at the collocation points: \n') ;
fprintf(' max(abs(u.n)) = %8.4e \n\n',max(abs(vNi))) ;

% Pressure coefficient
cP = 1 - vTi.^2 ./Vinf^2;
% P - P_infty: cp = ( P - P_infty ) / ( 0.5 * rho * V_infy^2 )
dP = 0.5 * rho * Vinf^2 .* cP;
% Lift and lift coefficient (pressure integration)
L  = sum (dP .* (-nvers(2,:)'*cos(alpha)+nvers(1,:)'*sin(alpha)) .* len');
cL = L ./ (0.5 * rho * Vinf^2 * c);
% Drag and drag coefficient (pressure integration)
D  = sum (dP .* ( nvers(1,:)'*cos(alpha)+nvers(2,:)'*sin(alpha)) .* len');
cD = D ./ (0.5 * rho * Vinf^2 * c);

% Some output
fprintf(' AOA        : alpha = %6.2f °\n',aldeg*180/pi)
%fprintf(' aerodynamic coefficients (pressure integration) \n')
fprintf(' lift coeff : cl    = %6.3f \n',cL)
fprintf(' drag coeff : cd    = %6.3f \n',cD)
% Dimensional output ?
% fprintf(...)
% fprintf(...)

fprintf('\n')

% Some plots
% pressure coefficient
figure(1)
plot(xc(1:n)./c,cP(1:n),'r',xc(n:end)./c,cP(n:end),'b','LineWidth',1), axis ij
legend('ventre','dorso')
hold on
plot(x./c,-y./c,'k','LineWidth',1)
ylabel('c_P'),xlabel('x/c'),axis([-0.1 1.1 -1.5 1.1]), axis equal

