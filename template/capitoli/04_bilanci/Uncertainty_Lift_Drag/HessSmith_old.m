% % Hess-Smith Method

close all
clear all
clc

% Dati della corrente esterna +++++++++++++++++++++++++++++++++++ 
rho   = 1.225;  Vinf  = 30;  alpha = 0.0*pi/180;

% Dati del profilo: NACA 4-digit 5-digit airfoil ++++++++++++++++
nDigit = 4;

% NACA 2412  : MPSS -------
M  = 2; P  = 4; SS = 12;
% NACA 23012 : ddd SS -----
ddd = 230;

c = 0.3;        % corda
nplot = 100;    % numero di punti per i grafici
n = 30 ;        % metà del numero di pannelli (nPan_dorso=nPan_ventre)

% Punti per il metodo di Hess-Smith
if nDigit == 4
    [x,y] = naca4digit(M,P,SS,c,n);
elseif nDigit == 5
    [x,y] = naca5digit(230,12,c,n);
end

% Save coordinates ------------------------
% M = [x',y'];
% save('naca2412bis','M','-ascii')

% Costruzione matrici per il metodo di Hess-Smith +++++++++++++++

% Grandezze caratteristiche dei pannelli
% Dimensioni pannelli
len = sqrt((x(2:2*n+1)-x(1:2*n)).^2 + (y(2:2*n+1)-y(1:2*n)).^2);
% Versori normale e tangente
sintheta = (y(2:2*n+1)-y(1:2*n))./len;
costheta = (x(2:2*n+1)-x(1:2*n))./len;
theta    = atan2(sintheta,costheta);
nvers = [-sintheta; costheta];
tvers = [ costheta; sintheta];
% Centri dei pannelli
xc = (x(2:2*n+1)+x(1:2*n))./2;
yc = (y(2:2*n+1)+y(1:2*n))./2;

% Calcolo di Matrici ausiliari
x1  = zeros(2*n);
y1  = zeros(2*n);
x2  = zeros(2*n);
y2  = zeros(2*n);
sij = zeros(2*n);
sinbij = zeros(2*n);
cosbij = zeros(2*n);
bij = zeros(2*n);
usijstar = zeros(2*n);
vsijstar = zeros(2*n);
uvijstar = zeros(2*n);
vvijstar = zeros(2*n);
usij = zeros(2*n);
vsij = zeros(2*n);
uvij = zeros(2*n);
vvij = zeros(2*n);
uvi  = zeros(2*n,1);
vvi  = zeros(2*n,1);

A = zeros(2*n+1);
b = zeros(2*n+1,1);

Au = zeros(2*n,2*n+1);
Av = zeros(2*n,2*n+1);

for ii = 1 : 2*n
    for jj = 1 : 2*n
        x1(ii,jj) = xc(ii) - x(jj);    y1(ii,jj) = yc(ii) - y(jj);    r1 = [x1(ii,jj);y1(ii,jj)];
        x2(ii,jj) = xc(ii) - x(jj+1);  y2(ii,jj) = yc(ii) - y(jj+1);  r2 = [x2(ii,jj);y2(ii,jj)];
        
        if ii == jj
            sij(ii,jj) = 1;
            bij(ii,jj) = pi;
        else
            sij(ii,jj) = norm(r2)/norm(r1);
            sinbij(ii,jj) = (x1(ii,jj)*y2(ii,jj) - x2(ii,jj)*y1(ii,jj))/ norm(r1) / norm(r2);
            cosbij(ii,jj) = (x1(ii,jj)*x2(ii,jj) + y1(ii,jj)*y2(ii,jj))/ norm(r1) / norm(r2);
            bij(ii,jj) = atan2(sinbij(ii,jj),cosbij(ii,jj));
        end
        
        usijstar(ii,jj) = -log(sij(ii,jj))/(2*pi);
        vsijstar(ii,jj) =  bij(ii,jj)/(2*pi);
        
        uvijstar(ii,jj) = -bij(ii,jj)/(2*pi);
        vvijstar(ii,jj) = -log(sij(ii,jj))/(2*pi);
        
        usij(ii,jj) = usijstar(ii,jj)*costheta(jj) - vsijstar(ii,jj)*sintheta(jj);
        vsij(ii,jj) = usijstar(ii,jj)*sintheta(jj) + vsijstar(ii,jj)*costheta(jj);
        
        uvij(ii,jj) = uvijstar(ii,jj)*costheta(jj) - vvijstar(ii,jj)*sintheta(jj);
        vvij(ii,jj) = uvijstar(ii,jj)*sintheta(jj) + vvijstar(ii,jj)*costheta(jj);
        
        A (ii,jj)   = -sintheta(ii)*usij(ii,jj) + costheta(ii)*vsij(ii,jj);
        Au(ii,jj)   = usij(ii,jj);
        Av(ii,jj)   = vsij(ii,jj);
    end
    uvi(ii) = sum(uvij(ii,:));
    vvi(ii) = sum(vvij(ii,:));
    
    A(ii,2*n+1) = -sintheta(ii)* uvi(ii) + costheta(ii)* vvi(ii);
    b(ii)       =  sintheta(ii)* Vinf*cos(alpha) - costheta(ii)* Vinf*sin(alpha);
    
    Au(ii,2*n+1) = uvi(ii);
    Av(ii,2*n+1) = vvi(ii);
end
A(2*n+1 , 1:2*n) = costheta(1)   .*usij(1,1:2*n)   + sintheta(1)   .*vsij(1,1:2*n) + ...
                   costheta(2*n) .*usij(2*n,1:2*n) + sintheta(2*n) .*vsij(2*n,1:2*n) ;
A(2*n+1 , 2*n+1) = costheta(1)   .*uvi(1)   + sintheta(1)   .*vvi(1) + ...
                   costheta(2*n) .*uvi(2*n) + sintheta(2*n) .*vvi(2*n) ;
b(2*n+1)         = -(costheta(1)   .*Vinf * cos(alpha)  + sintheta(1)   .*Vinf * sin(alpha) + ...
                     costheta(2*n) .*Vinf * cos(alpha)  + sintheta(2*n) .*Vinf * sin(alpha));

% Solve the linear system : unk = [ q_1 ; ... ; q_2n ; gamma ]
unk = A\b;

% Compute velocity components con the pannels from unk
u = Au * unk + Vinf*cos(alpha);   % x-component
v = Av * unk + Vinf*sin(alpha);   % y-component
modu = sqrt(u.^2 + v.^2);
angolo = atan2(v,u);

% Projection of the velocity parallel and normal to the panels
vTi =  u .* costheta' + v .* sintheta';
vNi = -u .* sintheta' + v .* costheta'; 

% Pressure coefficients, pressure, lift and drag
cP = 1 - vTi.^2 ./Vinf^2;
dP = 0.5 * rho * Vinf^2 .* cP;

L  = sum (dP .* (-nvers(2,:)'*cos(alpha)+nvers(1,:)'*sin(alpha)) .* len');
cL = L ./ (0.5 * rho * Vinf^2 * c);

D  = sum (dP .* (nvers(1,:)'*cos(alpha)+nvers(2,:)'*sin(alpha)) .* len');
cD = D ./ (0.5 * rho * Vinf^2 * c);

% Useless arrays -----
% uvct = [((1:2*n))' ,u, v, modu, angolo, theta'];
% ang  = [sintheta' costheta' theta'.*180/pi];
% Result = [xc',yc',dP./1000];

figure(1)
plot(xc(1:n)./c,cP(1:n),'r',xc(n:end)./c,cP(n:end),'b','LineWidth',1), axis ij
legend('ventre','dorso')
hold on
plot(x./c,-y./c,'k','LineWidth',1)
ylabel('c_P'),xlabel('x/c'),axis([-0.1 1.1 -1.5 1.1]), axis equal


% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Lift, drag and uncertainty from box
[ x_ctr , y_ctr , u_ctr , v_ctr ] = control_volume ( x , y , Vinf , alpha , unk ) ;

Um_ctr = ( u_ctr.^2 + v_ctr.^2 ) .^ 0.5 ;
P0_ctr = 0.0 ; % 1e5 ;
P_ctr  = P0_ctr + 0.5 * rho * ( Vinf^2 -  Um_ctr.^2 ) ; 
% Lift and drag from integral balance on the control volume surface
% Augment box vectors for the integrals
x_c = [ x_ctr' ; x_ctr(1) ];
y_c = [ y_ctr' ; y_ctr(1) ];
u_c = [ u_ctr  ; u_ctr(1) ];
v_c = [ v_ctr  ; v_ctr(1) ];
p_c = [ P_ctr  ; P_ctr(1) ];

% vector \hat{n}*dS
ns_c = [ - ( y_c(2:end) - y_c(1:end-1) ) , ( x_c(2:end) - x_c(1:end-1) ) ] ;
uu_c = [ u_c , v_c ] ;

% No pressure terms : L = - int_S{ rho u u.n }
dR_app = - 0.5*rho*(uu_c(2:end  ,:)'.*(ones(2,1)*sum(uu_c(2:end  ,:)'.*ns_c'))+ ...
                    uu_c(1:end-1,:)'.*(ones(2,1)*sum(uu_c(1:end-1,:)'.*ns_c'))) ;

R_app = sum(dR_app,2) ;

% With pressure terms : L = - int_S{ rho u u.n } - int_S{Pn}
dR = -0.5 * rho*(uu_c(2:end  ,:)'.*(ones(2,1)*sum(uu_c(2:end  ,:)'.*ns_c'))+ ...
                 uu_c(1:end-1,:)'.*(ones(2,1)*sum(uu_c(1:end-1,:)'.*ns_c'))) ...
    - 0.5 *( ones(2,1)*p_c(2:end,:)' + ones(2,1)*p_c(1:end-1)' ) .* ns_c' ;

R = sum(dR,2) ;

% Uncertainty +++++++++++
i_sigP = 0.0  ;
dsigUx = 0.05 ;
dsigUy = 0.05 ;
dsigD = 0.5 * rho.^2 * ( ...
        ( sum(uu_c(2:end,:)'.*ns_c') + uu_c(2:end,1)'.*ns_c(:,1)' + ...
                             i_sigP .* uu_c(2:end,1)'.*ns_c(:,1)' ).^2     .* dsigUx .^2 + ...
        ( uu_c(2:end,2)'.*ns_c(:,1)'                              + ... 
                             i_sigP .* uu_c(2:end,2)'.*ns_c(:,1)' ).^2     .* dsigUy .^2 + ...
        ( sum(uu_c(1:end-1,:)'.*ns_c') + uu_c(1:end-1,1)'.*ns_c(:,1)' + ...
                             i_sigP .* uu_c(1:end-1,1)'.*ns_c(:,1)' ).^2   .* dsigUx .^2 + ...
        ( uu_c(1:end-1,2)'.*ns_c(:,1)'                                + ...
                             i_sigP .* uu_c(1:end-1,2)'.*ns_c(:,1)' ).^2   .* dsigUy .^2 ) ;
        
dsigL = 0.5 * rho.^2 * ( ...
        ( sum(uu_c(2:end,:)'.*ns_c') + uu_c(2:end,2)'.*ns_c(:,2)' + ...
                             i_sigP .* uu_c(2:end,2)'.*ns_c(:,2)').^2      .* dsigUy .^2 + ...
        ( uu_c(2:end,1)'.*ns_c(:,2)'                              + ...
                             i_sigP .* uu_c(2:end,1)'.*ns_c(:,2)' ).^2     .* dsigUx .^2 + ...
        ( sum(uu_c(1:end-1,:)'.*ns_c') + uu_c(1:end-1,2)'.*ns_c(:,2)' + ...
                             i_sigP .* uu_c(1:end-1,2)'.*ns_c(:,2)' ).^2   .* dsigUy .^2 + ...
        ( uu_c(1:end-1,1)'.*ns_c(:,2)'                                + ...
                             i_sigP .* uu_c(1:end-1,1)'.*ns_c(:,2)' ).^2   .* dsigUx .^2 ) ;
        

sig2L = sum(dsigL);
sig2D = sum(dsigD);

rel_err_L = abs( sig2L^0.5 / R(2) ) ;
rel_err_D = abs( sig2D^0.5 / R(1) ) ;

