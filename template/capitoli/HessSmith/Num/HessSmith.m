% % Hess-Smith Method

close all
clear all
%clc

nplot = 100;    % numero di punti per i grafici

% % Dati della corrente esterna 
rho   = 1.225;
Vinf  = 1 ;
alpha = 0.0*pi/180;

% % Dati del profilo: NACA 4-digit 5-digit airfoil
nDigit = 4;

% NACA 2412  MPSS
M  = 0;
P  = 4;
SS = 12;

% NACA 23012    ddd SS
ddd = 230;

nny = 1 ;
y_c = 100*(1:nny);
cL1v = zeros(1,nny);
cL2v = zeros(1,nny);

for iy = 1 : nny
% Profilo 1 : naca 4/5 digit
c = 1.0;        % corda
theta1= 0 * pi/180;
n = 30;         % metà del numero di pannelli (nPannelli_dorso=nPanelli_ventre)
% Profilo 2 : naca 4 digit
nAirfoil = 2;
c2  = 1.00;
Dx  = 0.0;
Dy  = y_c(iy)*c;
cAC = 0.25;
Dal = 0 * pi/180;

% Coordinate del primo profilo
if nDigit == 4
%    [x,y] = naca4digit(M,P,SS,c,n);
    mirror = 0;
    [x,y] = setAirfoil4(M,P,SS,c,n,0.0,0.0,cAC,theta1,mirror);
elseif nDigit == 5
    [x,y] = naca5digit(230,12,c,n);
end
Mat = [x',y'];
save('naca2412bis','Mat','-ascii')

% Coordinate del secondo profilo
mirror = 0;
[xa,ya] = setAirfoil4(M,P,SS,c2,n,Dx,Dy,cAC,Dal,mirror);


% Costruzione delle matrici contenenti la geometria del problema
xx = zeros(2,2*nAirfoil*n);
yy = zeros(2,2*nAirfoil*n);
xx(1,:) = [x(1:end-1) , xa(1:end-1)];
xx(2,:) = [x(2:end)   , xa(2:end)  ];
yy(1,:) = [y(1:end-1) , ya(1:end-1)];
yy(2,:) = [y(2:end)   , ya(2:end)  ];

figure
plot(x,y,xa,ya)
axis equal, grid on

% Grandezze caratteristiche dei pannelli
% Dimensioni pannelli
  len = sqrt((xx(2,:)-xx(1,:)).^2 + (yy(2,:)-yy(1,:)).^2);
% Versori normale e tangente
  sintheta = (yy(2,:)-yy(1,:))./len;
  costheta = (xx(2,:)-xx(1,:))./len;
  theta    = atan2(sintheta,costheta);
  nvers = [-sintheta; costheta];
  tvers = [ costheta; sintheta];
% Centri dei pannelli
  xc = (xx(2,:)+xx(1,:))./2;
  yc = (yy(2,:)+yy(1,:))./2;


% Calcolo di Matrici ausiliari
x1       = zeros(2*nAirfoil*n);
y1       = zeros(2*nAirfoil*n);
x2       = zeros(2*nAirfoil*n);
y2       = zeros(2*nAirfoil*n);
sij      = zeros(2*nAirfoil*n);
sinbij   = zeros(2*nAirfoil*n);
cosbij   = zeros(2*nAirfoil*n);
bij      = zeros(2*nAirfoil*n);
usijstar = zeros(2*nAirfoil*n);
vsijstar = zeros(2*nAirfoil*n);
uvijstar = zeros(2*nAirfoil*n);
vvijstar = zeros(2*nAirfoil*n);
usij     = zeros(2*nAirfoil*n);
vsij     = zeros(2*nAirfoil*n);
uvij     = zeros(2*nAirfoil*n);
vvij     = zeros(2*nAirfoil*n);
uv1i     = zeros(2*nAirfoil*n,1);
vv1i     = zeros(2*nAirfoil*n,1);
uv2i     = zeros(2*nAirfoil*n,1);
vv2i     = zeros(2*nAirfoil*n,1);

A = zeros((2*n+1)*nAirfoil);
b = zeros((2*n+1)*nAirfoil,1);

Au = zeros(2*n*nAirfoil,(2*n+1)*nAirfoil);
Av = zeros(2*n*nAirfoil,(2*n+1)*nAirfoil);

for ii = 1 : 2*nAirfoil*n
    for jj = 1 : 2*nAirfoil*n
        x1(ii,jj) = xc(ii) - xx(1,jj);    y1(ii,jj) = yc(ii) - yy(1,jj);    r1 = [x1(ii,jj);y1(ii,jj)];
        x2(ii,jj) = xc(ii) - xx(2,jj);    y2(ii,jj) = yc(ii) - yy(2,jj);    r2 = [x2(ii,jj);y2(ii,jj)];
        
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
    uv1i(ii) = sum(uvij(ii,1:2*n));
    vv1i(ii) = sum(vvij(ii,1:2*n));
    uv2i(ii) = sum(uvij(ii,1+2*n:4*n));
    vv2i(ii) = sum(vvij(ii,1+2*n:4*n));
    
    A(ii,4*n+1) = -sintheta(ii)* uv1i(ii) + costheta(ii)* vv1i(ii);
    A(ii,4*n+2) = -sintheta(ii)* uv2i(ii) + costheta(ii)* vv2i(ii);
    b(ii)       =  sintheta(ii)* Vinf*cos(alpha) - costheta(ii)* Vinf*sin(alpha);
    Au(ii,4*n+1) = uv1i(ii);
    Au(ii,4*n+2) = uv2i(ii);
    Av(ii,4*n+1) = vv1i(ii);
    Av(ii,4*n+2) = vv2i(ii);
end

A(4*n+1 , 1:2*n)     = costheta(1)        .*usij(1,1:2*n)           + sintheta(1)       .*vsij(1,1:2*n) + ...
                       costheta(2*n)      .*usij(2*n,1:2*n)         + sintheta(2*n)     .*vsij(2*n,1:2*n) ;
A(4*n+1 , 1+2*n:4*n) = costheta(1)        .*usij(1  ,1+2*n:4*n)     + sintheta(1)       .*vsij(1,1+2*n:4*n) + ...
                       costheta(2*n)      .*usij(2*n,1+2*n:4*n)     + sintheta(2*n)     .*vsij(2*n,1+2*n:4*n) ;
A(4*n+2 , 1:2*n)     = costheta(2*n+1)    .*usij(2*n+1,1:2*n)       + sintheta(2*n+1)   .*vsij(2*n+1,1:2*n) + ...
                       costheta(2*n+2*n)  .*usij(2*n+2*n,1:2*n)     + sintheta(2*n+2*n) .*vsij(2*n+2*n,1:2*n) ;
A(4*n+2 , 1+2*n:4*n) = costheta(2*n+1)    .*usij(2*n+1  ,1+2*n:4*n) + sintheta(2*n+1)   .*vsij(2*n+1,1+2*n:4*n) + ...
                       costheta(2*n+2*n)  .*usij(2*n+2*n,1+2*n:4*n) + sintheta(2*n+2*n) .*vsij(2*n+2*n,1+2*n:4*n) ;
A(4*n+1 , 4*n+1) = costheta(1)   .*uv1i(1)   + sintheta(1)   .*vv1i(1) + ...
                   costheta(2*n) .*uv1i(2*n) + sintheta(2*n) .*vv1i(2*n) ;
A(4*n+1 , 4*n+2) = costheta(1)   .*uv2i(1)   + sintheta(1)   .*vv2i(1) + ...
                   costheta(2*n) .*uv2i(2*n) + sintheta(2*n) .*vv2i(2*n) ;
A(4*n+2 , 4*n+1) = costheta(2*n+1)   .*uv1i(1+2*n)   + sintheta(2*n+1)   .*vv1i(1+2*n) + ...
                   costheta(2*n+2*n) .*uv1i(2*n+2*n) + sintheta(2*n+2*n) .*vv1i(2*n+2*n) ;
A(4*n+2 , 4*n+2) = costheta(2*n+1)   .*uv2i(2*n+1)   + sintheta(2*n+1)   .*vv2i(2*n+1) + ...
                   costheta(2*n+2*n) .*uv2i(2*n+2*n) + sintheta(2*n+2*n) .*vv2i(2*n+2*n) ;
b(4*n+1)         = -(costheta(1)   .*Vinf * cos(alpha)  + sintheta(1)   .*Vinf * sin(alpha) + ...
                     costheta(2*n) .*Vinf * cos(alpha)  + sintheta(2*n) .*Vinf * sin(alpha));
b(4*n+2)         = -(costheta(1+2*n)   .*Vinf * cos(alpha)  + sintheta(1+2*n)   .*Vinf * sin(alpha) + ...
                     costheta(4*n)     .*Vinf * cos(alpha)  + sintheta(4*n)     .*Vinf * sin(alpha));

vv = A\b;


u = Au * vv + Vinf*cos(alpha);
v = Av * vv + Vinf*sin(alpha);
modu = sqrt(u.^2 + v.^2);
angolo = atan2(v,u);

uvct = [((1:4*n))' ,u, v, modu, angolo, theta'];

ang = [sintheta' costheta' theta'.*180/pi];

vTi =  u .* costheta' + v .* sintheta';
vNi = -u .* sintheta' + v .* costheta'; 

% % Ricostruzione velocit� tangenziale
% vTi = zeros(2*n,1);
% 
% for ii = 1 : 2*n
%     vT = Vinf * cos(theta(ii) - alpha);
%     for jj = 1 : 2*n
%         vT = vT + v(jj)/(2*pi)*(sin(theta(ii)-theta(jj))*bij(ii,jj) - ...
%             cos(theta(ii)-theta(jj))*log(sij(ii,jj))) - ...
%             v(2*n+1)/2*pi * (cos(theta(ii)-theta(jj))*bij(ii,jj) + ...
%             sin(theta(ii)-theta(jj))*log(sij(ii,jj)));
%     end
%     
%     vTi(ii) = vT;
% end
% 
cP = 1 - vTi.^2 ./Vinf^2;

dP = 0.5 * rho * Vinf^2 .* cP;

Result = [xc',yc',dP./1000];

L  = sum (dP .* (-nvers(2,:)'*cos(alpha)+nvers(1,:)'*sin(alpha)) .* len');
L1 = sum (dP(1:2*n) .* (-nvers(2,1:2*n)'*cos(alpha)+nvers(1,1:2*n)'*sin(alpha)) .* len(1:2*n)' ); 
L2 = sum (dP(1+2*n:2*n+2*n) .* (-nvers(2,1+2*n:2*n+2*n)'*cos(alpha)+ ...
                                  nvers(1,1+2*n:2*n+2*n)'*sin(alpha)) .* len(1+2*n:2*n+2*n)' ); 
cL1 = L1 ./ (0.5 * rho * Vinf^2 * c );
cL2 = L2 ./ (0.5 * rho * Vinf^2 * c2);

  cL1v(iy) = cL1;
  cL2v(iy) = cL2;


fprintf('++++++++++++++++++++++++++++++++++++++++++++++++++++++++ \n')
fprintf('Posizione relativa del profilo 2 rispetto al profilo 1:  \n')
fprintf(' x0/c = %6.3f   y0/c = %6.3f \n', Dx/c,Dy/c)
fprintf('Coefficienti di portanza dei due profili:                \n')
fprintf(' a1 = %6.3f     cL1 = %6.3f \n',(theta1+alpha)*180/pi,cL1)
fprintf(' a2 = %6.3f     cL2 = %6.3f \n',(Dal   +alpha)*180/pi,cL2)
fprintf('++++++++++++++++++++++++++++++++++++++++++++++++++++++++ \n\n')


end
D  = sum (dP .* (nvers(1,:)'*cos(alpha)+nvers(2,:)'*sin(alpha)) .* len');
cD = D ./ (0.5 * rho * Vinf^2 * c);

rangx = max(max(xx)) -min(min(xx)) ;
rangy = max(cP)-min(cP);
xmin  = min(min(xx)) - 0.1*rangx;
xmax  = max(max(xx)) + 0.1*rangx;
ymin  = min(cP)- 0.1*rangy;
ymax  = max(cP)+ 0.1*rangy;
ScaleF= (ymax-ymin)/(xmax-xmin);


figure
h = plot(y_c,cL2v,'-o',y_c,cL1v,'-o',[0 15],[0.594 0.594]), grid on
set(h(1),'LineWidth',2)
set(h(2),'LineWidth',2)
set(h(3),'LineWidth',2)
xlabel('y/c'), ylabel('cL')
axis([0 15 0.48 0.60])
legend('upper','lower','free','Location','SouthEast')

figure(2)
plot(xc(1:n),cP(1:n),'r',xc(n:2*n),cP(n:2*n),'b','LineWidth',1), axis ij
ylabel('c_P'),xlabel('x/c'),axis([xmin xmax ymin ymax]), axis square, grid on
legend('ventre','dorso')
hold on
plot(xc(2*n+1:3*n),cP(2*n+1:3*n),'r',xc(3*n:4*n),cP(3*n:4*n),'b','LineWidth',1), axis ij
plot(x ,-y .*ScaleF,'k','LineWidth',1)
plot(xa,-ya.*ScaleF,'k','LineWidth',1)

%, axis equal
% hold on
% plot(x,y,'k')

%figure
%plot(xc,dP)
