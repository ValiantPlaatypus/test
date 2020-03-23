function [x,y] = naca5digit(ddd,SS,c,n)

t = SS / 100;

xv = linspace(0,c,n+1);
xv = c/2 .*(1-cos(pi.*xv./c));

switch ddd

   case 210
     q = 0.0580;
     k = 361.4;
     
   case 220
     q = 0.1260;
     k = 51.64;
     
   case 230
     q = 0.2025;
     k = 15.957;
     
   case 240
     q = 0.2900;
     k = 6.643; 
     
   case 250
     q = 0.3910;
     k = 3.230;
     
   otherwise
     error('Errore nella chiamata a linea_media_5d: la linea media richiesta non esiste.')
end


% Spessore
ytfcn = @(x) 5.*t.*c.*(0.2969.*(x./c).^0.5 - 0.1260.*(x./c) ...
    - 0.3516.*(x./c).^2 + 0.2843.*(x./c).^3 - 0.1036.*(x./c).^4);

yt = ytfcn(xv);

% Linea Media
yc = zeros(size(xv));

for ii = 1 : n+1
    if xv(ii) <= q*c
        yc(ii) = c*(k/6 * ((xv(ii)/c)^3 - 3*q*(xv(ii)/c)^2 + q^2*(3-q)*(xv(ii)/c)));
    else
        yc(ii) = c*(k/6 * q^3 * (1 - (xv(ii)/c)));
    end
end

% Derivata della Linea Media
dyc = zeros(size(xv));

for ii = 1 : n+1
    if xv(ii) <= q*c
        dyc(ii) = k/6 * (3*(xv(ii)/c)^2 - 6*q*(xv(ii)/c) + q^2*(3-q));
    else
        dyc(ii) = -k/6 * q^3;
    end
end

% Ventre e Dorso
th = atan2(dyc,1);
xU = xv - yt.*sin(th);
yU = yc + yt.*cos(th);
xL = xv + yt.*sin(th);
yL = yc - yt.*cos(th);

x = zeros(1,2*n+1);
y = zeros(1,2*n+1);
for ii = 1 : n
   x(ii) = xL(n+2-ii);
   y(ii) = yL(n+2-ii);
end

x(n+1:2*n+1) = xU;
y(n+1:2*n+1) = yU;