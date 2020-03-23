% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ %
%                                                                              %
% Hess-Smith method                                                            %
%  for 2-dimensional steady incompressible irrotational flows around airfoils  %
%                                                                              %
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ %

close all ; clear all ; clc

deg2rad = pi/180 ;

% Free-stream conditions
freeStream.rho   = 1.225 ;          % free-stream density
freeStream.v     = 1.0 ;            %             velocity absolute value
freeStream.alpha = 0.0 * deg2rad ;  %             velocity direction

% Airfoil definition
airfoil(1).airfoil_str  = 'NACA0012' ;
airfoil(1).chord        = 1.0 ;
airfoil(1).theta        = 0.0 * deg2rad ;
airfoil(1).xcRefPoint   = 0.25 ;
airfoil(1).refPoint     = [ 0.0 ; 0.0 ] ;
airfoil(1).nChordPanels = 40 ;

airfoil(2).airfoil_str  = 'NACA0012' ;
airfoil(2).chord        = 0.5 ;
airfoil(2).theta        = 7.0 * deg2rad ;
airfoil(2).xcRefPoint   = 0.25 ;
airfoil(2).refPoint     = [ 1.00 ; -0.1 ] ;
airfoil(2).nChordPanels = 40 ;

airfoil(3).airfoil_str  = 'NACA0012' ;
airfoil(3).chord        = 0.5 ;
airfoil(3).theta        = 0.0 * deg2rad ;
airfoil(3).xcRefPoint   = 0.25 ;
airfoil(3).refPoint     = [ 1.0 ; 1.0 ] ;
airfoil(3).nChordPanels = 30 ;

nAirfoils = length(airfoil) ;
nelems_te = nAirfoils ;

nelems = 0 ; 
npoints = 0 ;
rr = [] ; ee = [] ;
ee_te = [] ;

for ia = 1 : nAirfoils

  % Check that the first 4 digits are 'NACA'
  if ( strcmp(airfoil(ia).airfoil_str(1:4),'NACA') == 0 )
    error([' Only NACA airfoils are implemented. ', ...
            'The airfoil string must beging with NACA. STOP.' ])
  end
  
  if (     length(airfoil(ia).airfoil_str) == 8 ) % NACA 4-digit airfoil
    M = str2num(airfoil(ia).airfoil_str(5)  ) ;  
    P = str2num(airfoil(ia).airfoil_str(6)  ) ;
    SS= str2num(airfoil(ia).airfoil_str(7:8)) ;
    mirror = 0 ;
    [ x1 , y1 ]  = setAirfoil4( M , P , SS , ...
           airfoil(ia).chord       , airfoil(ia).nChordPanels , ...
           airfoil(ia).refPoint(1) , airfoil(ia).refPoint(2)  , ...
           airfoil(ia).xcRefPoint  , airfoil(ia).theta , mirror ) ;
                                   
  elseif ( length(airfoil(ia).airfoil_str) == 9 ) % NACA 5-digit airfoil
    error([' NACA 5-digit airfoils not implemented yet ' ] ) 
  else
    error([' Only NACA 4 and 5-digit airfoils are implemented. ', ...
            'airfoil(',ia,').airfoil_str must be 8 or 9-character string'])
  end

  rr = [ rr ; [ x1' , y1' ] ] ;
  ee = [ ee ; ...
         npoints + [ (1:length(x1)-1)' , (2:length(x1))' ] ]  ;

  ee_te = [ ee_te ; ...
            [ 1 , 2*airfoil(ia).nChordPanels ] + nelems ] ;

  npoints = npoints + length(x1) ; 
  nelems = nelems + 2 * airfoil(ia).nChordPanels ;

end

% Geometrical quantities of the elements
len = sqrt( sum( ( rr(ee(:,1),:) - rr(ee(:,2),:) ).^2 , 2 ) ) ;
sintheta = ( rr(ee(:,2),2) - rr(ee(:,1),2) )./len;
costheta = ( rr(ee(:,2),1) - rr(ee(:,1),1) )./len;
theta    = atan2(sintheta,costheta);
nvers = [-sintheta; costheta];
tvers = [ costheta; sintheta];
% Centri dei pannelli
rrc = ( rr(ee(:,1),:) + rr(ee(:,2),:) ) ./ 2.0 ;

% Build the linear system
A  = zeros(nelems+nelems_te   ) ; 
b  = zeros(nelems+nelems_te,1 ) ;
Au = zeros(nelems,nelems+nelems_te) ; Av = Au ;

usij     = zeros(nelems);
vsij     = zeros(nelems);
uvij     = zeros(nelems);
vvij     = zeros(nelems);
uv1i     = zeros(nelems,1);
vv1i     = zeros(nelems,1);

for ii = 1 : nelems

    % sources ------------- 
    for jj = 1 : nelems 
        x1 = rrc(ii,1) - rr(ee(jj,1),1);    y1 = rrc(ii,2) - rr(ee(jj,1),2);    r1 = [x1;y1];
        x2 = rrc(ii,1) - rr(ee(jj,2),1);    y2 = rrc(ii,2) - rr(ee(jj,2),2);    r2 = [x2;y2];
        
        if ii == jj
            sij = 1;
            bij = pi;
        else
            
            sij    = norm(r2)/norm(r1) ;
            sinbij = (x1*y2 - x2*y1) / norm(r1) / norm(r2) ;
            cosbij = (x1*x2 + y1*y2) / norm(r1) / norm(r2) ;
            bij    = atan2(sinbij,cosbij);
        end
        
        usijstar = -log(sij)/(2*pi) ;
        vsijstar =      bij /(2*pi) ;
        
        uvijstar =     -bij /(2*pi) ;
        vvijstar = -log(sij)/(2*pi) ;
        
        usij(ii,jj) = usijstar*costheta(jj) - vsijstar*sintheta(jj) ;
        vsij(ii,jj) = usijstar*sintheta(jj) + vsijstar*costheta(jj) ;
                                                                   
        uvij(ii,jj) = uvijstar*costheta(jj) - vvijstar*sintheta(jj) ;
        vvij(ii,jj) = uvijstar*sintheta(jj) + vvijstar*costheta(jj) ;
        
        A (ii,jj)   = -sintheta(ii)*usij(ii,jj) + costheta(ii)*vsij(ii,jj);
        Au(ii,jj)   = usij(ii,jj);
        Av(ii,jj)   = vsij(ii,jj);
    end

    % vortices ---------------- 
    nelems_tmp = 0 ;
    for ia = 1 : nAirfoils
      uv1i(ii) = sum(uvij(ii,nelems_tmp+1:nelems_tmp+2*airfoil(ia).nChordPanels));
      vv1i(ii) = sum(vvij(ii,nelems_tmp+1:nelems_tmp+2*airfoil(ia).nChordPanels));

      A( ii,nelems+ia) = -sintheta(ii)* uv1i(ii) + costheta(ii)* vv1i(ii) ;
      Au(ii,nelems+ia) = uv1i(ii) ;
      Au(ii,nelems+ia) = vv1i(ii) ;
      nelems_tmp = nelems_tmp + airfoil(ia).nChordPanels * 2 ;

    end

    b(ii)       =  sintheta(ii)* freeStream.v * cos(freeStream.alpha) ...
                 - costheta(ii)* freeStream.v * sin(freeStream.alpha)     ;

end

% Kutta conditions ------
nelems_tmp = 0 ;
for ia = 1 : nAirfoils
  A(nelems+ia, 1:nelems ) = costheta(ee_te(ia,1)) .* usij(ee_te(ia,1),:) + ...
                            sintheta(ee_te(ia,1)) .* vsij(ee_te(ia,1),:) + ...
                            costheta(ee_te(ia,2)) .* usij(ee_te(ia,2),:) + ...
                            sintheta(ee_te(ia,2)) .* vsij(ee_te(ia,2),:)       ;
  for ja = 1 : nAirfoils
    A(nelems+ia,nelems+ja) = costheta(ee_te(ia,1)) * uv1i(ee_te(ja,1)) + ... 
                             sintheta(ee_te(ia,1)) * vv1i(ee_te(ja,1)) + ... 
                             costheta(ee_te(ia,2)) * uv1i(ee_te(ja,2)) + ... 
                             sintheta(ee_te(ia,2)) * vv1i(ee_te(ja,2))       ;
  end

  b(nelems+ia) = - ( costheta(ee_te(ia,1)) .* freeStream.v * cos(freeStream.alpha) + ...
                     sintheta(ee_te(ia,1)) .* freeStream.v * sin(freeStream.alpha) + ...
                     costheta(ee_te(ia,2)) .* freeStream.v * cos(freeStream.alpha) + ...
                     sintheta(ee_te(ia,2)) .* freeStream.v * sin(freeStream.alpha)       ) ;

end


vv = A\b;

u = Au * vv + freeStream.v*cos(freeStream.alpha);
v = Av * vv + freeStream.v*sin(freeStream.alpha);
modu = sqrt(u.^2 + v.^2);
angolo = atan2(v,u);

% uvct = [((1:4*n))' ,u, v, modu, angolo, theta'];
% 
% ang = [sintheta' costheta' theta'.*180/pi];

vTi =  u .* costheta + v .* sintheta ;
vNi = -u .* sintheta + v .* costheta ; 

cP = 1 - vTi.^2 ./freeStream.v^2;

figure ; hold on
plot(rr(:,1),rr(:,2),'-o') , axis equal 
plot(rrc(:,1),cP,'-x')
hold off
