clear all ; close all ;

% Integral loads and their uncertainty ---------------- %
%                                                       % 
% - computed from momentum integral equation            %
% - RSS for the uncertainty propagation                 %
%                                                       % 
% ----------------------------------------------------- %
% assume u ~= piecewise linear function

% physical quantities ---------------------------------
rho = 1.225 ; 

% load nodal variables (x,y,u,v) ----------------------
data = load('box1000-drag.dat') ;
x = data(end:-1:1,1) ; % reverse order for anti-clockwise direction
y = data(end:-1:1,2) ;  
u = data(end:-1:1,3) ;  
v = data(end:-1:1,4) ;
nNodes = length(x) ;
nElems = nNodes    ;

% connectivity ----------------------------------------
rr = [ x' ; y' ] ; % nodes coordinates
uu = [ u' ; v' ] ; % velocity
ee = [ 1:nNodes ; [2:nNodes,1] ] ; % node-elem connectivity
sur= sqrt( (rr(1,ee(2,:))-rr(1,ee(1,:))).^2 + ...
           (rr(2,ee(2,:))-rr(2,ee(1,:))).^2 ) ; % elem surface (length in 2d)
nor= [ (  rr(2,ee(2,:))-rr(2,ee(1,:)) ) ./ sur ; ...
       ( -rr(1,ee(2,:))+rr(1,ee(1,:)) ) ./ sur  ]  ; % elem unit normal vector

% integral loads from momentum equation ---------------
% R = -\oint_S {rho u u.n} - \oint_S {p n},
% using Bernoulli equation for p,
%   p = p_inf + 1/2 rho V_inf^2 - 1/2 rho |u|^2
%
% Incremental computation of the aerodynamic force ----
R = [ 0 ; 0 ] ;
mref = [ 1/3 , 1/6 ; 1/6 , 1/3 ] ;
for is = 1 : nElems

  R1 = [ 0 ; 0 ] ;

  for i1 = 1 : 2
    for i2 = 1 : 2
    
      R1 = R1 + mref(i1,i2) * (-uu(:,ee(i1,is)) * uu(:,ee(i2,is))' * nor(:,is) ...
                          +0.5* uu(:,ee(i1,is))'* uu(:,ee(i2,is))  * nor(:,is) ) ;
    
    end
  end

  R1 = R1 * sur(is) ;
  R = R + R1 ;

end
R = R*rho ;
fprintf('\n Rx : %f   ',R(1))
fprintf(  ' Ry : %f \n',R(2))

% Uncertainty -----------------------------------------
sigma_u  = 0.05 ;
sigma_v  = sigma_u  ;
sigma_uu = [ sigma_u *ones(1,nNodes) ; sigma_v *ones(1,nNodes) ] ;

% mass matrix (sparse format)
irow = [ ee(1,:) , ee(1,:) , ee(2,:) , ee(2,:) ] ;
icol = [ ee(1,:) , ee(2,:) , ee(1,:) , ee(2,:) ] ;
emas = [ 1/3*sur, 1/6*sur, 1/6*sur, 1/3*sur ] ; 
Mass = sparse(irow,icol,emas) ;
% SensL = [ zeros(1,nNodes) ; -rho*(nor(1,:).*uu(2,:) - nor(2,:).*uu(1,:) ) ] ;
% SensD = [ rho*(nor(1,:).*uu(2,:) - nor(2,:).*uu(1,:) ) ; zeros(1,nNodes)  ] ;
SensL = [-rho*(nor(1,:).*uu(1,:) + nor(2,:).*uu(2,:) ) ; -rho*(nor(1,:).*uu(2,:) - nor(2,:).*uu(1,:) ) ] ;
SensD = [ rho*(nor(1,:).*uu(2,:) - nor(2,:).*uu(1,:) ) ; -rho*(nor(1,:).*uu(1,:) + nor(2,:).*uu(2,:) ) ] ;

sigma2_l = SensL(1,:) * Mass * diag(sigma_uu(1,:).^2) * Mass' * SensL(1,:)' + ...
           SensL(2,:) * Mass * diag(sigma_uu(2,:).^2) * Mass' * SensL(2,:)' ;
sigma2_d = SensD(1,:) * Mass * diag(sigma_uu(1,:).^2) * Mass' * SensD(1,:)' + ...
           SensD(2,:) * Mass * diag(sigma_uu(2,:).^2) * Mass' * SensD(2,:)' ;
fprintf('\n n.nodes  : %d ', nNodes)
fprintf(  ' sigma2_l : %f '    ,sigma2_l)
fprintf(  ' sigma2_d : %f \n\n',sigma2_d)

% % Incremental computation of sigma_l2
% sigma2_l = 0 ;
% Sens_l = -rho*[ nn(1,:).*uu(2,:) - nn(2,:).*uu(1,:) ] ;
% Sens_d =  rho*[ nn(1,:).*uu(2,:) - nn(2,:).*uu(1,:) ] ;
% for is = 1 : nElems
% 
%   sigma1 = 0 ;
% 
%   for i2 = 1 : 2
%     for i1 = 1 : 2
%       
%       sigma1 = sigma1 + ( sur(is) .* ( ...
%              Sens_l(ee(i1,is)) * mref(i1,i2)  ) ) ; %* sigma2_uu(2,ee(i2,is))  ) ; 
%     
%     end
% 
%     sigma2_l = sigma2_l + sigma1 ;
% 
%   end
% 
% 
% end



figure ; hold on
plot(rr(1,:),rr(2,:),'-o') 
plot(rr(1,1),rr(2,1),'o','MarkerFaceColor','red')
quiver(rr(1,:),rr(2,:),uu(1,:),uu(2,:))
axis equal , grid on


% convergenza
nnodes = [40,80,120,400,1000] ; 
sigma2 = [5.7179,2.9431,1.9808,0.6021,0.2416] ;
figure
loglog(nnodes,sigma2) , grid on , axis([10 1000 0.1 10.0])
