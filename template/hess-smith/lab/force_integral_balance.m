
% xmin = -2.0; xmax =  3.0; nxp = 30;
% ymin = -2.5; ymax =  2.5; nyp = 30;
% xmin = -4.0; xmax =  5.0; nxp = 40;
% ymin = -4.5; ymax =  4.5; nyp = 40;
% xmin = -1.0; xmax =  2.0; nxp = 30;
% ymin = -1.5; ymax =  1.5; nyp = 30;
xmin = -9.0; xmax = 10.0; nxp = 60;
ymin = -9.5; ymax =  9.5; nyp = 60;

xv = linspace(xmin,xmax,nxp);
yv = linspace(ymin,ymax,nyp);

%> === Box ===
%> Points
rr_box = [ xv(end:-1:2); ymin*ones(1,nxp-1) ];
rr_box = [ rr_box , ...
       [ xmin*ones(1,nyp-1); yv(1:end-1) ] ];
rr_box = [ rr_box , ...
       [ xv(1:end-1); ymax*ones(1,nxp-1) ] ] ;
rr_box = [ rr_box , ...
       [ xmax*ones(1,nyp-1); yv(end:-1:2) ] ] ;

npoints = 2*(nxp+nyp-2);

%> Connectivity, and normal and tangent unit vectors
ee_box = [ 1:2*(nxp+nyp-2); 2:2*(nxp+nyp-2)+1 ]; ee_box(2,end) = 1;
nelem_box = size(ee_box,2);
elem_len = zeros(1,nelem_box);
elem_tan = zeros(2,nelem_box);
elem_nor = zeros(2,nelem_box);
for k = 1 : nelem_box
  elem_tan(:,k) = rr_box(:,ee_box(2,k)) - rr_box(:,ee_box(1,k));
  elem_len(  k) = norm(elem_tan(:,k));
  elem_tan(:,k) = elem_tan(:,k) / elem_len(k);
  elem_nor(:,k) = [ -elem_tan(2,k); elem_tan(1,k)];
end

%> Evaluate velocity in the nodes
node_vel = zeros(2,npoints);
for ip = 1 : npoints % loop over box points

  for ie = 1 : length(elems) % loop over airfoil elements
    node_vel(:,ip) = node_vel(:,ip) + ...
      sing_intensity(ie) * compute_velocity_source_p( elems(ie), rr_box(:,ip) ) + ...
      sing_intensity(end)* compute_velocity_vortex_p( elems(ie), rr_box(:,ip) );
  end

end
node_vel(1,:) = node_vel(1,:) + freeStream.vvec(1);
node_vel(2,:) = node_vel(2,:) + freeStream.vvec(2);

figure;
quiver(rr_box(1,:),rr_box(2,:),node_vel(1,:),node_vel(2,:), ...
       'LineWidth',2), axis equal, grid on

node_p = - 0.5 * freeStream.rho * ( node_vel(1,:).^2 + node_vel(2,:).^2 );

%> Compute the integral forces w/ integral balance of the box
I_vel = zeros(2,1);
I_p   = zeros(2,1);
for k = 1 : nelem_box
  I_p = I_p - 0.5* ( node_p( ee_box(1,k) ) + node_p( ee_box(2,k) ) ) * ...
        elem_len(k) * elem_nor(:,k);
  I_vel = I_vel - 0.5* ( ...
        node_vel(:, ee_box(1,k))'*elem_nor(:,k) * node_vel(:,ee_box(1,k)) + ...
        node_vel(:, ee_box(2,k))'*elem_nor(:,k) * node_vel(:,ee_box(2,k)) ) * ...
        elem_len(k);
end
I = I_vel + I_p;

cf_p   = I_p   / ( 0.5 * freeStream.rho * freeStream.v^2 * airfoil.chord );
cf_vel = I_vel / ( 0.5 * freeStream.rho * freeStream.v^2 * airfoil.chord );
cf_box = cf_p + cf_vel;
fprintf('Integral balance \n')
fprintf(' cf_box   = %12.6f, %12.6f \n', cf_box(1), cf_box(2))
fprintf(' cf_box_p = %12.6f, %12.6f \n', cf_p  (1), cf_p  (2))
fprintf(' cf_box_v = %12.6f, %12.6f \n', cf_vel(1), cf_vel(2))

