
%> Test element
elem_test.ver1 = [ -0.5; 0.0 ];
elem_test.ver2 = [  0.5; 0.0 ];
elem_test.cen  = 0.5 * ( elem_test.ver1 + elem_test.ver2 ) ;
elem_test.len  =   norm( elem_test.ver2 - elem_test.ver1 ) ;
elem_test.tver = (  elem_test.ver2 - elem_test.ver1 ) / elem_test.len ;
elem_test.nver = [ -elem_test.tver(2) ; elem_test.tver(1) ] ;

xvec1 = [ -1.55: 0.1 : 1.55 ]; nx = length(xvec1);
yvec1 = [ -1.55: 0.1 : 1.55 ]; ny = length(yvec1);
xvec   = zeros(nx*ny,1);  yvec   = zeros(nx*ny,1);
uvec_s = zeros(nx*ny,1);  uvec_v = zeros(nx*ny,1);
vvec_s = zeros(nx*ny,1);  vvec_v = zeros(nx*ny,1);
for iy = 1 : ny
  for ix = 1 : nx
    xvec(ix+(iy-1)*nx) = xvec1(ix);
    yvec(ix+(iy-1)*nx) = yvec1(iy);
    vv = compute_velocity_vortex_p( elem_test, [ xvec1(ix); yvec1(iy) ] );
    vs = compute_velocity_source_p( elem_test, [ xvec1(ix); yvec1(iy) ] );
    uvec_v(ix+(iy-1)*nx) = vv(1); vvec_v(ix+(iy-1)*nx) = vv(2);
    uvec_s(ix+(iy-1)*nx) = vs(1); vvec_s(ix+(iy-1)*nx) = vs(2);
  end
end

figure('Position',[ 50 50 1200 600 ])
subplot(1,2,1), quiver(xvec,...
                       yvec,uvec_s, ...
                            vvec_s,'b','LineWidth',2'   , ...
                                       'MaxHeadSize',0.2), title('Source'), hold on
plot([-0.5,0.5],[0,0],'k','LineWidth',3)
axis equal, hold off
subplot(1,2,2), quiver(xvec,...
                       yvec,uvec_v, ...
                            vvec_v,'b','LineWidth',2'   , ...
                                       'MaxHeadSize',0.2), title('Vortex'), hold on
plot([-0.5,0.5],[0,0],'k','LineWidth',3)
axis equal, hold off
