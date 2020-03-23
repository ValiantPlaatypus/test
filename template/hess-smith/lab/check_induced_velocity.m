close all ; clear all ; clc

% -> check induction
x_v = [ -2.9 : 0.2 : 2.9 ] ; nx = length(x_v) ;
y_v = [ -2.9 : 0.2 : 2.9 ] ; ny = length(y_v) ;
elem_test.ver1 = [ -1.0 ; 0.0 ] ; elem_test.ver2 = [  1.0 ; 0.0 ] ;
elem_test.cen = 0.5 * elem_test.ver1 + ...
                      elem_test.ver2 ;
elem_test.nver = [ 0.0 ; 1.0 ] ;
elem_test.tver = [ 1.0 ; 0.0 ] ;
u_sou = [] ; v_sou = [] ; x_plot = [] ;
u_vor = [] ; v_vor = [] ; y_plot = [] ;
for j = 1 : ny
    for i = 1 : nx

        vs = compute_velocity_source( elem_test , [ x_v(i) ; y_v(j) ] ) ;        
        vv = compute_velocity_vortex( elem_test , [ x_v(i) ; y_v(j) ] ) ;        
        
        u_sou = [ u_sou , vs(1) ] ; v_sou = [ v_sou  , vs(2) ] ; 
        u_vor = [ u_vor , vv(1) ] ; v_vor = [ v_vor  , vv(2) ] ;
        x_plot = [ x_plot , x_v(i) ] ;
        y_plot = [ y_plot , y_v(j) ] ;

  end
end

figure
quiver(x_plot,y_plot,u_sou,v_sou,'LineWidth',2) ; grid on ; axis equal ; title('line source')
figure
quiver(x_plot,y_plot,u_vor,v_vor,'LineWidth',2) ; grid on ; axis equal ; title('line vortex')


