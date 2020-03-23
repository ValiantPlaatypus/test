close all ; clear all ; clc


data = load('box1000.dat') ;

x_tol = 0.01 ;
x_out = 3.4 - x_tol ;
y0 = 0.0 ;
sigma = 0.1 ;
du = 0.2 ;

figure ; hold on
quiver(data(:,1),data(:,2),data(:,3),data(:,4))

for i1 = 1 : size(data,1)

  if ( data(i1,1) > x_out ) 

    data(i1,3:4) = data(i1,3:4) * ( 1.0 - du*exp(-(data(i1,2)-y0).^2 ./ sigma^2) ) ;

  end

end

quiver(data(:,1),data(:,2),data(:,3),data(:,4))
grid on

save('box1000-drag.dat','data','-ascii')
