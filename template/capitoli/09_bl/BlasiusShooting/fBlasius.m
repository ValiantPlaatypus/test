function dy = fBlasius(t,y)
% Equazione di Blasius:
% 0.5 g * g'' + g''' = 0    + b.c.
%
% y = [ g ; g' ; g'' ]
% 
%
%  d    / g  \      /  g'            \     % dg/dx = g'
% ____ |  g'  | =  |   g''            |    % dg'/dx = g''
%  dx   \ g''/      \  -0.5 * g * g''/     % dg''/dx = -0.5*g*g'' da Blasius

dy = zeros(3,1);
dy(1) = y(2);
dy(2) = y(3);
dy(3) = -0.5 .* y(1) .* y(3);