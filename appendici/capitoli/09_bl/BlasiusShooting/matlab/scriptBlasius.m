% === Shooting Method: equazione di Blasius dello strato limite ==========

close all
clear all
clc

% -- Descrizione del metodo ----------------------------------------------
% Invece di risolvere l'equazione di Blasius nell'intervallo [0,inf) con 
% le condizioni al contorno g(0) = 0 , g'(0) = 0 , g'(inf) = 1 ...
% 
% 1. si tronca l'intervallo: [0,inf) viene sostituito da [0,xinfty]
% 2. si sostituisce la condizione g'(inf) = 1 con la condizione g''(0) = a:
%    il problema diventa un problema ai valori iniziali, che può essere 
%    scritto nella forma: dy/dx = f(y(x)) , y(0) = y0.
%    L'equazione del primo ordine è scritta per il vettore y, definito 
%    come y = [g ; g' ; g'']
%    Si integra numericamente (con il metodo che preferite; nel programma
%    viene usata la ode45 di matlab) per ottenere y(x) per x in [0,xinfty].
%    La soluzione ottenuta dipende dal valore di g''(0) = a
% 3. si costruisce un metodo numerico "di Newton" per determinare il valore
%    di a.
%    La condizione che si vuole imporre è un'approssimazione di g'(inf) = 1
%    In particolare la condizione è g'(xinfty) = dginfty = 1
%    g' è la seconda componente di y; chiamiamola y2; essa dipende da a
%    
%    0 = y2(xinfty; a) - dginfty = F(a)
%
%    Il metodo di Newton per trovare a t.c F(a) = 0, si scrive come:
%
%    0 = F(a+da) = F(a) + [dF/da (a)] * da  ---> da = -[dF/da (a)]^-1 * F(a)
%    aggiornamento valore di a: a = a + da
%
%    Si itera fino a quando il modulo del residuo F(a) non è minore della
%    tolleranza o si supera un numero di iterazioni massimo.
%
% == OSS: come caloclare dF/da ????? ======================================
% Il calcolo di dF/da avviene per via numerica, poichè non è nota un'espres-
% sione analitica. Viene calcolata con differenze finite centrate.
% dF / da = ( F(xinfty; a+h) -  F(xinfty; a-h)) / (2*h) = 
%         = (y2(xinfty; a+h) - y2(xinfty; a-h)) / (2*h)
% dove si è indicata con F(xinfty,a+h) = yx(xinfty; a+h) - dginfty. Inoltre
%   . la dginfty vale 1 e quindi la sua derivata è nulla
%   . y2(xinfty,a+h) è il valore della derivata g' calcolata in xinfty, con 
%     condizioni iniziali g(0)=0;g'(0)=0;g''(0)=a+h
% =========================================================================

% -- Definizione dell'intervallo -----------------------------------------
%   L'intervallo [0 , inf) viene troncato; l'equazione viene risolta
%   numericamente sull'intervallo [0,xinfty]
x0 = 0.0;
xinfty = 30.0;

% -- Valore g'(inf) ------------------------------------------------------
dginfty = 1.0;

% -- Parametri del metodo numerico ---------------------------------------
tol = 1e-9;     % tolleranza per il metodo di Newton
h   = 1e-3;     % viene usata per calcolare le derivate numeriche 
a   = 0.0;      % guess iniziale per il valore di g''(0)

maxIter = 100;  % numero massimo di iterazioni

nIter = 0;      % inizializzazione numero di iterazioni e res > tol
res = tol + 1;

fprintf('nIter    g"(0)    res \n')

% -- Metodo di Newton ----------------------------------------------------
while ( nIter < maxIter && res > tol )
    
    % Calcolo di y(x;a)
    y0 = [0.0 ; 0.0 ; a];
    [x,y] = ode45(@fBlasius,[x0 xinfty],y0);
    
    % Calcolo di y(x;a+h) e y(x;a-h) da utilizzare nella derivata numerica
    y0_1 = [0.0 ; 0.0 ; a-h];
    [x_1,y_1] = ode45(@fBlasius,[x0 xinfty],y0_1);
    
    y01  = [0.0 ; 0.0 ; a+h];
    [x1,y1] = ode45(@fBlasius,[x0 xinfty],y01);
    
    % 0 = F(a+da) = F(a)+F'(a)*da ---> da = -F(a)/F'(a)
    F  = y(end,2) - dginfty;                 % F(a)
    dF = (y1(end,2) - y_1(end,2))./(2*h);    % F'(a) con diff.finite centrate

    % Calcolo residuo e iterazioni
    res = abs(F);
    nIter = nIter + 1;
    fprintf('  %i     %6.4f  %6.3e \n', nIter, a, res)

    % incremento da e aggiornamento di a
    da = - F / dF ;
    a  = a + da;
    
end

figure(1)
plot(y,x)
legend('g(\eta)','g''(\eta)','g"(\eta)')
xlabel('g,g'',g"')
ylabel('\eta = y/\delta(x)')


% approx d retta asintotica a g(\eta) per \eta \to \infty
m = (y(end,1) - y(end-1,1)) / (x(end) - x(end-1));
yr = y(end,1) + m * (x - x(end));

figure(2)
hFig = figure(2);
set(hFig, 'Position', [300.0 300.0 900 300])
subplot(1,3,1),plot(y(:,1),x,yr,x),xlabel('g(\eta)')  ,ylabel('\eta = y/\delta(x)'),axis equal,axis([-2 8 0 10])
subplot(1,3,2),plot(y(:,2),x)     ,xlabel('g''(\eta)'),ylabel('\eta = y/\delta(x)'),axis square%,axis equal,axis([-1 3 0 4])
subplot(1,3,3),plot(y(:,3),x)     ,xlabel('g"(\eta)') ,ylabel('\eta = y/\delta(x)'),axis square%,axis equal,axis([-1 3 0 4])
