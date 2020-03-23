[ A , b , Au , Av ] = build_linsys_test( freeStream , elems , ii_te )

nelems = length(elems) ;

% Inizializzazione
A = zeros(nelems+1,nelems+1) ;
b = zeros(nelems+1,1) ;

Au = zeros(nelemes,nelems+1) ;
Av = zeros(nelemes,nelems+1) ;

% Costruzione del sistema
% condizioni al contorno
for i = 1 : nelems

  for j = 1 : nelems
  
    

      A(i,j) =  ############
      A(i,nelems+1) = A(i,nelems+1) + ##########
  
  end

  b(i) = ########### ;

end

% condizione di kutta













