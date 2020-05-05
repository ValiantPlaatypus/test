function [ A , b , Au , Av ] = build_linsys( freeStream , elems , ee_te )

nelems    = length(elems) ;   % n. of panels
n_te      = size(ee_te,1) ;   % n. of bodies = n. of te = n. of Kutta conditions to be assigned

% Initialise matrices 
A  = zeros(nelems+n_te   ) ; 
b  = zeros(nelems+n_te,1 ) ;
Au = zeros(nelems,nelems+n_te) ; Av = Au ;

% === Fill A matrix and RHS b vector ===
% 1. assign the non-penetration b.c. ( u.n = 0 ) -------------------------------
% rows = 1:nelems, cols = 1:nelems+nTe 
for ii = 1 : nelems    % elems where velocity is induced

    for jj = 1 : nelems    % inducing elems

      % Compute the elements of the matrix of the linear system ----------------

      % Velocity induced by unitary intensity singularities of the j-th elem on
      %  the centre of the i-th elem
      vs(ii,jj).v = compute_velocity_source( elems(jj) , elems(ii) ) ;  % sources
      vv(ii,jj).v = compute_velocity_vortex( elems(jj) , elems(ii) ) ;  % vortices

      % Fill source to elems matrix elements
      A (ii,jj)   = elems(ii).nver' * vs(ii,jj).v ;
      % Update vortex to elems matrix elements
      A( ii,nelems+elems(jj).airfoilId ) = A( ii,nelems+elems(jj).airfoilId ) + elems(ii).nver' * vv(ii,jj).v ; 

      % Fill matrices to retrieve the velocity field
      Au(ii,jj)   = vs(ii,jj).v(1) ;    % sources
      Av(ii,jj)   = vs(ii,jj).v(2) ;
      Au(ii,nelems+elems(jj).airfoilId ) = Au(ii,nelems+elems(jj).airfoilId ) + vv(ii,jj).v(1) ; % vortices
      Av(ii,nelems+elems(jj).airfoilId ) = Av(ii,nelems+elems(jj).airfoilId ) + vv(ii,jj).v(2) ;

    end

    % 
    b(ii) = - elems(ii).nver' * freeStream.vvec ;

end

% 2. assign Kutta condition at all the trailing edges -------------------------
% Kutta condition is approximated as:
%
%  U_TE_upper . tTE_upper + U_TE_lower . tTE_lower = 0
%
% rows = nelems+1:nelems+nTe, cols = 1:nelems+nTe 

for ia = 1 : n_te

  for jj = 1  : nelems

  A(nelems+ia, jj ) = elems(ee_te(ia,1)).tver' * vs(ee_te(ia,1),jj).v + ...
                      elems(ee_te(ia,2)).tver' * vs(ee_te(ia,2),jj).v ;
  A(nelems+ia,nelems+elems(jj).airfoilId) = A(nelems+ia,nelems+elems(jj).airfoilId) + ...
                             elems(ee_te(ia,1)).tver' * vv(ee_te(ia,1),jj).v +  ...
                             elems(ee_te(ia,2)).tver' * vv(ee_te(ia,2),jj).v ; 
  end

  b(nelems+ia) = - ( elems( ee_te(ia,1) ).tver + elems( ee_te(ia,2) ).tver )' * freeStream.vvec ;

end

