function [ A , b , Au , Av ] = build_linsys( freeStream , elems , ii_te )

nelems    = length(elems) ;   % n. of panels
n_te      = size(ii_te,1) ;   % n. of bodies = n. of te = n. of Kutta conditions to be assigned

% === Initialise matrices to zero ===
% --- matrix and rhs vector of the linsys, A*x = b ----------------------------- 
A  = zeros(nelems+n_te   ) ; 
b  = zeros(nelems+n_te,1 ) ;
% --- matrices to retrieve velocity vectors at the control points with       ---
% ---  mat*vec multiplication:   u = Au*x , v = Av*x                         ---
Au = zeros(nelems,nelems+n_te) ; Av = Au ;

% === Fill A matrix and RHS b vector ===
% 1. assign the non-penetration b.c. ( u.n = 0 ) -------------------------------
% rows = 1:nelems, cols = 1:nelems+nTe 
%
for ii = 1 : nelems    % elems where velocity is induced

    for jj = 1 : nelems    % inducing elems

      % Velocity induced by unitary intensity singularities of the j-th elem on
      %  the centre of the i-th elem
      % ( save the velocity AIC to "arrays of vectors" vs(:,:).v , vv(:,:).v, so
      %   that the vectors vs(ii,jj).v and vv(ii,jj).v are the velocity induced
      %   by the jj-th source and vortex of unitary intensity on the control pt
      %   of the ii-th elem, elems(ii).cen )
      vs(ii,jj).v = compute_velocity_source( elems(jj) , elems(ii).cen ) ;  % sources
      vv(ii,jj).v = compute_velocity_vortex( elems(jj) , elems(ii).cen ) ;  % vortices

      % === Fill b.c. block of the matrix A with source AIC ===
      A( ii,jj)   = elems(ii).nver' * vs(ii,jj).v ;
      % === Accumulate vortex AIC in the b.c. block of the matrix A === 
      A( ii,nelems+elems(jj).airfoilId ) =  ...
          A( ii,nelems+elems(jj).airfoilId ) + elems(ii).nver' * vv(ii,jj).v ; 


      % === Fill matrices to retrieve the velocity field ===
      % as before, fill sources block, accumulate vortex contributions
      % -> component x of the velocity field
      Au(ii,jj)   = vs(ii,jj).v(1) ;                            % sources
      Au(ii,nelems+elems(jj).airfoilId ) = ...
          Au(ii,nelems+elems(jj).airfoilId ) + vv(ii,jj).v(1) ; % vortices
      % -> component y of the velocity field
      Av(ii,jj)   = vs(ii,jj).v(2) ;                            % sources
      Av(ii,nelems+elems(jj).airfoilId ) = ...
          Av(ii,nelems+elems(jj).airfoilId ) + vv(ii,jj).v(2) ; % vortices


    end

    % === Fill the b.c. block of the rhs vector b ===
    b(ii) = - elems(ii).nver' * freeStream.vvec ;

end

% 2. assign Kutta condition at all the trailing edges -------------------------
% Kutta condition is approximated as:
%
%  U_TE_upper . tTE_upper + U_TE_lower . tTE_lower = 0
%
% rows = nelems+1:nelems+nTe, cols = 1:nelems+nTe 
%
for ia = 1 : n_te  % n_te = 1 here

  % indices of the elems at the te
  i_te_1 = ii_te(ia,1) ;
  i_te_N = ii_te(ia,2) ;

  for jj = 1  : nelems

  % fill the Kutta condition row with the source contributions 
  A(nelems+ia, jj ) = elems(i_te_1).tver' * vs(i_te_1,jj).v + ...
                      elems(i_te_N).tver' * vs(i_te_N,jj).v ;
  % accumulate the contribution of the vortex in the last elem of the Kutta
  % condition row
  A(nelems+ia,nelems+elems(jj).airfoilId) = ...
         A(nelems+ia,nelems+elems(jj).airfoilId) + ...
            elems(i_te_1).tver' * vv(i_te_1,jj).v +  ...
            elems(i_te_N).tver' * vv(i_te_N,jj).v ; 
  end

  % fill the last component of the rhs ( -(t1+tN).U_inf )
  b(nelems+ia) = - ( elems( i_te_1 ).tver + elems( i_te_N ).tver )' * ...
                     freeStream.vvec ;

end

