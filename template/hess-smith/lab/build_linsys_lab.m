function [ A , b , Au , Av ] = build_linsys_lab( freeStream , elems , ee_te )

nelems    = length(elems) ;   % n. of panels
n_te      = 1 ;               % n. of te = n. of Kutta conditions to be assigned

% === Initialize matrices ===
%> Linear system
A  = zeros(nelems+n_te  ) ; % nelems+1, nelems+1 
b  = zeros(nelems+n_te,1) ; % nelems+1
%> Auxiliary matrices for computing on-body velocity
Au = zeros(nelems,nelems+n_te) ; Av = Au ;

% === Fill A matrix and RHS b vector ===
%> 1. assign non-penetration b.c. ( u.n = 0 )
% rows = 1:nelems, cols = 1:nelems+1 
for ii = 1 : nelems    % "passive" elem

  for jj = 1 : nelems    % "active" elem

    %> Induced velocity by elems j on the center of the element i
    % is evaluated as v_xxx = compute_veloctiy_xxx( elems(j), elems(i) )
    vs_ij = compute_velocity_source( elems(jj), elems(ii) );
    vv_ij = compute_velocity_vortex( elems(jj), elems(ii) );

    %> A matrix
    A(ii,jj) = elems(ii).nver' * vs_ij;  % OK
    A(ii,nelems+1) = A(ii,nelems+1) + ...
               elems(ii).nver' * vv_ij ;

    %> Auxiliary matrices ( aerodynamic influence coefficients to 
    % evaluate on-body velocity )
    Au(ii,jj) = vs_ij(1);
    Av(ii,jj) = vs_ij(2);
    Au(ii,nelems+1) = Au(ii,nelems+1) + vv_ij(1);
    Av(ii,nelems+1) = Av(ii,nelems+1) + vv_ij(2);

  end

  %> RHS
  b(ii) = - elems(ii).nver' * freeStream.vvec ;

end

% Kutta condition ????
A(end,:) = 0.0; A(end,end) = 1.0;
b(end,:) = 0.0;
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
%     % Velocity induced by unitary intensity singularities of
%     % the j-th elem on the centre of the i-th elem
%     vs(ii,jj).v =  % ##### 
%     vv(ii,jj).v =  % ##### 
% 
%     %> === A matrix of the linear system ===
%     % Fill source to elems matrix elements
%     %  A(i,j) = ni' * vs_j(r_i)
%     A(ii,jj) = % ##### 
%     % Update vortex to elems matrix elements (accumulation)
%     %  A(i,j) += ni . vv_j(r_i)
%     A(ii,nelems+1) = ...
%        A(ii,nelems+1) + % #####
% 
%     %> === Auxiliary matrices ===
%     % Matrices to retrieve on-body velocity
%     Au(ii,jj) =  % ##### ; % x-comp    % sources
%     Av(ii,jj) =  % ##### ; % y-comp
%     Au(ii,nelems+1 ) = Au(ii,nelems+1 ) + % ##### ; % vortices
%     Av(ii,nelems+1 ) = Av(ii,nelems+1 ) + % ##### ;
% 
%   end
% 
%   %> === RHS, vector b === 
%   b(ii) = - elems(ii).nver' * freeStream.vvec ;
% 
% end
%  
% %> Kutta condition
% % HOMEWORK: fill the last row of the matrix A, and the last element
% % of the vector b
% A(end,:) = 0.0; A(end,end) = 1.0;
% b(end)   = 0.0;
% 
% 
% 
% % % === Fill A matrix and RHS b vector ===
% % %> 1. assign non-penetration b.c. ( u.n = 0 )
% % % rows = 1:nelems, cols = 1:nelems+nTe 
% % for ii = 1 : nelems    % "passive" elem
% % 
% %     for jj = 1 : nelems    % "active" elem
% % 
% %       % Velocity induced by unitary intensity singularities of
% %       % the j-th elem on the centre of the i-th elem
% %       vs(ii,jj).v =  % ##### 
% %       vv(ii,jj).v =  % ##### 
% % 
% %       %> === A matrix of the linear system ===
% %       % Fill source to elems matrix elements
% %       %  A(i,j) = ni' * vs_j(r_i)
% %       A(ii,jj) = % ##### 
% %       % Update vortex to elems matrix elements (accumulation)
% %       %  A(i,j) += ni . vv_j(r_i)
% %       A(ii,nelems+1) = ...
% %          A(ii,nelems+1) + % #####
% % 
% %       %> === Auxiliary matrices ===
% %       % Matrices to retrieve on-body velocity
% %       Au(ii,jj) =  % ##### ; % x-comp    % sources
% %       Av(ii,jj) =  % ##### ; % y-comp
% %       Au(ii,nelems+1 ) = ...
% %          Au(ii,nelems+1 ) + % ##### ; % vortices
% %       Av(ii,nelems+1 ) = ...
% %          Av(ii,nelems+1 ) + % ##### ;
% % 
% %     end
% % 
% %     % 
% %     b(ii) = - elems(ii).nver' * freeStream.vvec ;
% % 
% % end
% % 
% % % w/o Kutta condition
% % A(end,:) = 0.0; A(end,end) = 1.0;
% % b(end)   = 0.0;
% 
% 
% % %> 2. assign Kutta condition at all the trailing edges
% % % Kutta condition is approximated as:
% % %
% % %  U_TE_upper . tTE_upper + U_TE_lower . tTE_lower = 0
% % %
% % % rows = nelems+1:nelems+nTe, cols = 1:nelems+nTe 
% % 
% % for ia = 1 : n_te
% % 
% %   for jj = 1  : nelems
% % 
% %   A(nelems+ia, jj ) = elems(ee_te(ia,1)).tver' * vs(ee_te(ia,1),jj).v + ...
% %                       elems(ee_te(ia,2)).tver' * vs(ee_te(ia,2),jj).v ;
% %   A(nelems+ia,nelems+elems(jj).airfoilId) = A(nelems+ia,nelems+elems(jj).airfoilId) + ...
% %                              elems(ee_te(ia,1)).tver' * vv(ee_te(ia,1),jj).v +  ...
% %                              elems(ee_te(ia,2)).tver' * vv(ee_te(ia,2),jj).v ; 
% %   end
% % 
% %   b(nelems+ia) = - ( elems( ee_te(ia,1) ).tver + elems( ee_te(ia,2) ).tver )' * freeStream.vvec ;
% % 
% % end
% 
