% TODO:
%  complete the function, filling the lines indicated with <<<

% compute the velocity induced by a unitary linear source elem_j on the point rr_i
% Inputs :
% - elem_j: structure containing the all the relevant info of the line elems,
%     as described in build_geometry.m
% - rr_i(2,1): point where the velocity is computed 
% Outputs:
% - v(2,1): velocity induced in the point rr_i by the source distribution on elem_j
% 

function v = compute_velocity_source( elem_j , rr_i )

tol = 1.0e-6 ; % tolerance to avoid singularities

% Vectors connecting the vertices of elem_j to the point rr_i
r1 = % <<<
r2 = % <<<

% First, the velocity is computed in the normal and tangential directions w.r.t.
% the linear element elem_j; in the reference frame {elem_j.tver,elem_j.nver},
% the expression of the velocity reads:
% v = -1/(2*pi)*ln(sij) * elems_j.nver + ...
%     +1/(2*pi)*bij     * elems_j.tver , 
% where:
% - sij is the ratio |r2|/|r1|
sij = % <<<

% - bij is the angle [rad!] between the vectors r1 and r2
%   ( regularisation is needed, as follows )
if ( norm(rr_i - elem_j.cen) >= tol )   % "regular" case: formula of the exercise
 
  % cross and dot product of r1,r2 to get sin(bij) and cos(bij) ...
  vcross = cross([r1;0.0],[r2;0.0]) ;
  sinbij = vcross(3) / norm(r1) / norm(r2) ;
  cosbij = r1' * r2 / norm(r1) / norm(r2) ;

  % ... and thus bij as atan2(sin(bij),cos(bij))
  bij = atan2(sinbij,cosbij) ;

else                                   % limiting cases for rr_i ---> elem_j.cen
  if ( ( rr_i - elem_j.cen )' * elem_j.nver >= 0.0 ) ;  bij = pi ;
  else  bij = -pi ; % (rr_i - elem_j.cen)'*elem_j.nor < 0.0 );  
  end
end

% v = vt * elem_j.tver + vn * elem_j.nver ;
vt = % <<<
vn = % <<<
vstar = [ vt , vn ] ;

% Then, rotation to get the components of the velocity in the global ref.frame from
% its components in the tangential and normal directions w.r.t. the linear elem
% elem_j
% [ vx ] = [  Rot_mat  ][ vt ] 
% [ vy ]   [           ][ vn ] 

% v = vx * \hat{x} + vy * \hat{y} ;
v = [ elem_j.tver , elem_j.nver ] * vstar ; 

