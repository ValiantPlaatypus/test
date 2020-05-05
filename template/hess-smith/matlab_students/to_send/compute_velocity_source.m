% compute the velocity induced by a linear source elem_j on the centre of the elem_i

function v = compute_velocity_source(elem_j,elem_i)

% Vectors connecting the nodes of the "active" j-th elem, with 
% the centre of the "passive" i-th elem
r1 = elem_i.cen - elem_j.ver1 ;
r2 = elem_i.cen - elem_j.ver2 ;

if ( elem_i.id == elem_j.id ) % self induction
  sij = 1.0 ;
  bij = pi ;
else % 
  % ...
  sij = % ...
  bij = % ...
end

% Velocity in the local reference frame of the j-th elem
vstar = [ -log(sij)/(2*pi) ; bij /(2*pi) ] ;

% Velocity in the global reference frame. Rotation
v = [ elem_j.tver , elem_j.nver ] * vstar ; 

