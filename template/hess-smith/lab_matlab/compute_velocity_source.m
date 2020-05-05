% compute the velocity induced by a linear source elem_j on the centre of the elem_i

function v = compute_velocity_source(elem_j,elem_i)

r1 = elem_i.cen - elem_j.ver1 ;
r2 = elem_i.cen - elem_j.ver2 ;

if ( elem_i.id == elem_j.id )
  sij = 1.0 ;
  bij = pi ;
else
  sij = norm(r2)/norm(r1) ;
  vcross = cross([r1;0.0],[r2;0.0]) ;
  sinbij = vcross(3) / norm(r1) / norm(r2) ;
  cosbij = r1' * r2 / norm(r1) / norm(r2) ;
  bij = atan2(sinbij,cosbij) ;
end

vstar = [ -log(sij)/(2*pi) ; bij /(2*pi) ] ;

v = [ elem_j.tver , elem_j.nver ] * vstar ; 

