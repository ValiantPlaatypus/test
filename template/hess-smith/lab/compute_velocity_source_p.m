% compute the velocity induced by a linear source elem_j on the centre of the elem_i

function v = compute_velocity_source_p(elem_j, r_i)

r1 = r_i - elem_j.ver1 ;
r2 = r_i - elem_j.ver2 ;

sij = norm(r2)/norm(r1) ;
vcross = cross([r1;0.0],[r2;0.0]) ;
sinbij = vcross(3) / norm(r1) / norm(r2) ;
cosbij = r1' * r2 / norm(r1) / norm(r2) ;
bij = atan2(sinbij,cosbij) ;

vstar = [ -log(sij)/(2*pi) ; bij/(2*pi) ] ;

v = [ elem_j.tver , elem_j.nver ] * vstar ; 

