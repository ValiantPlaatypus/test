% Update geometry appending new airfoils to ee,rr,ee_te arrays and incrementing nelems,npoints

function [ ee , rr , ee_te , elems , nelems , npoints ] = ...
            build_geometry( eei , rri , ee_tei , elemsi , nelemsi , npointsi , airfoil )

% Initialise output equal to the input (to be updated within this function)
nelems  = nelemsi ;
npoints = npointsi ;

% === Read NACA airfoil and build the geometry ===
if ( strcmp(airfoil.airfoil_str(1:4),'NACA') == 0 )
  % Check that the first 4 digits are 'NACA'
  error([' Only NACA airfoils are implemented. ', ...
          'The airfoil string must beging with NACA. STOP.' ])
end

if (     length(airfoil.airfoil_str) == 8 ) % NACA 4-digit airfoil
  M = str2num(airfoil.airfoil_str(5)  ) ;  
  P = str2num(airfoil.airfoil_str(6)  ) ;
  SS= str2num(airfoil.airfoil_str(7:8)) ;
  mirror = 0 ;
  [ x1 , y1 ]  = setAirfoil4( M , P , SS , ...
         airfoil.chord       , airfoil.nChordPanels , ...
         airfoil.refPoint(1) , airfoil.refPoint(2)  , ...
         airfoil.xcRefPoint  , airfoil.theta , mirror ) ;
                                 
elseif ( length(airfoil.airfoil_str) == 9 ) % NACA 5-digit airfoil
  error([' NACA 5-digit airfoils not implemented yet ' ] ) 
else
  error([' Only NACA 4 and 5-digit airfoils are implemented. ', ...
          'airfoil(',ia,').airfoil_str must be 8 or 9-character string'])
end

% === Append new points and elements to node and connectivity arrays ===
rr = [ rri ; [ x1' , y1' ] ] ;
ee = [ eei ; ...
       npointsi + [ (1:length(x1)-1)' , (2:length(x1))' ] ]  ;
ee_te = [ ee_tei ; ...
          [ 1 , 2*airfoil.nChordPanels ] + nelemsi ] ;

% === Add new elems to the elems structure ===
% elems is an array of "elements type", objects containing all the relevant info
for ie = nelems + 1 : nelems + 2*airfoil.nChordPanels
   elemsi(ie).airfoilId = airfoil.id ;
   elemsi(ie).id   = ie ;
   elemsi(ie).ver1 = rr(ee(ie,1),:)' ;
   elemsi(ie).ver2 = rr(ee(ie,2),:)' ;
   elemsi(ie).cen  = 0.5 * ( elemsi(ie).ver1 + elemsi(ie).ver2 ) ;
   elemsi(ie).len  =   norm( elemsi(ie).ver2 - elemsi(ie).ver1 ) ;
   elemsi(ie).tver = (  elemsi(ie).ver2 - elemsi(ie).ver1 ) / elemsi(ie).len ;
   elemsi(ie).nver = [ -elemsi(ie).tver(2) ; elemsi(ie).tver(1) ] ;
end

elems = elemsi ; % Update elems structure 

% === Update number of points and elems ===
npoints = npoints + length(x1) ; 
nelems = nelems + 2 * airfoil.nChordPanels ;

