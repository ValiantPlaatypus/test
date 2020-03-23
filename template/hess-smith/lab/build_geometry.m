% TODO:
% 1. build the required geometrical quantities
% 2. uncomment last lines to build the elems structure

% Input:
% - airfoil (structure)
%   .id            (integer): id. number of the airfoil (almost useless for just one airfoil)
%   .airfoil_str   (string) : identifying the airfoil, e.g 'NACA0012'
%   .chord         (real)   : chord
%   .theta         (real)   : angle of attack [rad]
%   .xcRefPoint    (real)   : position along the chord (in %) of the reference point, e.g. 0.25
%   .refPoint(2,1) (real)   : vector of size (2,1) containing x0,y0 coord. of the refPoint in 
%                             in the "global" reference system 
%   .nChordPanels  (integer): n. of surface elements for the discretisation of the upper and lower
%                             surface of the airfoil: the total n. of surface elements will be 2* 
% Output:
% - npoints        (integer): n. of points on the surface of the airfoil, used for discretisation
% - nelems         (integer): n. of "surface" line elems of the airfoil, used for discretisation
% - rr(2,npoints)  (real)   : array containing the x,y coords of the points.
%                             rr(1,i), rr(2,i) contain the x,y coord of the i-th point
% - ee(2,nelems)   (integer): connectivity matrix.
%                             ee(1,ie), ee(2,i2) contain the id of the first and second node of the 
%                             elems
% - ii_te(2,1)     (integer): vector containing the id of the elems at the trailing edge
% - elems(nelems) (structure): structure containing all the relevant informations of the elements.
%                              see last lines of the function.
%                              Uncomment and fill these lines, after having computed all the other
%                              outputs 

function [ rr , ee , ii_te , elems , nelems , npoints ] = build_geometry( airfoil )

% === Read NACA airfoil and build the geometry ===
if ( strcmp(airfoil.airfoil_str(1:4),'NACA') == 0 )
  % Check that the first 4 digits are 'NACA'
  error([' Only NACA airfoils are implemented. ', ...
          'The airfoil string must begin with NACA. STOP.' ])
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

npoints = length(x1) ;
nelems  = npoints-1 ;

% === Define the arrays of points, node-to-elem connectivity and id of the   ===
% === elems at the TE (used for Kutta condition)                             ===
rr = [  x1 ; y1 ] ;
ee = [ (1:npoints-1) ; (2:npoints) ]  ;
ii_te = [ 1 , nelems ] ;

% === Create elems structure ===
% elems is an array of "elements type", objects containing all the relevant info
for ie = 1 : nelems 
   %> id. of the airfoil (quite useless with only one airfoil)
   elems(ie).airfoilId = airfoil.id ;
   %> id. number of the elem
   elems(ie).id   = ie ;                                      
   %> coordinates of the 1st node of the ie-th elem 
   elems(ie).ver1 = rr(:,ee(1,ie))  ;                         
   %> coordinates of the 2nd node of the ie-th elem 
   elems(ie).ver2 = rr(:,ee(2,ie))  ;                         
   %> coordinates of the center of the ie-th elem 
   elems(ie).cen  = 0.5 * ( elems(ie).ver1 + elems(ie).ver2 ) ;                
   %> length of the ie-th elem
   elems(ie).len  = norm( elems(ie).ver2 - elems(ie).ver1 ) ;
   %> tangent unit vector
   elems(ie).tver = ( elems(ie).ver2 - elems(ie).ver1 ) / elems(ie).len ;                         
   % outward pointing normal unit vector
   elems(ie).nver = [ -elems(ie).tver(2) ; elems(ie).tver(1) ] ;                   
end


