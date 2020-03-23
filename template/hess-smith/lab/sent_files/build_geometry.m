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
% ...
% ...
% ...
% ...
% ...
% ...
% ...


% % === Create elems structure ===
% % elems is an array of "elements type", objects containing all the relevant info
% for ie = 1 : nelems 
%    elems(ie).airfoilId = airfoil.id ;
%    elems(ie).id   = % ie ; % id. number of the elem
%    elems(ie).ver1 = % coordinates of the 1st node of the ie-th elem 
%    elems(ie).ver2 = % coordinates of the 2nd node of the ie-th elem 
%    elems(ie).cen  = % coordinates of the center of the ie-th elem 
%    elems(ie).len  = % length of the ie-th elem
%    elems(ie).tver = % tangent unit vector
%    elems(ie).nver = % outward pointing normal unit vector
% end


