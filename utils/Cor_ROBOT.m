%==========================================================================
% Convert a transformation matrix into the translation and
% the rotations coordinates expressed in the EQUIPMENT ROBOT format
%
% INPUT:
%		T: the 4x4 transformation matrix
%
% OUTPUT:
%		Coord: vector with the translation and rotation
%				[X Y Z rot roll pitch]
%               angles in degree
%                   rot [0, 360]
%                   roll [-90, 90] (more correctly [0,90] and [270, 360])
%                   pitch [0, 360]
%                    In the undeterminate case, roll = +-90, we define
%                    pitch =0
%               distance in cm
%	NB: the parameters of the ROBOToffset function must be correstly configured for this function
%		to work correctly
%
% Author: rla
%	Date: 13/3/06
% MODIFICATION
%   25/2/09: Deal with undeterminate case cos(Y) = 0 (see MID 19674)
%   6/11/09: Use atan2 instead of getangle and use asin(roll) instead of
%               acos(roll) and directly extract translation from matrix to improve precision
%=====================================================================================================

function coord = Cor_ROBOT(T)

offset = ROBOToffset;

%Get the angles
%===============
sinY = -T(3,1);
roll = asin(sinY);

%We explicitely state that the cosY is POSITIVE. Indeed, by convention, we
%define -pi/2 < Ry < pi/2 (see MID 19674)
%cosY = sqrt(1-T(3,1).^2);
%roll = getangle(cosY,sinY);

if(sinY == 1)
    %We are in the undetermined case (see MID 19674)

    %Define the pitch equal to zero
    cosX = 1;
    sinX = 0;
    pitch =0;
    
    cosZ = T(2,2);
    sinZ = -T(1,2);
    rot = atan2(sinZ,cosZ);
    
else
    %We are in the general case
    pitch= atan2(T(3,2) , T(3,3));

    rot = atan2(T(2,1),T(1,1));
end

%Get the translations
%====================
vec = [T(1,4) , T(2,4) , T(3,4)];

%Prepare final vector
%====================
coord = [vec(1) vec(2) vec(3) rot*180/pi roll*180/pi pitch*180/pi];

return