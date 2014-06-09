%=====================================================
% Convert a transformation matrix into the translation and
% the rotations expressed in the ISOCENTRIC format
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
%               distance in mm
%
% Author: rla
%	Date: 13/3/06
% MODIFICATION
%   25/2/09: Deal with undeterminate case cos(Y) = 0 (see MID 19674)
%===========================================================

function coord = Cor_ISO(T)

%Get the angles
% The angles are the same as for the equipment mode
%===================================================
sinY = -T(3,1);

%We explicitely state that the cosY is POSITIVE. Indeed, by convention, we
%define -pi/2 < Ry < pi/2 (see MID 19674)
%cosY = sqrt(1-T(3,1).^2);
%roll = getangle(cosY,sinY);
roll = asin(sinY);

if(sinY == 1)
    %We are in the undetermined case (see MID 19674)

    %Define the pitch equal to zero
    cosX = 1;
    sinX = 0;
    pitch =0;
    
    %cosZ = T(2,2);
    %sinZ = -T(1,2);
    %rot = getangle(cosZ,sinZ);
	rot = atan2(-T(1,2),T(2,2));
    
else
    %We are in the general case
    %cosX = T(3,3) / cosY;
    %sinX = T(3,2) / cosY;
    %pitch= getangle(cosX,sinX);
	pitch = atan2(T(3,2),T(3,3));

    %cosZ = T(1,1) / cosY;
    %sinZ = T(2,1) / cosY;
    %rot = getangle(cosZ,sinZ);
	rot = atan2(T(2,1),T(1,1));
end

%Get the translations
%====================
IsoFRS = [0 0 0 1]'; %Coordinate of isocente in the FRS

%Get the coordinate of the isocentre expressed in the TTCS
% Solve the system of equations
%IsoFRS = T * vec

vec = T\IsoFRS;

%Prepare final vector
%====================
% The position in isocentric mode is equal to minus the coordinate of the
% isocentre in the TTCS
%The angles in isocentric mode are equal to the angles computed in the
%equipment mode
coord = [-vec(1) -vec(2) -vec(3) rot*180/pi roll*180/pi pitch*180/pi];

return