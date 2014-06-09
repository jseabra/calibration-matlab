%=====================================================
% Convert a transformation matrix into the translation and
% the rotations coordinates expressed in the EQUIPMENT PPS format
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
%	NB: the parameters of the PPSoffset function must be correstly configured for this function
%		to work correctly
%
% Author: rla
%	Date: 13/3/06
% MODIFICATION
%   25/2/09: Deal with undeterminate case cos(Y) = 0 (see MID 19674)
%===========================================================

function coord = Cor_EQUI(T)

offset = PPSoffset;

%Get the angles
%===============
sinY = -T(3,1);
roll = asin(sinY);
%We explicitely state that the cosY is POSITIVE. Indeed, by convention, we
%define -pi/2 < Ry < pi/2 (see MID 19674)

cosY = sqrt(1-T(3,1).^2);
%roll = getangle(cosY,sinY);


if(cosY == 0)
    %We are in the undetermined case (see MID 19674)

    %Define the pitch equal to zero
    cosX = 1;
    sinX = 0;
    pitch =0;
    
    cosZ = T(2,2);
    sinZ = -T(1,2);
    %rot = getangle(cosZ,sinZ);
	rot = atan2(sinZ,cosZ);
    
else
    %We are in the general case
    cosX = T(3,3) / cosY;
    sinX = T(3,2) / cosY;
    %pitch= getangle(cosX,sinX);
	pitch= atan2(T(3,2) , T(3,3));

    cosZ = T(1,1) / cosY;
    sinZ = T(2,1) / cosY;
    %rot = getangle(cosZ,sinZ);
	rot = atan2(T(2,1),T(1,1));
end

%Get the translations
%====================
B = zeros(3,4);

%These are the coefficient forthe following system of equations
% B11 X + B12 Y + B13 Z + B14 =0
% B21 X + B22 Y + B23 Z + B24 =0
% B31 X + B32 Y + B33 Z + B34 =0
B(1,1)= 0;
B(1,2)= 0;
B(1,3)= 1;
B(1,4)= -sinY * offset.XOti + cosY * sinX * offset.YOti ...
   + cosY * cosX * offset.ZOti ...
   		+ cosY * offset.ZOio + offset.ZOoe - T(3,4);

B(2,1)= 0;
B(2,2)= 1;
B(2,3)= 0;
B(2,4)= sinZ * offset.XOti*cosY + sinZ * offset.XOoe + sinZ * sinY * sinX * offset.YOti...
   + sinZ * sinY * cosX * offset.ZOti+ sinZ * sinY * offset.ZOio ...
   + cosZ * cosX * offset.YOti - cosZ * sinX * offset.ZOti + cosZ * offset.YOio ...
   + offset.YOes - T(2,4);

B(3,1)= 1;
B(3,2)= 0;
B(3,3)= 0;
B(3,4)= cosZ * offset.XOti * cosY + cosZ * offset.XOoe + cosZ * sinY * sinX * offset.YOti...
   + cosZ * sinY * cosX * offset.ZOti + cosZ * sinY * offset.ZOio - sinZ * cosX * offset.YOti...
   + sinZ * sinX * offset.ZOti - sinZ * offset.YOio + offset.XOes  - T(1,4);

vec = gettrans(B);

%Prepare final vector
%====================
coord = [vec(1) vec(2) vec(3) rot*180/pi roll*180/pi pitch*180/pi];

return