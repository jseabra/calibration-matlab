%=====================================================
% Convert a transformation matrix into the translation and
% the rotations expressed in the IEC-61217 format
%
% INPUT:
%		T: the 4x4 transformation matrix
%
% OUTPUT:
%		Coord: vector with the translation and rotation
%				[X Y Z rot roll pitch]
%                   rot [0, 360]
%                   pitch [-90, 90] (more correctly [0,90] and [270, 360])
%                   roll [0, 360]
%                   In the undeterminate case, pitch = +-90, we define
%                    roll =0
%               distance in cm
%
% Author: rla
%	Date: 17/7/09
%===========================================================

function coord = Cor_IECV(T)

offset.Ydelt = 0;
offset.Zdelt = 0;

%Get the angles
%===============
sinX = T(3,2);
pitch = asin(sinX);

%We explicitely state that the cosY is POSITIVE. Indeed, by convention, we
%define -pi/2 < Rx < pi/2 (see MID 19674)
cosX = sqrt(1-T(3,2).^2);
%pitch = getangle(cosX,sinX);

if(cosX == 0)
    %We are in the undetermined case (see MID 19674)

    %Define the roll equal to zero
    cosY = 1;
    sinY = 0;
    roll =0;
    
    cosZ = T(1,1);
    sinZ = T(2,1);
    %rot = getangle(cosZ,sinZ);
    rot = atan2(T(2,1),T(1,1));
    
else
    %We are in the general case

    cosY = T(3,3) / cosX;
    sinY = -T(3,1) / cosX;
    roll = atan2(-T(3,1), T(3,3));
    %roll = getangle(cosY,sinY);

    cosZ = T(2,2) / cosX;
    sinZ = -T(1,2) / cosX;
    rot = atan2(-T(1,2),T(2,2));
    %rot = getangle(cosZ,sinZ);
end

%Get the translations
%====================
B = zeros(3,4);

%These are the coefficient forthe following system of equations
% B11 X + B12 Y + B13 Z + B14 =0
% B21 X + B22 Y + B23 Z + B24 =0 
% B31 X + B32 Y + B33 Z + B34  =0

B(1,1)= cosZ;
B(1,2)= -sinZ;
B(1,3)= 0;
B(1,4)= - T(1,4) ...
        - offset.Ydelt * sinZ * cosX ...
        + offset.Zdelt * (cosZ*sinY+sinZ*sinX*cosY);

B(2,1)= sinZ ;
B(2,2)= cosZ ;
B(2,3)= 0;
B(2,4)=  - T(2,4)...
        + offset.Ydelt *cosZ * cosX ...
        + offset.Zdelt * (sinZ * sinY - sinX * cosY * cosZ);

B(3,1)= 0;
B(3,2)= 0;
B(3,3)= 1;
B(3,4)= -T(3,4)...
        + offset.Ydelt * sinX...
        + offset.Zdelt * cosX*cosY;

%vec = gettrans(B);

%Prepare final vector
%====================
coord = [0 0 0 rot*180/pi roll*180/pi pitch*180/pi];


return