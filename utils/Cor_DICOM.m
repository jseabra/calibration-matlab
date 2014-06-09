%=====================================================
% Convert a transformation matrix into the translation and
% the rotations expressed in the DICOM format
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
%	NB: the parameters of the DICOMoffset function must be correstly configured for this function
%		to work correctly
%
% Author: rla
%	Date: 13/3/06
% MODIFICATION
%   25/2/09: Deal with undeterminate case cos(Y) = 0 (see MID 19674)
%===========================================================

function coord = Cor_DICOM(T)

offset = DICOMoffset;

%Get the angles
%===============
sinX = T(3,2);

%We explicitely state that the cosY is POSITIVE. Indeed, by convention, we
%define -pi/2 < Rx < pi/2 (see MID 19674)
cosX = sqrt(1-T(3,2).^2);

%pitch = getangle(cosX,sinX);
pitch = asin(sinX);

if(cosX == 0)
    %We are in the undetermined case (see MID 19674)

    %Define the roll equal to zero
    cosY = 1;
    sinY = 0;
    roll =0;
    
    cosZ = -T(2,3)/T(3,2);
    sinZ = T(2,1);
    %rot = getangle(cosZ,sinZ);
	rot = atan2(sinZ ,cosZ );
    
else
    %We are in the general case

    cosY = T(3,3) / cosX;
    sinY = -T(3,1) / cosX;
    %roll = getangle(cosY,sinY);
	roll=atan2(-T(3,1) ,T(3,3));

    cosZ = T(2,2) / cosX;
    sinZ = -T(1,2) / cosX;
    %rot = getangle(cosZ,sinZ);
	rot = atan2(-T(1,2),T(2,2));
end

%Get the translations
%====================
B = zeros(3,4);

%These are the coefficient forthe following system of equations
% B11 X + B12 Y + B13 Z + B14 =0
% B21 X + B22 Y + B23 Z + B24 =0 
% B31 X + B32 Y + B33 Z + B34  =0

B(1,1)= -cosX * sinY;
B(1,2)= sinX;
B(1,3)= cosX * cosY;
B(1,4)= sinX * offset.Ydelt + cosX * cosY * offset.Zdelt - T(3,4);

B(2,1)= sinZ * cosY + cosZ * sinX * sinY;
B(2,2)= cosZ * cosX;
B(2,3)= sinZ * sinY - cosZ * sinX * cosY;
B(2,4)= cosZ * cosX * offset.Ydelt - cosZ * sinX * cosY * offset.Zdelt ...
   		+ sinZ * sinY *offset.Zdelt - T(2,4);

B(3,1)= cosZ * cosY - sinZ * sinX * sinY;
B(3,2)= -sinZ * cosX;
B(3,3)= cosZ * sinY + sinZ * sinX * cosY;
B(3,4)= cosZ * sinY * offset.Zdelt + sinZ * sinX * cosY * offset.Zdelt...
   		-sinZ*cosX*offset.Ydelt - T(1,4);

vec = gettrans(B);

%Prepare final vector
%====================
coord = [vec(1) vec(2) vec(3) rot*180/pi roll*180/pi pitch*180/pi];

return