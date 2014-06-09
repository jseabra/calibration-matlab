function [o1, pt_, vec_] = getBackprojectionRay(pt, angle, center, P, T, R, SAD, AID)
% GETBACKPROJECTIONGEOMETRY: computes the ray(line) that starts at a given location on the flat panel plane and ends at the tube focal spot
%
% INPUTS:
% pt = position of point in FP RAD Coordinates [X,Y pixels]
% angle = projection angle [deg]
% center = center of the machine (FRS) [x,y,z mm]
% P = Vector (Px,Py,Pz) with translations [mm] of FP in RAD CS
% T = Vector (Tx,Ty,Tz) with translations [mm] of Source in RAD CS
% R = Vector (Rx,Ry,Rz) with rotations of FP [deg] in RAD CS
% SAD = source to axis distance [mm]
% AID = axis to image distance [mm]
%
% OUTPUTS:
% o1: source position in FRS [mm]
% pt_: detector position in FRS [mm]
% vec_: backprojection vector [mm]

if(AID<0) 
    AID=-AID;
end
SDD = SAD+AID;

P(3) = P(3) -(SDD - SAD);

rotX = [1 0 0
       0 cosd(R(1)) -sind(R(1))
       0 sind(R(1)) cosd(R(1))];  

rotY = [cosd(R(2))     0   sind(R(2))
        0               1   0          
        -sind(R(2))    0   cosd(R(2)) ];   

rotZ = [cosd(R(3)) -sind(R(3)) 0
        sind(R(3)) cosd(R(3)) 0
        0 0 1]; 
   
% The ideal position of X ray tube:
osr = [0 0 SAD];

vector = ones(size(pt,1),1);

pt3 = [pt zeros(size(pt,1),1)];

center_vector = vector * center;
P_vector = vector * P;

pt3 = pt3 - center_vector;

ptc3 = rotZ * ((rotY * rotX * pt3') + P_vector');
% Apply the X-ray tube offset
osc = osr + T;

%Transformation matrices FRS = R * RAD
Mrot = [cosd(angle)     0   sind(angle);
        0               1   0          ;
        -sind(angle)    0   cosd(angle) ];    

% Convert the IEC RAD coordinates into IEC FRS coordinates
ptt = Mrot * ptc3;

o1 = osc  * Mrot';
pt_ = ptt';
vec_ = pt_ - o1;

end