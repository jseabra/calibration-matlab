function M = ieCtotransform(pitch, roll, yaw, tx, ty, tz)
% from pitch yaw roll to 4x4 transform matrix

% !! Angles are given in degrees
pitch = pitch*pi/180;
roll = roll*pi/180;
yaw = yaw*pi/180;

Rz = [  cos(yaw)    -sin(yaw)   0   0
        sin(yaw)    cos(yaw)    0   0
        0           0           1   0
        0           0           0   1];

Rx = [  1   0           0           0
        0   cos(pitch)  -sin(pitch) 0
        0   sin(pitch)  cos(pitch)  0
        0   0           0           1];
    
Ry = [  cos(roll)  0   sin(roll)  0
        0           1   0           0
        -sin(roll) 0   cos(roll)  0
        0           0   0           1];
    
Ti = [  1   0   0   tx
        0   1   0   ty
        0   0   1   tz
        0   0   0   1];

M = Rz * Ti * Rx * Ry;