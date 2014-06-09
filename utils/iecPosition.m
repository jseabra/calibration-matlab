% Gives the 6 values of a transform in the iec order given the transform
% matrix.
% Return in degrees
function [pitch, roll, yaw, tx, ty, tz] = iecPosition(transform)

almostOne = 0.999;

sinPitch = transform(3, 2);
pitch = asin(sinPitch);

if (abs(sinPitch) < almostOne) 
roll = atan2(-transform(3, 1), transform(3, 3));
yaw = atan2(-transform(1, 2), transform(2, 2));
else
        roll = 0.0;
        if (sinPitch < 0.0)
            yaw = atan2(-transform(1, 3), transform(1, 1));
        else
            yaw = atan2(transform(2, 1), transform(1, 1));
        end
end

tz = transform(3, 4);
tx = transform(1, 4) * cos(yaw) + transform(2, 4) * sin(yaw);
ty = transform(2, 4) * cos(yaw) - transform(1, 4) * sin(yaw);

pitch = pitch*180/pi;
roll = roll*180/pi;
yaw = yaw*180/pi;
