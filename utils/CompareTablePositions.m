function [diffPitch, diffRoll, diffYaw, diffTx, diffTy, diffTz] = CompareTablePositions(pitchPosInit, rollPosInit, yawPosInit, txPosInit, tyPosInit, tzPosInit,...
    pitch, roll, yaw, tx, ty, tz)

MInit = ieCtotransform(pitchPosInit, rollPosInit, yawPosInit, txPosInit, tyPosInit, tzPosInit);

M = ieCtotransform(pitch, roll, yaw, tx, ty, tz);

[diffPitch, diffRoll, diffYaw, diffTx, diffTy, diffTz] = iecPosition (inv(M) * inv(inv(MInit)));