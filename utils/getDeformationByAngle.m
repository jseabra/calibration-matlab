function [ P,T,R ] = getDeformationByAngle(gantryAngle, coeffs_Px, coeffs_Py, coeffs_Pz, coeffs_Rx, coeffs_Ry, coeffs_Rz, coeffs_Sx, coeffs_Sy, coeffs_Sz)
%GETDEFORMATIONBYANGLE Computes the deformation for a given gantry angle
%(use gantryAngle = 0 for fixed or 'free' geometry.
%
% INPUTS:
% gantryAngle: gantry angle in degrees [use =0 for free geometry]
% coeffs_{Px,Py,Pz,Rx,Ry,Rz,Sx,Sy,Sz} must be a 5-element vector containing
% the flexmap coefficients (a0,...a4)
%
parameters = {'Sx', 'Sy', 'Sz', 'Px', 'Py', 'Pz', 'Rx', 'Ry', 'Rz'};
for i = 1:length(parameters),
    eval( [parameters{i},' = generateModel(gantryAngle , coeffs_', parameters{i} ');'] );
end
P = [Px, Py, Pz];
T = [Sx, Sy, Sz];
R = [Rx, Ry, Rz];

end