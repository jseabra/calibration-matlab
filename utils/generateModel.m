function y = generateModel(t,p)
%GENERATEMODEL sinusoidal model described by parameters p, given as a
%function of t
%
% INPUTS:
% t: variable [gantry angle, for example, given in deg]
% p: 5-element vector

y = p(1) + p(2).*t + p(3)*cos( p(4).*(pi/180)*t + p(5) );

end

