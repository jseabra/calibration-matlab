function [ CP ] = getCylinderPhantomPositions(M, phantomModel)
% GETCYLINDERPHANTOMPOSITIONS Retrieves positions of calibration spheres (cylinder axis along Z) and rotates
% according to a given rotation matrix M
% Cylinder axis should be placed along Yfrs using transformation M
% 
% Inputs: 
%     M - transformation matrix
%     phantomModel - <0=Medcom, 1=iMagX>
%
% Outputs:
%     CP - N x 3 matrix of sphere positions given in FRS
% -----------------------------------------

% Phantom positions in FRS:
%A0 Sphere 0:
CP(1,:) = [0.000, 0.000, -90.000];

%B0 Sphere 0 (0.0,:):   
CP(2,:)=[0.000, -65.000, -75.000];
%B0 Sphere 1 (72.0,:): 
CP(3,:)=[61.819, -20.086, -75.000];
%B0 Sphere 2 (144.0,:): 
CP(4,:)=[38.206, 52.586, -75.000];
%B0 Sphere 3 (216.0,:): 
CP(5,:)=[-38.206, 52.586, -75.000];
%B0 Sphere 4 (288.0,:): 
CP(6,:)=[-61.819, -20.086, -75.000];

%C0 Sphere 0 (0.0,:): 
CP(7,:)=[0.000, -32.500, -60.000];
%C0 Sphere 1 (120.0,:): 
CP(8,:)=[28.146, 16.250, -60.000];
%C0 Sphere 2 (240.0,:): 
CP(9,:)=[-28.146, 16.250, -60.000];

%A1 Sphere 0: 
CP(10,:)=[0.000, 0.000, -45.000];

%B1 Sphere 0 (12.0,:): 
CP(11,:)=[13.514, -63.580, -30.000];
%B1 Sphere 1 (84.0,:): 
CP(12,:)=[64.644, -6.794, -30.000];
%B1 Sphere 2 (156.0,:): 
CP(13,:)=[26.438, 59.380, -30.000];
%B1 Sphere 3 (228.0,:): 
CP(14,:)=[-48.304, 43.493, -30.000];
%B1 Sphere 4 (300.0,:): 
CP(15,:)=[-56.292, -32.500, -30.000];

%C1 Sphere 0 (12.0,:): 
CP(16,:)=[6.757, -31.790, -15.000];
%C1 Sphere 1 (132.0,:): 
CP(17,:)=[24.152, 21.747, -15.000];
%C1 Sphere 2 (252.0,:): 
CP(18,:)=[-30.909, 10.043, -15.000];

%A2 Sphere 0: 
CP(19,:)=[0.000, 0.000, 0.000];

%B2 Sphere 0 (24.0,:): 
CP(20,:)=[26.438, -59.380, 15.000];
%B2 Sphere 1 (96.0,:): 
CP(21,:)=[64.644, 6.794, 15.000];
%%B2 Sphere 2 (168.0,:): 
CP(22,:)=[13.514, 63.580, 15.000];
%B2 Sphere 3 (240.0,:): 
CP(23,:)=[-56.292, 32.500, 15.000];
%B2 Sphere 4 (312.0,:): 
CP(24,:)=[-48.304, -43.493, 15.000];

%C2 Sphere 0 (24.0,:): 
CP(25,:)=[13.219, -29.690, 30.000];
%C2 Sphere 1 (144.0,:): 
CP(26,:)=[19.103, 26.293, 30.000];
%C2 Sphere 2 (264.0,:): 
CP(27,:)=[-32.322, 3.397, 30.000];

%A3 Sphere 0: 
CP(28,:)=[0.000, 0.000, 45.000];

%B3 Sphere 0 (36.0,:): 
CP(29,:)=[38.206, -52.586, 60.000];
%B3 Sphere 1 (108.0,:): 
CP(30,:)=[61.819, 20.086, 60.000];
%B3 Sphere 2 (180.0,:): 
CP(31,:)=[0.000, 65.000, 60.000];
%B3 Sphere 3 (252.0,:): 
CP(32,:)=[-61.819, 20.086, 60.000];
%B3 Sphere 4 (324.0,:): 
CP(33,:)=[-38.206, -52.586, 60.000];

%C3 Sphere 0 (36.0,:): 
CP(34,:)=[19.103, -26.293, 75.000];
%C3 Sphere 1 (156.0,:): 
CP(35,:)=[13.219, 29.690, 75.000];
%C3 Sphere 2 (276.0,:): 
CP(36,:)=[-32.322, -3.397, 75.000];

%A4 Sphere 0: 
if phantomModel,
    CP(37,:) = [NaN, NaN, NaN];
else
    CP(37,:)=[0.000, 0.000, 90.000]; 
end

if size(M,1) ~= 4,
    if size(M,2) ~= 4,
        M = eye(4); 
        disp('Transformation matrix must be 4x4');
        disp('Using identity transformation ...');
    end
end

for i=1:size(CP,1),
    vec = [CP(i,:) 1]';
    vec2 = M* vec;
    CP(i,:)=vec2(1:3)';
end
   
CP( find(isnan(CP(:,1))), : ) = [];

end