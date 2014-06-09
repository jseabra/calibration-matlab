function transformToRad = createTransform2RadFromPanelProps(rotationType, PhysDimX, PhysDimY)
% CREATETRANSFORM2RADFROMPANELPROPS creates a transformation matrix 
% 'transformToRad' used in the imagx room setup file
% from the flat panel properties

% Remark: rotation names correspond to the transformation needed to put the
% flat panel in IEC (Beam's eye view). E.g.: if the flat panel is upside
% down, this would be a 'rotationX'.
% The transformations are the ones read by imagx. They do not match their
% rotation names because imagx does a View transform by dafault.

if strcmp(rotationType,'rotationX'),
    transformToRad = [1  0  0  -PhysDimX;
                      0 -1  0   PhysDimY;
                      0  0 -1   0
                      0  0  0   1];
                  
elseif strcmp(rotationType,'rotationY'),
    transformToRad = [-1  0  0  PhysDimX;
                       0  1  0  -PhysDimY;
                       0  0 -1   0
                       0  0  0   1];
                  
elseif strcmp(rotationType,'rotationZ'),  
    transformToRad = [-1  0   0   PhysDimX;
                      0   1   0   -PhysDimY;
                      0   0  -1    0
                      0   0   0    1];

elseif strcmp(rotationType,'identity'),
    transformToRad = [1  0   0   -PhysDimX;
                      0  -1  0   PhysDimY;
                      0  0   1   0
                      0  0   0   1];    

else
    disp('Unknown rotation type.');
    return;
end