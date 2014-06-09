function [ detectorPos ] = computeDetectorPosDeformationFRS( detTfrs, TVecIEC, RVecIEC )
%COMPUTEDETECTORPOSDEFORMATIONFRS computes position matrix of flat panel in FRS CS
%from transformation matrix and deformation vector T and R (X,Y,Z) given in local
%IEC CS. T is for flat panel translation, R is for flat panel rotation

detectorPos = detTfrs*rot(RVecIEC(3),[0,0,0])*... % first Rz 
    trans(TVecIEC(1),TVecIEC(2),TVecIEC(3))*...   % translation
    roll(RVecIEC(2),[0,0,0])*pitch(RVecIEC(1),[0,0,0]); % Ry and Rx

end

% -----------------------------------------
% Rotation around Y
function R = roll(ang, off)

a = pi *ang / 180;
ca = cos(a);
sa = sin(a);

R = [ca , 0 , sa , off(1) ;
    0      , 1 , 0      , off(2) ;
    -sa, 0 , ca , off(3) ;
    0      , 0 , 0      , 1 ];

end

% -----------------------------------------
function T = trans(X , Y , Z)

T = [1  ,   0   ,   0   ,   X;
    0  ,   1   ,   0   ,   Y;
    0  ,   0   ,   1   ,   Z;
    0  ,   0   ,   0   ,   1 ];

end

% -----------------------------------------
% Rotation around Z
function R = rot(ang , off)

a = pi *ang / 180;
ca = cos(a);
sa = sin(a);

R = [ca , -sa , 0 , off(1) ;
    sa , ca  , 0 , off(2) ;
    0      , 0       , 1 , off(3) ;
    0      , 0       , 0 , 1 ];
end

% -----------------------------------------
function R = pitch(ang, off)
% rotation around X

a = pi *ang / 180;
ca = cos(a);
sa = sin(a);

R = [1  , 0     , 0         ,off(1) ;
    0  ,ca , -sa   ,off(2) ;
    0  ,sa , ca    ,off(3) ;
    0  , 0     ,    0      , 1 ];

end
