function pos3D = reconstructPoint3DFromImgPos(imagePosA, imagePosB, panelA, panelB, geometryA, geometryB)
% RECONSTRUCT3DFROMIMGPOS reconstructs point in 3D from
% backprojections of two points given in image coordinates
% Note: if backprojections do not intersect, closest point to the given
% lines is estimated.
% Note: points must be given in image coordinate system (pixels) and
% oriented in beam's eye view
%
% INPUTS:
% imagePosA: position from first DR in image coordinates (pixels)
% imagePosB: position from first DR in image coordinates (pixels)
% panelA/B: a structure with the following attributes:
%     Dimension
%         Pixel
%     Phys2RadT
%      PhysDimX
%      PhysDimY
%    invertGray
%          Trad
%        factor
% geometryA/B: a structure with the following attributes:
%              axis (free/gantry)
%            config (free/gantry)
%               Vec (free/gantry)
%       srcPosition (free)
%         detectorT (free)
%               SAD (free/gantry)
%               AID (free/gantry)
%                rx (free)
%                ry (free)
%                rz (free)
%                tx (free)
%                ty (free)
%                tz (free)
%       gantryAngle (gantry)
% gantryAngleOffset (gantry)


center = [0,0,0]; % if there is an offset wrt isocenter, change this.

X = imagePosA(1);
Y = imagePosA(2);
corner = [panelA.Dimension(1)/2 panelA.Dimension(2)/2];
ptRADA = [(X - corner(1))*panelA.Pixel, (-Y + corner(2))*panelA.Pixel];

if strcmp(geometryA.config,'gantry')
    [srcFRSA, ptFRSA, dirVec] = getBackprojectionRay(ptRADA, geometryA.gantryAngle+geometryA.gantryAngleOffset, center, geometryA.Vec(4:6), geometryA.Vec(1:3), geometryA.Vec(7:9), geometryA.SAD, geometryA.AID);
    
elseif strcmp(geometry.config,'fixed')
    ptFRSA = geometryA.detectorT * ptRADA;
    ptFRSA = ptFRSA(1:3,4);
    srcFRSA = geometryA.srcPosition;
    %vecA = ptFRSA-srcFRSA;
    
end
    
X = imagePosB(1);
Y = imagePosB(2);
corner = [panelB.Dimension(1)/2 panelB.Dimension(2)/2];
ptRADB = [(X - corner(1))*panelB.Pixel, (-Y + corner(2))*panelB.Pixel];

if strcmp(geometryB.config,'gantry')
    [srcFRSB, ptFRSB, dirVec] = getBackprojectionRay(ptRADB, geometryB.gantryAngle+geometryB.gantryAngleOffset, center, geometryB.Vec(4:6), geometryB.Vec(1:3), geometryB.Vec(7:9), geometryB.SAD, geometryB.AID);
    
elseif strcmp(geometry.config,'fixed')
    ptFRSB = geometryB.detectorT * ptRADB;
    ptFRSB = ptFRSB(1:3,4);
    srcFRSB = geometryB.srcPosition;
    vecB = ptFRSB-srcFRSB;
    %vecB = vecB./norm(vecB);
    
end

sourcePositions = [srcFRSA; srcFRSB];
markerPositions = [ptFRSA; ptFRSB];

[pos3D, distances] = lineIntersect3D(sourcePositions, markerPositions);

end