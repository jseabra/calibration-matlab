function [measurements,geometry] = static_geometric_calibration()
% STATIC_GEOMETRIC_CALIBRATION: Estimates deformation based on a single DR
% defined for a given gantry angle
clc;
close all;
clear;
addpath(genpath('../iOptim/')); % import optimization library
addpath(genpath('../utils/')); % import utils

%% ---| INPUTS |---
% display options:
displayMode = 1; %0=release, 1=debug display

% data management:
inputDir = 'C:\imagx_data\20091222_T07Essen_GTR2'; % root directory with calibration data
imagesDir = [inputDir,'\270']; % directory where images is stored
imageFname = '270_cw_B.raw';
anglesFromHeaders = 0; % =1 if angle is read from image header, =0 is angle is given explicitly 
tempDir = [inputDir,'\results-28042014'];

% load configurations:
addpath(inputDir); % import utils
[panel, geometry, phan, optimization] = calibration_config_essengtr2()
spAttrb = sphereDetectionAttributes_essengtr2(geometry.SAD, geometry.AID, panel.Pixel)
%% end of ---| INPUTS |---

imageTag = imageFname(1:end-4);
imageFullPath = [imagesDir,'\',imageFname];
    
if anglesFromHeaders
    gantryAngle=getValueFromXMLByTag([imageFullPath(1:end-4),'.xml'], 'angle');
else
    gantryAngle=geometry.gantryAngle;
end

[ok, Img] = readRawImage(imageFullPath, panel.Dimension, panel.Phys2RadT, panel.invertGray);

if (~ok)
    disp('Skipping image >> ', imageFname, ' >> with gantryAngle = ', gantryAngle, ' from processing ...');
    return
end

tic
disp('Finding calibration markers ... ');
% phantom:
phantom = getCylinderPhantomPositions(phan.transform, phan.model);
[measuredSpheres, isValid] = sphereDetection(Img, phantom, phan.offset, panel, geometry, spAttrb, displayMode);

if (~isValid)
    disp('Skipping image >> ', imageFname, ' >> with gantryAngle = ', gantryAngle, ' from processing ...');
    return
end
    
disp('Estimating geometrical deformation ... ');
OptVec = computeDeformation(phantom, phan.offset, panel, measuredSpheres, geometry, optimization, displayMode);
% % % gof = sum of distances (not squared!) between measured markers
% % % and estimated projection
% % % (in pixels) as given by the optimizer
% % % MINI = distance (not squared!) between measured markers and estimated projection
% % % (in pixels) as calculated doing mean(M - Fout)
% % % Remark: if there are outliers, the gof would be strongly affected.
% % % Suggestion: check MINI (if there is a distance clearly higher than
% % % others, it means the sphere detection detected some wrong point that has
% % % no matching)
% % % remove the point: MINI(idx_of_outlier) = [];
% % % then, calculate average distance mean(MINI)
    
%% optimize on axial distance only:
geometry.Vec = OptVec;
optimization.dof = logical([0 0 1 0 0 1 0 0 0]);
[OptVec, Fout, gof, idx, MINI] = computeDeformation(phantom, phan.offset, panel, measuredSpheres, geometry, optimization, displayMode);
%%
toc

% write list of measured spheres to file:
csvwrite([tempDir,'\',imageTag,'_SPHERES.csv'], measuredSpheres);
% deformation measurements are written according to the gantry angle:
currMeas = [geometry.gantryAngle, OptVec, gof/size(Fout,1), size(measuredSpheres,1)];
% write estimated parameters to file:
fid = fopen([tempDir,'\',imageTag,'_FLEX.csv'], 'w');
fprintf(fid, '%4.3f, %4.3f, %4.3f, %4.3f, %4.3f, %4.3f, %4.3f, %4.3f, %4.3f, %4.3f, %4.3f, %f \n',currMeas);
fclose(fid);
    
end % end(main)