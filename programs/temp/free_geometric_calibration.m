function free_geometric_calibration(room, config)
%FUNCTION_NAME - One line description of what the function or script performs (H1 line)
%Optional file header info (to give more details about the function than in the H1 line)
%Optional file header info (to give more details about the function than in the H1 line)
%Optional file header info (to give more details about the function than in the H1 line)
%
% Syntax:  [output1,output2] = function_name(input1,input2,input3)
%
% Inputs:
%    input1 - Description
%    input2 - Description
%    input3 - Description
%
% Outputs:
%    output1 - Description
%    output2 - Description
%
% Example: 
%    Line 1 of example
%    Line 2 of example
%    Line 3 of example
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2

% Author: Denis Gilbert, Ph.D., physical oceanography
% Maurice Lamontagne Institute, Dept. of Fisheries and Oceans Canada
% email address: gilbertd@dfo-mpo.gc.ca  


% FREE_GEOMETRIC_CALIBRATION: Estimates free geometry deformation flex
clc;
close all;
clear;
addpath(genpath('../iOptim/')); % import optimization library
addpath(genpath('../utils/')); % import utils



%% INPUTS:
%# CASEMATE
% measurements JSE (10/04/2014)
% % % [geometryOb1, geometryOb2, panelOb1, panelOb2] = roomSetup_testbench_measJSE(),
% % % measurements GEO (23/04/2014)
[geometryOb1, geometryOb2, panelOb1, panelOb2] = roomSetup_testbench_measGEOTOPV2()
fnameOb1 = 'C:\imagx_data\testbench-calibration-oblique\cylinder\RADA Apr 10 2014 160925 kVp 80 mA 100 ms 80.raw';
fnameOb2 = 'C:\imagx_data\testbench-calibration-oblique\cylinder\RADB Apr 10 2014 160931 kVp 80 mA 100 ms 80.raw';
objOffset = [-6.3,4.9,44.94]; % cylinder center (A2) is 45mm above the room center (given by the laser)
frsToPhanTransform = pitch(-90,[0,0,0]); disp('Phantom C3 (3spheres) toward +Yfrs');
phantomExtraRot = rot(-90,[0,0,0])*pitch(90,[0,0,0]); % this is required because the phantom was placed with its top toward the casemate ceiling/ cylinder ref. is pointing toward -Xfrs and +Zfrs

%# SHREVEPORT
% [geometryOb1, geometryOb2, panelOb1, panelOb2] = roomSetup_shreveport();
% % % fnameOb1 = 'C:\imagx_data\shreveport-calibration-24032014\raw_images\RADA Feb 14 2014 105746 kVp 70 mA 200 ms 125.raw';
% % % fnameOb2 =
% 'C:\imagx_data\shreveport-calibration-24032014\raw_images\RADB Feb 14
% 2014 105753 kVp 70 mA 200 ms 125.raw';
% fnameOb1 = 'C:\imagx_data\shreveport-calibration-26052014\raw_images\180\RADA\1.raw';
% fnameOb2 = 'C:\imagx_data\shreveport-calibration-26052014\raw_images\180\RADB\1.raw';

objOffset = [0,0,0];
%frsToPhanTransform = pitch(-90,[0,0,0]), disp('Phantom C3 (3spheres) toward +Yfrs');
frsToPhanTransform = pitch(+90,[0,0,0]), disp('Phantom B0 (5spheres) toward +Yfrs');
phantomExtraRot = eye(4);

% set cylinder orientation
phantom = getCylinderPhantomPositions(frsToPhanTransform, 1);
phantom = phantomExtraRot*[phantom,ones(size(phantom,1),1)]';
phantom = phantom'; phantom = phantom(:,1:3); 

% fnameOb1, geometryOb1, panelOb1, fnameOb2, geometryOb2, panelOb2, objOffset);
spAttrbOb1 = sphereDetectionAttributes(geometryOb1.SAD, geometryOb1.AID, panelOb1.Pixel);
spAttrbOb2 = sphereDetectionAttributes(geometryOb2.SAD, geometryOb2.AID, panelOb2.Pixel);


%% calibrate:
optimizationOb1.dof = logical([1 1 0 1 1 0 1 1 1]);
optimizationOb1.method = 'cmaes';

[okOb1, ImgOb1] = readRawImage(fnameOb1, panelOb1.Dimension, panelOb1.Phys2RadT, panelOb1.invertGray);
if (~okOb1)
    disp('Skipping image >> ', fnameOb1, ' >> from processing.');
    return;
end

% detected markers:
[measuredSpheres, isValid] = sphereDetection(ImgOb1, phantom, objOffset, panelOb1, geometryOb1, spAttrbOb1, 1);
if (~isValid)
    return;
end

% compute deformation:
[OptVec, Fout, gof, idx, MINI] = computeDeformation(phantom, objOffset, panelOb1, measuredSpheres, geometryOb1, optimizationOb1, 1);

%% optimize on axial distance only:
optimizationOb1Tmp = optimizationOb1;
geometryOb1Tmp = geometryOb1;
geometryOb1Tmp.Vec = OptVec;
optimizationOb1Tmp.dof = logical([0 0 1 0 0 1 0 0 0]);
[OptVec, Fout, gof, idx, MINI] = computeDeformation(phantom, objOffset, panelOb1, measuredSpheres, geometryOb1Tmp, optimizationOb1Tmp, 1);    
%%
geometryOb1Tmp.Vec = OptVec;
geometryOb1 = geometryOb1Tmp;
% write list of measured spheres to file:
imageTag = [fnameOb1(1:end-4),'_',datestr(now,30)];
csvwrite([imageTag,'_spheres.csv'], measuredSpheres);
results = [OptVec, gof/size(Fout,1)];
% write estimated parameters to file:
fid = fopen([imageTag,'_deformation.csv'], 'w');
fprintf(fid, '%4.3f, %4.3f, %4.3f, %4.3f, %4.3f, %4.3f, %4.3f, %4.3f, %4.3f, %4.3f \n', results);
fclose(fid);
fid = fopen([imageTag,'_MINI.csv'], 'w');
fprintf(fid, '%4.3f\n', MINI);
fclose(fid);

pause;

%% Ob-2:
%optimization:
optimizationOb2.dof = logical([1 1 0 1 1 0 1 1 1]);
optimizationOb2.method = 'cmaes';

[okOb2, ImgOb2] = readRawImage(fnameOb2, panelOb2.Dimension, panelOb2.Phys2RadT, panelOb2.invertGray);
if (~okOb2)
    disp('Skipping image >> ', fnameOb2, ' >> from processing.');
    return;
end
% detected markers:
[measuredSpheres, isValid] = sphereDetection(ImgOb2, phantom, objOffset, panelOb2, geometryOb2, spAttrbOb2, 1);
if (~isValid)
    return;
end
% compute deformation:
[OptVec, Fout, gof, idx, MINI] = computeDeformation(phantom, objOffset, panelOb2, measuredSpheres, geometryOb2, optimizationOb2, 1);


%% optimize on axial distance only:
geometryOb2Tmp = geometryOb2;
geometryOb2Tmp.Vec = OptVec;
optimizationOb2Tmp = optimizationOb2;
optimizationOb2Tmp.dof = logical([0 0 1 0 0 1 0 0 0]);
[OptVec, Fout, gof, idx, MINI] = computeDeformation(phantom, objOffset, panelOb2, measuredSpheres, geometryOb2Tmp, optimizationOb2Tmp, 1);        
%%
geometryOb2Tmp.Vec = OptVec;
geometryOb2 = geometryOb2Tmp;
% write list of measured spheres to file:
imageTag = [fnameOb2(1:end-4),'_',datestr(now,30)];
csvwrite([imageTag,'_spheres.csv'], measuredSpheres);
results = [OptVec, gof/size(Fout,1)];
% write estimated parameters to file:
fid = fopen([imageTag,'_deformation.csv'], 'w');
fprintf(fid, '%4.3f, %4.3f, %4.3f, %4.3f, %4.3f, %4.3f, %4.3f, %4.3f, %4.3f, %4.3f \n', results);
fclose(fid);
fid = fopen([imageTag,'_MINI.csv'], 'w');
fprintf(fid, '%4.3f\n', MINI);
fclose(fid);

% % write all the results into file:
eval(['save obliqueCalibration_results_',datestr(now,30),'.mat;']);


end