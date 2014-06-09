function geometric_calibration(geometry, panel, phantom, data_mngt, config)
%========================================================================== 
%GEOMETRIC_CALIBRATION - Performs calibration of free/gantry geometry
%
% Syntax:  geometric_calibration(geometry, panel, phantom, data_mngt,
% config)
%
% Inputs (!!! Check room_configuration folder for inputs !!!):
%    panel - panel attributes. Structure as given in (*)
%    geometry - geometry parameters. Structure as given in (**)
%    phantom - phantom specifications. Structure as given in (***)
%    data - data/results handler. Structure as given in (****)  
%    config - algorithm configuration. Structure as givne in (*****)
%
% (*) Structure of input 'panel':
%    panel.Name = <flat_panel_model>; e.g.: THALES2630CS
%    panel.Dimension = [<DIM_X_PIXELS> <DIM_Y_PIXELS>];
%    panel.Pixel = <PIXEL_SIZE_MM>;
%    panel.Phys2RadT = <'identity','rotationX','rotationY','rotationZ'>; 
%    panel.PhysDimX = panel.Dimension(1)/2*panel.Pixel;  
%    panel.PhysDimY = panel.Dimension(2)/2*panel.Pixel;
%    panel.invertGray = <0,1>;
%    panel.factorLow = <ADJUST_LUT_TO_LOW_FACTOR> 
%    panel.factorHigh = <ADJUST_LUT_TO_HIGH_FACTOR>
%
% (**) Structure of input 'geometry':
%    geometry.roomName = <ROOM_NAME>
%    geometry.roomSystem = <ROOM_SYSTEM>; E.g.: "Proteus235", "ObliqueSetup"
%    geometry.roomPlace
%    geometry.date = <DATE_OF_INSTALLATION>; % Format: ddMMyy
%    geometry.time = <TIME_OF_INSTALLATION>; % Format: hhMMss
%    geometry.axis = <'RADA','RADB','RADC',...>; % axis to calibrate
%    geometry.config = <'gantry','free'>
%    geometry.Vec = [Sx, Sy, Sz, Dx, Dy, Dz, Rx, Ry, Rz]; % starting guess
%    geometry.gantryAngleOffset = <PROJECTION_TO_GANTRY_ANGLE_DEGREES> ['gantry' geometry only]
%    geometry.SAD = <SOURCE_TO_AXIS_DISTANCE_MM>                 
%    geometry.AID = <AXIS_TO_DETECTOR_DISTANCE_MM>    
%    geometry.srcPosition = <(X_frs,Y_frs,Z_frs)> ['free' geometry only]
%    geometry.detectorT = <DETECTOR_MATRIX_frs>; 4x4 ['free' geometry only] 
%    geometry.rx = <PITCH_frs_DEG> [used in IMAGX only] ['free' geometry only] 
%    geometry.ry = <ROLL_frs_DEG> [used in IMAGX only] ['free' geometry only] 
%    geometry.rz = <YAW_frs_DEG> [used in IMAGX only] ['free' geometry only] 
%    geometry.tx = <trans_X_Frs> [used in IMAGX only] ['free' geometry only] 
%    geometry.ty = <trans_Y_Frs> [used in IMAGX only] ['free' geometry only] 
%    geometry.tz = <trans_Z_Frs> [used in IMAGX only] ['free' geometry only]            
%
% (***) Structure of input 'phantom':
%    phantom.spheres = <list_of_spheres>; % spheres(id,:) = (X_frs,Y_frs,Z_frs)
%    phantom.offset = [X_frs,Y_frs,Z_frs]; % offset central sphere to
%    isocenter
%
% (****) Structure of input 'data_mngt':
%    data_mngt.input = <PATH_TO_DIRECTORY (for 'gantry'), PATH_TO_IMAGE (for 'free')>
%    data_mngt.imagesRegExp = 'RADB'; % regular expression ['gantry' only]
%    data_mngt.sortBy = <'date','none'> % ['gantry' only]
%    data_mngt.anglesFromHeaders = <0,1>; % get gantry angles from xml
%    data_mngt.output = <PATH_TO_DIRECTORY>;
%    data_mngt.anglesFromHeaders = <0,1>;
%    data_mngt.anglesFname = <FILENAME_ANGLES> % ['gantry' only]
%
% (*****) Structure of input 'config':
%    config.dof = [Sx, Sy, Sz, Dx, Dy, Dz, Rx, Ry, Rz]; % <0=constant,
%    1=optimize)
%    config.spAttrb (see ******)
%    config.optimizer = <'cmaes','simplex','lsqnonlin','levmar'>
%    config.debugMode = <0,1>
%    config.dof = [Sx, Sy, Sz, Dx, Dy, Dz, Rx, Ry, Rz]; %  <0=computes
%    average of flexmap, 1=fits model to flexmap>
%    config.remoteOutliersFromModel = <PARAMETER_BETWEEN_0_AND_1>
%
% (******) Structure of input 'config.spAttrb':
%    spAttrb.beadSize = <SPHERE_DIAMETER_MM>
%    spAttrb.projDiam = <PROJECTED_SPHERE_DIAM_PIXELS>
%    spAttrb.rSize = <SEARCH_WINDOW_WIDTH_PIXELS>
%    spAttrb.sMin  = <MIN_NR_OBJ_PIXELS>
%    spAttrb.sMax  = <MAX_NR_OBJ_PIXELS>
%    spAttrb.dMin  = <MIN_DIAM_OBJ_PIXELS>
%    spAttrb.dMax  = <MAX_DIAM_OBJ_PIXELS>
%    spAttrb.ratio = <RATIO_OBJ_H/V>
%    spAttrb.min2diff = <MAX_DIST_BETWEEN_OBJ_PIXELS>
%    spAttrb.threshold = <THRESHOLD_BETWEEN_-1_1>
%    spAttrb.do_filter = <0,1> % median filter
%    spAttrb.ecc = <ECCENTRICITY_0_1> % 1 for ideal symmetry
%    spAttrb.imopen_param = <REMOVE_CLUTTER_0_1>
%    spAttrb.detectionMode = <'auto','manual'>
%
% Other m-files required: see /iOptim and /utils folder
% MAT-files required: none
%
% Author: J Seabra, Ph.D.
% Universite Catholique de Louvain, LLN, Belgique
% email address: mail2jseabra@gmail.com
%==========================================================================
clc;
close all;

% import optimization library:
addpath(genpath('../iOptim/'));
% import utils:
addpath(genpath('../utils/')); 

% check if inputs are valid:
if(~checkPanelIntegrity(panel))
    disp('Panel parameters missing or wrong.');
    return;
end
if(~checkGeometryIntegrity(geometry))
    disp('Geometry parameters missing or wrong.');
    return;
end
if(~checkPhantomIntegrity(phantom))
    disp('Phantom parameters missing or wrong.');
    return;
end
if(~checkConfigIntegrity(config))
    disp('Config parameters missing or wrong.');
    return;
end

% current date and time to use as uid:
DT = now;
% do 'gantry geometry' workflow:
if(strcmp(geometry.config,'gantry'))
    flexmap = gantry_calibration(panel, geometry, phantom, data_mngt, config);
    % write room setup XML:
    writeRoomSetupXML(roomFout, geometry.roomName, geometry.roomPlace, geometry.roomSystem, panel, geometry.axis,...
    geometry.config, geometry.SAD, geometry.AID, geometry.gantryAngleOffset);
else
% do 'free geometry' workflow:
    flexmap = free_calibration(panel, geometry, phantom, data_mngt, config);
    roomFout = [data_mngt.output,'/room_',datestr(DT,30),'_',geometry.axis,'.xml'];
    % write room XML:
    writeRoomSetupXML(roomFout, geometry.roomName, geometry.roomPlace, geometry.roomSystem, panel, geometry.axis,...
        geometry.config, geometry.SAD, geometry.AID, geometry.rx, geometry.ry, geometry.rz, geometry.tx, geometry.ty, geometry.tz);
end

% fit model to flexmap:
 [coeffsS, coeffsD, coeffsR] = fit_model_to_data(geometry.config, flexmap, config.dof, config.remoteOutliersFromModel, config.debugMode);
calibrationFout = [data_mngt.output,'/calibration_',datestr(DT,30),'_',geometry.axis,'.xml'];
% write calibration XML:
writeCalibrationXML(calibrationFout, geometry.roomName, geometry.roomPlace, DT, geometry.roomSystem, panel.Name, geometry.axis, coeffsD, coeffsR, coeffsS)

end 
% end of geometric_calibration() ========================================


% auxiliar functions:
function flexmap = free_calibration(panel, geometry, phantom, data_mngt, config)
% FREE_CALIBRATION - Calibration for free geometry:
flexmap = [];
% load X-ray image:
[ok, Img] = readRawImage(data_mngt.input, panel.Dimension, panel.Phys2RadT, panel.invertGray);
if (~ok)
    disp(['Skipping image >> ', data_mngt.input, ' >> from processing.']);
    return;
end
% detected markers:
[measuredSpheres, isValid] = sphereDetection(Img, phantom.spheres, phantom.offset, panel, geometry, config.spAttrb, config.debugMode);
if (~isValid)
    return;
end
% compute deformation:
[OptVec, Fout, gof, idx, MINI] = computeDeformation(phantom.spheres, phantom.offset, panel, measuredSpheres, geometry, config.optimizer, config.dof, config.debugMode);
% optimize on axial distances only:
configTmp = config;
geometryTmp = geometry;
geometryTmp.Vec = OptVec;
configTmp.dof = logical([0 0 1 0 0 1 0 0 0]);
[OptVec, Fout, gof, idx, MINI] = computeDeformation(phantom.spheres, phantom.offset, panel, measuredSpheres, geometryTmp, configTmp.optimizer, configTmp.dof, configTmp.debugMode);
geometryTmp.Vec = OptVec;
% write list of measured spheres to file:
imageTag = [data_mngt.input(1:end-4),'_',datestr(now,30)];
csvwrite([imageTag,'_SPHERES.csv'], measuredSpheres);
results = [OptVec, gof/size(Fout,1)];
% write estimated parameters to file:
fid = fopen([imageTag,'_FLEX.csv'], 'w');
fprintf(fid, '%4.3f, %4.3f, %4.3f, %4.3f, %4.3f, %4.3f, %4.3f, %4.3f, %4.3f, %4.3f \n', results); fclose(fid);
fid = fopen([imageTag,'_MINI.csv'], 'w');
fprintf(fid, '%4.3f\n', MINI); fclose(fid);
flexmap = OptVec;

end

function flexmap = gantry_calibration(panel, geometry, phantom, data_mngt, config)
% GANTRY_CALIBRATION - Calibration for gantry geometry:
flexmap = [];
% load directory with X-ray images:
[ok,fnames] = loadData(data_mngt.input, data_mngt.imagesRegExp, data_mngt.sortBy);
if (~ok)
    disp('Images were not loaded.');
    disp('Aborting program ...');
    return;
end

if (~data_mngt.anglesFromHeaders)
    gantryAngles = csvread(data_mngt.anglesFname);
    if (isempty(gantryAngles))
        disp('Gantry angles file is empty');
        disp('Aborting program ...');
        return;
    end
    if (numel(gantryAngles)~=numel(fnames))
        disp('Number of gantry angles and images does not match.');
        disp(['Number of gantry angles = ', num2str(numel(gantryAngles))]);
        disp(['Number of images = ', num2str(numel(fnames))]);
        disp('Aborting program ...');
        return;
    end
end

% ///// on DEBUG mode
if (config.debugOn)
    disp(' Entering DEBUG ...  ');
    idx = input('Give image index: ');
    
    imageFname = fnames{idx};
    imageFullPath = [data_mngt.imagesDir,'\',imageFname];
    imageTag = imageFname(1:end-4);
    
    if data_mngt.anglesFromHeaders
        gantryAngle = getValueFromXMLByTag([imageFullPath(1:end-4),'.xml'], 'angle');
    else
        gantryAngle = gantryAngles(idx);
    end
    
    % read current projection:
    [ok, Img] = readRawImage(imageFullPath, panel.Dimension, panel.Phys2RadT, panel.invertGray);
    
    % update geometry struct with current gantry angle:
    geometryTmp = geometry;
    geometryTmp.gantryAngle = gantryAngle;

    % detected markers:
    [measuredSpheres, isValid] = sphereDetection(Img, phantom.spheres, phantom.offset, panel, geometryTmp, config.spAttrb, config.debugMode);
    if isValid,
        % compute deformation:
        [OptVec, Fout, gof, idx, MINI] = computeDeformation(phantom.spheres, phantom.offset, panel, measuredSpheres, geometry, config.optimizer, config.dof, config.debugMode);
        geometryTmp.Vec = OptVec;
        %         optimizationTmp = optimization;
        %         optimizationTmp.dof = logical([0 0 1 0 0 1 0 0 0]);
        %         [OptVec, Fout, gof, idx, MINI] = computeDeformation(phantom, objOffset, panel, measuredSpheres, geometryTmp, optimizationTmp, displayMode);
    end
    return
    
end

% ///// on RELEASE mode
stopIdx = numel(fnames);
measurements = [];

for k = startIdx : step : stopIdx,
    
    imageFname = fnames{k};
    imageTag = imageFname(1:end-4);
    imageFullPath = [data_mngt.imagesDir,'\',imageFname];
    
    if data.anglesFromHeaders
        gantryAngle=getValueFromXMLByTag([imageFullPath(1:end-4),'.xml'], 'angle');
    else
        gantryAngle=gantryAngles(k);
    end
    
    % read current projection:
    [ok, Img] = readRawImage(imageFullPath, panel.Dimension, panel.Phys2RadT, panel.invertGray);
    
    if (~ok)
        disp(['Skipping image >> ', imageFname, ' >> with gantryAngle = ', gantryAngle, ' from processing ...']);
        continue
    end
    
    % update geometry struct with current gantry angle:
    geometryTmp = geometry;
    geometryTmp.gantryAngle = gantryAngle;
    
    % detected markers:
    [measuredSpheres, isValid] = sphereDetection(Img, phantom.spheres, phantom.offset, panel, geometryTmp, config.spAttrb, config.debugMode);
    if (~isValid)
        disp(['Skipping image >> ', imageFname, ' >> with gantryAngle = ', gantryAngle, ' from processing ...']);
        continue
    end
    % compute deformation:
    [OptVec, Fout, gof, idx, MINI] = computeDeformation(phantom, objOffset, panel, measuredSpheres, geometryTmp, optimization, displayMode);
    %     %% optimize on axial distance only:
    %     geometryTmp.Vec = OptVec;
    %     optimizationTmp = optimization;
    %     optimizationTmp.dof = logical([0 0 1 0 0 1 0 0 0]);
    %     [OptVec, Fout, gof, idx, MINI] = computeDeformation(phantom, objOffset, panel, measuredSpheres, geometryTmp, optimizationTmp, displayMode);
    %%
    toc
    disp(['[Progress = ', num2str(k/stopIdx * 100), '%]']);
    % write list of measured spheres to file:
    csvwrite([data.tempDir,'\',imageTag,'_SPHERES.csv'], measuredSpheres);
    
    % deformation measurements are written according to the gantry angle:
    currMeas = [geometryTmp.gantryAngle, OptVec, gof/size(Fout,1), size(measuredSpheres,1)];
    % write estimated parameters to file:
    measurements = [measurements; currMeas];
    fid = fopen([data.tempDir,'\',imageTag,'_FLEX.csv'], 'w');
    fprintf(fid, '%4.3f, %4.3f, %4.3f, %4.3f, %4.3f, %4.3f, %4.3f, %4.3f, %4.3f, %4.3f, %4.3f, %f \n',currMeas); fclose(fid); 
    fid = fopen([data.tempDir,'\',imageTag,'_MINI.csv'], 'w');
    fprintf(fid, '%4.3f\n', MINI); fclose(fid);
    csvwrite([data.tempDir,'\FLEXMAP.csv'], measurements);
    close;
    flexmap = measurements(:,1:10); % gantry angles + 9 dof parameters

end

end

% -----------------------------------------
% Lists the files from a specific directory
function [ok,fnames] = loadData(imagesDir, imagesRegExp, sortBy)
data = dir([imagesDir,'\*',imagesRegExp,'*.raw']);

if (isempty(data))
    disp(['No files exist in >> ', imagesDir, ' >> for regExp = >> ', imagesRegExp]);
    ok = false;
    
else
    if (strcmp(sortBy,'date'))
        S = [data(:).datenum]';
        [S,S] = sort(S);
        fnames = {data(S).name};
    else
        fnames = {data.name};
    end
    ok = true;
    
end

end