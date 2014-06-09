function [measurements,geometry] = cbct_geometric_calibration(panel, geometry, phantom, objOffset, optimization, data, spAttrb)
% CBCT_GEOMETRIC_CALIBRATION: Estimates gantry deformation flex
clc;
close all;
addpath(genpath('../iOptim/')); % import optimization library
addpath(genpath('../utils/')); % import utils

%% ---| INPUTS |---

% display options:
displayMode = 0; %0=release, 1=debug display
debugOn = 0; %1=debug/test image orientation/start deformation, 0=skip
startIdx = 1;
step = 2;

%% end of ---| INPUTS |---

measurements = []; % store calibration measurements here
[ok,fnames] = loadData(data.imagesDir, data.imagesRegExp, data.sortBy);

if (~ok)
    disp('Images were not loaded.');
    disp('Aborting program ...');
    return;
end

if ~data.anglesFromHeaders
   
    gantryAngles = csvread(data.anglesFname);
    
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
if (debugOn)    
    disp(' Entering DEBUG ...  ');
    idx = input('Give image index: ');
 
    imageFname = fnames{idx};    
    imageFullPath = [data.imagesDir,'\',imageFname];
    imageTag = imageFname(1:end-4);
    
    if data.anglesFromHeaders
        gantryAngle=getValueFromXMLByTag([imageFullPath(1:end-4),'.xml'], 'angle')
    else
        gantryAngle=gantryAngles(idx);
    end
    
    [ok, Img] = readRawImage(imageFullPath, panel.Dimension, panel.Phys2RadT, panel.invertGray);
    
    geometryTmp = geometry;
    geometryTmp.gantryAngle = gantryAngle;
    
    tic
    [measuredSpheres, isValid, matchpair] = sphereDetection(Img, phantom, objOffset, panel, geometryTmp, spAttrb, displayMode);
 
    if isValid,
        [OptVec, Fout, gof, idx, MINI] = computeDeformation(phantom, objOffset, panel, measuredSpheres, geometryTmp, optimization, displayMode);
        
        geometryTmp.Vec = OptVec;
%         optimizationTmp = optimization;
%         optimizationTmp.dof = logical([0 0 1 0 0 1 0 0 0]);
%         [OptVec, Fout, gof, idx, MINI] = computeDeformation(phantom, objOffset, panel, measuredSpheres, geometryTmp, optimizationTmp, displayMode);        
    end
    toc
    
    csvwrite([data.tempDir,'\',imageTag,'_SPHERES.csv'], measuredSpheres);    
    % deformation measurements are written according to the gantry angle:
    currMeas = [geometryTmp.gantryAngle, OptVec, gof/size(Fout,1), size(measuredSpheres,1)];
    % write estimated parameters to file:
    fid = fopen([data.tempDir,'\',imageTag,'_FLEX.csv'], 'w');
    fprintf(fid, '%4.3f, %4.3f, %4.3f, %4.3f, %4.3f, %4.3f, %4.3f, %4.3f, %4.3f, %4.3f, %4.3f, %f \n',currMeas);
    fclose(fid);
    
    fid = fopen([data.tempDir,'\',imageTag,'_MINI.csv'], 'w');
    fprintf(fid, '%4.3f\n', MINI);
    fclose(fid);
  
    keyboard

end

% ///// on RELEASE mode
stopIdx = numel(fnames);

for k = startIdx : step : stopIdx,
    
    imageFname = fnames{k};
    imageTag = imageFname(1:end-4);
    imageFullPath = [data.imagesDir,'\',imageFname];
    
    if data.anglesFromHeaders
        gantryAngle=getValueFromXMLByTag([imageFullPath(1:end-4),'.xml'], 'angle');
    else
        gantryAngle=gantryAngles(k);
    end

    [ok, Img] = readRawImage(imageFullPath, panel.Dimension, panel.Phys2RadT, panel.invertGray);
    
    if (~ok)
        disp('Skipping image >> ', imageFname, ' >> with gantryAngle = ', gantryAngle, ' from processing ...');
        continue
    end
    
    % update geometry struct with current gantry angle:
    geometryTmp = geometry;
    geometryTmp.gantryAngle = gantryAngle,
    
    tic
    disp('Finding calibration markers ... ');
    [measuredSpheres, isValid, matchpair] = sphereDetection(Img, phantom, objOffset, panel, geometryTmp, spAttrb, displayMode);
    
    if (~isValid)
        disp(['Skipping image >> ', imageFname, ' >> with gantryAngle = ', gantryAngle, ' from processing ...']);
        continue
    end
    
    disp('Estimating geometrical deformation ... ');
    [OptVec, Fout, gof, idx, MINI] = computeDeformation(phantom, objOffset, panel, measuredSpheres, geometryTmp, optimization, displayMode);
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

    
%     %% optimize on axial distance only:
%     geometryTmp.Vec = OptVec;
%     optimizationTmp = optimization;
%     optimizationTmp.dof = logical([0 0 1 0 0 1 0 0 0]);
%     [OptVec, Fout, gof, idx, MINI] = computeDeformation(phantom, objOffset, panel, measuredSpheres, geometryTmp, optimizationTmp, displayMode);
    %%
    toc
    pause(3);
    disp(['[Progress = ', num2str(k/stopIdx * 100), '%]']);
    % write list of measured spheres to file:
    csvwrite([data.tempDir,'\',imageTag,'_SPHERES.csv'], measuredSpheres);

    % deformation measurements are written according to the gantry angle:
    currMeas = [geometryTmp.gantryAngle, OptVec, gof/size(Fout,1), size(measuredSpheres,1)];
    % write estimated parameters to file:
    measurements = [measurements; currMeas];
    fid = fopen([data.tempDir,'\',imageTag,'_FLEX.csv'], 'w');
    fprintf(fid, '%4.3f, %4.3f, %4.3f, %4.3f, %4.3f, %4.3f, %4.3f, %4.3f, %4.3f, %4.3f, %4.3f, %f \n',currMeas);
    fclose(fid);
    
    fid = fopen([data.tempDir,'\',imageTag,'_MINI.csv'], 'w');
    fprintf(fid, '%4.3f\n', MINI);
    fclose(fid);
    
    csvwrite([data.tempDir,'\FLEXMAP.csv'], measurements);

    close;
end
    
end % end(main)

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

