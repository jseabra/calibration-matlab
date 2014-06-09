function [panel, geometry, phantom, objOffset, optimization, data, spAttrb] = calibration_setup_upenngtr4_RADB()
clc;
close all;
clear;
addpath(genpath('../iOptim/')); % import optimization library
addpath(genpath('../utils/')); % import utils

% flat panel attributes:
panel.Dimension = [2304 3200]; % in PIXELS (X,Y)
panel.Pixel = 0.127; % pixel size [mm]
panel.Phys2RadT = 'rot0'; % rotation needed to put flat panel in RAD CS
panel.PhysDimX = panel.Dimension(1)/2*panel.Pixel;  
panel.PhysDimY = panel.Dimension(2)/2*panel.Pixel;
panel.invertGray = 1;
panel.factorLow = 50000; % non-normalized 
panel.factorHigh = 0; % non-normalized

% geometry:
geometry.axis = 'RADB'; % axis to calibrate
geometry.config = 'gantry'; % 'gantry','free'
geometry.Vec = [0 0 0 0 0 0 0 0 0];   % (Sx, Sy, Sz, Dx, Dy, Dz, Rx, Ry, Rz)
geometry.gantryAngleOffset = -90; % offset projection to gantry angle
geometry.SAD = 2875; % source to axis distance [mm]                 
geometry.AID = 3470-2875; % axis to detector distance [mm]                 
geometry.gantryAngle = 0; % = no longer used (read from image header)

% cylinder phantom:
frsToPhanTransform = pitch(-90,[0,0,0]); %disp('Phantom C3 (3spheres) toward +Yfrs');
cylinderModel = 1; %=1 casemate, 0 medcom
phantom = getCylinderPhantomPositions(frsToPhanTransform, cylinderModel);
objOffset = [0,0,0]; % isocenter offset (x,y,z) [mm]

% optimization:
optimization.dof = logical([0 0 0 1 1 1 1 1 1]);
optimization.method = 'cmaes'; % 'simplex', 'cmaes', 'lsqnonlin', 'levmar', 'cmb'      

% data management:
data.inputDir = 'C:\imagx_data\Calib_AdaPT_P4'; % root directory with calibration data
data.imagesDir = [data.inputDir,'\RADB']; % directory where images were stored
data.imagesRegExp = 'RADB'; % regular expression to identify images in directory
data.sortBy = 'none'; % 'date'or 'none'
data.anglesFromHeaders = 1;
data.tempDir = [data.inputDir,'\calibration\RADB\6DOF'];
% define if data.anglesFromHeaders=0, otherwise set to ""
data.anglesFname = [data.inputDir,'\gantryAngles_CW.txt']; % file must contain an angle for each image

% sphere attributes:
dratio_min  = 0.8; % relative lower tol ratio of sphere diam
dratio_max  = 2; % relative upper tol ratio of sphere diam
mf = (geometry.SAD  + geometry.AID )/geometry.SAD; % magnification factor
searchSize = 20; % search by given distance around theoretical marker [mm]
spAttrb.beadSize = 1.5; % diameter of metallic sphere [mm]
projDiam = (mf * spAttrb.beadSize)/panel.Pixel; % projected sphere diameter [pixels]
dMin = dratio_min*projDiam; dMax = dratio_max*projDiam;
sMin = dMin*dMin; sMax = dMax*dMax;
spAttrb.projDiam = projDiam;
spAttrb.rSize = floor(searchSize/panel.Pixel + 0.5); % region size [pixels]
spAttrb.sMin  = sMin; % min nr elements [pixels]
spAttrb.sMax  = sMax; % max nr elements [pixels]
spAttrb.dMin  = dMin; % min diameter [pixels]
spAttrb.dMax  = dMax; % max diameter [pixels]
spAttrb.ratio = 0.8; % ratio vertical/horizontal axis
spAttrb.min2diff = projDiam/2; % sphere max resolution (pixels)
spAttrb.threshold = -0.3; % threshold for local binarization
spAttrb.do_filter = 1; % median filter: 0=off, 1=on
spAttrb.ecc = 0.5; % eccentricity
spAttrb.imopen_param = 0.3; % image open parameter (to remove clutter)
spAttrb.detectionMode = 'auto'; % detect spheres: auto, manual

end