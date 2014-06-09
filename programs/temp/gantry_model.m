% RUN FIT MODEL creates calibration file and room setup
% (example for gantry geometry)
clc
close all
warning off
addpath(genpath('../utils')); % import libraries

% roomFout = 'C:\imagx_data\cylinder-static-11042014\tmp-noRot\roomSetupTestbenchRADB.xml';
% modelFout = 'C:\imagx_data\cylinder-static-11042014\tmp-noRot\geometricModelTestbenchRADB.xml';
% % calibration measurements or filename (either one or the other):
% flexFname = 'C:\imagx_data\cylinder-static-11042014\tmp-noRot\FLEXMAP.csv';
% flexMeas = [];
% 
% roomName = 'CasemateUCL';
% roomPlace = 'LLN';
% roomSystem = 'testbench';
% date = datestr(now,'ddmmyyyy');
% time = datestr(now,13);
% fpName = 'THALES4343RF';
% 
% displayMode = 1; % it means some plots are displayed at the end of computation
% optionPlotFlexmaps = 'cartesian';
% fitVec = [1 1 0 1 1 1 1 1 1]; % (Sx,Sy,Sz,Dx,Dy,Dz,Rx,Ry,Rz) fit =1; take average =0
% 
% % flat panel attributes:
% panel.Dimension = [2880 2881]; % in PIXELS (X,Y)
% panel.Pixel = 0.148; % pixel size [mm]
% panel.Phys2RadT = 'rotX'; % rotation needed to put flat panel in RAD CS
% panel.PhysDimX = panel.Dimension(1)/2*panel.Pixel;  
% panel.PhysDimY = panel.Dimension(2)/2*panel.Pixel;
% T = createTransform2RadFromPanelProps(panel);
% % geometry:
% geometry.axis = 'RADB'; % axis to calibrate
% geometry.config = 'gantry'; % 'gantry','free'
% geometry.gantryAngleOffset = -90; % offset projection to gantry angle
% geometry.SAD = 2834; % source to axis distance [mm]                 
% geometry.AID = 559; % axis to detector distance [mm]                 

% roomFout = 'C:\imagx_data\cylinder-static-11042014\tmp-noRot\roomSetupTestbenchRADA.xml';
% modelFout = 'C:\imagx_data\cylinder-static-11042014\tmp-noRot\geometricModelTestbenchRADA.xml';
% % calibration measurements or filename (either one or the other):
% flexFname = 'C:\imagx_data\cylinder-static-11042014\tmp-noRot\FLEXMAP.csv';
% flexMeas = [];
% 
% roomName = 'CasemateUCL';
% roomPlace = 'LLN';
% roomSystem = 'testbench';
% date = datestr(now,'ddmmyyyy');
% time = datestr(now,13);
% fpName = 'THALES4343RF';
% 
% displayMode = 1; % it means some plots are displayed at the end of computation
% optionPlotFlexmaps = 'cartesian';
% fitVec = [1 1 0 1 1 1 1 1 1]; % (Sx,Sy,Sz,Dx,Dy,Dz,Rx,Ry,Rz) fit =1; take average =0
% 
% % flat panel attributes:
% panel.Dimension = [2880 2881]; % in PIXELS (X,Y)
% panel.Pixel = 0.148; % pixel size [mm]
% panel.Phys2RadT = 'rotX'; % rotation needed to put flat panel in RAD CS
% panel.PhysDimX = panel.Dimension(1)/2*panel.Pixel;  
% panel.PhysDimY = panel.Dimension(2)/2*panel.Pixel;
% T = createTransform2RadFromPanelProps(panel);
% % geometry:
% geometry.axis = 'RADA'; % axis to calibrate
% geometry.config = 'gantry'; % 'gantry','free'
% geometry.gantryAngleOffset = 0; % offset projection to gantry angle
% geometry.SAD = 2834; % source to axis distance [mm]                 
% geometry.AID = 559; % axis to detector distance [mm]  


% roomFout = 'C:\imagx_data\20091222_T07Essen_GTR2\results-28042014\roomSetupEssenGTR2_RADB.xml';
% modelFout = 'C:\imagx_data\20091222_T07Essen_GTR2\results-28042014\geometricalCalibrationEssenGTR2_RADB.xml';
% % calibration measurements or filename (either one or the other):
% flexFname = 'C:\imagx_data\20091222_T07Essen_GTR2\results-28042014\270_cw_B_FLEX.csv';
% flexMeas = [];
% 
% roomName = 'GTR2';
% roomPlace = 'Essen';
% roomSystem = 'PROTEUS235';
% date = datestr(now,'ddmmyyyy');
% time = datestr(now,13);
% fpName = 'Varian4030A';
% 
% displayMode = 1; % it means some plots are displayed at the end of computation
% optionPlotFlexmaps = 'cartesian';
% fitVec = [1 1 0 1 1 1 1 1 1]; % (Sx,Sy,Sz,Dx,Dy,Dz,Rx,Ry,Rz) fit =1; take average =0
% 
% % flat panel attributes:
% panel.Dimension = [2304 3200]; % in PIXELS (X,Y)
% panel.Pixel = 0.127; % pixel size [mm]
% panel.Phys2RadT = 'rotY'; % rotation needed to put flat panel in RAD CS
% panel.PhysDimX = panel.Dimension(1)/2*panel.Pixel;  
% panel.PhysDimY = panel.Dimension(2)/2*panel.Pixel;
% T = createTransform2RadFromPanelProps(panel);
% % geometry:
% geometry.axis = 'RADB'; % axis to calibrate
% geometry.config = 'gantry'; % 'gantry','free'
% geometry.gantryAngleOffset = -90; % offset projection to gantry angle
% geometry.SAD = 2875; % source to axis distance [mm]
% geometry.AID = 3470 - geometry.SAD; % axis to detector distance [mm]


roomFout = 'C:\imagx_data\Calib_AdaPT_P4\calibration\RADA\roomSetup.xml';
modelFout = 'C:\imagx_data\Calib_AdaPT_P4\calibration\RADA\calibration.xml';
% calibration measurements or filename (either one or the other):
flexFname = 'C:\imagx_data\Calib_AdaPT_P4\calibration\RADA\6DOF\FLEXMAP.csv';
flexMeas = [];

roomName = 'GTR4';
roomPlace = 'UPenn';
roomSystem = 'PROTEUS235';
date = datestr(now,'ddmmyyyy');
time = datestr(now,13);
fpName = 'Varian4030A';

displayMode = 1; % it means some plots are displayed at the end of computation
optionPlotFlexmaps = 'cartesian';
fitVec = [1 1 0 1 1 1 1 1 1]; % (Sx,Sy,Sz,Dx,Dy,Dz,Rx,Ry,Rz) fit =1; take average =0

% flat panel attributes:
panel.Dimension = [2304 3200]; % in PIXELS (X,Y)
panel.Pixel = 0.127; % pixel size [mm]
panel.Phys2RadT = 'rotY'; % rotation needed to put flat panel in RAD CS
panel.PhysDimX = panel.Dimension(1)/2*panel.Pixel;  
panel.PhysDimY = panel.Dimension(2)/2*panel.Pixel;
T = createTransform2RadFromPanelProps(panel);
% geometry:
geometry.axis = 'RADA'; % axis to calibrate
geometry.config = 'gantry'; % 'gantry','free'
geometry.Vec = [0 0 0 0 0 0 0 0 0];   % (Sx, Sy, Sz, Dx, Dy, Dz, Rx, Ry, Rz)
geometry.gantryAngleOffset = 0; % offset projection to gantry angle
geometry.SAD = 1511; % source to axis distance [mm]                 
geometry.AID = 2108-1511; % axis to detector distance [mm]                 
geometry.gantryAngle = 0; % = no longer used (read from image header)





roomFout = 'C:\imagx_data\Calib_AdaPT_P4\calibration\RADB\roomSetup.xml';
modelFout = 'C:\imagx_data\Calib_AdaPT_P4\calibration\RADB\calibration.xml';
% calibration measurements or filename (either one or the other):
flexFname = 'C:\imagx_data\Calib_AdaPT_P4\calibration\RADB\6DOF\FLEXMAP.csv';
flexMeas = [];

roomName = 'GTR4';
roomPlace = 'UPenn';
roomSystem = 'PROTEUS235';
date = datestr(now,'ddmmyyyy');
time = datestr(now,13);
fpName = 'Varian4030A';

displayMode = 1; % it means some plots are displayed at the end of computation
optionPlotFlexmaps = 'cartesian';
fitVec = [1 1 0 1 1 1 1 1 1]; % (Sx,Sy,Sz,Dx,Dy,Dz,Rx,Ry,Rz) fit =1; take average =0

% flat panel attributes:
panel.Dimension = [2304 3200]; % in PIXELS (X,Y)
panel.Pixel = 0.127; % pixel size [mm]
panel.Phys2RadT = 'rotY'; % rotation needed to put flat panel in RAD CS
panel.PhysDimX = panel.Dimension(1)/2*panel.Pixel;  
panel.PhysDimY = panel.Dimension(2)/2*panel.Pixel;
T = createTransform2RadFromPanelProps(panel);

geometry.axis = 'RADB'; % axis to calibrate
geometry.config = 'gantry'; % 'gantry','free'
geometry.Vec = [0 0 0 0 0 0 0 0 0];   % (Sx, Sy, Sz, Dx, Dy, Dz, Rx, Ry, Rz)
geometry.gantryAngleOffset = -90; % offset projection to gantry angle
geometry.SAD = 2875; % source to axis distance [mm]                 
geometry.AID = 3470-2875; % axis to detector distance [mm]                 
geometry.gantryAngle = 0; % = no longer used (read from image header)


























% write room setup xml file:
%writeRoomSetupXML(roomFout, roomName, roomPlace, roomSystem, fpName, geometry.axis, geometry.config, panel.Dimension(1), panel.Dimension(2), panel.Pixel, panel.Pixel, geometry.SAD, geometry.AID, T, geometry.gantryAngleOffset) 

% fit model and write calibration xml file:
fit_model(modelFout, geometry.config, flexMeas, flexFname, roomName, roomPlace, date, time, roomSystem, fpName, geometry.axis, fitVec, displayMode, optionPlotFlexmaps);

