% RUN FIT MODEL creates calibration file and room setup
% (example for free geometry)
clc
close all
warning off
addpath(genpath('../utils')); % import libraries              

%# INPUTS
% % % roomFout = 'C:\imagx_data\testbench-oblique-10042014\calibration_results\measurementsGetop-06052014\roomSetupTestbenchRADB.xml';
% % % modelFout = 'C:\imagx_data\testbench-oblique-10042014\calibration_results\measurementsGetop-06052014\geometricModelTestbenchRADB.xml';
% % % % calibration measurements or filename (either one or the other):
% % % flexFname = 'C:\imagx_data\testbench-calibration-oblique\calibration_results\measurementsGetop-06052014\deformation\RADB Apr 10 2014 160931 kVp 80 mA 100 ms 80_20140507T142433_deformation.csv';
% % % flexMeas = [];
% % % roomName = 'CasemateUCL';
% % % roomPlace = 'LLN';
% % % roomSystem = 'testbench_oblique';
% % % date = datestr(now,'ddmmyyyy');
% % % time = datestr(now,13);
% % % fpName = 'THALES2630';
% % % displayMode = 1; % it means some plots are displayed at the end of computation
% % % optionPlotFlexmaps = 'cartesian';
% % % fitVec = [1 1 0 1 1 1 1 1 1]; % (Sx,Sy,Sz,Dx,Dy,Dz,Rx,Ry,Rz) fit =1; take average =0
% % % % load geometry and panel attributes
% % % load C:\imagx_data\testbench-oblique-10042014\calibration_results\measurementsGetop\cmaes1.mat;

roomFout = 'C:\imagx_data\shreveport-calibration-26052014\roomSetupSHRRADB.xml';
modelFout = 'C:\imagx_data\shreveport-calibration-26052014\modelSHRRADB.xml';
% calibration measurements or filename (either one or the other):
flexFname = 'C:\imagx_data\shreveport-calibration-26052014\raw_images\180\RADB\1_20140526T125536_deformation.csv';
flexMeas = [];
roomName = 'PGTR';
roomPlace = 'Shreveport';
roomSystem = 'PGTR';
date = datestr(now,'ddmmyyyy');
time = datestr(now,13);
fpName = 'THALES2630';
displayMode = 1; % it means some plots are displayed at the end of computation
optionPlotFlexmaps = 'cartesian';
fitVec = [1 1 0 1 1 1 1 1 1]; % (Sx,Sy,Sz,Dx,Dy,Dz,Rx,Ry,Rz) fit =1; take average =0
% load geometry and panel attributes
load C:\imagx_data\shreveport-calibration-26052014\results\obliqueCalibration_results_20140526T125536.mat;
%# end of INPUTS

%# Ob-1
%# % write room setup xml file:
% % % writeRoomSetupXML(roomFout, roomName, roomPlace, roomSystem, fpName, geometryOb1.axis, geometryOb1.config, panelOb1.Dimension(1), panelOb1.Dimension(2), panelOb1.Pixel, panelOb1.Pixel, geometryOb1.SAD, geometryOb1.AID, panelOb1.Trad, geometryOb1.rx, geometryOb1.ry, geometryOb1.rz, geometryOb1.tx, geometryOb1.ty, geometryOb1.tz);
% % % %# % fit model and write calibration xml file:
% % % fit_model(modelFout, geometryOb1.config, flexMeas, flexFname, roomName, roomPlace, date, time, roomSystem, fpName, geometryOb1.axis, fitVec, displayMode, optionPlotFlexmaps);

% % Ob-2
% write room setup xml file:
writeRoomSetupXML(roomFout, roomName, roomPlace, roomSystem, fpName, geometryOb2.axis, geometryOb2.config, panelOb2.Dimension(1), panelOb2.Dimension(2), panelOb2.Pixel, panelOb2.Pixel, geometryOb2.SAD, geometryOb2.AID, panelOb2.Trad, geometryOb2.rx, geometryOb2.ry, geometryOb2.rz, geometryOb2.tx, geometryOb2.ty, geometryOb2.tz);
% fit model and write calibration xml file:
fit_model(modelFout, geometryOb2.config, flexMeas, flexFname, roomName, roomPlace, date, time, roomSystem, fpName, geometryOb2.axis, fitVec, displayMode, optionPlotFlexmaps);