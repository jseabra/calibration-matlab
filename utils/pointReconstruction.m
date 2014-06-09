function [absErr, relErr] = pointReconstruction()
%POINTRECONSTRUCTION Enables to reconstruct a point in 3D from 2D points
%selected in two or more DRs.
clc;
close all;
clear;

%% ---| INPUTS |---
% Cylinder phantom:
frsToPhanTransform = pitch(-90,[0,0,0]); %disp('Phantom C3 (3spheres) toward +Yfrs');
cylinderModel = 1; %=1 casemate, 0 otherwise
phantom = getCylinderPhantomPositions(frsToPhanTransform, cylinderModel);
objOffset = [0,0,0]; % isocenter offset (x,y,z) [mm]

% flat panel attributes:
% RadOne:
panelRadOne.Dimension = [2880 2881]; % in PIXELS (X,Y)
panelRadOne.Pixel = 0.148; % pixel size [mm]
panelRadOne.Phys2RadT = 'rotX'; % rotation needed to put flat panel in RAD CS
panelRadOne.PhysDimX = panelRadOne.Dimension(1)/2*panelRadOne.Pixel;  
panelRadOne.PhysDimY = panelRadOne.Dimension(2)/2*panelRadOne.Pixel;
panelRadOne.invertGray = 1;
% RadTwo:
panelRadTwo.Dimension = [2880 2881]; % in PIXELS (X,Y)
panelRadTwo.Pixel = 0.148; % pixel size [mm]
panelRadTwo.Phys2RadT = 'rotX'; % rotation needed to put flat panel in RAD CS
panelRadTwo.PhysDimX = panelRadTwo.Dimension(1)/2*panelRadTwo.Pixel;  
panelRadTwo.PhysDimY = panelRadTwo.Dimension(2)/2*panelRadTwo.Pixel;
panelRadTwo.invertGray = 1;

% Geometry struct:
% RadOne:
geometryRadOne.axis = 'RADA'; % axis to calibrate
geometryRadOne.config = 'gantry'; % 'gantry','free'
geometryRadOne.gantryAngle = 0; 
geometryRadOne.gantryAngleOffset = 0; % offset projection to gantry angle
geometryRadOne.SAD = 2834; % source to axis distance [mm]                 
geometryRadOne.AID = 559; % axis to detector distance [mm]                 
geometryRadOne.model = 'C:\imagx_data\cylinder-static-11042014\tmp-noRot\model.xml';
% RadTwo:
geometryRadTwo.axis = 'RADB'; % axis to calibrate
geometryRadTwo.config = 'gantry'; % 'gantry','free'
geometryRadTwo.gantryAngle = 0;
geometryRadTwo.gantryAngleOffset = -90; % offset projection to gantry angle
geometryRadTwo.SAD = 2834; % source to axis distance [mm]                 
geometryRadTwo.AID = 559; % axis to detector distance [mm]                 
geometryRadTwo.model = 'C:\imagx_data\cylinder-static-11042014\tmp-noRot\model.xml';

% Drs:
% RadOne:
% image at gantry angle = 90deg [for testbench image at RadOne must be gantryAngle+90]
fnameRadOne = 'C:\imagx_data\cylinder-static-11042014\drs\RADB Apr 11 2014 174352 kVp 100 mA 250 ms 125_001.raw';
% RadTwo:
% image at gantry angle = 0deg
fnameRadTwo = 'C:\imagx_data\cylinder-static-11042014\drs\RADB Apr 11 2014 174052 kVp 100 mA 250 ms 125_001.raw';

% label of cylinder sphere to reconstruct:
pointGT = phantom(2,:);

%% end of ---| INPUTS |---

% raw is read and must be displayed in beam's eye view:
[okRadOne, drRadOne] = readRawImage(fnameRadOne, panelRadOne.Dimension, panelRadOne.Phys2RadT, panelRadOne.invertGray);
if (~okRadOne)
    disp('Cannot read >> ', fnameRadOne, '. Aborting.');
    return;
end

[okRadTwo, drRadTwo] = readRawImage(fnameRadTwo, panelRadTwo.Dimension, panelRadTwo.Phys2RadT, panelRadTwo.invertGray);
if (~okRadTwo)
    disp('Cannot read >> ', fnameRadTwo, '. Aborting.');
    return;
end

% beam must consist of: geometry, panel and image
beamOne.geometry = geometryRadOne;
beamOne.panel = panelRadOne;
beamOne.dr = drRadOne;

beamTwo.geometry = geometryRadTwo;
beamTwo.panel = panelRadTwo;
beamTwo.dr = drRadTwo;

[absErr, relErr] = reconstructPoint3D(phantom, objOffset, pointGT, beamOne, beamTwo);

end

