function [geometry, panel, phantom, data_mngt, config] = configuration_testbench_oblique_ob2()
% This is a configuration for oblique (free) room geometry.
% Copy onto new file and edit to create new room geometry.
clc
clear
addpath(genpath('../../utils/')); % import utils

% panel:
panel.Name = 'THALES2630';
panel.Dimension = [1440 1440];
panel.Pixel = 0.184;
panel.Phys2RadT = 'identity';
panel.PhysDimX = panel.Dimension(1)/2*panel.Pixel;
panel.PhysDimY = panel.Dimension(2)/2*panel.Pixel;
panel.invertGray = 1;
panel.factorLow = 0;
panel.factorHigh = 1

% phantom:
frsToPhanTransform = pitch(-90,[0,0,0]), disp('Phantom C3 (3spheres) toward +Yfrs');
phantom.spheres = getCylinderPhantomPositions(frsToPhanTransform, 1);
phantom.offset = [-6.3,4.9,44.94]; % cylinder center (A2) is 45mm above the room center (given by the laser)

% required because the phantom was placed with its top toward the casemate
% ceiling/ cylinder ref. is pointing toward -Xfrs and +Zfrs
phantomExtraRot = rot(-90,[0,0,0])*pitch(90,[0,0,0]); 
phantom.spheres = phantomExtraRot*[phantom.spheres,ones(size(phantom.spheres,1),1)]';
phantom.spheres = phantom.spheres'; phantom.spheres = phantom.spheres(:,1:3); 

% data mngt:
data_mngt.input = 'C:\imagx_data\testbench-calibration-oblique\cylinder\RADB Apr 10 2014 160931 kVp 80 mA 100 ms 80.raw';
data_mngt.output = 'C:\imagx_data\testbench-calibration-oblique\results\Ob2\'

%% =====  Measurements Geotop
% points sur surface à la sortie du tube Ob-2
TB1 = [148.8174 48.3123 3.8496];
TB2 = [148.7776 48.3409 3.8494];
TB3 = [148.8058 48.3754 3.7723];
TB4 = [148.8488 48.3441 3.7724];
TBF = [148.7478 48.2554 3.8300];
TBC = mean([TB1;TB2;TB3;TB4],1);
% points sur flat panel B
PB1 = [150.5722 50.4319 5.4044];			
PB2	= [150.4896 50.6259 5.2137];			
PB3	= [150.3019 50.4966 5.5826];			
PB4	= [150.2215 50.6858 5.3953];
PBC = [150.4072 50.5784 5.4130];			
PBC_ = mean([PB1;PB2;PB3;PB4],1);		
% isocenter (fil a plomb)
ISO = [150.0000, 50.000, 5.0000];		
% table inclination:
TAB1 = [149.7242 49.6115 5.0201];			
TAB2 = [150.2490 49.6023 5.0206];		
TAB3 = [150.4253 49.8683 4.9731];		
TAB4 = [150.4328 50.2974 4.8974];		
TAB5 = [149.5546 49.8434 4.9796];		
TAB6 = [149.5650 50.4344 4.8755];

vTAB64 = TAB6 - TAB4; vTAB64 = vTAB64./norm(vTAB64,2);
vTAB24 = TAB2 - TAB4; vTAB24 = vTAB24./norm(vTAB24,2);
vTABn = cross(vTAB24,vTAB64);
vZ = [0 0 -1];
%tablePitch =acosd( vTABn*vZ')/(norm(vTABn,2)*norm(vZ,2)) //11 deg?
tablePitch = 10; % measured with telemeter

%% Calculate positions, vectors and transformation matrices:
% direction vectors flat panel B
% direction vector y
vB34 = PB3-PB4;
vB34 = vB34./norm(vB34,2);
% direction vector x
vB24 = PB2-PB4;
vB24 = vB24./norm(vB24,2);
% direction vector z
vBn = cross(vB24,vB34);
% check angle betweetn two vectors
%acosd((vB24*vB34')/(norm(vB24,2)*norm(vB34,2)))
% plane B transformation matrix
TB = [[vB24,0]' [vB34,0]' [vBn,0]' [0 0 0 1]'];

PBC_ = (PBC - ISO)*10^3 % frs (origin ISO)

TB(1:3,4) = PBC_';

% normal vector of tube window B
vBT12 = TB1-TB2; vBT12 = vBT12./norm(vBT12,2);
vBT13 = TB1-TB3; vBT13 = vBT13./norm(vBT13,2);
vBTn = cross(vBT12,vBT13)

%# Determine angle between detectorA-tubeA:
deltaB = acosd(vBn*vBTn'/(norm(vBn)*norm(vBTn)));

% calculate shortest point along vector to the focal spot:
TBI = project_point_to_line_segment(TBC,TBC+vBTn,TBF);
TBI_ = (TBI - ISO)*10^3; % frs (origin ISO)

% SDD, SID, AID:
SDDB = norm(TBI-PBC,2)*10^3;
SADB = norm(TBI-ISO,2)*10^3;
AIDB = SDDB-SADB;

TBs = TB * trans(0,0,SADB); TBs = TBs(1:3,4)';

% diff between GEOTOP measured source positions and those obtained with
% simplied model
abs(TBI_-TBs);

%% plot:
figure,
grid on, hold on,
plot3(TAB1(1),TAB1(2),TAB1(3),'r*'); hold on, text(TAB1(1),TAB1(2),TAB1(3),'TAB1'); hold on, 
plot3(TAB2(1),TAB2(2),TAB2(3),'b*'); hold on, text(TAB2(1),TAB2(2),TAB2(3),'TAB2'); hold on, 
plot3(TAB3(1),TAB3(2),TAB3(3),'m*'); hold on, text(TAB3(1),TAB3(2),TAB3(3),'TAB3'); hold on, 
plot3(TAB4(1),TAB4(2),TAB4(3),'g*'); hold on, text(TAB4(1),TAB4(2),TAB4(3),'TAB4'); hold on, 
plot3(TAB5(1),TAB5(2),TAB5(3),'k*'); hold on, text(TAB5(1),TAB5(2),TAB5(3),'TAB5'); hold on, 
plot3(TAB6(1),TAB6(2),TAB6(3),'c*'); hold on, text(TAB6(1),TAB6(2),TAB6(3),'TAB6'); hold on, 
hold on, xlabel('X'), ylabel('Y'), zlabel('Z'); hold on,
plot3(PB1(1),PB1(2),PB1(3),'r*'); hold on, text(PB1(1),PB1(2),PB1(3),'PB1'); hold on, 
plot3(PB2(1),PB2(2),PB2(3),'b*'); hold on, text(PB2(1),PB2(2),PB2(3),'PB2'); hold on, 
plot3(PB3(1),PB3(2),PB3(3),'m*'); hold on, text(PB3(1),PB3(2),PB3(3),'PB3'); hold on, 
plot3(PB4(1),PB4(2),PB4(3),'g*'); hold on, text(PB4(1),PB4(2),PB4(3),'PB4');
plot3(PBC(1),PBC(2),PBC(3),'k+'); hold on, text(PBC(1),PBC(2),PBC(3),'PBC');
hold on, 
quiver3(PBC(1),PBC(2),PBC(3), vB24(1)*panel.PhysDimX*10^-3, vB24(2)*panel.PhysDimX*10^-3, vB24(3)*panel.PhysDimX*10^-3);
hold on, quiver3(PBC(1),PBC(2),PBC(3), vB34(1)*panel.PhysDimY*10^-3, vB34(2)*panel.PhysDimY*10^-3, vB34(3)*panel.PhysDimY*10^-3);
hold on, quiver3(PBC(1),PBC(2),PBC(3), vBn(1)*SDDB*10^-3, vBn(2)*SDDB*10^-3, vBn(3)*SDDB*10^-3);
hold on, plot3(TBI(1), TBI(2), TBI(3), 'o', 'MarkerSize',6); hold on, text(TB1(1),TB1(2),TB1(3),'FSB'); hold on, 


%% ... add extra table inclination of 10deg:
TtablePitch = pitch(tablePitch,[0,0,0]);
TB_ = TtablePitch * TB;
srcB_ = TtablePitch*[TBI_,1]'; srcBPos = srcB_(1:3);

%# use if no pitch is to be added:
% TA_=TA;
% fpAPos_ = TA_(1:3,4);
% srcAPos = TAI_';

%% compute rotation angles to use in iMagX:
betaB = asind(TB_(1,3));
alphaB = acosd(TB_(3,3)/cosd(betaB));
gammaB = acosd(TB_(1,1)/cosd(betaB));
gammaB = - gammaB;

geometry.roomName = 'Testbench';
geometry.roomSystem = 'ObliqueSetup';
geometry.roomPlace = 'Casemate';
geometry.date = 'date_of_install';
geometry.time = 'time_of_install';
geometry.axis = 'RADB';
geometry.config = 'free';
geometry.Vec = [0 0 0 0 0 0 0 0 0];
geometry.SAD = SADB;
geometry.AID = AIDB;
geometry.srcPosition = srcBPos;
geometry.detectorT = TB_;
geometry.rx = alphaB;
geometry.ry = betaB;
geometry.rz = gammaB;
geometry.tx = TB_(1,4);
geometry.ty = TB_(2,4);
geometry.tz = TB_(3,4)  

% algorithm config:
config.dof = logical([1, 1, 0, 1, 1, 0, 1, 1, 1]);
%    config.spAttrb (see ******)
config.optimizer = 'cmaes';
config.debugMode = 1;
config.remoteOutliersFromModel = 1;
%# default:
dratio_min  = 0.5; % relative lower tol ratio of sphere diam
dratio_max  = 1.7; % relative upper tol ratio of sphere diam
mf = (geometry.SAD + geometry.AID)/geometry.SAD; % magnification factor
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
spAttrb.threshold = -0.5; % threshold for local binarization
spAttrb.do_filter = 1; % median filter: 0=off, 1=on
spAttrb.ecc = 0.5; % eccentricity
spAttrb.imopen_param = 0.3; % image open parameter (to remove clutter)
spAttrb.detectionMode = 'auto'; % detect spheres: auto, manual
config.spAttrb = spAttrb;

end