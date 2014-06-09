function [geometry, panel, phantom, data_mngt, config] = configuration_testbench_oblique_ob1()
% This is a configuration for oblique (free) room geometry.
% Copy onto new file and edit to create new room geometry.
clc
clear
addpath(genpath('../../utils/')); % import utils

% panel:
panel.Name = 'THALES2630';
panel.Dimension = [1440 1440];
panel.Pixel = 0.184;
panel.Phys2RadT = 'rotationZ';
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
data_mngt.input = 'C:\imagx_data\testbench-calibration-oblique\cylinder\RADA Apr 10 2014 160925 kVp 80 mA 100 ms 80.raw';
data_mngt.output = 'C:\imagx_data\testbench-calibration-oblique\results\Ob1\'

%% =====  Measurements Geotop
% points sur surface à la sortie du tube Ob-1
TA1 = [151.1583 48.3074 3.8840];
TA2 = [151.2069 48.3407 3.8848];
TA3 = [151.1772 48.3748 3.8058];
TA4 = [151.1423 48.3503 3.8062];
TAF = [151.2339 48.2554 3.8650];			
TAC = mean([TA1;TA2;TA3;TA4],1);
% points sur flat panel A
PA1 = [149.4339 50.4355 5.4043];			
PA2	= [149.5108 50.6224 5.2164];			
PA3 = [149.6985 50.4999 5.5793];			
PA4 = [149.7864 50.6838 5.4069];
PAC = [149.5898 50.5775 5.4110];			
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
% direction vectors flat panel A
% direction vector y
vA12 = PA1-PA2;
vA12 = vA12./norm(vA12,2);
% direction vector x
vA31 = PA3-PA1;
vA31 = vA31./norm(vA31,2);
% direction vector z
vAn = cross(vA31,vA12);
% check angle betweetn two vectors
%acosd((vA31*vA12')/(norm(vA31,2)*norm(vA12,2)))
% plane A transformation matrix:
TA = [[vA31,0]' [vA12,0]' [vAn,0]' [0 0 0 1]'];

PAC_ = (PAC - ISO)*10^3 % frs (origin ISO)

TA(1:3,4) = PAC_';

% normal vector of tube window A
vAT12 = TA1-TA2; vAT12 = vAT12./norm(vAT12,2);
vAT13 = TA1-TA3; vAT13 = vAT13./norm(vAT13,2);
vATn = cross(vAT13,vAT12);

%# Determine angle between detectorA-tubeA:
deltaA = acosd(vAn*vATn'/(norm(vAn)*norm(vATn)));

% calculate shortest point along vector to the focal spot:
TAI = project_point_to_line_segment(TAC,TAC+vATn,TAF);
TAI_ = (TAI - ISO)*10^3; % frs (origin ISO)

% SDD, SID, AID:
SDDA = norm(TAI-PAC,2)*10^3; % m to mm
SADA = norm(TAI-ISO,2)*10^3;
AIDA = SDDA-SADA;

TAs = TA * trans(0,0,SADA); TAs = TAs(1:3,4)';

% diff between GEOTOP measured source positions and those obtained with
% simplied model
abs(TAI_-TAs);

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
plot3(PA1(1),PA1(2),PA1(3),'r*'); hold on, text(PA1(1),PA1(2),PA1(3),'PA1'); hold on, 
plot3(PA2(1),PA2(2),PA2(3),'b*'); hold on, text(PA2(1),PA2(2),PA2(3),'PA2'); hold on, 
plot3(PA3(1),PA3(2),PA3(3),'m*'); hold on, text(PA3(1),PA3(2),PA3(3),'PA3'); hold on, 
plot3(PA4(1),PA4(2),PA4(3),'g*'); hold on, text(PA4(1),PA4(2),PA4(3),'PA4'); hold on, 
plot3(PAC(1),PAC(2),PAC(3),'k+'); hold on, text(PAC(1),PAC(2),PAC(3),'PAC'); hold on, 
quiver3(PAC(1),PAC(2),PAC(3), vA31(1)*panel.PhysDimX*10^-3, vA31(2)*panel.PhysDimX*10^-3, vA31(3)*panel.PhysDimX*10^-3); hold on,
quiver3(PAC(1),PAC(2),PAC(3), vA12(1)*panel.PhysDimY*10^-3, vA12(2)*panel.PhysDimY*10^-3, vA12(3)*panel.PhysDimY*10^-3); hold on, 
quiver3(PAC(1),PAC(2),PAC(3), vAn(1)*SDDA*10^-3, vAn(2)*SDDA*10^-3, vAn(3)*SDDA*10^-3);

%% ... add extra table inclination of 10deg:
TtablePitch = pitch(tablePitch,[0,0,0]);
TA_ = TtablePitch * TA;
srcA_ = TtablePitch*[TAI_,1]'; srcAPos = srcA_(1:3);

%# use if no pitch is to be added:
% TA_=TA;
% fpAPos_ = TA_(1:3,4);
% srcAPos = TAI_';

%% compute rotation angles to use in iMagX:
betaA = asind(TA_(1,3));
alphaA = acosd(TA_(3,3)/cosd(betaA));
gammaA = acosd(TA_(1,1)/cosd(betaA));

geometry.roomName = 'Testbench';
geometry.roomSystem = 'ObliqueSetup';
geometry.roomPlace = 'Casemate';
geometry.date = 112233;
geometry.time = 112233;
geometry.axis = 'RADA';
geometry.config = 'free';
geometry.Vec = [0 0 0 0 0 0 0 0 0];
geometry.SAD = SADA;
geometry.AID = AIDA;
geometry.srcPosition = srcAPos;
geometry.detectorT = TA_;
geometry.rx = alphaA;
geometry.ry = betaA;
geometry.rz = gammaA;
geometry.tx = TA_(1,4);
geometry.ty = TA_(2,4);
geometry.tz = TA_(3,4)  

% algorithm config:
config.dof = logical([1, 1, 0, 1, 1, 0, 1, 1, 1]);
%    config.spAttrb (see ******)
config.optimizer = 'cmaes';
config.debugMode = 1;
config.remoteOutliersFromModel = 1;
%# default:
dratio_min  = 0.7; % relative lower tol ratio of sphere diam
dratio_max  = 1.5; % relative upper tol ratio of sphere diam
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
spAttrb.threshold = 0; % threshold for local binarization
spAttrb.do_filter = 1; % median filter: 0=off, 1=on
spAttrb.ecc = 0.5; % eccentricity
spAttrb.imopen_param = 0.3; % image open parameter (to remove clutter)
spAttrb.detectionMode = 'auto'; % detect spheres: auto, manual
config.spAttrb = spAttrb;

end