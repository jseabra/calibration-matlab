function [geometryOb1, geometryOb2, panelOb1, panelOb2] = roomSetup_testbench_measGEOTOPV2()
clear all
clc
close all

%% Specify flat panel attributes:
panelOb1.Dimension = [1440 1440]; % in PIXELS (X,Y)
panelOb1.Pixel = 0.184; % pixel size (mm)
panelOb1.Phys2RadT = 'rotZ'; % rotation needed to put flat panel in RAD CS
panelOb1.PhysDimX = panelOb1.Dimension(1)/2*panelOb1.Pixel;  
panelOb1.PhysDimY = panelOb1.Dimension(2)/2*panelOb1.Pixel;
panelOb1.invertGray = 1;
T = createTransform2RadFromPanelProps(panelOb1);
panelOb1.Trad = T;
panelOb1.factor = 0; % put zero for no effect

panelOb2.Dimension = [1440 1440];
panelOb2.Pixel = 0.184;
panelOb2.Phys2RadT = 'rot0';
panelOb2.PhysDimX = panelOb2.Dimension(1)/2*panelOb2.Pixel;  
panelOb2.PhysDimY = panelOb2.Dimension(2)/2*panelOb2.Pixel;
panelOb2.invertGray = 1;
T = createTransform2RadFromPanelProps(panelOb2);
panelOb2.Trad = 0;
panelOb2.factor = 0; % put zero for no effect


%% Measurements Geotop:
% points sur surface à la sortie du tube Ob-1
TA1 = [151.1583 48.3074 3.8840];
TA2 = [151.2069 48.3407 3.8848];
TA3 = [151.1772 48.3748 3.8058];
TA4 = [151.1423 48.3503 3.8062];
TAF = [151.2339 48.2554 3.8650];			
TAC = mean([TA1;TA2;TA3;TA4],1);

% points sur surface à la sortie du tube Ob-2
TB1 = [148.8174 48.3123 3.8496];
TB2 = [148.7776 48.3409 3.8494];
TB3 = [148.8058 48.3754 3.7723];
TB4 = [148.8488 48.3441 3.7724];
TBF = [148.7478 48.2554 3.8300];
TBC = mean([TB1;TB2;TB3;TB4],1);

% points sur flat panel A
PA1 = [149.4339 50.4355 5.4043];			
PA2	= [149.5108 50.6224 5.2164];			
PA3 = [149.6985 50.4999 5.5793];			
PA4 = [149.7864 50.6838 5.4069];
PAC = [149.5898 50.5775 5.4110];			
PAC_ = mean([PA1;PA2;PA3;PA4],1);
%PAC-PAC_

% points sur flat panel B
PB1 = [150.5722 50.4319 5.4044];			
PB2	= [150.4896 50.6259 5.2137];			
PB3	= [150.3019 50.4966 5.5826];			
PB4	= [150.2215 50.6858 5.3953];
PBC = [150.4072 50.5784 5.4130];			
PBC_ = mean([PB1;PB2;PB3;PB4],1);
%PBC-PBC_

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
vTABn = cross(vTAB24,vTAB64),
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
TA = [[vA31,0]' [vA12,0]' [vAn,0]' [0 0 0 1]']

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
TB = [[vB24,0]' [vB34,0]' [vBn,0]' [0 0 0 1]']

PAC_ = (PAC - ISO)*10^3 % frs (origin ISO)
PBC_ = (PBC - ISO)*10^3 % frs (origin ISO)

TA(1:3,4) = PAC_';
TB(1:3,4) = PBC_';

% normal vector of tube window A
vAT12 = TA1-TA2; vAT12 = vAT12./norm(vAT12,2);
vAT13 = TA1-TA3; vAT13 = vAT13./norm(vAT13,2);
vATn = cross(vAT13,vAT12)
% normal vector of tube window B
vBT12 = TB1-TB2; vBT12 = vBT12./norm(vBT12,2);
vBT13 = TB1-TB3; vBT13 = vBT13./norm(vBT13,2);
vBTn = cross(vBT12,vBT13)

%# Determine angle between detectorA-tubeA and detectorB-tubeB:
deltaA = acosd(vAn*vATn'/(norm(vAn)*norm(vATn)))
deltaB = acosd(vBn*vBTn'/(norm(vBn)*norm(vBTn)))

% calculate shortest point along vector to the focal spot 
TAI = project_point_to_line_segment(TAC,TAC+vATn,TAF);
TBI = project_point_to_line_segment(TBC,TBC+vBTn,TBF);

TAI_ = (TAI - ISO)*10^3 % frs (origin ISO)
TBI_ = (TBI - ISO)*10^3 % frs (origin ISO)

% SDD, SID, AID:
SDDA = norm(TAI-PAC,2)*10^3 % m to mm
SDDB = norm(TBI-PBC,2)*10^3

SADA = norm(TAI-ISO,2)*10^3
SADB = norm(TBI-ISO,2)*10^3

AIDA = SDDA-SADA;
AIDB = SDDB-SADB;

TAs = TA * trans(0,0,SADA); TAs = TAs(1:3,4)';
TBs = TB * trans(0,0,SADB); TBs = TBs(1:3,4)';

% diff between GEOTOP measured source positions and those obtained with
% simplied model
abs(TAI_-TAs)
abs(TBI_-TBs)


%% plot:
figure,
grid on, hold on,
plot3(TAB1(1),TAB1(2),TAB1(3),'r*'); hold on, text(TAB1(1),TAB1(2),TAB1(3),'TAB1'); hold on, 
plot3(TAB2(1),TAB2(2),TAB2(3),'b*'); hold on, text(TAB2(1),TAB2(2),TAB2(3),'TAB2'); hold on, 
plot3(TAB3(1),TAB3(2),TAB3(3),'m*'); hold on, text(TAB3(1),TAB3(2),TAB3(3),'TAB3'); hold on, 
plot3(TAB4(1),TAB4(2),TAB4(3),'g*'); hold on, text(TAB4(1),TAB4(2),TAB4(3),'TAB4'); hold on, 
plot3(TAB5(1),TAB5(2),TAB5(3),'k*'); hold on, text(TAB5(1),TAB5(2),TAB5(3),'TAB5'); hold on, 
plot3(TAB6(1),TAB6(2),TAB6(3),'c*'); hold on, text(TAB6(1),TAB6(2),TAB6(3),'TAB6'); hold on, 
hold on, xlabel('X'), ylabel('Y'), zlabel('Z');

hold on,
plot3(PA1(1),PA1(2),PA1(3),'r*'); hold on, text(PA1(1),PA1(2),PA1(3),'PA1'); hold on, 
plot3(PA2(1),PA2(2),PA2(3),'b*'); hold on, text(PA2(1),PA2(2),PA2(3),'PA2'); hold on, 
plot3(PA3(1),PA3(2),PA3(3),'m*'); hold on, text(PA3(1),PA3(2),PA3(3),'PA3'); hold on, 
plot3(PA4(1),PA4(2),PA4(3),'g*'); hold on, text(PA4(1),PA4(2),PA4(3),'PA4'); hold on, 
plot3(PAC(1),PAC(2),PAC(3),'k+'); hold on, text(PAC(1),PAC(2),PAC(3),'PAC'); hold on, 
plot3(PB1(1),PB1(2),PB1(3),'r*'); hold on, text(PB1(1),PB1(2),PB1(3),'PB1'); hold on, 
plot3(PB2(1),PB2(2),PB2(3),'b*'); hold on, text(PB2(1),PB2(2),PB2(3),'PB2'); hold on, 
plot3(PB3(1),PB3(2),PB3(3),'m*'); hold on, text(PB3(1),PB3(2),PB3(3),'PB3'); hold on, 
plot3(PB4(1),PB4(2),PB4(3),'g*'); hold on, text(PB4(1),PB4(2),PB4(3),'PB4');
plot3(PBC(1),PBC(2),PBC(3),'k+'); hold on, text(PBC(1),PBC(2),PAC(3),'PBC'); hold on, 

hold on,
quiver3(PAC(1),PAC(2),PAC(3), vA31(1)*panelOb1.PhysDimX*10^-3, vA31(2)*panelOb1.PhysDimX*10^-3, vA31(3)*panelOb1.PhysDimX*10^-3);
hold on, quiver3(PAC(1),PAC(2),PAC(3), vA12(1)*panelOb1.PhysDimY*10^-3, vA12(2)*panelOb1.PhysDimY*10^-3, vA12(3)*panelOb1.PhysDimY*10^-3);
hold on, quiver3(PAC(1),PAC(2),PAC(3), vAn(1)*SDDA*10^-3, vAn(2)*SDDA*10^-3, vAn(3)*SDDA*10^-3);
hold on,
quiver3(PBC(1),PBC(2),PBC(3), vB24(1)*panelOb2.PhysDimX*10^-3, vB24(2)*panelOb2.PhysDimX*10^-3, vB24(3)*panelOb2.PhysDimX*10^-3);
hold on, quiver3(PBC(1),PBC(2),PBC(3), vB34(1)*panelOb2.PhysDimY*10^-3, vB34(2)*panelOb2.PhysDimY*10^-3, vB34(3)*panelOb2.PhysDimY*10^-3);
hold on, quiver3(PBC(1),PBC(2),PBC(3), vBn(1)*SDDB*10^-3, vBn(2)*SDDB*10^-3, vBn(3)*SDDB*10^-3);

% hold on,
% hold on, quiver3(TAC(1),TAC(2),TAC(3), vATn(1)*SDDA*10^-3, vATn(2)*SDDA*10^-3, vATn(3)*SDDA*10^-3);
% hold on, plot3(TAF(1), TAF(2), TAF(3), '*', 'MarkerSize',6);
% hold on, quiver3(TBC(1),TBC(2),TBC(3), vBTn(1)*SDDB*10^-3, vBTn(2)*SDDB*10^-3, vBTn(3)*SDDB*10^-3);
% hold on, xlabel('X'), ylabel('Y'), zlabel('Z');
% hold on, plot3(TBF(1), TBF(2), TBF(3), '*', 'MarkerSize',6);
% 
hold on, plot3(TAI(1), TAI(2), TAI(3), 'o', 'MarkerSize',6); hold on, text(TA1(1),TA1(2),TA1(3),'FSA'); hold on, 
hold on, plot3(TBI(1), TBI(2), TBI(3), 'o', 'MarkerSize',6); hold on, text(TB1(1),TB1(2),TB1(3),'FSB'); hold on, 
% 

%% ... add extra table inclination of 10deg:
TtablePitch = pitch(tablePitch,[0,0,0]);
TA_ = TtablePitch * TA;
TB_ = TtablePitch * TB;

srcA_ = TtablePitch*[TAI_,1]'; srcAPos = srcA_(1:3);
srcB_ = TtablePitch*[TBI_,1]'; srcBPos = srcB_(1:3);

%# use if no pitch is to be added:
% TA_=TA;
% TB_=TB;
% fpAPos_ = TA_(1:3,4);
% fpBPos_ = TB_(1:3,4);
% srcAPos = TAI_';
% srcBPos = TBI_';


%% compute rotation angles to use in iMagX:
betaA = asind(TA_(1,3));
alphaA = acosd(TA_(3,3)/cosd(betaA));
gammaA = acosd(TA_(1,1)/cosd(betaA));

betaB = asind(TB_(1,3));
alphaB = acosd(TB_(3,3)/cosd(betaB));
gammaB = acosd(TB_(1,1)/cosd(betaB));
gammaB = - gammaB;


%% generate 'geometry' outputs:
% output ob-1:
geometryOb1.axis = 'OB-1';               % axis to calibrate
geometryOb1.config = 'free';             % 'gantry','free'
geometryOb1.Vec = [0 0 0 0 0 0 0 0 0];   % (Sx, Sy, Sz, Dx, Dy, Dz, Rx, Ry, Rz)
%geometryOb1.srcT = srcA_;                % 4x4 source transformation matrix into FRS
geometryOb1.srcPosition = srcAPos;
geometryOb1.detectorT = TA_; 
geometryOb1.SAD = SADA;                 
geometryOb1.AID = AIDA;                 
geometryOb1.rx = alphaA; 
geometryOb1.ry = betaA;
geometryOb1.rz = gammaA;
geometryOb1.tx = 0; 
geometryOb1.ty = 0;
geometryOb1.tz = 0;

% output ob-2:
geometryOb2.axis = 'OB-2';               % axis to calibrate
geometryOb2.config = 'free';             % 'gantry','free'
geometryOb2.Vec = [0 0 0 0 0 0 0 0 0];   % (Sx, Sy, Sz, Dx, Dy, Dz, Rx, Ry, Rz)
%geometryOb2.srcT = srcB_;                % 4x4 source transformation matrix into FRS
geometryOb2.srcPosition = srcBPos;
geometryOb2.detectorT = TB_; 
geometryOb2.SAD = SADB;
geometryOb2.AID = AIDB;
geometryOb2.rx = alphaB; 
geometryOb2.ry = betaB;
geometryOb2.rz = gammaB;
geometryOb2.tx = 0; 
geometryOb2.ty = 0;
geometryOb2.tz = 0;
