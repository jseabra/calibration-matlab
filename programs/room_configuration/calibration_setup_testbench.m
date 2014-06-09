function configuration_testbench_oblique_ob1();


panel.Dimension = [1440 1440];
panel.Pixel = 0.184;
panel.Phys2RadT = 'rot0';
panel.PhysDimX = panelOb2.Dimension(1)/2*panelOb2.Pixel;
panel.PhysDimY = panelOb2.Dimension(2)/2*panelOb2.Pixel;
panel.invertGray = 1;
panel.factorLow = 0;
panel.factorHigh = 1;

%% Measurements Geotop:
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
PAC_ = mean([PA1;PA2;PA3;PA4],1);
%PAC-PAC_
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
plot3(PB1(1),PB1(2),PB1(3),'r*'); hold on, text(PB1(1),PB1(2),PB1(3),'PB1'); hold on, 
plot3(PB2(1),PB2(2),PB2(3),'b*'); hold on, text(PB2(1),PB2(2),PB2(3),'PB2'); hold on, 
plot3(PB3(1),PB3(2),PB3(3),'m*'); hold on, text(PB3(1),PB3(2),PB3(3),'PB3'); hold on, 
plot3(PB4(1),PB4(2),PB4(3),'g*'); hold on, text(PB4(1),PB4(2),PB4(3),'PB4');
plot3(PBC(1),PBC(2),PBC(3),'k+'); hold on, text(PBC(1),PBC(2),PAC(3),'PBC'); hold on, 
hold on,
quiver3(PBC(1),PBC(2),PBC(3), vB24(1)*panelOb2.PhysDimX*10^-3, vB24(2)*panelOb2.PhysDimX*10^-3, vB24(3)*panelOb2.PhysDimX*10^-3);
hold on, quiver3(PBC(1),PBC(2),PBC(3), vB34(1)*panelOb2.PhysDimY*10^-3, vB34(2)*panelOb2.PhysDimY*10^-3, vB34(3)*panelOb2.PhysDimY*10^-3);
hold on, quiver3(PBC(1),PBC(2),PBC(3), vBn(1)*SDDB*10^-3, vBn(2)*SDDB*10^-3, vBn(3)*SDDB*10^-3);
hold on, plot3(TBI(1), TBI(2), TBI(3), 'o', 'MarkerSize',6); hold on, text(TB1(1),TB1(2),TB1(3),'FSB'); hold on, 

% ... add extra table inclination of 10deg:
TtablePitch = pitch(tablePitch,[0,0,0]);
TB_ = TtablePitch * TB;

srcB_ = TtablePitch*[TBI_,1]'; srcBPos = srcB_(1:3);

% compute rotation angles to use in iMagX:
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



end

% Inputs:
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
%    config.optimizer = <'cmaes',...>
%    config.debugMode = <0,1>
%    config.fitModel = [Sx, Sy, Sz, Dx, Dy, Dz, Rx, Ry, Rz]; %  <0=computes
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








% flat panel attributes:
panel.Dimension = [2880 2881]; % in PIXELS (X,Y)
panel.Pixel = 0.148; % pixel size [mm]
panel.Phys2RadT = 'rotX'; % rotation needed to put flat panel in RAD CS
panel.PhysDimX = panel.Dimension(1)/2*panel.Pixel;  
panel.PhysDimY = panel.Dimension(2)/2*panel.Pixel;
panel.invertGray = 1;

% geometry:
geometry.axis = 'RADB'; % axis to calibrate
geometry.config = 'gantry'; % 'gantry','free'
geometry.Vec = [0 0 100 0 0 0 0 0 0];   % (Sx, Sy, Sz, Dx, Dy, Dz, Rx, Ry, Rz)
geometry.gantryAngleOffset = -90; % offset projection to gantry angle
geometry.SAD = 2834; % source to axis distance [mm]                 
geometry.AID = 559; % axis to detector distance [mm]                 
geometry.gantryAngle = 0; % = no longer used (read from image header)

% cylinder phantom:
frsToPhanTransform = pitch(-90,[0,0,0]); %disp('Phantom C3 (3spheres) toward +Yfrs');
cylinderModel = 1; %=1 casemate, 0 medcom
phantom = getCylinderPhantomPositions(frsToPhanTransform, cylinderModel);
objOffset = [0,0,0]; % isocenter offset (x,y,z) [mm]

% optimization:
optimization.dof = logical([1 1 0 1 1 0 0 0 1]);
optimization.method = 'cmaes'; % 'simplex', 'cmaes', 'lsqnonlin', 'levmar', 'cmb'      

% data management:
data.inputDir = 'C:\imagx_data\cylinder-static-11042014'; % root directory with calibration data
data.imagesDir = [data.inputDir,'\drs']; % directory where images were stored
data.imagesRegExp = 'RADB'; % regular expression to identify images in directory
data.sortBy = 'date'; % 'date'or 'none'
data.anglesFromHeaders = 1;
data.tempDir = [data.inputDir,'\tmp-noRot'];
% define if data.anglesFromHeaders=0, otherwise set to ""
data.anglesFname = [data.inputDir,'\gantryAngles_CW.txt']; % file must contain an angle for each image

end