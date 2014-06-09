function [geometryOb1, geometryOb2, panelOb1, panelOb2] = roomSetup_shreveport()
%ROOMSETUP_SHREVEPORT retrieves the geometry and flat panel structures
%needed to perform calibration.
% Oblique geometry only for shreveport installation.
% Room setup was calculated based on the IBA specs and mechanical drawings
clc
clear
close all

%% Specify flat panel attributes:
panelOb1.Dimension = [1440 1440]; % in PIXELS (X,Y)
panelOb1.Pixel = 0.184; % pixel size (mm)
panelOb1.Phys2RadT = 'rotZ'; % rotation needed to put flat panel in RAD CS
panelOb1.PhysDimX = panelOb1.Dimension(1)/2*panelOb1.Pixel;  
panelOb1.PhysDimY = panelOb1.Dimension(2)/2*panelOb1.Pixel;
T = createTransform2RadFromPanelProps(panelOb1);
panelOb1.Trad = T;
panelOb1.invertGray=1;
panelOb1.factorLow = 20000; % non-normalized 
panelOb1.factorHigh = 0; % non-normalized

panelOb2.Dimension = [1440 1440];
panelOb2.Pixel = 0.184;
panelOb2.Phys2RadT = 'rot0';
panelOb2.PhysDimX = panelOb2.Dimension(1)/2*panelOb2.Pixel;  
panelOb2.PhysDimY = panelOb2.Dimension(2)/2*panelOb2.Pixel;
T = createTransform2RadFromPanelProps(panelOb2);
panelOb2.Trad = T;
panelOb2.invertGray=1;
panelOb2.factorLow = 20000; % non-normalized 
panelOb2.factorHigh = 0; % non-normalized

SAD_Ob1 = 2500; SAD_Ob2 = SAD_Ob1; 
AID_Ob1 = 820; AID_Ob2 = AID_Ob1;

% Points in FP Ob-2 plane FRS: 
p2H = [293.009, 409.426, 662.958];
p2V = [448.204, 296.864, 632.813];
p2C = [410, 410, 579.83];

% Center of FP Ob-1 FRS:
p1C = [-410, 410, 579.83];

v2H = p2C-p2H; v2Hn = v2H./sqrt(sum(v2H.^2));
v2V = p2V-p2C; v2Vn = v2V./sqrt(sum(v2V.^2));

v1H = -v2H; v1H(1) = -v1H(1);
v1V = v2V; v1V(1) = -v1V(1);
v1Hn = -v2Hn; v1Hn(1) = -v1Hn(1);
v1Vn = v2Vn; v1Vn(1) = -v1Vn(1);

% FP normal vectors:
vn1 = cross(v1Hn,v1Vn);
vn2 = cross(v2Hn,v2Vn);

% Transformation matrices:
Tob1 = [[v1Hn,0]' [v1Vn,0]' [vn1,0]' [0 0 0 1]'];
Tob2 = [[v2Hn,0]' [v2Vn,0]' [vn2,0]' [0 0 0 1]'];

FPTob1 = Tob1*trans(0,0,-AID_Ob1)
SRCTob1 = Tob1*trans(0,0,SAD_Ob1)

FPTob2 = Tob2*trans(0,0,-AID_Ob2)
SRCTob2 = Tob2*trans(0,0,SAD_Ob2)


figure, 
quiver3(p2C(1),p2C(2),p2C(3),v2H(1),v2H(2),v2H(3)); hold on,
quiver3(p2C(1),p2C(2),p2C(3),v2V(1),v2V(2),v2V(3)); hold on,
quiver3(p2C(1),p2C(2),p2C(3),vn2(1).*100,vn2(2).*100,vn2(3).*100); hold on,

quiver3(p1C(1),p1C(2),p1C(3),v1H(1),v1H(2),v1H(3)); hold on,
quiver3(p1C(1),p1C(2),p1C(3),v1V(1),v1V(2),v1V(3)); hold on,
quiver3(p1C(1),p1C(2),p1C(3),vn1(1).*100,vn1(2).*100,vn1(3).*100); hold on,

xlabel('X'); ylabel('Y'); zlabel('Z'); 
hold on,
box on;

% flat panel dimensions:
Dimension = [1440 1440]; % in PIXELS (X,Y)
Pixel = 0.184; % pixel size (mm)
dimx = Dimension(1)/2*Pixel;  
dimy = Dimension(2)/2*Pixel;
% fp corners in local CS:
bl = [-dimx, -dimy, 0, 1];
br = [dimx, -dimy, 0, 1];
tl = [-dimx, dimy, 0, 1];
tr = [dimx, dimy, 0, 1];

detectorPosblOb2 = FPTob2 * bl';
detectorPosbrOb2 = FPTob2 * br';
detectorPostlOb2 = FPTob2 * tl';
detectorPostrOb2 = FPTob2 * tr';

hold on,
plot3([detectorPosblOb2(1),detectorPosbrOb2(1)], [detectorPosblOb2(2),detectorPosbrOb2(2)], [detectorPosblOb2(3), detectorPosbrOb2(3)],'k', 'MarkerSize', 2, 'LineWidth',2);hold on,
plot3([detectorPostlOb2(1),detectorPostrOb2(1)], [detectorPostlOb2(2),detectorPostrOb2(2)], [detectorPostlOb2(3), detectorPostrOb2(3)],'k', 'MarkerSize', 2, 'LineWidth',2);hold on,
plot3([detectorPosblOb2(1),detectorPostlOb2(1)], [detectorPosblOb2(2),detectorPostlOb2(2)], [detectorPosblOb2(3), detectorPostlOb2(3)],'k', 'MarkerSize', 2, 'LineWidth',2);hold on,
plot3([detectorPosbrOb2(1),detectorPostrOb2(1)], [detectorPosbrOb2(2),detectorPostrOb2(2)], [detectorPosbrOb2(3), detectorPostrOb2(3)],'k', 'MarkerSize', 2, 'LineWidth',2);hold on,

detectorPosblOb1 = FPTob1 * bl';
detectorPosbrOb1 = FPTob1 * br';
detectorPostlOb1 = FPTob1 * tl';
detectorPostrOb1 = FPTob1 * tr';

hold on,
plot3([detectorPosblOb1(1),detectorPosbrOb1(1)], [detectorPosblOb1(2),detectorPosbrOb1(2)], [detectorPosblOb1(3), detectorPosbrOb1(3)],'k', 'MarkerSize', 2, 'LineWidth',2);hold on,
plot3([detectorPostlOb1(1),detectorPostrOb1(1)], [detectorPostlOb1(2),detectorPostrOb1(2)], [detectorPostlOb1(3), detectorPostrOb1(3)],'k', 'MarkerSize', 2, 'LineWidth',2);hold on,
plot3([detectorPosblOb1(1),detectorPostlOb1(1)], [detectorPosblOb1(2),detectorPostlOb1(2)], [detectorPosblOb1(3), detectorPostlOb1(3)],'k', 'MarkerSize', 2, 'LineWidth',2);hold on,
plot3([detectorPosbrOb1(1),detectorPostrOb1(1)], [detectorPosbrOb1(2),detectorPostrOb1(2)], [detectorPosbrOb1(3), detectorPostrOb1(3)],'k', 'MarkerSize', 2, 'LineWidth',2);hold on,

% Source positions FRS:
src1 = vn1.*SAD_Ob1
src2 = vn2.*SAD_Ob2

vb1 = p1C - src1; vb1n = vb1./sqrt(sum(vb1.^2)); sqrt(sum(vb1.^2))
vb2 = p2C - src2; vb2n = vb2./sqrt(sum(vb2.^2)); sqrt(sum(vb2.^2))

% test orthogonality: dot product
sum( vb2n.*v2Vn )
sum( vb2n.*v2Hn )

sum( vb1n.*v1Vn )
sum( vb1n.*v1Hn )

alpha = 144.7360;
beta = 30;
gamma = 19.7270

% create geometry outputs:
geometryOb1.axis = 'OB-1';               % axis to calibrate
geometryOb1.config = 'free';             % 'gantry','free'
geometryOb1.Vec = [0 0 0 0 0 0 0 0 0];   % (Sx, Sy, Sz, Dx, Dy, Dz, Rx, Ry, Rz)
geometryOb1.srcT = SRCTob1;              % 4x4 source transformation matrix into FRS
geometryOb1.detectorT = FPTob1; 
geometryOb1.SAD = SAD_Ob1;               % rotation around Zfrs  
geometryOb1.AID = AID_Ob1;               % rotation around Zfrs  
geometryOb1.rx = alpha; 
geometryOb1.ry = beta;
geometryOb1.rz = gamma;
geometryOb1.tx = 0;
geometryOb1.ty = 0;
geometryOb1.tz = 0;

geometryOb2.axis = 'OB-2';               % axis to calibrate
geometryOb2.config = 'free';             % 'gantry','free'
geometryOb2.Vec = [0 0 0 0 0 0 0 0 0];   % (Sx, Sy, Sz, Dx, Dy, Dz, Rx, Ry, Rz)
geometryOb2.srcT = SRCTob2;              % 4x4 source transformation matrix into FRS
geometryOb2.detectorT = FPTob2; 
geometryOb2.SAD = SAD_Ob2;               % rotation around Zfrs  
geometryOb2.AID = AID_Ob2;               % rotation around Zfrs  
geometryOb2.rx = alpha; 
geometryOb2.ry = -beta;
geometryOb2.rz = -gamma;
geometryOb2.tx = 0;
geometryOb2.ty = 0;
geometryOb2.tz = 0;


end