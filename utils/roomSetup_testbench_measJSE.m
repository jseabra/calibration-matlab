function [geometryOb1, geometryOb2, panelOb1, panelOb2] = roomSetup_testbench_measJSE()
%ROOMSETUP_TESTBENCH retrieves the geometry and flat panel structures
%needed to perform calibration.
% Oblique geometry only.
% Based on measurements made in 10/04/2014
clc
clear
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

panelOb2.Dimension = [1440 1440];
panelOb2.Pixel = 0.184;
panelOb2.Phys2RadT = 'rot0';
panelOb2.PhysDimX = panelOb2.Dimension(1)/2*panelOb2.Pixel;  
panelOb2.PhysDimY = panelOb2.Dimension(2)/2*panelOb2.Pixel;
panelOb2.invertGray = 1;
T = createTransform2RadFromPanelProps(panelOb2);
panelOb2.Trad = T;

%% OBLIQUE 1
% flat panel positions in "casemate" CS (origin is back left corner, when
% facing the wall from the corridor)
p11 = [4640, -1103, 2246];
p21 = [4725, -887, 2041];
p22 = [5021, -814, 2240];

v1121 = p11 - p21; % vertical vector
v2221 = p22 - p21; % horizontal vector
% check angle betweetn two vectors (they should be perpendicular)
acosd((v1121*v2221')/(sqrt(sum(v1121.^2))*sqrt(sum(v2221.^2))));

% centroid
pC = p21+v1121*.5+v2221*.5;

% direction vectors:
v1121 = v1121./sqrt(sum(v1121.^2));
v2221 = v2221./sqrt(sum(v2221.^2));
vN = cross(v2221,v1121);
Tob1 = [[v2221,0]' [v1121,0]' [vN,0]' [0 0 0 1]'];

% source position:
srcPos = [6500,-3242,612];
sdd = 3209+25+80;
sddp = 2860; % projection of beam axis on the floor
%acosd(sddp/sdd)
figure(1),
hold on, quiver3(pC(1),pC(2),pC(3), v2221(1)*panelOb1.PhysDimX, v2221(2)*panelOb1.PhysDimX, v2221(3));
hold on, quiver3(pC(1),pC(2),pC(3), v1121(1)*panelOb1.PhysDimY, v1121(2)*panelOb1.PhysDimY, v1121(3));
hold on, quiver3(pC(1),pC(2),pC(3), vN(1)*sdd, vN(2)*sdd, vN(3)*sdd);
hold on, xlabel('X'), ylabel('Y'), zlabel('Z');
plot3([srcPos(1),-vN(1)*sdd+srcPos(1)],[srcPos(2),-vN(2)*sdd+srcPos(2)], [srcPos(3),-vN(3)*sdd+srcPos(3)], 'k');

% pC(1)+vN(1)*sdd - srcPos(1)
% pC(2)+vN(2)*sdd - srcPos(2)
% pC(3)+vN(3)*sdd - srcPos(3)


%% OBLIQUE 2
% flat panel positions in "casemate" CS (origin is back left corner, when
% facing the wall from the corridor)
p12 = [5853, -1017, 2250];
p21 = [5452, -792,  2230];
p22 = [5757, -871,  2037];

v2221 = p22 - p21;
v1222 = p12 - p22;
% check angle betweetn two vectors (they should be perpendicular)
acosd((v2221*v1222')/(sqrt(sum(v2221.^2))*sqrt(sum(v1222.^2))));

% centroid
pC = p21+v1222*.5+v2221*.5;

% direction vectors:
v2221 = v2221./sqrt(sum(v2221.^2));
v1222 = v1222./sqrt(sum(v1222.^2));
%vN = cross(v2221,v1222);

% USE OB-2 PANEL SYMMETRIC TO OB-1 --------------------
p11 = [4640, -1103, 2246];
p21 = [4725, -887, 2041];
p22 = [5021, -814, 2240];

v1121 = p11 - p21; % vertical vector
v2221 = p22 - p21; % horizontal vector

v1121 = v1121./sqrt(sum(v1121.^2));
v2221 = v2221./sqrt(sum(v2221.^2));

v1121(1)=-v1121(1);
v2221(2)=-v2221(2);
v2221(3)=-v2221(3);
vN = cross(v2221,v1121);
%------------------------------------------------------

Tob2 = [[v2221,0]' [v1121,0]' [vN,0]' [0 0 0 1]'];

% source position:
srcPos = [4020,-3293,584];
sdd = sqrt(sum((srcPos - pC).^2));
sdd = 3209+25+80;

sddp = 2860; % projection of beam axis on the floor
%acosd(sddp/sdd)
figure(1), 
hold on, quiver3(pC(1),pC(2),pC(3), v2221(1)*panelOb1.PhysDimX, v2221(2)*panelOb1.PhysDimX, v2221(3)); 
hold on, quiver3(pC(1),pC(2),pC(3), v1222(1)*panelOb1.PhysDimY, v1222(2)*panelOb1.PhysDimY, v1222(3)); 
hold on, quiver3(pC(1),pC(2),pC(3), vN(1)*sdd, vN(2)*sdd, vN(3)*sdd); 
hold on, xlabel('X'), ylabel('Y'), zlabel('Z');
plot3([srcPos(1),-vN(1)*sdd+srcPos(1)],[srcPos(2),-vN(2)*sdd+srcPos(2)], [srcPos(3),-vN(3)*sdd+srcPos(3)], 'k');
box on, legend('X1', 'Y1', 'FPtoSrc1', 'SrctoFP1', 'X2', 'Y2', 'FPtoSrc2', 'SrctoFP2')

% pC(1)+vN(1)*sdd - srcPos(1)
% pC(2)+vN(2)*sdd - srcPos(2)
% pC(3)+vN(3)*sdd - srcPos(3)


%% Compute flat panel and tube positions...
aid1 = -820;
sid1 = 2500;

% transformations
src1 = Tob1 * trans(0,0,sid1);
fp1 = Tob1 * trans(0,0,aid1);
% position vectors
src1Pos = src1(1:3,4);
fp1Pos = fp1(1:3,4);

aid2 = -820;
sid2 = sdd+aid2;

% transformations
src2 = Tob2 * trans(0,0,sid2);
fp2 = Tob2 * trans(0,0,aid2);
% position vectors
src2Pos = src2(1:3,4);
fp2Pos = fp2(1:3,4);


%% ... add extra table inclination of 10deg:
tablePitch = 10; 
TtablePitch = pitch(tablePitch,[0,0,0]);
Tob1_ = TtablePitch * Tob1;
Tob2_ = TtablePitch * Tob2;

% transformations
src1r = Tob1_ * trans(0,0,sid1);
fp1r = Tob1_ * trans(0,0,aid1);
% position vectors
src1Posr = src1r(1:3,4);
fp1Posr = fp1r(1:3,4);

% transformations
src2r = Tob2_ * trans(0,0,sid2);
fp2r = Tob2_ * trans(0,0,aid2);
% position vectors
src2Posr = src2r(1:3,4);
fp2Posr = fp2r(1:3,4);


%% compute rotation angles to use in iMagX:
betaOb1 = asind(Tob1_(1,3));
alphaOb1 = acosd(Tob1_(3,3)/cosd(betaOb1));
gammaOb1 = acosd(Tob1_(1,1)/cosd(betaOb1));

betaOb2 = asind(Tob2_(1,3));
alphaOb2 = acosd(Tob2_(3,3)/cosd(betaOb2));
gammaOb2 = acosd(Tob2_(1,1)/cosd(betaOb2));
gammaOb2 = - gammaOb2;


%% generate 'geometry' outputs:
% output ob-1:
geometryOb1.axis = 'OB-1';               % axis to calibrate
geometryOb1.config = 'free';             % 'gantry','free'
geometryOb1.Vec = [0 0 0 0 0 0 0 0 0];   % (Sx, Sy, Sz, Dx, Dy, Dz, Rx, Ry, Rz)
geometryOb1.srcT = src1r                 % 4x4 source transformation matrix into FRS
geometryOb1.detectorT = fp1r 
geometryOb1.SAD = sid1;                 
geometryOb1.AID = -aid1;                 
geometryOb1.rx = alphaOb1; 
geometryOb1.ry = betaOb1;
geometryOb1.rz = gammaOb1;
geometryOb1.tx = 0; 
geometryOb1.ty = 0;
geometryOb1.tz = 0;

% output ob-2:
geometryOb2.axis = 'OB-2';               % axis to calibrate
geometryOb2.config = 'free';             % 'gantry','free'
geometryOb2.Vec = [0 0 0 0 0 0 0 0 0];   % (Sx, Sy, Sz, Dx, Dy, Dz, Rx, Ry, Rz)
geometryOb2.srcT = src2r                 % 4x4 source transformation matrix into FRS
geometryOb2.detectorT = fp2r 
geometryOb2.SAD = sid2;                  
geometryOb2.AID = -aid2;                   
geometryOb2.rx = alphaOb2; 
geometryOb2.ry = betaOb2;
geometryOb2.rz = gammaOb2;
geometryOb2.tx = 0; 
geometryOb2.ty = 0;
geometryOb2.tz = 0;


%% Now plot this geometry:
figure(100),
grid on, hold on,
iso = [0, 0, 0];
u = [500, 0, 0, 1];
v = [0, 500, 0, 1];
w = [0, 0, 500, 1];
ur = TtablePitch * u';
vr = TtablePitch * v';
wr = TtablePitch * w';

% FRS CS:
plot3([iso(1),u(1)],[iso(2),u(2)],[iso(3),u(3)],'k','LineWidth',2); hold on;
plot3([iso(1),v(1)],[iso(2),v(2)],[iso(3),v(3)],'k','LineWidth',2); hold on;
plot3([iso(1),w(1)],[iso(2),w(2)],[iso(3),w(3)],'k','LineWidth',2); hold on;
% Rotated CS:
plot3([iso(1),ur(1)],[iso(2),ur(2)],[iso(3),ur(3)],'m'); hold on;
plot3([iso(1),vr(1)],[iso(2),vr(2)],[iso(3),vr(3)],'m'); hold on;
plot3([iso(1),wr(1)],[iso(2),wr(2)],[iso(3),wr(3)],'m'); hold on;

% fp dimensions [mm]
dimx = 143.5;
dimy = 130.5;
bl = [-dimx, -dimy, 0, 1];
br = [dimx, -dimy, 0, 1];
tl = [-dimx, dimy, 0, 1];
tr = [dimx, dimy, 0, 1];

fpbl1 = fp1 * bl';
fpbr1 = fp1 * br';
fptl1 = fp1 * tl';
fptr1 = fp1 * tr';
fpbl2 = fp2 * bl';
fpbr2 = fp2 * br';
fptl2 = fp2 * tl';
fptr2 = fp2 * tr';

fpbl1r = fp1r * bl';
fpbr1r = fp1r * br';
fptl1r = fp1r * tl';
fptr1r = fp1r * tr';
fpbl2r = fp2r * bl';
fpbr2r = fp2r * br';
fptl2r = fp2r * tl';
fptr2r = fp2r * tr';

% tube:
plot3(src1Pos(1), src1Pos(2), src1Pos(3), 'b*'); hold on;
plot3(src1Posr(1), src1Posr(2), src1Posr(3), 'r*'); hold on;
plot3(src2Pos(1), src2Pos(2), src2Pos(3), 'b*'); hold on;
plot3(src2Posr(1), src2Posr(2), src2Posr(3), 'r*'); hold on;

% flat panel:
plot3(fp1Pos(1), fp1Pos(2), fp1Pos(3),'bd', 'MarkerSize', 4, 'LineWidth',2); hold on,
plot3([fpbl1(1),fpbr1(1)], [fpbl1(2),fpbr1(2)], [fpbl1(3), fpbr1(3)],'b', 'MarkerSize', 2, 'LineWidth',2);hold on,
plot3([fptl1(1),fptr1(1)], [fptl1(2),fptr1(2)], [fptl1(3), fptr1(3)],'b', 'MarkerSize', 2, 'LineWidth',2);hold on,
plot3([fpbl1(1),fptl1(1)], [fpbl1(2),fptl1(2)], [fpbl1(3), fptl1(3)],'b', 'MarkerSize', 2, 'LineWidth',2);hold on,
plot3([fpbr1(1),fptr1(1)], [fpbr1(2),fptr1(2)], [fpbr1(3), fptr1(3)],'b', 'MarkerSize', 2, 'LineWidth',2);hold on,
text(src1Pos(1),src1Pos(2),src1Pos(3),'OB-1'); hold on;
plot3(fp1Posr(1), fp1Posr(2), fp1Posr(3),'rd', 'MarkerSize', 4, 'LineWidth',2); hold on,
plot3([fpbl1r(1),fpbr1r(1)], [fpbl1r(2),fpbr1r(2)], [fpbl1r(3), fpbr1r(3)],'r', 'MarkerSize', 2, 'LineWidth',2);hold on,
plot3([fptl1r(1),fptr1r(1)], [fptl1r(2),fptr1r(2)], [fptl1r(3), fptr1r(3)],'r', 'MarkerSize', 2, 'LineWidth',2);hold on,
plot3([fpbl1r(1),fptl1r(1)], [fpbl1r(2),fptl1r(2)], [fpbl1r(3), fptl1r(3)],'r', 'MarkerSize', 2, 'LineWidth',2);hold on,
plot3([fpbr1r(1),fptr1r(1)], [fpbr1r(2),fptr1r(2)], [fpbr1r(3), fptr1r(3)],'r', 'MarkerSize', 2, 'LineWidth',2);hold on,
text(src1Posr(1),src1Posr(2),src1Posr(3),'OB-1_{pitch}'); hold on;

plot3(fp2Pos(1), fp2Pos(2), fp2Pos(3),'bd', 'MarkerSize', 4, 'LineWidth',2); hold on,
plot3([fpbl2(1),fpbr2(1)], [fpbl2(2),fpbr2(2)], [fpbl2(3), fpbr2(3)],'b', 'MarkerSize', 2, 'LineWidth',2);hold on,
plot3([fptl2(1),fptr2(1)], [fptl2(2),fptr2(2)], [fptl2(3), fptr2(3)],'b', 'MarkerSize', 2, 'LineWidth',2);hold on,
plot3([fpbl2(1),fptl2(1)], [fpbl2(2),fptl2(2)], [fpbl2(3), fptl2(3)],'b', 'MarkerSize', 2, 'LineWidth',2);hold on,
plot3([fpbr2(1),fptr2(1)], [fpbr2(2),fptr2(2)], [fpbr2(3), fptr2(3)],'b', 'MarkerSize', 2, 'LineWidth',2);hold on,
text(src2Pos(1),src2Pos(2),src2Pos(3),'OB-2'); hold on;
plot3(fp2Posr(1), fp2Posr(2), fp2Posr(3),'rd', 'MarkerSize', 4, 'LineWidth',2); hold on,
plot3([fpbl2r(1),fpbr2r(1)], [fpbl2r(2),fpbr2r(2)], [fpbl2r(3), fpbr2r(3)],'r', 'MarkerSize', 2, 'LineWidth',2);hold on,
plot3([fptl2r(1),fptr2r(1)], [fptl2r(2),fptr2r(2)], [fptl2r(3), fptr2r(3)],'r', 'MarkerSize', 2, 'LineWidth',2);hold on,
plot3([fpbl2r(1),fptl2r(1)], [fpbl2r(2),fptl2r(2)], [fpbl2r(3), fptl2r(3)],'r', 'MarkerSize', 2, 'LineWidth',2);hold on,
plot3([fpbr2r(1),fptr2r(1)], [fpbr2r(2),fptr2r(2)], [fpbr2r(3), fptr2r(3)],'r', 'MarkerSize', 2, 'LineWidth',2);hold on,
text(src2Posr(1),src2Posr(2),src2Posr(3),'OB-2_{pitch}'); hold on;

hold on, xlabel('X'), ylabel('Y'), zlabel('Z');


end