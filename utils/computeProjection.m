function F = computeProjection(phantom, geometry, panel, objOffset)
% COMPUTEPROJECTION calculates the projection of object 'phantom' onto the
% flat panel 'panel' given a known geometry 'geometry'

if isempty(phantom),
    disp('Please check: empty phantom.');
    return
end
if isempty(geometry),
    disp('Please check: empty geometry.');
    return
end
if isempty(panel),
    disp('Please check: empty panel.');
    return
end

F = [];

    Vec = geometry.Vec;

if (strcmp(geometry.config,'gantry'))
 
    gantryAngle = geometry.gantryAngle;
    angleOffset = geometry.gantryAngleOffset;
    SAD = geometry.SAD;
    AID = geometry.AID;
    
    if AID<0,
        disp('Axis to Isocenter distance must be positive.');
        return
    end
    
    if SAD<0,
        disp('Source to Axis distance must be positive.');
        return
    end
         
    M = roll(gantryAngle+angleOffset,[0,0,0])*trans(Vec(1), Vec(2), Vec(3)+SAD);
    % Position of the source in the FRS. The source is located at the origin of the radiographic CS
    X4 = M * [0,0,0,1]';
    sourcePosition = X4(1:3)';
    
    % Detector position in FRS:    
    detectorPosition = roll(gantryAngle+angleOffset,[0,0,0])*rot(Vec(9),[0,0,0])*trans(Vec(4),Vec(5),Vec(6)-AID)*roll(Vec(8),[0,0,0])*pitch(Vec(7),[0,0,0]);
 
elseif strcmp(geometry.config,'free')
    
    if(isfield(geometry,'srcPosition'))

        %sourcePosition = trans(Vec(1),Vec(2),Vec(3)) *
        %[geometry.srcPosition; 1]; % translation parameters given in FRS 
        %sourcePosition = trans(geometry.srcPosition(1),geometry.srcPosition(2),geometry.srcPosition(3)) * trans(Vec(1),Vec(2),Vec(3));
        %sourcePosition = sourcePosition(1:3,4)';
        
        % express in FRS CS:
        sourcePosition = trans(Vec(1),Vec(2),Vec(3))*[geometry.srcPosition;1];
        sourcePosition = sourcePosition(1:3)';      
        
    else
        srcT = geometry.srcT;
        sourcePosition = srcT*trans(Vec(1),Vec(2),Vec(3))*[0,0,0,1]'; sourcePosition = sourcePosition(1:3)';        
    end
    
    detT = geometry.detectorT;  
    % Usual way:
    detectorPosition = detT*rot(Vec(9),[0,0,0])*trans(Vec(4),Vec(5),Vec(6))*roll(Vec(8),[0,0,0])*pitch(Vec(7),[0,0,0]);
   
    % New way: do all translations first, then rotations
    %detectorPosition =
    %detT*trans(Vec(4),Vec(5),Vec(6))*rot(Vec(9),[0,0,0])*roll(Vec(8),[0,0,
    %0])*pitch(Vec(7),[0,0,0]);
    
end

F  = ProjectDRR(phantom, sourcePosition, detectorPosition , panel, objOffset);

end

% -----------------------------------------
% Project object onto the detector.
% F is the list of projected points in image coordinates (pixels)
function F = ProjectDRR(phantom, sourcePosition, detectorPosition, panel, ObjPos)

% Position of the flat panel in the FRS. The source is located at the
% origin of the radiographic CS:
O4 = detectorPosition * [0,0,0,1]';
% Component of the rad axis vectors:
u4 = detectorPosition * [1,0,0,1]';
v4 = detectorPosition * [0,1,0,1]';
O = O4(1:3);
u = u4(1:3)-O;
v = v4(1:3)-O;
O = O';
u = u';
v = v';

pxl = panel.Pixel;
Dim = panel.Dimension;

M = trans(ObjPos(1),ObjPos(2),ObjPos(3));
sPhan = size(phantom);
Corner = [Dim(1)/2, Dim(2)/2]; %Coordinates of the corner of the flat panel

for i = 1:sPhan(1)
    Obj =	M*[phantom(i,:) , 1]';
    P = Obj(1:3)';
    t = P-sourcePosition;
    Pt = InterPlaneLine(u,v,O,t,P);
    tmp = Pt./pxl+Corner;
    F(i,:) =  tmp;
    % Note: Y is flipped to be plotted on image oriented X:left->right,
    % Y:bottom->up F(:,i) = (X,Y)
    F(i,:) =  [tmp(1), Dim(2)-tmp(2)];
end

end

% -----------------------------------------
function [ I , k ] = InterPlaneLine(u,v,O,t,P)

if(CheckLine(u) ~= 1)
    error('Vector u is not a line vector');
end
if(CheckLine(v) ~= 1)
    error('Vector v is not a line vector');
end
if(CheckLine(O) ~= 1)
    error('Vector O is not a line vector');
end

% Let's define a system of equation A.x = B
B = P - O; % ok!
B = B' ;
A = [u'  v'  -t'];

%Now we solve the system
x = A\B;
I = [x(1) x(2)]; %Coordinate of the intersection point in the plane coordiantes
k = x(3); % Index along the line

end

% -----------------------------------------
function bool = CheckLine(vec)

svec = size(vec);
if (svec(1) ~= 1)
    bool = 0;
    return
end
bool = 1;
end