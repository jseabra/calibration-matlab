function test_geometry3D2D()

srcA = [1212.71, -1503.73, -1517.89];

srcB = [-1365.24, -1470.98, -1608.73];

fpA = [0.80789, -0.383467, 0.447512, -351.889;...
    0.0747818, -0.68651, -0.723265, 441.325;...
    0.58457, 0.617784, -0.525948, 459.114;...
    0, 0, 0, 1];

fpB = [0.804613, -0.126276, -0.580217, -5.26379;...
    0.35309, -0.683865, 0.638479, -153.487;...
    -0.477415, -0.718598, -0.50566, 742.858;...
    0, 0, 0, 1];


vec_a = [];
pt0 = [10,10,1,1]; % frs (x,y,z)

for a=-1:0.1:1,


pt1 = pitch(a,[0,0,0]) * pt0';
pt1 = pt1';

disp(['pt0: ', num2str(pt0)]);
disp(['after rotation: ', num2str(pt1)]);

pt0_A = projectFRSOntoPanel(srcA, fpA, pt0);
pt1_A = projectFRSOntoPanel(srcA, fpA, pt1);

curA = acosd((pt0_A*pt1_A')/(norm(pt0_A)*norm(pt1_A))); % angle between two vectors
vec_a = [vec_a, curA];

end

hold on;
plot(-1:0.1:1,vec_a);


%disp(['projected pt0_A: ', num2str(pt0_A)]);
%disp(['after rotation: ', num2str(pt1_A)]);


end

function Pt=projectFRSOntoPanel(src, fp, ptFRS)

% Position of the flat panel in FRS. The source is located at the
% origin of the radiographic CS:
O4 = fp * [0,0,0,1]';
% Component of the rad axis vectors:
u4 = fp * [1,0,0,1]';
v4 = fp * [0,1,0,1]';
O = O4(1:3);
u = u4(1:3)-O;
v = v4(1:3)-O;
O = O';
u = u';
v = v';

if size(ptFRS,1)==1,
    ptFRS=ptFRS';
end

P = ptFRS(1:3)';
t = P-src;
Pt = InterPlaneLine(u,v,O,t,P);

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

end