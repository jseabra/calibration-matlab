function writeRoomSetupXML(xmlFileName, name, place, system, panel, axis, emplacement, sad, aid, varargin) 
% WRITEROOMSETUPXML Writes the room setup in xml format
% --------------------------------------------

% constants:
NRMANDARGS=9;
NROPTARGSFREE=6;
NROPTARGSGANTRY=1;
EMPLACEMENTFREE = 'free';
EMPLACEMENTGTR = 'gantry';

fpName = panel.Name;
dims = panel.Dimension;
pixel = panel.Pixel;
radT = panel.Phys2RadT;

%fprintf('Number of inputs = %d\n', nargin)
optargin = size(varargin,2);
stdargin = nargin - optargin;

if stdargin~=NRMANDARGS,
   disp(['Number of mandatory inputs must be ', num2str(NRMANDARGS), '. Aborting']);
   return
end

if strcmp(emplacement, EMPLACEMENTFREE),
    if optargin~=NROPTARGSFREE
        disp(['Number of optional inputs must be ', num2str(NROPTARGSFREE), ' for free emplacement. Aborting.']);
        return
    else
        pitchAng = varargin{1};
        rollAng = varargin{2};
        yawAng = varargin{3};
        tx = varargin{4};
        ty = varargin{5};
        tz = varargin{6};
    end
    
elseif strcmp(emplacement, EMPLACEMENTGTR),
    if optargin~=NROPTARGSGANTRY,
        disp(['Number of optional inputs must be ', num2str(NROPTARGSGANTRY), ' for gantry emplacement. Aborting.']);
        return
    else
        angleOffset = varargin{1};
        
    end
end

dimX = dims(1);
dimY = dims(2);
pixelX = pixel;
pixelY = pixel;

docNode = com.mathworks.xml.XMLUtils.createDocument('iMagXRoomSetup');
docRootNode = docNode.getDocumentElement;
docRootNode.setAttribute('roomName', name);
docRootNode.setAttribute('system', system);
docRootNode.setAttribute('place', place);

% node 'FlatPanel': 
node = docNode.createElement('FlatPanel');
node.setAttribute('name', fpName);
node.setAttribute('axis', axis);
node.setAttribute('emplacement', emplacement);
docRootNode.appendChild(node);

% child 'pixelSize'
param = docNode.createElement('pixelSize');
node.appendChild(param);
param.setAttribute('x', sprintf('%g',pixelX));
param.setAttribute('y', sprintf('%g',pixelY));

% child 'dimensions'
param = docNode.createElement('dimensions');
node.appendChild(param);
param.setAttribute('x', sprintf('%i',dimX));
param.setAttribute('y', sprintf('%i',dimY));

% child 'Distances'
param = docNode.createElement('Distances');
node.appendChild(param);
param.setAttribute('sid', sprintf('%i',sad));
param.setAttribute('sdd', sprintf('%i',sad+aid));

% child 'transformToRad'
T = createTransform2RadFromPanelProps(radT,panel.PhysDimX,panel.PhysDimY);
param = docNode.createElement('transformToRad');
node.appendChild(param);
for i=1:size(T,1)
    for j=1:size(T,2)
        str = sprintf('a%i%i',i-1,j-1);
        param.setAttribute(str,  sprintf('%d',T(i,j)));
    end
end

if strcmp(emplacement, EMPLACEMENTGTR)
    % child 'ProjectionAxis'
    param = docNode.createElement('offset');
    node.appendChild(param);
    param.setAttribute('angle', sprintf('%d',angleOffset));
    
elseif strcmp(emplacement, EMPLACEMENTFREE)
    % child 'transformToFrs'
    param = docNode.createElement('transformToFrs');
    node.appendChild(param);
    param.setAttribute('pitch', sprintf('%g',pitchAng));
    param.setAttribute('roll', sprintf('%g',rollAng));
    param.setAttribute('yaw', sprintf('%g',yawAng));
    param.setAttribute('tx', sprintf('%g',tx));
    param.setAttribute('ty', sprintf('%g',ty));
    param.setAttribute('tz', sprintf('%g',tz));
    
else
    disp('Emplacement must either be gantry or fixed. Aborting.');
    return
    
end

% write model into file:
xmlwrite(xmlFileName,docNode);
type(xmlFileName);
end