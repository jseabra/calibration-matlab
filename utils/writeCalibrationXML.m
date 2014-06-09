function writeCalibrationXML(xmlFileName, name, place, datetime, system, fpName, axis, P, R, T)
% WRITECALIBRATIONXML Writes the geometric model on xml format, with the following structure:
% Room:
%     roomName   place    system    date
% FlatPanel:
%     name    axis
%     DetectorTranslation
%     DetectorRotation
%     SourceTranslation
% --------------------------------------------

% geometrical model tags:
tags={'a0', 'a1', 'a2', 'a3', 'a4', 'MSE'};

docNode = com.mathworks.xml.XMLUtils.createDocument('iMagXGeometricalCalibrationModel');
docRootNode = docNode.getDocumentElement;
docRootNode.setAttribute('roomName', name);
docRootNode.setAttribute('place', place);
docRootNode.setAttribute('system', system);
docRootNode.setAttribute('time', datestr(datetime,13));
docRootNode.setAttribute('date', datestr(datetime,23));
docRootNode.setAttribute('uid', datestr(datetime,30));

% node 'FlatPanel': 
node = docNode.createElement('FlatPanel');
node.setAttribute('name', fpName);
node.setAttribute('axis', axis);
docRootNode.appendChild(node);

% child 'DetectorTranslation'
param = docNode.createElement('DetectorTranslation');
node.appendChild(param);
param.appendChild(docNode.createComment('Detector translation parameters'));
ptags={'Px','Py','Pz'};
for pt_id=1:length(ptags)
    thisElement = docNode.createElement(ptags(pt_id)); 
    for t_id=1:length(tags)
        thisElement.setAttribute(tags(t_id),sprintf('%i',P(t_id,pt_id)));
    end
    param.appendChild(thisElement);
end

% child 'DetectorRotation'
param = docNode.createElement('DetectorRotation');
node.appendChild(param);
param.appendChild(docNode.createComment('Detector rotation parameters'));
ptags={'Rx','Ry','Rz'};
for pt_id=1:length(ptags)
    thisElement = docNode.createElement(ptags(pt_id)); 
    for t_id=1:length(tags)
        thisElement.setAttribute(tags(t_id),sprintf('%i',R(t_id,pt_id)));
    end
    param.appendChild(thisElement);
end

% child 'SourceTranslation'
param = docNode.createElement('SourceTranslation');
node.appendChild(param);
param.appendChild(docNode.createComment('Source translation parameters'));
ptags={'Tx','Ty','Tz'};
for pt_id=1:length(ptags)
    thisElement = docNode.createElement(ptags(pt_id)); 
    for t_id=1:length(tags)
        thisElement.setAttribute(tags(t_id),sprintf('%i',T(t_id,pt_id)));
    end
    param.appendChild(thisElement);
end

xmlwrite(xmlFileName,docNode);
type(xmlFileName);
end