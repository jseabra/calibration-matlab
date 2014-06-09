function [coeffs_Px, coeffs_Py, coeffs_Pz, coeffs_Rx, coeffs_Ry, coeffs_Rz, coeffs_Sx, coeffs_Sy, coeffs_Sz] = getFlexmapCoeffsFromModelXML(fname_model, calibrationAxis)
% GETFLEXMAPCOEFFSFROMMODELXML: get deformation coefficients
% (a0,a1,a2,a3,a4) for each radiographic deformation component from a
% geometric model given in XML format
% 
% INPUTS:
% xml file containing a flexmap model with the following structure (given as example, do not use):
%
%    <FlatPanel name="Varian4030" axis="RADB">
%            <DetectorTranslation><!-- Model coefficients for detector translation + Mean Squared Error Model vs. Raw meas. -->
%             <Px MSE="2.48519" a0="13.5039" a1="-2.57e-05" a2="-1.64171" a3="0.996129" a4="0.141973" />
%             <Py MSE="1.24685" a0="-13.0596" a1="-6.76125e-05" a2="0.858309" a3="1.00626" a4="0.968553" />
%             <Pz MSE="-1" a0="1.81243" a1="0" a2="0" a3="0" a4="0" />
%        </DetectorTranslation>
%        <DetectorRotation><!-- Model coefficients for detector rotation + Mean Squared Error Model vs. Raw meas. -->
%             <Rx MSE="-1" a0="0.190462" a1="0" a2="0" a3="0" a4="0" />
%             <Ry MSE="-1" a0="1.35318" a1="0" a2="0" a3="0" a4="0" />
%             <Rz MSE="0.0943452" a0="0.019252" a1="-3.92204e-07" a2="0.075318" a3="0.998621" a4="-0.349217" />
%        </DetectorRotation>
%        <SourceTranslation><!-- Model coefficients for source translation + Mean Squared Error Model vs. Raw meas. -->
%             <Tx MSE="-1" a0="-2.0083" a1="0" a2="0" a3="0" a4="0" />
%             <Ty MSE="4.25088" a0="0.601343" a1="8.55825e-05" a2="5.21253" a3="1.0028" a4="4.46338" />
%             <Tz MSE="-1" a0="27.8726" a1="0" a2="0" a3="0" a4="0" />
%        </SourceTranslation>
%        <offset angle="-90" />
%        <Distances sid="2875" sdd="3470" />
%    </FlatPanel>
%
% OUTPUS:
% coeffs_Px, coeffs_Py, coeffs_Pz: coefficients (x,y,z) [mm] of flat panel
% translation
% coeffs_Rx, coeffs_Ry, coeffs_Rz: coefficients (x,y,z) [mm] of flat panel
% rotation
% coeffs_Sx, coeffs_Sy, coeffs_Sz: coefficients (x,y,z) [mm] of source
% translation 
% -----------------------------------------

struct = parseXML(fname_model);
[coeffs_Px, coeffs_Py, coeffs_Pz,...
    coeffs_Rx, coeffs_Ry, coeffs_Rz,...
    coeffs_Sx, coeffs_Sy, coeffs_Sz] = getCoefficientsFromStruct(struct, calibrationAxis);
end

% ----- ROUTINE TO PARSE MODEL COEFFICIENTS FROM STRUCT OBTAINED FROM XML PARSER -----
function [coeffs_Px, coeffs_Py, coeffs_Pz, coeffs_Rx, coeffs_Ry, coeffs_Rz, coeffs_Sx, coeffs_Sy, coeffs_Sz, sid, sdd, angleOffset] = getCoefficientsFromStruct(struct, calibrationAxis)

coeffs_Px = [];
coeffs_Py = [];
coeffs_Pz = [];
coeffs_Rx = [];
coeffs_Ry = [];
coeffs_Rz = [];
coeffs_Sx = [];
coeffs_Sy = [];
coeffs_Sz = [];
sid = [];
sdd = [];
angleOffset = [];

tagFlatPanel = 'FlatPanel';
tagAxis = 'axis';
tagDetectorTranslation = 'DetectorTranslation';
tagDetectorRotation = 'DetectorRotation';
tagSourceTranslation = 'SourceTranslation';
tagPx = 'Px'; tagPy = 'Py'; tagPz = 'Pz';
tagRx = 'Rx'; tagRy = 'Ry'; tagRz = 'Rz';
tagSx = 'Tx'; tagSy = 'Ty'; tagSz = 'Tz';

NrElements = numel(struct.Children);

for i=1:NrElements,
    
    newElement = struct.Children(i);
    if(strcmp(newElement.Name, tagFlatPanel))
        
        for j=1:numel(newElement.Attributes)
            
            if(strcmp(newElement.Attributes(j).Name, tagAxis))
                
                if(strcmp(newElement.Attributes(j).Value, calibrationAxis))
                     
                    NrNodes = numel(newElement.Children);
                    
                    % parse nodes and children:
                    for k=1:NrNodes,
                        
                        newNode = newElement.Children(k);
                        if( strcmp(newNode.Name,tagDetectorTranslation) )
                            
                            for l=1:numel(newNode.Children),
                                newChild = newNode.Children(l);
                                if( strcmp(newChild.Name,tagPx) )
                                    coeffs_Px = parse_coefficients(newChild);
                                elseif ( strcmp(newChild.Name,tagPy) )
                                    coeffs_Py = parse_coefficients(newChild);
                                elseif ( strcmp(newChild.Name,tagPz) )
                                    coeffs_Pz = parse_coefficients(newChild);
                                end
                            end
                            
                        elseif ( strcmp(newNode.Name,tagDetectorRotation) )
                            
                            for l=1:numel(newNode.Children),
                                newChild = newNode.Children(l);
                                if( strcmp(newChild.Name,tagRx) )
                                    coeffs_Rx = parse_coefficients(newChild);
                                elseif ( strcmp(newChild.Name,tagRy) )
                                    coeffs_Ry = parse_coefficients(newChild);
                                elseif ( strcmp(newChild.Name,tagRz) )
                                    coeffs_Rz = parse_coefficients(newChild);
                                end
                            end
                            
                        elseif ( strcmp(newNode.Name,tagSourceTranslation) )
                            
                            for l=1:numel(newNode.Children),
                                newChild = newNode.Children(l);
                                if( strcmp(newChild.Name,tagSx) )
                                    coeffs_Sx = parse_coefficients(newChild);
                                elseif ( strcmp(newChild.Name,tagSy) )
                                    coeffs_Sy = parse_coefficients(newChild);
                                elseif ( strcmp(newChild.Name,tagSz) )
                                    coeffs_Sz = parse_coefficients(newChild);
                                end
                            end
                            
                        end
                        
                    end
                    % end parsing nodes
                    
                end
                
            end
            
        end
        
    end
    
end %endfor

if (isempty(coeffs_Px) || isempty(coeffs_Py) || isempty(coeffs_Pz) ||...
        isempty(coeffs_Sx) || isempty(coeffs_Sy) || isempty(coeffs_Sz) ||...
        isempty(coeffs_Rx) || isempty(coeffs_Ry) || isempty(coeffs_Rz))
    fprintf('Unable to parse coefficients. Check tags and calibration axis\n');
    return 
end

end

% -----------------------------------------
function theStruct = parseXML(filename)
% PARSEXML Convert XML file to a MATLAB structure.
try
   tree = xmlread(filename);
catch
   error('Failed to read XML file %s.',filename);
end

% Recurse over child nodes. This could run into problems 
% with very deeply nested trees.
try
   theStruct = parseChildNodes(tree);
catch
   error('Unable to parse XML file %s.',filename);
end

end


% ----- Subfunction PARSECHILDNODES -----
function children = parseChildNodes(theNode)
% Recurse over node children.
children = [];
if theNode.hasChildNodes
   childNodes = theNode.getChildNodes;
   numChildNodes = childNodes.getLength;
   allocCell = cell(1, numChildNodes);

   children = struct(             ...
      'Name', allocCell, 'Attributes', allocCell,    ...
      'Data', allocCell, 'Children', allocCell);

    for count = 1:numChildNodes
        theChild = childNodes.item(count-1);
        children(count) = makeStructFromNode(theChild);
    end
end
end


% ----- Subfunction MAKESTRUCTFROMNODE -----
function nodeStruct = makeStructFromNode(theNode)
% Create structure of node info.

nodeStruct = struct(                        ...
   'Name', char(theNode.getNodeName),       ...
   'Attributes', parseAttributes(theNode),  ...
   'Data', '',                              ...
   'Children', parseChildNodes(theNode));

if any(strcmp(methods(theNode), 'getData'))
   nodeStruct.Data = char(theNode.getData); 
else
   nodeStruct.Data = '';
end
end


% ----- Subfunction PARSEATTRIBUTES -----
function attributes = parseAttributes(theNode)
% Create attributes structure.

attributes = [];
if theNode.hasAttributes
   theAttributes = theNode.getAttributes;
   numAttributes = theAttributes.getLength;
   allocCell = cell(1, numAttributes);
   attributes = struct('Name', allocCell, 'Value', ...
                       allocCell);

   for count = 1:numAttributes
      attrib = theAttributes.item(count-1);
      attributes(count).Name = char(attrib.getName);
      attributes(count).Value = char(attrib.getValue);
   end
end
end

% -----------------------------------------
function coeffs = parse_coefficients(newChild)

coeffs = -1.*ones(1,5);

for k=1:numel(newChild.Attributes)
    newAttribute = newChild.Attributes(k);
    
    if (strcmp(newAttribute.Name,'a0'))
        coeffs(1) = str2double(newAttribute.Value);
    elseif (strcmp(newAttribute.Name,'a1'))
        coeffs(2) = str2double(newAttribute.Value);
    elseif (strcmp(newAttribute.Name,'a2'))
        coeffs(3) = str2double(newAttribute.Value);
    elseif (strcmp(newAttribute.Name,'a3'))
        coeffs(4) = str2double(newAttribute.Value);
    elseif (strcmp(newAttribute.Name,'a4'))
        coeffs(5) = str2double(newAttribute.Value);
    end
    
end

if ~isempty(coeffs==-1) return; end
end
