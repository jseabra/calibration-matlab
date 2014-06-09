function value=getValueFromXMLByTag(filename, tagName)
%GETVALUEFROMXMLBYTAG parses an xml file and returns the value associated
%to a tag given as argument. Returns empty if no value associated to the
%tagName is found
theStruct = parseXML(filename);

value = str2num(getStructValueByTag(theStruct, tagName));

end

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

function value = getStructValueByTag(struct, tagName)
value = [];

for j=1:numel(struct.Attributes)
%     struct.Attributes(j).Name
    if(strcmp(struct.Attributes(j).Name, tagName))
        value = struct.Attributes(j).Value;
        return
    end
end

NrElements = numel(struct.Children);

for i=1:NrElements
    newElement = struct.Children(i);
    
    for j=1:numel(newElement.Attributes)
%         newElement.Attributes(j).Name
        if(strcmp(newElement.Attributes(j).Name, tagName))
            value = newElement.Attributes(j).Value;
            return
        end
    end
    
    for i=1:numel(newElement.Children)
        newNode = newElement.Children(i);
        
        for j=1:numel(newNode.Attributes)
%             newNode.Attributes(j).Name
            if(strcmp(newNode.Attributes(j).Name, tagName))
                value = newNode.Attributes(j).Value;
                return
            end
        end
        
        for k=1:numel(newNode.Children)
            newChild = newNode.Children(i);
            
            for w=1:numel(newChild.Attributes)
%                 newChild.Attributes(w).Name
                if(strcmp(newChild.Attributes(w).Name,tagName))
                    value = newChild.Attributes(w).Value;
                    return
                end
            end
        end
        
    end
    
end

end

   
