function ok = checkDataIntegrity(data_mngt)
% (****) Structure of input 'data_mngt':
%    data_mngt.input = <PATH_TO_DIRECTORY (for 'gantry'), PATH_TO_IMAGE (for 'free')>
%    data_mngt.imagesRegExp = 'RADB'; % regular expression ['gantry' only]
%    data_mngt.sortBy = <'date','none'> % ['gantry' only]
%    data_mngt.anglesFromHeaders = <0,1>; % get gantry angles from xml
%    data_mngt.output = <PATH_TO_DIRECTORY>;
%    data_mngt.anglesFromHeaders = <0,1>;
%    data_mngt.anglesFname = <FILENAME_ANGLES> % ['gantry' only]
ok = 0;
if(~isstruct(data_mngt))
    return;
end

if(~isfield(data_mngt,'input'))
    return;
end
if(~isfield(data_mngt,'imagesRegExp'))
    return;
end
if(~isfield(data_mngt,'sortBy'))
    return;
end
if(~isfield(data_mngt,'anglesFromHeaders'))
    return;
end
if(~isfield(data_mngt,'output'))
    return;
end
if(~isfield(data_mngt,'anglesFname'))
    return;
end
ok = 1;
end