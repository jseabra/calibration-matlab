function ok = checkGeometryIntegrity(geometry)
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
ok = 0;
if(~isstruct(geometry))
    return;
end

if(~isfield(geometry,'roomName'))
    return;
end

if(~isfield(geometry,'roomSystem'))
    return;
end

if(~isfield(geometry,'roomPlace'))
    return;
end

if(~isfield(geometry,'axis'))
    return;
end

if(isfield(geometry,'config'))
    if(strcmp(geometry.config,'gantry')==0 && strcmp(geometry.config,'free')==0)
        return;
    end
else
    return;
end

if(isfield(geometry,'Vec'))
    if(numel(geometry.Vec)~=9)
       return; 
    end
else
    return;
end

if(isfield(geometry,'SAD'))
    if(geometry.SAD<0)
       return; 
    end
else
    return;
end

if(isfield(geometry,'AID'))
    if(geometry.AID<0)
       return; 
    end
else
    return;
end

if(geometry.AID>geometry.SAD)
    return;
end

if(strcmp(geometry.config,'gantry'))
    if(~isfield(geometry,'gantryAngleOffset'))
        return;
    end
end

if(strcmp(geometry.config,'free'))
    if(isfield(geometry,'srcPosition'))
        if(numel(geometry.srcPosition)~=3)
            return;
        end
    else
        return;
    end
    if(isfield(geometry,'detectorT'))
        if(size(geometry.detectorT,1)~=4)
            return;
        end
        if(size(geometry.detectorT,2)~=4)
            errorNotify();
            return;
        end
    else
        return;
    end
    if(~isfield(geometry,'rx'))
        return;
    end
    if(~isfield(geometry,'ry'))
        return;
    end
    if(~isfield(geometry,'rz'))
        return;
    end
    if(~isfield(geometry,'tx'))
        return;
    end
    if(~isfield(geometry,'ty'))
        return;
    end
    if(~isfield(geometry,'tz'))
        return;
    end
end
ok=1;
end