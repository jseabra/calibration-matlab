function ok = checkConfigIntegrity(config)
% (*****) Structure of input 'config':
%    config.dof = [Sx, Sy, Sz, Dx, Dy, Dz, Rx, Ry, Rz]; % <0=constant,
%    1=optimize)
%    config.spAttrb (see ******)
%    config.optimizer = <'cmaes',...>
%    config.debugMode = <0,1>
%
% (******) Structure of input 'config.spAttrb':
%    spAttrb.beadSize = <SPHERE_DIAMETER_MM>
%    spAttrb.projDiam = <PROJECTED_SPHERE_DIAM_PIXELS>
%    spAttrb.rSize = <SEARCH_WINDOW_WIDTH_PIXELS>
%    spAttrb.sMin  = <MIN_NR_OBJ_PIXELS>
%    spAttrb.sMax  = <MAX_NR_OBJ_PIXELS>
%    spAttrb.dMin  = <MIN_DIAM_OBJ_PIXELS>
%    spAttrb.dMax  = <MAX_DIAM_OBJ_PIXELS>
%    spAttrb.ratio = <RATIO_OBJ_H/V>
%    spAttrb.min2diff = <MAX_DIST_BETWEEN_OBJ_PIXELS>
%    spAttrb.threshold = <THRESHOLD_BETWEEN_-1_1>
%    spAttrb.do_filter = <0,1> % median filter
%    spAttrb.ecc = <ECCENTRICITY_0_1> % 1 for ideal symmetry
%    spAttrb.imopen_param = <REMOVE_CLUTTER_0_1>
%    spAttrb.detectionMode = <'auto','manual'>
ok=0;
if(~isstruct(config))
    return;
end

if(~isfield(config,'dof'))
    return;
end
if(~isfield(config,'optimizer'))
    return;
end
if(~isfield(config,'debugMode'))
    return;
end
if(~isfield(config,'spAttrb'))
    return;
end
if(~isfield(config.spAttrb,'beadSize'))
    return;
end
if(~isfield(config.spAttrb,'projDiam'))
    return;
end
if(~isfield(config.spAttrb,'rSize'))
    return;
end
if(~isfield(config.spAttrb,'sMin'))
    return;
end
if(~isfield(config.spAttrb,'sMax'))
    return;
end
if(~isfield(config.spAttrb,'dMin'))
    return;
end
if(~isfield(config.spAttrb,'dMax'))
    return;
end
if(~isfield(config.spAttrb,'ratio'))
    return;
end
if(~isfield(config.spAttrb,'min2diff'))
    return;
end
if(~isfield(config.spAttrb,'threshold'))
    return;
end
if(~isfield(config.spAttrb,'do_filter'))
    return;
end
if(~isfield(config.spAttrb,'ecc'))
    return;
end
if(~isfield(config.spAttrb,'detectionMode'))
    return;
end
ok=1;
end