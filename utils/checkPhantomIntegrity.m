function ok = checkPhantomIntegrity(phantom)
% (***) Structure of input 'phantom':
%    phantom.spheres = <list_of_spheres>; % spheres(id,:) = (X_frs,Y_frs,Z_frs)
%    phantom.offset = [X_frs,Y_frs,Z_frs]; % offset central sphere to
%    isocenter
ok = 0;
if(~isstruct(phantom))
    return;
end

if(isfield(phantom,'spheres'))
    if(size(phantom.spheres,1)==0)
        return;
    end
else
    return;
end

if(isfield(phantom,'offset'))
    if(numel(phantom.offset)~=3)
        return;
    end
else
    return;
end
ok = 1;
end