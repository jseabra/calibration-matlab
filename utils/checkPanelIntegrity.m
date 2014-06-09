function ok = checkPanelIntegrity(panel)
% (*) Structure of input 'panel':
%    panel.Name = <flat_panel_model>; e.g.: THALES2630CS
%    panel.Dimension = [<DIM_X_PIXELS> <DIM_Y_PIXELS>];
%    panel.Pixel = <PIXEL_SIZE_MM>;
%    panel.Phys2RadT = <'identity','rotationX','rotationY','rotationZ'>;
%    panel.PhysDimX = panel.Dimension(1)/2*panel.Pixel;
%    panel.PhysDimY = panel.Dimension(2)/2*panel.Pixel;
%    panel.invertGray = <0,1>;
%    panel.factorLow = <ADJUST_LUT_TO_LOW_FACTOR>
%    panel.factorHigh = <ADJUST_LUT_TO_HIGH_FACTOR>
ok = 0;
if(~isstruct(panel))
    return;
end

if(~isfield(panel,'Name'))
    return;
end

if(isfield(panel,'Dimension'))
    if(panel.Dimension(1)<=0)
        return;
    end
    if(panel.Dimension(2)<=0)
        return;
    end
else
    return;
end

if(isfield(panel,'Pixel'))
    if(panel.Pixel<=0)
        return;
    end
else
    return;
end

if(isfield(panel,'Phys2RadT'))
    if(strcmp(panel.Phys2RadT,'identity')==0 &&...
            strcmp(panel.Phys2RadT,'rotationX')==0 &&...
            strcmp(panel.Phys2RadT,'rotationY')==0 &&...
            strcmp(panel.Phys2RadT,'rotationZ')==0)
        return;
    end
else
    return;
end

if(isfield(panel,'PhysDimX'))
    DimX = panel.Dimension(1)/2*panel.Pixel;
    if(panel.PhysDimX~=DimX)
        return;
    end
else
    return;
end

if(isfield(panel,'PhysDimY'))
    DimY = panel.Dimension(2)/2*panel.Pixel;
    if(panel.PhysDimY~=DimY)
        return;
    end
else
    return;
end

if(isfield(panel,'invertGray'))
    if(panel.invertGray~=0 && panel.invertGray~=1)
        return;
    end
else
    return;
end

ok = 1;
end