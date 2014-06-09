function [absError, relError,reconstructed_point] = reconstructPoint3D(phantom, objOffset, p_gt, beam1, beam2)
% RECONSTRUCTPOINT3D Reconstructs a point in 3D space from two points
% manually selected in the DRs.
%
% INPUTS:
% phantom (list N x 3 of sphere positions in FRS x,y,z)
% p_gt (ground truth sphere position in FRS x,y,z)
% beam1 (structure containing: panel, geometry and dr (See
% pointReconstruction for details))
% beam2 (same structure as for beam1)
%
% OUTPUT (reconstruction errors):
% absError: absolute error corresponds to the L2-distance [mm] between 
% reconstructed point and grount truth position
% relError: relative error corresponds to the minimum distance between the
% two backprojections (where rays should cross)
% -------------------------------------------------------------------------
sourcePositions = [];
markerPositions = [];

getColorStyle = [{'r'},{'b'},{'c'},{'m'}];


for index=1:2, % process for both axes
    
    if index==1,
        % Rad One:
        panel = beam1.panel;
        dr = beam1.dr;
        geometry = beam1.geometry;
        corner = [panel.Dimension(1)/2 panel.Dimension(2)/2];
        center = [0,0,0];
        axisName = 'RadOne';
    else
        % Rad Two:
        panel = beam2.panel;
        dr = beam2.dr;
        geometry = beam2.geometry;
        corner = [panel.Dimension(1)/2 panel.Dimension(2)/2];
        center = [0,0,0];
        axisName = 'RadTwo';
    end
    
    % get deformation parameters for a given axis:
    
    if(isfield(geometry,'model'))
        
        [coeffs_Px, coeffs_Py, coeffs_Pz,...
            coeffs_Rx, coeffs_Ry, coeffs_Rz,...
            coeffs_Sx, coeffs_Sy, coeffs_Sz] = getFlexmapCoeffsFromModelXML(geometry.model, geometry.axis);
        
        if strcmp(geometry.config,'gantry')
            [P, T, R] = getDeformationByAngle(geometry.gantryAngle, coeffs_Px, coeffs_Py, coeffs_Pz, coeffs_Rx, coeffs_Ry, coeffs_Rz, coeffs_Sx, coeffs_Sy, coeffs_Sz);
        else
            [P, T, R] = getDeformationByAngle(0, coeffs_Px, coeffs_Py, coeffs_Pz, coeffs_Rx, coeffs_Ry, coeffs_Rz, coeffs_Sx, coeffs_Sy, coeffs_Sz);
        end
        
    else
        
        if(~isfield(geometry,'Vec'))
            disp('Geometry must contain field named Vec');
            return
        else
            T = geometry.Vec(1:3);
            P = geometry.Vec(4:6);
            R = geometry.Vec(7:9);            
        end
        
    end
    
    % marker selection:
    F = computeProjection(phantom, geometry, panel, objOffset);
    figure, imshow(dr,[], 'InitialMagnification', 500), hold on, zoom(1.5)
    set(gcf, 'units','normalized','outerposition',[0, 0, 1, 1]);
    for idx = 1:size(F,1),
        plot (F(idx,1), F(idx,2) , 'c.');
    end
    for idx = 1:size(F,1),
        text (F(idx,1), F(idx,2) , num2str(idx),'color','c');
    end
    hold on;
    
    % draw ray projection on 'second' dr:
    if index==2,

        F2 = computeProjection(backRay, geometry, panel, [0,0,0]);
        for idx = 1:size(F2,1),
            plot(F2(idx,1), F2(idx,2),'m.','LineWidth',1);
        end
        hold on;
    end
    
    % first coordinate is X (along column), second is Y (along row). X is
    % pointing right, Y is pointing up
    [X,Y] = ginput(1);
    close
    
    % compute backprojection (with geometry correction):
    detectorPointRAD = [(X - corner(1))*panel.Pixel, (-Y + corner(2))*panel.Pixel];
    
    if strcmp(geometry.config,'gantry')
        [sourcePosFRS, markerPosFRS, dirVec] = getBackprojectionRay(detectorPointRAD, geometry.gantryAngle+geometry.gantryAngleOffset, center, P, T, R, geometry.SAD, geometry.AID);
    elseif strcmp(geometry.config,'fixed') || strcmp(geometry.config,'free')
        %[sourcePosFRS, markerPosFRS, dirVec] = getBackprojectionRay(detectorPointRAD, 0, center, P, T, R, geometry.SAD, geometry.AID); 
        markerPosFRS = geometry.detectorT * [detectorPointRAD';0;1];
        markerPosFRS = markerPosFRS(1:3)';
        sourcePosFRS = geometry.srcPosition';
        dirVec = (markerPosFRS-sourcePosFRS);
    end
    
    backRay = repmat(sourcePosFRS,1e3,1) + (linspace(0,1,1e3))'*dirVec;
    
    disp('* * * geometry * * *');
    disp(sprintf('source pos frs [mm]: %6.4f %6.4f %6.4f', sourcePosFRS));
    disp(sprintf('marker pos frs [mm]: %6.4f %6.4f %6.4f', markerPosFRS));
    
    sourcePositions = [sourcePositions; sourcePosFRS];
    markerPositions = [markerPositions; markerPosFRS];
    
    if index==1,
    hf = figure; hold on,
    xlabel('X_{FRS}');
    ylabel('Y_{FRS}');
    zlabel('Z_{FRS}');
    else 
        figure(hf), hold on,
    end
    plot3([markerPosFRS(:,1)'; sourcePosFRS(1)], ...
        [markerPosFRS(:,2)'; sourcePosFRS(2)],...
        [markerPosFRS(:,3)'; sourcePosFRS(3)], getColorStyle{index}); hold on,
    plot3(sourcePosFRS(1), sourcePosFRS(2), sourcePosFRS(3),'kd', 'MarkerSize', 4, 'LineWidth',2);
    text (sourcePosFRS(1), sourcePosFRS(2), sourcePosFRS(3), axisName);
    
    %Place the image in a 3D plot:
    showImage3d(dr(end:-1:1,:), geometry, panel);
    
end

if size(sourcePositions,1)==2 && size(markerPositions,1)==2,
    [reconstructed_point, distances] = lineIntersect3D(sourcePositions, markerPositions);
    disp('Pt reconstruction errors [L2-norm, mm] after geometrical correction: ');
    absError = sqrt(sum((reconstructed_point - p_gt).*(reconstructed_point - p_gt))),
    relError = sum(distances)/2,
    
    figure(hf), hold on,
    plot3( reconstructed_point(:,1), reconstructed_point(:,2), reconstructed_point(:,3), 'g+');
    box on, grid on, axis equal;
else
    disp('2 source and 2 marker positions must be given. Aborting.');
    return
end

end

% -----------------------------------------
function showImage3d(image, geometry, panel)

halfSizeX = panel.Dimension(1)/2*panel.Pixel;
halfSizeY = panel.Dimension(2)/2*panel.Pixel;

c1 = [-halfSizeX, -halfSizeY, 0, 1]';
c2 = [halfSizeX, -halfSizeY, 0, 1]';
c3 = [-halfSizeX, halfSizeY, 0, 1]';
c4 = [halfSizeX, halfSizeY, 0, 1]';


if strcmp(geometry.config,'gantry')
    
    angle = geometry.gantryAngle + geometry.gantryAngleOffset;
    aid = geometry.AID;
    Vec = geometry.Vec;
    M = roll(angle,[0,0,0])*rot(Vec(9),[0,0,0])*trans(Vec(4),Vec(5),Vec(6)-aid)*roll(Vec(8),[0,0,0])*pitch(Vec(7),[0,0,0]);
    
else
    
    M = geometry.detectorT;
    
end

c1 = M*c1;
c2 = M*c2;
c3 = M*c3;
c4 = M*c4;

hold on,
surf([c1(1), c2(1); c3(1), c4(1)],...
    [c1(2), c2(2); c3(2), c4(2)],...
    [c1(3), c2(3); c3(3), c4(3)],...
    'CData',image(1:2:end, 1:2:end),...
    'FaceColor','texturemap');
colormap(gray);

end

