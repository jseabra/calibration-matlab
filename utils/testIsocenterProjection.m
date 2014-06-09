function testIsocenterProjection(fnameRadOne, panelRadOne, geometryRadOne, fnameRadTwo, panelRadTwo, geometryRadTwo, phantom, objOffset)

% raw is read and must be displayed in beam's eye view:
[okRadOne, drRadOne] = readRawImage(fnameRadOne, panelRadOne.Dimension, panelRadOne.Phys2RadT, panelRadOne.invertGray);
if (~okRadOne)
    disp('Cannot read >> ', fnameRadOne, '. Aborting.');
    return;
end

[okRadTwo, drRadTwo] = readRawImage(fnameRadTwo, panelRadTwo.Dimension, panelRadTwo.Phys2RadT, panelRadTwo.invertGray);
if (~okRadTwo)
    disp('Cannot read >> ', fnameRadTwo, '. Aborting.');
    return;
end

% beam must consist of: geometry, panel and image
beamOne.geometry = geometryRadOne;
beamOne.panel = panelRadOne;
beamOne.dr = drRadOne;

beamTwo.geometry = geometryRadTwo;
beamTwo.panel = panelRadTwo;
beamTwo.dr = drRadTwo;

% % compute point onto isocenter onto flat panel
% F = computeProjection(phantom, geometry, panel, objOffset);
% 
% hfig = figure;
% I = adjustLUT( I, panel.factor); % panel factor between 0 and 1
% imshow(I, []);
% set(hfig, 'units','normalized','outerposition',[0, 0, 1, 1]);
% grid on;
% xlabel('U_{RAD} (pixels)')
% ylabel('V_{RAD} (pixels)')
% hold on,
% plot(F(1),F(2),'+g','MarkerSize',12);

% phantom = [0,0,100];
% F1 = computeProjection(phantom, geometry, panel, objOffset);
% hold on
% plot([F(1) F1(1)], [F(2) F1(2)], 'g');
pointGT = [0,0,0];


[absErr, relErr,reconstructed_point] = reconstructPoint3D(phantom, objOffset, pointGT, beamOne, beamTwo)




end