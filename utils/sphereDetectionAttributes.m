function spAttrb = sphereDetectionAttributes(sad, aid, pixelSize)
%SPHEREDETECTIONATTRIBUTES Retrieve sphere attributes

%# default:
% % % dratio_min  = 0.7; % relative lower tol ratio of sphere diam
% % % dratio_max  = 1.5; % relative upper tol ratio of sphere diam
% % % mf = (sad + aid)/sad; % magnification factor
% % % searchSize = 20; % search by given distance around theoretical marker [mm]
% % % spAttrb.beadSize = 1.5; % diameter of metallic sphere [mm]
% % % projDiam = (mf * spAttrb.beadSize)/pixelSize; % projected sphere diameter [pixels]
% % % dMin = dratio_min*projDiam; dMax = dratio_max*projDiam;
% % % sMin = dMin*dMin; sMax = dMax*dMax;
% % % spAttrb.projDiam = projDiam;
% % % spAttrb.rSize = floor(searchSize/pixelSize + 0.5); % region size [pixels]
% % % spAttrb.sMin  = sMin; % min nr elements [pixels]
% % % spAttrb.sMax  = sMax; % max nr elements [pixels]
% % % spAttrb.dMin  = dMin; % min diameter [pixels]
% % % spAttrb.dMax  = dMax; % max diameter [pixels]
% % % spAttrb.ratio = 0.8; % ratio vertical/horizontal axis
% % % spAttrb.min2diff = projDiam/2; % sphere max resolution (pixels)
% % % spAttrb.threshold = 0; % threshold for local binarization
% % % spAttrb.do_filter = 1; % median filter: 0=off, 1=on
% % % spAttrb.ecc = 0.5; % eccentricity
% % % spAttrb.imopen_param = 0.3; % image open parameter (to remove clutter)
% % % spAttrb.detectionMode = 'auto'; % detect spheres: auto, manual


% %# shreveport:
dratio_min  = 0.8;                                             % relative lower tol ratio of sphere diam
dratio_max  = 1.3;                                             % relative upper tol ratio of sphere diam
mf = (sad + aid)/sad; % magnification factor
searchSize = 20;                                               % search by given distance (mm) around theoretical marker
spAttrb.beadSize = 1.5;                                     % diameter of metallic sphere (mm)
projDiam = (mf * spAttrb.beadSize)/pixelSize; % projected sphere diameter [pixels]
dMin = dratio_min*projDiam; dMax = dratio_max*projDiam;
sMin = dMin*dMin; sMax = dMax*dMax;
spAttrb.projDiam = projDiam;
spAttrb.rSize = floor(searchSize/pixelSize + 0.5); % region size [pixels]
spAttrb.sMin  = sMin;                                       % min nr elements (pixels)
spAttrb.sMax  = sMax;                                       % max nr elements (pixels)
spAttrb.dMin  = dMin;                                       % min diameter (pixels)
spAttrb.dMax  = dMax;                                       % max diameter (pixels)
spAttrb.ratio = 0.8;                                        % ratio vertical/horizontal axis
spAttrb.min2diff = projDiam/2;                              % sphere max resolution (pixels)
spAttrb.threshold = -1.5;                                     % threshold for local binarization
spAttrb.do_filter = 1;                                      % median filter: 0=off, 1=on
spAttrb.ecc = 0.8;                                          % eccentricity
spAttrb.imopen_param = 0.3;                                 % image open parameter (to remove clutter)
spAttrb.detectionMode = 'auto',                             % detect spheres: auto, manual

end

