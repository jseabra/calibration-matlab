function [OptVec, Fout, resnorm, idx, MINI] = computeDeformation(phantom, objOffset, panel, measuredSpheres, geometry, optimizerType, bool, displayMode)
% COMPUTEDEFORMATION Estimates gantry deformation parameters using inverse
% projector optimization method
% INPUTS
%    
% OUTPUTS
%     OptVec = deformation vector (Sx*,Sy*,Sz*, Px*,Py*,Pz*, Rx!,Ry!,Rz!) *= in mm; != in deg 
%     Fout = projection at solution (OptVec) in pixels (X,Y)
%     resnorm = sum of distances (non squared) between measured markers and estimated projection
%              (in pixels) as given by the optimizer
%     idx = indices of detected spheres matching nominal projection; use Fout(idx) 
%     MINI = distance (not squared!) between measured markers and estimated projection
%           (in pixels) as calculated doing mean(M - Fout)
% 
% Remark: if there are outliers, the gof would be strongly affected. 
% Suggestion: check MINI (if there is a distance clearly higher than
% others, it means the sphere detection detected some wrong point that has
% no matching)
% remove the point: MINI(idx_of_outlier) = [];
% then, calculate average distance mean(MINI)
% 
% AUTHORS
% JSE
% ========================================================================
disp('Estimating geometrical deformation ... ');
addpath(genpath('./iOptim/')); % import libraries

Vec = geometry.Vec;

% compute position of fiducials:
fiducials = computeProjection(phantom, geometry, panel, objOffset);

%% Do rough estimation of FP displacement X/Y
% if the detector displacement is to be estimated we perform this step:
if bool(4)==1 || bool(5)==1,
    medInit = detection_analysis(measuredSpheres, fiducials);
    % update vector:
    Vec(4) = Vec(4) - medInit(1)*panel.Pixel;
    Vec(5) = Vec(5) - medInit(2)*panel.Pixel;
    
end

%% Set Upper and Lower Constraints:
tolDispl = 100;
tolRot = 10;
tmp = Vec;      
tmp_ub = [tmp(1:6)+tolDispl tmp(7:9)+tolRot];
tmp_lb = [tmp(1:6)-tolDispl tmp(7:9)-tolRot];
UB = tmp_ub(bool);
LB = tmp_lb(bool);

%% Set the goodness of fit function:
[Variable , Constant] = SplitVec(Vec , bool);
myfun = @(x)computeGof(x, Constant, geometry , bool, phantom , objOffset, panel , measuredSpheres);
myfun2 = @(x)computeGof2(x, Constant, geometry , bool, phantom , objOffset, panel , measuredSpheres);
myfun3 = @(x)computeGofMinDis(x , Constant, geometry , bool, phantom , objOffset, panel, measuredSpheres);

switch optimizerType,
    case 'simplex'
        
        options = fminsearch('defaults');
        options = optimset(options, 'MaxFunEvals', 1e5, 'MaxIter',1e6, 'TolFun',1e-14, 'TolX',1e-14);
        [OptVecVar, resnorm, flg, output] = fminsearch(myfun2, Variable, options); % Invoke optimizer
        
    case 'cmaes'
        
        options = fmincmaes('defaults');
        options.PopulationSize = 1e2;
        options = optimset(options, 'MaxFunEvals', 1e5, 'MaxIter',1e5, 'TolFun',1e-14, 'TolX',1e-14);
        constraints.min = LB;
        constraints.max = UB;
        [OptVecVar, resnorm, exitflg, output] = fmincmaes(myfun, Variable, options, constraints); % Invoke optimizer

    case 'lsqnonlin'
        
        options = lsqnonlin('defaults');
        options = optimset(options, 'MaxFunEvals', 1e4, 'MaxIter',1e4, 'TolFun',1e-14, 'TolX',1e-14,'PlotFcns', @optimplotresnorm);
        [OptVecVar, resnorm] = lsqnonlin(myfun3, Variable, LB, UB, options); % Invoke optimizer
        
    case 'levmar'
        
        options = lsqnonlin('defaults');
        options = optimset('Algorithm', 'levenberg-marquardt');
        options = optimset(options, 'MaxFunEvals', 1e4, 'MaxIter',1e4, 'TolFun',1e-13, 'TolX',1e-13);
        [OptVecVar, resnorm] = lsqnonlin(myfun, Variable, [], [], options); % Invoke optimizer
        
    case 'cmb'
        
        options = fmincmaes('defaults');
        options.PopulationSize = 1e2;
        options = optimset(options, 'MaxFunEvals', 1e5, 'MaxIter',1e5, 'TolFun',1e-14, 'TolX',1e-14);
        constraints.min = LB;
        constraints.max = UB;
        [tmpVec, resnorm] = fmincmaes(myfun, Variable, options, constraints); % Invoke optimizer
        
        options = lsqnonlin('defaults');
        options = optimset(options, 'MaxFunEvals', 1e6, 'MaxIter',1e6, 'TolFun',1e-15, 'TolX',1e-15);%,'PlotFcns', @optimplotresnorm);
        
        [OptVecVar1, resnorm1] = lsqnonlin(myfun, tmpVec, LB, UB, options); % Invoke optimizer
        
        % if we optimize Rx and Ry, we may want to check -Rx and -Ry
        if bool(7) == 1 && bool(8) == 1,
            
            if sum(bool) >= 6,         tmpVec(5) = -tmpVec(5); tmpVec(6) = -tmpVec(6);
            elseif sum(bool) >= 4,     tmpVec(3) = -tmpVec(3); tmpVec(4) = -tmpVec(4);
            end
            [OptVecVar2, resnorm2] = lsqnonlin(myfun, tmpVec, LB, UB, options);   % Invoke optimizer
            
            if resnorm1 < resnorm2,
                OptVecVar = OptVecVar1;
                resnorm = resnorm1;
            else
                OptVecVar = OptVecVar2;
                resnorm = resnorm2;
            end
            
        else
            
            resnorm = resnorm1;
            OptVecVar = OptVecVar1;
            
        end
        
end

% Build projection from solution:
OptVec = BuildVec(bool, OptVecVar, Constant);

% indices for labeling:
[MINI, idx] = computeGofMinDis(OptVecVar , Constant , geometry , bool , phantom , objOffset, panel, measuredSpheres);

% update fiducials with estimated deformation vector:
geometry.Vec = OptVec;
Fout = computeProjection(phantom, geometry, panel, objOffset);

% plot only matched points (colormap by score (L2 norm in pixels):
Ftmp = Fout(idx,:);
cmap = [];
caux = [];
for i=1:size(Ftmp,1)
    
     caux = [caux, {['id',num2str(idx(i)),' ',num2str(MINI(i)*panel.Pixel)]}];
    
    if MINI(i)<=1
        cmap = [cmap, {'+g'}];
    elseif MINI(i)>1 && MINI(i)<=2
        cmap = [cmap, {'+c'}];
    elseif MINI(i)>2 && MINI(i)<=4
        cmap = [cmap, {'+b'}];
    elseif MINI(i)>4 && MINI(i)<=7
        cmap = [cmap, {'+y'}];
    elseif MINI(i)>7 && MINI(i)<=10
        cmap = [cmap, {'+m'}];
    else
        cmap = [cmap, {'+r'}];
    end
end
hold on,
for i = 1:size(Ftmp,1)
    scatter(Ftmp(i,1), Ftmp(i,2) ,cmap{i});
end

hold on, legend(caux, 'Location', 'EastOutside');
fprintf('Done. Norm of residuals = %3.2fmm. \n', resnorm/size(Fout,1)*panel.Pixel);


end

% -----------------------------------------
function [med, da] = detection_analysis(M_in, F_in)

[da d idx dis] = associate_points(M_in, F_in);

mask=(da<1e3);
idx=idx(mask);
d=d(mask,:);
da=da(mask);
da = da';
Fq=F_in(mask,:);

med(:,1) = median(d(:,1));
dstd(:,1) = std(d(:,1));

med(:,2) = median(d(:,2));
dstd(:,2) = std(d(:,2));

end

% -----------------------------------------
function [da d idx gof] = associate_points(M_in, F_in)
% =========================================================================
% ASSOCIATE_POINTS:  associates closest pairs of measured and theoretical
% points.
% Note: the algorithm does not consider a unique associate between a given
% pair of points. More than one theoretical point can be associated with
% the same measured one.
%
% INPUTS
% M_in: list of measured points (X,Y)
% F_in: list of theoretical points (X,y)
%
% OUTPUTS
% da: vector with absolute distances between associated points
% d: vector (X,Y) with distances between associated points
% idx: indices of the measured sphere idx(i) associated with each theoretical
% point i
% gof:
% AUTHOR: 
% JSEA
% 30/08/2012: created
% 03/09/2012: updated such that there can be less measured points than
% theoretical ones. The opposite is not possible.
% =========================================================================
sM = size(M_in,1);
sF = size(F_in, 1);

Mx = repmat(M_in(:,1),[1 sF]);
Fx = repmat(F_in(:,1)',[sM 1]);
D(:,:,1) = Mx-Fx;

My = repmat(M_in(:,2),[1 sF]);
Fy = repmat(F_in(:,2)',[sM 1]);
D(:,:,2) = My-Fy;

Da = sqrt(D(:,:,1).^2+D(:,:,2).^2);

[da idx] = min(Da,[],1);

if sF > sM,
    
    [dd,cols] = sort(da);
    for i=1:sF - (sF-sM),
        d(i,:) = D(idx(cols(i)),cols(i),:);
    end
    da = da(cols(1:end-(sF-sM)));
    idx = idx(cols(1:end-(sF-sM)));
    
else
    
    for i=1:sF
        d(i,:) = D(idx(i),i,:);
    end
end

gof = da;                  % works for CMAES, LSQNONLIN, ...
%gof = d(:,1) + d(:,2);    %does not work with CMAES!
%gof = [d(:,1); d(:,2)];   %does not work with CMAES!
end

% -----------------------------------------
function [gof, d3, d1] = computeGof(Variable, Constant, geometry, bool, phantom, objOffset, panel, measuredSpheres)

Vec = BuildVec(bool , Variable , Constant);
geometryTmp = geometry;

% check that the SDD is kept within a certain range
% while sum of Tz and Pz is more than 100mm, reduce such parameters to half
while(abs(Vec(3)-Vec(6))>100)
    %disp('Reducing current SDD...');
    Vec(3) = .75*Vec(3);
    Vec(6) = .75*Vec(6);
end
    
geometryTmp.Vec = Vec;
fiducials = computeProjection(phantom, geometryTmp, panel, objOffset);
[d1, d2, d3, gof] = associate_points(measuredSpheres, fiducials);

end

% -----------------------------------------
function gof = computeGof2(Variable, Constant, geometry, bool, phantom, objOffset, panel, measuredSpheres)

Vec = BuildVec(bool , Variable , Constant);
geometryTmp = geometry;

% check that the SDD is kept within a certain range
% while sum of Tz and Pz is more than 100mm, reduce such parameters to half
while(abs(Vec(3)-Vec(6))>100)
    %disp('Reducing current SDD...');
    Vec(3) = .75*Vec(3);
    Vec(6) = .75*Vec(6);
end

geometryTmp.Vec = Vec;
fiducials = computeProjection(phantom, geometryTmp, panel, objOffset);
[d1, d2, d3, dis] = associate_points(measuredSpheres, fiducials);
gof = sum(dis);

end

% -----------------------------------------
function [MINI, matchpair] = computeGofMinDis(Variable , Constant , geometry , bool , phantom , objOffset, panel , measuredSpheres)

Vec = BuildVec(bool , Variable , Constant);
geometryTmp = geometry;
geometryTmp.Vec = Vec;
fiducials = computeProjection(phantom, geometryTmp, panel, objOffset);

sF = size(fiducials);
sM = size(measuredSpheres);

MINI = [];
matchpair = [];

%Build the distance table:
for idx = 1:sM(1)
    m = measuredSpheres(idx,:);
    m = repmat(m , sF(1) , 1);
    disTmp = sum((m - fiducials).^2,2);
    dis(idx,:) = disTmp';
end %for

MINI = [];
matchpair = zeros(sM(1),1);

maxdis = Inf;
% create the right association between closest pairs fiducial-marker
while sum( matchpair~=0 ) ~= min([sF(1) sM(1)]),    
   
    mindis = min(min(dis));   
    [wminX , wminY] = find(dis == mindis);
    for k = 1:numel(wminX)
    [dis , MINI , matchpair] = AssociateMinDist(dis , MINI , matchpair , wminX(k) , wminY(k) );
    end
    maxdis = max(max(dis));
    
end %for
MINI = sqrt(MINI);
end

% -----------------------------------------
function [dis , MINI , matchpair] = AssociateMinDist(dis , MINI , matchpair , wminX , wminY )
MaxDis = max(max(dis));
%Record the association
matchpair(wminX) = wminY;
MINI(wminX) = dis(wminX,wminY);
%Remove these two points from the list for next loop
dis(wminX,:) = MaxDis;
dis(:,wminY) = MaxDis;
end

% -----------------------------------------
function [Variable , Constant] = SplitVec(Vec , bool)

Constant = Vec(find(~bool));
Variable = Vec(find(bool));

end

% -----------------------------------------
function Vec = BuildVec(bool, Variable , Constant)

lbool = length(bool);
Vec = zeros(1,lbool);

wZero = find(bool == 0);
wOnes = find(bool == 1);

Vec(wZero) = Constant;
Vec(wOnes) = Variable;
end