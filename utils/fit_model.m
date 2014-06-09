function fit_model(fout, geometryCfg, measurements, filename, roomName, place, date, time, system, fpName, axis, fitVec, displayMode, optionPlotFlexmaps)
%========================================================================== 
% FIT_MODEL: 
% Creates a model describing the gantry deformation as a function of the
% gantry angle.
% The model of deformation for each radiographic element is described
% as follows:
%         M = (a0 + a1 t) + a2 cos(a3*t + a4) with {a0, a2} in [mm], a1 in
%         [deg] and a4 in [deg].
% 
% INPUTS
% geometryCfg = 'gantry' or 'free'
% measurements structure: angle, sx,sy,sz, dx,dy,dz, rx,ry,rz (one angle
% per row) or a text file containing same structure
% OUTPUTS
%   xml file whose structure is saved as xmlFilename
%   a csv file containing extrap vals from each model
%
% JSEA 
% 08/08/11: created
% 03/05/12: output given to xml file
% 19/07/12: changed so that only the calibrated parameters are fitted
% instead of all parameters
% 23/08/12: parameters parsed in different way
% 09/09/12: changed parameters parsing again
% 13/09/12: add filter before fitting to remove outliers
% 14/06/13: generate model regardless of geometry (gantry/free)
%========================================================================== 
clc
close all
warning off
addpath(genpath('../utils')); % import libraries

if ~strcmp(geometryCfg, 'gantry') && ~strcmp(geometryCfg, 'free')
    disp('Unknown geometry');
    disp('Aborting the program ...');
    return
end

% least squares model fitting:
options = optimset('TolX',1e-15,'TolFun',1e-15);
% Parametric model:
fitFun = @(x,xdata)(x(1) + x(2)*xdata + x(3)*cos( x(4)*xdata*(pi/180) + x(5)));
param_remove_outliers = 1;
% (0-1) if =1 the thresholds are less strict,
% <1 the thresholds are more and more strict
% and so more pts are considered as outliers
x0 = [0.5 0.5 1 1 0.5];
%x0 = [1 1 1 1 1]; % starting coefficients for model fitting

% optimize =1; take average =0
oSx = fitVec(1); oSy = fitVec(2); oSz = fitVec(3);
oDx = fitVec(4); oDy = fitVec(5); oDz = fitVec(6);
oRx = fitVec(7); oRy = fitVec(8); oRz = fitVec(9);

if isempty(measurements) || size(measurements,2)<9
    measurements = csvread(filename);
    if isempty(measurements) || size(measurements,2)<9
        disp('Could not get calibration measurements from file');
        disp('Aborting the program ...');
        return;
    end
end

if strcmp(geometryCfg, 'free') % free geometry --> take average
  
    const = [0,0,0,0,-1];
    displayMode = 0; % do not display plots at end
    
    % source translation coefficients:
    coeffsSx = [mean(measurements(:,1)), const];
    coeffsSy = [mean(measurements(:,2)), const];
    coeffsSz = [mean(measurements(:,3)), const];
    % detector translation coefficients:
    coeffsDx = [mean(measurements(:,4)), const];
    coeffsDy = [mean(measurements(:,5)), const];
    coeffsDz = [mean(measurements(:,6)), const];
    % detector rotation coefficients:        
    coeffsRx = [mean(measurements(:,7)), const];
    coeffsRy = [mean(measurements(:,8)), const];
    coeffsRz = [mean(measurements(:,9)), const];
 
else 
    
    if size(measurements,1) < 9, % not enough measurements --> take average
        
        const = [0,0,0,0,-1];
        displayMode = 0; % do not display plots at end
        
        % source translation coefficients:
        coeffsSx = [mean(measurements(:,2)), const];
        coeffsSy = [mean(measurements(:,3)), const];
        coeffsSz = [mean(measurements(:,4)), const];
        % detector translation coefficients:
        coeffsDx = [mean(measurements(:,5)), const];
        coeffsDy = [mean(measurements(:,6)), const];
        coeffsDz = [mean(measurements(:,7)), const];
        % detector rotation coefficients:
        coeffsRx = [mean(measurements(:,8)), const];
        coeffsRy = [mean(measurements(:,9)), const];
        coeffsRz = [mean(measurements(:,10)), const];
        
    else

        % duplicate measurements [0:360,0+360:360+360] for continuity:
        angles = measurements(:,1);
        angles = mod(angles+360, 360); % we put the angles between 0 and 360 deg
        params = {'Sx', 'Sy', 'Sz', 'Dx', 'Dy', 'Dz', 'Rx', 'Ry', 'Rz'};
        angles2 = [angles; angles+360];
        [d,idx] = sort(angles2);
        angles3 = angles2(idx);
        for i=1:numel(params),
            eval(['c',params{i},' = [measurements(:,',num2str(i+1),'); measurements(:,',num2str(i+1),')];']);
            eval(['c',params{i},' = c',params{i},'(idx);']); % order for optimization
        end
        
        % remove outliers before fitting:
        for i=1:numel(params),
            eval(['ub = medfilt1(c',params{i},',13) + param_remove_outliers*std(c',params{i},');']);
            eval(['lb = medfilt1(c',params{i},',13) - param_remove_outliers*std(c',params{i},');']);
            eval(['inc',params{i},' = c',params{i},'(c',params{i},' <= ub & c',params{i},' >= lb);']);
            eval(['inAng',params{i},' = angles3(c',params{i},' <= ub & c',params{i},' >= lb);']);
            
        end
        
        for i=1:numel(params),
            
            str = ['if ~o', params{i}, ', disp(''taking avg...''); m', params{i}, '= [mean(measurements(:,',num2str(i+1),')) 0 0 0 0]; resnorm', params{i}, '=0; end'];            
            eval(str);
            str = ['if c', params{i}, '(1)==0 && o',params{i},', m', params{i}, '= [0 0 0 0 0]; resnorm', params{i}, '=0; end'];
            eval(str);
            str = ['if c', params{i}, '(1)~=0 && o',params{i},', [m', params{i}, ', resnorm',params{i},'] = lsqcurvefit(fitFun, x0, inAng',params{i},', inc',params{i},', [], [], options); end'];
            eval(str);
             
        end
        
        for i=1:numel(params),
            eval(['coeffs',params{i},'=[m',params{i},', resnorm',params{i},'];']);
        end
        
    end
        
end

coeffsS = [coeffsSx', coeffsSy', coeffsSz'];
coeffsD = [coeffsDx', coeffsDy', coeffsDz'];
coeffsR = [coeffsRx', coeffsRy', coeffsRz'];

% write model:
writeCalibrationXML(fout, roomName, place, date, time, system, fpName, axis, coeffsD, coeffsR, coeffsS)

if displayMode == 1,
    
    % display:
    xx = 0:10:360;
    % create LUT from model:
    for i = 1:length(params),
        eval(['yy',params{i},' = fitFun( m',params{i},', xx);']);
    end
    
    if strcmp(optionPlotFlexmaps,'cartesian')
        
        figure, title(['X-ray tube ', 'Axis = ',axis]);
        hold on, h1 = plot(inAngSx, incSx, 'bo','MarkerFaceColor','b','MarkerSize',4);
        hold on, h2 = plot(inAngSy, incSy, 'ro','MarkerFaceColor','r','MarkerSize',4);
        hold on, h3 = plot(inAngSz, incSz, 'ko','MarkerFaceColor','k','MarkerSize',4);
        hold on, h4 = plot(xx, yySx, 'b');
        hold on, h5 = plot(xx, yySy, 'r');
        hold on, h6 = plot(xx, yySz, 'k');
        l=legend([h1 h2 h3], {'Sx', 'Sy', 'Sz', 'Location', 'NorthEastOutside'});
        xl=xlabel('Gantry angles [^o]'); yl=ylabel('Displacement [mm]');
        set(gca,'XLim', [0 360]);
        set(l,'Interpreter','latex')
        set(gca,'fontname','Calibri');
        set(xl, 'fontname', 'Calibri','fontsize',14);
        set(yl, 'fontname', 'Calibri','fontsize',14);
        print(gcf,'-dtiff','-r600','srcTransl.tif')
        
        figure, title(['Flat panel ', 'Axis = ',axis]);
        hold on, h1 = plot(inAngDx, incDx, 'bo','MarkerFaceColor','b','MarkerSize',4);
        hold on, h2 = plot(inAngDy, incDy, 'ro','MarkerFaceColor','r','MarkerSize',4);
        hold on, h3 = plot(inAngDz, incDz, 'ko','MarkerFaceColor','k','MarkerSize',4);
        hold on, h4 = plot(xx, yyDx, 'b');
        hold on, h5 = plot(xx, yyDy, 'r');
        hold on, h6 = plot(xx, yyDz, 'k');
        l=legend([h1 h2 h3], {'Dx', 'Dy', 'Dz', 'Location', 'NorthEastOutside'});
        xl=xlabel('Gantry angles [^o]'); yl=ylabel('Displacement [mm]');
        set(gca,'XLim', [0 360]);
        set(l,'Interpreter','latex')
        set(gca,'fontname','Calibri');
        set(xl, 'fontname', 'Calibri','fontsize',14);
        set(yl, 'fontname', 'Calibri','fontsize',14);
        print(gcf,'-dtiff','-r600','fpTransl.tif')
        
        figure, title(['Flat panel ', 'Axis = ',axis]);
        hold on, h1 = plot(inAngRx, incRx, 'bo','MarkerFaceColor','b','MarkerSize',4);
        hold on, h2 = plot(inAngRy, incRy, 'ro','MarkerFaceColor','r','MarkerSize',4);
        hold on, h3 = plot(inAngRz, incRz, 'ko','MarkerFaceColor','k','MarkerSize',4);
        hold on, h4 = plot(xx, yyRx, 'b');
        hold on, h5 = plot(xx, yyRy, 'r');
        hold on, h6 = plot(xx, yyRz, 'k');
        l=legend([h1 h2 h3], {'Rx', 'Ry', 'Rz', 'Location', 'NorthEastOutside'});
        xl=xlabel('Gantry angles [^o]'); yl=ylabel('Rotation [^o]');
        set(gca,'XLim', [0 360]);
        set(l,'Interpreter','latex')
        set(gca,'fontname','Calibri');
        set(xl, 'fontname', 'Calibri','fontsize',14);
        set(yl, 'fontname', 'Calibri','fontsize',14);
        print(gcf,'-dtiff','-r600','fpRot.tif')
        
    elseif strcmp(optionPlotFlexmaps,'polar')
        
        % display:
        xx = 0:10:360;
        % create LUT from model:
        for i = 1:length(params),
            eval(['yy',params{i},' = fitFun( m',params{i},', xx);']);
        end
        
        M_PI_2 = pi/180;
        xx = xx .* M_PI_2;
        %inAngSx = mod(inAngSx+360, 360);
        
        yySx(end) = yySx(1);
        yySy(end) = yySy(1);
        yySz(end) = yySz(1);
        
        MSx = min(yySx); MSx = abs(MSx(1));
        MSy = min(yySy); MSy = abs(MSy(1));
        MSz = min(yySz); MSz = abs(MSz(1));
        [M,idx] = max([MSx, MSy, MSz]);
        
        %figure, %title(['X-ray tube ', 'Axis = ',axis]);
        %hold on, polar(inAngSx.*M_PI_2, incSx, 'bo');%,'MarkerFaceColor','b','MarkerSize',4);
        %hold on, polar(inAngSy.*M_PI_2, incSy, 'ro');%,'MarkerFaceColor','r','MarkerSize',4);
        %hold on, polar(inAngSz.*M_PI_2, incSz, 'ko');%,'MarkerFaceColor','k','MarkerSize',4);
        %polar(xx, yySx + M);
        %hold on, polar(xx, yySy + M);
        %hold on, polar(xx, yySz + M);
        %l=legend([h1 h2 h3], {'Sx', 'Sy', 'Sz', 'Location', 'NorthEastOutside'});
        %xl=xlabel('Gantry angles [^o]'); yl=ylabel('Displacement [mm]');
        %set(l,'Interpreter','latex')
        %set(gca,'fontname','Calibri');
        %set(xl, 'fontname', 'Calibri','fontsize',14);
        %set(yl, 'fontname', 'Calibri','fontsize',14);
        %
        %     figure, title(['Flat panel ', 'Axis = ',axis]);
        %     hold on, h1 = polar(inAngDx.*M_PI_2, incDx, 'bo','MarkerFaceColor','b','MarkerSize',4);
        %     hold on, h2 = polar(inAngDy.*M_PI_2, incDy, 'ro','MarkerFaceColor','r','MarkerSize',4);
        %     hold on, h3 = polar(inAngDz.*M_PI_2, incDz, 'ko','MarkerFaceColor','k','MarkerSize',4);
        %     hold on, h4 = plot(xx.*M_PI_2, yyDx, 'b');
        %     hold on, h5 = plot(xx.*M_PI_2, yyDy, 'r');
        %     hold on, h6 = plot(xx.*M_PI_2, yyDz, 'k');
        %     l=legend([h1 h2 h3], {'Dx', 'Dy', 'Dz', 'Location', 'NorthEastOutside'});
        %     xl=xlabel('Gantry angles [^o]'); yl=ylabel('Displacement [mm]');
        %     set(gca,'XLim', [0 360]);
        %     set(l,'Interpreter','latex')
        %     set(gca,'fontname','Calibri');
        %     set(xl, 'fontname', 'Calibri','fontsize',14);
        %     set(yl, 'fontname', 'Calibri','fontsize',14);
        %
        %     figure, title(['Flat panel ', 'Axis = ',axis]);
        %     hold on, h1 = polar(inAngRx.*M_PI_2, incRx, 'bo','MarkerFaceColor','b','MarkerSize',4);
        %     hold on, h2 = polar(inAngRy.*M_PI_2, incRy, 'ro','MarkerFaceColor','r','MarkerSize',4);
        %     hold on, h3 = polar(inAngRz.*M_PI_2, incRz, 'ko','MarkerFaceColor','k','MarkerSize',4);
        %     hold on, h4 = plot(xx.*M_PI_2, yyRx, 'b');
        %     hold on, h5 = plot(xx.*M_PI_2, yyRy, 'r');
        %     hold on, h6 = plot(xx.*M_PI_2, yyRz, 'k');
        %     l=legend([h1 h2 h3], {'Rx', 'Ry', 'Rz', 'Location', 'NorthEastOutside'});
        %     xl=xlabel('Gantry angles [^o]'); yl=ylabel('Rotation [^o]');
        %     set(gca,'XLim', [0 360]);
        %     set(l,'Interpreter','latex')
        %     set(gca,'fontname','Calibri');
        %     set(xl, 'fontname', 'Calibri','fontsize',14);
        %     set(yl, 'fontname', 'Calibri','fontsize',14);
        %
    end

end

end