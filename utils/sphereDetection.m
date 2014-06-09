function [Mak, isValid, matchpair, deltaAID] = sphereDetection(I, phantom, objOffset, panel, geometry, spAttrb, display_mode)
% SPHEREDETECTION detects spheres (auto/manual) in DRs of the calibration cylinder
Mak = [];
isValid = 0;
matchpair = [];
deltaAID = [];

if( ~isSpAttrbValid(spAttrb) )
    isValid = 0;
    return;
end

% compute position of fiducials:
F = computeProjection(phantom, geometry, panel, objOffset);
SF = size(F);

% if (display_mode>=1), %plot image
hfig = figure;
I = adjustLUT( I, panel.factorLow, panel.factorHigh); % panel factor between 0 and 1
imshow(I, []);
set(hfig, 'units','normalized','outerposition',[0, 0, 1, 1]);
grid on;
xlabel('U_{RAD} (pixels)')
ylabel('V_{RAD} (pixels)')
hold on,
for i = 1:SF(1)
    plot(F(i,1), F(i,2) ,'k+'); 
    strTxt = num2str(i);
    text (F(i,1)+10, F(i,2)+10,strTxt,'color','k');
end
drawnow;
hold on;
% end

disp('Finding calibration markers ... ');

% determines the position of the centroid of each object:
MakTMP = [];
SphTMP = [];
index = 1;
for M = 1:SF(1)
    if display_mode>1,
        fprintf('Sphere %d \n',M)
        fprintf('==========\n');
    end
    [PosX, PosY, sphDiam] = getObjectsCentroidPosition(I, F(M,:), spAttrb, display_mode);

    SphTMP = [SphTMP, sphDiam];
    for j = 1:length(PosX)
        MakTMP(index,1) = PosX(j);
        MakTMP(index,2) = PosY(j);
        index = index + 1;
    end
end %M

% Remove the points that are identified twice:
if ~isempty(MakTMP),
    Mak = RemoveDoublet(MakTMP,spAttrb.min2diff);
    sMak = size(Mak);
else
    Mak = [];
end

if display_mode>1,
    fprintf('No. spheres before removing doublets: %d \n', size(MakTMP,1));
    fprintf('No. spheres after removing doublets: %d \n', size(Mak,1));
end

if size(MakTMP,1)==0,
    disp('Algorithm failed to detect markers. Aborting.');
    return
end

% associate measured pts with num proj:
sF = size(F);
sM = size(Mak);
%Build the distance table:
for idx = 1:sM(1)
    m = Mak(idx,:);
    m = repmat(m , sF(1) , 1);
    disTmp = sum((m - F).^2,2);
    dis(idx,:) = disTmp';
end %for
MINI = [];
matchpair = zeros(sM(1),1);
% create the right association between closest pairs fiducial-marker
while sum( matchpair~=0 ) ~= min([sF(1) sM(1)]),
    mindis = min(min(dis));
    [wminX , wminY] = find(dis == mindis);
    for k = 1:numel(wminX)
        [dis , MINI , matchpair] = AssociateMinDist(dis , MINI , matchpair , wminX(k) , wminY(k) );
    end
end %for

Mak_ = Mak;

% if display_mode>=1, %plot image
figure(hfig),
hold on,
for i = 1:sMak(1)
    scatter(Mak_(i,1), Mak_(i,2) ,'w');
    strTxt = num2str(matchpair(i));
    text (Mak_(i,1)+10, Mak_(i,2)+10, strTxt, 'color','w');
end
drawnow;
axis image;
% end

% Remove false positives manually:
if display_mode==1,
    figure(hfig),
    but=1;
    tol = 10;
    title('Detection ok? Remove outliers by pressing right mouse button. left button to exit.');
    while(but~=3)
        [x,y,but] = ginput(1);
        
        pt=[x y];
        dif=repmat(pt,[size(Mak_,1),1])-Mak_;
        dis=Inf;
        for i=1:size(dif,1),
            if(norm(dif(i,:),2)<dis)
                idx=i;
                dis=norm(dif(i,:),2);
            end
        end
        
        if(dis<tol)
            disp('Delete point...');
            Mak_(idx,:)=[NaN NaN];
            
            cla;
            imshow(I, []);
            hold on,
            for i = 1:SF(1)
                plot(F(i,1), F(i,2) ,'+k'); %white circle
                strTxt = num2str(i);
                text (F(i,1)+10, F(i,2)+10,strTxt,'color','k');
            end
            hold on;
            
            for i = 1:size(Mak_,1)
                if ~isnan(Mak_(i,1))
                    scatter(Mak_(i,1), Mak_(i,2) ,'w');
                    strTxt = num2str(matchpair(i));
                    text (Mak_(i,1)+10, Mak_(i,2)+10, strTxt, 'color','w');
                end
            end
            axis image;
            
        end
        
    end
    title(''); drawnow;
end

% Eliminate all outliers from list:
Mak_(find(isnan(Mak_(:,1))),:) = [];
Mak = Mak_;

% Check output and set isValid flag:
if (~isempty(Mak) && size(Mak,1)>10)
    isValid = 1;
    fprintf('Done. Detection rate = %i/%i = %3.1f%% \n', size(Mak,1), size(F,1), size(Mak,1)/size(F,1)*100);
else
    isValid = 0;
    fprintf('Check number of detected spheres! Discarding image for calibration...\n');
end

% if the algorithm detects outliers, we simply remove some spheres from the
% end of the list
if size(Mak,1)>size(F,1),
   
    nel = size(Mak,1)-size(F,1);
    Mak = Mak(1:end-nel,:);
    matchpair = matchpair(1:end-nel,:);
    
end

if ~isempty(SphTMP),    
    deltaAID = (1/spAttrb.beadSize)*(mean(SphTMP)-2)*panel.Pixel*geometry.SAD-geometry.SAD-geometry.AID;
end

end

% -----------------------------------------
function bool = isSpAttrbValid(spAttrb)

if ((spAttrb.rSize <= 0) || (spAttrb.sMin >= spAttrb.sMax)...
        || (spAttrb.dMin >= spAttrb.dMax) || (isempty(spAttrb.ecc))...
        || (spAttrb.do_filter~=0 && spAttrb.do_filter~=1)...
        || (spAttrb.imopen_param<0)...
        || (~strcmp(spAttrb.detectionMode,'auto') && ~strcmp(spAttrb.detectionMode,'manual')))
    
    bool = false;
    disp('Sphere detection attributes are invalid');
else
    bool = true;
end

end

% -----------------------------------------
function LocFin = RemoveDoublet(Location,Thres)

LocFin = [];
LocInd = 1;

sLocation = size(Location);
NumLocation = sLocation(1);

for i = 1:NumLocation
	b = 0;
	if(i+1 <= NumLocation)
		% Loop every position following the current one. Seach for double of the current position
		for j = i+1:NumLocation
			if (SamePoint(Location(i,:), Location(j,:),Thres));
				%fprintf('Found double \n');
				b = 1;
			end %if
		end %j
	end %if
	if (b==0)
		%There is no double => copy this point in the lost
		LocFin(LocInd,:) =  Location(i,:);
		LocInd = LocInd + 1;
	end %if
end %i

end

% -----------------------------------------
function b = SamePoint(itemR, itemT , Thres)
D = norm(itemR-itemT); %Distance between the two points
if(D <= Thres)
	b = 1;
	return;
end %if

b = 0;
end

% -----------------------------------------
function [PosX, PosY, SphDiam] = getObjectsCentroidPosition(I, F, spAttrb, display_mode)

PosX = [];
PosY = [];
SphDiam = [];

% get sphere searching attributes:
detectionMode = spAttrb.detectionMode;
zone  = spAttrb.rSize; % search window size (pixels)
sMin  = spAttrb.sMin; % minimum sphere size (pixels)
sMax  = spAttrb.sMax; % maximum sphere size (pixels)
dMin  = spAttrb.dMin; % minimum sphere diameter (pixels)
dMax  = spAttrb.dMax; % maximum sphere diameter (pixels)
ratio = spAttrb.ratio; % ratio between ax-trans axes
ecc = spAttrb.ecc;
projDiam = spAttrb.projDiam;
threshold = spAttrb.threshold;
do_filter = spAttrb.do_filter;
imopen_param = spAttrb.imopen_param;

if strcmp(detectionMode,'manual') 
    % if detection is manual, force image regions display for point
    % selection
    display_mode = 2;
end

% check that the search window is inside the image:
sI = size(I);
while((F(2)-zone) < 1 || (F(2)+zone) > sI(1) || (F(1)-zone) < 1 || (F(1)+zone) > sI(2))
	%fprintf('Zone out of original image. Chopping. Zone = %d \n',zone)
	zone = zone - 1;
	if zone < 1
		PosX =[];
		PosY =[];
		props = [];
		return;
	end %if
end %while
 
tmpY = floor(F(2)); 
tmpX = floor(F(1)); 
cornerX =tmpX - zone;
cornerY = tmpY - zone;

Icut = I(tmpY-zone:tmpY+zone , tmpX-zone:tmpX+zone);
if(display_mode>1)
    figure(10), subplot(221), imshow(Icut, []); title('ROI centered at fiducial');
end

% filter image and remove background:
It0 = Icut;
if do_filter,
    mv = mean(Icut, 1);
%     It0 = Icut - repmat(mv, size(mv,2), size(mv,1));
%     It1 = It0;
    It1 = It0 - medfilt2(It0,[15 15], 'symmetric');
else
    It1 = It0;
end
It2 = (It1 - mean(It1(:)))/std(It1(:));
if(display_mode>1)   
	figure(10), subplot(222), imshow(It2,[]), title('Filtered'); 
end %endif

% binarize the image by considering objects different than the background:
BCor = It2 < threshold;
se = strel('disk', floor(projDiam*imopen_param), 0);
BCor = imopen(BCor,se);
if(display_mode>1),
   figure(10), subplot(223), imshow(BCor,[]), title('Binarized')
end
L = bwlabel(BCor);
if(display_mode>1),
   figure(10), subplot(224), imshow(label2rgb(L, @jet, [.7 .7 .7])); title('Labeled objects') 
end

if strcmp(detectionMode,'auto')
    
    % get properties from labeled regions:
    props = regionprops(BCor,'Eccentricity', 'Area', 'Centroid', 'BoundingBox');
    num = length(props);
    
    % Compute the centroid of the different objects:

    index = 1;
    % loop for every object:
    for Ob = 1:num
        if(isfield(props,'Area'))
            NPixel = props(Ob).Area; % count number of pixels in the objects
            
            if props(Ob).BoundingBox(4) >= props(Ob).BoundingBox(3),
                Objratio = props(Ob).BoundingBox(3)/props(Ob).BoundingBox(4);
            else
                Objratio = props(Ob).BoundingBox(4)/props(Ob).BoundingBox(3);
            end
            if display_mode>1,
                fprintf('Number of pixels) = %d \n',NPixel);
                fprintf('Dimension u (pixels) = %d \n', props(Ob).BoundingBox(3));
                fprintf('Dimension v (pixels) = %d \n',props(Ob).BoundingBox(4));
                fprintf('Ratio = %d \n', Objratio);
                fprintf('Eccentricity = %d \n', props(Ob).Eccentricity);
            end
            
            % Check that the object matches the sphere criteria:
            if((NPixel <= sMax) && (NPixel >= sMin) && (props(Ob).BoundingBox(3) >= dMin ) &&...
                    (props(Ob).BoundingBox(3) <= dMax) && (props(Ob).BoundingBox(4) >= dMin ) &&...
                    (props(Ob).BoundingBox(4) <= dMax) && (Objratio >= ratio) && props(Ob).Eccentricity<=ecc)
                
                Cxm = props(Ob).Centroid(1); % Horizontal coordinate (X) in PIXELS
                Cym = props(Ob).Centroid(2); % Vertical coordinate (Y) in PIXELS
                
                SphDiam = [SphDiam, props(Ob).BoundingBox(3)];
                    
                if ((Cxm-dMin) < 1 || (Cxm+dMin) > (2*zone) ||...
                        (Cym-dMin) < 1 || (Cym+dMin) > (2*zone))
                    boolC = 0;
                else
                    boolC = 1;
                end %if
                
                if(boolC),
                    
                    PosX(index) = cornerX + Cxm - 1;
                    PosY(index) = cornerY + Cym - 1;
                    
                    if(display_mode>1),
                        figure(10), subplot(221), hold on, scatter(Cxm, Cym, 'r+');
                        fprintf('PosX = %4.3f PosY = %4.3f \n', PosX(index), PosY(index));
                        pause;
                    end
                    index = index + 1;
                else
                    if display_mode>1, fprintf('Background noise. Rejected \n'); end
                end %if
            else
                if display_mode>1, fprintf('Object too big or too small. Rejected \n'); end
            end %if
        else
            if display_mode>1, fprintf('Empty structure. Rejected \n'); end
        end %if
    end	%Ob
    
else % we go on manual selection:
    
    figure(10), subplot(222), hold on,
    [Cxm, Cym] = ginput; % Horizontal and Vertical coordinates (X,Y) in PIXELS
    
    PosX = cornerX + Cxm - 1;
    PosY = cornerY + Cym - 1;
    
    for index=1:numel(PosX),  
        figure(10), subplot(221), hold on, scatter(Cxm(index), Cym(index), 'r+');
        fprintf('PosX = %4.3f PosY = %4.3f \n', PosX(index), PosY(index));
    end
    hold off,
    
end
                   
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
