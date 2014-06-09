function [ok, Img] = readRawImage(filename, size, rotType, invert)
%READRAWIMAGE reads raw image and displays it in beam's eye view
%orientation, i.e as if the object is projected onto the flat panel and
%seen from the tube
% INPUTS:
% filename: name of raw file
% size: image size in pixels (X,Y)
% rotType: how to rotate around the IEC 3-axes to put the image in beam's eye
% invert: invert grayscale
% view orientation
% % % IEC 3-axes CS:
% % % X+ left->right
% % % Y+ bot->top
% % % Z+ toward tube
% -----------------------------------------

ok = false;
Img = [];
arch = 'ieee-le';
[fid , msg] = fopen (filename, 'r', arch);
[I, count] = fread (fid, size, 'uint16', 0, arch);
fclose (fid);

if invert
    I = 65535-I;
end

if fid == -1
    disp('Failed to read image');
else
    if (strcmp(rotType, 'rotationX')) % CPO (rotation matrix X 180 deg) / Use for casemate
        Img = rot90(I, 1);
        Img = Img(end:-1:1, :);
        ok = true;
    elseif (strcmp(rotType,'rotationY'))
        Img = rot90(I, -1); % ESSEN (rotation matrix Y 180 deg)
        ok = true;
    elseif (strcmp(rotType,'identity'))
        Img = rot90(I, 1); % counterclockwise rotation (UPENN)
        %Img = Img(end:-1:1, :);
        ok = true;
    elseif (strcmp(rotType,'rotationZ'))
        Img = rot90(I, -1); % clockwise rotation
        ok = true;
    else
        disp('Missing orientation info.');
        ok = false;
    end
    
end

if (ok)
    Img = Img(end:-1:1, :); % flip image in Y, Y positive pointing upwards
end

end