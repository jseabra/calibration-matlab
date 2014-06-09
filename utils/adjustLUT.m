function adjustedImage = adjustLUT( grayImage, val_lowN, val_highN)
% using =0 has no effect
if val_highN==0,
    adjustedImage=grayImage;
else
    
    grayImageN = grayImage./max(grayImage(:));
    adjustedImage = imadjust(grayImageN,[val_lowN; val_highN],[0 1]).*max(grayImage(:));
    
end
end
