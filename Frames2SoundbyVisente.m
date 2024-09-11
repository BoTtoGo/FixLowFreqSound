function [pos, dpos, L, corrValues, corrImages] = Frames2SoundbyVisente(vidImages, bkg_removal, varargin)


% 
% [pos, dpos, corrValues, corrImages] = Frames2Sound(vid, bkg_removal)
%
% [pos, dpos, corrValues, corrImages] = Frames2Sound(vid, bkg_removal, noFarting)
%
% [pos, dpos, corrValues, corrImages] = Frames2Sound(vid, bkg_removal, noFarting, subpixelExponential)
%
% [pos, dpos, corrValues, corrImages] = Frames2Sound(vid, bkg_removal, noFarting, subpixelExponential, 'nofigs')
%
% INPUT:
%  vid = 3D matrix, 1 and 2 dimension are the image, and 3 number of frames
%        this can be also a struct array of images contained in 'cdata'
%
%  bkg_removal = 'FALSE' uses images as they are
%                'TRUE' removes a moving average background
% 
%  noFarting = 'FALSE' (default), search correlation in full image
%              'TRUE' , search correlation near the previous peak
%
%  subpixelExponential = if 1 subpixel is calculated with an exponent of
%                           correlation.  (Default = 1)
%
% OUTPUT:
%  pos = 2D vector containing x and y movements
%  dpos = same as pos, but without low freq movements
%  corrValues = max correlation values
%  corrImages = correlation images

%
% LAST MODIF. DATE: 16/10/2015 by Martin Sanz and Javier Garcia
%

TRUE = 1;
FALSE = 0;
DRAW_FIGS = TRUE;
noFarting = FALSE;
subpixelExponential = TRUE;

TrackMaskDefaultRadius=7;
BigTrackMaskDefaultRadius=TrackMaskDefaultRadius+2;


if ischar(vidImages) == 1
    vidImages = VideoReader(vidImages);
end

classType = class(vidImages);
if length(classType) == length('VideoReader')
    if classType == 'VideoReader'
        vidImages = ReadVideoReader(vidImages);
    end
end

if isstruct(vidImages) == 1
    [t1, t2, ~] = size(vidImages(1).cdata);
    aux = zeros(t1,t2,length(vidImages));
    for p=1:length(vidImages)
        aux(:,:,p) = vidImages(p).cdata(:,:,1);
    end
    
    vidImages = aux;
    clear t1;
    clear t2;
    clear aux;
end



vidImages = double(vidImages);


if (ischar(bkg_removal) == TRUE) && (strcmp(bkg_removal,'TRUE') == 1)
    noFarting = TRUE;
end
if isempty(varargin) == 0
    if (length(varargin) > 0)
        if (ischar(varargin{1}) == TRUE) && (strcmp(varargin{1},'TRUE') == 1)
            noFarting = TRUE;
        end
        if (length(varargin) > 1)
            if (ischar(varargin{2}) == TRUE) && (strcmp(varargin{2},'FALSE') == 1)
                subpixelExponential = FALSE;
            end
            if (length(varargin) > 2)
                if (ischar(varargin{3}) == TRUE) && (strcmp(varargin{3},'nofigs') == 0)
                    DRAW_FIGS = FALSE;
                end
            end
        end
    end
end


% Getting size
s1 = size(vidImages,1);
s2 = size(vidImages,2);
Nframes = size(vidImages,3);
L = size(vidImages,3);
range1=(1:s1)-(s1/2+1);
range2=(1:s2)-(s2/2+1);
gw=exp(-(range1/(s1/4)).^2)'*exp(-(range2/(s2/4)).^2); 




% Initialization
pos = zeros(2, Nframes);
dpos = zeros(2, Nframes);
corrValues = zeros(1, Nframes);
corrImages = zeros(s1,s2,Nframes);

cen1= floor(s1/2)+1;
cen2= floor(s2/2)+1;

[X,Y] = meshgrid(-s1/2:s1/2-1,-s2/2:s2/2-1);
R=sqrt(X.^2+Y.^2);



trackMask = double(R<TrackMaskDefaultRadius);

if noFarting == TRUE
    borderTrackMask =trackMask-double(R<(TrackMaskDefaultRadius-1));
    bigTrackMask = double(R<BigTrackMaskDefaultRadius);
    borderBigTrackMask = bigTrackMask - double(R<(BigTrackMaskDefaultRadius-1));
else
    % the process is the same but with all 1s or 0s masks
    trackMask = ones(s1,s2);
    borderTrackMask = zeros(s1,s2);
    bigTrackMask = ones(s1,s2);
    borderBigTrackMask = zeros(s1,s2);
end
                                        
bkgpool = 128;
backgroundImages = BackgroundCalculation(vidImages, bkgpool);

if (bkg_removal==1)
    previousImage=vidImages(:,:,1)-backgroundImages(:,:,1);

else
    previousImage=vidImages(:,:,1);

end;

previousImage = previousImage-mean(mean(previousImage));
previousImage = previousImage/sqrt(sum(sum(previousImage.^2)));
previousImageFT=fftshift(fft2(fftshift(previousImage)));

previousMax = [cen1 cen2];

texto = '';
% --------------------------- MAIN LOOP  -----------------------------------------
for p=2:Nframes
    
    parabolicCalculation = FALSE;
    
    % running display info
    if mod(p,500)==0
        for pp=1:length(texto)
            fprintf('\b');
        end
        texto = sprintf('%i%%',round(p/Nframes*100));
        fprintf('%s',texto);
    end;

    % getting current image
    currentImage=double(vidImages(:,:,p));
    %currentImage(end,end)=currentImage(end-1,end-1);
    
    if (bkg_removal==1)
        currentImage= currentImage - backgroundImages(:,:,p);
    end;
    
    % zero mean. energy normalized
    currentImage = currentImage-mean(mean(currentImage));
    currentImage = currentImage/sqrt(sum(sum(currentImage.^2)));
    
    currentImageFT=fftshift(fft2(fftshift(currentImage)));
    
    % filter applied
   
    fil=conj(previousImageFT);

    if (subpixelExponential == TRUE)
        % This gaussian is a soft low pass that opens the correlation peak
        % speciually important if correlation peak is subpixel size
        fil=fil.*gw;
    end;

    
    previousImageFT = currentImageFT;
    
    % Corr Positive and Negative correlation
    corr=real(fftshift(ifft2(fftshift(currentImageFT.*fil))));
    corrImages(:,:,p) = corr;
      
    % Apply track mask
    corrModif = corr.*circshift(trackMask,[previousMax(1)-cen1 previousMax(2)-cen2]);
    
    % Max point
    [C ,I ] = max(corrModif);
    [~,I1] = max(C);
    
    p1 = I(I1);
    p2 = I1;
    
    % check if the peak moved more than 1/2 window (disaster)
    if (p1-previousMax(1)+cen1 > 0) && (p1-previousMax(1)+cen1 <= s1) && ...
            (p2-previousMax(2)+cen2 > 0) && (p2-previousMax(2)+cen2 <= s2)
        % if the peak is in the border of the track mask we cannot know if
        % it is a local max. BAD. so we default to center of image
        if borderTrackMask(p1-previousMax(1)+cen1, p2-previousMax(2)+cen2) == 1
            % search near the center
            corrModif = corr.*bigTrackMask;
            % Max point
            [C ,I ] = max(corrModif,[],1);
            [~,I1] = max(C);
            
            p1 = I(I1);
            p2 = I1;
            
            % check if it is a inner peak (local max)
            if borderBigTrackMask(p1,p2) == 1
                % default value is zero, if it is not a local max
                previousMax=[cen1 cen2];
                p1 = 0;
                p2 = 0;
                parabolicCalculation = FALSE;
                
            else % good peak
                parabolicCalculation = TRUE;
            end
            
        else % good peak
            parabolicCalculation = TRUE;
        end
    else
        % disaster case. put null shifts
        previousMax=[cen1 cen2];
        p1 = 0;
        p2 = 0;
        parabolicCalculation = FALSE;
    end
    
    if parabolicCalculation == TRUE
        previousMax=[p1 p2];
        
        if (subpixelExponential == TRUE) && (p1 > 1) && (p1 < s1 -1) && (p2 > 1) && (p2 < s2 -1)
            corr(p1-1:p1+1,p2-1:p2+1) = abs(corr(p1-1:p1+1,p2-1:p2+1))- mean(mean(abs(corr)));
            corr(p1-1:p1+1,p2-1:p2+1) = (log(abs(corr(p1-1:p1+1,p2-1:p2+1))));

        end
        
        dp1 = 0;
        dp2 = 0;
        CalcCorr = 1;
        if (p1 > 1) && (p1 < s1 -1)
            % Parabolic max point (1st derivative)
            y1=corr(p1-1,p2);
            y2=corr(p1,p2);
            y3=corr(p1+1,p2);

            dp1=(y3-y1)/2/(2*y2-y1-y3);
        else
            CalcCorr = 0;
        end;
        if (p2 > 1) && (p2 < s2 -1)
            % Parabolic max point (1st derivative)
            y1=corr(p1,p2-1);
            y2=corr(p1,p2);
            y3=corr(p1,p2+1);

            dp2=(y3-y1)/2/(2*y2-y1-y3);
        else
            CalcCorr = 0;
        end;
        if CalcCorr == 1
            corrValues(1,p) = (1-abs(dp1))*(1-abs(dp2))*corr(p1,p2) + abs(dp1)*abs(dp2)*corr(p1+floor(dp1)+ceil(dp1),p2+floor(dp2)+ceil(dp2));
        else
            corrValues(1,p) = 0;
        end

        
        p1=p1-cen1+dp1;
        p2=p2-cen2+dp2;
        
    end
    
    pos(1,p) = p1;
    pos(2,p) = p2;
    
end %for
            
lenconv=51;
lenconv=min([lenconv Nframes-2]);
gw=gausswin(lenconv,3);gw=gw/sum(gw);
poscum=cumsum(pos,2);
posc=conv2(poscum, gw','same');
%HighPass filter

dpos=poscum-posc;
dpos(:,end-lenconv-1:end)=0;

pos = -pos;
dpos = -dpos;

%Plot results
if DRAW_FIGS == TRUE
t=1:Nframes;
    figure(1);plot(t,pos(1,:), t, pos(2,:));title('pos');
    figure(2);plot(t,dpos(1,:), t, dpos(2,:));title('dpos');
end

fprintf('\n\n');

end


function backgroundImages = BackgroundCalculation(vidImages, bkgpool)

% -
% 
% backgroundImages = BackgroundCalculation(vidImages, bkgpool)
%
% INPUT:
%  vid = 3D matrix, 1 and 2 dimension are the image, and 3 number of frames
%  bkgpool = number of images taken to calculate background
%
% OUTPUT:
%  backgroundImages = 3D matrix, 1 and 2 dimension are the image, and 3 number of frames
%

% Getting size
s1 = size(vidImages,1);
s2 = size(vidImages,2);
Nframes = size(vidImages,3);

% Initialization
backgroundImages = zeros(s1, s2, Nframes);
bkgpool = min([bkgpool Nframes-1]);

for p=2:bkgpool+1
    backgroundImages(:,:,p) = sum(vidImages(:,:,1:p-1),3)/length(1:p-1);
    backgroundImages(end,end,p)=backgroundImages(end-1,end-1,p);
end

for p=bkgpool+2:Nframes
    backgroundImages(:,:,p) = backgroundImages(:,:,p-1) + double(vidImages(:,:,p-1))/bkgpool - double(vidImages(:,:,p-bkgpool-1))/bkgpool;
    backgroundImages(end,end,p)=backgroundImages(end-1,end-1,p);
end
 
end

function images = ReadVideoReader(vrObj)
    images = read(vrObj);
    images = squeeze(images(:,:,1,:));
end
