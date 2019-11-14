function [shift, cc, fa1, fa2] = ImRegFFt3D_FindSpots(o,Im1, Im2, CorrThresh, MinSize)
% [shift, cc, f1, ft2] = ImRegFft2(Im1, Im2, CorrThresh)
%
% do image registration via fft convolution, finding match as point of 
% maximum correlation in the unwhitened images
%
% If Im1 and Im2 both come from the same global image, shift
% is the position of Im2's origin - Im1's origin.
%
% Equivalently, shift the vector such that Im2(x-shift) = Im1(x)
% and Im2(x) = Im1(x+shift) [approximately]
%
% No match if correl<CorrThresh, returns [nan nan].
%
% Correlation returned as cc. NOTE if you pass a 2-element
% vector to CorrThresh, the second entry is an extra-stringent threshold it
% uses for offsets of exactly [0 0], which is often obtained spuriously in
% microscope images. Default is [.2 .6]
%
% MinSize is a number of pixels you need to have matching before it can
% give you a good score (used to regularize the correlation)
%
% Direction specifes the direction in which the images overlap. Can be
% 'South' or 'East'.
%
% If instead of a matrix you pass a 2-element cell array for Im1 or Im2,
% this contains the fft and energy arrays, to save time. These are
% optionally returned as fa1 and fa2.
%
% This differs from ImRegFFt3D in that it can only be used for the
% registration step of the pipeline. It uses know direction and approximate
% amount of overlap of the images to restrict the final shift value. 
% Also, had to change single to double for memory purposes.
% 
% Kenneth D. Harris, 9/8/17
% GPL 3.0 https://www.gnu.org/licenses/gpl-3.0.en.html
 

% not tapering images yet but could


if nargin<4; CorrThresh = [.2 .6]; end
if length(CorrThresh)<2; CorrThresh = CorrThresh*[1, 1]; end

if nargin<5
    MinSize = 100;
end

nTries = 13; % how many local maxima to try for CorrThresh before giving up
Graphics = 1;

%if iscell(Im1)
    %sz = size(Im1{1}, :)/2;% /2 because it was zero-padded
%else
    %sz = size(Im1);
%end
sz = size(Im1);
%%
if ~iscell(Im1)
    % convert to double because matlab has all sorts of problems with integer data types
    Im1d = single(Im1);

    % create arrays of z-scored original images
    Im1z = (Im1d - mean(Im1d(:)))/std(Im1d(:));

    % compute total energy in sub-images of different sizes
    % first make indefinite integrals of energy, starting with a zero:
    Cum1 = zeros(sz+1);
    Cum1(2:sz(1)+1,2:sz(2)+1,2:sz(3)+1) = cumsum(cumsum(cumsum(Im1z.^2,1),2),3);
    Cum1 = single(Cum1);
    
    % next find box edges (inclusive), as a function of dy and dx. 0 or sz+1 means
    % no overlap
    Box1Top = [1:sz(1), ones(1,sz(1))]';
    Box1Bot = [sz(1)*ones(1,sz(1)) , 0:(sz(1)-1)]';
    Box1Left = [1:sz(2), ones(1,sz(2))];
    Box1Right = [sz(2)*ones(1,sz(2)) , (0:sz(2)-1)];
    %Up/Down refers to z direction
    Box1Up = [1:sz(3), ones(1,sz(3))];
    Box1Down = [sz(3)*ones(1,sz(3)) , 0:(sz(3)-1)];
    

    % finally, doing the 2d definite integral means a difference of a
    % difference NOT SURE ABOUT THIS AT ALL
    Energy1 = Cum1(Box1Top,Box1Left,Box1Up) + Cum1(Box1Bot+1,Box1Right+1,Box1Up)...
            - Cum1(Box1Top,Box1Right+1,Box1Up) - Cum1(Box1Bot+1,Box1Left,Box1Up)...
            + Cum1(Box1Top,Box1Left,Box1Down+1) + Cum1(Box1Bot+1,Box1Right+1,Box1Down+1)...
            - Cum1(Box1Top,Box1Right+1,Box1Down+1) - Cum1(Box1Bot+1,Box1Left,Box1Down+1);
    
        
    % zero pad them
    Im1zp = zeros(sz*2);
    Im1zp(1:sz(1),1:sz(2),1:sz(3)) = Im1z;
    Im1zp = single(Im1zp);
    
    % Fourier 
    f1 = fftn(Im1zp);
    
    %Clear variables for memory
    clear Box1Top Box1Bot Box1Left Box1Right Box1Up Box1Down Im1 Im1d Im1z Im1zp Cum1
else 
    f1 = Im1{1};
    Energy1 = Im1{2};
end


% now for image 2 - note box computation is different.
if ~iscell(Im2)
    Im2d = single(Im2);
    Im2z = (Im2d - mean(Im2d(:)))/std(Im2d(:));

    Cum2 = zeros(sz+1);
    Cum2(2:sz(1)+1,2:sz(2)+1,2:sz(3)+1) = cumsum(cumsum(cumsum(Im2z.^2,1),2),3);
    Cum2 = single(Cum2);
    
    Box2Top = [ones(1,sz(1)), (sz(1)+1):-1:2]';
    Box2Bot = [sz(1):-1:1, sz(1)*ones(1,sz(1))]';
    Box2Left = [ones(1,sz(2)), (sz(2)+1):-1:2];
    Box2Right = [sz(2):-1:1, sz(2)*ones(1,sz(2))];
    Box2Up = [ones(1,sz(3)), (sz(3)+1):-1:2];
    Box2Down = [sz(3):-1:1, sz(3)*ones(1,sz(3))];

    Energy2 = Cum2(Box2Bot+1,Box2Right+1,Box2Up) + Cum2(Box2Top,Box2Left,Box2Up) ...
            - Cum2(Box2Top,Box2Right+1,Box2Up) - Cum2(Box2Bot+1,Box2Left,Box2Up)...
            + Cum2(Box2Bot+1,Box2Right+1,Box2Down+1) + Cum2(Box2Top,Box2Left,Box2Down+1) ...
            - Cum2(Box2Top,Box2Right+1,Box2Down+1) - Cum2(Box2Bot+1,Box2Left,Box2Down+1);
        
    Im2zp = zeros(sz*2);
    Im2zp(1:sz(1),1:sz(2),1:sz(3)) = Im2z;
    Im2zp = single(Im2zp);
    
    clear Box2Top Box2Bot Box2Left Box2Right Box2Up Box2Down Im2 Im2d Im2z Cum2
    
    f2 = fftn(Im2zp);
    
    clear Im2zp    
    
else
    f2 = Im2{1};
    Energy2 = Im2{2};
end

% convolve
Conv = ifftn(f1 .* conj(f2));
clear f1 f2

% compute correlation for each shift
Correl = (Conv./(MinSize + sqrt(Energy1.*Energy2)));
%clear Conv
%Take real parts if Correl is complex
if ~isreal(Correl)
    Correl = real(Correl);
end

[cc, MaxShift] = max(Correl(:));

[dy0, dx0,dz0] = ind2sub(size(Conv), MaxShift);
ShiftTry = mod([dy0, dx0, dz0] +sz, sz*2) - sz - 1;


%Use known constrains of overlap to find shift
%NEED TO PUT THE DIVIDING VALUES INTO ISS OBJECT
if abs(ShiftTry(1))<o.MaxRoundShift(1) && abs(ShiftTry(2))<o.MaxRoundShift(1) && abs(ShiftTry(3))<o.MaxRoundShift(2)
    shift = ShiftTry;
else
    shift = [NaN, NaN, NaN];
    [sorted, order] = sort(Correl(:), 'descend');
    i=1;
    while i<=floor(o.ToTest*size(sorted,1))
        [dy0, dx0,dz0] = ind2sub(size(Conv), order(i));
        ShiftTry = mod([dy0, dx0, dz0] +sz, sz*2) - sz - 1;
        if abs(ShiftTry(1))<o.MaxRoundShift(1) && abs(ShiftTry(2))<o.MaxRoundShift(1) && abs(ShiftTry(3))<o.MaxRoundShift(2)
            cc = sorted(i);
            if cc>CorrThresh(1)
                shift = ShiftTry;
                i = size(sorted,1)+1;
                clear sorted order
                break
            end
        end
        i=i+1;
    end
    clear sorted order
end
    


if o.Graphics == 2
    plotCorrelation(abs(Correl),shift,'Correlation');
end

% optional pre-computation outputs:
%if nargout>=3
%    fa1 = {f1; Energy1};
%end
%
%if nargout>=4
%    fa2 = {f2; Energy2};
end