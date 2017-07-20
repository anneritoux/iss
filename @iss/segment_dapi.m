function [o, CellMap] = segment_dapi(o, Dapi0)
% CellMap = o.segment_dapi(DapiIm)
%
% segments a DAPI image and assigns each pixel to a cell. Input is
% optional, otherwise it will load from o.BigDapiFile
%
% only works within region outlined by o.CellCallRegionYX
%
% Output CellMap is same size as input DapiIm, with integer entries for
% each pixel, saying which cell it belongs to. (Zero if unassigned)
%
% also saved to o.CellMapFile

if nargin<2
    Dapi0 = imread(o.BigDapiFile);
end

%% find Cell Calling Region
y0 = min(o.CellCallRegionYX(:,1));
x0 = min(o.CellCallRegionYX(:,2));
y1 = max(o.CellCallRegionYX(:,1));
x1 = max(o.CellCallRegionYX(:,2));

Mask = poly2mask(o.CellCallRegionYX(:,2)-x0+1, o.CellCallRegionYX(:,1)-y0+1, y1-y0+1, x1-x0+1);
Dapi = Dapi0(y0:y1, x0:x1) .*uint16(Mask);

%%
Dapi = imadjust(Dapi); % contrast enhancement
ImSz = size(Dapi);
Debug = 0;
%% threshold the map
ThreshVal = prctile(Dapi(:), o.DapiThresh);

bwDapi = imerode(imfill(Dapi>ThreshVal, 'holes'), strel('disk', 2));

if Debug
    figure(300)
    subplot(2,1,1);
    imagesc(Dapi); 
    subplot(2,1,2);
    imagesc(bwDapi);
    colormap bone
end
%% find local maxima 
dist = bwdist(~bwDapi);
dist0 = dist;
dist0(dist<o.DapiMinSize)=0;
ddist = imdilate(dist0, strel('disk', o.DapiMinSep));
clear dist 
impim = imimposemin(-dist0, imregionalmax(ddist));
clear dist0
if Debug
    figure(301);
    subplot(2,2,1)

    imagesc(dist);
    subplot(2,2,2)
    subplot(2,2,3)
    imagesc(maxxes);
    subplot(2,2,4)
    imagesc(impim);
end
%% segment
% remove pixels at watershed boundaries
bwDapi0 = bwDapi;
bwDapi0(watershed(impim)==0)=0;

% assign all pixels a label
labels = uint32(bwlabel(bwDapi0));
[d, idx] = bwdist(bwDapi0);

% now expand the regions by a margin
CellMap = zeros(ImSz, 'uint32');
Expansions = (d<o.DapiMargin);
CellMap(Expansions) = labels(idx(Expansions));

if Debug
    figure(302)
    subplot(2,2,1);
    image(label2rgb(ws, 'jet', 'w', 'shuffle'));

    subplot(2,2,2);
    imagesc(bwDapiSep);

    subplot(2,2,3);
    imagesc(labels);

    % now give every pixel to its nearest neighbor
    subplot(2,2,4);
    image(colors);
end

o.CellMapFile = fullfile(o.OutputDirectory, 'CellMap.mat');
save(o.CellMapFile, 'CellMap', 'y0', 'y1', 'x0', 'x1');
end

