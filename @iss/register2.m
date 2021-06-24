function o=register2(o)
% o=iss_register2(o)
%
% register images based on tile files
% creates arrays o.RefPos(y, x): origin of tile (y,x) in pixels on
% reference round relative to global coordinate frame
% 
% o.RelativePos(r, 1:2, t1, t2): origin of tile t2 on 
% reference round minus origin of tile t1 round r. In other words,
% Im2(x;rr) = Im1(x + RelativePos; r). Nan if not defined
% t1 and t2 are linear indices for the tile (y,x)
%
% also makes a global DAPI image 
%
% This finds shifts between tiles using point cloud not by finding the max
% correlation between images
% 
% Kenneth D. Harris, 29/3/17
% GPL 3.0 https://www.gnu.org/licenses/gpl-3.0.en.html
 
%% Logging

if o.LogToFile
    diary(o.LogFile);
    cleanup = onCleanup(@()diary('off'));
end

%% basic variables
rr = o.ReferenceRound;
NonemptyTiles = find(~o.EmptyTiles)';
if size(NonemptyTiles,2)==1
    NonemptyTiles = NonemptyTiles';
end

[nY, nX] = size(o.EmptyTiles);
nTiles = nY*nX;

%% loop through all tiles, finding spots in anchor channel on ref round
o.RawLocalYXZ = cell(nTiles,1);  % cell array, giving spots in local coordinates
o.RawIsolated = cell(nTiles,1);
%SE = fspecial3('ellipsoid',[o.SmoothSize*o.Zpixelsize/o.XYpixelsize, o.SmoothSize*o.Zpixelsize/o.XYpixelsize, o.SmoothSize]);
%SE = fspecial3('ellipsoid',o.SmoothSize);

fprintf('\nDetecting reference spots in tile   ');
for t=NonemptyTiles
    if t<10
        fprintf('\b%d',t);
    else
        fprintf('\b\b%d',t);
    end 
    
    [y,x] = ind2sub([nY nX], t);
    AnchorIm = o.load_3D(rr,y,x,o.AnchorChannel)-o.TilePixelValueShift;
    if o.SmoothSize   
        SE = fspecial3('ellipsoid',o.SmoothSize);
        AnchorImSm = imfilter(AnchorIm, SE);
    else
        AnchorImSm = AnchorIm;
    end
    [tRawLocalYXZ, tRawIsolated] = o.detect_spots(AnchorImSm,t,o.AnchorChannel,rr);
    %Exclude first z-plane as sometimes contains loads of spots.
    Use = tRawLocalYXZ(:,3)~=1;
    o.RawLocalYXZ{t} = tRawLocalYXZ(Use,:);
    o.RawIsolated{t} = tRawIsolated(Use);
end
fprintf('\n');
%% get arrays ready


% WithinTileShift(t,:,r) is origin of tile t round r relative to origin of 
% tile t ref round
%WithinTileShift = nan(nTiles,2,o.nRounds);

% VerticalPairs: n x 2 array of tile IDs
% vShifts: n x 3 array of YXZ shifts
% vScore: n x 1 array of finding shift scores (larger is better)
% HorizontalPairs, hShifts, hScore: similar
VerticalPairs = zeros(0,2);
HorizontalPairs = zeros(0,2);
vShifts = zeros(0,3);
hShifts = zeros(0,3);
vScore = zeros(0,1);
hScore = zeros(0,1);
vChangedSearch = 0;
hChangedSearch = 0;

%% now do the alignments
for t=NonemptyTiles
    [y,x] = ind2sub([nY nX], t);    
    
    % can I align ref round to south neighbor?
    if y<nY && ~o.EmptyTiles(t+1)
        tic
        [shift, score, ChangedSearch] = o.get_initial_shift2(o.RawLocalYXZ{t}, o.RawLocalYXZ{t+1}, o.RegSearch.South,'Register');
        toc
        if all(isfinite(shift))
            VerticalPairs = [VerticalPairs; t, t+1];
            vShifts = [vShifts; shift];
            vScore = [vScore; score];
            vChangedSearch = vChangedSearch + ChangedSearch;
        end
        %ShowPos(o, y, x, y+1, x, rr, shift);
        fprintf('Tile %d (%d, %d), down: shift %d %d %d, score %f\n', t, y, x, shift, score);
        
        %Change search range after 3 tiles or if search has had to be widened twice (This is for speed).
        if size(vShifts,1) == 3 || (mod(vChangedSearch,2) == 0) && (vChangedSearch>0)
            o = o.GetNewSearchRange_Register(vShifts,'South');
        end

    end
    
    % can I align to east neighbor
    if x<nX && ~o.EmptyTiles(t+nY)
        tic
        [shift, score, ChangedSearch] = o.get_initial_shift2(o.RawLocalYXZ{t}, o.RawLocalYXZ{t+nY}, o.RegSearch.East,'Register');
        toc
        if all(isfinite(shift))
            HorizontalPairs = [HorizontalPairs; t, t+nY];
            hShifts = [hShifts; shift];
            hScore = [hScore; score];
            hChangedSearch = hChangedSearch + ChangedSearch;
        end        
        %ShowPos(o, y, x, y, x+1, rr, shift);
        fprintf('Tile %d (%d, %d), right: shift %d %d %d, score %f\n', t, y, x, shift, score);
        
        %Change search range after 3 tiles or if search has had to be widened twice (This is for speed).
        if size(hShifts,1) == 3 || (mod(hChangedSearch,2) == 0) && (hChangedSearch>0)
            o = o.GetNewSearchRange_Register(hShifts,'East');
        end

    end
                    
end

%Convert z shift back to z pixel units
vShifts = vShifts.*[1,1,o.XYpixelsize/o.Zpixelsize];
hShifts = hShifts.*[1,1,o.XYpixelsize/o.Zpixelsize];


%Set any anomalous shifts to average of all other shifts
%Anomalous if awful score or either shift is an outlier
[vShifts, vOutlier] = o.AmendShifts(vShifts,vScore,'Register');
[hShifts, hOutlier] = o.AmendShifts(hShifts,hScore,'Register');


%Save registration info for debugging
o.RegInfo.VerticalPairs = VerticalPairs;
o.RegInfo.vShifts = vShifts;
o.RegInfo.Scorev = vScore;
o.RegInfo.vChangedSearch = vChangedSearch;
o.RegInfo.vOutlier = vOutlier;
o.RegInfo.HorizontalPairs = HorizontalPairs;
o.RegInfo.hShifts = hShifts;
o.RegInfo.Scoreh = hScore;
o.RegInfo.hChangedSearch = hChangedSearch;
o.RegInfo.hOutlier = hOutlier;
%save(fullfile(o.OutputDirectory, 'o2.mat'), 'o');

%% now we need to solve a set of linear equations for each shift,
% This will be of the form M*x = c, where x and c are both of length 
% nTiles=nY*nX. The t'th row is the equation for tile t. 
% c has columns for y and x coordinates

M = zeros(nTiles, nTiles);
c = zeros(nTiles, 3);
for i=1:size(VerticalPairs,1)
    if isnan(vShifts(i,1)); continue; end
    t1 = VerticalPairs(i,1);
    t2 = VerticalPairs(i,2);
    M(t1,t1) = M(t1,t1)+1;
    M(t1,t2) = M(t1,t2)-1;
    c(t1,:) = c(t1,:) - vShifts(i,:);
    M(t2,t2) = M(t2,t2)+1;
    M(t2,t1) = M(t2,t1)-1;
    c(t2,:) = c(t2,:) + vShifts(i,:);
end

for i=1:size(HorizontalPairs,1)
    if isnan(hShifts(i,1)); continue; end
    t1 = HorizontalPairs(i,1);
    t2 = HorizontalPairs(i,2);
    M(t1,t1) = M(t1,t1)+1;
    M(t1,t2) = M(t1,t2)-1;
    c(t1,:) = c(t1,:) - hShifts(i,:);
    M(t2,t2) = M(t2,t2)+1;
    M(t2,t1) = M(t2,t1)-1;
    c(t2,:) = c(t2,:) + hShifts(i,:);
end

% now we want to anchor one of the tiles to a fixed coordinate. We do this
% for a home tile in the middle, because it is going to be connected; and we set
% its coordinate to a large value, so any non-connected ones can be
% detected. (BTW this is why spectral clustering works!!)
Huge = 1e6;
TileDistFromCenter = abs(mod(0:nTiles-1, nY)-nY/2) + ...
    abs(floor((0:nTiles-1)/nY)-nX/2);
[~, HomeTile] = min(TileDistFromCenter(:)./~o.EmptyTiles(:));
%sub2ind([nY nX], ceil(nY/2), ceil(nX/2));
M(nTiles+1,HomeTile) = 1;
c(nTiles+1,:) = [Huge, Huge, Huge];

Tiny = 1e-4; % for regularization
TileOffset0 = (M+Tiny*eye(nTiles+1, nTiles))\c;

% find tiles that are connected to the home tile 
AlignedOK = (TileOffset0(:,1)>Huge/2);
TileOffset1 = nan(nTiles, 3);
TileOffset1(AlignedOK,:) = TileOffset0(AlignedOK,:)-Huge;

% RefPos(t,1:2) is origin of reference tile
RefPos = bsxfun(@minus,TileOffset1, nanmin(TileOffset1))+1;

% tile origin(t,1:2,r)
o.TileOrigin = zeros(nTiles,3,o.nRounds+o.nExtraRounds);
o.TileOrigin(:,:,rr) =  RefPos;

%%

%save(fullfile(o.OutputDirectory, 'o1.mat'), 'o');



%% now make background image
AnchorOrigin = round(o.TileOrigin(:,1:2,rr));           %Only consider YX coordinates
ZOrigin = round(o.TileOrigin(:,3,rr));                  %To align between Z planes if necessary
MaxTileLoc = max(AnchorOrigin);
MaxZ = ceil((max(ZOrigin) + o.nZ));

o.BigDapiFile = fullfile(o.OutputDirectory, 'background_image.tif');
o.BigAnchorFile = fullfile(o.OutputDirectory, 'anchor_image.tif');
if exist(o.BigDapiFile, 'file')
    SaveDapi = false;
else
    SaveDapi = true;
end
if exist(o.BigAnchorFile, 'file')
    SaveAnchor = false;
else
    SaveAnchor = true;
end

%saving options
options.big=true;
options.append=true;
options.message=false;
options.compress='lzw';

fprintf('\nLoading Dapi Image, Z Plane   ');
for z = 1:MaxZ
    BigDapiIm = zeros(ceil((MaxTileLoc + o.TileSz)), 'uint16');
    BigAnchorIm = zeros(ceil((MaxTileLoc + o.TileSz)), 'uint16');
    if z<10
        fprintf('\b%d',z);
    else
        fprintf('\b\b%d',z);
    end 

    for t=NonemptyTiles
        [y,x] = ind2sub([nY nX], t);
        MyOrigin = AnchorOrigin(t,:);
        FileZ = z-ZOrigin(t)+1;
        
        if ~isfinite(MyOrigin(1)); continue; end
        if FileZ < 1 || FileZ > o.nZ
            %Set tile to 0 if currently outside its area
            LocalDapiIm = zeros(o.TileSz,o.TileSz);
            LocalAnchorIm = zeros(o.TileSz,o.TileSz);
        else
            DapiTiff = Tiff(o.TileFiles{o.ReferenceRound,y,x,o.DapiChannel},'r');
            setDirectory(DapiTiff,FileZ);
            LocalDapiIm = read(DapiTiff);
            close(DapiTiff);
            AnchorTiff = Tiff(o.TileFiles{o.ReferenceRound,y,x,o.AnchorChannel},'r');
            setDirectory(AnchorTiff,FileZ);
            LocalAnchorIm = read(AnchorTiff);
            close(AnchorTiff);
        end
        BigDapiIm(floor(MyOrigin(1))+(1:o.TileSz), ...
            floor(MyOrigin(2))+(1:o.TileSz)) ...
            = imresize(LocalDapiIm, 1);        
        BigAnchorIm(floor(MyOrigin(1))+(1:o.TileSz), ...
            floor(MyOrigin(2))+(1:o.TileSz)) ...
            = LocalAnchorIm;
    end
    if SaveDapi
        saveastiff(uint16(BigDapiIm), o.BigDapiFile, options);
    end
    if SaveAnchor
        saveastiff(uint16(BigAnchorIm), o.BigAnchorFile, options);
    end
end
fprintf('\n');

return
end
%%
function ShowPos(o, y, x, y1, x1, r, shift)
	if all(isfinite(shift))
		Color = 'b';
	else
		Color = 'r';
	end
    %figure(239856); 
    clf; hold on
    plot(o.TilePosYXC(:,2), o.TilePosYXC(:,1), 'k.');
    plot([x x1], [y y1], Color);
    plot(x, y, [Color 'o'], 'markersize', r*3);
    set(gca, 'ydir', 'reverse');
    title(sprintf('Round %d', r));
    drawnow;
end