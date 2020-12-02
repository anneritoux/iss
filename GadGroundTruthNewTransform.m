%% Stuff to change
%o.TileFiles
%o.nExtraRounds
%o.FileBase
%o.AutoThresh
%Need AllBaseLocalYX in workspace

o.GadChannel = 4;
o.GadRound = 9;
o.GcampChannel = 6;
o.GcampRound = 9;
o.GadReferenceChannel = o.GadChannel;       %Anchor is channel 7 but want best to Gad channel
o.UseRounds = o.GadRound;
%o.UseChannels = [o.GadChannel,o.GcampChannel];

%% Get Gad spots
[nY, nX] = size(o.EmptyTiles);
nTiles = nY*nX;
o.TileCentre = 0.5*[o.TileSz+1,o.TileSz+1];
NonemptyTiles = find(~o.EmptyTiles)';

o.GadRawLocalYX = cell(nTiles,o.nBP, o.nRounds+o.nExtraRounds);
o.GadRawIsolated = cell(nTiles,o.nBP, o.nRounds+o.nExtraRounds);

for t=NonemptyTiles(:)'
    if mod(t,10)==0; fprintf('Loading tile %d reference image\n', t); end
    FileName = o.TileFiles{o.GadRound,t};
    TifObj = Tiff(FileName);
    for b=o.UseChannels     
        TifObj.setDirectory(o.FirstBaseChannel + b - 1);
        ReferenceIm = int32(TifObj.read())-o.TilePixelValueShift;            
        if o.SmoothSize
            SE = fspecial('disk', o.SmoothSize);
            ReferenceImSm = imfilter(ReferenceIm ,SE);
        else
            ReferenceImSm = ReferenceIm;
        end
        [o.GadRawLocalYX{t,b,o.GadRound}, o.GadRawIsolated{t,b,o.GadRound}] = o.detect_spots(ReferenceImSm,t,b,o.GadRound);
    end
end


%% Align to anchor round
o.GadInfo.Scores = zeros(nTiles,1);
o.GadInfo.ChangedSearch = 0;

for t=NonemptyTiles
    for r = o.UseRounds
        tic
        [o.D0(t,:,r), o.GadInfo.Scores(t),tChangedSearch] = o.get_initial_shift2(o.GadRawLocalYX{t,o.GadReferenceChannel,o.GadRound},...
            vertcat(o.RawLocalYX{t,:}), o.FindSpotsSearch{3},'FindSpots');
        o.GadInfo.ChangedSearch = o.GadInfo.ChangedSearch+tChangedSearch;
        
        fprintf('Tile %d, shift from reference round (%d) to round %d: [%d %d], score %d \n',...
            t, o.ReferenceRound,r, o.D0(t,:,r),o.GadInfo.Scores(t));
        
        %Change search range after 3 tiles or if search has had to be widened twice (This is for speed).
        if t == 3 || (mod(o.GadInfo.ChangedSearch,2) == 0) && (o.GadInfo.ChangedSearch>0)
            o = o.GetNewSearchRange_FindSpots(t,3);
        end
    end
end

[o,x] = o.PointCloudRegister6_Gad(o.GadRawLocalYX, o.RawLocalYX, nTiles);
%Reference round coordinates are adjusted as to take account of chromatic
%aberration.
o.RawLocalYX = x;
clearvars x;


%% Save background images
%Compute approx new shifts from D matrices
YXGadRoundTileShifts = permute(squeeze(o.D(3,:,:,o.GadRound,o.GadReferenceChannel)),[2,1,3]);
o.TileOrigin(:,:,o.GadRound) =  o.TileOrigin(:,:,o.ReferenceRound) - YXGadRoundTileShifts;  

AnchorOrigin = round(o.TileOrigin(:,:,o.GadRound));
MaxTileLoc = max(AnchorOrigin);
NegOriginShift = abs(min(min(AnchorOrigin(:)),1)-1);        %To correct for negative origin - remove later
MissedSpotsShift = abs(min(min(ceil((MaxTileLoc + o.TileSz) - max(o.pxSpotGlobalYX))),0));    %To correct for spots being outside image bounds i.e. zero pad.
BigGadIm = zeros(ceil((MaxTileLoc + o.TileSz + NegOriginShift+MissedSpotsShift)), 'uint16');
BigGcampIm = zeros(ceil((MaxTileLoc + o.TileSz + NegOriginShift+MissedSpotsShift)), 'uint16');

for t=NonemptyTiles
    MyOrigin = AnchorOrigin(t,:)+NegOriginShift;
    if mod(t,10)==0; fprintf('Loading tile %d DAPI image\n', t); end
    if ~isfinite(MyOrigin(1)); continue; end
    LocalGadIm = imread(o.TileFiles{o.GadRound,t}, o.GadChannel);
    BigGadIm(floor(MyOrigin(1))+(1:o.TileSz), ...
        floor(MyOrigin(2))+(1:o.TileSz)) ...
        = imresize(LocalGadIm, 1);
    LocalGcampIm = imread(o.TileFiles{o.GadRound,t}, o.GcampChannel);
    BigGcampIm(floor(MyOrigin(1))+(1:o.TileSz), ...
        floor(MyOrigin(2))+(1:o.TileSz)) ...
        = LocalGcampIm;
end

o.BigGadFile = fullfile(o.OutputDirectory, 'gad_image.tif');
imwrite(BigGadIm(1+NegOriginShift:end,1+NegOriginShift:end), o.BigGadFile);
o.BigGcampFile = fullfile(o.OutputDirectory, 'gcamp_image.tif');
imwrite(BigGcampIm(1+NegOriginShift:end,1+NegOriginShift:end), o.BigGcampFile);


%% Find intensity of all pixel based spots in the Gad round
o.UseRounds = o.GadRound;
[nY, nX] = size(o.EmptyTiles);
nTiles = nY*nX;
o.TileCentre = 0.5*[o.TileSz+1,o.TileSz+1];
NonemptyTiles = find(~o.EmptyTiles)';

nSpots = length(o.pxSpotGlobalYX);
GlobalYX = reshape([o.pxSpotGlobalYX(1:nSpots,2),o.pxSpotGlobalYX(1:nSpots,1)]',1,2,nSpots);
GlobalYX = flip(GlobalYX,2);
% Dist2Tiles = GlobalYX - repmat(o.TileOrigin(:,:,o.ReferenceRound),1,1,nSpots);
% Dist2Tiles(Dist2Tiles<0) = Inf;
% [Dist2Tile,LocalTile] = min(sum(Dist2Tiles,2));
% o.pxLocalTile = squeeze(LocalTile);
pxSpotLocalYX = o.pxSpotGlobalYX-o.TileOrigin(o.pxLocalTile,:,o.ReferenceRound);


ndRoundTile = nan(nSpots,o.nRounds+o.nExtraRounds);
ndRoundYX = nan(nSpots,2,o.nRounds+o.nExtraRounds);
PossNeighbs = [-1 -nY 1 nY 0]; % NWSE then same tile - same will have priority by being last

for r=o.UseRounds
    fprintf('Finding appropriate tiles for round %d\n', r);
    
    for n = PossNeighbs
        % find origins of each tile's neighbor, NaN if not there
        NeighbTile = (1:nTiles)+n;
        NeighbOK = (NeighbTile>=1 & NeighbTile<=nTiles);
        NeighbOrigins = nan(nTiles,2);
        NeighbOrigins(NeighbOK,:) = round(o.TileOrigin(NeighbTile(NeighbOK),:,r));
        
        % now for each spot see if it is inside neighbor's tile area
        SpotsNeighbOrigin = NeighbOrigins(o.pxLocalTile,:);
        SpotsInNeighbTile = all(o.pxSpotGlobalYX>=SpotsNeighbOrigin+1+o.ExpectedAberration...
            & o.pxSpotGlobalYX<=SpotsNeighbOrigin+o.TileSz-o.ExpectedAberration, 2);
        
        % for those that were in set this to be its neighbor
        ndRoundTile(SpotsInNeighbTile,r) = NeighbTile(o.pxLocalTile(SpotsInNeighbTile));    
    end
    
    % compute YX coord
    HasTile = isfinite(ndRoundTile(:,r));
    ndRoundYX(HasTile,:,r) = o.pxSpotGlobalYX(HasTile,:) - round(o.TileOrigin(ndRoundTile(HasTile,r),:,r));
    
end


% loop through all tiles, finding spot colors
ndLocalYX = [pxSpotLocalYX-o.TileCentre,ones(nSpots,1)];
ndSpotColors = nan(nSpots, o.nBP);
o.UseRounds = o.GadRound;
for t=NonemptyTiles
    [y, x] = ind2sub([nY nX], t);
   
    for r=o.UseRounds         
        % find spots whose home tile on round r is t
        MySpots = (ndRoundTile(:,r)==t);
        if ~any(MySpots); continue; end
        
        % open file for this tile/round
        FileName = o.TileFiles{r,t};
        TifObj = Tiff(FileName);
        
        % find the home tile for all current spots in the ref round
        RefRoundHomeTiles = o.pxLocalTile(ndRoundTile(:,r)==t);
        MyRefTiles = unique(RefRoundHomeTiles);
        fprintf('\nRef round home tiles for spots in t%d at (%2d, %2d), r%d: ', t, y, x, r);
        for i=MyRefTiles(:)'
            fprintf('t%d, %d spots; ', i, sum(RefRoundHomeTiles==i));
        end
        fprintf('\n');        
        
        
        % now read in images for each base
        for b=o.UseChannels               %No 0 as trying without using anchor

            
            TifObj.setDirectory(o.FirstBaseChannel + b - 1);
            BaseIm = int32(TifObj.read())-o.TilePixelValueShift;
            
            if o.SmoothSize
                BaseImSm = imfilter(double(BaseIm), fspecial('disk', o.SmoothSize));
            else
                BaseImSm = BaseIm;
            end
            
            for t2 = MyRefTiles(:)'
                MyBaseSpots = (ndRoundTile(:,r)==t & o.pxLocalTile==t2);
                MyLocalYX = ndLocalYX(MyBaseSpots,:);
                
                if t == t2
                    fprintf('Point cloud: ref round tile %d -> tile %d round %d base %d, %d matches, error %f\n', ...
                        t, t2, r, b,  o.GadInfo.nMatches(t,b,r), o.GadInfo.Error(t,b,r));
                    MyPointCorrectedYX = MyLocalYX*o.D(:,:,t,r,b)+o.TileCentre;
                    MyPointCorrectedYX = round(MyPointCorrectedYX);
                    ndSpotColors(MyBaseSpots,b) = IndexArrayNan(BaseImSm, MyPointCorrectedYX');
                else
                    [MyPointCorrectedYX, Error, nMatches] = o.different_tile_transform(o.GadRawLocalYX,o.RawLocalYX, ...
                        MyLocalYX,t,t2,r,b);
                    fprintf('Point cloud: ref round tile %d -> tile %d round %d base %d, %d matches, error %f\n', ...
                        t, t2, r, b,  nMatches, Error);
                    ndSpotColors(MyBaseSpots,b) = IndexArrayNan(BaseImSm, MyPointCorrectedYX');
                end
               
            end    
        end
        TifObj.close();       
    end
end
fprintf('\n');

o.pxGadColor = ndSpotColors;


%% Take Gad channel local maxima as the spots - find SpotColors for these
%Get Anchor round coordinates of Gad Peaks
o.UseRounds = 1:o.nRounds;
o.UseChannels = 1:o.nBP;
rr = o.ReferenceRound;
nAll = sum(sum(cellfun(@numel, o.GadRawIsolated(:,o.GadChannel,o.GadRound))));

AllGlobalYX = zeros(nAll,2);
AllLocalYX = zeros(nAll,2);
AllIsolated = zeros(nAll,1);
OriginalTile = zeros(nAll,1);
AllGadIntensity = zeros(nAll,1);

ind = 1;
for t=NonemptyTiles
    nMySpots = sum(sum(cellfun(@numel, o.GadRawIsolated(t,o.GadChannel,o.GadRound))));
    FileName = o.TileFiles{o.GadRound,t};
    TifObj = Tiff(FileName);
    TifObj.setDirectory(o.FirstBaseChannel + o.GadChannel - 1);
    BaseIm = int32(TifObj.read())-o.TilePixelValueShift;
    AllGadIntensity(ind:ind+nMySpots-1,:) = ...
        IndexArrayNan(BaseIm, o.GadRawLocalYX{t,o.GadChannel,o.GadRound}');
       
    AllLocalYX(ind:ind+nMySpots-1,:) = (o.GadRawLocalYX{t,o.GadChannel,o.GadRound}-...
        o.TileCentre-o.D(3,:,t,o.GadRound,o.GadChannel))/...
        o.D(1:2,:,t,o.GadRound,o.GadChannel)+o.TileCentre;
    AllGlobalYX(ind:ind+nMySpots-1,:) = AllLocalYX(ind:ind+nMySpots-1,:)+o.TileOrigin(t,:,rr);
    OriginalTile(ind:ind+nMySpots-1) = t;
    ind = ind+nMySpots;
end

%Remove duplciates
[AllLocalTile, ~] = which_tile(AllGlobalYX, o.TileOrigin(:,:,rr), o.TileSz);
NotDuplicate = (AllLocalTile==OriginalTile);
ndGlobalYX = AllGlobalYX(NotDuplicate,:);
ndLocalYX = AllLocalYX(NotDuplicate,:);
ndLocalTile = AllLocalTile(NotDuplicate,:);
ndGadIntensity = AllGadIntensity(NotDuplicate);
nnd = sum(NotDuplicate);

% decide which tile to read each spot off in each round. 
% They are read of home tile if possible (always possible in ref round)
% in other rounds might have to be a NWSE neighbor - but never a diagonal
% neighbor
% ndRoundTile(s,r) stores appropriate tile for spot s on round r
% ndRoundYX(s,:,r) stores YX coord on this tile

%Compute approx new shifts from D matrices
YXRoundTileShifts = permute(squeeze(mean(o.D(3,:,:,1:o.nRounds,:),5)),[2,1,3]);
o.TileOrigin(:,:,1:o.nRounds) =  o.TileOrigin(:,:,rr) - YXRoundTileShifts;  

ndRoundTile = nan(nnd,o.nRounds);
ndRoundYX = nan(nnd,2,o.nRounds);

PossNeighbs = [-1 -nY 1 nY 0]; % NWSE then same tile - same will have priority by being last

for r=o.UseRounds
    fprintf('Finding appropriate tiles for round %d\n', r);
    
    for n = PossNeighbs
        % find origins of each tile's neighbor, NaN if not there
        NeighbTile = (1:nTiles)+n;
        NeighbOK = (NeighbTile>=1 & NeighbTile<=nTiles);
        NeighbOrigins = nan(nTiles,2);
        NeighbOrigins(NeighbOK,:) = round(o.TileOrigin(NeighbTile(NeighbOK),:,r));
        
        % now for each spot see if it is inside neighbor's tile area
        SpotsNeighbOrigin = NeighbOrigins(ndLocalTile,:);
        SpotsInNeighbTile = all(ndGlobalYX>=SpotsNeighbOrigin+1+o.ExpectedAberration...
            & ndGlobalYX<=SpotsNeighbOrigin+o.TileSz-o.ExpectedAberration, 2);
        
        % for those that were in set this to be its neighbor
        ndRoundTile(SpotsInNeighbTile,r) = NeighbTile(ndLocalTile(SpotsInNeighbTile));    
    end
    
    % compute YX coord
    HasTile = isfinite(ndRoundTile(:,r));
    ndRoundYX(HasTile,:,r) = ndGlobalYX(HasTile,:) - round(o.TileOrigin(ndRoundTile(HasTile,r),:,r));
    
end

% loop through all tiles, finding spot colors
ndLocalYX = [ndLocalYX-o.TileCentre,ones(nnd,1)];
ndSpotColors = nan(nnd, o.nBP, o.nRounds);
ndPointCorrectedLocalYX = nan(nnd, 2, o.nRounds, o.nBP);

for t=NonemptyTiles
    [y, x] = ind2sub([nY nX], t);
   
    for r=o.UseRounds         
        % find spots whose home tile on round r is t
        MySpots = (ndRoundTile(:,r)==t);
        if ~any(MySpots); continue; end
        
        % open file for this tile/round
        FileName = o.TileFiles{r,t};
        TifObj = Tiff(FileName);
        
        % find the home tile for all current spots in the ref round
        RefRoundHomeTiles = ndLocalTile(ndRoundTile(:,r)==t);
        MyRefTiles = unique(RefRoundHomeTiles);
        fprintf('\nRef round home tiles for spots in t%d at (%2d, %2d), r%d: ', t, y, x, r);
        for i=MyRefTiles(:)'
            fprintf('t%d, %d spots; ', i, sum(RefRoundHomeTiles==i));
        end
        fprintf('\n');        
        
        
        % now read in images for each base
        for b=o.UseChannels               %No 0 as trying without using anchor

            
            TifObj.setDirectory(o.FirstBaseChannel + b - 1);
            BaseIm = int32(TifObj.read())-o.TilePixelValueShift;
            
            if o.SmoothSize
                BaseImSm = imfilter(double(BaseIm), fspecial('disk', o.SmoothSize));
            else
                BaseImSm = BaseIm;
            end
            
            for t2 = MyRefTiles(:)'
                MyBaseSpots = (ndRoundTile(:,r)==t & ndLocalTile==t2);
                MyLocalYX = ndLocalYX(MyBaseSpots,:);
                
                if t == t2
                    fprintf('Point cloud: ref round tile %d -> tile %d round %d base %d, %d/%d matches, error %f\n', ...
                        t, t2, r, b,  o.nMatches(t,b,r), o.RawLocalNo(t2), o.Error(t,b,r));
                    if o.nMatches(t,b,r)<o.MinPCMatchFract*o.AllBaseSpotNo(t,b,r) || isempty(o.nMatches(t,b,r))
                        warning('Tile %d, channel %d, round %d has %d point cloud matches, which is below the threshold of %d.',...
                            t,b,r,o.nMatches(t,b,r),o.MinPCMatchFract*o.AllBaseSpotNo(t,b,r));
                    end
                    MyPointCorrectedYX = MyLocalYX*o.D(:,:,t,r,b)+o.TileCentre;
                    MyPointCorrectedYX = round(MyPointCorrectedYX);
                    ndPointCorrectedLocalYX(MyBaseSpots,:,r,b) = MyPointCorrectedYX;
                    ndSpotColors(MyBaseSpots,b,r) = IndexArrayNan(BaseImSm, MyPointCorrectedYX');
                else
                    [MyPointCorrectedYX, Error, nMatches] = o.different_tile_transform(AllBaseLocalYX,o.RawLocalYX, ...
                        MyLocalYX,t,t2,r,b);
                    fprintf('Point cloud: ref round tile %d -> tile %d round %d base %d, %d/%d matches, error %f\n', ...
                        t, t2, r, b,  nMatches, o.RawLocalNo(t2), Error);
                    if nMatches<o.MinPCMatchFract*o.AllBaseSpotNo(t,b,r) || isempty(nMatches)
                        continue;
                    end
                    ndPointCorrectedLocalYX(MyBaseSpots,:,r,b) = MyPointCorrectedYX;
                    ndSpotColors(MyBaseSpots,b,r) = IndexArrayNan(BaseImSm, MyPointCorrectedYX');
                end
               
            end    
        end
        TifObj.close();       
    end
end
fprintf('\n');

ndSpotColorsToUse = ndSpotColors(:,o.UseChannels,o.UseRounds);
Good = all(isfinite(ndSpotColorsToUse(:,:)),2);
GoodGlobalYX = ndGlobalYX(Good,:);
GoodSpotColors = ndSpotColors(Good,:,:);
GoodLocalTile = ndLocalTile(Good);

%Append to DotProduct/Prob method method
%Append twice as one will be Gad score
o.cSpotColors = [o.cSpotColors;GoodSpotColors;GoodSpotColors];
o.cSpotIsolated = [o.cSpotIsolated;false(size(GoodSpotColors,1),1);false(size(GoodSpotColors,1),1)];
o.SpotCombi = [o.SpotCombi;true(size(GoodSpotColors,1),1);true(size(GoodSpotColors,1),1)];
o.SpotGlobalYX = [o.SpotGlobalYX;GoodGlobalYX;GoodGlobalYX];
o.GadPeakSpots = zeros(size(o.cSpotIsolated));
o.GadPeakSpots(end-2*size(GoodSpotColors,1)+1:end-size(GoodSpotColors,1))=1;
o.GadPeakSpots(end-size(GoodSpotColors,1)+1:end)=2;

o.GadPeakColor = zeros(size(o.cSpotIsolated));
o.GadPeakColor(o.GadPeakSpots==1) = ndGadIntensity(Good);
o.GadPeakColor(o.GadPeakSpots==2) = ndGadIntensity(Good);


%Now need to do call_spots stuff
%call_spots but:
%o.SpotScore(o.GadPeakSpots==2) = SpotScores(o.GadPeakSpots==2,22);
%o.SpotCodeNo(o.GadPeakSpots==2)==22;
%call_spots_prob but:
%LogProbNoSort = LogProb;
%o.pSpotCodeNo(o.GadPeakSpots==2)=22;
%o.pLogProbOverBackground(o.GadPeakSpots==2) = ...
%LogProbNoSort(o.GadPeakSpots==2,22)-BackgroundLogProb(o.GadPeakSpots==2);
%o.pSpotScore(o.GadPeakSpots==2) = ...
%LogProbNoSort(o.GadPeakSpots==2,22)-LogProb(o.GadPeakSpots==2,2);


%% Find intensity of all dotproduct/prob based spots in the Gad round
o.UseRounds = o.GadRound;
nSpots = length(o.SpotGlobalYX);
GlobalYX = reshape([o.SpotGlobalYX(1:nSpots,2),o.SpotGlobalYX(1:nSpots,1)]',1,2,nSpots);
GlobalYX = flip(GlobalYX,2);
Dist2Tiles = GlobalYX - repmat(o.TileOrigin(:,:,o.ReferenceRound),1,1,nSpots);
Dist2Tiles(Dist2Tiles<0) = Inf;
[Dist2Tile,LocalTile] = min(sum(Dist2Tiles,2));
o.pLocalTile = squeeze(LocalTile);
SpotLocalYX = o.SpotGlobalYX-o.TileOrigin(LocalTile,:,o.ReferenceRound);


ndRoundTile = nan(nSpots,o.nRounds+o.nExtraRounds);
ndRoundYX = nan(nSpots,2,o.nRounds+o.nExtraRounds);
PossNeighbs = [-1 -nY 1 nY 0]; % NWSE then same tile - same will have priority by being last

for r=o.UseRounds
    fprintf('Finding appropriate tiles for round %d\n', r);
    
    for n = PossNeighbs
        % find origins of each tile's neighbor, NaN if not there
        NeighbTile = (1:nTiles)+n;
        NeighbOK = (NeighbTile>=1 & NeighbTile<=nTiles);
        NeighbOrigins = nan(nTiles,2);
        NeighbOrigins(NeighbOK,:) = round(o.TileOrigin(NeighbTile(NeighbOK),:,r));
        
        % now for each spot see if it is inside neighbor's tile area
        SpotsNeighbOrigin = NeighbOrigins(o.pLocalTile,:);
        SpotsInNeighbTile = all(o.SpotGlobalYX>=SpotsNeighbOrigin+1+o.ExpectedAberration...
            & o.SpotGlobalYX<=SpotsNeighbOrigin+o.TileSz-o.ExpectedAberration, 2);
        
        % for those that were in set this to be its neighbor
        ndRoundTile(SpotsInNeighbTile,r) = NeighbTile(o.pLocalTile(SpotsInNeighbTile));    
    end
    
    % compute YX coord
    HasTile = isfinite(ndRoundTile(:,r));
    ndRoundYX(HasTile,:,r) = o.SpotGlobalYX(HasTile,:) - round(o.TileOrigin(ndRoundTile(HasTile,r),:,r));
    
end


% loop through all tiles, finding spot colors
ndLocalYX = [SpotLocalYX-o.TileCentre,ones(nSpots,1)];
ndSpotColors = nan(nSpots, o.nBP);

for t=NonemptyTiles
    [y, x] = ind2sub([nY nX], t);
   
    for r=o.UseRounds         
        % find spots whose home tile on round r is t
        MySpots = (ndRoundTile(:,r)==t);
        if ~any(MySpots); continue; end
        
        % open file for this tile/round
        FileName = o.TileFiles{r,t};
        TifObj = Tiff(FileName);
        
        % find the home tile for all current spots in the ref round
        RefRoundHomeTiles = o.pLocalTile(ndRoundTile(:,r)==t);
        MyRefTiles = unique(RefRoundHomeTiles);
        fprintf('\nRef round home tiles for spots in t%d at (%2d, %2d), r%d: ', t, y, x, r);
        for i=MyRefTiles(:)'
            fprintf('t%d, %d spots; ', i, sum(RefRoundHomeTiles==i));
        end
        fprintf('\n');        
        
        
        % now read in images for each base
        for b=o.UseChannels               %No 0 as trying without using anchor

            
            TifObj.setDirectory(o.FirstBaseChannel + b - 1);
            BaseIm = int32(TifObj.read())-o.TilePixelValueShift;
            
            if o.SmoothSize
                BaseImSm = imfilter(double(BaseIm), fspecial('disk', o.SmoothSize));
            else
                BaseImSm = BaseIm;
            end
            
            for t2 = MyRefTiles(:)'
                MyBaseSpots = (ndRoundTile(:,r)==t & o.pLocalTile==t2);
                MyLocalYX = ndLocalYX(MyBaseSpots,:);
                
                if t == t2
                    fprintf('Point cloud: ref round tile %d -> tile %d round %d base %d, %d matches, error %f\n', ...
                        t, t2, r, b,  o.GadInfo.nMatches(t,b,r), o.GadInfo.Error(t,b,r));
                    MyPointCorrectedYX = MyLocalYX*o.D(:,:,t,r,b)+o.TileCentre;
                    MyPointCorrectedYX = round(MyPointCorrectedYX);
                    ndSpotColors(MyBaseSpots,b) = IndexArrayNan(BaseImSm, MyPointCorrectedYX');
                else
                    [MyPointCorrectedYX, Error, nMatches] = o.different_tile_transform(o.GadRawLocalYX,o.RawLocalYX, ...
                        MyLocalYX,t,t2,r,b);
                    fprintf('Point cloud: ref round tile %d -> tile %d round %d base %d, %d matches, error %f\n', ...
                        t, t2, r, b,  nMatches, Error);
                    ndSpotColors(MyBaseSpots,b) = IndexArrayNan(BaseImSm, MyPointCorrectedYX');
                end
               
            end    
        end
        TifObj.close();       
    end
end
fprintf('\n');

o.pGadColor = ndSpotColors;
o.UseRounds = 1:o.nRounds;
%save(fullfile(o.OutputDirectory, 'oGad'), 'o', '-v7.3');


