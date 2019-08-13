function o = find_spots(o)          %ADDING t2 BIT BACK IN
% o = o.find_spots;
%
% finds spots in all tiles using the reference channel, removes
% duplicates in overlap regions and returns nSpots x 2 array o.SpotGlobalYX of
% coordinates in global frame
% 
% Looks up colors from apporpriate files and makes nSpots x nBP x nRounds
% array o.SpotColors
%
% o.Isolated is a nSpots x 1 binary array giving 1 for
% well-isolated spots
%
% NB spots that can't be read in all rounds are discarded
% 
% Kenneth D. Harris, 29/3/17
% GPL 3.0 https://www.gnu.org/licenses/gpl-3.0.en.html

%% variable naming conventions:
% spot subgroups:
% All: Any spot included in any tile (includes duplicates)
% nd: only spots whose anchor coordinate is in its home tile (no duplicates)
% Good: Spots for which could be read for all rounds

% coordinate frames or other info
% LocalYX: relative to home tile on the reference round
% LocalTile: number of home tile on the reference round
% GlobalYX: relative to the stitched image on the reference round
% RoundYX: relative to home tile after registration on each round
% RoundTile: number of home tile after registration on each round
% Isolated: binary number, saying if it is isolated
% SpotColors: the answer:

%% basic variables
rr = o.ReferenceRound;
Tiles = find(~o.EmptyTiles)';

[nY, nX] = size(o.EmptyTiles);
nTiles = nY*nX;

%% loop through all tiles, finding spots in anchor channel on ref round
RawLocalYXZ = cell(nTiles,1);  % cell array, giving spots in local coordinates
RawIsolated = cell(nTiles,1);
%SE = fspecial3('ellipsoid',[o.SmoothSize*o.Zpixelsize/o.XYpixelsize, o.SmoothSize*o.Zpixelsize/o.XYpixelsize, o.SmoothSize]);
SE = fspecial3('ellipsoid',o.SmoothSize);
for t=Tiles
    if mod(t,10)==0; fprintf('Detecting reference spots in tile %d\n', t); end
    [y,x] = ind2sub([nY nX], t);
    AnchorIm = o.load_3D(rr,y,x,o.AnchorChannel);
    if o.SmoothSize
        %NOT SURE HOW TO DO THIS SMOOTHING YET        
        AnchorImSm = imfilter(AnchorIm, SE);
    else
        AnchorImSm = AnchorIm;
    end
    [RawLocalYXZ{t}, RawIsolated{t}] = o.detect_spots(AnchorImSm);
end
    
%% now make array of global coordinates
AllIsolated = logical(vertcat(RawIsolated{:})); % I HATE MATLAB - for converting logical to doubles for no reason
nAll = length(AllIsolated);

AllGlobalYXZ = zeros(nAll,3);
AllLocalYXZ = zeros(nAll,3);
OriginalTile = zeros(nAll,1);

ind = 1;
for t=Tiles
    MySpots = RawLocalYXZ{t};
    nMySpots = size(MySpots, 1);
    AllGlobalYXZ(ind:ind+nMySpots-1,:) = bsxfun(@plus, MySpots, o.TileOrigin(t,:,rr));
    AllLocalYXZ(ind:ind+nMySpots-1,:) = MySpots;
    OriginalTile(ind:ind+nMySpots-1) = t;
    ind = ind+nMySpots;
end
if o.Graphics
    plotAnchorSpotsGlobal(AllGlobalYXZ,o.nZ,'All global coords including duplicates')
end

%% now remove duplicates by keeping only spots detected on their home tile

[AllLocalTile, ~] = which_tile(AllGlobalYXZ, o.TileOrigin(:,:,rr), o.TileSz);
NotDuplicate = (AllLocalTile==OriginalTile);
ndGlobalYXZ = AllGlobalYXZ(NotDuplicate,:);
ndLocalYXZ = AllLocalYXZ(NotDuplicate,:);
ndIsolated = AllIsolated(NotDuplicate,:);
ndLocalTile = AllLocalTile(NotDuplicate,:);

nnd = sum(NotDuplicate);

if o.Graphics
    plotAnchorSpotsGlobal(ndGlobalYXZ,o.nZ,'Global coords without duplicates')
end




%% get spot local coordinates in all colour channels and run PCR
AllBaseLocalYXZ = cell(nTiles,o.nBP, o.nRounds);



%Specify which rounds/colour channels to use (default is all)
if isempty(o.UseChannels)
    o.UseChannels = 1:o.nBP;
end

if isempty(o.UseRounds)
    o.UseRounds = 1:o.nRounds;
end

% loop through all tiles, finding PCR outputs
fprintf('\nLocating spots in each colour channel of tile   ');

%For scaling need to be centered about 0 hence subtract this
o.CentreCorrection = [1+(o.TileSz-1)/2,1+(o.TileSz-1)/2,1+(o.nZ-1)/2];

for t=1:nTiles
    if o.EmptyTiles(t); continue; end
    
    if t<10
        fprintf('\b%d',t);
    else
        fprintf('\b\b%d',t);
    end 
    
    [y, x] = ind2sub([nY nX], t);
    %Reload AnchorIm to do ImRegFFT3
    %AnchorIm = o.load_3D(rr,y,x,o.AnchorChannel);
    %if o.SmoothSize
    %    AnchorImSm = single(imfilter(AnchorIm, SE));
    %    clear AnchorIm
    %else
    %    AnchorImSm = AnchorIm;
    %end

    for r = o.UseRounds
        % find spots whose home tile on round r is t      
        % open file for this tile/round       
        % now read in images for each base
        
        FinalBaseIm = zeros(o.TileSz,o.TileSz,o.nZ);        %Max projected image of round to find initial shift
        for b = o.UseChannels              

            BaseIm = o.load_3D(r,y,x,o.FirstBaseChannel + b - 1);
            %SE = fspecial3('ellipsoid',[o.SmoothSize*o.Zpixelsize/o.XYpixelsize, o.SmoothSize*o.Zpixelsize/o.XYpixelsize, o.SmoothSize]);
            BaseIm = imfilter(BaseIm, SE);
            %o.MinThresh = max(mean(mean(BaseIm)));
            %o.DetectionThresh = 1.5*max(mean(mean(BaseIm)));
            % find spots for base b on tile t - we will use this for point
            % cloud registration only, we don't use these detections to
            % detect colors, we read the colors off the
            % pointcloud-corrected positions of the spots detected in the reference round home tiles  
            CenteredSpots = o.detect_spots(BaseIm) - o.CentreCorrection;
            %Scale so all in terms of XY pixel size. Import for PCR as find
            %nearest neighbours
            AllBaseLocalYXZ(t,b,r) = {CenteredSpots.*[1,1,o.Zpixelsize/o.XYpixelsize]};
            FinalBaseIm = max(FinalBaseIm,BaseIm); %INSTEAD OF THIS, JUST USE IMAGE WITH MOST SPOTS
        end
        [o.D0(t,:,r), cc] = o.ImRegFFt3D_FindSpots(FinalBaseIm,AnchorImSm, 0, o.RegMinSize); 
        %cc
    end      
end
fprintf('\n');

%PCR initial shifts
o.D0 = o.D0.*[1,1,o.Zpixelsize/o.XYpixelsize];       %Convert shift so units are XY pixels. 
%o = o.get_initial_shift(AllBaseLocalYXZ, RawLocalYXZ, nTiles, o.RegMinSize*1000);
%A=o.D0(:,:,1);
%A(:,3)=0;
%o.D0(:,:,1)=A;
%A=o.D0(:,:,2);
%A(:,3)=0;
%o.D0(:,:,2)=A;
%o.D0(7,:,1) = o.TileOrigin(1,:,rr)-o.TileOrigin(1,:,7);
%o.D0(14,:,3) = [10,45]; %This makes sense as there is a major artifact in
%this tile/round

o = o.PointCloudRegister_NoAnchor3DNoCA(AllBaseLocalYXZ, RawLocalYXZ, nTiles);


%% decide which tile to read each spot off in each round. 
% They are read of home tile if possible (always possible in ref round)
% in other rounds might have to be a NWSE neighbor - but never a diagonal
% neighbor
% ndRoundTile(s,r) stores appropriate tile for spot s on round r
% ndRoundYX(s,:,r) stores YX coord on this tile

%Compute approx new shifts in XY pixels, by taking the bottom row of the
%transform R. Then convert z shift back to units of z pixels for origin
XYPixelShifts = permute(squeeze(o.R(4,:,:,1:o.nRounds).*[1,1,o.XYpixelsize/o.Zpixelsize]),[2 1 3]);
o.TileOrigin(:,:,1:o.nRounds) =  o.TileOrigin(:,:,rr) - XYPixelShifts(:,:,1:o.nRounds);     

ndRoundTile = nan(nnd,o.nRounds);
ndRoundYXZ = nan(nnd,3,o.nRounds);

PossNeighbs = [-1 -nY 1 nY 0]; % NWSE then same tile - same will have priority by being last

for r=o.UseRounds
    fprintf('Finding appropriate tiles for round %d\n', r);
    
    for n = PossNeighbs
        % find origins of each tile's neighbor, NaN if not there
        NeighbTile = (1:nTiles)+n;
        NeighbOK = (NeighbTile>=1 & NeighbTile<=nTiles);
        NeighbOrigins = nan(nTiles,3);
        NeighbOrigins(NeighbOK,:) = round(o.TileOrigin(NeighbTile(NeighbOK),:,r));
        
        % now for each spot see if it is inside neighbor's tile area
        SpotsNeighbOrigin = NeighbOrigins(ndLocalTile,:);
        SpotsInNeighbTile = all(ndGlobalYXZ>=SpotsNeighbOrigin+1+o.ExpectedAberration...
            & ndGlobalYXZ<=SpotsNeighbOrigin+o.TileSz-o.ExpectedAberration, 2);
        
        % for those that were in set this to be its neighbor
        ndRoundTile(SpotsInNeighbTile,r) = NeighbTile(ndLocalTile(SpotsInNeighbTile));    
    end
    
    % compute YXZ coord
    HasTile = isfinite(ndRoundTile(:,r));
    ndRoundYXZ(HasTile,:,r) = ndGlobalYXZ(HasTile,:) - round(o.TileOrigin(ndRoundTile(HasTile,r),:,r));
    
end

%% loop through all tiles, finding spot colors

ndSpotColors = nan(nnd, o.nBP, o.nRounds);
ndPointCorrectedLocalYXZ = nan(nnd, 3, o.nRounds, o.nBP);

for t=1:nTiles
    if o.EmptyTiles(t); continue; end
    [y, x] = ind2sub([nY nX], t);
   
    for r=o.UseRounds      
        % find spots whose home tile on round r is t
        MySpots = (ndRoundTile(:,r)==t);
        if ~any(MySpots); continue; end
        
        % find the home tile for all current spots in the ref round
        RefRoundHomeTiles = ndLocalTile(ndRoundTile(:,r)==t);
        MyRefTiles = unique(RefRoundHomeTiles);
        fprintf('\nRef round home tiles for spots in t%d at (%2d, %2d), r%d: ', t, y, x, r);
        for i=MyRefTiles(:)'
            fprintf('t%d, %d spots; ', i, sum(RefRoundHomeTiles==i));
        end
        fprintf('\n');        
        
        
        % now read in images for each base
        for b=o.UseChannels              %No 0 as trying without using anchor

            
            BaseIm = o.load_3D(r,y,x,o.FirstBaseChannel + b - 1);
            
            if o.SmoothSize
                %BaseImSm = imfilter(double(BaseIm), fspecial('disk', o.SmoothSize));
                BaseImSm = imfilter(BaseIm, SE);
            else
                BaseImSm = BaseIm;
            end
            
            for t2 = MyRefTiles(:)'
                MyBaseSpots = (ndRoundTile(:,r)==t & ndLocalTile==t2);
                CenteredScaledMyLocalYXZ = [(ndLocalYXZ(MyBaseSpots,:) - o.CentreCorrection).*[1,1,o.Zpixelsize/o.XYpixelsize],...
                    ones(size(ndLocalYXZ(MyBaseSpots,:),1),1)];
                
                if t == t2
                    fprintf('Point cloud: ref round tile %d -> tile %d round %d base %d, %d/%d matches, error %f\n', ...
                        t, t2, r, b,  o.nMatches(t,b,r), size(RawLocalYXZ{t2},1), o.Error(t,b,r));
                    if o.nMatches(t,b,r)<o.MinPCMatches || isempty(o.nMatches(t,b,r))
                        continue;
                    end
                    CenteredMyPointCorrectedYXZ = (CenteredScaledMyLocalYXZ*o.R(:,:,t,r));
                    MyPointCorrectedYXZ = round(CenteredMyPointCorrectedYXZ.*[1,1,o.XYpixelsize/o.Zpixelsize] + o.CentreCorrection);
                    ndPointCorrectedLocalYXZ(MyBaseSpots,:,r,b) = MyPointCorrectedYXZ;
                    ndSpotColors(MyBaseSpots,b,r) = IndexArrayNan(BaseImSm, MyPointCorrectedYXZ');
                else
                    [MyPointCorrectedYXZ, error, nMatches] = o.different_tile_transform(AllBaseLocalYXZ,RawLocalYXZ, ...
                        CenteredScaledMyLocalYXZ,t,t2,r,b);
                    fprintf('Point cloud: ref round tile %d -> tile %d round %d base %d, %d/%d matches, error %f\n', ...
                        t, t2, r, b,  nMatches, size(RawLocalYXZ{t2},1), error);
                    if nMatches<o.MinPCMatches || isempty(nMatches)
                        continue;
                    end
                    ndPointCorrectedLocalYXZ(MyBaseSpots,:,r,b) = MyPointCorrectedYXZ;
                    ndSpotColors(MyBaseSpots,b,r) = IndexArrayNan(BaseImSm, MyPointCorrectedYXZ');
                end
               
            end    
        end      
    end
end
fprintf('\n');

%% now find those that were detected in all tiles
ndSpotColorsToUse = ndSpotColors(:,o.UseChannels,o.UseRounds);
Good = all(isfinite(ndSpotColorsToUse(:,:)),2);
GoodGlobalYXZ = ndGlobalYXZ(Good,:);
GoodSpotColors = ndSpotColors(Good,:,:);
GoodLocalTile = ndLocalTile(Good);
GoodIsolated = ndIsolated(Good);

save(fullfile(o.OutputDirectory, 'Intensities_NoAnchor.mat'), 'Good', 'ndGlobalYXZ', 'ndSpotColors', 'ndLocalTile');

%% plot those that were found and those that weren't
if o.Graphics
    plotSpotsResolved(o,ndGlobalYXZ,Good,Tiles,'Resolved Spots')
end
       

%% sanity check
plsz = 7;
if o.Graphics ==2
    GoodRoundYXZ = ndRoundYXZ(Good,:,:);
    GoodRoundTile = ndRoundTile(Good,:);
    GoodCorrectedYXZ = ndPointCorrectedLocalYXZ(Good,:,:,:);

    roi = o.FindSpotsRoi;
    PlotSpots = find(GoodGlobalYXZ(:,1)>roi(1) & GoodGlobalYXZ(:,1)<roi(2) & GoodGlobalYXZ(:,2)>roi(3) & GoodGlobalYXZ(:,2)<roi(4)...
        & round(GoodGlobalYXZ(:,3))>roi(5)& round(GoodGlobalYXZ(:,3))<roi(6));
    
    for s=(PlotSpots(:))' %PlotSpots(randperm(length(PlotSpots)))'
        figure(s); clf
        for r=o.UseRounds
            t=GoodRoundTile(s,r);
            [yTile,xTile] = ind2sub([nY nX], t);
            fprintf('Spot %d, round %d, tile %d: y=%d, x=%d, z=%d\n', s, r, t, GoodRoundYXZ(s,1,r), GoodRoundYXZ(s,2,r), GoodRoundYXZ(s,3,r));

            Ylegends = {o.bpLabels{:}};
            for b=o.UseChannels
                
                      
%                 if b==0                    
%                     y0 = GoodRoundYX(s,1,r);
%                     x0 = GoodRoundYX(s,2,r);
%                 else
%                     y0 = GoodCorrectedYX(s,1,r,b);
%                     x0 = GoodCorrectedYX(s,2,r,b);
%                 end
                y0 = GoodCorrectedYXZ(s,1,r,b);
                x0 = GoodCorrectedYXZ(s,2,r,b);
                z = round(GoodCorrectedYXZ(s,3,r,b));
                if ~isfinite(x0) || ~isfinite(y0)
                    continue;
                end
                y1 = max(1,y0 - plsz);
                y2 = min(o.TileSz,y0 + plsz);
                x1 = max(1,x0 - plsz);
                x2 = min(o.TileSz,x0 + plsz);
           
                BaseIm = imread(o.TileFiles{r,yTile,xTile,o.FirstBaseChannel + b - 1}, z, 'PixelRegion', {[y1 y2], [x1 x2]});
                if o.SmoothSize
                    %BaseImSm = imfilter(double(BaseIm), fspecial('disk', o.SmoothSize));
                    BaseImSm = imfilter(BaseIm, SE);
                else
                    BaseImSm = BaseIm;
                end

                subplot(o.nBP, o.nRounds, (b-1)*o.nRounds + r)
                imagesc([x1 x2], [y1 y2], BaseImSm); hold on
                axis([x0-plsz, x0+plsz, y0-plsz, y0+plsz]);
                plot(xlim, [y0 y0], 'w'); plot([x0 x0], ylim, 'w');
                caxis([0 o.DetectionThresh*2]);
                if r==1; ylabel(Ylegends{b+1}); end
                colorbar;
                
                title(sprintf('Round %d, Base %d, Tile %d', r, b, t));
                drawnow
            end
        end
        fprintf('\n');
        %figure(92); clf
        %imagesc(sq(GoodSpotColors(s,:,:)));
        %set(gca, 'ytick', 1:5); set(gca, 'yticklabel', {'Anchor', o.bpLabels{:}});
        %caxis([0 o.DetectionThresh*2]);
%         fprintf('local YX = (%f, %f) screen YX = (%f, %f) Called as %s, %s, quality %f\n', ...
%             GoodRoundYX(s,1), GoodRoundYX(s,2), GoodGlobalYX(s,1)/4, GoodGlobalYX(s,2)/4, ...
%             GoodCodes{s}, GoodGenes{s}, GoodMaxScore(s));
        %figure(1003); hold on
        %squarex = [-1 1 1 -1 -1]*plsz; squarey = [-1 -1 1 1 -1]*plsz;
        %h = plot(GoodGlobalYXZ(s,2)+squarex, GoodGlobalYXZ(s,1)+squarey, 'g');
        %pause;
        %delete(h);
    end
end



%%
o.SpotGlobalYXZ = GoodGlobalYXZ;
o.cSpotColors = GoodSpotColors;          
%o.cAnchorIntensities = squeeze(GoodSpotColors(:,1,:));
o.cSpotIsolated = GoodIsolated;
