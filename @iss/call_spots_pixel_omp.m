function o = call_spots_pixel_omp(o)
%% o = o.call_spots_pixel_omp
% 
% This is another probability method for gene calling.
% The difference here, is that an image is built up for each gene and then
% the spots are the local maxima on each gene image. This allows for
% multiple gene matches at each location i.e. overlapping spots. 
%
% o: iss object
% LookupTable: should be returned from call_spots_prob. It
% just gives the the probabilities that each spot score is explained by each
% gene. It saves calculating the probabilities explicitly each time.
%
% produces 
% pxSpotColors(Spot,b,r): intensity of Spot in channel b, round r
% pxSpotGlobalYX(Spot,:): global yx coordinate for each spot
% pxSpotCodeNo(Spot): gene index for each spot
% pxLogProbOverBackground(Spot): log of probability spot can be explained
% by gene relative to probability it can be explained by background.
% pxSpotScore(Spot): pxLogProbOverBackground of best gene match relative to
% second best gene match at that location.
% pxSpotScoreDev(Spot): standard deviation in spot scores across all genes
% at that location.
% pxSpotIntensity(Spot): intensity of the spot. Takes into account
% pxSpotCodeNo. Calculated by get_spot_intensity.
% 
%% Logging
if o.LogToFile
    diary(o.LogFile);
    cleanup = onCleanup(@()diary('off'));
end

% %% Specify version of python
% pcPythonExe = '/Users/joshduffield/Documents/UCL/ISS/VAE/BasicNeuralNetwork/venv/bin/python';
% PyDetails = pyenv;
% if PyDetails.Executable ~= pcPythonExe
%     [ver, exec, loaded]	= pyversion(pcPythonExe);
%     %pyenv('ExecutionMode','/Users/joshduffield/Documents/UCL/ISS/VAE/BasicNeuralNetwork/venv/bin/python');
%     % Add folder with desired script to python path
% end
% pyLibraryFolder = '/Users/joshduffield/Documents/UCL/ISS/VAE/BasicNeuralNetwork/';
% insert(py.sys.path, int64(0), pyLibraryFolder);
%% get scales for rounds and channels
nPixels = sum(o.HistCounts(:,1,1));
nCodes = length(o.CharCodes);
o.SHIFT = sum(o.HistValues'.*o.HistCounts)/nPixels;  %Want 0 mean
MeanOfSquare = sum(o.HistValues.^2'.*o.HistCounts)/nPixels;
o.SCALE = sqrt(MeanOfSquare - o.SHIFT.^2);
% Scale bled codes
o.ScaledBM = o.pBleedMatrix(:,:,1)./squeeze(o.SCALE);
UnbledCodes = reshape(o.UnbledCodes,[nCodes,o.nBP,o.nRounds]);
o.ScaledBledCodes = zeros(nCodes+o.nBP,o.nBP,o.nRounds);
for i=1:nCodes
    o.ScaledBledCodes(i,:,:) = o.ScaledBM * squeeze(UnbledCodes(i,:,:));
end
%Add background codes
for b=1:o.nBP
    UnbledBackground = zeros(o.nBP,o.nRounds);
    UnbledBackground(b,:) = 1;
    o.ScaledBledCodes(b+nCodes,:,:) = o.ScaledBM * UnbledBackground;
end
%All BledCodes need same norm for OMP to be fair.
o.ScaledBledCodes = o.ScaledBledCodes./vecnorm(o.ScaledBledCodes(:,:),2,2);
nBackground = size(o.ScaledBledCodes,1)-nCodes;
% ScaledBledCodes = mat2np(ScaledBledCodes);
% UnbledCodes = mat2np(UnbledCodes);
% ScaledBM = mat2np(ScaledBM);
%% Load in images, make images for each gene, find local maxima in images
rr = o.ReferenceRound;      %Round to base coordinate system on
[nY, nX] = size(o.EmptyTiles);
nTiles = nY*nX;
NonemptyTiles = find(~o.EmptyTiles)';

%Spots on achor round cover whole range of coordinates, same for each tile
AnchorLocalYX = zeros(o.TileSz^2,2);
AnchorLocalYX(:,1) = repelem(1:o.TileSz,1,o.TileSz);
AnchorLocalYX(:,2) = repmat(1:o.TileSz,1,o.TileSz);
nPixels = size(AnchorLocalYX,1);
AnchorLocalYX = [AnchorLocalYX,ones(nPixels,1)];

PeakSpotColors = cell(nCodes,1);
PeakLocalYX = cell(nCodes,1);
PeakCoefs = cell(nCodes,1);
PeakNeighbNonZeros = cell(nCodes,1);
OriginalTile = cell(nCodes,1);

for t=1:length(NonemptyTiles)  
    tile_no = NonemptyTiles(t);
    %Get pixel colors
    [GoodAnchorLocalYX,GoodSpotColors] = o.get_spot_colors(tile_no,AnchorLocalYX,nPixels);
    GoodSpotColors = (double(GoodSpotColors)-o.SHIFT)./o.SCALE;
    %OMP stops when reduction in residual drops below ResidualThresh.
    %Prctile bit gets 2nd largest intensity for each spot.
    ResidualThresh = o.ResidualThreshParam*prctile(abs(GoodSpotColors(:,:))',47.5*100/49.0)';
    ResidualThresh(ResidualThresh<o.ResidualThreshMin) = o.ResidualThreshMin;
    ResidualThresh(ResidualThresh>o.ResidualThreshMax) = o.ResidualThreshMax;
    nSpots = length(GoodSpotColors);
    coefs = zeros(nSpots,nCodes+nBackground);
    fprintf('Tile %d: Percentage of spot coefs found:       ',tile_no);
    for s=1:nSpots
        coefs(s,:) = omp_free_background(o.ScaledBledCodes(:,:)',GoodSpotColors(s,:)',...
            o.ompMaxGenes,ResidualThresh(s),nCodes+1:nCodes+nBackground)';
        if mod(s,round(nSpots/100))==0
            Percent = sprintf('%.6f', round(s*100/nSpots));
            fprintf('\b\b\b\b\b%s%%',Percent(1:4));
        end
    end 
    fprintf('\n');
    
    %Get local maxima log probabilities for each gene
    [tPeakLocalYX,tPeakSpotColors,tPeakCoefs,tPeakNeighbNonZeros,tOriginalTile] = ...
    o.detect_peak_genes_omp(coefs,GoodSpotColors,GoodAnchorLocalYX,tile_no);
    clearvars GoodSpotColors GoodAnchorLocalYX;
    
    %Keep data for all tiles together
    PeakSpotColors = cellfun( @(x,y) [x;y], PeakSpotColors, tPeakSpotColors, 'UniformOutput', false );
    PeakLocalYX = cellfun( @(x,y) [x;y], PeakLocalYX, tPeakLocalYX, 'UniformOutput', false );
    clearvars tPeakSpotColors tPeakLocalYX;
    PeakCoefs = cellfun( @(x,y) [x;y], PeakCoefs, tPeakCoefs, 'UniformOutput', false );
    PeakNeighbNonZeros = cellfun( @(x,y) [x;y], PeakNeighbNonZeros, tPeakNeighbNonZeros, 'UniformOutput', false );
    OriginalTile = cellfun( @(x,y) [x;y], OriginalTile, tOriginalTile, 'UniformOutput', false );
    clearvars tPeakCoefs tPeakNeighbourhoodNonZeros tOriginalTile;    
end


%% Deal with each file one by one
PeakGlobalYX = cell(nCodes,1);
for GeneNo = 1:nCodes
    PeakGlobalYX{GeneNo} = bsxfun(@plus,double(PeakLocalYX{GeneNo}),o.TileOrigin(OriginalTile{GeneNo},:,rr));
end 

% Remove duplicates by keeping only spots detected on their home tile
ndSpotColors = cell(nCodes,1);
ndGlobalYX = cell(nCodes,1);
ndCoefs = cell(nCodes,1);
ndNeighbNonZeros = cell(nCodes,1);
ndOriginalTile = cell(nCodes,1);

for GeneNo = 1:nCodes
    if o.Graphics==2
        figure(1001)
        plot(PeakGlobalYX{GeneNo}(:,2), PeakGlobalYX{GeneNo}(:,1), '.', 'markersize', 1);
        title('All global coords including duplicates');
        %set(gca, 'YDir', 'reverse');
    end

    [AllLocalTile, ~] = which_tile(PeakGlobalYX{GeneNo}, o.TileOrigin(:,:,rr), o.TileSz);
    NotDuplicate = (AllLocalTile==OriginalTile{GeneNo});
    ndSpotColors{GeneNo} = PeakSpotColors{GeneNo}(NotDuplicate,:,:);
    ndGlobalYX{GeneNo} = PeakGlobalYX{GeneNo}(NotDuplicate,:);
    ndCoefs{GeneNo} = PeakCoefs{GeneNo}(NotDuplicate,:);
    ndNeighbNonZeros{GeneNo} = PeakNeighbNonZeros{GeneNo}(NotDuplicate);
    ndOriginalTile{GeneNo} = OriginalTile{GeneNo}(NotDuplicate);

    if o.Graphics==2
        figure(1002); clf
        plot(ndGlobalYX{GeneNo}(:,2), ndGlobalYX{GeneNo}(:,1), '.', 'markersize', 1);
        title('Global coords without duplicates');
        drawnow;
        %set(gca, 'YDir', 'reverse');
    end
end

%Free up memory
clearvars PeakSpotColors PeakGlobalYX PeakCoefs PeakNeighbNonZeros ...
    OriginalTile AllLocalTile NotDuplicate

% Get final results
SpotCodeNo = cell(1,1);
SpotColors = cell(1,1);
GlobalYX = cell(1,1);
Coefs = cell(1,1);
NeighbNonZeros = cell(1,1);
AllOriginalTile = cell(1,1);
nGeneSpots = cell2mat(cellfun(@length,ndNeighbNonZeros,'uni',false));    
for GeneNo = 1:nCodes
    SpotCodeNo{1} = [SpotCodeNo{1};ones(nGeneSpots(GeneNo),1)*GeneNo];
    SpotColors{1} = [SpotColors{1};ndSpotColors{GeneNo}];
    GlobalYX{1} = [GlobalYX{1};ndGlobalYX{GeneNo}];
    Coefs{1} = [Coefs{1};ndCoefs{GeneNo}];
    NeighbNonZeros{1} = [NeighbNonZeros{1};ndNeighbNonZeros{GeneNo}];
    AllOriginalTile{1} = [AllOriginalTile{1};ndOriginalTile{GeneNo}];
end
clearvars ndSpotColors ndGlobalYX ndCoefs ndNeighbNonZeros ndOriginalTile


%% Add results to iss object

o.pxSpotColors = cell2mat(SpotColors);
o.pxSpotColors=o.pxSpotColors.*o.SCALE+o.SHIFT;
o.ompSpotCodeNo = cell2mat(SpotCodeNo);
o.pxSpotGlobalYX = cell2mat(GlobalYX);
o.pxLocalTile = cell2mat(AllOriginalTile);  
o.ompCoefs = cell2mat(Coefs);
o.ompNeighbNonZeros = cell2mat(NeighbNonZeros);
o.ompSpotIntensity = o.get_spot_intensity(o.ompSpotCodeNo,o.pxSpotColors);
end

