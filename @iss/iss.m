%% iss: code for processing of in situ sequencing
% Kenneth D. Harris and Xiaoyan Qian
% GPL 3.0 https://www.gnu.org/licenses/gpl-3.0.en.html
%
% to use:
% o = iss; % create structure, default parameters
% % change any parameters you want, and set o.FileN
% o = o.extract_and_filter; % create top-hat filtered tiffs for each tile
% o = o.find_spots; % find spot positions in global coordinates
% o = o.call_spots; % allocate each spot to a gene
% o = o.call_cells; % identify cells
%
% this current file iss.m just contains default parameter values

classdef iss
    properties
        %% interactive graphics mode. 1 means some, 2 means a lot.
        Graphics = 1;
        
        % region to show during cell calling
        CellCallShowCenter = [385 1100];
        CellCallShowRad = 200;
        
        % example cell to diagnose calling, plus class pair to compare (two
        % strings)
        ExampleCellCenter = [1670 250];
        CellCallDiagnosisPair = [];
        
        %% parameters: file locations
        % where the input czi files are kept
        InputDirectory;
        
        % where the top-hat image files are kept (one for each tile and round, all channels)
        TileDirectory;
        
        % where .mat file outputs and intermediates are kept
        OutputDirectory;
        
        % code file for spot calling
        CodeFile;
        
        %% parameters: stuff that might vary between experiments
        
        % one more round after the main rounds for Anchor
        nExtraRounds = 1; 
        
        % For single-gene rounds, a n by 3 structure array of {GeneName, round, channel, threshold}
        ExtraCodes = {'Sst', 6, 3, 500; 'Npy', 6, 5, 400};
        
        % BasePair labels
        bpLabels = {'A', 'C', 'G', 'T'};
        
        %% parameters: extract and filter
        
        % To extract spots, a filter is used with inner radius ExtractR1 within 
        % which it is positive and outer radius ExtractR2 so the annulus between 
        % R1 and R2 is negative. Overall sums up to 0. R1 should be the
        % approximate radius of spot. R2 should be such that between 
        % R1 and R2 there is a dark region - good guess is double R1.
        % RYX has units of YX pixels, RZ has units of Z pixels.
        ExtractR1YX = 5;
        ExtractR1Z = 5;
        ExtractR2YX = 10;
        ExtractR2Z = 10;
        
        % DapiR is the same as above but is applied to the Dapi image so
        % should be the size of the Dapi spots. HAVN'T CHECKED THE FILTER
        % FOR THE DAPI IMAGE YET
        DapiR1YX = 20;
        DapiR1Z = 9;
        DapiR2YX = 40;
        DapiR2Z = 18;      
        
        % Each filtered image is multiplied by ExtractScale. This is
        % because the image is saved as uint16 so to gain information from
        % the decimal points, should multiply image so max pixel number is
        % in the 10,000s (less than 65,536). If 'auto', it sets to
        % 10,000/max(Tile 1 round 1 colour channel 1).
        ExtractScale = 'auto';
        %ExtractScale = 5*10^7;
        
        % Dapi has different scaling as different filter used.
        DapiScale = 'auto';
        
        %% parameters: registration and alignment
        
        % correlation threshold for image alignment. Can be low since 
        % chance correls are miniscule. Second value is for shifts of 0 
        % which can occur spuriously due to hot pixels so need a stronger
        % threshold
        RegCorrThresh = [.3 .6]; 
        
        % smoothing before registration (size of disk filter)
        RegSmooth = 1;
                
        % minimum size overlap for tile matching (pixels - so side squared)
        RegMinSize = 200^2; 
        
        % distance scale for point cloud registration (pixels)
        PcDist = 3; 
        
        % PcIter is the max number of iterations done by the PCR algorithm.
        PcIter = 50;
        
        % number of point cloud matches needed to count an overlap
        MinPCMatches = 50; 
        
        % ToPlot: [r,c,t], plot of round r, colour channel c, tile t
        % is shown for debugging
        ToPlot;
        
        % InitialShiftChannel is the chosen channel to use to register all
        % rounds to the anchor. Choose the one with clear spots, recommend 5,6 or 7.
        InitialShiftChannel = 5;
        
        % MaxRoundShift(1) and MaxRoundShift(2) are the maximum absolute
        % shifts in the XY and Z directions respectively
        % of a tile in any round relative to the anchor round.
        MaxRoundShift = [250,10];
        
        %Below are parameters for ImRegFFt3D_Register. They shouldn't need
        %changing.
        
        %Registration is forced to find an overlap smaller than MaxOverlapFract*TileSz
        %between neighbouring Tiles
        MaxOverlapFract = 0.2;
        
        %MaxRegShift is the maximum shift of a tile in the non overlapping
        %direction i.e. South/North if looking at its east neighbour.
        MaxRegShift = 20;
        
        %ToTest is the fraction of shifts to look at before returning a NaN shift
        ToTest = 0.25;
        
        %The Score used to find the best shift in get_initial_shift2 is 
        %sum(exp(-Dist.^2/(2*o.ShiftScoreThresh^2)))
        ShiftScoreThresh = 2;
        
        %RegSearch.Direction.Y,RegSearch.Direction.X,RegSearch.Direction.Z 
        %(.Z in Z pixel units) are the ranges values of
        %shifts to check during the registration to the neighbour in the corresponding
        %Direction (South or East)
        RegSearch;
        
        %RegStep specifies the step size to use in the Y,X,Z direction of
        %RegSearch (Z in Z pixel units)
        RegStep = [5,5,2];
        
        %if the score mentioned above is below RegMinScore, the search
        %range will be enlarged.
        RegMinScore = 30;
        
        %RegWidenSearch specifies how much to widen the search range in the
        %Y,X,Z directions respectively if the score is below MinRegScore.
        %The search will go from min(RegSearch):RegStep:max(RegSearch) to
        %min(RegSearch)-RegWidenSearch:RegStep:max(RegSearch)+RegWidenSearch
        RegWidenSearch = [50,50,9];
        
        %After the initial broad search, get_initial_shift will do a
        %refined search with a range of
        %BestShift-RegRefinedSearch:BestShift+RegRefinedSearch inclusive.
        RegRefinedSearch = [12,12,1];
        
        %RegisterRefinedStep is the step size to use in this refined
        %search
        RegisterRefinedStep = [3,3,1];
        
        %if the score is below RegAbsoluteMinScore, the shift found will be
        %set to the average of all the other shifts
        RegAbsoluteMinScore = 4;
        
        %OutlierThresh is the number of scaled MAD away from the median
        %that constitutes an outlier when considering the shifts in the
        %registration or find_steps methods. These outliers are set to the
        %average of all acceptable shifts if the score is low enough
        OutlierThresh = 5;
        
        
        %% parameters: spot detection
        
        % smooth images before reading out fluorescence with a disk of this radius:
        SmoothSize = 1;
        
        % to detect spot, pixel needs to be above dilation with this radius
        DetectionRadius = 1;
        
        % and pixel raw fluorescence needs to be above this value:
        DetectionThresh = 300;
        
        % if o.DetectionThresh is 'auto' you get m* the pth percentile of
        % each frame
        AutoThreshPercentile = 99.95;
        AutoThreshMultiplier = .25;
        
        %If find less than o.minPeaks spots then DectionThresh is lowered by
        %o.ThreshParam        
        minPeaks = 1000;
        ThreshParam = 5;
        
        % find isolated spots by annular filtering with these radii
        IsolationRadius1 = 2;
        IsolationRadius2 = 7;
        
        % annular filtered value needs to be less than this:
        IsolationThresh = 60;
        
        % for visualization during spot detection
        FindSpotsRoi = [1742 1755 213 227];
        
        % Each spot will be allocated to home tile if possible - but not if
        % it is this close to the edge, because of aberrations
        ExpectedAberration = 3;
        
        %Centre Correction is the value to be subtracted from spot
        %positions so they are centered about 0
        CentreCorrection;
        
        %MinThresh is the smallest value DectionThresh can get to       
        MinThresh;
        
        %FindSpotsSearch.Y,FindSpotsSearch.X,FindSpotsSearch.Z 
        %(.Z in Z pixel units) are the ranges values of
        %shifts to check when looking for the initial shifts between rounds
        %for each tile. Can also set to cell(o.nRounds,1) and give a
        %different search range for each round.
        FindSpotsSearch;
        
        %FindSpotsStep specifies the step size to use in the Y,X,Z direction of
        %FindSpotsSearch (Z in Z pixel units)
        FindSpotsStep = [5,5,2];
        
        %if the score mentioned above is below FindSpotsMinScore, the search
        %range will be enlarged.
        FindSpotsMinScore = 90;
        
        %RegWidenSearch specifies how much to widen the search range in the
        %Y,X,Z directions respectively if the score is below MinRegScore.
        %The search will go from min(FindSpotsSearch):RegStep:max(FindSpotsSearch) to
        %min(FindSpotsSearch)-FindSpotsWidenSearch:FindSpotsStep:max(FindSpotsSearch)+FindSpotsWidenSearch
        FindSpotsWidenSearch = [50,50,9];
        
        %After the initial broad search, get_initial_shift will do a
        %refined search with a range of
        %BestShift-FindSpotsRefinedSearch:BestShift+FindSpotsRefinedSearch
        %inclusive.
        FindSpotsRefinedSearch = [12,12,1];
        
        %FindSpotsRefinedStep is the step size to use in this refined
        %search
        FindSpotsRefinedStep = [3,3,1];
        
        %if the score is below FindSpotsAbsoluteMinScore, the shift found will be
        %set to the average of all the other shifts in a particular round
        FindSpotsAbsoluteMinScore = 10;
        
        
        %% parameters: spot calling
        % normalizes spot fluorescence so this percentile = 1
        SpotNormPrctile = 98;
        
        % if BleedMatrixType == 'Separate', then a bleed matrix will be
        % computed for each round. If BleedMatrixType == 'Single', a single
        % bleed matrix will be computed, combining spot colours from all
        % rounds.
        BleedMatrixType = 'Single';
        
        % score and intensity thresholds to plot a spot (combi codes)
        CombiQualThresh = .8;         
        CombiIntensityThresh = .1;
        CombiAnchorsReq = 4; % need at least this many anchor chans above threshold
        
        nRedundantRounds = 0;
        RedundantPseudobases = {'AC', 'GT'}; 
%         RedundantCodes = {...
%                 '.11..', '.12..', '.22..', '.21..';...
%                 '..11.', '..12.', '..22.', '..21.';...
%                 '1..1.', '2..1.', '2..2.', '1..2.';...
%                 '11...', '12...', '22...', '21...';...
%             };

        RedundantCodes = {...
                '.11..', '.12..', '.21..', '.22..';...
                '..11.', '..12.', '..21.', '..22.';...
                '1..1.', '2..1.', '1..2.', '2..2.';...
                '11...', '12...', '21...', '22...';...
            };

        
        %% parameters: cell segmentation
        
        % percentile threshold for dapi image (after imadjust)
        DapiThresh = 90;
        
        % how close local maxima can be together (in pixels)
        DapiMinSep = 7;
        
        % how big cell radius needs to be to be detected (in pixels)
        DapiMinSize = 5;
        
        % how much margin around each cell to assign to it (in pixels)
        DapiMargin = 10;
        
        % minimum area in pixels for a cell not to be dropped
        MinCellArea = 200;
        
        %% parameters: cell calling
        
        % CellCallRegionXY(:,1:2): XY coords of polygod outlining region to
        % look for cells in
        CellCallRegionYX;
       
        % how many neighboring cells to consider for each spot
        nNeighbors = 3; 
        
        % gamma dist prior for each spot
        rSpot = 2; 
        
        % gamma dist shape for scaling whole gene
        rGene = 20; 
        
        % probability of misread per pixel
        MisreadDensity = 1e-5; 
        
        % prior on how much to scale scRNAseq down by to match iss
        Inefficiency = .2; 
        
        % how much to add to expression so as not to get infities
        SpotReg = .1; 
        
        % additional log likelihood for a gene inside the segmented region
        InsideCellBonus = 2; 
        
        % converges when no probabilities have changed more than this
        CellCallTolerance = .02; 
        
        % converges when no probabilities have changed more than this
        CellCallMaxIter = 100; 
        
        % for pie-plots: don't both showing probs below this
        MinPieProb = .1;
        
        % size of pies
        PieSize = 25;
        
        % for plotting: any cell classes starting with these names will be
        % collapsed and given the specified color
        ClassCollapse = {{'Astro', 'Endo', 'Oligo', 'Eryth', 'Vsmc'}, 'NonNeuron', [.5 .5 .5] ; ...
                         {'PC.CA1'}, 'PC CA1', [1 .8 .8] ; ...
                         {'PC.CA2', 'PC.CA3'}, 'PC subtypes', [.8 1 .8] ; ...
%                          {'PC.CA3'}, 'PC CA3', [.8 1 .8] ; ...
                         {'Zero'}, 'Zero', [0 0 0]};


        
        %% parameters: stuff that should be the same between experiments
        
        % for dapi images: scale is downsampling factor for final .fig file
        DapiChannel = 1;  
        
        % which channel of each file is anchor images
        AnchorChannel = 2;
        
        % which channel has first round of sequencing images
        FirstBaseChannel = 3;
        
        % which sequencing round to align all others to
        ReferenceRound = 2;
        
        % how many combinatorial sequencing rounds excluding anchor round
        nRounds = 5;
        
        % Number of possible basepairs (always 4 for life as we know it but
        % hardcoding is bad programming style!)
        nBP = 4;
        
        % tile size (assumed square tiles)
        TileSz = 2048;
        
        % expected step size (in coordinates returned by bioformats)
        MicroscopeStepSize = 2048;
        
        RawFileExtension = '.czi';
        
        % when decoding spots, will use all colour channels in
        % UseChannels (Array of numbers in range 1 to o.nBP)
        UseChannels;
        
        % when decoding spots, will use all rounds in
        % UseRounds (Array of numbers in range 1 to o.nRounds)
        UseRounds;
        
        %nZ is the number of Z planes per tile
        nZ;
        
        %pixel sizes are different in the 2 directions, need to consider
        %when filtering.
        XYpixelsize;
        Zpixelsize;
        

        
        %% variables: filenames
        
        % TileFiles{r, y, x, z} is full pathname of tile (y,x,z) on round r
        % (mutlicolor tiff file). Empty if no file there
        TileFiles; 
        
        % EmptyTiles{y,x} is zero if that tile location wasn't scanned
        EmptyTiles;
        
        % FileBase{r} is local filename of round r .czi file; other
        % extensions added in output files
        FileBase;
        
        % BigDapiFile is full path of stitched DAPI image (tiff format)
        BigDapiFile;
        
        % CellMapFile is full path of cell map (.mat format)
        CellMapFile;
        
        %% variables: registration
        
        % TilePosYXZ(t,:): grid position of tile t in integers. t is what you
        % find in the file name - so it only counts non-empty tiles. This
        % is confusing and ought to be refactored. Manana. YX and Colour
        % channel position of Tiff
        TilePosYXC;
        
        % TileOriginYX(t,:,r): YX coordinate of origin of tile t on round
        % r. To get the global coordinate of a point, add this to the local
        % coordinate within the tile (counting from 1)
        TileOrigin;
        
        % TileInitialPosXY(t,:): coordinate of tile t in integers. This is
        % the derived coordinates, before using the fact we know the answer
        % to get TilePosYXC.
        TileInitialPosYX;
        
        %RawLocalYXZ{t} stores the YXZ coordinates of spots found in the
        %anchor round of tile t
        RawLocalYXZ;
        
        %RawIsolated{t} labels each spot in the anchor round as isolated or not
        RawIsolated;
        
        %RegInfo saves debugging information for the registration section.
        %h = horizontal, v = vertical
        %Score(i) is the score for the Shifts(i) found between tiles given by Pairs(i)
        %ChangedSearch is number of times the search range had to be changed.
        RegInfo;    
        
        % D0(t,2,r) stores the initial shift to use as a starting point for
        % the PCR on round r tile t.
        D0;
        
        %FindSpotsInfo saves debugging information for finding the initial
        %shifts. 
        %Score(t,r) is the score for the shift D0(t,r) found for tile t 
        %between the anchor round and round r.
        %ChangedSearch is number of times the search range had to be changed.
        %Outlier is a shift found that was wrong and subsequently changed
        %to the average of all other shifts.
        FindSpotsInfo;
        
        % A(2,2,c): stores the scaling correction for chromatic aberration
        % found by point cloud registration for color channel c
        A;
        
        % D(t,2,r): stores the final shift found by point cloud registration
        % on round r tile t.
        D;
        
        %R(4,3,t,r) stores the final transformation found by point cloud
        %registration on rounf r tile t.
        R;
        
        % nMatches(t,c,r): stores number of matches found by point cloud
        % registration for tile t, color channel c, round r
        nMatches;
        
        % error(t,c,r): stores error found by point cloud registration
        % for tile t, color channel c, round r
        Error;

        %% variables: spot calling outputs
       
        % cSpotColors(Spot, Base, Round) contains spot color on each base
        % and round. only for combinatorial splots
        cSpotColors;             
        
        % cSpotIsolated(Spot) is a binary array saying if the spot is well isolated
        % again for combinatorial spots only
        cSpotIsolated; 
        
        % SpotCombi: binary array which is one for combinatorial spots,
        % zero for extra spots
        SpotCombi;
        
        % SpotGlobalYX(Spot,1:2) contains y,x coordinates of every spot in
        % global coordiate system. Both combinatorial and extra spots
        SpotGlobalYXZ;
        
        % SpotCodeNo(Spot): code number for each spot (uint16). Both combinatorial and extra spots
        SpotCodeNo;        
     
        % SpotScore(Spot): score saying how well the code fits (0...1).
        % 1 for extras
        SpotScore;
        
        % SpotIntensity(Spot): RMS intensity of the spot. Zero for
        % extras!
        SpotIntensity;
        
        % CharCodes{Code}: text representation of each code in the
        % codebook. For extra spots, this just says "EXTRA"
        CharCodes;
        
        % GeneNames{Code}: gene name for each code. Note that some genes
        % can have multiple entries, because they have multiple codes or
        % because they have extra rounds
        GeneNames;
        
        % UnbledCodes(nCodes, nBP*nRounds): binary code vectors
        UnbledCodes;
        
        % BledCodes(nCodes, nBP*nRounds): code vectors after modeling
        % crosstalk
        BledCodes;
        
        % Normalised Spot Scores
        NormBledCodes;
        cNormSpotColors;
        
        % BleedMatrix used to estimate BledCodes
        BleedMatrix;
		
        %% variables: cell calling outputs
        % pCellClass(cell, class); % prob each cell goes to each class: last class is zero expression
        pCellClass;
        
        % ClassNames(class): name of each class, taken from gSet file
        ClassNames;
        
        % pSpotCell(spot, cell): sparse array containing prob of each spot
        % to belong to each cell. Last cell is background
        pSpotCell;
        
        % position of each cell centroid
        CellYX;
        
        
    end
    
end

 





