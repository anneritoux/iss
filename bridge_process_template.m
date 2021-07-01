%% Template of a pipeline that links together all the different processes 
%in the correct order.
%The first two sections below will be different for each experiment so
%should be checked before each run.

%% Parameters that should be checked before each run
%CHECK BEFORE EACH RUN
o = iss;
o.ReferenceRound = 8;           %Round that contains Dapi image (Starting from 1).
o.AnchorChannel =  ;            %Channel that has most spots in o.AnchorRound (Starting from 1).
o.DapiChannel = 1;              %Channel in o.AnchorRound that contains Dapi images (Starting from 1).
o.InitialShiftChannel = o.AnchorChannel;      %Channel to use to find initial shifts between rounds (Starting from 1).
o.ExtractScaleChannel = o.AnchorChannel;	  %Channel which is used to get scaling for imaging rounds (Starting from 1).
o.RawFileExtension = '.nd2';    %Format of raw data
o.LogToFile = 1;                %Set to 1 if you want to save command window to txt file, else set to 0.
%o.EmptyTiles = [1,5,7];        %Can specify tiles to work with if don't want to get tiles for all data.

%% File Names
%CHECK BEFORE EACH RUN
o.InputDirectory = '...\Experiment1\raw_data';     %Folder path of raw data

%FileBase{r} is the file name of the raw data of round r in o.InputDirectory
o.nRounds = 7;
o.nExtraRounds = 1;         %Treat Anchor channel as extra round
o.FileBase = cell(o.nRounds+o.nExtraRounds,1);
o.FileBase{1} = 'round0';
o.FileBase{2} = 'round1';
o.FileBase{3} = 'round2';
o.FileBase{4} = 'round3';
o.FileBase{5} = 'round4';
o.FileBase{6} = 'round5';
o.FileBase{7} = 'round6';
o.FileBase{8} = 'anchor';    %Make sure the last round is the anchor

o.TileDirectory = '...\Experiment1\tiles'; 
o.OutputDirectory = '...\Experiment1\output';  
%Codebook is a text file containing 2 columns - 1st is the gene name. 2nd is
%the code, length o.nRounds and containing numbers in the range from 0 to o.nBP-1.
o.CodeFile = '\\zserver\Data\ISS\codebook_73gene_6channels_2col.txt';

%% Logging
if o.LogToFile
    if isempty(o.LogFile)
        o.LogFile = fullfile(o.OutputDirectory,'Log.txt');
    end
end

%% extract and filter

%parameters
o.TileSz = 2048;
o.bpLabels = {'0', '1', '2', '3','4','5','6'}; %order of bases

%These specify the dimensions of the filter. R1 should be approximately the
%size of the spot in the respective direction and R2 should be double this.
o.ExtractR1 = 'auto';
o.ExtractR2 = 'auto';           %NOT SURE ABOUT THIS. IN 2D MUCH BETTER WITH 6
o.DapiR = 'auto';
%o.DapiR1YX = 20;
%o.DapiR1Z = 9;
%o.DapiR2YX = 40;
%o.DapiR2Z = 18;

o.ExtractScale = 'auto';
o.TilePixelValueShift = 15000;

%1:o.nZ will be the z planes read in for each imaging round.
o.nZ = 'auto';

%Max time (seconds) to wait for raw .nd2 files to be obtained
o.MaxWaitTime1 = 60;      %Less time for round 1 incase name is wrong
o.MaxWaitTime = 21600;  

%run code
% try
%     o = o.extract_and_filter;
% catch
%     o = o.extract_and_filter_NoGPU;
% end
save(fullfile(o.OutputDirectory, 'oExtract'), 'o', '-v7.3');
o = o.extract_and_filter_NoGPU;
save(fullfile(o.OutputDirectory, 'oExtract'), 'o', '-v7.3');

%% register

%parameters

%Anchor spots are detected in register2
o.DetectionRadius=[2,1];    %[YX_radius, Z_radius]
o.SmoothSize = 0;     %[Y_radius,X_radius,Z_radius]
o.IsolationRadius1 = 4;
o.IsolationRadius2 = [14,1];

o.DetectionThresh = 'auto';
o.MinThresh = 10;
o.minPeaks = 1;

%Registration Search parameters
o.RegMinScore = 'auto';     
o.RegStep = [5,5,2];
o.RegSearch.South.Y = -1900:o.RegStep(1):-1700;
o.RegSearch.South.X = -50:o.RegStep(2):50;
o.RegSearch.South.Z = -2:o.RegStep(3):2;

o.RegSearch.East.Y = -50:o.RegStep(1):50;
o.RegSearch.East.X = -1900:o.RegStep(2):-1700;
o.RegSearch.East.Z = -2:o.RegStep(3):2;

%run code
o = o.register2;
save(fullfile(o.OutputDirectory, 'oRegister'), 'o', '-v7.3');

%% find spots

%parameters
o.nBP = 7;

%If a channel or round is faulty, you can ignore it by selecting only the
%good ones in o.UseChannels and o.UseRounds.
o.UseChannels = 1:o.nBP;
o.UseRounds = 1:o.nRounds;

o.FirstBaseChannel = 1;

%Search parameters
o.FindSpotsMinScore = 'auto';     
o.FindSpotsStep = [5,5,2];
%FindSpotsSearch can either be a 1x1 struct or a o.nRounds x 1 cell of
%structs - have a different range for each round: 
%o.FindSpotsSearch = cell(o.nRounds,1);
o.FindSpotsSearch.Y = -100:o.FindSpotsStep(1):100;
o.FindSpotsSearch.X = -100:o.FindSpotsStep(2):100;
o.FindSpotsSearch.Z = -1:o.FindSpotsStep(3):1;

o.PcDist = 5; 
o.MinPCMatches = 1; 


%run code
o = o.find_spots2;
save(fullfile(o.OutputDirectory, 'oFind_spots'), 'o', '-v7.3');

%% call spots

%run code
o.CallSpotsCodeNorm = 'WholeCode';      %Or 'Round'
o = o.call_spots;
o = o.call_spots_prob;

save(fullfile(o.OutputDirectory, 'oCall_spots'), 'o', '-v7.3');

%% plot results

o.CombiQualThresh = 0.7;

Roi = round([1, max(o.SpotGlobalYXZ(:,2)), ...
    1, max(o.SpotGlobalYXZ(:,1)),...
    min(o.SpotGlobalYXZ(:,3)), max(o.SpotGlobalYXZ(:,3))]);
BackgroundImage = zeros(Roi(4),Roi(2),Roi(6)-Roi(5)+1,'uint16');
for z = Roi(5):Roi(6)
    BackgroundImage(:,:,z-Roi(5)+1) = imread(o.BigDapiFile, z,'PixelRegion', {Roi(3:4), Roi(1:2)});
end

o.PlotZThick=0;           %See all spots from o.PlotZThick+/- o.PlotZThick on single plane
o.plot(BackgroundImage);       %plot from call_spots

%Can interactively change plot by changing o or S.ZThick and then running
%iss_change_plot(o,'CallSpotsMethod')

%o.pScoreThresh = 10;
%iss_change_plot(o,'prob');             %plot from call_spots_prob

%Norm = 1;                      %Different normalisations for Norm=1,2,3
%iss_view_codes(o,93454,Norm);     %See result from call_spots
%iss_view_prob(o,93454,Norm);     %See result from call_spots_prob

%iss_change_plot(o,'DotProduct');       %change back to call_spots method
