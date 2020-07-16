%% Get gene calling data for Gad ground truth for different o.ScoreScale.
% Then put data into csv so can do thresholding later
ScoreScales = [0:0.1:1,2,10000];
nSpots = size(o.cSpotColors,1);
%ToUse = GadGroundTruthLogical(o,'Prob') & o.GadPeakSpots==1 & o.GadPeakColor>250;
ToUse = o.GadPeakSpots<2;
nUse = sum(ToUse);
nData = nUse*length(ScoreScales);
nColumns = 9;

nCodes = length(o.CharCodes);
gChannelIndex = repmat(1:o.nBP,1,o.nRounds);
gRoundIndex = repelem(1:o.nRounds,1,o.nBP);
ChannelIndex = repmat(gChannelIndex,1,nCodes);
RoundIndex = repmat(gRoundIndex,1,nCodes);
GeneIndex = repelem(1:nCodes,1,o.nRounds*o.nBP);

BackgroundLogProb = zeros(nSpots,49);
for s=1:nSpots
    BackgroundIndices = sub2ind(size(o.BackgroundProb),o.ZeroIndex-1+o.cSpotColors(s,:),gChannelIndex,gRoundIndex);
    BackgroundLogProb(s,:) = log(o.BackgroundProb(BackgroundIndices));
end


LogProbMultiplierDefault = zeros(o.nRounds*o.nBP,nCodes);
for g=1:nCodes
    GeneChannels = str2double(regexp(cell2mat(o.CharCodes(g)),'\d','match'))+1;  
    gUnbled = o.UnbledCodes(g,:);
    LogProbMultiplierDefault(:,g) = o.UnbledCodes(g,:);
end

LogProbMultiplier = zeros(o.nRounds*o.nBP,nCodes);
ScoreScaleCompAllScaleValuesData = zeros(nData,nColumns);
idx = 0;
for ScoreScale=ScoreScales
    o.ScoreScale = ScoreScale;
    NormFactor = 49.0/(7+o.ScoreScale*42);  %How many squares that contribute.
                                            %Normalise by this to allow valid comparison
    LogProbMultiplier(:) = LogProbMultiplierDefault;
    LogProbMultiplier(LogProbMultiplierDefault==0) = o.ScoreScale;
    LogProbMultiplier = LogProbMultiplier*NormFactor;   
    
    LogProbOverBackground = zeros(nSpots,nCodes);
    for s=1:nSpots
        SpotIndex = repmat(o.ZeroIndex-1+o.cSpotColors(s,:),1,nCodes); %-1 due to matlab indexing I think
        Indices = sub2ind(size(LookupTable),SpotIndex,GeneIndex,ChannelIndex,RoundIndex);
        LogProbMatrix = reshape(LookupTable(Indices),[o.nRounds*o.nBP,nCodes]);
        LogProbOverBackground(s,:) = sum((LogProbMatrix-BackgroundLogProb(s,:)').*LogProbMultiplier);
    end
    GadLogProbOverBackground = LogProbOverBackground(ToUse,22);
    [LogProbOverBackground,SpotCodeNo] = sort(LogProbOverBackground,2,'descend');
    BestLogProbOverBackground = LogProbOverBackground(ToUse,1);
    SecondBestLogProbOverBackground = LogProbOverBackground(ToUse,2);
    BestSpotCodeNo = SpotCodeNo(ToUse,1);
    SpotScore = LogProbOverBackground(ToUse,1)-LogProbOverBackground(ToUse,2);
    GadSpotScore = GadLogProbOverBackground - LogProbOverBackground(ToUse,2);
    SpotScoreDev = std(LogProbOverBackground(ToUse,:),[],2);
    AllSpotIntensity = o.get_spot_intensity(SpotCodeNo(:,1),o.cSpotColors);
    SpotIntensity = AllSpotIntensity(ToUse);
    [GadRowIdx,SpotNumbers] = find(SpotCodeNo'==22);
    NestToGadRowIdx = min(GadRowIdx+1,nCodes);
    NextToGadIdx = sub2ind(size(SpotCodeNo),SpotNumbers,NestToGadRowIdx);
    NextToGadLogProbOverBackground = LogProbOverBackground(NextToGadIdx);
    NextToGadLogProbOverBackground = NextToGadLogProbOverBackground(ToUse);
    GadOverNext = GadLogProbOverBackground-NextToGadLogProbOverBackground;
    GadRank = GadRowIdx(ToUse);
    
    ScoreScaleCompAllScaleValuesData(idx*nUse+1:(idx+1)*nUse,1) = o.ScoreScale;
    ScoreScaleCompAllScaleValuesData(idx*nUse+1:(idx+1)*nUse,2) = BestSpotCodeNo;
    ScoreScaleCompAllScaleValuesData(idx*nUse+1:(idx+1)*nUse,3) = BestLogProbOverBackground; 
    ScoreScaleCompAllScaleValuesData(idx*nUse+1:(idx+1)*nUse,4) = SecondBestLogProbOverBackground;
    ScoreScaleCompAllScaleValuesData(idx*nUse+1:(idx+1)*nUse,5) = GadLogProbOverBackground;
    ScoreScaleCompAllScaleValuesData(idx*nUse+1:(idx+1)*nUse,6) = NextToGadLogProbOverBackground;
    ScoreScaleCompAllScaleValuesData(idx*nUse+1:(idx+1)*nUse,7) = SpotIntensity;
    ScoreScaleCompAllScaleValuesData(idx*nUse+1:(idx+1)*nUse,8) = SpotScoreDev;
    ScoreScaleCompAllScaleValuesData(idx*nUse+1:(idx+1)*nUse,9) = GadRank;
    idx=idx+1;
end
csvwrite(fullfile(o.OutputDirectory, 'ScoreScaleCompAllScaleValuesData.csv'),...
    ScoreScaleCompAllScaleValuesData);


%% Read in csv
% ScoreScaleCompAllSpotData columns are:
% SpotNo | cSpotIsolated | GadPeakColor | GadSpotIntensity | GadPeakSpots
AllSpotData = dlmread(fullfile(o.OutputDirectory, 'ScoreScaleCompAllSpotData.csv'));
SpotNo = AllSpotData(:,1); cSpotIsolated = AllSpotData(:,2);
GadPeakColor = AllSpotData(:,3); GadSpotIntensity = AllSpotData(:,4);
GadPeakSpots = AllSpotData(:,5);
GadGroundTruth = GadGroundTruthLogical(o,'Prob');
ToUse = o.GadPeakSpots<2;
GadGroundTruth = GadGroundTruth(ToUse);
GadPeakColorThresh = 250;
% Set to look at to see if correctly getting Gad
FalseNegativesSet = GadGroundTruth &...
    GadPeakColor > GadPeakColorThresh & GadPeakSpots == 1;
nFalseNegativesSet = sum(FalseNegativesSet);
% Set to look at to see if we aren't just saying everything is Gad
FalsePositivesSet = ~GadGroundTruth &...
    GadPeakColor < GadPeakColorThresh & GadPeakSpots == 0;
nFalsePositivesSet = sum(FalsePositivesSet);
clear AllSpotData;

% ScoreScaleCompAllScaleValuesData columns are:
% ScoreScale | BestSpotCodeNo | BestLogProbOverBackground |
% SecondBestLogProbOverBackground | GadLogProbOverBackground | 
% NextToGadLogProbOverBackground | SpotIntensity | SpotScoreDev | GadRank
AllScaleValuesData = dlmread(fullfile(o.OutputDirectory, 'ScoreScaleCompAllScaleValuesData.csv'));
load('/Users/joshduffield/Documents/UCL/ISS/B9S4_Slice001-002/Slice001output/ByPixelAnchor/NewTransform-PCR6/Gamma/GadGroundTruthLogicals.mat');
%% Analysis
% ROC curve
% With default thresholds (tuned to Scale=1 but should apply to other
% scales???), see how true positive and false negatives vary with scale.
o.pScoreThresh = 80;
o.pScoreThresh2 = 0;
o.pLogProbThresh = 0;
o.pIntensityThresh = 100;
o.pIntensityThresh2 = 0;
o.pLogProbThresh2 = 0;
TruePositveRate = zeros(size(ScoreScales))';
FalsePositveRate = zeros(size(ScoreScales))';
idx=1;
for Scale = ScoreScales
    ToUse = AllScaleValuesData(:,1) == round(Scale,1);
    AllScaleValuesData2 = AllScaleValuesData(ToUse,:);
    ScoreScale = AllScaleValuesData2(:,1); BestSpotCodeNo = AllScaleValuesData2(:,2);
    BestLogProbOverBackground = AllScaleValuesData2(:,3);
    SecondBestLogProbOverBackground = AllScaleValuesData2(:,4);
    GadLogProbOverBackground = AllScaleValuesData2(:,5);
    NextToGadLogProbOverBackground = AllScaleValuesData2(:,6);
    SpotIntensity = AllScaleValuesData2(:,7);
    SpotScoreDev = AllScaleValuesData2(:,8); GadRank = AllScaleValuesData2(:,9);
    SpotScore = BestLogProbOverBackground - SecondBestLogProbOverBackground;
    GadSpotScore = GadLogProbOverBackground - SecondBestLogProbOverBackground;
    GadOverNext = GadLogProbOverBackground - NextToGadLogProbOverBackground;
    %QualOK = GadRank==1;
    %NoOverlap:
    QualOK = (GadSpotScore>o.pScoreThresh & GadSpotIntensity>o.pIntensityThresh2 | ...
        GadSpotScore>o.pScoreThresh2 & GadSpotScore+GadLogProbOverBackground>o.pLogProbThresh2 &...
        GadSpotIntensity>o.pIntensityThresh2);
        %| GadSpotIntensity>o.pIntensityThresh & GadLogProbOverBackground>o.pLogProbThresh & ...
        %GadSpotScore>o.pScoreThresh2);    
    %Overlapping:
%     QualOK = (GadSpotScore>o.pScoreThresh & GadSpotIntensity>o.pIntensityThresh2 | ...
%         GadSpotIntensity>o.pIntensityThresh & GadLogProbOverBackground>o.pLogProbThresh & GadSpotScore>o.pScoreThresh2 |...
%         GadRank==2 & GadOverNext>o.pScoreThresh & GadSpotIntensity>o.pIntensityThresh2 | ...
%         GadOverNext > BestLogProbOverBackground-GadLogProbOverBackground &...
%         GadOverNext > o.pScoreThresh & GadSpotIntensity>o.pIntensityThresh2);
    TruePositveRate(idx) = sum(QualOK&FalseNegativesSet);
    FalsePositveRate(idx) = sum(QualOK&FalsePositivesSet);
    idx=idx+1;
    
end
figure;
plot(FalsePositveRate(1:end-1),TruePositveRate(1:end-1));
title('ROC curve with ScoreScale as variable parameter');
ylabel(['True Positive Rate (/',num2str(nFalseNegativesSet),')']);
xlabel(['False Positive Rate (/',num2str(nFalsePositivesSet),')']);

%Image of how FalseNegative+FalsePositive varies with pScoreThresh and
%pIntensityThresh2 as these are the thresholding parameters of the most
%consequence. Want this to be minimised
Scale=0;
ToUse = AllScaleValuesData(:,1) == round(Scale,1);
AllScaleValuesData2 = AllScaleValuesData(ToUse,:);
ScoreScale = AllScaleValuesData2(:,1); BestSpotCodeNo = AllScaleValuesData2(:,2);
BestLogProbOverBackground = AllScaleValuesData2(:,3);
SecondBestLogProbOverBackground = AllScaleValuesData2(:,4);
GadLogProbOverBackground = AllScaleValuesData2(:,5);
NextToGadLogProbOverBackground = AllScaleValuesData2(:,6);
SpotIntensity = AllScaleValuesData2(:,7);
SpotScoreDev = AllScaleValuesData2(:,8); GadRank = AllScaleValuesData2(:,9);
SpotScore = BestLogProbOverBackground - SecondBestLogProbOverBackground;
GadSpotScore = GadLogProbOverBackground - SecondBestLogProbOverBackground;
GadOverNext = GadLogProbOverBackground - NextToGadLogProbOverBackground;

ScoreThreshArray = -50:50;
IntensityThresh2Array = -100:100;
ThresholdQuality = zeros(length(ScoreThreshArray),length(IntensityThresh2Array));
for i=1:length(ScoreThreshArray)
    for j=1:length(IntensityThresh2Array)
    o.pScoreThresh = ScoreThreshArray(i);
    pIntensityThresh2 = IntensityThresh2Array(j);
    %QualOK = GadRank==1;
    QualOK = (GadSpotScore>o.pScoreThresh & GadSpotIntensity>pIntensityThresh2 | ...
        GadSpotIntensity>o.pIntensityThresh & GadLogProbOverBackground>o.pLogProbThresh & GadSpotScore>o.pScoreThresh2 |...
        GadRank==2 & GadOverNext>o.pScoreThresh & GadSpotIntensity>pIntensityThresh2 | ...
        GadOverNext > BestLogProbOverBackground-GadLogProbOverBackground &...
        GadOverNext > o.pScoreThresh & GadSpotIntensity>pIntensityThresh2);
    ThresholdQuality(i,j) = 2*(nFalseNegativesSet-sum(QualOK&FalseNegativesSet))+sum(QualOK&FalsePositivesSet);
    end
end
figure; imagesc(IntensityThresh2Array,ScoreThreshArray,ThresholdQuality);
xlabel('pIntensityThreshold2');
ylabel('pScoreThresh');
        


