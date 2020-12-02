%% Saving raw data
Method = 'Pixel';
%Method = 'Pixel';
%Use spots which have GadGroundTruth data to test
nSpots = sum(o.pxGadTruePositiveSet)+sum(o.pxGadFalsePositiveSet);
SpotNumbers = zeros(nSpots,1);
SpotNumbers(1:sum(o.pxGadTruePositiveSet))=find(o.pxGadTruePositiveSet);
SpotNumbers(sum(o.pxGadTruePositiveSet)+1:end)=find(o.pxGadFalsePositiveSet);
SpotIdentity = SpotNumbers;
SpotIdentity(1:sum(o.pxGadTruePositiveSet))=1;
SpotIdentity(sum(o.pxGadTruePositiveSet)+1:end)=2;

%Need means and std to z score data
if strcmpi('OMP',Method)
    IntensityMean = mean(o.ompSpotIntensity);
    IntensityStd = std(o.ompSpotIntensity);
    Intensity2Mean = mean(o.ompSpotIntensity2);
    Intensity2Std = std(o.ompSpotIntensity2);
    NeighboursMean = mean(double(o.ompNeighbNonZeros));
    NeighboursStd = std(double(o.ompNeighbNonZeros));
    ScoreMean = mean(o.ompScore);
    ScoreStd = std(o.ompScore);
    
    nData = 4;
    SpotData = zeros(nSpots,nData);
    SpotData(:,1)=(o.ompSpotIntensity(SpotNumbers)-IntensityMean)/IntensityStd;
    SpotData(:,2)=(o.ompSpotIntensity2(SpotNumbers)-Intensity2Mean)/Intensity2Std;
    SpotData(:,3)=(double(o.ompNeighbNonZeros(SpotNumbers))-NeighboursMean)/NeighboursStd;
    SpotData(:,4)=(o.ompScore(SpotNumbers)-ScoreMean)/ScoreStd;  
    
    OutputDirectory = '/Users/joshduffield/Documents/UCL/ISS/B9S4_Slice001-002/Slice001output/ByPixelAnchor/NewTransform-PCR6/OMP/Without1stZPlane/PyTorch/OMP/';
    
elseif strcmpi('Pixel',Method)
    IntensityMean = mean(o.pxSpotIntensity);
    IntensityStd = std(o.pxSpotIntensity);
    Intensity2Mean = mean(o.pxSpotIntensity2);
    Intensity2Std = std(o.pxSpotIntensity2);  
    LogProbMean = mean(o.pxLogProbOverBackground);
    LogProbStd = std(o.pxLogProbOverBackground);
    ScoreMean = mean(o.pxSpotScore);
    ScoreStd = std(o.pxSpotScore);    
    ScoreDevMean = mean(o.pxSpotScoreDev);
    ScoreDevStd = std(o.pxSpotScoreDev);
    ScoreSum = o.pxLogProbOverBackground+o.pxSpotScore;
    ScoreSumMean = mean(ScoreSum);
    ScoreSumStd = std(ScoreSum);
    
    nData = 6;
    SpotData = zeros(nSpots,nData);
    SpotData(:,1)=(o.pxSpotIntensity(SpotNumbers)-IntensityMean)/IntensityStd;
    SpotData(:,2)=(o.pxSpotIntensity2(SpotNumbers)-Intensity2Mean)/Intensity2Std;
    SpotData(:,2)=(o.pxLogProbOverBackground(SpotNumbers)-LogProbMean)/LogProbStd;
    SpotData(:,4)=(o.pxSpotScore(SpotNumbers)-ScoreMean)/ScoreStd;  
    SpotData(:,5)=(o.pxSpotScoreDev(SpotNumbers)-ScoreDevMean)/ScoreDevStd; 
    SpotData(:,6)=(ScoreSum(SpotNumbers)-ScoreSumMean)/ScoreSumStd; 
    
    OutputDirectory = '/Users/joshduffield/Documents/UCL/ISS/B9S4_Slice001-002/Slice001output/ByPixelAnchor/NewTransform-PCR6/OMP/Without1stZPlane/PyTorch/PixelBased/';
end

%Save ZScored data
writeNPY(SpotNumbers,fullfile(OutputDirectory,'SpotNumbers.npy'));
writeNPY(SpotData,fullfile(OutputDirectory,'SpotData.npy'));
writeNPY(SpotIdentity,fullfile(OutputDirectory,'SpotIdentity.npy'));

%% Do Python bit
%/Users/joshduffield/Documents/UCL/ISS/B9S4_Slice001-002/Slice001output/ByPixelAnchor/NewTransform-PCR6/OMP/Without1stZPlane/PyTorch/ThresholdingNN.ipynb

%% Load in model
NN=load('/Users/joshduffield/Documents/UCL/ISS/B9S4_Slice001-002/Slice001output/ByPixelAnchor/NewTransform-PCR6/OMP/Without1stZPlane/PyTorch/PixelBased/PixelBased_NeuralNetwork_Train_Full_Set.mat');
nSpotsAll = size(o.pxSpotColors,1);
nData = size(NN.NN_layer_values{1},2);


SpotDataAll = zeros(nSpotsAll,nData);
if strcmpi('OMP',Method)
    SpotDataAll(:,1)=(o.ompSpotIntensity-IntensityMean)/IntensityStd;
    SpotDataAll(:,2)=(o.ompSpotIntensity2-Intensity2Mean)/Intensity2Std;
    SpotDataAll(:,3)=(double(o.ompNeighbNonZeros)-NeighboursMean)/NeighboursStd;
    SpotDataAll(:,4)=(o.ompScore-ScoreMean)/ScoreStd;
elseif strcmpi('Pixel',Method)
    SpotDataAll(:,1)=(o.pxSpotIntensity-IntensityMean)/IntensityStd;
    SpotDataAll(:,2)=(o.pxSpotIntensity2-Intensity2Mean)/Intensity2Std;
    SpotDataAll(:,2)=(o.pxLogProbOverBackground-LogProbMean)/LogProbStd;
    SpotDataAll(:,4)=(o.pxSpotScore-ScoreMean)/ScoreStd;  
    SpotDataAll(:,5)=(o.pxSpotScoreDev-ScoreDevMean)/ScoreDevStd; 
end
    

Layer1 = SpotDataAll*NN.NN_layer_values{1}'+NN.NN_layer_values{2};
Layer1(Layer1<0)=0;     %Relu
Layer2 = Layer1*NN.NN_layer_values{3}'+NN.NN_layer_values{4};
Layer2 = exp(Layer2)./sum(exp(Layer2),2);     %SoftMax

