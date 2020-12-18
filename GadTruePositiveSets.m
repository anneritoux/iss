%load('/Users/joshduffield/Documents/UCL/ISS/B9S4_Slice001-002/Slice001output/ByPixelAnchor/NewTransform-PCR6/OMP/oOMP.mat')
%Load in logicals that tell you coordinates of spots that are peaks in Gad
%channel image and are in true positive set
%load('/Users/joshduffield/Documents/UCL/ISS/B9S4_Slice001-002/Slice001output/ByPixelAnchor/NewTransform-PCR6/Gamma/GadGroundTruthLogicals.mat')
%Need to run GadGroundTruthNewTransform first
function o = GadTruePositiveSets(o,Method)
%% o = GadTruePositiveSets(o,Method)
% True positive if spot assigned as Gad and within o.ExtractR2 of Gad peak.
% False positive if spot as Gad but further from Gad peak that o.ExtractR2
% and has Gad Channel intensity below o.GadColorFalsePositiveThresh
% Method = 'OMP' or 'Pixel'

%Get Gad peak locations
GadGroundTruth = GadGroundTruthLogical(o,'Prob');       %High only in Gad Channel of Gad Round
%Use GadPeakColor not o.pGadColor(:,o.GadChannel) as deals with stitched
%regions better.
TruePositiveSet = GadGroundTruth &...         %Exceed thresh intensity and be a Gad gene peak
o.GadPeakColor > o.GadColorTruePositiveThresh & o.GadPeakSpots == 1;
%TruePositiveSet = GadGroundTruth &...         %Exceed thresh intensity and be a Gad gene peak
%o.pGadColor(:,o.GadChannel) > o.GadColorFalsePositiveThresh & o.GadPeakSpots == 1;
GadPeakGlobalYX_TP = o.SpotGlobalYX(TruePositiveSet,:);
fprintf('There are %d Gad peak spots\n', length(GadPeakGlobalYX_TP));

%For FalsePositiveSet, relax the constraint that gad peaks are only high in
%Gad Channel. 
FalsePositiveSet = ...
    o.GadPeakColor > o.GadColorFalsePositiveThresh & o.GadPeakSpots == 1;
GadPeakGlobalYX_FP = o.SpotGlobalYX(FalsePositiveSet,:);

%For each point GadPeakTruePositive, find nearest Gad spot found by pixel
%based method
GadCodeNo = 22;
if strcmpi('OMP',Method)
    SetSize = size(o.ompSpotCodeNo);
    pxGadSpotsIndex = find(o.ompSpotCodeNo==GadCodeNo);
elseif strcmpi('Pixel',Method)
    SetSize = size(o.pxSpotCodeNo);
    pxGadSpotsIndex = find(o.pxSpotCodeNo==GadCodeNo);
end
pxGadSpotsYX = o.pxSpotGlobalYX(pxGadSpotsIndex,:);
treeTruePositive = KDTreeSearcher(pxGadSpotsYX);
[Index2_TP,Dist_TP] = treeTruePositive.knnsearch(GadPeakGlobalYX_TP);

%Spot is true positive if within filter size of a peak in gad channel
%image.
%I.e. of the 2980 true positives in Gad channel image, how many spots
%identified as Gad using OMP are close to these. 
TruePositiveIndex = unique(pxGadSpotsIndex(Index2_TP(Dist_TP<o.GadTruePosMaxSep)));

%False positive if spot is further from GadChannel peak than filter size
%, has GadChannelColor below o.GadColorFalsePositiveThresh in below and assigned as Gad by OMP.
treeFalsePositive = KDTreeSearcher(GadPeakGlobalYX_FP);
[Index2_FP,Dist_FP] = treeFalsePositive.knnsearch(pxGadSpotsYX);
FalsePositiveIndex = unique(pxGadSpotsIndex(Dist_FP>o.GadFalsePosMinSep));

o.pxGadTruePositiveSet = false(SetSize);
o.pxGadTruePositiveSet(TruePositiveIndex) = true;
o.pxGadFalsePositiveSet = false(SetSize);
o.pxGadFalsePositiveSet(FalsePositiveIndex) = true;
o.pxGadFalsePositiveSet = o.pxGadFalsePositiveSet&o.pxGadColor(:,o.GadChannel)<o.GadColorFalsePositiveThresh;

%Save GadTruePositive that we missed
FOUNDTruePositive = Dist_TP<o.GadTruePosMaxSep;
MISSEDTruePositive = Dist_TP>=o.GadTruePosMaxSep;
TruePositiveSetIndex = find(TruePositiveSet);
o.GadPeakFound = zeros(size(o.GadPeakSpots));
o.GadPeakFound(TruePositiveSetIndex(FOUNDTruePositive))=1;
o.GadPeakFound(TruePositiveSetIndex(MISSEDTruePositive))=2;


end