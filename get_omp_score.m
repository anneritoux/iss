function NormExpScore = get_omp_score(o,SpotNumbers,SpotCodeNo)
%% ModScore = get_omp_score(o)
%omp score is reduction in error caused by introduction of gene in
%rounds/channels where there is not already a gene.

if nargin==1
    nSpots = size(o.ompCoefs,1);
    SpotNumbers = 1:nSpots;
else
    nSpots = length(SpotNumbers);
end
if nargin<3
    SpotCodeNo = o.ompSpotCodeNo(SpotNumbers); 
end

OriginalCoefs = o.ompCoefs(SpotNumbers,:);
SpotColors = o.pxSpotColors(SpotNumbers,:,:);
OriginalPredCodes = OriginalCoefs*o.ScaledBledCodes(:,:);
BestGeneIndex = sub2ind(size(OriginalCoefs),(1:nSpots)',SpotCodeNo);
RemoveGeneCoefs = OriginalCoefs;
RemoveGeneCoefs(BestGeneIndex) = 0;
RemoveGenePredCodes = RemoveGeneCoefs*o.ScaledBledCodes(:,:);
ScaledSpotColors = (double(SpotColors)-o.SHIFT)./o.SCALE;
OriginalError = abs(ScaledSpotColors(:,:)-OriginalPredCodes);
RemoveGeneError = abs(ScaledSpotColors(:,:)-RemoveGenePredCodes); 
%AbsErrorFactor so rounds/channels where error is greater contribute more
AbsErrorFactor = RemoveGeneError./prctile(RemoveGeneError',80)';
AbsErrorFactor(AbsErrorFactor<0.5)=0.5;
AbsErrorFactor(AbsErrorFactor>1.5)=1.5;
OverallScore = RemoveGeneError - OriginalError;

%Get multiplier that only takes into account rounds/channels in unbled code
nCodes = length(o.CharCodes);
GeneMultiplier = zeros(o.nRounds*o.nBP,nCodes);
for g=1:nCodes
    GeneMultiplier(:,g) = o.UnbledCodes(g,:);
end
GeneMultiplier(:,nCodes+1) = 0;     %Use when no overlapping spots
ScoreMultiplier = GeneMultiplier(:,SpotCodeNo);

%Now want to modify this so neglect rounds/channels which already have a
%gene in.
[CoefValues,SortedCoefs]=sort(OriginalCoefs(:,1:73)','descend');
SortedCoefs = SortedCoefs';
CoefValues=CoefValues';
SortedCoefs(CoefValues<=0) = nCodes+1;

%Deal with each overlapping genes in turn until no spot has any more genes
MinCodeNo = 1;
nOverlaps = 0;
ScoreMultiplierOverlaps = zeros(size(ScoreMultiplier));
while MinCodeNo<nCodes+1
    nOverlaps=nOverlaps+1;
    ScoreMultiplierOverlaps = ScoreMultiplierOverlaps+...
        GeneMultiplier(:,SortedCoefs(:,nOverlaps));
    MinCodeNo=min(SortedCoefs(:,nOverlaps));
end
ScoreMultiplier(ScoreMultiplierOverlaps>1&ScoreMultiplier==1)=0.15;
%nRoundChannelsNorm = o.nRounds./sum(ScoreMultiplier==1)';
nRoundChannelsNorm=1;

ModScore = OverallScore.*ScoreMultiplier';
NormModScore = ModScore./RemoveGeneError;
ExpScore = (exp(NormModScore)-1).*AbsErrorFactor;
if nSpots>1
    ExpScore = sum(ExpScore,2);
    NormExpScore = ExpScore.*nRoundChannelsNorm;
else
    NormExpScore = reshape(ExpScore.*nRoundChannelsNorm,o.nRounds,o.nBP);
end
%ModScore = sum(ModScore,2)./nRoundChannelsUsed;
%Fraction of error present without gene is removed due to gene
%ModScore = sum(ModScore,2)./RemoveGeneErrorSum; 
end