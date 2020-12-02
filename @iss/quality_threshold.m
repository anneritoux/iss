function QualOK = quality_threshold(o,Method)
% QualOK = o.quality_threshold
% quick function that returns a binary saying which spots are above quality
% threshold
%Method = 'DotProduct','Prob' or 'Pixel' to consider gene assignments given
%by o.SpotCodeNo, o.pSpotCodeNo and o.pxSpotCodeNo respectively.


if strcmpi('Prob',Method)
    %QualOK = (o.pSpotScore>o.pScoreThresh & o.pSpotIntensity>o.pIntensityThresh2 | ...
    %o.pSpotScore>o.pScoreThresh2 & o.pSpotScore+o.pLogProbOverBackground>o.pLogProbThresh2 &...
    %o.pSpotIntensity>o.pIntensityThresh); 
    %QualOK = QualOK & o.pSpotIntensity2 > o.pIntensity2Thresh;
    QualOK = o.pSpotScore>0 & o.pLogProbOverBackground+o.pQualParam1*o.pSpotScore>o.pQualThresh1 | ...
        o.pSpotScore==0 & o.pLogProbOverBackground+o.pQualParam2*o.pSpotScore>o.pQualThresh2;
elseif strcmpi('Pixel',Method)
    %QualOK = (o.pxSpotScore>o.pScoreThresh & o.pxSpotIntensity>o.pIntensityThresh2 | ...
    %o.pxSpotScore>o.pScoreThresh2 & o.pxSpotScore+o.pxLogProbOverBackground>o.pLogProbThresh2 &...
    %o.pxSpotIntensity>o.pIntensityThresh); 
    %QualOK = QualOK & o.pxSpotIntensity2 > o.pIntensity2Thresh;
    QualOK = o.pxSpotScore>0 & o.pxLogProbOverBackground+o.pQualParam1*o.pxSpotScore>o.pQualThresh1 | ...
        o.pxSpotScore==0 & o.pxLogProbOverBackground+o.pQualParam2*o.pxSpotScore>o.pQualThresh2;
% Extra last condition is for overlapping spots
%| o.pSpotIntensity>1000);
elseif strcmpi('DotProduct',Method)
    QualOK = (o.SpotCombi & o.SpotScore>o.CombiQualThresh & o.SpotIntensity>o.CombiIntensityThresh & o.SpotScoreDev>o.CombiDevThresh);
elseif strcmpi('OMP',Method)
     %Old method below
     %QualOK = o.ompNeighbNonZeros>o.ompNeighbThresh | (o.ompSpotIntensity>o.ompIntensityThresh & o.ompNeighbNonZeros>o.ompNeighbThresh2);
     %QualOK = QualOK & o.ompSpotIntensity2 > o.ompIntensity2Thresh;
     %New method, found using PyTorch
     QualOK = o.ompNeighbNonZeros>o.ompNeighbThresh | o.ompSpotIntensity>o.ompIntensityThresh | o.ompScore>o.ompScoreThresh;
else
    error('Method not valid, must be DotProduct, Prob or Pixel');
end


% % HACK ALERT
% QualOK = QualOK & o.cSpotIsolated;

nCombiCodes = sum(~strcmp(o.CharCodes, 'EXTRA'));

% now extras - they have their own thresholds, set manually for each type
for i=1:size(o.ExtraCodes,1)
    MySpots = (o.SpotCodeNo == nCombiCodes+i);
    QualOK(MySpots) = o.SpotIntensity(MySpots)>o.ExtraCodes{i,4};
end