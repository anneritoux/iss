function QualOK = quality_threshold_prob(o,Method)
% QualOK = o.quality_threshold
% quick function that returns a binary saying which spots are above quality
% threshold
% If Method == 'Prob', this changes the spots and gene assignments to those given by
% the probability method using reference round. If Method == 'Pixel', use
% spots and gene assignments using pixel based method

if strcmpi('Prob',Method)
    QualOK = (o.pSpotScore>o.pScoreThresh & o.pSpotIntensity>0 | ...
    o.pSpotIntensity>o.pIntensityThresh & o.pLogProbOverBackground>o.pLogProbThresh & o.pSpotScore+o.pSpotScoreDev>o.pDevThresh...
    & o.pSpotScore>0); 
elseif strcmpi('Pixel',Method)
    QualOK = (o.pxSpotScore>o.pScoreThresh & o.pxSpotIntensity>0 | ...
    o.pxSpotIntensity>o.pIntensityThresh & o.pxLogProbOverBackground>o.pLogProbThresh & o.pxSpotScore+o.pxSpotScoreDev>o.pDevThresh...
    & o.pxSpotScore>0); 
else
    error('Method not valid, must be Prob or Pixel');
end
%| o.pSpotIntensity>1000);

% % HACK ALERT
% QualOK = QualOK & o.cSpotIsolated;

nCombiCodes = sum(~strcmp(o.CharCodes, 'EXTRA'));

% now extras - they have their own thresholds, set manually for each type
for i=1:size(o.ExtraCodes,1)
    MySpots = (o.SpotCodeNo == nCombiCodes+i);
    QualOK(MySpots) = o.SpotIntensity(MySpots)>o.ExtraCodes{i,4};
end