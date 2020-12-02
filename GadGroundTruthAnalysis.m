%% GadGroundTruthNewTransform %Need to do this first
% Or if already done it with one:
% o.GadRawLocalYX = oOMP.GadRawLocalYX;
% o.GadRawIsolated = oOMP.GadRawIsolated;
% o.GadInfo = oOMP.GadInfo;
% o.D0 = oOMP.D0;
% o.D = oOMP.D;
% o.TileOrigin = oOMP.TileOrigin;
% o.BigGadFile = oOMP.BigGadFile;
% o.BigGcampFile = oOMP.BigGcampFile;
% Then need to do last three sections of GadGroundTruthNewTransform. Need
% to load FindSpotsWorkspace to do this.

%%
Method = 'Pixel'; %OMP or Pixel
%Method = 'OMP';
o = GadTruePositiveSets(o,Method);
fprintf('of which, we can achieve %d\n',sum(o.pxGadTruePositiveSet));
fprintf('False positive set has %d spots.\n',sum(o.pxGadFalsePositiveSet));
if strcmpi('OMP',Method)
    %Get primary and secondary sets 
    [~,SortedCoefs]=sort(o.ompCoefs(:,1:73)','descend');
    SortedCoefs = SortedCoefs';
    PrimarySet = o.ompSpotCodeNo==SortedCoefs(:,1);
    SecondarySet = o.ompSpotCodeNo==SortedCoefs(:,2);
    o.ompNeighbThresh = 13;
    o.ompIntensityThresh = 700;
elseif strcmpi('Pixel',Method)
     %Get primary and secondary sets 
    PrimarySet = o.pxSpotScore>0;
    SecondarySet = o.pxSpotScore==0;
    o.pScoreThresh = 0;
    o.pScoreThresh2 = -10;  
    o.pIntensityThresh = 600;
    o.pLogProbThresh2 = 0;
    o.pIntensityThresh2 = 0;
end
QualOK = o.quality_threshold(Method);
%print excel data
fprintf('Total Primary Spots: \t\t\t\t%d\n',sum(QualOK&PrimarySet));
fprintf('Total Primary or Secondary Spots: \t\t%d\n',sum(QualOK&(PrimarySet|SecondarySet)));
fprintf('Total Spots: \t\t\t\t\t%d\n',sum(QualOK));
fprintf('Total Primary True Positives: \t\t\t%d\n',sum(QualOK&PrimarySet&o.pxGadTruePositiveSet));
fprintf('Total Primary or Secondary True Positives: \t%d\n',sum(QualOK&(PrimarySet|SecondarySet)&o.pxGadTruePositiveSet));
fprintf('Total True Positives: \t\t\t\t%d\n',sum(QualOK&o.pxGadTruePositiveSet));
fprintf('Total Primary False Positives: \t\t\t%d\n',sum(QualOK&PrimarySet&o.pxGadFalsePositiveSet));
fprintf('Total Primary or Secondary False Positives: \t%d\n',sum(QualOK&(PrimarySet|SecondarySet)&o.pxGadFalsePositiveSet));
fprintf('Total False Positives: \t\t\t\t%d\n',sum(QualOK&o.pxGadFalsePositiveSet));


%% %Visualisation
% IntensityThresh = 100:100:2000;
% NeighbThresh = 3:35;
% ScoreImage = zeros(length(IntensityThresh),length(NeighbThresh));
% for i=1:length(IntensityThresh)
%     o.ompIntensityThresh = IntensityThresh(i);
%     for j=1:length(NeighbThresh)
%         o.ompNeighbThresh = NeighbThresh(j);
%         QualOK = o.quality_threshold('OMP');
%         ScoreImage(i,j) = sum(QualOK&o.pxGadTruePositiveSet)-0.65*sum(QualOK&o.pxGadFalsePositiveSet);
%     end
% end
% figure; imagesc(ScoreImage);

% ScoreThresh = 0:5:80;
% IntensityThresh = 0:50:1000;
% o.pScoreThresh2 = -0.000001;
% ScoreImage = zeros(length(ScoreThresh),length(IntensityThresh));
% for i=1:length(ScoreThresh)
%     o.pScoreThresh = ScoreThresh(i);
%     for j=1:length(IntensityThresh)
%         o.pIntensityThresh = IntensityThresh(j);
%         QualOK = o.quality_threshold('Pixel');
%         ScoreImage(i,j) = sum(QualOK&o.pxGadTruePositiveSet)-0.65*sum(QualOK&o.pxGadFalsePositiveSet);
%     end
% end
% figure; imagesc(ScoreImage);

%% Testing OMP thresholds
NeighbThresh = 25:40;
IntensityThresh = 2000:50:4000;
ScoreThresh = 5:0.25:15;

TestArray = zeros(length(NeighbThresh),length(IntensityThresh),length(ScoreThresh));

for n=1:length(NeighbThresh)
    for i=1:length(IntensityThresh)
        for s=1:length(ScoreThresh)
            o.ompNeighbThresh=NeighbThresh(n);
            o.ompIntensityThresh=IntensityThresh(i);
            o.ompScoreThresh=ScoreThresh(s);
            QualOK = o.quality_threshold('OMP');
            TestArray(n,i,s) = sum(QualOK&o.pxGadTruePositiveSet)+sum(~QualOK&o.pxGadFalsePositiveSet);
        end
    end
end
[a,b] = max(TestArray(:));
[a,b,c]=ind2sub(size(TestArray),b);
o.ompNeighbThresh=NeighbThresh(a);
o.ompIntensityThresh=IntensityThresh(b);
o.ompScoreThresh = ScoreThresh(c);
NeighbThresh(a)
IntensityThresh(b)
ScoreThresh(c)