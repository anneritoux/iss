function o = call_spots_prob(o)
%% o = o.call_spots
% calls spots to codes for in-situ sequencing. Run this after find_spots
% 
% produces SpotGene{Spot}: name of gene for each spot
% SpotCode{Spot}: text representation of code for each spot 
% SpotScore(Spot): score saying how well the code fits (0...1)
% SpotIntensity(Spot): RMS intensity of the spot
% 
% Using o.UseChannels and o.UseRounds, you can do spot calling
% without using certain rounds and colour channels.
%
% Kenneth D. Harris, 29/3/17
% GPL 3.0 https://www.gnu.org/licenses/gpl-3.0.en.html

%% Logging
if o.LogToFile
    diary(o.LogFile);
    cleanup = onCleanup(@()diary('off'));
end

%% Make Bleed Matrices
%Only using channels and rounds given by o.UseChannels and o.UseRounds
if isempty(o.UseChannels)
    o.UseChannels = 1:o.nBP;
end
    
if isempty(o.UseRounds)
    o.UseRounds = 1:o.nRounds;
end

nChans = size(o.UseChannels,2);
nRounds = size(o.UseRounds,2);

%Normalise each colour channel by a percentile as to correct for weaker
%colour channels
if strcmpi(o.BleedMatrixType,'Separate')
    p = prctile(o.cSpotColors, o.SpotNormPrctile);
elseif strcmpi(o.BleedMatrixType,'Single')
    p = zeros(1,o.nBP,o.nRounds);
    for b = 1:o.nBP
        bSpotColors = o.cSpotColors(:,b,:);
        p(:,b,:) = prctile(bSpotColors(:), o.SpotNormPrctile);
    end    
else
    warning('Wrong o.BleedMatrixType entry, should be either Separate or Single')
end

NormBleedMatrix = o.BleedMatrix;
%Unnormalise Bleed matrices by multiplying rows by percentiles
BleedMatrix = zeros(o.nBP,o.nBP,nRounds);
for r=1:nRounds
    for b=1:o.nBP
        BleedMatrix(b,:,r) = p(:,b,o.UseRounds(r))*NormBleedMatrix(b,:,r);
    end
end

if o.Graphics
    figure(98043764); clf
    for i=1:nRounds
        subplot(ceil(nRounds/3),3,i); 
        imagesc(BleedMatrix(:,:,i)); 
        %caxis([0 1]); 
        title(sprintf('Cycle %d', o.UseRounds(i))); 
        set(gca, 'xtick', 1:nChans);
        set(gca, 'XTickLabel', o.bpLabels(o.UseChannels));
        set(gca, 'ytick', 1:nChans);
        set(gca, 'yTickLabel', o.bpLabels(o.UseChannels));
        if i==4
            xlabel('Actual')
            ylabel('Measured');
        end
    end
%     subplot(2,3,6);
%     caxis([0 1]); 
%     axis off
%     colormap hot
% %     colorbar
end
o.pBleedMatrix = BleedMatrix;

%% Get BledCodes for each gene
% now load in the code book and apply bleeds to it
%codebook_raw = importdata(o.CodeFile);
%CharCode = codebook_raw.textdata(2:end,5);
%GeneName = codebook_raw.textdata(2:end,3);
GeneName = {};
CharCode = {};
fp = fopen(o.CodeFile, 'r');
tmp = textscan(fp, '%s %s', inf);
GeneName=tmp{1};
CharCode=tmp{2};
fclose(fp);

% bit of a hack to get rid of Sst and Npy (assume always in the end)
nCodes = size(CharCode,1) - nnz(cellfun(@(v) strcmp(v(1:2),'SW'), CharCode));

% put them into object o but without the extras
o.CharCodes=CharCode(1:nCodes);
o.GeneNames=GeneName(1:nCodes);

% create numerical code (e.g. 33244 for CCGAA)
NumericalCode = zeros(nCodes, o.nRounds);
for r=1:o.nRounds
    if r<=o.nRounds-o.nRedundantRounds
        for c=1:nCodes
            [~, NumericalCode(c,r)] = ismember(CharCode{c}(r), o.bpLabels);
        end
    else
        % redundant round - compute codes automatically
        % find pseudobases for this code
        for c=1:nCodes
            PseudoCode = repmat('0',1,o.nRounds-o.nRedundantRounds);
            for p = 1:length(o.RedundantPseudobases)
                PseudoCode(1,ismember(CharCode{c}, o.RedundantPseudobases{p}))=('0'+p);
            end
            % now match them to the redundant codes
            for cc=1:nChans
                rrn = r-o.nRounds+o.nRedundantRounds;
                if ~isempty(regexp(PseudoCode, o.RedundantCodes{rrn,cc}, 'once'))
                    NumericalCode(c,r)=cc;
                end
            end
        end
    end
end

BledCodes = nan(nCodes, o.nBP*o.nRounds);
UnbledCodes = zeros(nCodes, o.nBP*o.nRounds);
% make starting point using bleed vectors (means for each base on each day)
for i=1:nCodes
    for r=1:nRounds
        if any(o.UseChannels == NumericalCode(i,o.UseRounds(r))) == 0 continue; end
        BledCodes(i,o.UseChannels+o.nBP*(r-1)) = BleedMatrix(o.UseChannels, NumericalCode(i,o.UseRounds(r)), r);
        UnbledCodes(i,NumericalCode(i,o.UseRounds(r))+o.nBP*(r-1)) = 1;
    end
end

o.pBledCodes = BledCodes;
o.UnbledCodes = UnbledCodes;

%% Find probability of each spot to each gene
%Comes from P(f) = P_lambda(lambda)P_hist(f-lambda*g) as f = lambda*g+background

%Load histogram data - background prob distribution
%Need to make hist data symmetric and include all data - 0 in middle
%This assumes -NewValMax < min(o.HistValues).
[NonZeroValIdx,~,~] = ind2sub(size(o.HistCounts),find(o.HistCounts>0));
NewValMax = max(max(o.cSpotColors(:)),o.HistValues(max((NonZeroValIdx))));
o.SymmHistValues = -NewValMax:NewValMax;
nBins = length(o.SymmHistValues);
nPixels = sum(o.HistCounts(:,1,1));
SymmHistCounts = zeros(nBins,o.nBP,o.nRounds);
LastIdx = find(o.HistValues==o.SymmHistValues(end));
if o.SymmHistValues(1)<o.HistValues(1)
    FirstIdx = find(o.SymmHistValues==o.HistValues(1));
    SymmHistCounts(FirstIdx:end,:,:) = o.HistCounts(1:LastIdx,:,:);
else
    FirstIdx = find(o.HistValues==o.SymmHistValues(1));
    SymmHistCounts(:,:,:) = o.HistCounts(FirstIdx:LastIdx,:,:);
end
HistProbs = SymmHistCounts/nPixels;
o.HistProbs = (HistProbs+o.alpha)./(1+nBins*o.alpha);

%Get Lambda probability distribution for all genes
x = min(o.cSpotColors(:))-1:max(o.cSpotColors(:))-1;    %subsitiution x=lambda*g, -1 due to matlab indexing
o.ZeroIndex = find(x==0);     %need when looking up conv values
o.LambdaDist = zeros(length(x),nCodes,o.nBP,o.nRounds);

for GeneNo = 1:nCodes
    BledCode = reshape(o.pBledCodes(GeneNo,:),[o.nBP,o.nRounds]);
    numCharCode = str2double(regexp(cell2mat(o.CharCodes(GeneNo)),'\d','match'))+1;

    for b=o.UseChannels
        for r=1:o.nRounds
            g = BledCode(b,r);
            if numCharCode(r)==b
                %for b/r in CharCodes, expect non zero lambda.
                %g always >0 in this case
                o.LambdaDist(x>0,GeneNo,b,r) = raylpdf(x(x>0)/g,o.RaylConst)/g;
                o.LambdaDist(:,GeneNo,b,r) = o.LambdaDist(:,GeneNo,b,r)/sum(o.LambdaDist(:,GeneNo,b,r));
            else
                %for b/r not in CharCodes, expect approx zero lambda.
                o.LambdaDist(:,GeneNo,b,r) = (o.ExpConst/2)*exp(-o.ExpConst*abs(x/g))/abs(g);
                o.LambdaDist(:,GeneNo,b,r) = o.LambdaDist(:,GeneNo,b,r)/sum(o.LambdaDist(:,GeneNo,b,r));
            end
        end
    end
end

%Store convolution results as look up table
LookupTable = zeros(length(x),nCodes,o.nBP,o.nRounds);
fprintf('\nDoing convolutions for Channel  ');
for b=o.UseChannels
    fprintf('\b%d',b);
    for r=1:o.nRounds
        LookupTable(:,:,b,r)=log(conv2(o.LambdaDist(:,:,b,r),o.HistProbs(:,b,r),'same'));
    end
end
fprintf('\n');

%Get log probs for each spot 
nSpots = size(o.cSpotColors,1);
LogProb = zeros(nSpots,nCodes);
BackgroundLogProb = zeros(nSpots,1);

gChannelIndex = repmat(o.UseChannels,1,nRounds);
gRoundIndex = repelem(o.UseRounds,1,nChans);
ChannelIndex = repmat(gChannelIndex,1,nCodes);
RoundIndex = repmat(gRoundIndex,1,nCodes);
GeneIndex = repelem(1:nCodes,1,nRounds*nChans);
HistZeroIndex = find(o.SymmHistValues == 0);            %As HistProbs is of different length to x
SpotColors = o.cSpotColors(:,o.UseChannels,o.UseRounds);
for s=1:nSpots
    SpotIndex = repmat(o.ZeroIndex-1+SpotColors(s,:),1,nCodes); %-1 due to matlab indexing I think
    Indices = sub2ind(size(LookupTable),SpotIndex,GeneIndex,ChannelIndex,RoundIndex);
    LogProb(s,:)=sum(reshape(LookupTable(Indices),[nRounds*nChans,nCodes]));
    BackgroundIndices = sub2ind(size(o.HistProbs),HistZeroIndex+SpotColors(s,:),gChannelIndex,gRoundIndex);
    BackgroundLogProb(s) = sum(log(o.HistProbs(BackgroundIndices)));
end
[LogProb,SpotCodeNo] = sort(LogProb,2,'descend');

o.pLogProbOverBackground = LogProb(:,1)-BackgroundLogProb;
o.pSpotCodeNo = SpotCodeNo(:,1);
o.pSpotScore = LogProb(:,1)-LogProb(:,2);
%Store deviation in spot scores - can rule out matches based on a low
%deviation.
o.pSpotScoreDev = std(LogProb,[],2);
o.pSpotIntensity = o.get_spot_intensity(o.pSpotCodeNo);
end
