function iss_view_codes(o, FigNo)
% iss_view_codes(o, FigNo)
% 
% run this after iss.plot, then click on a gene read to see its actual 
% code value and the bled template for this code.
%
% note it is NOT a member function, you need to provide o as an argument.

    if nargin>=2
        figure(FigNo);
    end
    
    % use normed SpotColors that are actually used
    % to determine spot scores
    SpotColors = bsxfun(@rdivide, o.cSpotColors, prctile(o.cSpotColors, o.SpotNormPrctile));
    FlatSpotColors = SpotColors(:,:);
    o.SpotIntensity = sqrt(sum(FlatSpotColors.^2,2));
    NormFlatSpotColors = bsxfun(@rdivide, FlatSpotColors, o.SpotIntensity);
    cSpotColors = reshape(NormFlatSpotColors,size(o.cSpotColors));
    
    
    set(gca, 'Color', [1 1 1]*.2);
    CrossHairColor = [1,1,1];   %Make white as black background
    xy = ginput_modified(1,CrossHairColor);
    set(gca, 'color', 'k');
    [~,SpotNo] = min(sum(abs(o.SpotGlobalYX-[xy(2),xy(1)]),2));
    CodeNo = o.SpotCodeNo(SpotNo);
    
    MeasuredCode = squeeze(cSpotColors(SpotNo,:,:));
    CodeShape = size(MeasuredCode);
    
    figure(930476530)
    subplot(2,1,1);
    imagesc(MeasuredCode); colorbar
    caxis([0 max(MeasuredCode(:))]);
    title(sprintf('Measured code: match %.3f to %s', o.SpotScore(SpotNo), o.GeneNames{CodeNo}));
    
    set(gca, 'ytick', 1:o.nBP);
    set(gca, 'YTickLabel', o.bpLabels);
    ylabel('Color Channel');
    
    subplot(2,1,2)
    cBledCode = o.NormBledCodes(CodeNo,:);
    imagesc(reshape(cBledCode, CodeShape)); colorbar
    caxis([0 max(cBledCode(:))]);

    title(sprintf('Predicted Code for %s, code #%d (%s)', o.GeneNames{CodeNo}, CodeNo, o.CharCodes{CodeNo}));
    
    
    set(gca, 'ytick', 1:o.nBP);
    set(gca, 'YTickLabel', o.bpLabels);
    ylabel('Color Channel');
    xlabel('Round');

    fprintf('Spot %d at yx=(%d,%d): code %d, %s\n', ...
        SpotNo, o.SpotGlobalYX(SpotNo,1),o.SpotGlobalYX(SpotNo,2), CodeNo, o.GeneNames{CodeNo});

    
end
    