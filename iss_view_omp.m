function SpotNo = iss_view_omp(o, FigNo, Norm, SpotNum)
%% SpotNo = iss_view_prob(o, FigNo, Norm, Method, SpotNum)
%
% This function lets you view the spot code, gene code
% and ln(Prob) for the best match for a chosen spot.
%
% o: iss object
% FigNo: figure number (default, current figure)
% Norm: normalization option (described below)
% Method: 'Prob' or 'Pixel' to consider gene assignments given
% by o.pSpotCodeNo or o.pxSpotCodeNo respectively.
% SpotNum: number of spot to analyze (default is to click)
% SpotNo: returns the number of the spot analyzed.
%
% Norm = 1: Raw Colors
% Norm = 2: Z-score each round and colour channel. (SpotColors were in this
% form to find OMP coefficients).
% Norm = 3: Normalised by percentile for each color channel across all
% rounds


%%
if nargin>=4
    SpotNo = SpotNum;
    if nargin>=2 && ishandle(FigNo)
        figure(FigNo);
        xlim([o.pxSpotGlobalYX(SpotNo,2)-20,o.pxSpotGlobalYX(SpotNo,2)+20]);
        ylim([o.pxSpotGlobalYX(SpotNo,1)-20,o.pxSpotGlobalYX(SpotNo,1)+20]);
        rectangle('Position',[o.pxSpotGlobalYX(SpotNo,2)-5,o.pxSpotGlobalYX(SpotNo,1)-5,10,10],'EdgeColor','g','LineWidth',4);
        drawnow;
    end
else
    if nargin>=2
        figure(FigNo);
    end
    CrossHairColor = [1,1,1];   %Make white as black background
    xy = ginput_modified(1,CrossHairColor);
    S = evalin('base', 'issPlot2DObject');
    InRoi = all(int64(round(S.SpotYX))>=S.Roi([3 1]) & round(S.SpotYX)<=S.Roi([4 2]),2);
    PlotSpots = find(InRoi & S.QualOK);         %Only consider spots that can be seen in current plot
    [~,SpotIdx] = min(sum(abs(S.SpotYX(PlotSpots,:)-[xy(2),xy(1)]),2));
    SpotNo = PlotSpots(SpotIdx);      
end

%Different parameters for different methods
CodeNo = o.ompSpotCodeNo(SpotNo);
%SpotColor is raw values (Norm 1)
SpotColor = o.pxSpotColors(SpotNo,:,:);
CodeShape = size(SpotColor);
SpotGlobalYX = o.pxSpotGlobalYX(SpotNo,:);
SpotCoefs = o.ompCoefs(SpotNo,:);
SpotNeighbNonZero = o.ompNeighbNonZeros(SpotNo);
SpotIntensity = o.ompSpotIntensity(SpotNo);
%PredCode is scaled and shifted (Norm 2)
PredCode = SpotCoefs*o.ScaledBledCodes(:,:);
PredCode = reshape(PredCode,CodeShape);
TotalError = sqrt(sum(sum((PredCode - (double(SpotColor)-o.SHIFT)./o.SCALE).^2)));

if nargin<3 || isempty(Norm)
    Norm = 2;
end

%Different Normalisations
if Norm == 1
    cSpotColor = double(SpotColor);
    cPredCode = PredCode.*o.SCALE+o.SHIFT;
elseif Norm == 2    
    cSpotColor = (double(SpotColor)-o.SHIFT)./o.SCALE;
    cPredCode = PredCode;
elseif Norm == 3
    cSpotColor = double(SpotColor);
    cPredCode = PredCode.*o.SCALE+o.SHIFT;
    for b = 1:o.nBP
        bSpotColors = o.pxSpotColors(:,b,:);
        p = double(prctile(bSpotColors(:), o.SpotNormPrctile));
        cSpotColor(:,b,:) = cSpotColor(:,b,:)/p;
        cPredCode(:,b,:) = cPredCode(:,b,:)/p;
    end
end

cSpotColor = squeeze(cSpotColor);
CodeShape = size(cSpotColor);
cPredCode = squeeze(cPredCode);
Error = cSpotColor - cPredCode;
v0 = min(-0.2,min(cSpotColor(:)));
v1 = max(0.2,max(cSpotColor(:)));
caxis_lims = [v0 v1];
error_caxis_lims = [min(v0,min(Error(:))), max(v1,max(Error(:)))];

%Get square outlining unbled code
gUnbled = reshape(o.UnbledCodes(CodeNo(1),:,:),CodeShape);
gSquares = zeros(o.nRounds,4);
for r=1:o.nRounds
    try
        gSquares(r,:) = [r-0.5,find(gUnbled(:,r,:)==1)-0.5,1,1];
    end
end

try
    clf(430476599)
    figure(430476599)
catch
    figure(430476599)
end
subplot(4,1,1);
imagesc(cSpotColor); colorbar('Color','w');
caxis(caxis_lims);
if Norm>1
    sq_color = 'g';
    colormap(gca,bluewhitered);
else
    sq_color = 'r';
end
title('Spot Code','Color','w');
set(gca, 'ytick', 1:o.nBP);
set(gca, 'YTickLabel', o.bpLabels);
ylabel('Color Channel','Color','w');
hold on
for r=1:o.nRounds
    rectangle('Position',gSquares(r,:),'EdgeColor',sq_color,'LineWidth',2,'LineStyle',':')
end
hold off
set(gca,'XColor','w');
set(gca,'YColor','w');

subplot(4,1,2)
imagesc(cPredCode); colorbar('Color','w');
caxis(caxis_lims);
if Norm>1
    colormap(gca,bluewhitered);
end
title(sprintf('Predicted Code for %s, code #%d', o.GeneNames{CodeNo}, CodeNo),'Color','w');
set(gca, 'ytick', 1:o.nBP);
set(gca, 'YTickLabel', o.bpLabels);
ylabel('Color Channel','Color','w');
hold on
for r=1:o.nRounds
    rectangle('Position',gSquares(r,:),'EdgeColor',sq_color,'LineWidth',2,'LineStyle',':')
end
hold off
set(gca,'XColor','w');
set(gca,'YColor','w');

subplot(4,1,3);
imagesc(Error); colorbar('Color','w');
caxis(error_caxis_lims);
colormap(gca,bluewhitered);
title('Error','Color','w');
set(gca, 'ytick', 1:o.nBP);
set(gca, 'YTickLabel', o.bpLabels);
ylabel('Color Channel','Color','w');
xlabel('Round','Color','w');
%title('$\ln\left({\frac{P(spot\,\mid \,gene\,\, and\,\, background)}{P(spot\,\mid \,background)}}\right)$','interpreter','latex','FontSize',15)
%title(sprintf('Log Probability the spot can be explained by gene - Log Probability it can be explained by background alone'));
%set(ClickPlot,'ButtonDownFcn',{@getCoord,o,SpotNo,CodeNo,SpotColor});
hold on
for r=1:o.nRounds
    rectangle('Position',gSquares(r,:),'EdgeColor','g','LineWidth',2,'LineStyle',':')
end
hold off
set(gca,'XColor','w');
set(gca,'YColor','w');

nCodes = size(o.CharCodes,1);
ClickPlot = subplot(4,1,4);
%scatter(1:nCodes,SpotCoefs(1:nCodes),46,[0, 0.4470, 0.7410],'.');
hold on
for i=1:nCodes
    S.h(i) = plot(i, SpotCoefs(i), '.');
end
legend(S.h,o.GeneNames);
legend off
change_gene_symbols(0);
%set(gcf, 'color', 'w');
set(gca, 'color', 'k');
h=findobj(gcf,'type','axes');
delete(h(1));
%scatter(InitialCodeNo,0.8,46,[0.8500, 0.3250, 0.0980],'.');
%scatter(CodeNo,SpotCoefs(CodeNo),66,[0.8500, 0.3250, 0.0980],'.');
scatter(nCodes+1:nCodes+7,SpotCoefs(nCodes+1:nCodes+7),46,'cyan','.');
%ylim([min([backWt-0.05,-0.05,gWt-0.05]),max([backWt+0.05,gWt+0.05,1])]);
ylabel('Coefficient','Color','w');
set(gca, 'xtick', [1:nCodes,nCodes+int8((o.nBP+1)/2)]);
Labels = o.GeneNames;
Labels{end+1} = 'Background';
set(gca, 'XTickLabel', Labels);
set(gca,'XColor','w');
set(gca,'YColor','w');
ClickPlot.XTickLabel{CodeNo} = ['\color{red}',Labels{CodeNo}];
xtickangle(90);
hold off
set(ClickPlot,'ButtonDownFcn',{@getCoord,o,CodeShape,sq_color,'g',ClickPlot.XTickLabel,SpotCoefs});

%Color different parameters depending if over threshold
% if SpotScore>o.pScoreThresh
%     c1 = [0,0.7,0]; else; c1 = [0,0,0];end
% if LogProbOverBackground<o.pLogProbThresh
%     c2 = [1,0,0]; else; c2 = [0,0,0];end
% if SpotScore+SpotScoreDev<o.pDevThresh
%     c3 = [1,0,0]; else; c3 = [0,0,0];end
% if SpotIntensity<o.pIntensityThresh
%     c4 = [1,0,0]; else; c4 = [0,0,0];end
% 
set(gcf,'Position',[409,92,793,707])
figtitle = sgtitle('', 'interpreter', 'tex','Color','w');   %'tex' required for colors
figtitle.String = sprintf('Spot %d: code %d, %s. Coef = %.2f, NeighbNonZero = %d, Error = %.2f, Intensity = %.0f',...
    SpotNo, CodeNo, o.GeneNames{CodeNo}, SpotCoefs(CodeNo), SpotNeighbNonZero,TotalError,SpotIntensity);
%figtitle.Color='red';
%drawnow

fprintf('Spot %d at yx=(%d,%d): code %d, %s\n', ...
    SpotNo, SpotGlobalYX(1),SpotGlobalYX(2),...
    CodeNo, o.GeneNames{CodeNo});

if nargin>=4 && ishandle(FigNo)
    figure(FigNo);
    rect = findall(gcf,'Type', 'Rectangle');
    delete(rect);
end
end

function getCoord(aH,evnt,o,CodeShape,code_sq_color,error_sq_color,original_xtick_labels,Coefs)
%This plots the unbled codes
nCodes = size(o.CharCodes,1);
fig = ancestor(aH,'figure');
click_type = get(fig,'SelectionType');
ClickLoc = evnt.IntersectionPoint(1:2);
CodeNo = round(ClickLoc(1));
if CodeNo>nCodes
    CodeNo = nCodes;
elseif CodeNo < 1
    CodeNo = 1;
end


rect = findall(gcf,'Type', 'Rectangle');
if strcmp(click_type,'normal')        
    %Delete current rectangles    
    delete(rect);
    
    %Get square outlining unbled code
    gUnbled = reshape(o.UnbledCodes(CodeNo,:),CodeShape);
    gSquares = zeros(o.nRounds,4);
    for r=1:o.nRounds
        try
            gSquares(r,:) = [r-0.5,find(gUnbled(:,r,:)==1)-0.5,1,1];
        end
    end
    
    h=findobj(gcf,'type','axes');
    for subplot = 1:3
        set(gcf,'CurrentAxes',h(5-subplot));
        if subplot==3
            sq_color = error_sq_color; else; sq_color = code_sq_color;
        end
        if subplot==2
            title(sprintf('Predicted Code for %s, code #%d. Coef = %.2f', o.GeneNames{CodeNo}, CodeNo, Coefs(CodeNo)));
        end
        hold on
        for r=1:o.nRounds
            rectangle('Position',gSquares(r,:),'EdgeColor',sq_color,'LineWidth',2,'LineStyle',':')
        end
        hold off
    end  
    h(1).XTickLabel = original_xtick_labels;
    h(1).XTickLabel{CodeNo} = ['\color{green}',o.GeneNames{CodeNo}];
    drawnow;    
elseif strcmp(click_type,'alt')    
    sq_colors = {'yellow','gray','cyan','magenta','orange','blue'};
    sq_colors2 = {'yellow',[0.5,0.5,0.5],'cyan','magenta',[0.8500, 0.3250, 0.0980],'blue'};
    color_idx = min(length(rect)/(3*7),length(sq_colors));
    %Get square outlining unbled code
    gUnbled = reshape(o.UnbledCodes(CodeNo,:),CodeShape);
    gSquares = zeros(o.nRounds,4);
    
    for r=1:o.nRounds
        try
            gSquares(r,:) = [r-0.5,find(gUnbled(:,r,:)==1)-0.5,1,1];
        end
    end    
    
    h=findobj(gcf,'type','axes');
    for subplot = 1:3
        set(gcf,'CurrentAxes',h(5-subplot));
        if subplot==2
            title(sprintf('Predicted Code for %s, code #%d. Coef = %.2f', o.GeneNames{CodeNo}, CodeNo, Coefs(CodeNo)));
        end
        hold on
        for r=1:o.nRounds
            rectangle('Position',gSquares(r,:),'EdgeColor',sq_colors2{color_idx},'LineWidth',2,'LineStyle',':');
        end
        hold off
    end  
    h(1).XTickLabel{CodeNo} = ['\color{',sq_colors{color_idx},'}',o.GeneNames{CodeNo}];
    drawnow;        
end
end

