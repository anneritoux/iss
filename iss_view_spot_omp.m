function iss_view_spot_omp(o, FigNo, ImSz, SpotLocation,ScoreMethod, SpotNum)
%% iss_view_spot(o, FigNo, ImSz, SpotLocation,ScoreMethod, SpotNum)
%
% Check PCR by plotting location of spot in each round and color channel
%
% FigNo: figure number (default, current figure)
% ImSz: radius of image that is plotted for each round and channel.
% Default value is 7 pixels.
% SpotLocation: logical,  if true, will use location of spot closest to
% crosshair, otherwise will use actual position of crosshair. Default is false.
% ScoreMethod: 'DotProduct means use o.SpotCodeNo and ScoreMethod='Prob' means use
% o.pSpotCodeNo to highlight Gene squares. Set to 'Prob' by default
% SpotNum: index of spot that you want to look at.


%%
%If Dist to nearest spot is less than MaxDist, then the crosshair will
%be red on the squares corresponding to the gene the nearest spot
%matches to.
MaxDist = 10;

if nargin<3 || isempty(ImSz)
    ImSz = 7;
end
if ImSz>100
    warning('ImSz too large, setting to 7');
    ImSz = 7;
end

if nargin<4 || isempty(SpotLocation)
    SpotLocation = false;
end

if nargin>=6
    SpotLocation = true;
    SpotNo = SpotNum;
    Dist = 0;
    if strcmpi('Pixel',ScoreMethod)
        xy = o.pxSpotGlobalYX(SpotNo,[2,1]);
    else
        xy = o.SpotGlobalYX(SpotNo,[2,1]);
    end
else
    if nargin>=2
        figure(FigNo);
    end
    CrossHairColor = [1,1,1];   %Make white as black background
    xy = ginput_modified(1,CrossHairColor);
    S = evalin('base', 'issPlot2DObject');
    if nargin<5 || isempty(ScoreMethod)
        ScoreMethod = S.CallMethod;
    elseif ~strcmpi(S.CallMethod,ScoreMethod)
        if strcmpi('Prob',ScoreMethod) || strcmpi('DotProduct',ScoreMethod)
            S.SpotYX = o.SpotGlobalYX;
        elseif strcmpi('Pixel',ScoreMethod)
            S.SpotYX = o.pxSpotGlobalYX;
        end
        S.QualOK = 1;
    end
    InRoi = all(int64(round(S.SpotYX))>=S.Roi([3 1]) & round(S.SpotYX)<=S.Roi([4 2]),2);
    PlotSpots = find(InRoi & S.QualOK);         %Only consider spots that can be seen in current plot
    [Dist,SpotIdx] = min(sum(abs(S.SpotYX(PlotSpots,:)-[xy(2),xy(1)]),2));
    SpotNo = PlotSpots(SpotIdx);
    if SpotLocation || round(Dist)==0
        SpotLocation = true;
        Dist = 0;
        xy = S.SpotYX(SpotNo,[2,1]);
    end
    
end


fprintf('loading channel/round images...');
%Find tile that the point is on and local centered coordinates in reference round
Dist2Tiles = [xy(2),xy(1)]-o.TileOrigin(:,:,o.ReferenceRound);
Dist2Tile = min(sum(Dist2Tiles(Dist2Tiles(:,1)>=0 & Dist2Tiles(:,2)>=0,:),2));
t = find(sum(Dist2Tiles,2)==Dist2Tile);
LocalYX = [xy(2),xy(1)]-o.TileOrigin(t,:,o.ReferenceRound)-o.TileCentre;

S.NonZeroCoefs = [];
if strcmpi(ScoreMethod,'DotProduct')
    numCharCode = str2double(regexp(cell2mat(o.CharCodes(o.SpotCodeNo(SpotNo))),'\d','match'))+1;
elseif strcmpi(ScoreMethod,'Prob')
    numCharCode = str2double(regexp(cell2mat(o.CharCodes(o.pSpotCodeNo(SpotNo))),'\d','match'))+1;   
elseif strcmpi(ScoreMethod,'Pixel')
    numCharCode = str2double(regexp(cell2mat(o.CharCodes(o.pxSpotCodeNo(SpotNo))),'\d','match'))+1;
    S.NonZeroCoefs = find(o.ompCoefs(SpotNo,:)~=0);
    S.SpotCoefs = o.ompCoefs(SpotNo,:);
end
 

%Get spot colors for whole grid
SpotColors = zeros((ImSz*2+1)^2,o.nBP,o.nRounds);
for r=1:o.nRounds  
    for b=1:o.nBP        
        rbYX = round([LocalYX,1]*o.D(:,:,t,r,b)+o.TileCentre);
        y0 = rbYX(1);
        x0 = rbYX(2);
        if y0>o.TileSz || y0<1 || x0>o.TileSz || x0<1
            continue;
        end
        y1 = max(1,y0 - ImSz);
        y2 = min(o.TileSz,y0 + ImSz);
        x1 = max(1,x0 - ImSz);
        x2 = min(o.TileSz,x0 + ImSz);
        BaseIm = int32(imread(o.TileFiles{r,t}, b, 'PixelRegion', {[y1 y2], [x1 x2]}))-o.TilePixelValueShift;
        if o.SmoothSize
            SE = fspecial3('ellipsoid',o.SmoothSize);
            BaseImSm = imfilter(BaseIm, SE);
        else
            BaseImSm = BaseIm;
        end
        SpotColors(:,b,r) = BaseImSm(:);               
    end
end

%Get dot product with the gene bled codes for each pixel in image
S.ImSz = ImSz;
S.ImShape = [2*S.ImSz+1,2*S.ImSz+1];
S.NormSpotColors = (double(SpotColors)-o.SHIFT)./o.SCALE;
S.NormSpotColors = S.NormSpotColors(:,:);
S.NormBledCodes = o.ScaledBledCodes(:,:)./vecnorm(o.ScaledBledCodes(:,:),2,2);
S.DotProduct = S.NormSpotColors*S.NormBledCodes';
S.GeneOrder = [];
S.ResidualOrder = [];
S.nCodes = 80;
S.GeneNames = o.GeneNames;
for g=length(o.CharCodes)+1:S.nCodes
    S.GeneNames{g}='Bckgrnd';
end

try
    clf(276465)
    figure(276465);
catch
    figure(276465);
end
set(gcf,'Position',[164,108,1621,805]);
S.x0 = ImSz+1;
S.y0 = ImSz+1;
%S.climits = [min(S.DotProduct(:)),max(S.DotProduct(:))];
plot_dot_product(S);

fprintf('done\n');
assignin('base','issViewSpotOMPObject',S);
%set(ClickPlot,'ButtonDownFcn',{@getCoord})
end

function plot_dot_product(S)
DotSum = sum(abs(S.DotProduct),1);
DotSum(S.GeneOrder) = -inf;
[~,NextBestGene] = max(DotSum);
S.climits = [min(S.DotProduct(:)),max(S.DotProduct(:))];
for g=1:S.nCodes  
    %Highlight gene with the largest dot product

    
    subplot(10, 8, g);
    DotProductIm = reshape(S.DotProduct(:,g),S.ImShape);
    imagesc(DotProductIm, 'ButtonDownFcn', {@getCoord, g}); hold on;
    caxis(S.climits);
    colormap(gca,bluewhitered);
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);
    

    if ismember(g,S.NonZeroCoefs)
        ax = gca;
        ax.XColor = 'r';
        ax.YColor = 'r';
        plot(xlim, [S.y0 S.y0], 'g'); plot([S.x0 S.x0], ylim, 'g'); 
        if ismember(g,S.GeneOrder)
            title(sprintf('%d: %s, %.2f', g, S.GeneNames{g}, S.SpotCoefs(g)),'FontWeight','normal','Color','r');
        elseif g==NextBestGene
            title(sprintf('%d: %s, %.2f', g, S.GeneNames{g}, S.SpotCoefs(g)),'FontWeight','normal','Color','g');
        else
            title(sprintf('%d: %s, %.2f', g, S.GeneNames{g}, S.SpotCoefs(g)),'FontWeight','normal');
        end
    else
        plot(xlim, [S.y0 S.y0], 'k'); plot([S.x0 S.x0], ylim, 'k');
        if ismember(g,S.GeneOrder)
            title([num2str(g), ': ',S.GeneNames{g}],'FontWeight','normal','Color','r');
        elseif g==NextBestGene
            title([num2str(g), ': ',S.GeneNames{g}],'FontWeight','normal','Color','g');
        else
            title([num2str(g), ': ',S.GeneNames{g}],'FontWeight','normal');
        end
    end
    %drawnow
    set(gca, 'YDir', 'normal');
    hold off

end
cbh = colorbar;
cbh.Position(3) = cbh.Position(3)*2;
cbh.Position(4) = cbh.Position(4)*10;
cbh.Position(1) = .95-cbh.Position(3)/2;
cbh.Position(2) = 0.5-cbh.Position(4)/2;
end

function plot_gene_coefs(S,xFull)
try
    clf(276466)
    figure(276466);
catch
    figure(276466);
end
climits = [min(xFull(:)),max(xFull(:))];
nRows = floor((length(S.GeneOrder)-0.000001)/5)+1;
nCols = min(length(S.GeneOrder),5);
for g=1:length(S.GeneOrder)
    subplot(nRows, nCols, g);
    xIm = reshape(xFull(:,S.GeneOrder(g)),S.ImShape);
    imagesc(xIm);
    hold on
    caxis(climits);
    colormap(gca,bluewhitered);
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);
    plot(xlim, [S.y0 S.y0], 'k'); plot([S.x0 S.x0], ylim, 'k');
    if length(S.ResidualOrder)==length(S.GeneOrder)
        title({sprintf('%d: %s, %.2f', S.GeneOrder(g), S.GeneNames{S.GeneOrder(g)}, xIm(S.x0,S.y0)),...
            sprintf('Residual = %.2f',S.ResidualOrder(g))},'FontWeight','normal');
    else
        title(sprintf('%d: %s, %.2f', S.GeneOrder(g), S.GeneNames{S.GeneOrder(g)}, xIm(S.x0,S.y0)),...
            'FontWeight','normal');
    end
    set(gca, 'YDir', 'normal');
    hold off
end
cbh = colorbar;
cbh.Position(3) = cbh.Position(3)*1;
cbh.Position(4) = cbh.Position(4)*nRows;
cbh.Position(1) = .95-cbh.Position(3)/2;
cbh.Position(2) = 0.5-cbh.Position(4)/2;
set(gcf,'Position',[50,500,290+90*nCols,250+100*nRows]);
end


function S = getCoord(aH,evnt,idx)
fig = ancestor(aH,'figure');
click_type = get(fig,'SelectionType');
S = evalin('base', 'issViewSpotOMPObject');
if strcmp(click_type,'normal')
    S.GeneOrder = idx;    
elseif strcmp(click_type,'alt')
    if ~ismember(idx,S.GeneOrder)
        S.GeneOrder = [S.GeneOrder, idx];
    end
end
NormSpotColorsNew = S.NormSpotColors;
xFull = zeros(S.ImShape(1)*S.ImShape(2),80);
residual = zeros(S.ImShape(1)*S.ImShape(2),1);
for b=1:S.ImShape(1)*S.ImShape(2)
    [x,r,residual(b)] = omp_specify_atoms(S.NormBledCodes',S.NormSpotColors(b,:)',S.GeneOrder);
    NormSpotColorsNew(b,:) = r';
    xFull(b,:) = x';
end
if strcmp(click_type,'normal')
    S.ResidualOrder = mean(residual);
elseif strcmp(click_type,'alt')
    if length(S.ResidualOrder) == length(S.GeneOrder)-1
        S.ResidualOrder = [S.ResidualOrder,mean(residual)];
    end
end
S.DotProduct = NormSpotColorsNew*S.NormBledCodes';
plot_dot_product(S);
plot_gene_coefs(S,xFull);

assignin('base','issViewSpotOMPObject',S);
%evnt.Source.Position
%fig = ancestor(aH,'figure');
%click_type = get(fig,'SelectionType');
%ClickLoc = evnt.IntersectionPoint(1:2);
end

