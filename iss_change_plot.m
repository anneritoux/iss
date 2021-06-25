function iss_change_plot(o,Method,GeneType,GenesToShow,UseSpots)
%Given issPlot3DObject, this function lets you change the details
%of the plot without closing the figure e.g. you can change
%o.CombiQualThresh or issPlot3DObject.ZThick to change the threshold value and the number
%of Z planes you see. 
%If Method == 'Prob', this changes the gene assignments to those given by
%the probability method rather than the dot product. In this case
%o.pScoreThresh is the threshold.
% GeneType: Neuron or Non-Neuron
%GenesToShow is a cell of gene names that you want to see e.g.
%[{'Npy'},{'Pvalb'}]. It is case sensitive.
% UseSpots: if you want to use your own thresholding, not
% o.quality_threshold. Logical array e.g. o.SpotScore>0.

S = evalin('base', 'issPlot3DObject');
S.ZThick = o.PlotZThick;
figure(S.FigNo);
h = findobj('type','line'); %KEY LINES: DELETE EXISTING SCATTER PLOTS SO CHANGE_SYMBOLS WORKS
delete(h);

if nargin<3 || isempty(GeneType)
    if ~isfield(S,'GeneType')
        GeneType = 'Neuron';
    else
        GeneType = S.GeneType;
    end
end
if ~strcmpi(GeneType,'Neuron') && ~strcmpi(GeneType,'NonNeuron')
    warning('Showing neuron type Genes');
    GeneType = 'Neuron';
end
S.GeneType = GeneType;

if nargin<4 || isempty(GenesToShow)
    GenesToShow = o.GeneNames;
    if ~isfield(S,'GeneNoToShow')
        %Only change if not previosuly given GenesToShow
        S.GeneNoToShow = find(ismember(o.GeneNames,GenesToShow));
    end
else
    if strcmpi(GenesToShow,'Neuron')
        GeneNames = change_gene_symbols(0);
        GenesToShow = GeneNames(:,1);
    elseif strcmpi(GenesToShow,'NonNeuron')
        GeneNames = change_gene_symbols_NonNeuron(0);
        GenesToShow = GeneNames(:,1);
    end
    S.GeneNoToShow = find(ismember(o.GeneNames,GenesToShow));
end

if nargin<2 || isempty(Method)  
    if strcmpi(S.CallMethod,'DotProduct')
        S.QualOK = o.quality_threshold & ismember(o.SpotCodeNo,S.GeneNoToShow);
    elseif strcmpi(S.CallMethod,'Prob')
        S.QualOK = o.quality_threshold_prob & ismember(o.pSpotCodeNo,S.GeneNoToShow);
    end
else
    if strcmpi('Prob',Method)
        S.SpotGeneName = o.GeneNames(o.pSpotCodeNo);
        S.uGenes = unique(S.SpotGeneName);
        % which ones pass quality threshold (combi first)
        S.QualOK = o.quality_threshold_prob & ismember(o.pSpotCodeNo,S.GeneNoToShow);
        S.CallMethod = 'Prob';
    else
        S.SpotGeneName = o.GeneNames(o.SpotCodeNo);
        S.uGenes = unique(S.SpotGeneName);
        % which ones pass quality threshold (combi first)
        S.QualOK = o.quality_threshold & ismember(o.SpotCodeNo,S.GeneNoToShow);
        S.CallMethod = 'DotProduct';
    end
end

if nargin>=5 && length(UseSpots)==length(o.SpotCodeNo) && islogical(UseSpots)
    S.QualOK = UseSpots & ismember(o.SpotCodeNo,S.GeneNoToShow);
else
    if nargin>=5; warning('UseSpots not valid, using quality_threshold');end
    if strcmpi('Prob',S.CallMethod)
        S.QualOK = o.quality_threshold_prob & ismember(o.SpotCodeNo,S.GeneNoToShow);
    else
        S.QualOK = o.quality_threshold & ismember(o.SpotCodeNo,S.GeneNoToShow);
    end
end

%S.SpotYXZ = o.SpotGlobalYXZ;
%S.Roi is the Roi for the current Z plane
S.Roi(5:6) = [S.CurrentZ-S.ZThick,S.CurrentZ+S.ZThick];
InRoi = all(int64(round(S.SpotYXZ))>=S.Roi([3 1 5]) & round(S.SpotYXZ)<=S.Roi([4 2 6]),2);
PlotSpots = find(InRoi & S.QualOK);
[~, S.GeneNo] = ismember(S.SpotGeneName(PlotSpots), S.uGenes);
S.h = zeros(size(S.uGenes));

%hold on; GET RID OF HOLD AND JUST DELETE PLOTS USING DELETE MEANS THAT THE
%ZOOM IS KEPT BETWEEN Z PLANES
for i=1:length(S.uGenes)
    MySpots = PlotSpots(S.GeneNo==i);
    if any(MySpots)
        S.h(i) = plot(S.SpotYXZ(MySpots,2), S.SpotYXZ(MySpots,1), '.', 'MarkerSize',...
            1,'Color',hsv2rgb([0 0 0.5]));
    end
end 
%hold off

legend(S.h(S.h~=0), S.uGenes(S.h~=0));
legend off;

set(gca, 'Clipping', 'off');

if ~isempty(PlotSpots)
    if strcmpi(GeneType,'Neuron')
        change_gene_symbols(0);
    elseif strcmpi(GeneType,'NonNeuron')
        change_gene_symbols_NonNeuron(0);
    end
else
    set(gcf, 'color', 'k');
    set(gcf, 'InvertHardcopy', 'off');    
end

S.CurrentZ = S.Roi(5)+S.ZThick;
assignin('base','issPlot3DObject',S)
