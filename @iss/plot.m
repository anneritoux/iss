function plot(o, BackgroundImageFile, Roi)
% o.plot(BackgroundImage, Roi)
%
% plot the results of in situ sequencing spot detection. 
% SpotYX gives coordinates; Gene is a list of strings; 
% 
% if BackgroundImage is empty or missing, loaded from file o.BigDapiFile
% and subsetted to Roi (if present)
% If it is a numerical array, that is plotted (not subsetted for ROI)
% If zero, nothing is plotted
%
% Roi = [xmin xmax ymin ymax zmin zmax] shows only this part. Whole thing
% shown if empty or missing. Must be integers, xmin and ymin must be 1
%
% All spots within +/- o.PlotZThick of the current Z plane will be shown on that
% Z plane. Default is 0
%
% sizes can be a vector or a scalar - only used for scatter, which isn't
% called anyway.
% 
% Kenneth D. Harris, 29/3/17
% GPL 3.0 https://www.gnu.org/licenses/gpl-3.0.en.html

if nargin<4 || isempty(Roi)
    Roi = round([1, max(o.SpotGlobalYXZ(:,2)), ...
    1, max(o.SpotGlobalYXZ(:,1)),...
    min(o.SpotGlobalYXZ(:,3)), max(o.SpotGlobalYXZ(:,3))]);
end

if Roi(1) ~= 1 || Roi(3) ~= 1
    %Causes bugs if doesn't start at 1
    warning('Set min Roi to 1');
    Roi(1) = 1;
    Roi(3) = 1;
end

if (nargin<2 || isempty(BackgroundImageFile)) && ~isempty(o.BigDapiFile) && ...
        ~isnumeric(o.BigDapiFile)
    if exist(o.BigDapiFile, 'file')
        fprintf('loading background image...');
        BackgroundImageFile = o.BigDapiFile;
        %Load in Dapi image
        Image3D = zeros(Roi(4),Roi(2),Roi(6)-Roi(5)+1,'uint16');
        for z = Roi(5):Roi(6)
        	Image3D(:,:,z-Roi(5)+1) = imread(BackgroundImageFile, z,'PixelRegion', {Roi(3:4), Roi(1:2)});
        end
        fprintf('done\n');
    else
        warning('not sure what to do with BackgroundImage, setting to off');
        Image3D = zeros(Roi(4),Roi(2),Roi(6)-Roi(5)+1,'uint16');
    end
    
elseif ~isempty(BackgroundImageFile) && ~isnumeric(BackgroundImageFile)
    if exist(BackgroundImageFile, 'file')
        %Load in Dapi image
        fprintf('loading background image...');
        Image3D = zeros(Roi(4),Roi(2),Roi(6)-Roi(5)+1,'uint16');
        for z = Roi(5):Roi(6)
            Image3D(:,:,z-Roi(5)+1) = imread(BackgroundImageFile, z,'PixelRegion', {Roi(3:4), Roi(1:2)});
        end
        fprintf('done\n');
    else
        warning('not sure what to do with BackgroundImage, setting to off');
        Image3D = zeros(Roi(4),Roi(2),Roi(6)-Roi(5)+1,'uint16');
    end
    
elseif isnumeric(BackgroundImageFile)
    try
        Image3D = BackgroundImageFile(Roi(3):Roi(4),Roi(1):Roi(2),Roi(5):Roi(6));
    catch
        warning('Background Image wrong size - setting to off');
        Image3D = zeros(Roi(4),Roi(2),Roi(6)-Roi(5)+1,'uint16');
    end
end


try
    S.FigNo = 93454;
    clf(S.FigNo)
catch
    S.FigNo = 93454;
end
S.fh = figure(S.FigNo);set(S.fh,'units','pixels','position',[500 200 800 600]);  %Left, Bottom, Width, Height
set(gcf, 'color', 'k');
set(gca, 'color', 'k');

S.Image = Image3D;
S.Background = imagesc(S.Image(:,:,1)); hold on; colormap bone;
%set(S.Background, 'XData', [Roi(3), Roi(4)]);
%set(S.Background, 'YData', [Roi(1), Roi(2)]);
xlim([Roi(1) Roi(2)]);
ylim([Roi(3) Roi(4)]);
S.MinZ = Roi(5);

title(['Z Plane ' num2str(S.MinZ)],'Color','w');


set(gca, 'YDir', 'normal');
%axis on

S.SpotGeneName = o.GeneNames(o.SpotCodeNo);
S.uGenes = unique(S.SpotGeneName);
% which ones pass quality threshold (combi first)
S.QualOK = o.quality_threshold;
S.SpotYXZ = o.SpotGlobalYXZ;
%S.Roi is the Roi for the current Z plane
S.ZThick = o.PlotZThick;
S.Roi = [Roi(1:4),S.MinZ-S.ZThick,S.MinZ+S.ZThick];
S.ZRange = Roi(6)-Roi(5);
InRoi = all(int64(round(S.SpotYXZ))>=S.Roi([3 1 5]) & round(S.SpotYXZ)<=S.Roi([4 2 6]),2);
PlotSpots = find(InRoi & S.QualOK);
[~, S.GeneNo] = ismember(S.SpotGeneName(PlotSpots), S.uGenes);
S.h = zeros(size(S.uGenes));

%hold on; GET RID OF HOLD AND JUST DELETE PLOTS USING DELETE MEANS THAT THE
%ZOOM IS KEPT BETWEEN Z PLANES
for i=1:length(S.uGenes)
    MySpots = PlotSpots(S.GeneNo==i);
    if any(MySpots)
        S.h(i) = plot(S.SpotYXZ(MySpots,2), S.SpotYXZ(MySpots,1), '.', 'MarkerSize', 1,'Color',hsv2rgb([0 0 0.5]));
    end
end 
%hold off

legend(S.h(S.h~=0), S.uGenes(S.h~=0));
legend off;

set(gca, 'Clipping', 'off');

if ~isempty(PlotSpots)
    change_gene_symbols(0);
else
    set(gcf, 'color', 'k');
    set(gcf, 'InvertHardcopy', 'off');    
end

%If anchor background or all zero then change climits.
if min(S.Background.CData(:))>20 || max(S.Background.CData(:))==0
    caxis([o.TilePixelValueShift,o.TilePixelValueShift+1500]);
end

S.CallMethod = 'DotProduct';
S.CurrentZ = S.MinZ;
assignin('base','issPlot3DObject',S)

S.sl = uicontrol('style','slide',...
                 'unit','pix',...
                 'position',[60 8 693 18],...
                 'min',1,'max',S.ZRange+1,'val',1,...
                 'sliderstep',[1/S.ZRange 1/S.ZRange],...
                 'callback',{@sl_call,S});  
set( findall( S.fh, '-property', 'Units' ), 'Units', 'Normalized' )





function [] = sl_call(varargin)
% Callback for the slider.
[h,~] = varargin{[1,3]};  % calling handle and data structure.
S = evalin('base', 'issPlot3DObject');  %Take S from workspace so can use iss_change_plot
ZPlane = round(get(h,'value'))+S.MinZ-1;
S.Background = imagesc(S.Image(:,:,ZPlane-S.MinZ+1)); hold on; colormap bone;
%set(S.Background, 'XData', [S.Roi(3), S.Roi(4)]);
%set(S.Background, 'YData', [S.Roi(1), S.Roi(2)]);
%xlim([S.Roi(1) S.Roi(2)]);
%ylim([S.Roi(3) S.Roi(4)]);

h = findobj('type','line'); %KEY LINES: DELETE EXISTING SCATTER PLOTS SO CHANGE_SYMBOLS WORKS
delete(h)

set(gca, 'YDir', 'normal');
%axis on
title(['Z Plane ' num2str(ZPlane)],'Color','w');
S.Roi(5:6) = [ZPlane-S.ZThick,ZPlane+S.ZThick];
InRoi = all(round(S.SpotYXZ)>=S.Roi([3 1 5]) & round(S.SpotYXZ)<=S.Roi([4 2 6]),2);
PlotSpots = find(InRoi & S.QualOK);
[~, S.GeneNo] = ismember(S.SpotGeneName(PlotSpots), S.uGenes);
S.h = zeros(size(S.uGenes));
%hold on;
for i=1:length(S.uGenes)
    MySpots = PlotSpots(S.GeneNo==i);
    if any(MySpots)
        S.h(i) = plot(S.SpotYXZ(MySpots,2), S.SpotYXZ(MySpots,1), '.', 'MarkerSize', 1,'Color',hsv2rgb([0 0 0.5]));
    end
end 
%hold off

legend(S.h(S.h~=0), S.uGenes(S.h~=0));
legend off;



set(gca, 'Clipping', 'off');

if ~isempty(PlotSpots)
    change_gene_symbols(0);
else
    set(gcf, 'color', 'k');
    set(gcf, 'InvertHardcopy', 'off');    
end

%Update current Z position in woprkspace so can use for iss_view_codes
S.CurrentZ = ZPlane;
assignin('base','issPlot3DObject',S)








