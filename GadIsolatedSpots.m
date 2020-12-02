%OriginalIsolated = o.cSpotIsolated;
NonemptyTiles = unique(o.pLocalTile)';
LocalYX = round(o.SpotGlobalYX-o.TileOrigin(o.pLocalTile,:,o.ReferenceRound));
o.cSpotIsolated = double(o.cSpotIsolated);

% annular filter
[xr, yr] = meshgrid(-o.IsolationRadius2:o.IsolationRadius2);
Annulus = (xr.^2 + yr.^2)<=o.IsolationRadius2.^2 & (xr.^2 + yr.^2)>o.IsolationRadius1.^2;
Annulus = double(Annulus)/sum(Annulus(:));
if ischar(o.IsolationThresh)
    IsolationThresh = -250;
else
    IsolationThresh = o.IsolationThresh;
end
for t=NonemptyTiles
    %Load reference round image
    FileName = o.TileFiles{o.ReferenceRound,t};
    TifObj = Tiff(FileName);
    TifObj.setDirectory(o.FirstBaseChannel + o.ReferenceChannel - 1);
    ReferenceIm = int32(TifObj.read())-o.TilePixelValueShift; 
    if o.SmoothSize
        SE = fspecial('disk', o.SmoothSize);
        ReferenceImSm = imfilter(ReferenceIm ,SE);
    else
        ReferenceImSm = ReferenceIm;
    end
    rng(1);     %So shift is always the same.
    RandSmallImShift = rand(o.TileSz,o.TileSz)/10;  %Shift so int value remains the same
    ReferenceImSm = double(ReferenceImSm) + RandSmallImShift;
    
    % filter the image
    AnnularFiltered = imfilter(ReferenceImSm, Annulus);
    % find isolated spots
    UseSpots = o.pLocalTile==t & LocalYX(:,1)>=1 & LocalYX(:,1)<=o.TileSz &...
               LocalYX(:,2)>=1 & LocalYX(:,2)<=o.TileSz;
    MaxPixels = int32(sub2ind(size(AnnularFiltered),LocalYX(UseSpots,1),LocalYX(UseSpots,2)));
    o.cSpotIsolated(UseSpots) = AnnularFiltered(MaxPixels);
    
    if o.Graphics==2
        figure(50965469); clf;
        %imagesc(ReferenceImSm); hold on; colormap hot
        imagesc(AnnularFiltered); hold on; colormap hot
        plot(LocalYX(UseSpots&o.cSpotIsolated<IsolationThresh,2),...
            LocalYX(UseSpots&o.cSpotIsolated<IsolationThresh,1), 'gx');
        plot(LocalYX(UseSpots&o.cSpotIsolated>IsolationThresh,2),...
            LocalYX(UseSpots&o.cSpotIsolated>IsolationThresh,1), 'wx');
        legend('Isolated', 'Not isolated');
        hold off
        drawnow;
    end
end
    
    
    