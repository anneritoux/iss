function o = get_TilePos(o, xypos, nTiles)
%% o = o = o.get_TilePos(xypos,nSeries);
% This function takes the tile coordinates and derives the Y,X index of
% each tile. Adds o.TileInitialPosYX and o.TilePosYX to o object.
% o: iss object
% xypos: [X,Y] coordinate of each tile as given from metadata.
% nTiles: number of tiles


% find x and y grid spacing as median of distances that are about
% right
dx = xypos(:,1)-xypos(:,1)'; % all pairs of x distances
xStep = median(dx(abs(1- dx(:)/o.MicroscopeStepSize)<.5));
dy = xypos(:,2)-xypos(:,2)'; % all pairs of y distances
yStep = median(dy(abs(1- dy(:)/o.MicroscopeStepSize)<.5));


% find index for each tile
if isempty(o.TileInitialPosYX)
    if nTiles==1
        o.TileInitialPosYX = [1,1];
    else
        o.TileInitialPosYX = fliplr(1+round((xypos - min(xypos))./[xStep yStep]));
    end
end

% only consider tiles with index less than number of tiles and only 1
% away from another tile index.
uniqueY = unique(o.TileInitialPosYX(:,1),'sorted');
uniqueX = unique(o.TileInitialPosYX(:,2),'sorted');
Y = uniqueY(diff([0;uniqueY])==1 & uniqueY<=nTiles);
X = uniqueX(diff([0;uniqueX])==1 & uniqueX<=nTiles);
MaxY = max(Y);
MaxX = max(X);

%Sometimes get Nan, if only one Nan, then check if all tiles 
%arranged along only one direction i.e. Nan should be 1.
if max(isnan([MaxX,MaxY])) && nanmax(MaxX,MaxY)==nTiles
    if isnan(MaxX); MaxX=1; else; MaxY=1; end
end

if MaxY*MaxX ~= nTiles
    %Hack if one of the dimensions is correct
    Ydiff = unique(diff(Y));
    Xdiff = unique(diff(X));
    if length(Ydiff)==1 && Ydiff==1 && mod(nTiles/MaxY,1)==0
        MaxX = nTiles/MaxY;
    elseif length(Xdiff)==1 && Xdiff==1 && mod(nTiles/MaxX,1)==0
        MaxY = nTiles/MaxX;
    else       
        % save workspace if error for debugging
        save(fullfile(o.OutputDirectory, 'ExtractError_TilePosErrorWorkspace'));
        error('Number of tiles (%d) is not equal to maximum Y position (%d) multiplied by maximum X position (%d)'...
            , nTiles, MaxY, MaxX);
    end
end

%Once found correct MaxY and MaxX, we know what the answer should be.
o.TilePosYX = zeros(size(o.TileInitialPosYX));
TilePosY = flip(repelem(1:MaxY,MaxX));
o.TilePosYX(:,1) = TilePosY;
TilePosX = repmat([flip(1:MaxX),1:MaxX],1,ceil(MaxY/2));
o.TilePosYX(1:nTiles,2) = TilePosX(1:nTiles);
end

