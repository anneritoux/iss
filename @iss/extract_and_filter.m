function o = extract_and_filter(o)
% create tiff files for each tile that are top-hat filtered versions of
% original czi files

    o.TileFiles = cell(o.nRounds+o.nExtraRounds,1,1); % 1,1 because we don't yet know how many tiles

    for r = 1:o.nRounds+o.nExtraRounds       
        imfile = fullfile(o.InputDirectory, [o.FileBase{r}, o.RawFileExtension]);

        % construct a Bio-Formats reader with the Memoizer wrapper
        bfreader = loci.formats.Memoizer(bfGetReader(), 0);
        % initiate reader
        bfreader.setId(imfile);

        % get some basic image metadata
        [nSeries, nSerieswPos, nChannels, nZstacks, xypos, pixelsize] = ...
            get_ome_tilepos(bfreader);
        if isempty(xypos) || size(xypos, 1)==1
            if r == 1
                warning('first round xypos empty - using values from initial manual input')
                assert(~isempty(o.TileInitialPosXY), 'xypos unavailable')
                xypos = o.TileInitialPosXY;
                xyposOld = xypos;
            else
                warning('xypos empty - using values from previous round')
                xypos = xyposOld;
            end
            nSerieswPos = size(xypos,1);
        else
            xyposOld = xypos;
        end
        
        scene = nSeries/nSerieswPos;

        bfreader.close();
        
        if r == 1
            % find x and y grid spacing as median of distances that are about
            % right
            dx = xypos(:,1)-xypos(:,1)'; % all pairs of x distances
            xStep = median(dx(abs(1- dx(:)/o.MicroscopeStepSize)<.5));
            dy = xypos(:,1)-xypos(:,1)'; % all pairs of y distances
            yStep = median(dy(abs(1- dy(:)/o.MicroscopeStepSize)<.5));
        
        
            % find coordinates for each tile
            o.TileInitialPosYX = fliplr(1+round((xypos - min(xypos))./[xStep yStep]));
            o.TilePosYX = o.TileInitialPosYX;
            %Below is a safeguard incase wrong positions found - can do
            %this as we knwo what the answer should be.
            MaxY = max(o.TileInitialPosYX(:,1));
            MaxX = max(o.TileInitialPosYX(:,2));
            if MaxY*MaxX ~= nSeries
                warning('Number of tiles (%d) is not equal to maximum Y position (%d) multiplied by maximum X position (%d)'...
                    , nSeries, MaxY, MaxX)
                break
            else
                TilePosY = flip(repelem(1:MaxY,MaxX));
                o.TilePosYX(:,1) = TilePosY;
                TilePosX = repmat([flip(1:MaxX),1:MaxX],1,ceil(MaxY/2));
                o.TilePosYX(1:nSeries,2) = TilePosX(1:nSeries);                
            end
        end

        % set up filename grid for this round
        fName = cell(nSerieswPos,1);
        
        %Set top hat structuring elements
        if strcmpi(o.ExtractR,'auto')
            o.ExtractR = round(1/pixelsize);
        end
        if strcmpi(o.DapiR,'auto')
            o.DapiR = round(8/pixelsize);
        end
        SE = strel('disk', o.ExtractR);
        DapiSE = strel('disk', o.DapiR);
        
        %parfor t = 1:nSerieswPos  
        for t = 1:nSerieswPos  
           
            fName{t} = fullfile(o.TileDirectory, ...
                    [o.FileBase{r}, '_t', num2str(t), '.tif']);
            
            if exist(fName{t}, 'file')
                fprintf('Round %d tile %d already done.\n', r, t);
                continue;
            end                   
            
                
            % a new reader per worker
            bfreader = javaObject('loci.formats.Memoizer', bfGetReader(), 0);
            % use the memo file cached before
            bfreader.setId(imfile);

            bfreader.setSeries(scene*t-1);

            for c = 1:nChannels
                % read z stacks
                I = cell(nZstacks,1);
                for z = 1:nZstacks
                    iPlane = bfreader.getIndex(z-1, c-1, 0)+1;
                    I{z} = bfGetPlane(bfreader, iPlane);
                end

                % focus stacking
                IFS = o.fstack_modified(I);

                % tophat
                if c == o.DapiChannel && r == o.ReferenceRound
                    IFS = imtophat(IFS, DapiSE);
                else
                    IFS = imtophat(IFS, SE);
                end

                % write stack image
                imwrite(IFS,...
                    fullfile(o.TileDirectory,...
                    [o.FileBase{r}, '_t', num2str(t), '.tif']),...
                    'tiff', 'writemode', 'append');
            end
            fprintf('Round %d tile %d finished.\n', r, t);
            bfreader.close();

        end        

        for t=1:nSerieswPos
            o.TileFiles{r,o.TilePosYX(t,1), o.TilePosYX(t,2)} = fName{t};
        end
    end
    
    o.EmptyTiles = cellfun(@isempty, squeeze(o.TileFiles(o.ReferenceRound,:,:)));

end
