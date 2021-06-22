function o = extract_and_filter(o)

%% So if no Parallel Computing Toolbox, fails straight away
GPU_test = gpuArray([1]);

%% Logging
if o.LogToFile
    diary(o.LogFile);
    cleanup = onCleanup(@()diary('off'));
end
%%
% create tiff files for each tile that are top-hat filtered versions of
% original czi files
% This uses a GPU for speed
    o.TileFiles = cell(o.nRounds+o.nExtraRounds,1,1,1); % 1,1,1 because we don't yet know how many tiles
    
    %New filter
    h = -hanning(o.ExtractR2*2+1);
    h = -h/sum(h);
    h(o.ExtractR2+1-o.ExtractR1:o.ExtractR2+1+o.ExtractR1) = ...
        h(o.ExtractR2+1-o.ExtractR1:o.ExtractR2+1+o.ExtractR1)+hanning(o.ExtractR1*2+1)/sum(hanning(o.ExtractR1*2+1));
    SE = ftrans2(h');
    SE = single(gpuArray(SE));
%     h2D = ftrans2(h');
%     hzdirection = hanning(3);
%     hzdirection = reshape(hzdirection,[1,1,3]);
%     SE = h2D.*hzdirection;

    AnchorChannelsToUse = [o.DapiChannel,o.AnchorChannel];
    UsedEmptyTiles = false;   %If running for only a few tiles, this will change to true.
    
    
    for r = 1:o.nRounds+o.nExtraRounds
        imfile = fullfile(o.InputDirectory, [o.FileBase{r}, o.RawFileExtension]);

        % construct a Bio-Formats reader with the Memoizer wrapper
        bfreader = loci.formats.Memoizer(bfGetReader(), 0);
        
        if exist(imfile)>0
            % initiate reader
            bfreader.setId(imfile);
        else
            %Wait for file to exist
            if r==1
                MaxTime = o.MaxWaitTime1;   %Don't wait long if first round
            else
                MaxTime = o.MaxWaitTime;   %Max amount of time to wait in seconds
            end
            count = 0;
            while exist(imfile)==0
                pause(1);
                count = count+1;
                if count >= MaxTime
                    error(sprintf(strcat('No file named:\n  ',imfile,'\ncreated in allowed time')));
                end
            end
            %Wait for file to stop loading
            OldBytes = 0;
            NewBytes = 0.00001;
            while NewBytes>OldBytes
                pause(5);
                fileinfo = dir(imfile);
                OldBytes = NewBytes;
                NewBytes = fileinfo.bytes;
            end
            % initiate reader
            bfreader.setId(imfile);
        end

        % get some basic image metadata
        [nSeries, nSerieswPos, nChannels, o.nZ, xypos, o.XYpixelsize,o.Zpixelsize] = ...
            get_ome_tilepos(bfreader);
        if isempty(xypos) || size(xypos, 1)==1
            if r == 1
                warning('first round xypos empty - using values from initial manual input')
                assert(~isempty(o.TileInitialPosYX), 'xypos unavailable')
                xypos = o.TileInitialPosYX;
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
            if isempty(o.AutoThresh)
                o.AutoThresh = zeros(nSerieswPos,nChannels,o.nRounds+o.nExtraRounds);  
            end
            if isempty(o.HistCounts)
                o.HistValues = -o.TilePixelValueShift:1:2^16-o.TilePixelValueShift;  %Entire range of possible pixel values
                o.HistCounts = zeros(length(o.HistValues),nChannels,o.nRounds);
            end
            
            % find x and y grid spacing as median of distances that are about
            % right
            dx = xypos(:,1)-xypos(:,1)'; % all pairs of x distances
            xStep = median(dx(abs(1- dx(:)/o.MicroscopeStepSize)<.5));
            dy = xypos(:,2)-xypos(:,2)'; % all pairs of y distances
            yStep = median(dy(abs(1- dy(:)/o.MicroscopeStepSize)<.5));
        
        
            % find coordinates for each tile
            if isempty(o.TileInitialPosYX)
                if nSeries==1
                    o.TileInitialPosYX = [1,1];
                else
                    o.TileInitialPosYX = fliplr(1+round((xypos - min(xypos))./[xStep yStep]));
                end
            end
            TilePosYX = o.TileInitialPosYX;
            %Below is a safeguard incase wrong positions found - can do
            %this as we know what the answer should be.
            MaxY = max(TilePosYX(:,1));
            MaxX = max(TilePosYX(:,2));
            %Sometimes get Nan, if only one Nan, then check if all tiles
            %arranged along only one direction i.e. Nan should be 1.
            if max(isnan([MaxX,MaxY])) && nanmax(MaxX,MaxY)==nSeries
                if isnan(MaxX); MaxX=1; else; MaxY=1; end
            end
            if MaxY*MaxX ~= nSeries
                warning('Number of tiles (%d) is not equal to maximum Y position (%d) multiplied by maximum X position (%d)'...
                    , nSeries, MaxY, MaxX)
                break
            else
                TilePosY = flip(repelem(1:MaxY,MaxX));
                TilePosYX(:,1) = TilePosY;
                TilePosX = repmat([flip(1:MaxX),1:MaxX],1,ceil(MaxY/2));
                TilePosYX(1:nSeries,2) = TilePosX(1:nSeries);
            end
        end
        
        %Tile index in nd2 file different to index in o.EmptyTiles
        t_save_value = sub2ind([MaxY,MaxX],TilePosYX(:,1),TilePosYX(:,2));
        
        o.TilePosYXC = zeros(nSerieswPos*nChannels,3);

        % set up filename grid for this round
        fName = cell(nSerieswPos*nChannels,1);
        
        Index = 1;
        %parfor t = 1:nSerieswPos  
        for t_index = 1:nSerieswPos  
            t = t_save_value(t_index);
                       
            % a new reader per worker
            bfreader = javaObject('loci.formats.Memoizer', bfGetReader(), 0);
            % use the memo file cached before
            bfreader.setId(imfile);

            bfreader.setSeries(scene*t_index-1);
            for c = 1:nChannels
                tic
                fName{Index} = fullfile(o.TileDirectory, ...
                    [o.FileBase{r}, '_t', num2str(t),'c', num2str(c), '.tif']);  
                
                if exist(fName{Index}, 'file')
                    fprintf('Round %d, tile %d, channel %d already done.\n', r, t, c);
                    o.TilePosYXC(Index,:) = [TilePosYX(t_index,:),c];          %Think first Z plane is the highest
                    o.TileFiles{r,o.TilePosYXC(Index,1), o.TilePosYXC(Index,2),o.TilePosYXC(Index,3)} = fName{Index};
                    if o.AutoThresh(t,c,r) == 0
                        if c == o.DapiChannel && r == o.ReferenceRound; continue; end
                        IFS = o.load_3D_GPU(r,o.TilePosYXC(Index,1),o.TilePosYXC(Index,2),c)-o.TilePixelValueShift;
                        o.AutoThresh(t,c,r) = gather(median(abs(IFS(:)))*o.AutoThreshMultiplier);
                        if r~= o.ReferenceRound
                            o.HistCounts(:,c,r) = o.HistCounts(:,c,r)+gather(histc(IFS(:),o.HistValues));
                        end
                    end
                    Index = Index+1;
                    continue;
                elseif r == o.ReferenceRound && ~ismember(c,AnchorChannelsToUse)
                    %Only need anchor and dapi tiles in reference round
                    continue;
                elseif min(size(o.EmptyTiles))==1 && ~ismember(t,o.EmptyTiles)
                    %If specify o.EmptyTiles, only run for tiles in o.EmptyTiles
                    o.TilePosYXC(Index,:) = [TilePosYX(t_index,:),c];          %Think first Z plane is the highest
                    o.TileFiles{r,o.TilePosYXC(Index,1), o.TilePosYXC(Index,2),o.TilePosYXC(Index,3)} = fName{Index};
                    UsedEmptyTiles = true;
                    EmptyTilesOrig = o.EmptyTiles;
                    Index = Index+1;
                    continue;
                end
                                                                        
                %TopHat SE
%                 if c == o.DapiChannel && r == o.ReferenceRound    
%                         %SE = strel3D_2(20,10);       % I.e. set to 8 microns for DAPI
%                         SE = get_3DSE(o.DapiR1YX,o.DapiR1Z,o.DapiR2YX,o.DapiR2Z);       
%                 else
%                         %SE = strel3D_2(3,3);    %I.e. Set to 1 micron
%                         SE = get_3DSE(o.ExtractR1YX,o.ExtractR1Z,o.ExtractR2YX,o.ExtractR2Z);
%                 end

                I = zeros(o.TileSz,o.TileSz,o.nZ,'gpuArray'); 
                for z = 1:o.nZ
                    iPlane = bfreader.getIndex(z-1, c-1, 0)+1;
                    I(:,:,z) = bfGetPlane(bfreader, iPlane);
                end       
                
                %Make noise white first by divding amplitude of FT
                
                %FT = fftn(I);
                %Norm_FT = FT ./ abs(FT);
                %filter = fspecial3('gaussian',size(I),2);%DESCRIBE BETTER!!!!
                %Shiftfilter = fftshift(filter);     %Shift for FT so centered on 0
                %FT_filter = fftn(Shiftfilter);          
                %NormFT_filter = FT_filter ./ abs(FT);
                %Final_FT = Norm_FT .* NormFT_filter;
                %IFS = ifftn(Final_FT);
                
                %I = ifftn(Norm_FT);
                I = single(padarray(I,(size(SE)-1)/2,'replicate','both'));
                IFS = convn(I,SE,'valid'); 
                clearvars I  %Free up GPU memory
                
                %Scaling so fills uint16 range.
                if c == o.DapiChannel && r == o.ReferenceRound  
                    if strcmpi(o.DapiScale, 'auto')
                        o.DapiScale = 10000/max(IFS(:));
                    end
                    IFS = IFS*o.DapiScale;
                else
                    %Finds o.ExtractScale from first image and uses this
                    %value for the rest
                    if strcmpi(o.ExtractScale, 'auto')
                        o.ExtractScale = 30000/max(IFS(:));
                    end
                    IFS = IFS*o.ExtractScale;
                    
                    %Determine auto thresholds
                    o.AutoThresh(t,c,r) = gather(median(abs(IFS(:)))*o.AutoThreshMultiplier);
                end
                
                if r ~= o.ReferenceRound  
                    %Get histogram data
                    IFS = int32(IFS);
                    %AbridgedBaseIm = IFS(:,:,8:15);
                    o.HistCounts(:,c,r) = o.HistCounts(:,c,r)+gather(histc(IFS(:),o.HistValues));
                end
                
                
                %Append each z plane to same tiff image
                %Add o.TilePixelValueShift so keep negative pixels for background analysis
                IFS = gather(uint16(IFS+o.TilePixelValueShift));   
                %IFS = uint16(IFS + o.TilePixelValueShift);
                for z = 1:o.nZ
                    imwrite(IFS(:,:,z),... 
                            fullfile(o.TileDirectory,...
                            [o.FileBase{r}, '_t', num2str(t),'c', num2str(c), '.tif']),...
                            'tiff', 'writemode', 'append');
                end
                

                o.TilePosYXC(Index,:) = [TilePosYX(t_index,:),c];          %Think first Z plane is the highest
                o.TileFiles{r,o.TilePosYXC(Index,1), o.TilePosYXC(Index,2),o.TilePosYXC(Index,3)} = fName{Index};
                fprintf('Round %d tile %d colour channel %d finished.\n', r, t, c);                                               
                Index = Index+1; 
                toc
            end
            bfreader.close();
            
        end
        

    end
    
    o.EmptyTiles = cellfun(@isempty, squeeze(o.TileFiles(o.ReferenceRound,:,:,1)))*0;
    if UsedEmptyTiles
        o.EmptyTiles(:) = true;
        o.EmptyTiles(EmptyTilesOrig) = false;
    end
    
    %Plot boxplots showing distribution af AutoThresholds
    if o.Graphics
        Thresholds = [];
        group = [];
        index = 1;
        for c=1:nChannels
            for r=1:o.nRounds
                Thresholds = [Thresholds;o.AutoThresh(:,c,r)];
                group = [group;index*ones(size(o.AutoThresh(:,1,1)))];
                index = index+1;
            end
        end
        %Add anchor
        AnchorChannels = AnchorChannelsToUse(AnchorChannelsToUse~=o.DapiChannel);
        AnchorLabels = cell(size(AnchorChannels));
        LabelIndex = 1;
        for c=AnchorChannels
            Thresholds = [Thresholds;o.AutoThresh(:,c,o.ReferenceRound)];
            group = [group;index*ones(size(o.AutoThresh(:,1,1)))];
            if c == o.AnchorChannel
                AnchorLabels(LabelIndex) = {'Full Anchor'};
            else
                AnchorLabels(LabelIndex) = {['Split Anchor ',num2str(LabelIndex-1)]};
            end
            LabelIndex = LabelIndex+1;
            index=index+1;
        end
        
        figure(43290);
        colors = colormap(lines(nChannels));
        Colors = repelem(colors,o.nRounds,1);
        Colors = [Colors;repelem([0,0,0],nChannels,1)];
        Labels = [string(repmat(1:o.nRounds,1,nChannels)),string(AnchorLabels)];
        boxplot(Thresholds,group,'Colors',Colors, 'plotstyle', 'compact','labels', Labels);
        set(gca,'TickLength',[0 0]);
        ylabel('AutoThreshold');
        xlabel('Round');
        hold on
        for c=1:nChannels
            plot(NaN,1,'color', colors(c,:), 'LineWidth', 4);       %For legend labels
        end
        leg = legend(o.bpLabels,'Location','northwest');
        title(leg,'Color Channel');
        hold off
    end
    
    %Plot histograms to make sure they are smooth
    nPixels = sum(o.HistCounts(:,1,1));
    if o.Graphics       
        figure(43291);
        index = 1;
        for r=1:o.nRounds
            for b=1:nChannels
                subplot(o.nRounds,nChannels,index)
                histogram('BinEdges',[o.HistValues-0.5,max(o.HistValues)+0.5],'BinCounts',o.HistCounts(:,b,r)/nPixels,'DisplayStyle','stairs');
                xlim([-1000,1000]);
                ylim([0,max(o.HistCounts(:,b,r))/nPixels]);
                if b==4
                    title(strcat('Round ',num2str(r)));
                end
                index = index+1;
            end
        end
    end
    
end
