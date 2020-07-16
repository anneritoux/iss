function [o,x] = PointCloudRegister5(o, y0, x0, A0, nTiles)     %MADE A THE SAME FOR ALL TILES
% o = o.PointCloudRegister5(y, x, A0, Options)
% 
% Perform point cloud registration to map points x onto points y by
% iterative closest point: repeatedly finding the best y for each x, 
% and doing linear regression to find the M that maps best maps x to y
%
% inputs:
% y0 is a cell containig the YX location of all spots in all rounds 
% and colour channels for all tiles
%
% x0{t,b} is a cell containing the YX location of spots in the 
% reference round for tile t, channel b
%
% ToPlot: array of form [t,b,r] of specific example case to show plot of
% for debugging purposes
%
% Output: x is new reference round YX local coordinates. They are different
% from x0 as they are adjusted as PCR proceeds to take account of chromatic
% aberration.
%
% PCR5 tries to provide a regularisation so for a particular tile, for a
% particular round, shifts are similar across all colour channels.
% For a particular tile, for a particular channel, scalings are the same
% across all rounds. 
%%
NonemptyTiles = find(~o.EmptyTiles)';
if isempty(o.TileCentre)
    o.TileCentre = 0.5*[o.TileSz+1,o.TileSz+1];
end

%Colour channels that aren't the RefChannel need adjusting as we go on, to
%account for chromatic aberration.
RefChannelsToAdjust = setdiff(o.ReferenceSpotChannels,o.ReferenceChannel);

%Centre SpotYX
x0(NonemptyTiles,o.ReferenceSpotChannels) = cellfun(@(x0) x0(:,1:2)-o.TileCentre,...
    x0(NonemptyTiles,o.ReferenceSpotChannels),'UniformOutput',false);
x = cell(nTiles,o.nBP);
for t=NonemptyTiles
    for b = o.ReferenceSpotChannels        
        %Append array of ones for translation
        x(t,b) = {[x0{t,b},ones(size(x0{t,b},1),1)]};
    end
end


if nargin<4 || isempty(A0)
    A0 = ones(o.nBP,2);
elseif max(size(A0))==1
    A = zeros(o.nBP,2);
    for b=1:o.nBP
        A(b,1:2) = A0;
    end
    A0 = A;    
end

if isempty(o.PcDist)
    o.PcDist = inf;
end

%Initialize variables
D = zeros(3,2,nTiles,o.nRounds,o.nBP);
for t=NonemptyTiles
    for r = o.UseRounds
        for b=o.UseChannels
            D(1:2,:,t,r,b) = eye(2)*A0(b);
            D(3,:,t,r,b) = o.D0(t,:,r);
        end
    end
end

fprintf('\nPCR - Finding well isolated points');
% find well isolated points as those whose second neighbor is far
y = y0;
for t=NonemptyTiles
    for r=o.UseRounds
        for b=o.UseChannels
            
            % make kd tree - default options!
            k0 = KDTreeSearcher(y0{t,b,r});
            [~, d2] = k0.knnsearch(y0{t,b,r}, 'k', 2);
            if isfinite(o.PcDist) && size(y0{t,b,r},1) > 1 
                y(t,b,r) = {y0{t,b,r}(d2(:,2)>o.PcDist*2,:)};
            end
            
        end
    end
end

fprintf('\nPCR - Making kd trees');
%Make kd trees out of these well isolated points
k = cell(nTiles,o.nBP,o.nRounds);
for t=NonemptyTiles
    for r=o.UseRounds
        for b=o.UseChannels
            k(t,b,r) = {KDTreeSearcher(y{t,b,r})};
        end
    end
end


%%
UseMe = cell(nTiles,o.nBP,o.nRounds);           %nP DIFFERENT FOR DIFFERENT TILES!!!
Neighbor = cell(nTiles,o.nBP,o.nRounds);
MyNeighb = cell(nTiles,o.nBP,o.nRounds);
xM = cell(nTiles,o.nBP,o.nRounds);
nMatches = zeros(nTiles,o.nBP,o.nRounds);
Error = zeros(nTiles,o.nBP,o.nRounds);
TotalNeighbMatches = length(NonemptyTiles)*length(o.UseChannels)*...
    length(o.UseRounds);
nRounds = length(o.UseRounds);
nChannels = length(o.UseChannels);

%C is regularisation paramater
C = zeros(nChannels*nRounds*3);
ScaleShiftRatio = 2000;     %To make getting scale correct equally as important as shift.
%This part is to make chromatic aberration scaling for a single channel
%across all rounds (Value-Mean = 0). Also makes rotation same for same
%colour channel - not sure how to get rid of this part.
for b=1:nChannels
    for r=1:nRounds
        for r2=1:nRounds
            if r==r2
                C(((nRounds*(b-1))+r-1)*3+1,((nRounds*(b-1))+r2-1)*3+1) = ScaleShiftRatio*(1-1/nRounds);
                C(((nRounds*(b-1))+r-1)*3+2,((nRounds*(b-1))+r2-1)*3+2) = ScaleShiftRatio*(1-1/nRounds);
            else
                C(((nRounds*(b-1))+r-1)*3+1,((nRounds*(b-1))+r2-1)*3+1) = ScaleShiftRatio*(-1/nRounds);
                C(((nRounds*(b-1))+r-1)*3+2,((nRounds*(b-1))+r2-1)*3+2) = ScaleShiftRatio*(-1/nRounds);
            end
        end
    end
end

%This part is to make the shift the same for a single round 
%across all channels (Value-Mean = 0).
for r=1:nRounds
    for b=1:nChannels
        for b2=1:nChannels
            if b==b2
                C(((nRounds*(b-1))+r-1)*3+3,((nRounds*(b2-1))+r-1)*3+3) = 1-1/nChannels;
            else
                C(((nRounds*(b-1))+r-1)*3+3,((nRounds*(b2-1))+r-1)*3+3) = -1/nChannels;
            end
        end
    end
end

for i=1:o.PcIter
    
    LastNeighbor = Neighbor;
    

    vertcat(o.RawLocalYX{:,b});
    
    for t=NonemptyTiles
        for b = RefChannelsToAdjust
            %Update position of reference round coordinates, based on new transform D.
            %I.e. (originalcoords-shift)*inv(RotationMatrix)
            x{t,b}(:,1:2) = (x0{t,b}-D(3,:,t,o.ReferenceRound,b))/D(1:2,:,t,o.ReferenceRound,b);
        end
        x_t = vertcat(x{t,:});
        for r=o.UseRounds
            for b=o.UseChannels                                
                xM(t,b,r) = {x_t*D(:,:,t,r,b)+o.TileCentre};   
            end
        end
    end
        
    %This part finds new neighbours and new estimates for A
    for t=NonemptyTiles
        x_t = vertcat(x{t,:});
        x_t_neighb0 = cell(o.nRounds,o.nBP);
        y_t_neighb = cell(o.nRounds,o.nBP);
        for b=o.UseChannels            
            for r=o.UseRounds              
                Neighbor(t,b,r) = {k{t,b,r}.knnsearch(xM{t,b,r})};
                [~,Dist] = k{t,b,r}.knnsearch(xM{t,b,r});
                UseMe(t,b,r) = {Dist<o.PcDist};                
                MyNeighb(t,b,r) = {Neighbor{t,b,r}(UseMe{t,b,r}>0)};
                nMatches(t,b,r) = sum(UseMe{t,b,r});
                Error(t,b,r) = sqrt(mean(Dist(UseMe{t,b,r}>0).^2));   
                x_t_neighb0{r,b} = x_t(UseMe{t,b,r}>0,:);
                y_t_neighb{r,b}=(y{t,b,r}(MyNeighb{t,b,r},:)-o.TileCentre);                                           
            end
        end
        y_t_neighb = cell2mat(y_t_neighb(:));
        nSpots = sum(sum(cell2mat(cellfun(@(x) size(x,1),x_t_neighb0,'uni',false))));
        x_t_neighb = zeros(nSpots,nChannels*nRounds);
        j1=1;
        j2=1;
        for b=o.UseChannels
            for r=o.UseRounds
                x_t_neighb(j1:j1+size(x_t_neighb0{r,b},1)-1,((j2-1)*3)+1:j2*3)=x_t_neighb0{r,b};
                j1=j1+size(x_t_neighb0{r,b},1);
                j2=j2+1;
            end
        end        
        x_t_neighb = [x_t_neighb;C];
        y_t_neighb = [y_t_neighb;zeros(nChannels*nRounds*3,2)];
        D_t = x_t_neighb\y_t_neighb;
        D_t = reshape(D_t',[2,3,7,7]);
        D(:,:,t,:,:) = permute(D_t,[2,1,3,4]);
    end
                
   
    if isempty(o.ToPlot) == 0
        t = o.ToPlot(1);
        b = o.ToPlot(2);
        r = o.ToPlot(3);
        if i == 1
            fprintf('\nPlotting tile %d, color channel %d, round %d', t, b,r);
        end
        figure(29387648);
        fprintf('\nIteration %d: %d matches, mean error %f', i, nMatches(t,b,r), Error(t,b,r));
        clf; hold on
        plot(y{t,b,r}(:,2), y{t,b,r}(:,1), 'g+');
        plot(xM{t,b,r}(:,2), xM{t,b,r}(:,1), 'r+');
        plot([xM{t,b,r}(UseMe{t,b,r}>0,2) y{t,b,r}(MyNeighb{t,b,r},2)]',...
            [xM{t,b,r}(UseMe{t,b,r}>0,1) y{t,b,r}(MyNeighb{t,b,r},1)]', 'k-', 'linewidth', 1);

        drawnow;
    end
    nNeighbMatches = sum(sum(sum(cellfun(@isequal, Neighbor(NonemptyTiles,o.UseChannels,o.UseRounds),...
        LastNeighbor(NonemptyTiles,o.UseChannels,o.UseRounds)))));
    fprintf('\nPCR - Iteration %d: Converged images = %d/%d',i,nNeighbMatches,TotalNeighbMatches);
    if min(min(min(cellfun(@isequal, Neighbor, LastNeighbor)))) == 1; break; end
    
end
fprintf('\n');
if nNeighbMatches<o.PcCovergedImgFrac*TotalNeighbMatches
    warning('\nPCR - Less than %d%% of images have converged',o.PcCovergedImgFrac*100);
end

%%

o.D = D;
o.nMatches = nMatches;
o.Error = Error;
o.PcFailed = nMatches<o.MinSpots;
o.nPcCovergedImg = nNeighbMatches/TotalNeighbMatches;
%Uncentre reference spot YX
x(NonemptyTiles,o.ReferenceSpotChannels) = cellfun(@(x) x(:,1:2)+o.TileCentre,...
    x(NonemptyTiles,o.ReferenceSpotChannels),'UniformOutput',false);

%Record chromatic aberration by taking mean of all transforms for which
%nMatches exceeded threshold
o.A = zeros(o.nBP,2);
CA_ToUse = ~permute(o.PcFailed,[1,3,2]);
for b=o.UseChannels
    CA_ToUse2 = false(size(CA_ToUse));
    CA_ToUse2(:,:,b)=true;
    o.A(b,1) = mean(squeeze(o.D(1,1,CA_ToUse&CA_ToUse2)));
    o.A(b,2) = mean(squeeze(o.D(2,2,CA_ToUse&CA_ToUse2)));
end

