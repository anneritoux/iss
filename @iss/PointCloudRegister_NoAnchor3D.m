function o = PointCloudRegister_NoAnchor3D(o, y, x0, A0, nTiles, Options)     %MADE A THE SAME FOR ALL TILES
% o = o.PointCloudRegister(y, x, A0, Options)
% 
% Perform point cloud registration to map points x onto points y by
% iterative closest point: repeatedly finding the best y for each x, 
% and doing linear regression to find the M that maps best maps x to y
%
% inputs:
% y is a cell containig the centered YX location of all spots in all rounds 
% and colour channels for all tiles
%
% x0 is a cell containing the non centered YX location of spots in the 
% anchor channel for all tiles
%
% A0 are the initial scaling matrices for each colour channel 
% taking account of chromatic aberration. All default to identity if not
% specified
%
% ToPlot: array of form [t,b,r] of specific example case to show plot of
% for debugging purposes
%
% Options: what type of fit. For now ignored, the only option is a general linear
% model where x gets an extra column of ones and M is 2x3.
%%
MaxIter = 100;
nD = 3;

%centre anchor channel spots
x = cell(nTiles,1);
for t=1:nTiles
    if o.EmptyTiles(t); continue; end
    %Center and scale z direction.
    x(t) = {(x0{t} - o.CentreCorrection).*[1,1,o.Zpixelsize/o.XYpixelsize]};
end


if nargin<4 || isempty(A0)
    A0 = zeros(nD,nD,o.nBP);
    for b=1:o.nBP
        A0(:,:,b) = eye(nD);
    end
elseif max(size(A0))==nD
    A = zeros(nD,nD,o.nBP);
    for b=1:o.nBP
        A(:,:,b) = A0;
    end
    A0 = A;    
end

if isempty(o.PcDist)
    o.PcDist = inf;
end

%Initialize variables
A = A0;
D = o.D0;

fprintf('\nPCR - Finding well isolated points');
% find well isolated points as those whose second neighbor is far
for t=1:nTiles
    if o.EmptyTiles(t); continue; end
    for r=o.UseRounds
        for b=o.UseChannels
            
            % make kd tree - default options!
            k0 = KDTreeSearcher(y{t,b,r});
            [~, d2] = k0.knnsearch(y{t,b,r}, 'k', 2);
            if isfinite(o.PcDist) && size(y{t,b,r},1) > 1 
                y(t,b,r) = {y{t,b,r}(d2(:,2)>o.PcDist*2,:)};
            end
            
        end
    end
end

fprintf('\nPCR - Making kd trees');
%Make kd trees out of these well isolated points
k = cell(nTiles,o.nBP,o.nRounds);
for t=1:nTiles
    if o.EmptyTiles(t); continue; end
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

if isempty(o.ToPlot)
    fprintf('\nPCR - Iteration   ');
end
for i=1:MaxIter
    if isempty(o.ToPlot)
        if i<10
            fprintf('\b%d', i);
        else
            fprintf('\b\b%d', i);
        end    
        if i ==MaxIter
            fprintf('\nPCR - Max number of iterations reached');
        end
    end
    
    LastNeighbor = Neighbor;
    
    for t=1:nTiles
        if o.EmptyTiles(t); continue; end
        for r=o.UseRounds 
            for b=o.UseChannels
                xM(t,b,r) = {(A(:,:,b)*(x{t} + D(t,:,r))')'};   
            end
        end
    end
        
    %This part finds new neighbours and new estimates for A
    for b=o.UseChannels
        xA = [];
        yA = [];
        for t=1:nTiles
            if o.EmptyTiles(t); continue; end
            for r=o.UseRounds          
                Neighbor(t,b,r) = {k{t,b,r}.knnsearch(xM{t,b,r})};
                [~,Dist] = k{t,b,r}.knnsearch(xM{t,b,r});
                UseMe(t,b,r) = {Dist<o.PcDist};                
                MyNeighb(t,b,r) = {Neighbor{t,b,r}(UseMe{t,b,r}>0)};
                nMatches(t,b,r) = sum(UseMe{t,b,r});
                Error(t,b,r) = sqrt(mean(Dist(UseMe{t,b,r}>0).^2));
            
                xShift = x{t}(UseMe{t,b,r}>0,:) + D(t,:,r);         %Add shift between rounds here
                xA = vertcat(xA, xShift);      
                yA = vertcat(yA, y{t,b,r}(MyNeighb{t,b,r},:));                      
            
            end
        end
        A(:,:,b) = xA\yA;
        %A_2D = xA(:,1:2)\yA(:,1:2);
        %A(1:2,1:2,b) = A_2D;
    end
    
    %This part finds new estimates for D
    for t=1:nTiles
        if o.EmptyTiles(t); continue; end
        for r=o.UseRounds
            xD = [];
            yD = [];
            for b=o.UseChannels
                xD = vertcat(xD,x{t}(UseMe{t,b,r}>0,:));
                yScaled = (A(:,:,b)\y{t,b,r}(MyNeighb{t,b,r},:)')';
                yD = vertcat(yD, yScaled);
            end
            D(t,:,r) = mean(yD - xD);
        end
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
        InZPlaney = y{t,b,r}(:,3)==(o.ZPlane-1-(o.nZ-1)/2)*o.Zpixelsize/o.XYpixelsize;
        plot(y{t,b,r}(InZPlaney,2), y{t,b,r}(InZPlaney,1), 'g+');
        InZPlanex = xM{t,b,r}(:,3)==(o.ZPlane-1-(o.nZ-1)/2)*o.Zpixelsize/o.XYpixelsize;
        Neighboury = y{t,b,r}(MyNeighb{t,b,r},:);
        InZPlaneNeighboury = Neighboury(:,3)==(o.ZPlane-1-(o.nZ-1)/2)*o.Zpixelsize/o.XYpixelsize;
        plot(xM{t,b,r}(InZPlanex,2), xM{t,b,r}(InZPlanex,1), 'r+');
        plot([xM{t,b,r}(UseMe{t,b,r}>0&InZPlanex,2) Neighboury(InZPlaneNeighboury,2)]',...
            [xM{t,b,r}(UseMe{t,b,r}>0&InZPlanex,1) Neighboury(InZPlaneNeighboury,1)]', 'k-', 'linewidth', 1);

        drawnow;
    end
    
    if min(min(min(cellfun(@isequal, Neighbor, LastNeighbor)))) == 1; break; end
    
end
fprintf('\n')
%%
o.A = A;
o.D = D.*[1,1,o.XYpixelsize/o.Zpixelsize];
o.nMatches = nMatches;
o.Error = Error;

