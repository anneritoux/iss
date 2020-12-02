%% read in pixel colors for the whole image
OutputDirectory = '/Volumes/G-DRIVE mobile USB 1/UCL/ISS/B9S4_Slice001-002/Slice001output/ByPixelAnchor/NewTransform-PCR6/OMP/Without1stZPlane/';
[nY, nX] = size(o.EmptyTiles);
nTiles = nY*nX;
NonemptyTiles = find(~o.EmptyTiles)';
%Spots on achor round cover whole range of coordinates, same for each tile
AnchorLocalYX = zeros(o.TileSz^2,2);
AnchorLocalYX(:,1) = repelem(1:o.TileSz,1,o.TileSz);
AnchorLocalYX(:,2) = repmat(1:o.TileSz,1,o.TileSz);
nPixels = size(AnchorLocalYX,1);
AnchorLocalYX = [AnchorLocalYX,ones(nPixels,1)];

%for t=1:length(NonemptyTiles) 
t=1;
tile_no = NonemptyTiles(t);
%Get pixel colors
[GoodAnchorLocalYX,GoodSpotColors] = o.get_spot_colors(tile_no,AnchorLocalYX,nPixels);
load([OutputDirectory,'tile',num2str(tile_no),'_all_pixel_coefs.mat'], 'coefs');
ScaledColors = (double(GoodSpotColors)-o.SHIFT)./o.SCALE;
%Only use pixels for which all gene coefficients were zero.
nZeroCoefs = sum(coefs(:,1:73)==0,2);
Use = nZeroCoefs==73;   
Tile1ZeroCoef = ScaledColors(Use,:);
%end

%Get eigenvalues
CovMatrix = cov([Tile1ZeroCoef;Tile2ZeroCoef;Tile3ZeroCoef;Tile4ZeroCoef;...
    Tile5ZeroCoef;Tile6ZeroCoef;Tile7ZeroCoef;Tile8ZeroCoef]);
[V,D] = eig(CovMatrix,'vector');
[D, ind] = sort(D);
V = V(:, ind);

%visualise eigenvalue
figure; imagesc(reshape(V(:,end),[7,7]));
colormap(gca,bluewhitered);
title(sprintf('Eigenvalue = %.2f', D(end)));
xlabel('Round');
set(gca, 'ytick', 1:o.nBP);
set(gca, 'YTickLabel', o.bpLabels);
ylabel('Color Channel');


%% Do OMP for all gene bled codes to see how well eigenvectors can explain them.
% Using Eigenvectors of covariance matrix
nCodes = length(o.CharCodes);
GeneBledCodes = o.ScaledBledCodes(1:nCodes,:);
GeneBledCodes = GeneBledCodes./vecnorm(GeneBledCodes,2,2);
nEigVector = 8;
GeneCodeEigResidue = zeros(nCodes,nEigVector);
GeneCodeEigCoefs = zeros(nCodes,nEigVector);
for i=1:nEigVector
    for g=1:nCodes
        [coef,~,GeneCodeEigResidue(g,i)] = omp_specify_atoms(V, GeneBledCodes(g,:)',50-i:49);
        if i== nEigVector
            GeneCodeEigCoefs(g,:) = coef(50-i:49)';
        end
    end
end

% Using colour channel strips
BackgroundVectors = o.ScaledBledCodes(nCodes+1:end,:);
BackgroundVectors = BackgroundVectors./vecnorm(BackgroundVectors,2,2);
BackgroundVectors = BackgroundVectors';
nBackground = 7;
GeneCodeBckgrndResidue = zeros(nCodes,nBackground);
GeneCodeBckgrndCoefs = zeros(nCodes,nBackground);
for i=1:nBackground
    for g=1:nCodes
        [coef,~,GeneCodeBckgrndResidue(g,i)] = omp_specify_atoms(BackgroundVectors, GeneBledCodes(g,:)',1:i);
        if i== nBackground
            GeneCodeBckgrndCoefs(g,:) = coef(1:i)';
        end
    end
end

%% Do OMP for some example spots with both old background and new eigenvectors.
%  See how many of each are within 3 best atoms for pixels.

AtomDict = o.ScaledBledCodes(:,:);
AtomDict(81:84,:) = V(:,end-3:end)';
AtomDict = AtomDict./vecnorm(AtomDict,2,2);

SpotColors = (double(o.pxSpotColors)-o.SHIFT)./o.SCALE;
SpotColors = SpotColors(:,:);
nSpots = size(SpotColors,1);
nAtoms = 3;     %number of atoms to select for each spot.

SelectedAtoms = zeros(nSpots,nAtoms);
for s=1:nSpots
    [~,SelectedAtoms(s,:)] = omp(AtomDict',SpotColors(s,:)',nAtoms);
end

%For spots with 1,2,3 background see how atoms distributed
nBackgroundSelected = sum(SelectedAtoms>73,2);
BackgroundCount = zeros(size(AtomDict,1)-73,3,3);
for n=1:3
    SelectedAtomsSet = SelectedAtoms(nBackgroundSelected==n,:);
    for b=1:size(AtomDict,1)-73
        for c=1:3
            BackgroundCount(b,c,n) =  sum(SelectedAtomsSet(:,c)==73+b);
        end
    end
end
        

BackgroundCount = zeros(11,4);
for b=1:11
    for c=1:3
        BackgroundCount(b,c) =  sum(SelectedAtoms(:,c)==73+b);
    end
end
BackgroundCount(:,4) = sum(BackgroundCount(:,1:3),2);


FinalBackgroundCodes = zeros(7,49);
%how to get background vectors with same intensity before z scoring
UnbledBackground = zeros(o.nBP,o.nRounds);
UnbledBackground(1,:) = 1;  %Color Channel 0
BledBackground = o.pBleedMatrix(:,:,1)* UnbledBackground;
BledBackground = BledBackground./squeeze(o.SCALE);
FinalBackgroundCodes(1,:) = BledBackground(:);
UnbledBackground = zeros(o.nBP,o.nRounds);
UnbledBackground(4,:) = 1;  %Color Channel 3
BledBackground = o.pBleedMatrix(:,:,1)* UnbledBackground;
BledBackground = BledBackground./squeeze(o.SCALE);
FinalBackgroundCodes(2,:) = BledBackground(:);
UnbledBackground = zeros(o.nBP,o.nRounds);
UnbledBackground(6,:) = 1;  %Color Channel 5
BledBackground = o.pBleedMatrix(:,:,1)* UnbledBackground;
BledBackground = BledBackground./squeeze(o.SCALE);
FinalBackgroundCodes(3,:) = BledBackground(:);
FinalBackgroundCodes(4,:) = V(:,end)';
FinalBackgroundCodes(5,:) = V(:,end-1)';
FinalBackgroundCodes(6,:) = V(:,end-2)';
FinalBackgroundCodes(7,:) = V(:,end-3)';
FinalBackgroundCodes = FinalBackgroundCodes./vecnorm(FinalBackgroundCodes,2,2);

%% Read in random 100,000 pixels from each tile.
[nY, nX] = size(o.EmptyTiles);
nTiles = nY*nX;
NonemptyTiles = find(~o.EmptyTiles)';
%Spots on achor round cover whole range of coordinates, same for each tile
AnchorLocalYX = zeros(o.TileSz^2,2);
AnchorLocalYX(:,1) = repelem(1:o.TileSz,1,o.TileSz);
AnchorLocalYX(:,2) = repmat(1:o.TileSz,1,o.TileSz);
nPixels = size(AnchorLocalYX,1);
AnchorLocalYX = [AnchorLocalYX,ones(nPixels,1)];

SpotColors = zeros(800000,7,7);
for t=1:nTiles
    [~,GoodSpotColors] = o.get_spot_colors(t,AnchorLocalYX,nPixels);
    SpotIndex = randi([1 length(SpotColors)],100000,1);
    SpotColors((t-1)*100000+1:t*100000,:,:) = GoodSpotColors(SpotIndex,:,:);
end
%ScaledColors = (double(SpotColors)-o.SHIFT)./o.SCALE;
ScaledColors = (double(SpotColorsToTestBackground)-o.SHIFT)./o.SCALE;
DotProduct =  ScaledColors(:,:)*o.SpatialBledCodes(:,:)';
[OrderedAbsValue,GeneOrder]=sort(abs(DotProduct),2,'descend');

%For each method, record number of spots whose best background vector is
%within the top 3.
%Data = cell(5,13);
%Names
Data{1,1} = 'Model';
Data{1,2} = '#Best Background = 1st';
Data{1,3} = 'Median Dot Product';
Data{1,4} = '#Best Background = 2nd';
Data{1,5} = 'Median Dot Product';
Data{1,6} = '#Best Background = 3rd';
Data{1,7} = 'Median Dot Product';
Data{1,8} = '#Best Background in top 3';
Data{1,9} = 'Median Dot Product';
Data{1,10} = 'Background Eigenvectors';
Data{1,11} = 'Background Eigenvalues';
Data{1,12} = 'SHIFT used';
Data{1,13} = 'SCALE used';

%Save data for each model
i=2;
Data{i,1} = 'MAD: ZScore: r and b, Eigenvectors: r and b';
Set1 = GeneOrder(:,1)>73;
Data{i,2} = sum(Set1);
Data{i,3} = median(OrderedAbsValue(Set1,1));
Set2 = GeneOrder(:,2)>73;
Data{i,4} = sum(~Set1&Set2);
Data{i,5} = median(OrderedAbsValue(~Set1&Set2,2));
Set3 = GeneOrder(:,3)>73;
Data{i,6} = sum(~Set1&~Set2&Set3);
Data{i,7} = median(OrderedAbsValue(~Set1&~Set2&Set3,3));
AllWithBackground = [OrderedAbsValue(Set1,1);OrderedAbsValue(~Set1&Set2,2);...
    OrderedAbsValue(~Set1&~Set2&Set3,3)];
Data{i,8} = length(AllWithBackground);
Data{i,9} = median(AllWithBackground);
Data{i,10} = o.SpatialBledCodes(74:end,:,:);
Data{i,11} = EigVal;
Data{i,12} = o.SHIFT;
Data{i,13} = o.SCALE;

%% Check how much background vectors explain gene bled codes
GeneCoefs = o.SpatialBledCodes(1:73,:)*o.SpatialBledCodes(74:end,:)';
ResidueLeft = zeros(73,1);
for GeneNo=1:73
    ExplainedByBackground = sum(o.SpatialBledCodes(74:end,:).*GeneCoefs(GeneNo,:)');
    ResidueLeft(GeneNo) = norm(o.SpatialBledCodes(GeneNo,:)-ExplainedByBackground);
end
Data{i,14} = ResidueLeft;