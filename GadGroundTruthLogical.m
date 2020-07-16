function ToUse = GadGroundTruthLogical(o,Method)
%% ToUse = GadGroundTruthLogical()
% Returns a logical array that is true for all spots that should be Gad
% i.e. that have high intensity in Gad round, Gad channel and not in any
% other channels of Gad round

%%
% Has to be closer to Line2 (y=0) than Line1 (y=x), where x is Gad colour
% channel intensity to be correct.
Line1 = [0,0;1,1];
Line2 = [0,0;1,0];

if strcmpi('Prob',Method)
    GadColor = o.pGadColor;
elseif strcmpi('Pixel',Method)
    GadColor = o.pxGadColor;
end
GadAnchorChannel = 7;
NonGadChannels = setdiff(1:o.nBP,[o.GadChannel,GadAnchorChannel]);
nChannels = length(NonGadChannels);
ToUseAll = false(length(GadColor),nChannels);
for b=1:nChannels
    pts = GadColor(:,[o.GadChannel,NonGadChannels(b)]);
    dist1 = point_to_line(pts,Line1(1,:),Line1(2,:));
    dist2 = point_to_line(pts,Line2(1,:),Line2(2,:));
    ToUseAll(:,b) = dist2<dist1;
end

ToUse = all(ToUseAll,2);

end


function d = point_to_line(pt, v1, v2)
% pt should be nx2 
% v1 and v2 are vertices on the line (each 1x2)
% d is a nx1 vector with the orthogonal distances
v1 = [v1,0];
v2 = [v2,0];
pt = [pt,zeros(size(pt,1),1)];
v1 = repmat(v1,size(pt,1),1);
v2 = repmat(v2,size(pt,1),1);
a = v1 - v2;
b = pt - v2;
d = sqrt(sum(cross(a,b,2).^2,2)) ./ sqrt(sum(a.^2,2));
end