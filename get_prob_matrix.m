function ProbMatrix = get_prob_matrix(o,SpotCode,GeneNo)
%Given the best gene, this finds the probability matrix for the match i.e.
%the breakdown of log(Prob) into the different rounds and channels

%x = min(o.cSpotColors(:))-1:max(o.cSpotColors(:))-1;    %subsitiution x=lambda*g, -1 due to matlab indexing
HistZeroIndex = find(o.SymmHistValues == 0);            %As HistProbs is of different length to x
ProbMatrix = zeros(7,7);
for b=1:7
    for r=1:7
        f = SpotCode(b,r);
        LogConvDist = log(conv(o.LambdaDist(:,GeneNo,b,r),o.HistProbs(:,b,r),'same'));
        %Need to relate it to something that is the same for all genes i.e.
        %probability spot explained by background alone
        ProbMatrix(b,r) = LogConvDist(o.ZeroIndex-1+f)-log(o.HistProbs(HistZeroIndex+f,b,r));
        
        %x2 = x(x<HistZeroIndex+f);      %So ensure indices>0
        %hIndices = HistZeroIndex+f-x2;
        %Use = hIndices<length(o.SymmHistValues);
        %HistDist = o.HistProbs(hIndices(Use),b,r);
        %LambdaIndices = find(x<HistZeroIndex+f);
        %ProbMatrix(b,r) = log(sum(HistDist.*o.LambdaDist(LambdaIndices(Use),GeneNo,b,r)));
    end
end