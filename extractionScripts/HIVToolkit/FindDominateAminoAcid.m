function PosSet = FindDominateAminoAcid( ProbMatrix,OccupyNum,Bound,Cutoff1,Cutoff2 )
%% find the main amino acid of each column according to the probability.
%  AASymbol='ARNDCQEGHILKMFPSTWYVBZX*-?';
if nargin == 2
   Cutoff2 = 0.1;
end
[Row,Col] = size( ProbMatrix );
if Row > 24 % we don't consider - and ? as dominate amino acids.
   ProbMatrix = ProbMatrix(1:24,:);
end
PosSet = cell(1,Col);
for p = 1:Col
    [ MaxValue,MaxPos ] = sort( ProbMatrix(:,p),'descend');
    if nargin > 2 && OccupyNum(p) > Bound
       if MaxValue(2) > Cutoff1
          PosSet{ p } = MaxPos(1:2);
       else
          PosSet{ p } = MaxPos(1);
       end
    else
       if MaxValue(2) > Cutoff2
          PosSet{ p } = MaxPos(1:2);
       else
          PosSet{ p } = MaxPos(1);
       end
    end 
end
end