function [ TransferedMatrix, ProbMatrix, FullSeqRegion ] = Nucleotide2Int( CurrentSeq, Threshold )
[ Row,Col ] = size( CurrentSeq );
if nargin < 2
   Threshold = ceil(Row*0.2);
end
FullSeqRegion = ones(1,Col);
TransferedMatrix = 16*ones(Row,Col);
for p = 1:Row
    for q = Col:-1:1 % A T C G N R W Y M K S H B V D -
        if CurrentSeq(p,q) == 'A' || CurrentSeq(p,q) == 'a'
           TransferedMatrix(p,q) = 1;
        elseif CurrentSeq(p,q) == 'T' || CurrentSeq(p,q) == 't' 
           TransferedMatrix(p,q) = 2;
        elseif CurrentSeq(p,q) == 'C' || CurrentSeq(p,q) == 'c'
           TransferedMatrix(p,q) = 3;
        elseif CurrentSeq(p,q) == 'G' || CurrentSeq(p,q) == 'g'
           TransferedMatrix(p,q) = 4;
        elseif CurrentSeq(p,q) == 'N' || CurrentSeq(p,q) == 'n'
           TransferedMatrix(p,q) = 5;
        elseif CurrentSeq(p,q) == 'R' || CurrentSeq(p,q) == 'r'
           TransferedMatrix(p,q) = 6;
        elseif CurrentSeq(p,q) == 'W' || CurrentSeq(p,q) == 'w'
           TransferedMatrix(p,q) = 7;
        elseif CurrentSeq(p,q) == 'Y' || CurrentSeq(p,q) == 'y'
           TransferedMatrix(p,q) = 8;
        elseif CurrentSeq(p,q) == 'M' || CurrentSeq(p,q) == 'm'
           TransferedMatrix(p,q) = 9;
        elseif CurrentSeq(p,q) == 'K' || CurrentSeq(p,q) == 'k'
           TransferedMatrix(p,q) = 10;
        elseif CurrentSeq(p,q) == 'S' || CurrentSeq(p,q) == 's'
           TransferedMatrix(p,q) = 11;
        elseif CurrentSeq(p,q) == 'H' || CurrentSeq(p,q) == 'h'
           TransferedMatrix(p,q) = 12;
        elseif CurrentSeq(p,q) == 'B' || CurrentSeq(p,q) == 'b'
           TransferedMatrix(p,q) = 13;           
        elseif CurrentSeq(p,q) == 'V' || CurrentSeq(p,q) == 'v'
           TransferedMatrix(p,q) = 14;
        elseif CurrentSeq(p,q) == 'D' || CurrentSeq(p,q) == 'd'
           TransferedMatrix(p,q) = 15;           
        end
    end
end

ProbMatrix = zeros(16,Col);
for p = 1:Col
    for q = 1:16
        ProbMatrix(q,p) = length( find( TransferedMatrix(:,p) == q ) );
    end
    if length( find( TransferedMatrix(:,p) == 16 ) ) > Threshold
       FullSeqRegion( p ) = 0;
    end
end
ProbMatrix = ProbMatrix/Row;


end

