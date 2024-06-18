function OrigSeqPop = MoveSeqRight2Left( OrigSeqPop,Lower,Upper,SelectRow )
% this function moves all the sequence in OrigSeqPop from right (Upper) to
% left side (Lower)

[Row,Col] = size( OrigSeqPop );

if nargin == 1
   Lower = 1; Upper = Col; SelectRow = 1:Row;
elseif nargin == 3
   SelectRow = 1:Row;
end

Len = Upper-Lower+1;
SelectRowLen = length( SelectRow );

AdjustSeq = cell( 1,SelectRowLen );
parfor p = 1:SelectRowLen
    AdjustSeq{ p } = OrigSeqPop( SelectRow(p),Lower:Upper );   
end

parfor n = 1: SelectRowLen % move sequence from right to left
    % All = OrigSeqPop(q,:);
    % local = OrigSeqPop(q,Lower:Upper);
    NonGapPos = find( OrigSeqPop( SelectRow(n),Lower:Upper ) ~='-' );
    if isempty( NonGapPos )==0
       NonGapPosLen = length( NonGapPos );
       Temp = AdjustSeq{ n }( NonGapPos );
       AdjustSeq{ n } = repmat('-',1,Len);
       AdjustSeq{ n }( 1:NonGapPosLen ) = Temp;
    end
end 
for p = 1:SelectRowLen
    OrigSeqPop( SelectRow(p),Lower:Upper ) = AdjustSeq{ p };
end
end