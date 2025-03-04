function OrigSeqPop = MoveSeqLeft2Right( OrigSeqPop,Lower,Upper,SelectRow )
% this function moves all the sequence in OrigSeqPop from left side (Lower) to
%  right side (Upper)

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
       AdjustSeq{ n }( Len-NonGapPosLen+1:Len ) = Temp;
    end
end
    
for p = 1:SelectRowLen
    OrigSeqPop( SelectRow(p),Lower:Upper ) = AdjustSeq{ p };
end
end