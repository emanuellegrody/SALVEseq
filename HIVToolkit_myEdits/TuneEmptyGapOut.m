function CurrentSeq = TuneEmptyGapOut( CurrentSeq,PosStart,PosEnd )
% move out gaps between the column PosStart and the column PosEnd

Test = 0;
Col = size( CurrentSeq,2 );

if nargin == 1
   PosStart = 1; PosEnd = Col;
end
SelectCol = ones(1,Col);
for p = PosStart:PosEnd
    if isempty( find(CurrentSeq(:,p)~='-',1) ) == 1
       SelectCol( p ) = 0;
    end
end
Start = find( SelectCol==1,1 );
p = find( SelectCol==0,1 );
if isempty(p) == 0
  while p <= PosEnd
    End = p+find( SelectCol(p+1:PosEnd)==1,1)-1;
    if isempty(End) == 1,break;end
    Num = 3 - mod( End-Start+1,3 );
    if Num ~= 3
       SelectCol( p:p+Num-1) = ones( 1,Num );
    end
    Start = End + 1;
    if Test
        Start
        End
        p
        Num
    end
    p = End + find( SelectCol(End+1:PosEnd)==0,1 );
    if isempty( p )==0,break;end
  end
  CurrentSeq = CurrentSeq(:,find(SelectCol==1));
end
end
