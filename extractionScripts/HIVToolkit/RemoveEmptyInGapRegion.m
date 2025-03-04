function CurrentSeq = RemoveEmptyInGapRegion(  CurrentSeq,NoEmptyRow )

Test = 0;


%% remove the gap regions without any nucleotides 3 by 3.
[Row,Col] = size(CurrentSeq);
if nargin == 1
    NoEmptyRow = Row;
    for p = 1:Row
        if isempty(find(CurrentSeq(p,:)~='-',1))==1, NoEmptyRow = NoEmptyRow - 1; end
    end
end

SelectCol = ones(1,Col);
parfor p = 1:Col
    if isempty(find(CurrentSeq(:,p)~='-',1)) == 1
        SelectCol( p ) = 0;
    end
end

Start = find( SelectCol==1,1 );
p = find( SelectCol==0,1 );
if isempty(p) == 0
  while p <= Col
    End = p+find( SelectCol(p+1:Col)==1,1)-1;
    if isempty(End) == 1,break;end
    Num = 3 - mod(End-Start+1,3);
    if Num ~= 3
       SelectCol( p:p+Num-1) = ones(1,Num);       
    end
    Start = End + 1;
    if Test == 1
        Start
        End
        p
        Num
    end
    p = End + find( SelectCol(End+1:Col)==0,1 );
    if isempty( p )==0,break;end
  end
  CurrentSeq = CurrentSeq(:,find(SelectCol==1));
end

%% remove regions where sequence is less than 5%, and not able to assemble as amino acid
[~,Col] = size(CurrentSeq);

if mod(Col,3)==0,return;end

[ OccupyCol,OccupyNum ] = SequenceColumnOccupation( CurrentSeq,0.3,NoEmptyRow,CurrentSeq(1,:) );
OccupyNum = OccupyNum/NoEmptyRow;
p = 1; SelectedColumn = ones(1,Col);
while p <= Col-2
    End = p+find(OccupyCol(p+1:Col)==0,1)-1;
    if Test == 1
        Occupy = p:p+2
        Nu = OccupyNum(p:p+2)
        End
        Val = mod( End-p,3)
    end
    if isempty(End)==1,End = Col;end
    if OccupyNum(p)<0.05 && sum(OccupyNum(p+1:p+2)>0.9)==2 && mod( End-p,3) ==0
       SelectedColumn( p ) = 0;
       p = p + 1;
    elseif sum(OccupyNum(p:p+1)<0.05)==2 && OccupyNum(p+2)<0.05 && mod( End-p,3) ==0
       SelectedColumn( p:p+1 ) = [0 0];
       p = p + 2;
    else
       p = p + 3;
    end
end
if Test == 1
   ShowCol2 = find(SelectedColumn==0)
end
CurrentSeq = CurrentSeq( :,find(SelectedColumn==1) );

end
