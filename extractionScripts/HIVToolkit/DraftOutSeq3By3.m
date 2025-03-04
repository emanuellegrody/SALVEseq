function CurrentSeq = DraftOutSeq3By3( CurrentSeq,NoEmptyRow,OccupyCutoff )

Test = 0;
[ Row,Col ] = size( CurrentSeq );

if nargin == 1
    NoEmptyRow = Row;
    OccupyCutoff = 0.8;
    Reference = CurrentSeq(1,:);
end

Reference = CurrentSeq(1,:); %I added this
OccupyCol = SequenceColumnOccupation( CurrentSeq,OccupyCutoff,NoEmptyRow,Reference );

p = 1;
while p < Col-2
    % if 3 adjacent columns are either occupied or empty. 
    if sum( OccupyCol( p:p+2 ) ) == 3 || sum( OccupyCol( p:p+2 ) ) == 0
       p = p + 3;
    % we have to adjust the next three occupied position
    else
       Pos = p + find( OccupyCol(p:Col),3,'first') - 1 ;% column positions of three nucleotides.
       if length(Pos)<3,break; end
       Lower = Pos(1);
       Upper = Pos(3);
       Difference = zeros(1,3);
       if Test == 1
          Lower
          Upper
          p
       end
       
       % a). converge to left side: --A---B----C--  => --ABC---------
       if Pos(1)+1 ~= Pos(2) % --A---B----C--  => --ABC---------
          OrigSeqPop1 = MoveSeqRight2Left( CurrentSeq,Pos(1),Pos(3) );
          Difference(1) = SeqPenalty( CurrentSeq( :,Lower+1:Lower+2 ),OrigSeqPop1( :,Lower+1:Lower+2 ) );
       else % --AB----C--  => --ABC---------
          OrigSeqPop1 = MoveSeqRight2Left( CurrentSeq,Pos(2),Pos(3) );
          Difference(1) = SeqPenalty( CurrentSeq( :,Lower+2 ),OrigSeqPop1( :,Lower+2 ) );
       end
       % NewScore = Difference(1)
       
       % b). converge to right side: --A---B----C--  => ---------ABC--
       if Pos(2)+1 ~= Pos(3) % --A---B----C--  => ---------ABC--
          OrigSeqPop2 = MoveSeqLeft2Right( CurrentSeq,Pos(1),Pos(3) );
          Difference(2) = SeqPenalty( CurrentSeq( :,Upper-2:Upper-1 ),OrigSeqPop2( :,Upper-2:Upper-1 ) ) ;
       else % --A---BC--  => ---------ABC--
          OrigSeqPop2 = MoveSeqLeft2Right( CurrentSeq,Pos(1),Pos(2) );
          Difference(2) = SeqPenalty( CurrentSeq( :,Upper-2 ),OrigSeqPop2( :,Upper-2 ) ) ;
       end
       
       % c). converge to middle: --A---B----C--  => -----ABC------
       OrigSeqPop3 = CurrentSeq;
       if Pos(1)+1 ~= Pos(2) % --A---B----C--  => -----AB----C--
          OrigSeqPop3 = MoveSeqLeft2Right( CurrentSeq,Pos(1),Pos(2) );
          Difference(3) = SeqPenalty( CurrentSeq( :,Pos(2)-1 ),OrigSeqPop3( :,Pos(2)-1 ) ) ;
       end
       
       if Pos(2)+1~=Pos(3) % -----AB----C-- => -----ABC------
          OrigSeqPop3 = MoveSeqRight2Left( OrigSeqPop3,Pos(2),Pos(3) );
          TempScore = SeqPenalty( CurrentSeq( :,Lower+2 ),OrigSeqPop3( :,Lower+2 ) );
          if Pos(1)+1 ~= Pos(2)
             Difference(3) = ( Difference(3) + TempScore)/2 ;
          else
             Difference(3) = TempScore;
          end
       end
       %Pos
       %Difference
       %p
       [ ~,Max ] = max( Difference );
       if Max == 1 % the best solution is to move 
          CurrentSeq = OrigSeqPop1;
       elseif Max == 2
          CurrentSeq = OrigSeqPop2;
       elseif Max == 3
          CurrentSeq = OrigSeqPop3;
       end
       
       OccupyCol = SequenceColumnOccupation( CurrentSeq, OccupyCutoff, NoEmptyRow );

       p = Pos(1);
    end
end
end