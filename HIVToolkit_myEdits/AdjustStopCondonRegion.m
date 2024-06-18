 function CurrentSeq = AdjustStopCondonRegion( CurrentSeq )
 % this function modify the end of sequences with regard to the stop codons.
Test = 0; Test1 = 0; Test2 = 0;

[ Row,Col ] = size( CurrentSeq ) ; 

End = find( CurrentSeq(1,:)~='-',1,'last' );
FirstSequenceEnd = 0;
if ~(isempty( find(CurrentSeq( 1,End-2:End ) =='-',1) ) && isequal( nt2aa( CurrentSeq( 1,End-2:End ),'ACGTOnly', false, 'AlternativeStartCodons', false),'*' ) )
    return;
end

%% if the first reference sequence has a stop codon at the end, look for stop codon positions at the end.
if Test
    Type = 'Stop condon *'
end
   CheckStep = 6;
   FirstStopCondonPos = Col; 
   for p = 1:Row
       NonGapPos = End - CheckStep + find( CurrentSeq(p,End-CheckStep+1:End) ~= '-',9,'last' );
       if isempty(NonGapPos)==1, continue; end
       %%%%%%%%%%%%%%%%%%% find the position of stop codon %%%%%%%%%%%%%
       q = length( NonGapPos )-2;
       while q > 0
            %Pos = NonGapPos( q:q+2 );
            %Lo = CurrentSeq(p,NonGapPos(q:q+2))
            if isequal( nt2aa(CurrentSeq(p,NonGapPos(q:q+2)),'ACGTOnly', false, 'AlternativeStartCodons', false),'*' )
               %SearchEnd = NonGapPos( q+2 )
               if NonGapPos( q+2 )+1 < Col
                  CurrentSeq(p,NonGapPos( q+2 )+1:Col) = repmat('-',1,Col-NonGapPos( q+2 ));
               end
               if FirstStopCondonPos > NonGapPos( q+2 )
                  FirstStopCondonPos = NonGapPos( q+2 );
               end
               if p == 1
                   FirstSequenceEnd = NonGapPos( q+2 );
               end
               break;
            end
            q = q - 3;
       end
       if Test == 1
          p
          End
          Col
          FirstSequenceEnd1 = FirstSequenceEnd
          FirstStopCondonPos
          if FirstStopCondonPos == 70
             FindRow = p
             d
          end          
       end
       
       if q <= 0 && End+1 < Col
          LocalEnd = End-11+find( CurrentSeq( p,End-10:End ) ~= '-',3,'last' );
          if isempty(LocalEnd)==1,continue;end
          NonGapPos = LocalEnd(1) - 1 + find( CurrentSeq( p,LocalEnd(1):Col ) ~='-',30 );
          q = 1;
          while q <= length( NonGapPos )-2
               %Pos = NonGapPos(q:q+2)
               %Lo = CurrentSeq(p,NonGapPos(q:q+2))
               if isequal( nt2aa(CurrentSeq(p,NonGapPos(q:q+2)),'ACGTOnly', false, 'AlternativeStartCodons', false),'*' )
                  if NonGapPos( q+2 )+1 < Col
                     CurrentSeq( p,NonGapPos( q+2 )+1:Col ) = repmat('-',1,Col-NonGapPos( q+2 ));
                  end
                  if FirstStopCondonPos > NonGapPos( q+2 )
                     FirstStopCondonPos = NonGapPos( q+2 );
                  end
                  if p == 1
                     FirstSequenceEnd = NonGapPos( q+2 );
                  end
                  break;
               end
               q = q + 3;
          end
          if Test == 1
             FirstSequenceEnd2 = FirstSequenceEnd
             FirstStopCondonPos
             if FirstStopCondonPos == 70
                FindRow = p
                d
             end
          end
          if q > length( NonGapPos )-2 % non * is found at the end of the sequences.
              CurrentSeq( p,FirstSequenceEnd+1:Col ) = repmat('-',1,Col-FirstSequenceEnd );
          end
       end
    end
   if 0
    OccupyCol = SequenceColumnOccupation( CurrentSeq,OccupyCutoff,NoEmptyRow );
    Start = find( OccupyCol(1:FirstStopCondonPos-3)==1,1,'last' ) + 1;
    if Test == 1
        %OccupyCol
        [ Row,Col ] = size(CurrentSeq)
        FirstStopCondonPos
        Start
    end
    if isempty( Start ) == 0 
        CurrentSeq = MoveSeqRight2Left( CurrentSeq,Start,Col );
    end
   end
if Test == 1
   ResultDir = '/home/lowie/HIV_DataSet/Full_Genome/AlignMafft/';
   WriteSequence2Fasta( CurrentSeq,[],[ResultDir 'Assemble1.fasta'] ); 
end


%% if there are so many gaps near to the stop condon, move columns from right to left
% example:  AAA-----AA*         AAAAA*-----         AAAAA*
%           AAA-----AA*  ===>   AAAAA*-----   ===>  AAAAA*
%           CCCAA*TTTTT         CCCAA*TTTTT         CCCAA*
%           CCCAA*TTTTT         CCCAA*TTTTT         CCCAA*
% we first identify the location of gap, then test whether the newly
% aligned columns can guarattee columns are correctly aligned.
[ ~,OccupyNum ] = SequenceColumnOccupation( CurrentSeq );
GapEnd = find( OccupyNum < 0.7*OccupyNum( End ),1,'last' );
if isempty( GapEnd ), return; end
GapStart = find( OccupyNum( 1:GapEnd ) > 0.7*OccupyNum( End ),1,'last' ) + 1;
if isempty( GapStart ), return; end

if Test1
    GapStart
    GapEnd
    End
end

if Col - GapStart < 20
    
Fail = 0;
for p = GapEnd:3:End-3
    [ ~,~,NonAAProb1 ] = Nucleotide2AA( CurrentSeq( :,p+1:p+3 ) );
    [ ~,~,NonAAProb2 ] = Nucleotide2AA( CurrentSeq( :,GapStart+p-GapEnd:GapStart+p-GapEnd+2 ) );
    if Test1
       Loc1 = p+1:p+3
       Loc2 = GapStart+p-GapEnd:GapStart+p-GapEnd+2
       Mut1 = find(NonAAProb1 > 0.1)
       Mut2 = find(NonAAProb2 > 0.1)
    end
    if isempty( intersect( NonAAProb1 > 0.1,NonAAProb2 > 0.1 ) )
        if Test1
           FailCol = [ p p+2 GapStart+p-GapEnd GapStart+p-GapEnd+2 ]
        end
        Fail = 1; break;
    end
end
if Fail == 0 %&& Col-GapStart < 90
    TempCurrentSeq = MoveSeqRight2Left( CurrentSeq,GapStart,End );
    [ ~,~,NonAAProb1 ] = Nucleotide2AA( TempCurrentSeq( :,GapStart+End-GapEnd-3:GapStart+End-GapEnd-1 ) );
    if Test1
       Case = 'Type 1: end condon'
       NonAAProb1
       P1 = GapStart+End-GapEnd-1
       End
       ResultDir = '/home/lowie/HIV_DataSet/Full_Genome/AlignMafft/';
       WriteSequence2Fasta( TempCurrentSeq,[],[ResultDir 'Assemble2.fasta'] ); 
    end
    if NonAAProb1( 24 ) > 0.6
       if Test1
          Case = 'Success'
       end
       CurrentSeq = TempCurrentSeq( :,1:GapStart+End-GapEnd-1 );
    end
end
end

if Test1 == 1 || Test2 == 1
   ResultDir = '/home/lowie/HIV_DataSet/Full_Genome/AlignMafft/';
   WriteSequence2Fasta( CurrentSeq,[],[ResultDir 'Assemble3.fasta'] ); 
   Step = 'The last step:'
end

%% improve the alignment of the stop codon, if exists.
End = find( CurrentSeq(1,:)~='-',1,'last' );
%[ Row,End ] = size( CurrentSeq );
for p = 1:Row
    LocalEnd = find( CurrentSeq(p,:)~='-',1,'last' );
    if LocalEnd < End - 18,continue; end
    Find = 0;
    for q = LocalEnd:-3:LocalEnd-9
        if mod( q,3 ) == 0 && isempty(find(CurrentSeq(p,q-2:q)=='-',1)) && isequal( nt2aa( CurrentSeq( p,q-2:q ),'ACGTOnly', false, 'AlternativeStartCodons', false),'*' )
           if Test2
              Replace = q-2:q
              End
           end
           Find = q;
        end
    end
    if Find > 0
        if Test2
            p
            Find
            End
        end
        if Find >= End
            CurrentSeq(p,Find+1:LocalEnd)=repmat('-',1,LocalEnd-Find);
            if Test2
               MoveNum = LocalEnd-Find
            end
        else
            Key= CurrentSeq(p,Find-2:Find);
            CurrentSeq(p,Find-2:LocalEnd)=repmat('-',1,LocalEnd-Find+3);
            CurrentSeq(p,End-2:End)=Key;
            if Test2
                Key
                Find-2:LocalEnd
                End-2:End
                MoveNum = LocalEnd-Find+3
            end            
        end
    end
end
end
