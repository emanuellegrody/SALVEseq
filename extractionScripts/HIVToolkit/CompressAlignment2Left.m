function CurrentSeq = CompressAlignment2Left( CurrentSeq,NoEmptyRow,OccupyCutoff )
%% this function moves the non gap region to gap region if both amino acid columns hold similar probability.
%MutationCutoff = 0.3;
Test = 0;

[ Row,Col ] = size( CurrentSeq );
if nargin == 1
    NoEmptyRow = Row;
    OccupyCutoff = 0.1;
elseif nargin <= 2
    OccupyCutoff = 0.1;
end
if Test
   ResultDir = '/home/lowie/';
end

MutationCutoff = 0.1;
%% move the columns from right 2 left
[ OccupyCol,OccupyNum ] = SequenceColumnOccupation( CurrentSeq,OccupyCutoff,NoEmptyRow );
for p =1:3:Col-2
    if CheckStopFlag( OccupyCol(p:Col) ) == 0, break; end
    p
    if Test == 1
       Currentdd = p:p+2
       NoticedNum = sum( OccupyCol( p:p+2 ) ) 
       Flag = p>1&& OccupyCol(p-1)==1 || p==1
    end    
    if sum( OccupyCol( p:p+2 ) ) == 0 &&( p>1&& OccupyCol(p-1)==1 || p==1)
       %% the gap region is found in between p:p+2
       Right = p+2 + find( OccupyCol(p+3:Col)==1,1 );
       if isempty(Right)==1 || Right + 2 > Col, break; end
       if sum( OccupyCol(Right:Right+2) ) < 3
          continue;
       end
       % Right is the start of non-gap region, make sure that Right:Right
       % are all in non-gap region, which we have to move from right to left.
       RightEnd = Right + 1 + find( OccupyCol(Right+3:Col)==1,1 );
       % RightEnd is the end of gap region after Right:Right. where
       % we must move Right:Right+2 to left and adjust the gap region between Right+3:RightEnd
       if isempty(RightEnd)==1,RightEnd=Col; end
       
       %% SelectedRow = find(CurrentSeq(:,p)~='-');
       [ AA,~,AAProbLeft ] = Nucleotide2AA( CurrentSeq(:,p:p+2) );
       MutationLeft = find(AAProbLeft>0.1);
       [ ~,TempNum ] = SequenceColumnOccupation( AA );
       MaxLeft = FindDominateAminoAcid( AAProbLeft,TempNum,50,MutationCutoff,0.5 );

       [ AA,~,AAProbRight ] = Nucleotide2AA( CurrentSeq( find(CurrentSeq(:,p)=='-'),Right:Right+2) );       
       [ ~,TempNum ] = SequenceColumnOccupation( AA );
       MaxRight = FindDominateAminoAcid( AAProbRight,TempNum,50,MutationCutoff,0.5 );
       MutationRight = find( AAProbRight > 0.1);
          
       SelectedRow = find(CurrentSeq(:,p)~='-');
       [ ~,~,AAProb1 ] = Nucleotide2AA( CurrentSeq(SelectedRow,p:p+2) ); 
       [ ~,~,AAProb2 ] = Nucleotide2AA( CurrentSeq(SelectedRow,Right:Right+2) ); 
       
       if Test == 1
          Col
          Local1= p:p+2
          Local2 = Right:Right+2
          MutationRight
          MutationLeft
          M1 = MaxLeft{1}
          M2 = MaxRight{1}
          S1 = sum( OccupyNum(p+3:RightEnd)>0.3*Row )
          Right
          RightEnd          
          T1 = MaxLeft{1}(1) == MaxRight{1}(1)
          T2 = isempty( intersect(MaxLeft{1},MaxRight{1}) ) == 0
          T3 = isempty( setdiff(MutationLeft,MutationRight))
          
          V1 = sum( OccupyNum(Right+3:RightEnd)>0.3*Row )
          V2 = OccupyNum(Right-1)<0.3*Row
          V3 = AAProbLeft( MaxLeft{1} )
          V4 = find(AAProb1~=0)
          V5 = find(AAProb2~=0)
          
          B1 = intersect( find(AAProb1~=0),find(AAProbRight>0.05) )
          B2 = intersect( find(AAProb2~=0),find(AAProbRight>0.05) )
          B3 = isempty( intersect(find(AAProb1~=0),find(AAProbRight>0.05)) ) && isempty( intersect(find(AAProb2~=0),find(AAProbRight>0.05)) )
       end
       ModifySucess = 0;
       if (MaxLeft{1}(1) == MaxRight{1}(1) || isempty( intersect(MaxLeft{1},MaxRight{1}) ) == 0) && isempty( setdiff(MutationLeft,MutationRight))==1 ...
               || isequal( find(AAProb1~=0),find(AAProb2~=0) ) || isempty( intersect(find(AAProb1~=0),find(AAProbRight>0.05)) ) && isempty( intersect(find(AAProb2~=0),find(AAProbRight>0.05)) )
         % if p > 1 && OccupyCol(p-1)==1
          if Test == 1
             WriteSequence2Fasta( CurrentSeq,[],[ResultDir 'Compress_.fasta']);
          end
          if RightEnd >= Right+3 && sum( OccupyNum(Right+3:RightEnd)>0.3*Row ) > 0           
              [ CurrentSeq, p ] = PerfectCurrentSeq( CurrentSeq,p,Right+2,RightEnd );
              if Test == 1
                 Go = 'Inside change'
                 WriteSequence2Fasta( CurrentSeq,[],[ResultDir 'Compress0.fasta']);
              end 
          elseif isequal( find(AAProb1~=0),find(AAProb2~=0) ) || OccupyNum(Right-1)<0.3*Row && sum(AAProbLeft(MaxLeft{1}))>0.8 || ...
                  isempty( intersect(find(AAProb1~=0),find(AAProbRight>0.05)) ) && isempty( intersect(find(AAProb2~=0),find(AAProbRight>0.05)) )
              if Test == 1
                 Go = 'Normal change'
                 WriteSequence2Fasta( CurrentSeq,[],[ResultDir 'Compress0.fasta']);
              end
              CurrentSeq = MoveSeqRight2Left( CurrentSeq,p,RightEnd );
          end
          if Test == 1
             WriteSequence2Fasta( CurrentSeq,[],[ResultDir 'Compress4.fasta']);
          end   
          ModifySucess = 1;
       elseif length(setdiff( find(CurrentSeq(:,p)=='-'),find(CurrentSeq(:,RightEnd )=='-')))<=0.05*Row && MaxLeft{1}(1) == MaxRight{1}(1)
           % consider two columns are prefectly fit with each other.
          CurrentSeq = MoveSeqRight2Left( CurrentSeq,p,RightEnd );
          ModifySucess = 1;
      end
       
      if ModifySucess == 1
          CurrentSeq = TuneEmptyGapOut( CurrentSeq,p,RightEnd );
          Cutoff = Col - size( CurrentSeq,2 ); 
          Col = size( CurrentSeq,2 );
          p = p - Cutoff;
          [ OccupyCol,OccupyNum ] = SequenceColumnOccupation( CurrentSeq,OccupyCutoff,NoEmptyRow );           
      end
    end
end
CurrentSeq = RemoveEmptyInGapRegion(  CurrentSeq );
end


function StopFlag = CheckStopFlag( OccupyCol )
%% check whether there are two regions of gap regions. 
% StopFlag = 1: yes,  StopFlag = 0: no

StopFlag = 0;
Len = length( OccupyCol );
Empty = find( OccupyCol==0,1 );
if isempty( Empty ), return; end
Empty = Empty + find( OccupyCol(Empty+1:Len)==1,1 );
if isempty( Empty ), return; end
Empty = Empty + find( OccupyCol(Empty+1:Len)==0,1 );
if isempty( Empty ), return; end
StopFlag = 1;
end


function [ CurrentSeq,NextIndex ] = PerfectCurrentSeq( CurrentSeq,LeftStart,RightEnd,GapEnd )
% move the large sequence area to right side, and rare area to left side.
Test = 0;
ResultDir = '/home/lowie/HIV_DataSet/Full_Genome/AlignMafft/';

NextIndex = LeftStart;
%% first, move current 3 columns in Right:Right+2 to p:p+2
CurrentSeq = MoveSeqRight2Left( CurrentSeq,LeftStart,RightEnd );
[ Row,Col ] = size( CurrentSeq );
CollectRow = zeros(1,Row);
for p = 1:Row
    if p == 1
       Neighbor = 2:4;
    elseif p == 2
       Neighbor = [1 3 4];
    elseif p == Row-1
       Neighbor = [Row-2 Row-3 Row];
    elseif p == Row
       Neighbor = Row-3:Row-1;
    else
       Neighbor = [p-2:p-1 p+1:p+2];
    end
    if isempty( find(CurrentSeq(Neighbor,RightEnd+1)~='-',1) ) == 0
       CollectRow( p ) = 1;
    end
end

%% CollectRow colloect the rows that connect with each other.
p = 1;
while p<=Row
   if CollectRow(p)==1
      End = p + find( CollectRow(p+1:Row)==0,1 );
      if isempty(End)==1,End=Col; end
      if End-p+1 < 10
         CollectRow(p:End) = zeros(1,End-p+1);
      end
      p = End + 1;
   elseif CollectRow(p)==0
      p = p + find( CollectRow(p+1:Row)==1,1);
      if isempty(p),break; end
   end
end

%% move the individual nucleotide from right to left for those unselected rows (so many gaps in such region)
CurrentCollectRowLeft = find( CollectRow==0 );
if Test == 1
    LeftStart
    RightEnd
    CollectRow1 = CurrentCollectRowLeft
end
if isempty( CurrentCollectRowLeft ) == 0 && length( CurrentCollectRowLeft ) < Row
   %FurtherEnd = RightEnd + 1 + find( OccupyCol( RightEnd+3:Col)==1,1);
   if Test == 1
      Go = 'right 2 left'
      LeftStart
      GapEnd
      WriteSequence2Fasta( CurrentSeq,[],[ResultDir 'Compress2.fasta']);
   end
    CurrentSeq = MoveSeqRight2Left( CurrentSeq,LeftStart,GapEnd,CurrentCollectRowLeft );
end

if Test == 1
   WriteSequence2Fasta( CurrentSeq,[],[ResultDir 'Compress3.fasta']);
end

%% move the aligned sequences from left to right
CurrentCollectRowRight = find( CollectRow==1 );
if isempty( CurrentCollectRowRight ) == 0
   if Test == 1
      Go = 'left 2 right'
      CollectRow2 = CurrentCollectRowRight
   end  
   
   CutGapRegion = CurrentSeq(CurrentCollectRowRight,RightEnd+1:GapEnd);
   Len = size( CutGapRegion,2 ); ColSelect = ones( 1,Len );
   for p = 1:Len
       if isempty(find(CutGapRegion(:,p)~='-',1))==1
          ColSelect( p ) = 0;
       end
   end
   CutGapRegion = CutGapRegion(:,find(ColSelect==1));
   for p = LeftStart:GapEnd
       if isempty( find(CurrentSeq(CurrentCollectRowLeft,p)~='-',1) ) == 1
          ModifiedEnd = p; break;
       end
   end
   Len = size( CutGapRegion,2 );
   if ModifiedEnd + Len <= GapEnd
      CurrentSeq( CurrentCollectRowRight,RightEnd+1:GapEnd ) = repmat('-',length(CurrentCollectRowRight),GapEnd-RightEnd);
      CurrentSeq( CurrentCollectRowRight,GapEnd-Len+1:GapEnd ) = CutGapRegion;
   else
      if Test
         Len
         ModifiedEnd
         GapEnd
      end
      CurrentSeq = MoveSeqRight2Left( CurrentSeq,LeftStart,GapEnd );
      
      if 0
      
      CurrentSeq = MoveSeqRight2Left( CurrentSeq,LeftStart,GapEnd,CurrentCollectRowLeft );
      
      AssembleLeft = repmat('-',Row,ModifiedEnd-LeftStart+1 );
      AssembleLeft( CurrentCollectRowLeft,: ) = CurrentSeq( CurrentCollectRowLeft,LeftStart:ModifiedEnd );
      AssembleRight = repmat( '-',Row,Len );
      AssembleRight( CurrentCollectRowRight,1:Len ) = CutGapRegion;
      
      CurrentSeq = [CurrentSeq(:,1:LeftStart-1), AssembleLeft, AssembleRight,CurrentSeq(:,RightEnd+1:Col)  ];
          if Test == 1
             WriteSequence2Fasta( CurrentSeq,[],[ResultDir 'Compress5.fasta']);
          end       
        d
      end
   end
   NextIndex = GapEnd;
end
   
end

