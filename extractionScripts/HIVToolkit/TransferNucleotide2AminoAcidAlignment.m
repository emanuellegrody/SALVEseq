function Final_CurrentSeq = TransferNucleotide2AminoAcidAlignment( Orig_CurrentSeq,CurrentTitle,ScoreMatrix )
% s = matlabpool('size');
% if s == 0, matlabpool open; end

WriteTempSequenceOut = 2;
if WriteTempSequenceOut == 2
   %currentFolder = pwd;
   %addpath( genpath( currentFolder ) )
   ResultDir = './TemporaryData/' ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% read the sequences where there is no empty sequences in our analysis.
if nargin==3 && length( CurrentTitle ) ~= size( 2,Orig_CurrentSeq )
    Error = 'title number is not equal to sequence number'
end

SelectedRow = zeros( 1,size(Orig_CurrentSeq,1) );
for p = 1:size(Orig_CurrentSeq,1)
    if length( find(Orig_CurrentSeq(p,:)~='-') ) > 2
   % if isempty( find(Orig_CurrentSeq(p,:)~='-',1) ) == 0
       SelectedRow( p ) = 1;
    end
end

SelectedRow = find( SelectedRow > 0 );
if isempty( SelectedRow ) == 1
   Final_CurrentSeq = repmat('-',size(Orig_CurrentSeq,1),1);
   return;
end
CurrentSeq = Orig_CurrentSeq( SelectedRow,: );

%% look at reference sequence, move columns to fit the transformation from nucleotide to AA (3 nucleotides to one AA)
[ Row,Col ] = size( CurrentSeq );
if nargin < 3
    ScoreMatrix = blosum62;
end
if nargin < 2
   CurrentTitle = cell( 1,Row );
   for p = 1:Row
       CurrentTitle{p} = ['>' num2str( p )];
   end
end
NoEmptyRow = Row;
for p = 1:Row
    if isempty(find(CurrentSeq(p,:)~='-',1))==1, NoEmptyRow = NoEmptyRow - 1; end
end

%% step 0, remove empty region at the end of sequences.
[ ~,Col ] = size( CurrentSeq );
for p = Col:-1:1
    if isempty( find(CurrentSeq(:,p)~='-',1)) == 0, break; end
end
CurrentSeq = CurrentSeq( :,1:p ); Col = p;

if WriteTempSequenceOut == 2
   WriteSequence2Fasta( CurrentSeq,CurrentTitle, [ResultDir 'Orig.fasta'] );    
end

%% step 1: put all nucleotides 3 by 3, to make AA sequences.
OccupyCutoff = 0.50; 
OccupyCol = SequenceColumnOccupation( CurrentSeq,OccupyCutoff,NoEmptyRow);
while 1
    if mod(sum(OccupyCol),3) == 0,break; end
    OccupyCutoff = OccupyCutoff + 0.01;
    if OccupyCutoff > 1
        Error = 'Check out the dataset, cannot aligned into amino acids.'
    end
    OccupyCol = SequenceColumnOccupation( CurrentSeq,OccupyCutoff,NoEmptyRow);
    %if OccupyCutoff>0.5,break;end
end
OccupyCutoff
CurrentSeq = CheckoutStopCondonFirstSequence( CurrentSeq );

CurrentSeq = RemoveSingleEmptyCol( CurrentSeq );
if 0
CurrentSeq = ModifyNucleotideNonGapRegion( CurrentSeq );
end
CurrentSeq = DraftOutSeq3By3( CurrentSeq,NoEmptyRow,OccupyCutoff );
if WriteTempSequenceOut == 1
   Procedure = 'Step 1'
elseif WriteTempSequenceOut == 2
   Procedure = 'Step 1'
   WriteSequence2Fasta( CurrentSeq,CurrentTitle, [ResultDir 'Tran1.fasta'] );    
end

%% Step 2: move the ---- AB----CCCC---- ==> -----ABC----CCC---- or ---------ABCCCC----
%CurrentSeq = LocalRemoveGapRegion( CurrentSeq );

%OccupyCol = SequenceColumnOccupation( CurrentSeq,0.1 );
CurrentSeq = Arrange1or2Nucleotide( CurrentSeq );
if 0
OccupyCol = SequenceColumnOccupation( CurrentSeq,OccupyCutoff );
Start = find( OccupyCol==1,1 )-1; KeepHeadSeq = [];
if isempty( Start ) == 0
    KeepHeadSeq = CurrentSeq( :,1:Start );
    CurrentSeq = CurrentSeq( :,Start+1:size(CurrentSeq,2) );
end

if WriteTempSequenceOut == 1
   Procedure = 'Step 2'   
elseif WriteTempSequenceOut == 2
   Procedure = 'Step 2'
   WriteSequence2Fasta( CurrentSeq,CurrentTitle, [ResultDir 'Tran2.fasta'] );
end

end

%% step 3: align the regions to left or right, if the occupy column is not 3 by 3, 
CurrentSeq = AdjustLocal3By3Nucleotide(CurrentSeq,NoEmptyRow );
if WriteTempSequenceOut == 1
   Procedure = 'Step 3a'   
elseif WriteTempSequenceOut == 2
   Procedure = 'Step 3a'
   WriteSequence2Fasta( CurrentSeq,CurrentTitle, [ResultDir 'Tran3.fasta'] );
end


CurrentSeq = Arrange1or2Nucleotide( CurrentSeq );
if WriteTempSequenceOut == 1
    Procedure = 'Step 3b'   
elseif WriteTempSequenceOut == 2
   Procedure = 'Step 3b'
   WriteSequence2Fasta( CurrentSeq,CurrentTitle, [ResultDir 'Tran4.fasta'] );
end

clear OrigSeqPop1 OrigSeqPop2 OrigSeqPop3 OrigSeqPop4

CurrentSeq = MoveGapRegionSeq2AlignedRegion( CurrentSeq,NoEmptyRow );
if WriteTempSequenceOut == 1
    Procedure = 'Step 3c'   
elseif WriteTempSequenceOut == 2
   Procedure = 'Step 3c'
   WriteSequence2Fasta( CurrentSeq,CurrentTitle, [ResultDir 'Tran5.fasta'] );
end

%% step 4: adjust the 3 nucleotides aligned with other sequences.
% example:  ---AAA---BBB---    ---AAA---BBB---
%           --------AAA----    --------BBB----
%    ===>   ---AAA------       ---------BBB---
% [ TransMatrix, ProbMatrix ] = Nucleotide2Int( CurrentSeq );
CurrentSeq = AdjustLocal3By3Nucleotide(CurrentSeq,NoEmptyRow );
if WriteTempSequenceOut == 1
   Procedure = 'Step 3d'
elseif WriteTempSequenceOut == 2
   Procedure = 'Step 3d'
   WriteSequence2Fasta( CurrentSeq,CurrentTitle, [ResultDir 'Tran6.fasta'] );
end


CurrentSeq = LocalRemoveGapRegion( CurrentSeq );
if WriteTempSequenceOut == 1
    Procedure = 'Step 3e'   
elseif WriteTempSequenceOut == 2
   Procedure = 'Step 3e'
   WriteSequence2Fasta( CurrentSeq,CurrentTitle, [ResultDir 'Tran7.fasta'] );
end



%% step 5, remove unknown amino acids in gap regions
CurrentSeq = FitAlignmentByCompress( CurrentSeq,NoEmptyRow,OccupyCutoff );
if WriteTempSequenceOut == 1
   Procedure = 'Step 4'   
elseif WriteTempSequenceOut == 2
   Procedure = 'Step 4'
   WriteSequence2Fasta( CurrentSeq,CurrentTitle,[ResultDir 'Tran8.fasta']);
end

if 0
CurrentSeq = OptimalSingleAminoAcidPosition( CurrentSeq,NoEmptyRow,OccupyCutoff );
if WriteTempSequenceOut == 1
   Procedure = 'Step 5'   
elseif WriteTempSequenceOut == 2
   Procedure = 'Step 5'
   WriteSequence2Fasta( CurrentSeq,CurrentTitle,[ResultDir 'Tran8_1.fasta']);
end
end

CurrentSeq = MoveGapRegionSeq2AlignedRegion( CurrentSeq,NoEmptyRow );
if WriteTempSequenceOut == 1
   Procedure = 'Step 6'   
elseif WriteTempSequenceOut == 2
   Procedure = 'Step 6'
   WriteSequence2Fasta( CurrentSeq,CurrentTitle,[ResultDir 'Tran9.fasta']);
end

CurrentSeq = CompressAlignment2Left( CurrentSeq,NoEmptyRow );
%CurrentSeq = CompressAlignment2Left( CurrentSeq,NoEmptyRow,OccupyCutoff );
if WriteTempSequenceOut == 1
   Procedure = 'Step 8'
elseif WriteTempSequenceOut == 2
   Procedure = 'Step 8'
   WriteSequence2Fasta( CurrentSeq,CurrentTitle,[ResultDir 'Tran10.fasta']);
end

CurrentSeq = CompressAlignment2Right( CurrentSeq,NoEmptyRow,OccupyCutoff );
if WriteTempSequenceOut == 1
   Procedure = 'Step 8'
elseif WriteTempSequenceOut == 2
   Procedure = 'Step 8'
   WriteSequence2Fasta( CurrentSeq,CurrentTitle,[ResultDir 'Tran11.fasta']);
end

CurrentSeq = Arrange1or2Nucleotide( CurrentSeq );
if WriteTempSequenceOut == 1
   Procedure = 'Step 9'
elseif WriteTempSequenceOut == 2
   Procedure = 'Step 9'
   WriteSequence2Fasta( CurrentSeq,CurrentTitle,[ResultDir 'Tran12.fasta']);
end

CurrentSeq = OptimizedAminAcidByScoreMatrix( CurrentSeq,NoEmptyRow,OccupyCutoff,ScoreMatrix );
if WriteTempSequenceOut == 1
   Procedure = 'Step 9'
elseif WriteTempSequenceOut == 2
   Procedure = 'Step 9'
   WriteSequence2Fasta( CurrentSeq,CurrentTitle,[ResultDir 'Tran13.fasta']);
end

CurrentSeq = CompressAlignment2Left( CurrentSeq,NoEmptyRow,OccupyCutoff );
if WriteTempSequenceOut == 1
   Procedure = 'Step 10'
elseif WriteTempSequenceOut == 2
   Procedure = 'Step 10'
   WriteSequence2Fasta( CurrentSeq,CurrentTitle,[ResultDir 'Tran18.fasta']);
end

CurrentSeq = AdjustLocal3By3Nucleotide(CurrentSeq,NoEmptyRow );
if WriteTempSequenceOut == 1
   Procedure = 'Step 10'
elseif WriteTempSequenceOut == 2
   Procedure = 'Step 10'
   WriteSequence2Fasta( CurrentSeq,CurrentTitle,[ResultDir 'Tran14.fasta']);
end


CurrentSeq = LocalRemoveGapRegion( CurrentSeq );
if WriteTempSequenceOut == 1
   Procedure = 'Step 11'
elseif WriteTempSequenceOut == 2
   Procedure = 'Step 11'
   WriteSequence2Fasta( CurrentSeq,CurrentTitle,[ResultDir 'Tran15.fasta']);
end

%CurrentSeq = CultivateSeqByScoreMatrix( CurrentSeq,ScoreMatrix );
if WriteTempSequenceOut == 1
   Procedure = 'Step 11'
elseif WriteTempSequenceOut == 2
   Procedure = 'Step 11'
   WriteSequence2Fasta( CurrentSeq,CurrentTitle,[ResultDir 'Tran16.fasta']);
end

%CurrentSeq = OptimalMultiAminoAcidPosition( CurrentSeq,OccupyCutoff );
if WriteTempSequenceOut == 1
   Procedure = 'Step 11'
elseif WriteTempSequenceOut == 2
   Procedure = 'Step 11'
   WriteSequence2Fasta( CurrentSeq,CurrentTitle,[ResultDir 'Tran17.fasta']);
end

%% step 6, adjust the stop condon region, if they are observed.
% in this step, all amino acids are aligned 3 by 3.
%CurrentSeq = OptimalSeperatedSeq( CurrentSeq );

CurrentSeq = AlignGapInStopCodonRegion( CurrentSeq );
%CurrentSeq = TuneEmptyGapOut( CurrentSeq);
CurrentSeq = LocalRemoveGapRegion( CurrentSeq );
%CurrentSeq = [ KeepHeadSeq CurrentSeq ];
%CurrentSeq = LocalRemoveGapRegion( CurrentSeq );
if WriteTempSequenceOut == 1
   Procedure = 'Step 12'
elseif  WriteTempSequenceOut == 2
   Procedure = 'Step 12'
   WriteSequence2Fasta( CurrentSeq,CurrentTitle,[ResultDir 'Tran18.fasta']);
end

%% recover back the orignal datasets by filling the empty rows if exist.
Final_CurrentSeq = repmat( '-',size(Orig_CurrentSeq,1),size(CurrentSeq,2) );
Final_CurrentSeq( SelectedRow,: ) = CurrentSeq;


%% remove the temporary files
if  WriteTempSequenceOut == 2
    delete( [ResultDir 'Orig.fasta'] );
    for p = 1:24
        delete( [ResultDir 'Tran' num2str(p) '.fasta'] );
    end
end
end

function CurrentSeq = CheckoutStopCondonFirstSequence( CurrentSeq )
% this function check the stop codon of the first sequence.
Test = 0;
Col = size( CurrentSeq,2 );
End = find( CurrentSeq(1,:)~='-',3,'last' );
if End(3)==Col || ~(isempty(End)==0 && isequal( nt2aa( CurrentSeq( 1,End ),'ACGTOnly', false, 'AlternativeStartCodons', false),'*'))
    return;
end
Count = 0;
while isempty(End)==0 && length(End)==3
    Location = Col-41 + strfind( CurrentSeq(2,Col-40:Col), CurrentSeq(1,End) );
    if Test
       Str1 = CurrentSeq(2,Col-40:Col)
       Str2 = CurrentSeq(1,End)
       Location
    end
    if isempty(Location), break; end
    Location = Location( 1 );
    if Location <= End(1),break; end

    Temp = CurrentSeq( 1,End );
    CurrentSeq( 1,End ) = '---';
    CurrentSeq( 1,Location:Location+2 ) = Temp;
    if Count == 0
        [ ~,~,NonAAProb ] = Nucleotide2AA( CurrentSeq(:,Location+3:Location+5) );
        if NonAAProb( 24 ) < 0.2
            CurrentSeq = CurrentSeq( :,1:Location+2 );
            Col = size( CurrentSeq,2 );
        end
        Count = 1;
    end
    End = find( CurrentSeq(1,1:End(1)-1)~='-',3,'last');
end
end

function CurrentSeq = ModifyNucleotideNonGapRegion( CurrentSeq )
% this function moves the segment of nucleotides from gap region to non-gap region.
Test = 0; Test1 = 0; Test2 = 0;
[ Row,Col ] = size( CurrentSeq );
OccupyCutoff = 0.5;
OccupyCol = SequenceColumnOccupation( CurrentSeq,OccupyCutoff );
GapStart = find( OccupyCol==0,1 );   % the beginning of gap

while isempty( GapStart ) == 0 && GapStart <= Col
    GapEnd = GapStart + find( OccupyCol(GapStart+1:Col)==1,1 ) - 1;
    if isempty(GapEnd), GapEnd = Col; end
    NonGapEnd_Left = find( OccupyCol( 1:GapStart-1)==1,1,'last' );
    if isempty(NonGapEnd_Left), NonGapEnd_Left = 0; end
    NonGapStart_Left = find( OccupyCol( 1:NonGapEnd_Left ) == 0,1,'last' ) + 1;
    if isempty( NonGapStart_Left ), NonGapStart_Left = 0; end
    
    NonGapStart_Right = GapEnd + find( OccupyCol( GapEnd+1:Col )==1,1 );
    if isempty(NonGapStart_Right), NonGapStart_Right = 0; end
    NonGapEnd_Right = NonGapStart_Right + find( OccupyCol( NonGapStart_Right+1:Col ) == 0,1 ) - 1;
    if isempty( NonGapEnd_Right ), NonGapEnd_Right = Col; end
    if Test
        Go = 'New'
        GapStart
        GapEnd
        NonGapStart_Left
        NonGapEnd_Left
        NonGapStart_Right
        NonGapEnd_Right
    end
    for q = 1:Row
        % move nucleotide sequences in the left side
        if NonGapStart_Left > 0
           LeftEnd = NonGapEnd_Left - 1 + find( CurrentSeq(q,NonGapEnd_Left:GapEnd) ~= '-',1,'last' );
           if isempty(LeftEnd), continue; end
           LeftStart = NonGapStart_Left + find( CurrentSeq(q,NonGapStart_Left:LeftEnd)=='-',1,'last' );
           if isempty( LeftStart ), continue; end
           LeftBound = NonGapStart_Left + find( CurrentSeq( q,NonGapStart_Left:LeftStart-1 )~='-',1,'last' );
           if isempty( LeftBound ), LeftBound = NonGapStart_Left; end

           if NonGapEnd_Left - LeftBound >= LeftEnd - LeftStart
              Segment = CurrentSeq( q,LeftStart:LeftEnd );
              ComparedRow = find( CurrentSeq(1:q-1,NonGapEnd_Left)~='-',1,'last' );
              if isempty( ComparedRow )
                 ComparedRow = q + find( CurrentSeq(q+1:Row,NonGapEnd_Left)=='-',1 );
                 if isempty(ComparedRow) == 1, continue; end
              end
           if Test
              LeftStart
              LeftEnd
              LeftBound
           end              
              ComparedSegment = CurrentSeq( ComparedRow,GapStart-length(Segment):GapStart-1 );
              if Test
                 Ro = q
                 ComparedRow
                 Segment
                 ComparedSegment
              end
             if SimilarityBetweenTwoSeq( Segment,ComparedSegment ) > 0.8
                Temp = CurrentSeq( q,LeftStart:LeftEnd );
                CurrentSeq( q,LeftStart:LeftEnd ) = repmat('-',1,LeftEnd-LeftStart+1 );
                CurrentSeq( q,GapStart-(LeftEnd-LeftStart+1):GapStart-1 ) = Temp;
                if Test
                   Cu = 'success'
                   Temp
                end
                continue;
             end
           end
        end  
        % move nucleotide sequences in the right side
        if NonGapStart_Right > 0
           RightStart = GapStart - 1 + find( CurrentSeq( q,GapStart:Col )~='-',1 );
           if RightStart >= NonGapStart_Right, continue; end
           RightEnd = RightStart + find( CurrentSeq( q,RightStart+1:Col ) == '-',1 )-1;
           if isempty(RightEnd), RightEnd = Col; end
           RightBound = RightEnd + find( CurrentSeq(q,RightEnd+1:Col)~='-',1 ) - 1;
           if isempty( RightBound ),RightBound = Col; end
           if Test1
               RightStart
               RightEnd
               RightBound
           end
           Segment = CurrentSeq( q,RightStart:RightEnd );
           ComparedRow = find( CurrentSeq(1:q-1,NonGapStart_Right)~='-',1,'last' );
           if isempty( ComparedRow )
               ComparedRow = q + find( CurrentSeq(q+1:Row,NonGapStart_Right)=='-',1 );
               if isempty(ComparedRow) == 1, continue; end
           end  
           if RightBound - NonGapStart_Right >= RightEnd - RightStart
              ComparedSegment = CurrentSeq( ComparedRow,GapEnd+1:GapEnd+RightEnd-RightStart+1 );
              if Test1
                 Ro = q
                 ComparedRow
                 Segment
                 ComparedSegment
              end
              if SimilarityBetweenTwoSeq( Segment,ComparedSegment ) > 0.7
                   Temp = CurrentSeq( q,RightStart:RightEnd );
                   CurrentSeq( q,RightStart:RightEnd ) = repmat('-',1,RightEnd-RightStart+1 );
                   CurrentSeq( q,GapEnd+1:GapEnd+(RightEnd-RightStart+1) ) = Temp;
                   if Test1
                      Cu = 'success_1'
                      Temp
                   end
              end
           elseif RightBound < NonGapEnd_Right
               Extra = RightEnd - RightStart - (RightBound - NonGapStart_Right);
               if Test2
                  Extra
                  S1 = CurrentSeq(q,RightEnd-Extra+1:RightEnd)
                  S2 = CurrentSeq(ComparedRow,RightBound+1:RightBound+Extra)
               end
               if RightBound+Extra<=Col && isequal( CurrentSeq(q,RightEnd-Extra+1:RightEnd),CurrentSeq(ComparedRow,RightBound+1:RightBound+Extra) )
                  for n = RightBound+1:Col
                      if n+Extra <= Col && CurrentSeq(q,n)~= CurrentSeq(ComparedRow,n+Extra)
                         break;
                      end
                  end
                  ComparedSegment = CurrentSeq( ComparedRow,GapEnd+1:GapEnd+RightEnd-RightStart+1 );
                  if Test2
                      n
                      Segment
                      ComparedSegment
                      SimilarityBetweenTwoSeq( Segment,ComparedSegment )
                  end                  
                  if SimilarityBetweenTwoSeq( Segment,ComparedSegment ) > 0.7
                     CurrentSeq( q,RightBound+1+Extra:n+Extra ) = CurrentSeq( q,RightBound+1:n );
                     Temp = CurrentSeq( q,RightStart:RightEnd );
                     CurrentSeq( q,RightStart:RightEnd ) = repmat('-',1,RightEnd-RightStart+1 );
                     CurrentSeq( q,GapEnd+1:GapEnd+(RightEnd-RightStart+1) ) = Temp;                  
                     if Test2
                        Cu = 'success_2'
                        Temp
                     end
                  end
               end
           end
        end
    end
    GapStart = GapEnd + find( OccupyCol(GapEnd+1:Col)==0,1 );
end
end

function Prob = CompareSeqDataNum( CurrentSeq,NucSeq,NonGap )
[ Row,Col ] = size( CurrentSeq );
Flag = zeros(1,Row);
parfor p = 1:Row
    if isequal(CurrentSeq(p,:),NucSeq)
       Flag(p) = 1;
    end
end
if nargin == 2
   Prob = sum(Flag)/Row;
elseif nargin == 3
    Gap = repmat('-',1,Col);
    GapFlag = ones(1,Row);
    parfor p = 1:Row
       if isequal(CurrentSeq(p,:),Gap)
          GapFlag( p ) = 0;
       end
    end
    if isempty(find(GapFlag==1,1))==0
        Prob = sum(Flag)/sum(GapFlag);
    else
        Prob = 0;
    end
end
end

function Similarity = SimilarityBetweenTwoSeq( Seq1,Seq2 )
% assume Seq1 and Seq2 take the same length.
if isempty(Seq1) || isempty(Seq2)
   Similarity = 0; return;
end
Len = min([length(Seq1), length(Seq2)]);
Similarity = sum( Seq1(1:Len)==Seq2(1:Len) )/max([length(Seq1), length(Seq2)]);
end

function CurrentSeq = CompressAlignment2Right( CurrentSeq,NoEmptyRow,OccupyCutoff )
MutationCutoff = 0.3;
if nargin == 2
   OccupyCutoff = 0.5;
end
Test = 0;
[ Row,Col ] = size( CurrentSeq );
%% move the columns from left 2 right
OccupyCol = SequenceColumnOccupation( CurrentSeq,OccupyCutoff,NoEmptyRow );
for p = 1:3:Col-2
    if sum( OccupyCol( p:p+2 ) ) == 0 &&( p<Col-2 && OccupyCol(p+3)==1 || p==Col-2)
       %% the gap region is found in between p:p+2
       LeftEnd = find( OccupyCol(1:p-1) == 1,1,'last' );
       if isempty(LeftEnd)==1 || LeftEnd <= 2, break; end
       if sum( OccupyCol(LeftEnd-2:LeftEnd) ) < 3, continue;  end
       
       [ ~,~,AAProbRight ] = Nucleotide2AA( CurrentSeq(:,p:p+2) );
       [ AA,~,AAProbLeft ] = Nucleotide2AA( CurrentSeq(:,LeftEnd-2:LeftEnd) );
       [ ~,OccupyNum ] = SequenceColumnOccupation( AA );
       MaxLeft = FindDominateAminoAcid( AAProbLeft,OccupyNum,50,MutationCutoff,0.5 );
       
       if Test
          GapRegion = p:p+2
          NonGap = LeftEnd-2:LeftEnd
          Mut1 = MaxLeft{1}
          Val = ismember( find(AAProbRight~=0),MaxLeft{1} )
       end
       if ismember( find(AAProbRight~=0),MaxLeft{1} )
          % it means that the rows in setdiff(1:Row,SelectedRow) can match better with region Right-2:Right
          if Test
             Go ='successfully move'
          end
          CurrentSeq = MoveSeqLeft2Right( CurrentSeq,LeftEnd-2,p+2 );
          OccupyCol( LeftEnd-2:LeftEnd ) = [0 0 0];
          if p > 6, p = p-6; end
       end
    end
end
end

function [SubSeq1,SubSeq2] = OptimalSeperatedSeq( CurrentSeq )
Test = 0; 
SubSeq1 = []; SubSeq2 = [];
[ Row,Col ] = size(CurrentSeq); 
MutationCutoff = 0.8; MutCol = floor(Col/3);
LongVector = zeros(1,MutCol); RowMajor = cell(1,MutCol);
parfor p = 1: MutCol
    if length( find(CurrentSeq(:,p)~='-') ) > 0.8*Row
       [ AAVector,Prob ] = Nucleotide2AA( CurrentSeq( :,3*(p-1)+1:3*p ) );
       if Test == 1
          Select = 3*(p-1) + 1:3*p
          Prob       
       end
       if sum( Prob<MutationCutoff ) == 26  % two major mutations,
          [ LongestValue,CollectRowMajor ] = FindMajorMutationSep( AAVector );
          LongVector( p ) = LongestValue;
          RowMajor{ p } = CollectRowMajor;
          if Test == 1
             p
             LongestValue
             CollectRowMajor
          end           
       end
    end
end

UniValue = unique(LongVector(find(LongVector~=0)));
Max = 0;
for p = 1:length( UniValue )
    Count = length(find(LongVector==UniValue(p)));
    if Count > Max
       Max = Count;
       MaxValue = UniValue(p);
    end
end
if Test == 1
   LongVector
   MaxValue
   Max
   Temp = 0.1*Col/3
   UniValue
end

if Max > 0.1*Col/3
   Local = find( LongVector == MaxValue,1 );
   SelectRow2 = RowMajor{ Local };
   SelectRow1 = setdiff(1:Row,RowMajor{ Local });
   SubSeq1 = CurrentSeq(SelectRow1,:);
   SubSeq2 = CurrentSeq(SelectRow2,:);
   ResultDir = './Full_Genome/AlignMafft/';
   WriteSequence2Fasta( SubSeq1,[],[ResultDir 'Test20.fasta'] );
   WriteSequence2Fasta( SubSeq2,[],[ResultDir 'Test21.fasta'] );  
   RealSeparate = 1;
end
end

function CurrentSeq = AlignGapInStopCodonRegion( CurrentSeq )
Test = 0;
[Row,Col] = size( CurrentSeq ) ; 
End = find( CurrentSeq(1,:)~='-',1,'last' ); Find = 0;
if isempty( find(CurrentSeq( 1,End-2:End ) =='-',1) ) && isequal( nt2aa( CurrentSeq( 1,End-2:End ),'ACGTOnly', false, 'AlternativeStartCodons', false),'*' ) 
   LeftEnd = find(CurrentSeq(1,1:End-3)~='-',1,'last');
   if LeftEnd == End-3,LeftEnd = End-2;end
   CurrentSeq = MoveSeqRight2Left( CurrentSeq,LeftEnd,Col );
   Find = 1;
end
if Test
    Find
    End
end

%% if the first sequence contains stop codon, then maybe also for others, however, remove double * at the end. only keep one.
StopNum = 0; Len = 11;
for p = 1:Row
    TempCol = find( CurrentSeq(p,:)~='-',1,'last' );
    if p>1 && Find == 1 && TempCol>=End+1 && isempty( find(CurrentSeq( p,End-2:End ) =='-',1) ) && isequal( nt2aa( CurrentSeq( p,End-2:End ),'ACGTOnly', false, 'AlternativeStartCodons', false),'*' ) 
       if Test 
          Remove = End+1:TempCol;
       end
       CurrentSeq( p, End+1:TempCol ) = repmat('-',1,TempCol-End);
       StopNum = StopNum + 1;
       continue;
    end
    for q = 0:3:Len
        if TempCol-q-2 >= 1
           LocalEnd = TempCol-q-3 + find( CurrentSeq(p,TempCol-q-2:TempCol-q)~='-',1,'last' );
           if isempty( find(CurrentSeq( p,LocalEnd-2:LocalEnd ) =='-',1) ) && isequal( nt2aa( CurrentSeq( p,LocalEnd-2:LocalEnd ),'ACGTOnly', false, 'AlternativeStartCodons', false),'*' ) 
               if LocalEnd + 3 <= TempCol
                  CurrentSeq( p, LocalEnd+1:LocalEnd+3 ) = '---';
               end
               if p == 1, End = LocalEnd; end
               StopNum = StopNum + 1;
           else
               break;
           end
        end
    end
   if StopNum<=3 && p > 10,break;end
end
end

function  [ Longest,CollectRowMajor,CollectRowRest,MajorInAA ] = FindMajorMutationSep( AAVector )
Test = 0;
Row = length( AAVector );
[ AAInt,AAProb ] = AminoAcid2Int(AAVector);
[~,MajorInAA] = max( AAProb );
CollectRowMajor = zeros(1,Row);
for p = 1:Row
    if p <= 2
       Neighbor = 1:3;
    elseif p >= Row-1
       Neighbor = Row-3:Row;
    else
       Neighbor = p-2:p+2;
    end
    if isempty( find(AAInt(Neighbor)==MajorInAA,1) ) == 0
       CollectRowMajor( p ) = 1;
    end
end
if Test == 1
   MajorAA = int2aa( MajorInAA )
   TempRow = find( CollectRowMajor==1 )
end
for p = 1:Row
    if CollectRowMajor( p )==1 && AAInt( p ) ~= MajorInAA
       CollectRowMajor( p ) = 0;
    else
       break;
    end
end
if Test == 1
   Temp1 = find( CollectRowMajor==0 )
end

%% CollectRow collect the rows that connect with each other.
p = find(CollectRowMajor==1,1);
while p<= Row
   if CollectRowMajor(p)==1
      End = p + find( CollectRowMajor(p+1:Row)==0,1 ) - 1;
      if isempty(End)==1,End=Row; end
      if End-p+1 <= Row*0.05
         CollectRowMajor( p:End ) = zeros(1,End-p+1);
      end
      p = End + 1;
   elseif CollectRowMajor( p ) == 0
      p = p + find( CollectRowMajor(p+1:Row)==1,1);
      if isempty(p),break; end
   end
end

if Test == 1
   Temp2 = find( CollectRowMajor==0 )
end

for p = 1:Row
    if CollectRowMajor( p )==1 && AAInt(p) ~= MajorInAA
       CollectRowMajor( p ) = 0;
    elseif CollectRowMajor( p )==1 && AAInt(p) == MajorInAA
        break;
    end
end
if Test == 1
   Temp3 = find( CollectRowMajor==0 )
end

Longest = 0; Count = 0; LongestPos = 0;
for p = 1:Row
    if CollectRowMajor( p ) == 1
       if Longest < Count
           Longest = Count;
           LongestPos = p;
       end
       Count = 0;
    else
        Count = Count +1;
    end
end

% CollectRowMajor = find(CollectRowMajor==1);
% CollectRowRest = setdiff(1:Row,CollectRowMajor);
CollectRowRest = LongestPos - Longest:LongestPos-1;
CollectRowMajor = setdiff(1:Row,CollectRowRest);
if Test == 1
    CollectRowMajor
    CollectRowRest
    Longest
end
end



function CurrentSeq = ShiftColumnRight2Left( CurrentSeq, StartCol,EndCol,MoveNum )
%% CurrentSeq: the original nucleotide sequences.
% StartCol : the start column
% EndColj : the end column
% MoveNum : the number of rows to move.
[Row,Col] = size(CurrentSeq);
if StartCol-MoveNum < 1 || EndCol>Col, return; end
SelectRow = [];
for p = StartCol-1:StartCol-MoveNum
    if isempty(SelectRow)
        SelectRow = find( CurrentSeq(:,p)~='-' );
    else
        SelectRow = union( SelectRow,find(CurrentSeq(:,p)~='-') );
    end
end
RestSelectRow = setdiff(1:Row,SelectRow);
if isempty(RestSelectRow),return;end
% adjust the selected row
TopSeq = CurrentSeq( RestSelectRow,StartCol:Col );
CurrentSeq( RestSelectRow,StartCol:Col ) = repmat('-',length(RestSelectRow),Col-StartCol+1);
CurrentSeq( RestSelectRow,StartCol-MoveNum:Col-MoveNum ) = TopSeq;

% adjust the non-selected row.
for p = EndCol+1:Col
    if isempty(find(CurrentSeq(SelectRow,p)~='-',1))
        TopSeq = CurrentSeq(SelectRow,p+1:Col);
        CurrentSeq(SelectRow,p+1:Col) = repmat('-',length(SelectRow),Col-p);
        CurrentSeq( SelectRow,p:Col-1) = TopSeq;
        break;
    end
end
end

function CurrentSeq = SeperateColumeBy2AA( CurrentSeq )
%% improve the region by shift if the columns' sequence (30% -- 50%) can be realigned.
[ Row,Col ] = size(CurrentSeq); MutationCutoff = 0.25;
p = 1; Test = 0;
while p <= Col-2
    if length(find(CurrentSeq(:,p)~='-')) > 0.8*Row
       [ AAVector,Prob ] = Nucleotide2AA( CurrentSeq( :,p:p+2 ) );
       if sum(Prob>MutationCutoff)==2 % two major mutations,
          if Test == 1
              p
              Col
          end
           [ Longest,CollectRowMajor,CollectRowRest,MajorInAA ] = FindMajorMutationSep( AAVector );
          if Longest > MutationCutoff * Row
             PartFirst = repmat('-',Row,3);
             PartFirst(CollectRowMajor,:) = CurrentSeq( CollectRowMajor,p:p+2 );
             PartSecond = repmat('-',Row,3);
             PartSecond(CollectRowRest,:) = CurrentSeq( CollectRowRest,p:p+2 );
             if Test == 1
                AA1 = CurrentSeq( CollectRowMajor(1),p:p+2 )
                AA2 = CurrentSeq( CollectRowMajor(1),p+3:p+5 )
                CompareRow = CollectRowMajor(1)
             end
             if isequal( CurrentSeq(CollectRowMajor(1),p:p+2),CurrentSeq(CollectRowMajor(1),p+3:p+5)) == 0
                %  [ FdRow,FdCol ] = size(CurrentSeq)
                %  Col
                %  p + 3
                CurrentSeq = [ CurrentSeq(:,1:p-1) PartFirst PartSecond CurrentSeq(:,p+3:Col) ];                                 
             else
                AAIntInMajorNum = 1; % how many columns share the same amino acids in major rows.
                for q = p+3:3:Col
                    % ComparedCol = q
                    [ ~,AAProb ] = AminoAcid2Int( Nucleotide2AA(CurrentSeq(CollectRowMajor,q:q+2)) );
                    [ ~,TempInAA ] = max( AAProb );
                    if TempInAA == MajorInAA
                       AAIntInMajorNum = AAIntInMajorNum + 1;
                    else
                        break;
                    end
                end
                AAIntInRestNum = 0; % how many columns share the same amino acids in minor rows.
                for q = p+3:3:Col
                    %ExamMinor = q
                    [ ~,AAProb ] = AminoAcid2Int( Nucleotide2AA(CurrentSeq(CollectRowRest,q:q+2) ));
                    [ ~,TempInAA ] = max( AAProb );
                    if TempInAA == MajorInAA
                       AAIntInRestNum = AAIntInRestNum + 1;
                       MinorEnd = q+2;
                    else
                        break;
                    end
                end
                %AAIntInMajorNum
                %AAIntInRestNum
                CurrentSeq = [ CurrentSeq(:,1:p-1) PartSecond PartFirst CurrentSeq(:,p+3:Col) ];
                if AAIntInMajorNum >= AAIntInRestNum && AAIntInRestNum > 0
                   CurrentSeq = MoveSeqRight2Left( CurrentSeq,p+3,MinorEnd+3,CollectRowRest );
                end
             end
             p = p + 3; Col = size(CurrentSeq,2);
          end
       end
    end
    p = p + 3;
end
end

function CurrentSeq = FitAlignmentByCompress( CurrentSeq,NoEmptyRow,OccupyCutoff )
% compress the column from seperated regions to fit 3 by 3 nucleotides.
Col = size( CurrentSeq,2 );
%% move the columns from right 2 left
OccupyCol = SequenceColumnOccupation( CurrentSeq,0.4,NoEmptyRow );
Test = 0; p = 1;
while p<=Col-2
    %OccupyCol = SequenceColumnOccupation( CurrentSeq,0.4,NoEmptyRow );
    if sum(OccupyCol( p:p+2 )) == 0 &&( p>1&& OccupyCol(p-1)==1 || p==1)
       Right = p+2 + find( OccupyCol(p+3:Col)==1,1 );
       if isempty( Right )==1 || Right + 2 > Col, break;end
       Left = find( OccupyCol( 1:p-1 )==0,1,'last') + 1;
       if isempty( Left ) == 1, Left = 1; end
       CutNum = 3 - mod( p - 1 - Left + 1,3);
       if Test == 1
          p 
          Right
          Left
          CutNum
       end
       if CutNum ~= 3
          % need to adjust the nucleotide column by fitting 3 * 3, either move the column to left, or to right.
          % first, consider move occupied column from right to left, if the
          % left side non-gap region is not 3 by 3 nucleotides.
          if Test == 1
             CutNum
             Right
             Left
          end 
          [ ~,ProbGapLeft ] = NucleotideProbability( CurrentSeq( :,p:p+CutNum-1 ) );
          GapLeftNucleotide = FindDominateAminoAcid( ProbGapLeft,OccupyCutoff );
          [ ~,ProbNonGapRight ] = NucleotideProbability( CurrentSeq(:,Right:Right+CutNum-1 ) );
          NonGapRightNucleotide = FindDominateAminoAcid( ProbNonGapRight,OccupyCutoff );
          ScoreLeft = 0;
          for n = 1:length( GapLeftNucleotide )
              if isempty( intersect(GapLeftNucleotide{n},NonGapRightNucleotide{n}) ) == 0, ScoreLeft = ScoreLeft + 1; end
          end
          
          ScoreLeft = ScoreLeft/length( GapLeftNucleotide );
          if ScoreLeft == 1 % the score is quite high, move nucleotide columns from right to left.
             if Test == 1
                 Go = 'Type 1'
             end
             CurrentSeq = MoveSeqRight2Left( CurrentSeq,p,Right+CutNum-1 );
             OccupyCol(p:Right+CutNum-1) = SequenceColumnOccupation( CurrentSeq(:,p:Right+CutNum-1),0.4,NoEmptyRow );
             %OccupyCol(p:p+CutNum-1) = ones(1,CutNum);
             %OccupyCol(Right:Right+CutNum-1) = zeros(1,CutNum);
             p = Right+CutNum;
             continue;
          end
          
          % score the procedure about move from right to left.
          CutNum = 3 - CutNum;
          [ ~,ProbNonGapLeft ] = NucleotideProbability( CurrentSeq(:,p-(CutNum-1):p ) );
          NonGapLeftNucleotide = FindDominateAminoAcid( ProbNonGapLeft,OccupyCutoff );
          [ ~,ProbGapRight ] = NucleotideProbability( CurrentSeq(:,Right-1-(CutNum-1):Right-1 ) );
          GapRightNucleotide = FindDominateAminoAcid( ProbGapRight,OccupyCutoff );
          ScoreRight = 0;
          for n = 1:length( GapRightNucleotide )
              if isempty( intersect(GapRightNucleotide{n},NonGapLeftNucleotide{n}) ) == 0, ScoreRight = ScoreRight + 1; end
          end
          ScoreRight = ScoreRight/length( GapRightNucleotide );
          if ScoreLeft >= ScoreRight % the score is quite high, move nucleotide columns from right to left.
             if Test == 1
                 Go = 'Type 2'
             end              
             CurrentSeq = MoveSeqRight2Left( CurrentSeq,p,Right+3-CutNum-1 );
             OccupyCol(p:Right+3-CutNum-1) = SequenceColumnOccupation( CurrentSeq(:,p:Right+3-CutNum-1),0.4,NoEmptyRow );
             %OccupyCol(p:p+3-CutNum-1) = ones(1,3-CutNum);
             %OccupyCol(Right:Right+3-CutNum-1) = zeros(1,3-CutNum);
             p = Right+3-CutNum;
          else % move nucleotide columns from left to right
             if Test == 1
                 Go = 'Type 3'
             end              
             CurrentSeq = MoveSeqLeft2Right( CurrentSeq,p-(CutNum-1),Right-1 );
             OccupyCol(p-(CutNum-1):Right-1) = SequenceColumnOccupation( CurrentSeq(:,p-(CutNum-1):Right-1),0.4,NoEmptyRow );
             %OccupyCol(p-(CutNum-1):p) = zeros( 1,CutNum );
             %OccupyCol(Right-1-(CutNum-1):Right-1) = ones( 1,CutNum );
             p = Right;
          end
       else
           p = p+1;
       end
    else
        p = p + 1;
    end
end
end


function CurrentSeq = AdjustNonGapRegion( CurrentSeq,NoEmptyRow )
%% removing columns in non gap regions which do not fit 3 amino acids
% ConsiderStopCodon = 1, consider the stop codon, so that the last part is not adjusted.
% ConsiderStopCodon = 1, all of gap regions are adjusted.
%% remove the column of reference sequence if all the rest are not aligned.
[ ~,Col ] = size( CurrentSeq );
[ OccupyCol,OccupyNum ] = SequenceColumnOccupation( CurrentSeq,0.8,NoEmptyRow,CurrentSeq(1,:) );
p = find( OccupyCol==1,1 ); Remove = [];
while p <= Col
    End = p+find( OccupyCol(p+1:Col)==0,1)-1;
    if isempty(End) == 1,break;end
    Num = 3 - mod(End-p+1,3);
    if Num ~= 3
       [ Value,Pos ] = sort( OccupyNum(p:End) );
       n = 1;
       while n <= Num
           if Value(n) < 0.1*NoEmptyRow
              Remove = [ Remove Pos( n ) + p-1 ];
           else
               break;
           end
           n = n + 1;
       end
       if n < Num
          Pos = End + find( OccupyCol(End+1:Col)==0,Num-n );
          CurrentSeq = MoveSeqRight2Left( CurrentSeq,End+1,Pos(length(Pos)) );
       end
    end
    p = End + find( OccupyCol(End+1:Col)==0,1 );
    if isempty( p )==0,break;end
end
if isempty(Remove)==0
   CurrentSeq = CurrentSeq(:,setdiff(1:Col,Remove));
end
end

function CurrentSeq = MoveGapRegionSeq2AlignedRegion( CurrentSeq,NoEmptyRow )
%% move those nucleotide sequences in gap region to aligned region which can improve the alignment
[Row,Col] = size(CurrentSeq);
OccupyCol = SequenceColumnOccupation( CurrentSeq,0.4,NoEmptyRow,CurrentSeq(1,:) );
[ TransMatrix, ProbMatrix ] = Nucleotide2Int( CurrentSeq );
Start = 0; 
for q = 4:Col-3
    if OccupyCol( q ) == 0 && Start == 0
       Start = q;
    elseif OccupyCol( q ) == 1 && Start ~= 0 
       End = q - 1;
       for p = 1:Row
         % if '---' appear in right side of the gap region
         if isequal( CurrentSeq( p,q:q+2),'---') && length( find(CurrentSeq( p,Start:End)~='-') ) >= 3
            LastPos = find( CurrentSeq(p,1:End)~='-',3,'last');
            LocalScore = 0.6*ProbMatrix( TransMatrix(p,LastPos(1)),q ) + 0.8*ProbMatrix( TransMatrix(p,LastPos(2)),q+1 ) + ProbMatrix( TransMatrix(p,LastPos(3)),q+2 );              
            if LocalScore > 1.2
                CurrentSeq( p,q:q+2 )= CurrentSeq( p,LastPos );
                CurrentSeq( p,LastPos ) = '---';
            end
         % if '--' appear in right side of the gap region
         elseif isequal( CurrentSeq( p,q:q+1),'--') && CurrentSeq( p,q+2) ~= '-' && length( find(CurrentSeq( p,Start:End)~='-') ) >= 2
            LastPos = find( CurrentSeq(p,1:End)~='-',2,'last');
            LocalScore = 0.8*ProbMatrix( TransMatrix(p,LastPos(1)),q ) + ProbMatrix( TransMatrix(p,LastPos(2)),q+1 );
            if LocalScore > 0.6
                CurrentSeq( p,q:q+1 )= CurrentSeq( p,LastPos );
                CurrentSeq( p,LastPos ) = '--';
            end
         % if '-' appear in right side of the gap region
         elseif CurrentSeq( p,q ) =='-' && CurrentSeq( p,q ) ~='-' && CurrentSeq( p,q ) ~='-' && length( find(CurrentSeq( p,Start:End)~='-') ) >= 2
            LastPos = find( CurrentSeq(p,1:End)~='-',1,'last');
            LocalScore = ProbMatrix( TransMatrix(p,LastPos),q );
            if LocalScore > 0.6
                CurrentSeq( p,q )= CurrentSeq( p,LastPos );
                CurrentSeq( p,LastPos ) = '-';
            end
         % if '---' appear in left side of the gap region
         elseif isequal( CurrentSeq( p,Start-3:Start-1),'---' ) && length( find(CurrentSeq( p,Start:End)~='-') ) >= 3
            LastPos = Start - 1 + find( CurrentSeq(p,Start:End)~='-',3,'first');
            LocalScore = 0.6*ProbMatrix( TransMatrix(p,LastPos(1)),Start-3 ) + 0.8*ProbMatrix( TransMatrix(p,LastPos(2)),Start-2 ) +...
                         ProbMatrix( TransMatrix(p,LastPos(3)),Start-1 );              
            if LocalScore > 1
                CurrentSeq( p,Start-3:Start-1 )= CurrentSeq( p,LastPos );
                CurrentSeq( p,LastPos ) = '---';
            end
         % if '--' appear in left side of the gap region
         elseif isequal( CurrentSeq( p,Start-2:Start-1),'--') && CurrentSeq( p,Start-3)~='-' && length( find(CurrentSeq( p,Start:End)~='-') ) >= 2
            LastPos = Start - 1 + find( CurrentSeq(p,Start:End)~='-',2,'first');
            LocalScore = 0.8*ProbMatrix( TransMatrix(p,LastPos(1)),Start-2 ) + ProbMatrix( TransMatrix(p,LastPos(2)),Start-1 );
            if LocalScore > 0.5
                CurrentSeq( p,Start-2:Start-1 )= CurrentSeq( p,LastPos );
                CurrentSeq( p,LastPos ) = '--';
            end
         % if '-' appear in left side of the gap region
         elseif CurrentSeq( p,Start-3)~='-' && CurrentSeq( p,Start-2)~='-' && CurrentSeq( p,Start-1)=='-' && length( find(CurrentSeq( p,Start:End)~='-') ) >= 1
            LastPos = Start - 1 + find( CurrentSeq(p,Start:End)~='-',1,'first');
            LocalScore = ProbMatrix( TransMatrix(p,LastPos),Start-1 );
            if LocalScore > 0.3
                CurrentSeq( p,Start-1 )= CurrentSeq( p,LastPos );
                CurrentSeq( p,LastPos ) = '-';
            end 
         end
       end
       Start = 0;
    end
end
end
