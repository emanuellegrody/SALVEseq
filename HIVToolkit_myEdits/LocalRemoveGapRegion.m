function CurrentSeq = LocalRemoveGapRegion( CurrentSeq,ProbCutoff )
% 1. remove unknown amino acids.
% 2. remove some single columns where few nucleotide appearing.
Test =0; 

if nargin == 1
   ProbCutoff = 0.2; 
end
DifferenceCutoff = 0.5;

CurrentSeq = SolveGapRegion( CurrentSeq,ProbCutoff,DifferenceCutoff );
if Test
    WriteSequence2Fasta( CurrentSeq,[], '/home/lowie/HIV_DataSet/Full_Genome/AlignMafft/Short1.fasta');    
    %Local = length(find(CurrentSeq(:,1290)~='-'))
end
if ValidateCurrentStructure( CurrentSeq,ProbCutoff )
    if Test
        Jump = 'Jump 1'
    end
    return; 
end


%% this function align all non-gap regions into 3 by 3 columns (without removing columns).    
[Row,Col] = size( CurrentSeq ); 
OccupyNum = zeros(1,Col);
parfor p = 1:Col, OccupyNum(p) = length( find(CurrentSeq(:,p) ~= '-' ) ); end
OccupyNum = OccupyNum/Row;

DifferenceCutoff1 = 0.5;
NonGapStart = find( OccupyNum > ProbCutoff ,1 );
while isempty( NonGapStart ) == 0 && NonGapStart <= Col % test the gap region.
    NonGapEnd = [];
    for p = NonGapStart+1:Col
        if Test && NonGapStart == 0
           p
           Temp1 = OccupyNum(NonGapStart)
           Term2 = OccupyNum(p)
           Temp3 = OccupyNum(p) > (1+DifferenceCutoff1)*OccupyNum(NonGapStart)
        end  
        %if abs( OccupyNum(p)-OccupyNum(p-1)) < 0.1, continue; end
        if OccupyNum(NonGapStart) > 0.5 && (OccupyNum(p) < ProbCutoff || OccupyNum(p) < DifferenceCutoff1*OccupyNum(NonGapStart) && OccupyNum(NonGapStart) > 0.7) ||...
           OccupyNum(NonGapStart) < 0.5 && OccupyNum(NonGapStart) >= 0.2 && (OccupyNum(p) > (1+DifferenceCutoff1)*OccupyNum(NonGapStart) || OccupyNum(p) <DifferenceCutoff1*OccupyNum(NonGapStart)) || ...
           OccupyNum(NonGapStart) < 0.2 && OccupyNum(NonGapStart) >= 0.1 && OccupyNum(p) >= 0.1&& abs( OccupyNum(NonGapStart)-OccupyNum(p))>0.1 ||...
           OccupyNum(p) <= 0.1
           NonGapEnd = p-1; break;
        end
    end
    
   % NonGapEnd = NonGapStart + find( OccupyNum(NonGapStart+1:Col)<DifferenceCutoff*OccupyNum(NonGapStart),1 );
   if isempty( NonGapEnd ), NonGapEnd = Col; end
    Len = 3 - mod( NonGapEnd-NonGapStart+1,3 );
    if Len ~= 3
       if Test
          Case = 'Type_0'
          NonGapStart
          NonGapEnd
          Len
       end  
       
        %% always move the short columns to the longer columns.
        %NextColumn = NonGapEnd + find( OccupyNum(NonGapEnd+1:Col)>DifferenceCutoff*OccupyNum(NonGapEnd),Len );
        Success = 0;
        AttempColumn = NonGapEnd + find( OccupyNum(NonGapEnd+1:Col)<ProbCutoff,1 );
        if isempty(AttempColumn)==0
           AttempColumn = AttempColumn + find( OccupyNum(AttempColumn+1:Col)>ProbCutoff,1 );
        end
        
        if Test && NonGapStart == 1501
            AttempColumn
           V_Test = isempty(AttempColumn) == 0 && AttempColumn>NonGapEnd+Len 
           V2 = OccupyNum(AttempColumn) < 0.5 && sum(OccupyNum(NonGapEnd+1:NonGapEnd+Len)>0.5)==0
        end
        if isempty(AttempColumn) == 0 && AttempColumn>NonGapEnd+Len && OccupyNum(AttempColumn) < 0.5 && sum(OccupyNum(NonGapEnd+1:NonGapEnd+Len)>0.5)==0
           Success = 1;
           for p = 1:Len
               if Test
                   
                  length( union( find(CurrentSeq(:,NonGapEnd+p)~='-'),find(CurrentSeq(:,AttempColumn+p-1)~='-'))) 
               end
               if length( union( find(CurrentSeq(:,NonGapEnd+p)~='-'),find(CurrentSeq(:,AttempColumn+p-1)~='-'))) < 0.5*Row
                  Success = 0;
               end
           end
        end
        if Test
            Final_Success = Success
            AttempColumn
            Len
            T = AttempColumn + 0:(Len-1)
        end
        if Success == 1
           NextColumn = AttempColumn + (0:(Len-1));
        else
           NextColumn = NonGapEnd + find( OccupyNum(NonGapEnd+1:Col)>DifferenceCutoff*OccupyNum(NonGapEnd),Len );
        end
        if Test
            NextColumn
        end
        if isempty(NextColumn)
            NonGapStart = NonGapEnd+Len+1;
            continue;
        end
       
       %% first, move from right to left.
       MaxLen = 0;
       for n = 1:Len
           Local = length(find(CurrentSeq(:,NextColumn(n))~='-'));
           if Local > MaxLen, MaxLen = Local; end
       end
       if length(find(CurrentSeq(:,NonGapEnd)~='-')) >= MaxLen
          %% pay attention, we might search less Len in NextColumn
          NextColumn = NextColumn(length(NextColumn));
          if Test
             Case = 'Type_1'
             NonGapStart
             NonGapEnd
             RMC = [ NonGapEnd+1,NextColumn ]
             L = OccupyNum(1501:1524)
          end
          
          CurrentSeq = MoveSeqRight2Left( CurrentSeq,NonGapEnd+1,NextColumn );
          [ CurrentSeq,OccupyNum ] = RemoveSingleEmptyCol( CurrentSeq,0.1,NonGapEnd+1,NextColumn,OccupyNum ); 
          
          ChangedColNum = Col - size( CurrentSeq,2 );
          [~,Col] = size( CurrentSeq );
          NonGapStart = NextColumn - ChangedColNum + find( OccupyNum( NextColumn+1-ChangedColNum:Col )>ProbCutoff,1 );
          if Test
             Case = 'Type_1: right 2 left'
             NextColumn
             StartSearch = NextColumn+1-ChangedColNum
             Update_NonGapStart = NonGapStart
             %if NonGapStart == 1501
                WriteSequence2Fasta( CurrentSeq,[], '/home/lowie/HIV_DataSet/Full_Genome/AlignMafft/Right1.fasta');
            %end
            % Case = 'close'
          end
          if OccupyNum( NonGapStart ) < 0.5 % if the sequence in this column is much less than 50%
             PreShortCol = find( OccupyNum(1:NonGapStart-1)>ProbCutoff,1,'last');
             if isempty( PreShortCol ) == 0 && PreShortCol > NonGapEnd
                CurrentSeq = MoveSeqLeft2Right( CurrentSeq,PreShortCol,NonGapStart);
                OccupyNum(NonGapStart) = length(find( CurrentSeq(:,NonGapStart)~='-' ))/Row;
%                WriteSequence2Fasta( CurrentSeq,[], '/home/lowie/HIV_DataSet/Full_Genome/AlignMafft/Right2.fasta');
             end
          end       
       else
           ForwardColumn = find( OccupyNum(1:NonGapEnd)>DifferenceCutoff*OccupyNum(NonGapEnd),3-Len,'last' );
           if length(ForwardColumn) < 3-Len, continue; end
           if Test && NonGapStart == 1300
              WriteSequence2Fasta( CurrentSeq,[], '/home/lowie/HIV_DataSet/Full_Genome/AlignMafft/Local1.fasta');
           end 
           if Test
               Case = 'Type_3'
               NonGapStart
               NonGapEnd
               RMC = [ForwardColumn(1),NextColumn(1)]
               Original_OccupyNum = OccupyNum( ForwardColumn(1)-3:NextColumn(1)+3 )
           end
           
           CurrentSeq = MoveSeqLeft2Right( CurrentSeq,ForwardColumn(1),NextColumn(1) );
           if Test && NonGapStart == 1300
              WriteSequence2Fasta( CurrentSeq,[], '/home/lowie/HIV_DataSet/Full_Genome/AlignMafft/Local2.fasta');
           end
           % save('/home/lowie/HIV_DataSet/Full_Genome/AlignMafft/Error.mat','CurrentSeq','OccupyNum','ForwardColumn','NextColumn');
           [ CurrentSeq,OccupyNum ] = RemoveSingleEmptyCol( CurrentSeq,0.1,ForwardColumn(1),NextColumn(1),OccupyNum );          
           if Test && NonGapStart == 1315
              WriteSequence2Fasta( CurrentSeq,[], '/home/lowie/HIV_DataSet/Full_Genome/AlignMafft/Local3.fasta');
              %Local = length( find(CurrentSeq(:,1288)~='-') )
              %Local = length( find(CurrentSeq(:,1289)~='-') )
           end           
           
           [ ~,Col ] = size( CurrentSeq );
           NonGapStart = ForwardColumn(1) - 1 + find( OccupyNum( ForwardColumn(1):Col ) > ProbCutoff,1 );
           if Test
              Case = 'Type_3: left 2 right'
              StartSearch = ForwardColumn(1)
             % ForwardColumn
             
              Update_NonGapStart = NonGapStart
             % N = OccupyNum(NonGapStart)
              Update_OccupyNum = OccupyNum( ForwardColumn(1):NextColumn(1)+3 )
            %  Case = 'close'
           end
           if Test && NonGapStart > 1279,
               
           end
       end
    else
       NonGapStart = NonGapEnd + find( OccupyNum( NonGapEnd+1:Col ) > ProbCutoff,1 );
    end
        
  %  WriteSequence2Fasta( CurrentSeq,[], '/home/lowie/HIV_DataSet/Full_Genome/AlignMafft/Last.fasta');        
   % if Test && isempty(NonGapStart)==0 && NonGapStart > 500,break;end
end

if Test
    WriteSequence2Fasta( CurrentSeq,[], '/home/lowie/HIV_DataSet/Full_Genome/AlignMafft/Short2.fasta');    
end

if ValidateCurrentStructure( CurrentSeq,ProbCutoff )
    if Test
        Jump = 'Jump 2'
    end
    return;
end

CurrentSeq = SolveGapRegion( CurrentSeq,ProbCutoff );
end


function Success = ValidateCurrentStructure( CurrentSeq,ProbCutoff )
% this function test whether all columns are fitting within amino acids
% according to the cutoff of ungap rows ProbCutoff.
Test = 0;
if nargin == 1, ProbCutoff = 0.2; end
Success = 1;
[Row,Col] = size( CurrentSeq ); 
OccupyCol = zeros(1,Col);
parfor p = 1:Col, OccupyCol(p) = length( find(CurrentSeq(:,p) ~= '-' ) ); end
OccupyCol = OccupyCol/Row > ProbCutoff;

Start = find( OccupyCol==1,1 );
while isempty(Start) == 0 && Start <= Col
    End = Start + find( OccupyCol(Start+1:Col)==0,1 ) - 1;
    if isempty( End ),End = Col;end
    if Test
       Start
       End
       Value =  ~(mod(Start,3)==1 && mod(End-Start+1,3)==0)
    end
    
    if ~(mod(Start,3)==1 && mod(End-Start+1,3)==0)
        Success = 0; break;
    end
    Start = End + find( OccupyCol(End+1:Col)==1,1 );
end
end


function CurrentSeq = SolveGapRegion( CurrentSeq,ProbCutoff,DifferenceCutoff )
%% process the gap region.
Test = 0; 

if nargin == 1
   ProbCutoff = 0.1; 
   DifferenceCutoff = 0.5;
elseif nargin == 2
   DifferenceCutoff = 0.3;
end


[ ~,Col ] = size( CurrentSeq );
%% remove completely empty columns in gap region.
SelectCol = ones(1,Col);
parfor p = 1:Col
    if isempty( find(CurrentSeq(:,p)~='-',1) ) == 1
       SelectCol( p ) = 0;
    end
end
CurrentSeq = CurrentSeq(:,find(SelectCol==1) );
if Test
    WriteSequence2Fasta( CurrentSeq,[], '/home/lowie/HIV_DataSet/Full_Genome/AlignMafft/GAp1.fasta');    
end

CurrentSeq = RemoveSingleEmptyCol( CurrentSeq );
if Test
    WriteSequence2Fasta( CurrentSeq,[], '/home/lowie/HIV_DataSet/Full_Genome/AlignMafft/GAp2.fasta');    
end


%% remove empty columns in gap region where 3 by 3 columns are satisfied.
[ Row,Col ] = size( CurrentSeq );
SelectedColumn = ones( 1,Col );
OccupyCol = zeros(1,Col);
parfor p = 1:Col, OccupyCol(p) = length( find(CurrentSeq(:,p) ~= '-' ) ); end

OccupyCol = OccupyCol/Row;
TempOccupyCol = OccupyCol > ProbCutoff;
while mod( sum(TempOccupyCol),3 ) ~= 0
    ProbCutoff = ProbCutoff + 0.02;
    TempOccupyCol = OccupyCol > ProbCutoff;
end

if Test
   Strange = find( OccupyCol > ProbCutoff == 1)
   %N = length(find( CurrentSeq(:,1327)~='-'))
end

OccupyCol = ModifyGapRegion( OccupyCol > ProbCutoff,OccupyCol,DifferenceCutoff );
if Test
    ProbCutoff
    EmptyCol = find(OccupyCol==1)
end

GapStart = find( OccupyCol==0,1 );
while isempty( GapStart ) == 0 && GapStart <= Col % test the gap region.
    GapEnd = GapStart + find( OccupyCol(GapStart+1:Col) == 1,1 ) - 1;
    if Test
       GapStart
       GapEnd
    end
    if isempty( GapEnd ), GapEnd = Col; end
    Len = mod( GapEnd-GapStart+1,3 );
    if Len ~= 0
       if Test
          Case = 'Type_1: remove columns'
          RMC = GapEnd-(Len-1):GapEnd
          Len
          
       end
       SelectedColumn( GapEnd-(Len-1):GapEnd ) = zeros(1,Len);
    end
    GapStart = GapEnd + find( OccupyCol(GapEnd+1:Col)==0,1 );
end
CurrentSeq = CurrentSeq( :,find(SelectedColumn==1) );
if Test
    WriteSequence2Fasta( CurrentSeq,[], '/home/lowie/HIV_DataSet/Full_Genome/AlignMafft/GAp3.fasta');
end


end


function OccupyCol = ModifyGapRegion( OccupyCol,OccupyNum,DifferenceCutoff )
Test = 0;
if Test
   Prepare = 'prepare the occupy vector'
   find( OccupyCol )
end
Col = length( OccupyCol ); Pos = Col;
if Col-find(OccupyNum > 0.8,1,'last') > 120, Pos = find(OccupyNum > 0.8,1,'last'); end
for p = 1:Col-1
    if p>=Pos
        Cutoff = 0.3;
    else
        Cutoff = 0.8;
    end
    if p~=Pos && p+1~=Pos && sum(OccupyCol(p:p+1))==2 && OccupyNum(p) > Cutoff && OccupyNum(p+1) < DifferenceCutoff*OccupyNum(p)
       OccupyCol(p+1)=0;
       if Test
          Type = 'Type 1'
          RmCol = p+1
       end
       for q = p+2:Col-1
           if OccupyCol(q)==1 && OccupyNum(q) < DifferenceCutoff*OccupyNum(p)
               OccupyCol(q)=0;
               if Test
                  RmCol = q
               end
           else
               break;
           end
       end
    elseif p~=Pos && p+1~=Pos && sum(OccupyCol(p:p+1))==2 && OccupyNum(p+1) > Cutoff && OccupyNum(p) < DifferenceCutoff*OccupyNum(p+1)
       OccupyCol(p)=0;
       if Test
          Type = 'Type 2'
          RmCol = p
       end       
       for q = p-1:-1:1
           if OccupyCol(q)==1 && OccupyNum(q) < DifferenceCutoff*OccupyNum(p+1)
               OccupyCol(q)=0;
               if Test
                  RmCol = q
               end               
           else
               break;
           end
       end
    elseif OccupyCol(p)==1 && OccupyCol(p+1)==0 && abs(OccupyNum(p+1) - OccupyNum(p)) < 0.1
        OccupyCol(p) = 0; 
       if Test
          Type = 'Type 3'
          RmCol = p
       end           
       for q = p-1:-1:1
           if OccupyCol(q)==1 && abs(OccupyNum(p+1) - OccupyNum(q)) < 0.1
               OccupyCol(q)=0;
               if Test
                  RmCol = q
               end               
           else
               break;
           end
       end
    end
end

end

