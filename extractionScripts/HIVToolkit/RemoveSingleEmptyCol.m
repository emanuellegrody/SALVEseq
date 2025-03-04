function [ CurrentSeq,OccupyNum ] = RemoveSingleEmptyCol( CurrentSeq,ProbCutoff,Start,End,OccupyNum )
%% remove empty columns in gap region
Test = 0;
[ Row,Col ] = size( CurrentSeq );
if nargin == 1
   ProbCutoff = 0.1; Start = 1; End = Col;
   OccupyNum = zeros( 1,Col );
elseif nargin == 2
   Start = 1; End = Col;
   OccupyNum = zeros( 1,Col );
elseif nargin < 5 
   OccupyNum = ones( 1,Col );
   if Start > End
      return;
   end   
end

parfor p = Start:End, OccupyNum(p) = length( find(CurrentSeq(:,p) ~= '-' ) ); end
if Test
    Start
    End
   Count = OccupyNum( Start-2:End+2 )
   Count = OccupyNum( Start-2:End+2 )/Row
end
OccupyNum(Start:End) = OccupyNum(Start:End)/Row;

SelectedColumn = ones( 1,Col ); 
for p = End:-1:Start % unstricted remove. 
    if p>1 && p+1<=Col && OccupyNum(p-1) >= 2*OccupyNum(p) && OccupyNum(p) < ProbCutoff && OccupyNum(p+1) > 2*OccupyNum(p)
       SelectedColumn(p)=0;
       if Test
           Remove = p
       end
    elseif p+2 <= Col
       MaxProb = max(OccupyNum(p:p+1));
       if  p>1 && OccupyNum(p-1) > 2*MaxProb && OccupyNum(p) < ProbCutoff && OccupyNum(p+1) < ProbCutoff && OccupyNum(p+2) > 2*MaxProb
           SelectedColumn(p:p+1)=[0 0];
           if Test
              Remove = p:p+1
           end 
       end
    end
end

SelectedCol = find(SelectedColumn==1);
CurrentSeq = CurrentSeq(:,SelectedCol);
OccupyNum = OccupyNum(SelectedCol);

if Test
    Removed = find(SelectedColumn==0)
    WriteSequence2Fasta( CurrentSeq,[], '/home/lowie/HIV_DataSet/Full_Genome/AlignMafft/Short.fasta');    
end
end


function StrictRemove()

SelectedColumn = ones( 1,Col ); 
for p = End-1:-1:Start+1
    if OccupyNum(p-1) > 5*ProbCutoff && OccupyNum(p) < ProbCutoff && OccupyNum(p+1) > 5*ProbCutoff
       SelectedColumn(p)=0;
       if Test
           Remove = p
       end
    elseif p+2 <= Col && OccupyNum(p-1) > 5*ProbCutoff && OccupyNum(p) < ProbCutoff && OccupyNum(p+1) < ProbCutoff && OccupyNum(p+2) > 5*ProbCutoff
       SelectedColumn(p:p+1)=[0 0];
       if Test
           Remove = p:p+1
       end       
    end
end
CurrentSeq = CurrentSeq(:,find(SelectedColumn==1));
OccupyNum = OccupyNum(find(SelectedColumn==1));

if Test
    WriteSequence2Fasta( CurrentSeq,[], '/home/lowie/HIV_DataSet/Full_Genome/AlignMafft/Short.fasta');    
end

end
