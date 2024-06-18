function CurrentSeq = AdjustLocal3By3Nucleotide(CurrentSeq,NoEmptyRow )
%% step 3: align the regions to left or right, if the occupy column is not 3 by 3,
% example:  ---AAA---BBB---    ---AAA---BBB---
%           --------AAA----    --------BBB----
%    ===>   ---AAA------       ---------BBB---
Test = 0;
[Row,Col] = size(CurrentSeq);
if nargin == 1, NoEmptyRow = Row; end
OccupyCol = SequenceColumnOccupation( CurrentSeq,0.4,NoEmptyRow );
% if amino acids are not aligned well with other column, with minor probablity.
for q = 5:Col-5
    if length( find(CurrentSeq(:,q)~='-')) < 0.1*NoEmptyRow
        for p = 1:Row
            if length( find(CurrentSeq(p,q:q+2)~='-')) == 3
                if q+3<=Col && sum( OccupyCol(q+1:q+3)) == 3 && CurrentSeq(p,q+3)=='-'
                   CurrentSeq( p,q+1:q+3) = CurrentSeq( p,q:q+2);
                   CurrentSeq( p,q)='-';
                elseif q+4<=Col && sum( OccupyCol(q+2:q+4)) == 3 && isequal( CurrentSeq(p,q+3:q+4),'--')
                   CurrentSeq( p,q+2:q+4) = CurrentSeq( p,q:q+2);
                   CurrentSeq( p,q:q+1)='--';
                elseif q-3 > 0 && sum( OccupyCol(q-3:q-1)) == 3 && CurrentSeq(p,q-3) == '-'
                   CurrentSeq( p,q-3:q-1) = CurrentSeq( p,q-2:q);
                   CurrentSeq( p,q)='-';
                elseif q-4 > 0 && sum( OccupyCol(q-4:q-2)) == 3 && isequal( CurrentSeq(p,q-4:q-3),'--')
                   CurrentSeq( p,q-4:q-2) = CurrentSeq( p,q-2:q);
                   CurrentSeq( p,q-1:q)='--';
                end
            end
        end
    end
end
if Test == 1
   ResultDir = '/home/lowie/HIV_DataSet/Full_Genome/AlignMafft/';
   WriteSequence2Fasta( CurrentSeq,[],[ResultDir 'Compress1.fasta']);
end

% the amino acid in non-gap regions are not formed 3 by 3.
for p = 1:Row
    if Test == 1
        p 
    end
    q = find( CurrentSeq(p,:)~='-',1 );
    while q<=Col
        End = q + find( CurrentSeq(p,q+1:Col)=='-',1 );
        if isempty(End)==1,End=Col+1; end
        Num = mod(q-1,3);
        if End - q == 3 && Num ~= 0 && isequal( CurrentSeq(p,q-Num:q-1),repmat('-',1,Num))
           if Test == 1
              Top = q-Num:q-1
              Org = q:q+2
              Loc = q-Num:q-Num+2
           end
           
           Temp = CurrentSeq(p,q:q+2);
           CurrentSeq(p,q:q+2) = '---';
           CurrentSeq(p,q-Num:q-Num+2) = Temp;
        end
        q = End + find( CurrentSeq(p,End+1:Col)~='-',1 );
        if isempty(q)==1,break;end
    end
end
end