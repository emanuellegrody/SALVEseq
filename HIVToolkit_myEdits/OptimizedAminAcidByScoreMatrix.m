function CurrentSeq = OptimizedAminAcidByScoreMatrix( CurrentSeq,NoEmptyRow,OccupyCutoff,ScoreMatrix )

Test = 0;
ProbCutOff = 0.3;
if nargin ==3,
    ScoreMatrix = blosum62;
end
[ Row,Col ] = size(CurrentSeq);
OccupyCol = SequenceColumnOccupation( CurrentSeq,OccupyCutoff,NoEmptyRow);
MajorAACol = zeros(1,Col); AAColProb = zeros(24,Col); 

for p = 1:Row
    q = find(CurrentSeq(p,:)~='-',1)-1;
    if isempty(q),continue;end
    q = q-mod(q,3)+1;
    while q < Col-3
        q = q + 3; q = q-mod(q,3)+1;
        %% move amino acids to left side.
        if q+2<=Col && CurrentSeq(p,q-1)=='-' && isempty(find(CurrentSeq(p,q:q+2)=='-',1)) && sum(OccupyCol(q:q+2))==3
            Start = find(CurrentSeq(p,1:q-1)~='-',1,'last')+1;
            if isempty(Start),continue;end
            Start = Start - 1 + find(OccupyCol(Start:Col)==1,1);
            if Start>=q || mod(Start,3)~=1 || sum(OccupyCol(Start:Start+2))~=3,continue; end
            if MajorAACol(Start)==0
               [ MajorAACol(Start),AAColProb(:,Start) ] = MajorAAInColumn( CurrentSeq(:,Start:Start+2),ProbCutOff );
            end
            if MajorAACol(q) == 0
               [ MajorAACol(q),AAColProb(:,q) ] = MajorAAInColumn( CurrentSeq(:,q:q+2),ProbCutOff );
            end
            if Test
               p
               Pos_Left = q:q+2
               Replace = Start:Start+2
               S1 = CurrentSeq( p,q:q+2)
               V1 = MajorAACol(Start)
            end
            
            KeyAAInt = aa2int(nt2aa(CurrentSeq( p,q:q+2),'ACGTOnly', false, 'AlternativeStartCodons', false));
            if MajorAACol(q) ~= KeyAAInt && ( MajorAACol(Start) > 0 &&  MajorAACol(q) > 0 && ScoreMatrix(KeyAAInt,MajorAACol(Start)) > ...
                ScoreMatrix(KeyAAInt,MajorAACol(q))||MajorAACol(Start)>0&&MajorAACol(q)<0&&ScoreMatrix(KeyAAInt,MajorAACol(Start)) >=0 ||...
                AAColProb(KeyAAInt,Start) > AAColProb(KeyAAInt,q) )
               if Test
                  Replace = 'to left'
               end
               Temp = CurrentSeq( p,q:q+2);
               CurrentSeq( p,q:q+2) = '---';
               CurrentSeq( p,Start:Start+2) = Temp;
               if isempty(q),break;end
               continue;
            end
        end

        %% move amino acids to right side
        if q+3 <=Col && CurrentSeq(p,q+3)=='-' && isempty(find(CurrentSeq(p,q:q+2)=='-',1)) && sum(OccupyCol(q:q+2))==3
            End = q + 2 + find(CurrentSeq(p,q+3:Col)~='-',1)-1;
            if isempty(End),continue;end
            End = find(OccupyCol(1:End)==1,1,'last')-2;
            
            if End<=p+2 || mod(End,3)~=1 || sum(OccupyCol(End:End+2))~=3,continue; end
            if MajorAACol(End)==0
               [ MajorAACol(End),AAColProb(:,End) ] = MajorAAInColumn( CurrentSeq(:,End:End+2),ProbCutOff );
            end
            if MajorAACol(q) == 0
               [ MajorAACol(q),AAColProb(:,q) ] = MajorAAInColumn( CurrentSeq(:,q:q+2),ProbCutOff );
            end
            if Test
               p
               Pos_Right = q:q+2
               Replace = End:End+2
               S1 = CurrentSeq( p,q:q+2)
            end            
            KeyAAInt = aa2int(nt2aa(CurrentSeq( p,q:q+2),'ACGTOnly', false, 'AlternativeStartCodons', false));
            if MajorAACol(q) ~= KeyAAInt && ( MajorAACol(End) > 0 &&  MajorAACol(q) > 0 && ScoreMatrix(KeyAAInt,MajorAACol(End)) > ...
                ScoreMatrix(KeyAAInt,MajorAACol(q))||MajorAACol(End)>0&&MajorAACol(q)<0&&ScoreMatrix(KeyAAInt,MajorAACol(End)) >=0 ||...
                AAColProb(KeyAAInt,End) > AAColProb(KeyAAInt,q) )
               if Test
                  Replace = 'to right'
               end
               
               Temp = CurrentSeq( p,q:q+2);
               CurrentSeq( p,q:q+2) = '---';
               CurrentSeq( p,End:End+2) = Temp;
               q = q-6;
            end
        end
       
    end
end
end

function [AAMajor,AAProbNonGap] = MajorAAInColumn( CurrentSeq,ProbCutOff )
% find the major amino acids  in each column, whose probabilities are
% larger than ProbCutoff.
[ ~,~,AAProbNonGap ] = Nucleotide2AA( CurrentSeq );
[~,AAMajor] = max(AAProbNonGap(1:24,:)); Find = 0;
if AAProbNonGap( AAMajor ) < ProbCutOff
   AAMajor = -1;
end
end