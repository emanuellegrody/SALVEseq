function CurrentSeq = Arrange1or2Nucleotide( CurrentSeq )
 %% step 3: arrange the 1 or 2, 3 nucleotide where the 3 nucleotides are not fit well.
 Test = 0;
 [ Row,Col ] = size(CurrentSeq);
 OccupyCol = SequenceColumnOccupation( CurrentSeq,0.8);
for p = 1:Row
    for q = 3:Col-3
        Find = 0; Len =0;
        if CurrentSeq(p,q)~= '-'&&isequal(CurrentSeq(p,q-1),'-') && isequal(CurrentSeq(p,q+1),'-')  % --A--
            Find = q; Len = 1;
        elseif CurrentSeq(p,q)~='-'&&CurrentSeq(p,q+1)~='-'&&isequal(CurrentSeq(p,q-1),'-')&&isequal(CurrentSeq(p,q+2),'-')  % --AA--
            Find = q; Len = 2;
        elseif isempty(find(CurrentSeq(p,q:q+2)=='-',1)) && isequal(CurrentSeq(p,q-1),'-') && isequal(CurrentSeq(p,q+3),'-')...
               && sum(OccupyCol(q:q+2)) < 3 % --AAA--
            Find = q; Len = 3;            
        end
        if Find == q
            % check the right side.
            RightStart = q + Len - 1 + find( CurrentSeq(p,q+Len:Col)~='-',1); 
            % check the left side.            
            LeftEnd = find( CurrentSeq(p,1:q-1) ~= '-',1,'last'); LeftStart = [];
            if isempty(LeftEnd)==0
               LeftStart = find( CurrentSeq(p,1:LeftEnd-1)=='-',1,'last')+1;
               if isempty(LeftStart),LeftStart = 1;end
            end
            if Test
               q
               RightStart
               LeftStart
               LeftEnd
            end
            Temp = CurrentSeq( p,q:q+Len-1 ); 
            if isempty(LeftEnd) == 0 && ( isempty(RightStart) == 1 || mod( LeftEnd-LeftStart+1,3 ) == 3 - Len)
                % retrain the movement of nucleotides by considering the probability.
              %  Prob = NucleotideProbability( CurrentSeq(:,LeftEnd+1:LeftEnd+Len),Temp );
                if Test
                    Move = 'to left'
                    Vector = LeftEnd+1:LeftEnd+Len
                end                
              %  if sum(Prob) > 0
                   
                   CurrentSeq( p,q:q+Len-1 ) = repmat('-',1,Len);
                   CurrentSeq( p,LeftEnd+1:LeftEnd+Len ) = Temp;
              %  end
            elseif isempty( RightStart ) == 0
               % if the right side of sequence exactly need 2 nucleotides to fit 3 by 3.
               % retrain the movement of nucleotides by considering the probability.
              %  Prob = NucleotideProbability( CurrentSeq(:,RightStart-Len:RightStart-1),Temp );
               if Test
                  Move = 'to right'
                  Vector = RightStart-Len:RightStart-1
               end                
             %   if sum(Prob) > 0
                   CurrentSeq( p,q:q+Len-1 ) = repmat('-',1,Len);
                   CurrentSeq( p,RightStart-Len:RightStart-1 ) = Temp;
             %   end
            end
        end
    end
end
end
