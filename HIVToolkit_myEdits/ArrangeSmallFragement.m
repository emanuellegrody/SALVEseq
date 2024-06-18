function CurrentSeq = ArrangeSmallFragement( CurrentSeq,NoEmptyRow,OccupyCutoff )
%% step 5: arrange the small fragement where the 3 nucleotides are not fit well.
Test = 0;
[ Row,Col ] = size( CurrentSeq );
[OccupyCol,OccupyNum] = SequenceColumnOccupation( CurrentSeq,OccupyCutoff,NoEmptyRow);%,CurrentSeq(1,:) );
for p = 1:Row 
    for q = 3:Col-3
        if CurrentSeq( p,q-1 ) == '-' && CurrentSeq( p,q ) ~= '-' %%% -AAAAA
           Start = q;
           End = q - 2 + find( CurrentSeq( p,q:Col ) =='-',1,'first');
           if isempty( End ) == 1, End=Col; end
            % if there is odd nucleotide region in the sequence
            % check the left side to see whether there is an odd region.
              LocalEnd = find( CurrentSeq(p,1:q-1)~='-',1,'last' );   % TTTTT----AAAA, check the position of TTTTT
              if isempty( LocalEnd ) == 0
                 LocalStart = find( CurrentSeq(p,1:LocalEnd-1)=='-',1,'last') + 1;
                 if isempty(LocalStart) == 1,LocalStart=1; end
                 if Test == 1
                    p
                    q
                    Start
                    End
                    LocalStart
                    LocalEnd
                    Cut1 = mod( LocalEnd-LocalStart+1,3 )
                    Cut2 = mod( End-Start+1,3 )
                    Off1 = mod( sum(OccupyCol(LocalStart:LocalEnd)),3 )
                    Off2 = mod( sum(OccupyCol(Start:End)),3 )
                    Val1 = End - Start + 1 - sum( OccupyCol(Start:End) )
                    Val2 = LocalEnd - LocalStart + 1 - sum( OccupyCol(LocalStart:LocalEnd) )
                    Val3 = mod( LocalEnd-LocalStart+1,3 ) == 0 && mod(LocalStart,3)~=0 && LocalStart-1>0 && mod( length(find(CurrentSeq(p,LocalStart-1)~='-')),3 ) ~= 0
                    if LocalStart-1 > 0
                       Val4 = length( find(CurrentSeq(p,1:LocalStart-1)~='-') )
                       Val5 =  mod( length(find(CurrentSeq(p,1:LocalStart-1)~='-')),3 ) ~= 0
                    end
                 end
                 if mod( End-Start+1,3 ) ~= 0 && (mod( LocalEnd-LocalStart+1,3 ) ~= 0 || mod( LocalEnd-LocalStart+1,3 ) == 0 && ...
                         mod(LocalStart,3)~=0 && LocalStart-1>0 && mod( length(find(CurrentSeq(p,1:LocalStart-1)~='-')),3 ) ~= 0)
                    % TTTT----AAAA or TTTTT----AAAA , the left side of AAAA is not the times of 3 nucleotides.
                    % check the amino acid probability and determine to move nucleotide from left side to right side, or right to left
                    
                    if mod( LocalEnd-LocalStart+1,3 ) == 0 && mod(LocalStart,3)~=0 && mod( length(find(CurrentSeq(p,1:LocalStart-1)~='-')),3 ) ~= 0
                        LeftGap = mod( length(find(CurrentSeq(p,1:LocalStart-1)~='-')),3 );
                    else
                        LeftGap = mod( LocalEnd-LocalStart+1,3 );
                    end
                    RightGap = mod( End-Start+1,3 );
                    if Test
                        Go  = 'Case 1'
                        Start
                        End
                        LocalStart
                        LocalEnd                        
                        LeftGap
                        RightGap
                        Lo1 = [ LocalEnd-LeftGap+1:LocalEnd,Start:Start+2-LeftGap ]
                        Replace1 = LocalEnd-LeftGap+1:LocalEnd+3-LeftGap
                        Replace2 = Start-3+RightGap:Start-1+RightGap
                        
                        Lo4 = [ LocalEnd+RightGap-2:LocalEnd,Start:Start+RightGap-1 ]
                    end
                    if isempty(find(CurrentSeq( p,[ LocalEnd-LeftGap+1:LocalEnd,Start:Start+2-LeftGap ] )=='-',1) )==0,continue;end                       
                    if p == 1 % the first sequence is considered as reference sequence, so let it formed directly.
                       TempSeq = CurrentSeq( p,[ LocalEnd-LeftGap+1:LocalEnd,Start:Start+2-LeftGap ] );
                       CurrentSeq( p,[ LocalEnd-LeftGap+1:LocalEnd,Start:Start+2-LeftGap ] ) = '---';
                       CurrentSeq( p,LocalEnd-LeftGap+1:LocalEnd+3-LeftGap) = TempSeq;
                       continue;
                    end
                    
                    Replace1 = LocalEnd-LeftGap+1:LocalEnd+3-LeftGap;
                    Replace2 = Start-3+RightGap:Start-1+RightGap;                  
                    if isequal( Replace1,Replace2),continue; end
                    
                    CurrentNucleotide = CurrentSeq( p,[ LocalEnd-LeftGap+1:LocalEnd,Start:Start+2-LeftGap ] );
                    CurrentAA = nt2aa(CurrentNucleotide,'ACGTOnly', false, 'AlternativeStartCodons', false);
                    [ ~,AAProb1 ] = Nucleotide2AA( CurrentSeq( :,Replace1 ),CurrentAA );
                    [ ~,AAProb2 ] = Nucleotide2AA( CurrentSeq( :,Replace2 ),CurrentAA );
                    [ ~,NucProb1 ] = NucleotideProbability( CurrentSeq( :,Replace1 ),CurrentNucleotide );
                    [ ~,NucProb2 ] = NucleotideProbability( CurrentSeq( :,Replace2 ),CurrentNucleotide );
                    SumNucProb = [ sum(NucProb1) sum(NucProb2) ]/3;
                    if Test
                        CurrentNucleotide
                        AAProb1
                        AAProb2
                        NucProb1
                        NucProb2
                        SumNucProb
                    end
                    if AAProb1 > AAProb2 || SumNucProb(1) > SumNucProb(2) && SumNucProb(1) > 0.7
                        % move the nucleotide from right to left side, if the probability is higher
                        LeftNucleotide = CurrentSeq( p,[ LocalEnd-LeftGap+1:LocalEnd,Start:Start+2-LeftGap ] );
                        CurrentSeq( p,[ LocalEnd-LeftGap+1:LocalEnd,Start:Start+2-LeftGap ] ) = repmat('-',1,3);
                        CurrentSeq( p,Replace1 ) = LeftNucleotide;
                        if OccupyNum( LocalEnd+3-LeftGap ) < 10
                           CurrentSeq( p,Start+2-LeftGap ) = CurrentSeq( p,LocalEnd+3-LeftGap );
                           CurrentSeq( p,LocalEnd+3-LeftGap ) = '-';
                        end
                        if Test == 1
                           Go = 'Case 1_1'
                        end
                    elseif AAProb2 > 0 || SumNucProb(1) < SumNucProb(2) && SumNucProb(2) > 0.1
                        % move the nucleotides from left side to right, if the probability is higher.
                        RightNucleotide = CurrentSeq( p,[ LocalEnd+RightGap-2:LocalEnd,Start:Start+RightGap-1 ] );
                        if isempty(find(RightNucleotide=='-',1)) && isequal('*',nt2aa(RightNucleotide,'ACGTOnly', false, 'AlternativeStartCodons', false)) == 0
                           CurrentSeq( p,[ LocalEnd+RightGap-2:LocalEnd,Start:Start+RightGap-1 ] ) = repmat('-',1,3);
                           CurrentSeq( p,Replace2 ) = RightNucleotide;
                           if Test == 1
                              Go = 'Case 1_2'
                           end
                        end
                    end
                 else
                     GapStart = find( CurrentSeq(p,1:LocalStart-1)~='-',1,'last' ) + 1;
                     if isempty( GapStart ) == 1
                         TempGapNum = 0;
                     else
                         TempGapNum = mod( LocalStart-GapStart,3 );
                     end
                     Value1 = 3 - mod( LocalEnd-LocalStart+1,3 );
                     if Test
                        Come = 'Else part'
                        Start
                        End
                        LocalStart
                        LocalEnd
                        GapStart
                        TempGapNum
                        Value1
                     end
                     
                     if Value1 == 3 && mod( End-Start+1,3 ) == 0 ...
                         &&  sum( OccupyCol( LocalEnd:Start ) ) == Start - LocalEnd + 1 && End + 3<Col && ...
                         CompareSeqDataNum( CurrentSeq(:,LocalEnd+1:LocalEnd+3),CurrentSeq(:,Start:Start+2)) > 0.1
                        % move  TTTTTT---AAA--- to TTTTAAA--- if part of AAA is in gap region
                         TempSeq = CurrentSeq( p,Start:End );
                         CurrentSeq( p,Start:End) = repmat('-',1,End-Start+1);
                         CurrentSeq( p,LocalEnd+1:LocalEnd+End-Start+1 ) = TempSeq;
                        if Test == 1
                           Go = 'Case 3_1'
                        end 
                     elseif Value1 ~= 3 && sum( OccupyCol(Start:Start+Value1-1) ) == 0 
                         % move to left.
                         TempSeq = CurrentSeq( p,Start:Start+Value1-1 );
                         CurrentSeq( p,Start:Start+Value1-1) = repmat('-',1,Value1);
                         CurrentSeq( p,LocalEnd+1:LocalEnd+Value1 ) = TempSeq;
                        if Test == 1
                           Go = 'Case 3_2'
                        end                          
                     elseif Value1 == 3 && mod( End-Start+1,3 ) == 0 && TempGapNum ~= 0 && sum( OccupyCol( LocalEnd+1:Start-1 ) ) > 0 ...
                            && sum( OccupyCol( LocalEnd+1:Start-1 ) ) < Start - (LocalEnd + 1)
                         TempSeq = CurrentSeq( p,LocalStart:LocalEnd );
                         CurrentSeq( p,LocalStart:LocalEnd ) = repmat('-',1,LocalEnd-LocalStart+1);
                         CurrentSeq( p,LocalStart-TempGapNum:LocalEnd-TempGapNum ) = TempSeq;                         
                        if Test == 1
                           Go = 'Case 3_3'
                        end                            
                     elseif  Value1 == 3 && mod( End-Start+1,3 )==1 && sum( OccupyCol( LocalEnd+1:Start-1 ) ) == Start - (LocalEnd + 1)
                         TempSeq = CurrentSeq( p,Start );
                         CurrentSeq( p,Start) = '-';
                         CurrentSeq( p,LocalEnd+1 ) = TempSeq;
                        if Test == 1
                           Go = 'Case 3_4'
                        end                                                   
                     end
                 end
              end
        end
    end
end
end