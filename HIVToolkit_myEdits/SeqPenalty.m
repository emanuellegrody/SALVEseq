function PenaltyScore = SeqPenalty( OrigSeq,UpdateSeq,Type )
% this function compare UpdataSeq with OrigSeq, and compute the penalty.
% OrigSeq and UpdateSeq must have the same dimension

if isequal( size( OrigSeq ), size( UpdateSeq ) ) ~= 1;
    Error = 'The dimension of OrigSeq is not equal to that of UpdateSeq.';
end

if nargin < 3, Type = 'Probability'; end
[ Row,Col ] = size( OrigSeq ); 
PenaltyScore = 0; 
Num = 0;
if isequal( Type,'Absolute' )
   for p = 1: Row
       for q = 1:Col
          if OrigSeq( p,q ) ~= '-'
             Num = Num + 1;
             if OrigSeq(p,q) == UpdateSeq(p,q)
                PenaltyScore = PenaltyScore + 1;
             end
          end        
       end
   end 
elseif isequal( Type , 'Probability' )
   OrigTran = Nucleotide2Int( OrigSeq );
   [ ~, Prob ] = Nucleotide2Int( UpdateSeq );
   for p = 1: Row
       for q = 1:Col
           if OrigSeq( p,q ) ~= '-'
              Num = Num + 1;
              PenaltyScore = PenaltyScore + Prob( OrigTran(p,q),q );           
           end        
       end
   end
end
% PenaltyScore
% Num
PenaltyScore = PenaltyScore/Num;

end