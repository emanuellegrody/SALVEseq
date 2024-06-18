function [ OccupyCol,OccupyNum,NonGapOccupyProb ] = SequenceColumnOccupation( CurrentSeq,Threshold,NoEmptyRow,Reference )
% this function test the columns of sequence is occupied, if the nongap is more than Threshold.
[ Row,Col ] = size( CurrentSeq );
OccupyCol = zeros(1,Col); OccupyNum = zeros(1,Col);
parfor p = 1:Col
    OccupyNum(p) = length( find(CurrentSeq(:,p) ~= '-' ) );
end
if nargin == 1
   NonGapOccupyProb = [];
elseif nargin == 2
   OccupyCol = OccupyNum >= Threshold*Row;
elseif nargin == 3
   OccupyCol = OccupyNum >= Threshold*NoEmptyRow;
elseif nargin == 4
   OccupyCol = OccupyNum >= Threshold*NoEmptyRow;
   parfor n = 1:Col
      if  Reference(n) ~= '-', OccupyCol( n ) = 1; end
   end
end

if nargout== 3
    NonGapOccupyProb = zeros( 1,Col );
    parfor p = 1:Col
        Temp = Row-length( find(CurrentSeq(:,p) == '-' ) );
        if Temp > 0
           NonGapOccupyProb( p ) = OccupyNum( p )/Temp;
        end
    end
end

end