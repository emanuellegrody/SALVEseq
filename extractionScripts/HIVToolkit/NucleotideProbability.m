function [ Prob,ProbNonGap ] = NucleotideProbability( Nucleotide,Current )
% AAProb  shows the probability of amino acids in each column of AA.
% AAProbNonGap is the probability of amino acids without considering unknown amino acids

Symbol = 'ACGTURYKMSWBDHVN-*';
[Row,Col] = size(Nucleotide);
Len = length(Symbol);
if nargout == 1
   if nargin == 1
      Prob = zeros( Len,Col );
      for q = 1:Len
         parfor p = 1:Col
            Prob( q,p ) = length( find( Nucleotide(:,p) == Symbol(q) ) );
         end
      end
   elseif nargin == 2
      Prob = zeros( 1,Col );
      parfor p = 1:Col
         Prob( p ) = length( find( Nucleotide(:,p) == Current(p) ) );
      end
   end
   Prob = Prob/Row;
elseif nargout == 2
    Prob = [];
   if nargin == 1
      ProbNonGap = zeros( Len-2,Col ); % in this case, '-' '*' are not interested.
      for q = 1:Len-2
         parfor p = 1:Col
            if Row -length( find( Nucleotide(:,p) == '-') ) > 0
               ProbNonGap( q,p ) = length( find( Nucleotide(:,p) == Symbol(q) ) )/( Row-length( find( Nucleotide(:,p) == '-') ));
            end
         end
      end
   elseif nargin == 2
      ProbNonGap = zeros( 1,Col ); 
      parfor p = 1:Col
         ProbNonGap( p ) = length(find(Nucleotide(:,p)==Current(p)))/( Row-length( find( Nucleotide(:,p) == '-') ));
      end
   end    
end
end


