function [ AA,AAProb,AAProbNonGap ] = Nucleotide2AA( Nucleotide,CurrentAA )
% AAProb  shows the probability of amino acids in each column of AA.
% AAProbNonGap is the probability of amino acids without considering unknown amino acids
%s = matlabpool('size');
%if s == 0, matlabpool open; end
%I had to update this because matlabpool is no longer supported past Matlab
%v2014

s = gcp('nocreate'); % Check if a parallel pool already exists
if isempty(s)
    % Get the number of available physical CPU cores
    numCores = feature('numcores');
    % Create a parallel pool with the default number of workers (numCores)
    parpool('local', numCores);
end


AASymbol='ARNDCQEGHILKMFPSTWYVBZX*-?';
[Row,Col] = size(Nucleotide);
Len = floor(Col/3);
if Len*3 ~= Col
   Len = Len -1;
end
AA = repmat( '-',Row,Len );
parfor n = 1:Row
    AA(n,:) = ComputeLocalAA( Nucleotide(n,:),Len );
end
if nargout == 2
   if nargin == 1
      AAProb = zeros( 26,Len );
      for q = 1:26
         parfor p = 1:Len
            AAProb( q,p ) = length( find( AA(:,p) == AASymbol(q) ) );
         end
      end
   elseif nargin == 2
        AAProb = zeros( 1,Len );
       if length(CurrentAA) == 1
         parfor p = 1:Len
             AAProb( p ) = length( find( AA(:,p) == CurrentAA ) );
         end
      else % there are more than one amino acid.
         parfor p = 1:Len
             AAProb( p ) = length( find( AA(:,p) == CurrentAA(p) ) );
         end
      end
   end
   %AANum = AAProb;
   %AAProb = AAProb/Row;
elseif nargout == 3
    AAProb = [];
   if nargin == 1
      AAProbNonGap = zeros( 24,Len ); % in this case, '-' and '?' is not interested.
      for q = 1:24
         parfor p = 1:Len
            if Row -length( find( AA(:,p) == '-')) > 0
               AAProbNonGap( q,p ) = length( find( AA(:,p) == AASymbol(q) ) )/( Row-length( find( AA(:,p) == '-') ));
            end
         end
      end
   elseif nargin == 2
      AAProbNonGap = zeros( 1,Len ); 
      if length( CurrentAA )==1
         parfor p = 1:Len
            AAProbNonGap( p ) = length(find(AA(:,p)==CurrentAA))/( Row-length( find( AA(:,p) == '-') ));
         end
      else % there are more than one amino acid.
         parfor p = 1:Len
            AAProbNonGap( p ) = length(find(AA(:,p)==CurrentAA(p)))/( Row-length( find( AA(:,p) == '-') ));
         end
      end
   end    
else
    
end
end

function AA = ComputeLocalAA( Nucleotide,Len )
AA = repmat( '-',1,Len );
for p = 1:Len
    Temp= Nucleotide( (p-1)*3+1:p*3);
    if isequal(Temp,'') == 0 && isempty( find( Temp =='-',1) ) == 1
       AA(p) = nt2aa( Temp,'ACGTOnly', false, 'AlternativeStartCodons', false);
    end
end
end

function Num = ComparableNucleotide( AASet, CompareAA )
Row = size(AASet,1);
Num = zeros( 1,Row );
parfor p = 1:Num
    if isequal( AASet(p,:),CompareAA )
       Num(p) = 1;
    end
end
Num = sum(Num);
end
