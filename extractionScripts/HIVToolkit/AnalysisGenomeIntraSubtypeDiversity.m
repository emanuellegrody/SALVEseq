function  [ GenomeIntraDiversity,AveDiversity,ProteinDiversityScore,LocalCount,LocalTotal ] = AnalysisGenomeIntraSubtypeDiversity( GenomeSeq, AlignProteinLen, GagpCutoff )
%% if Option = 1, calculate the intra-subtype diversity according to the reference strain in the first one.
[ RowA,Col ] = size( GenomeSeq );
MapA = 1:Col;

%% do not compare the gap regions
 if nargin == 3
  for p = 1:Col
      if length( find( GenomeSeq(:,p)=='-' ) ) >= GagpCutoff*RowA 
          MapA(p) = 0;
      end
  end
 end

%% intra-subtype diversity
   LocalCount = -1*ones( 1,Col );  LocalTotal = -1*ones( 1,Col );
   GenomeIntraDiversity = -1*ones( 1,Col );
   parfor i = 1:Col
       [ Count,Total ] = LocalPairwiseIntra( GenomeSeq( 2:RowA,: ),MapA,i );
       LocalCount( i ) = Count; LocalTotal(i) = Total;
       if Total > 0,  GenomeIntraDiversity( i ) = Count/Total;end
   end
   
 %% option 1: calculation of average diversity: based on all difference.
 AveDiversity = sum(LocalCount)/sum(LocalTotal);
 
%% map to each protein based on reference strain in the first sequence
ProteinDiversityScore = [];
if nargout > 2
   % SumCount = 0;  SumTotal = 0; rato = zeros(1,10000);Increase = 0;
    RefA = find( GenomeSeq(1,:)~='-' );
    ProteinNum = length( AlignProteinLen );
    ProteinDiversityScore = -1*ones( ProteinNum,1000 );
    SumProteinLen = zeros( 1,ProteinNum );  MaxNum = 0;
    for p = 1:ProteinNum, SumProteinLen(p) = sum( AlignProteinLen(1:p) ) ;  end  
    for n = 1:length( RefA )
        if MapA( RefA(n) ) > 0
           LocalProtein = find( n <= SumProteinLen,1 );  Start = 0;
           if LocalProtein > 1,Start = SumProteinLen( LocalProtein-1 );end
           if LocalTotal( RefA(n) ) > 0
              ProteinDiversityScore( LocalProtein,n-Start ) = LocalCount( RefA(n) )/LocalTotal( RefA(n) );
              if 0
              SumCount = SumCount + LocalCount( RefA(n) );
              SumTotal = SumTotal + LocalTotal( RefA(n) );
              Increase = Increase + 1;
              rato(Increase) = LocalCount( RefA(n) )/LocalTotal( RefA(n) );
              end
           end
           if MaxNum < n - Start
              MaxNum = n - Start;
           end
        end
    end
    ProteinDiversityScore = ProteinDiversityScore(:,1:MaxNum);
    %trymean = mean(rato(1:Increase))
    
%% option 2: calculation of average diversity: based on all difference at reference strains
%   AveDiversity = SumCount/SumTotal
    for e = 1:0
              total = [];
              for l = 1:17
                  local = ProteinDiversityScore(l,:); local=local(find(local>-1));
                  total = [total local];
              end
              Mout = mean(total)
              AveDiversity
              d
    end
 end
end

function [ Count,Total ] = LocalPairwiseIntra( SeqA,MapA,i )
Count = 0 ; Total = 0;
RowA = size( SeqA,1);
  for p = 1:RowA
      for q = p+1:RowA
          if MapA(i) > 0 && SeqA( p, MapA(i) ) ~='-' && SeqA( q, MapA(i) ) ~='-' && SeqA( p, MapA(i) ) ~='X' && SeqA( q, MapA(i) ) ~='X'
             Total = Total + 1;
             if SeqA( p, MapA(i) ) ~= SeqA( q, MapA(i) )
                Count = Count + 1;
             end
          end
      end
  end    
end
