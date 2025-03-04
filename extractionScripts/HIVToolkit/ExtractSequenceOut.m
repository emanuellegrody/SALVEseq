function [ OrigSeqPop,SequenceTitle,SequenceID ] = ExtractSequenceOut( AADir,InformationType )
%% read sequences from sequence file in the format of .fas, .fasta, .aln , .msf, .phy
% if InformationType == 0, default, output all variables
% if InformationType == 1, only SequenceTitle and SequenceID

if nargin == 1, InformationType = 0; end
FileName = AADir;
if isempty( find( FileName == '/',1 ) ) == 0
   End = find( FileName == '/',1,'last' );
   FileName = AADir( (End+1):length(AADir));
elseif isempty( find( FileName == '\',1 ) ) == 0
   End = find( FileName == '\',1,'last' ); 
   FileName = AADir( (End+1):length(AADir));
end

if isempty( findstr( FileName, '.fas') ) == 0
    [ SequenceTitle, Sequences ] = fastaread( AADir );
elseif isempty( findstr(FileName, '.aln') ) == 0 || isempty( findstr(FileName, '.msf') ) == 0 || isempty( findstr(FileName, '.phy') ) == 0
    [ SequenceTitle, Sequences ] = multialignread( AADir );
end

SeqNum = length( SequenceTitle );
if InformationType == 1
   SequenceID = cell( 1,SeqNum );
   parfor p = 1:SeqNum
       Pos = find( SequenceTitle{p}=='/',1 );
       if isempty(Pos)==0
          SequenceID{ p } = SequenceTitle{p}(1:Pos-1);
       else
          Pos = find( SequenceTitle{p}==' ',1 );
          if isempty(Pos)== 1
             TempSeq = SequenceTitle{p};
             SequenceID{ p } = TempSeq( 1:length(TempSeq) );
          else
             SequenceID{ p } = SequenceTitle{p}(1:(Pos-1));
          end
       end
   end
   OrigSeqPop = [];
   return;
end

%% first, estimate the length of sequences.
MaxLen = 0; 
for p = 1:SeqNum;
    if MaxLen < length( Sequences{ p } ), MaxLen = length( Sequences{ p } ); end
end

%% second, store all sequences back into matrix.
OrigSeqPop = repmat( '-',SeqNum,MaxLen );
for p = 1:SeqNum, OrigSeqPop( p, 1:length( Sequences{ p })  ) = Sequences{ p }; end

%% change it into upper case sequences.
if OrigSeqPop( 1,find( OrigSeqPop(1,:)~='-',1 ) ) > 96
   OrigSeqPop = upper( OrigSeqPop );
end

%% if the id is required from output.
if nargout > 2
   SequenceID = cell( 1,SeqNum );
   for p = 1:SeqNum
       Pos = find( SequenceTitle{p}=='/',1 );
       if isempty(Pos)==0
          SequenceID{ p } = SequenceTitle{p}(1:Pos-1);
       else
          Pos = find( SequenceTitle{p}==' ',1 );
          if isempty(Pos)== 1
             TempSeq = SequenceTitle{p};
             SequenceID{ p } = TempSeq(1:length(TempSeq));
          else
             SequenceID{ p } = SequenceTitle{p}(1:Pos-1);
          end
       end
   end
end
end
