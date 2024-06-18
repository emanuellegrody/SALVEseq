function  [Ok,Empty] = WriteSequence2Fasta( TotalSequence,TotalSeqTitle,FileDir, SelectRow )
% length(TotalSeqTitle) must be larger than Row of TotalSequence

Ok = 0;
[Row,Col] = size(TotalSequence);
if nargin < 4
   SelectRow = 1:size( TotalSequence,1 );
end


%%%%%%%% write the title of sequences into fasta file.%%%%%%%%%%
%if Row < length(TotalSeqTitle)
%   Error ='The number of sequence title is less than total seqence ';
%end

if isempty( TotalSeqTitle ) == 1
   TotalSeqTitle = cell(1,Row);
   parfor p = 1:Row
       TotalSeqTitle{ p } = ['>' num2str( p )];
   end
else
   for p = 1:Row
       if TotalSeqTitle{ p }(1) ~= '>'
          TotalSeqTitle{ p } = ['>' TotalSeqTitle{ p }];
       end
   end    
end

NonGap = Col*ones( 1,Row );
if nargout == 0
   parfor p = 1:Row
       Index = find( TotalSequence(p,:) ~= '-',1,'last' );
       if isempty( Index ) == 0
          NonGap( p ) = Index;
       end
   end
end


% NonGap = max( NonGap );
Newfid = fopen( FileDir,'w' );
for p = 1:length(SelectRow)
    
    %Temp = SelectRow(p)
    Sequence = TotalSequence( SelectRow(p),1:NonGap( SelectRow(p) ) );
   % Sequence = TotalSequence( SelectRow(p),1:NonGap );
   
   %% remove empty sequenecs.
   if nargout == 2 && isempty(find( TotalSequence(p,:) ~= '-',1 ))==1
       continue;
   end
   
   fprintf( Newfid,'%s\n',TotalSeqTitle{SelectRow(p)} );
    Len = length( Sequence );
    while Len > 80
        fprintf( Newfid,'%s\n',Sequence(1:80) );
        Sequence = Sequence(81:Len);
        Len = Len - 80;
    end
    fprintf( Newfid,'%s\n',Sequence);    
end

if nargout == 2
    Empty = [];
end
fclose(Newfid);

end
