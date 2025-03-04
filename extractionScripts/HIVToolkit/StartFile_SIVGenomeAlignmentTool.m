% This example document provides our protocol for HIV full-length genome alignment
%
% THE APPROACH  
% ********************
% our protocol involves the following procedural steps: 

%   (0) perform the nucleotide sequence alignment on the original sequences
%       Many software are available for nucleotide alignment such as Muscle (up to 1000 sequences) or Mafft (more than 1000 sequences).
%       Users can use any of these software to perform nucleotide
%       alignment. The prepared multiple sequence alignments (MSAs) of nucleotide
%       genomes are used for the following procedure which optimize the
%       amino acid alignment in the full-length genome. Please keep the
%       HXB2 referene as the first genome sequence in the prepared MSA.

%   (1) Concatenate each HIV protein region from the aligned nucleotide genome
%       HIV proteins are translated from three open reading frames. We thus
%       concatenate individual HIV protein regions based on the HXB2
%       reference genome, which is the first sequence in the nuleotide MSA.
%       Each of the protein sequences are stored into one Fasta file (see
%       the directory called "TemporaryData"
        
%   (2) perform the codon sequence alignments
%       For each protein sequence file, we improve the sequence alignment
%       using the function TransferNucleotide2AminoAcidAlignment.m. In this
%       function, the nucleotide alignment will be automatically
%       transformed into codon sequence alignment. The quality of local
%       alignments will be improved based on the heuristic algorithm.

%   (3) Manual inspection of the sequence quality
%       One needs to manually inspect the sequence quality of the codon
%       sequence alignment. We also provide several functions to assist
%       manual alignment. For instance, MoveSeqLeft2Right can move a
%       defined codon positions from the left side to the right side (see 
%       example data in step (2).

%   (4) Transform the sequences from nucleotide forms to amino acid forms.
%       Each of the protein codon alignment files is transformed into an amino
%       acid fasta file.

%   (5) Assemble all amino acid protein sequences into full-length genomes.
%       When the amino acid alignments of HIV-1 proteins are prepared, we
%       assemble them into a full-length genome alignment. Each protein
%       region is considered as an independent region in the concatenated
%       genome sequences.
%
% *************************************************
% DEPENDENCIES 
% *********************** 
% Add Ons: Bioinformatics Toolbox; Parallel Computing Toolbox
%
% ***********************
% INPUTS
% ***********************
% Nucleotide full-length SIV genome sequences, aligned using MUSCLE
%
% ***********************
% OUTPUTS 
% *********************** 
% The aligned amino acid sequences of individual SIV proteins and the
% concatenated genome sequences
%
% *********************** 

% Copyright for original implementation: 
%             Guangdi Li (lowieli@gmail.com):  2012-2014
% 
% Permission is granted for anyone to copy, use, or modify this
% software and accompanying documents for any uncommercial
% purposes, provided this copyright notice is retained, and note is
% made of any changes that have been made. This software and
% documents are distributed without any warranty, express or
% implied. All use is entirely at the user's own risk.

% cd /home/lowie/HIVGenomeAnalysisToolbox
%currentFolder = pwd
%addpath( genpath( currentFolder ) )

%%---------------------------------------SIV full-length genome information ----------------------------------------%%
%% protein information of SIV genome, based on Mac239 reference strain (see website: https://www.hiv.lanl.gov/content/sequence/HIV/MAP/annotation.html )
LocPronRegion =        {  'P17'    'P24',  'P2', 'P7',         'P1','P6',   'Protease','P51', 'P31',    'Vif', 'Vpr',   'Tat', 'Rev', 'Vpx',  'GP120','GP41','Nef'  };
SIVProteinRegionName = { 'Matrix','Capsid','p2','Nucleocapsid','p1','p6',   'Protease','RT','Integrase','Vif', 'Vpr',   'Tat', 'Rev', 'Vpx',  'GP120','GP41','Nef'  };
ProteinNumber = length( SIVProteinRegionName );

% stop codon information 
StopCodonProtein = [      0        0       0         0          0     1       0          0        1      1       1       1        1      1       0      1     1    ];
% gene information, whether the translated protein has two genes in different open reading frames
Protein_Gene  =    [      1         1       1        1          1     1       1          1        1      1       1       2        2      1       1      1     1    ];

% Genome locations (relative to the aligned genomes), length unchanged
 %SIVProteinRegionName ={'Matrix','Capsid','p2','Nucleocapsid','p1',  'p6','Protease','RT','Integrase','Vif',  'Vpr', 'Tat',          'Rev',        'Vpx', 'GP120', 'GP41','Nef'};
 %LocationLowBound  =[    1053     1458    2145     2196       2352   2394   2555      2852   4529      5340    6151   6302    8806   6528   8806    5812    6670   8179   9077 ];
 %LocationUpperBound=[    1457     2144    2195     2351       2393   2585   2851      4528   5407      5984    6456   6597    8899   6597   9056    6150    8178   9240   9865 ]; 
 %ProteinLen       = [     135      229     17       52         14     64      99       559    293       215      102   131            108            113    503    354     264   ];
 LocationLowBound  =[     1238     1765    2467     2524       2684   2726   2888      3186   4529      5678    6498   6648    9226   6874   9226    6157    7016   8600   9497 ];
 LocationUpperBound=[     1764     2466    2523     2683       2725   2918   3185      4867   5407      6330    6803   6943    9323   6943   9480    6497    8599   9664   10305 ]; 
 ProteinLen       = [      169      229     17       52         14     64      99       559    293       215     102    131            108            113    503     354    264   ];

%%-----------------------------------HIV full-length genome information ----------------------------------------%%


%% the output fasta files of aligned protein sequences will be saved in the directory called: TemporaryData
OutputDir = './TemporaryData/'; 

%% step 1: concatenate each HIV-1 protein region in the full-length genome based on the reference strain M33262.1 (Mac239)
if 1 %% if you choose to run this step, set the value 1 here, otherwise 0.
CurrentStep = 'Step 1'
% make sure the first sequence of fasta input is the reference.
InputDir  = './alignment_SIV-nucleotide-complete.fasta' ; 

% extract out fasta sequences of each protein region, return sequences, sequence titles, sequence IDs
[ Seq,SeqTitle,SeqID ] = ExtractSequenceOut( InputDir );
    [ Row,Col ]= size( Seq );
    FullGenomeSequence = repmat( '-',Row,10000); Index = 0;
    FirstSeq = Seq(1,:); Pos = find(FirstSeq~='-');  No = 0;
    for p = 1:ProteinNumber
        IndexSet =[];
        for q = 0:Protein_Gene(p)-1 %% define the protein region "[Start - End]"
            Start = Pos( LocationLowBound(p+q+No) );  End = Pos( LocationUpperBound(p+q+No) );
            IndexSet = [IndexSet Start:End];
        end
        No = No + Protein_Gene(p) - 1;
        LocalSeq = Seq( :,IndexSet );
        if StopCodonProtein(p) == 0
           FullGenomeSequence( :,1+Index:Index+length(IndexSet) ) = LocalSeq ;
           Index = Index + length(IndexSet);
        else
           FullGenomeSequence( :,(1+Index):(Index+length(IndexSet)-3) ) = LocalSeq( :,1:(length(IndexSet)-3) ) ;
           Index = Index + length(IndexSet) - 3;
        end
        % write out the fasta files into the output directory called: "TemporaryData".
        [a,b] = WriteSequence2Fasta( LocalSeq,SeqTitle, [OutputDir SIVProteinRegionName{p} '.fasta' ] );
        fclose( 'all' );
    end
end



%% step 2: multiple sequence alignment of HIV-1 subtype B genome. 
% improve the alignment of nucleotide sequences by amino acid sequence alignment
if 1 %% if you choose to run this step, set the value 1 here, otherwise 0.
    CurrentStep = 'Step 2'
    for p = 1:17
        % for each protein, we transform nucleotide sequences to codon sequences
        LocalProtein = [ OutputDir SIVProteinRegionName{p} '.fasta' ]
        [ Seq,SeqTitle,SeqID ] = ExtractSequenceOut( LocalProtein );
        
        %% if there are many gaps and insertions in the sequences, the program may meet difficult to identify the ideal codon residues
        %% please check the temporary fasta outputs, named: "Tran1.fasta" to "Tran18.fasta"
        Seq = TransferNucleotide2AminoAcidAlignment( Seq,SeqTitle );
        
        % output the aligned sequences into the directory called "TemporaryData".
        [a,b] = WriteSequence2Fasta( Seq, SeqTitle, [ OutputDir SIVProteinRegionName{p} '_CodonAlignment.fasta' ] );
    end
    %d
end




%% step 3: Seaview manual examination: the aligned sequences are available in the fold: TemporaryData
if 0 %% if you choose to run this step, set the value 1 here, otherwise 0.
CurrentStep = 'Step 3'
% move sequences in a small region from left to right
[ CurrentSeq,SeqTitle,SeqID ] = ExtractSequenceOut( './TemporaryData/p6_CodonAlignment.fasta' );
CurrentSeq = MoveSeqLeft2Right( CurrentSeq,40,57 ); %% move the regions between the position 40 and 57 from left to right side

% move sequences in a small region from right to left
CurrentSeq = MoveSeqRight2Left( CurrentSeq,133,138 ); %% move the regions between the position 133 and 138 from right to left side

[ a,b ] = WriteSequence2Fasta( CurrentSeq, SeqTitle, './TemporaryData/p6_CodonAlignment_Improve.fasta'  );
fclose( 'all' );
end




%% step 4: transform the sequences from nucleotide forms to amino acid forms after the manual inspectation
if 1 %% if you choose to run this step, set the value 1 here, otherwise 0. Needs Step 2 or 3
 CurrentStep = 'Step 4'
 for p = 1:ProteinNumber
     %% transform the nucleotide to amino acid sequences 
     if fopen([ OutputDir SIVProteinRegionName{p} '_CodonAlignment_Improve.fasta' ],'r')>0
        LocalProtein = [ OutputDir HSIVProteinRegionName{p} '_CodonAlignment_Improve.fasta' ]
     else
        LocalProtein = [ OutputDir SIVProteinRegionName{p} '_CodonAlignment.fasta' ]
     end
     [ Seq,SeqTitle,SeqID ] = ExtractSequenceOut( LocalProtein );
     [a,b] = WriteSequence2Fasta( Nucleotide2AA(Seq), SeqTitle, [ OutputDir SIVProteinRegionName{p} '_AminoAcidAlignment.fasta' ] );
 end
end




%% step 5: assemble all amino acid protein sequences into full-length genome
%Needs Step 4
if 1
    CurrentStep = 'Step 5'
    [ Seq,SeqTitle,SeqID ] = ExtractSequenceOut( [ OutputDir SIVProteinRegionName{1} '_AminoAcidAlignment.fasta' ] );
    Row = 356; %added
    Col = size( Seq, 2 ); %updated
    FullGenome = repmat('-',Row,10000); Increase = Col;
    FullGenome(:,1:Col) = Seq(1:Row,:); %update
    for p = 1:ProteinNumber
        [ Seq, SeqTitle, SeqID ] = ExtractSequenceOut( [ OutputDir SIVProteinRegionName{ p } '_AminoAcidAlignment.fasta' ] );
        Col = size( Seq, 2 ); %updated
        FullGenome( :,(Increase+1):(Increase+Col) ) = Seq(1:Row,:); %update
        Increase = Increase + Col;
    end
    FullGenome = FullGenome( 1:Row , 1:Increase ); %update
    [a,b] = WriteSequence2Fasta( FullGenome, SeqTitle, [ OutputDir 'Mac239Genome_AminoAcidAlignment.fasta' ] ); 
end


