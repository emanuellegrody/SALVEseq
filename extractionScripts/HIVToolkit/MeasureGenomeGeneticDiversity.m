% function MeasureGenomeGeneticDiversity()
% In this function, the genetic diversity will be calculated at each position of HIV genome
% Input: Full-length genome sequences (see example data in the directory called ExampleData: "A1_GenomeAA.fasta", "B_GenomeAA.fasta")
% Output: A list of genetic diversity values for all positions in HIV genome (see output data in the directory called TemporaryData: "GeneticDiversity_InterGenome.txt", "GeneticDiversity_BGenome.txt" );

%cd /home/lowie/HIVGenomeAnalysisToolbox
%currentFolder = pwd
%addpath( genpath( currentFolder ) )

%% the length of amino acid sequences
HIVProteinRegionName = { 'Matrix','Capsid','p2','Nucleocapsid','p1','p6',      'Protease','RT','Integrase','Vif','Vpr',    'Tat', 'Rev', 'Vpu','GP120','GP41','Nef'  };
ProteinNumber = length( HIVProteinRegionName );
ProteinLen           = [  132       231    14        55         16   52          99        560     288      192    96       101    116    82     481    345    206   ];
StopCodonProtein     = [   0         0      0         0          0    1          0          0       1        1      1        1      1      1      0      1      1    ];

%% step 1: calculation of the intra-clade genome diversity
for e = 1:0
CurrentStep = 'Step 1'
% multiple sequence alignment of amino acid genomes in HIV-1 subtype B
[ B_GenomeSeq, B_SeqTitle, B_SeqID ] = ExtractSequenceOut( './ExampleData/B_GenomeAA.fasta' );
[ GenomeIntraDiversity, AveDiversity, ProteinDiversityScore ] = AnalysisGenomeIntraSubtypeDiversity( B_GenomeSeq, ProteinLen(1,1:17) );

%% output the postional genetic diversity for subtype B genome:
fid = fopen( './GeneticDiversity_BGenome.txt', 'w' );
fprintf(fid, 'Average genetic diversity: %.4f\n\n', AveDiversity );
for p = 1:ProteinNumber
    for q = 1 : ( ProteinLen(p) - StopCodonProtein(p) )
        fprintf(fid, '%s,%d,%.4f\n',HIVProteinRegionName{p},q,ProteinDiversityScore(p,q));
    end
    fprintf(fid, '\n');
end
fclose('all');
end

%% step 2: calculate amino acid genome diversity between HIV-1 subtype A1 and B genomes
for e = 1:1
CurrentStep = 'Step 2'
% multiple sequence alignment of amino acid genomes in HIV-1 subtype A1
[ A1_GenomeSeq,A1_SeqTitle,A1_SeqID ] = ExtractSequenceOut( '/ExampleData/A1_GenomeAA.fasta' );
% multiple sequence alignment of amino acid genomes in HIV-1 subtype B
[ B_GenomeSeq, B_SeqTitle, B_SeqID ] = ExtractSequenceOut( currentFolder '/ExampleData/B_GenomeAA.fasta' );

[ GenomeIntraDiversity, AveDiversity, PairwiseProteinDiversityScore ] = AnalysisGenomeInterSubtypeDiversity( A1_GenomeSeq,B_GenomeSeq, 1, ProteinLen(1,1:17) );

%% output the postional genetic diversity for subtype B genome:
fid = fopen( './TemporaryData/GeneticDiversity_InterGenome.txt', 'w' );
fprintf(fid, 'Average inter-clade genome diversity: %.4f\n\n', AveDiversity );
for p = 1:ProteinNumber
    for q = 1:ProteinLen(p) - StopCodonProtein(p)
        fprintf( fid, '%s,%d,%.4f\n',HIVProteinRegionName{p},q,PairwiseProteinDiversityScore(p,q) );
    end
    fprintf(fid, '\n');
end
fclose('all');
d
end