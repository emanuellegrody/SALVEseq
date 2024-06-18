% function MeasureGenomeGeneticDiversity()
% In this function, the genetic diversity will be calculated at each position of SIV genome
% Input: Full-length genome sequences (see example data in the directory called ExampleData: "A1_GenomeAA.fasta", "B_GenomeAA.fasta")
% Output: A list of genetic diversity values for all positions in SIV genome (see output data: "GeneticDiversity_BGenome.txt" );

%cd /MATLAB Drive
%currentFolder = pwd
%addpath( genpath( currentFolder ) )

%% the length of amino acid sequences

SIVProteinRegionName = { 'Matrix','Capsid','p2','Nucleocapsid','p1',  'p6', 'Protease', 'RT','Integrase','Vif', 'Vpr',   'Tat', 'Rev', 'Vpx',  'GP120','GP41','Nef'  };
SIVProteinNumber = length( SIVProteinRegionName );
SIVProteinLen       = [    169      229     17        52        14     64      99        559      293     215     102      131    108    113     503     354    264  ];
SIVStopCodonProtein     = [ 0        0       0         0         0     1       0          0        1       1       1        1      1      1       0       1      1   ];

%% step 1: calculate nucleotide genome diversity
for e = 1:1 %updated to 1:1
CurrentStep = 'Step 1'
% multiple sequence alignment of amino acid genomes in SIVmac239
[ GenomeSeq, SeqTitle, SeqID ] = ExtractSequenceOut( './alignment_SIV-nucleotide-complete.fasta' );
[ GenomeIntraDiversity, AveDiversity, ProteinDiversityScore ] = AnalysisGenomeIntraSubtypeDiversity( GenomeSeq, SIVProteinLen(1,1:17) );

%% output the postional genetic diversity for SIVmac239 genome:
fid = fopen('./GeneticDiversity_CompleteGenomes.txt', 'w' );
fprintf(fid, 'Average genetic diversity: %.4f\n\n', AveDiversity );
for p = 1:SIVProteinNumber
    for q = 1 : ( SIVProteinLen(p) - SIVStopCodonProtein(p) )
        fprintf(fid, '%s,%d,%.4f\n',SIVProteinRegionName{p},q,ProteinDiversityScore(p,q));
        endhe s
    fprintf(fid, '\n');
end
fclose('all');
end

%% step 2: calculate amino acid genome diversity
for e = 1:1 %updated to 1:1
CurrentStep = 'Step 2'
% multiple sequence alignment of amino acid genomes in SIVmac239
[ GenomeSeq, SeqTitle, SeqID ] = ExtractSequenceOut( './TemporaryData/Mac239Genome_AminoAcidAlignment.fasta' );
[ GenomeIntraDiversity, AveDiversity, ProteinDiversityScore ] = AnalysisGenomeIntraSubtypeDiversity( GenomeSeq, SIVProteinLen(1,1:17) );

%% output the postional genetic diversity for SIVmac239 genome:
fid = fopen('./GeneticDiversity_AminoAcids.txt', 'w' );
fprintf(fid, 'Average genetic diversity: %.4f\n\n', AveDiversity );
for p = 1:SIVProteinNumber
    for q = 1 : ( SIVProteinLen(p) - SIVStopCodonProtein(p) )
        fprintf(fid, '%s,%d,%.4f\n',SIVProteinRegionName{p},q,ProteinDiversityScore(p,q));
    end
    fprintf(fid, '\n');
end
fclose('all');
end