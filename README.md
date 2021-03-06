# IRFinder
--------------
|  Introduction  |
--------------
This package will do the following steps:

Step 1. Run mugsy to get multiple sequence alignment results.

Step 2. Get scaffolds for each strain with some parameters.

Step 3. Do pairwise core-genome scaffold comparison between any two strains and find inverted repeats at the two ends of an inversion.

Step 4: Search the sequence of the inverted repeats in all the strains¡¯ genomes and then mark down its relative positions in all the strains¡¯ scaffolds.

The package works under Linux system.
 
After download and unzip the folder named SourceCode, you will find five sub-Folders and a README.txt in folder SourceCode.
-Folder Musgy contains code by Angiuoli SV and Salzberg SL. The code was downloaded from http://mugsy.sourceforge.net/
-Folder GRIMM-synteny and GRIMM contains code by Glenn Tesler. The code was downloaded from http://grimm.ucsd.edu/DIST/
-Folder blast-2.2.26 contains code by NCBI. The code was downloaded from ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.26/ 
-Folder scr contains code developed by me (Author: Dan WANG, danwang5-c@my.cityu.edu.hk).

--------------
|  EXAMPLE  |
--------------
Before running the code, please put all the five folders (Musgy,blast-2.2.26,scr,GRIMM-synteny and GRIMM) in your home directory.
Here is an example that you can run the following four steps separately.
All the inputs and outputs of each step are in the ~/src/Example directory.
Step 1 may take a long time. For example, for aligning 25 bacterial chromosomes, it will take almost 1 day. But for aligning the 3 bacterial chromosomes in our example, it will only take several minutes.
Step 2~Step 4 are fast and normally will only cost a few minutes.
---------------------------------------------------------------------------------------------------------
Step 1. Run mugsy to get multiple sequence alignment results:
% cd ~/Mugsy/mugsy_x86-64-v1r2.3
% mugsy --directory ~/scr/Example/1.MSA_results --prefix MSA_Result ~/scr/Example/Genomes/*.fna
where
	--directory ~/scr/Example/1.MSA_results: the output file to direcotry ~/scr/Example/1.MSA_results
	--prefix MSA_Result: the prefix of the output file is MSA_Result
	~/scr/Example/Genomes/*.fna: all the input chromosomes which are in the fna format
	(For the input format of mugsy, you can refer to the mugsy's website:)

This step outputs
  	~/scr/Example/1.MSA_results/MSA_Result.maf: multiple sequence alignment results

Step 2. Get scaffolds for each strain with some parameters
% cd ~/scr
% sh getScaffolds.sh 500 3000
where 
	500: is the minimum block size
	3000: is the maximum gap threshold
This step outputs
	~/scr/Example/2.Scaffolds/anchors_500_3000_c/mgr_macro.txt: the scaffolds for each strain. (Each integer in the scaffolds stands for a core-genome block whose 	length is larger than 500. And consecutive blocks in all the strains are merged into one block.)
	~/scr/Example/2.Scaffolds/anchors_500_3000_c/blocks.txt: keeps the position of each block in the mgr_macro.txt file
	~/scr/Example/2.Scaffolds/cordis.maf: the multiple sequence alignment result without the strain-specific segments
	~/scr/Example/2.Scaffolds/core_coords.txt: is the positions of core-genome blocks in each strain without filtering short blocks and merging consecutive blocks.
 
Step 3. Do pairwise core-genome scaffold comparison between any two strains and find inverted repeats at the two ends of an inversion:
#Firstly, merget the consecutive blocks in both the two strains;
#Secondly, delete the independent block-interchange and transposition;
#Thirdly, calculate the inversions between the two genomes by using grimm;
#Lastly, check whether a pair of inverted repeats exist at the two ends of an inversion and store the results in report.txt
% cd ~/scr
% sh reptA.sh ./Example/2.Scaffolds/anchors_500_3000_c/mgr_macro.txt ./Example/2.Scaffolds/anchors_500_3000_c/blocks.txt ./Example/2.Scaffolds/cordis.maf ./Example/Genomes ./Example/3.IRreport 20
where
	./Example/2.Scaffolds/anchors_500_3000_c/mgr_macro.txt: the scaffolds for each strain.
	./Example/2.Scaffolds/anchors_500_3000_c/blocks.txt: the position of each block in the ./Example/2.Scaffolds/anchors_500_3000_c/mgr_macro.txt file
	./Example/2.Scaffolds/cordis.maf: the multiple sequence alignment result without the strain-specific segments
	./Example/Genomes: the directory where input chromosome sequences are stored.
	./Example/3.IRreport: the output files to directory ./Example/3.IRreport
	20: The maximum inversion steps threshold. If the inversion steps between any two strains are larger than 20, we will skip this pair of strains.
	
This step outputs
	~/scr/Example/3.IRreport/reportA.txt: provides the inversion scenario between a pair of strains with inversion steps less than the threshold and whether a pair of Inverted Repeats (we use A/-A to represent a pair of IRs) exists at the two ends of an inversion. For example: The following paragraph means that to transform Strain 2 (source) to Strain 3 (destination), there are 7 inversions after deleting the independent transpositions and block interchanges between them. The first inversion is to inverse Block 23 trough Block 25 in strain 2 and no pair of inverted repeats (A/-A) are found at the two ends of this inversion in both source and destination's genomes. The third inversion is to inverse Block 7 trough Block 21 in strain 2. At the two ends of this inversion, a pair of IRs are found in Strain 2 (source), the length of this IR(A/-A) is 5001 bp and the similarity betWeen +A and -A is 99%. Also, a pair of IRs are found in Strain 3 (destination), the length of this IR(A/-A) is 5000 bp and the similarity betWeen +A and -A is 99%.
--------------------------------------------------------------------------------------------------------------------------------------------
	>Genome 2 to 3: 7 steps
	Step 1: 23 through 25 Reversal: Source:No A/-A Destination:No A/-A
	Step 2: -8 through -23 Reversal: Source:No A/-A Destination:No A/-A
	Step 3: 7 through 21 Reversal: Source:Found A/-A,Length = 5001,Similarity = 99% Destination:Found A/-A,Length = 5000,Similarity = 99%
	Step 4: -16 through 22 Reversal: Source:No A/-A Destination:No A/-A
	Step 5: 9 through 25 Reversal: Source:No A/-A Destination:No A/-A
	Step 6: -28 through 8 Reversal: Source:No A/-A Destination:No A/-A
	Step 7: -22 through 18 Reversal: Source:No A/-A Destination:No A/-A
---------------------------------------------------------------------------------------------------------------------------------------------

Step 4.Search the sequence of the inverted repeats in all the strains¡¯ genomes and then mark down its relative positions in all the strains¡¯ scaffolds.
% cd ~/scr
% sh plotScafA.sh A_Sq.fa faNamesPre.txt ./Example/Genomes A ./Example/4.PlotIR ./Example/2.Scaffolds/anchors_500_3000_c/blocks.txt 5000 0.5
where
	A_Sq.fa: the sequence of inverted repeats A/-A. It is in FASTA format.
	faNamesPre.txt: keeps the chromosome file names of each strain
	./Example/Genomes: the directory for storing our input chromosome files
	A: the name of the inverted repeats
	./Example/4.PlotIR: the output files to directory ./Example/4.PlotIR
	./Example/2.Scaffolds/anchors_500_3000_c/blocks.txt: the position of each block in ./Example/2.Scaffolds/anchors_500_3000_c/mgr_macro.txt file
	5000: the length of repeats A
	0.5: the threshold for short IR. It means that if the length of the found segment is less than 0.5 * IR's length, it will be defined as short IR and we add an s after the label. 
             (For example, if the IR name is A, then the short IR is labeled as As.)

This step outputs:
	~/scr/Example/4.PlotIR/scafMapWithA.txt: Similar to mgr_macro.txt file but the relative positions of IR are also shown in this file. 
	~/scr/Example/4.PlotIR/blks_A.txt: keeps the positions of blocks and inverted repeats in each strain.
	~/scr/Example/4.PlotIR/N_A_Sq.famap.txt: the blast result when comparing the sequence of inverted repeat A with the chromosome sequence of Strain N. In our example, N=1,2,3

