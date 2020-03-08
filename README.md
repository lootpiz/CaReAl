# CaReAl

CaReAl (Capturing Read Alignments) is a high-performance alignment capturing tool for visualizing the read-alignment status of nucleotide sequences and associated genome features.

 * Features
   - Visualizing full-depth of aligned reads.
   - Capturing read-level data rapidly and conveniently.
   - Optimized for the systematic exploration of regions of interest.
   - Useful for evaluating variant calls and detecting technical biases.

 * Prerequisite
   - Python (version >= 2.7) with packages _tabix, multiprocessing, optparse_
   - R (version >= 3.3) with library _graphics_
   - SAMTools (version >= 0.1.19)
   
 * Synopsis  
 ```./careal -b <BAM or BAMs in TXT> -t <target or targets in TXT> -r <reference genome> [OPTIONS]```
    * Required fields  
        -b, --bamfile <FILE> : input a BAM or BAMs with directory paths in TXT, (i.e., sample_01.bam, samples.txt)  
        -t, --target <FILE> : chromosome and position or chromosome and range, or targets in TXT, (i.e., chr10:53933206, chr10:53933106-53933306, targets.txt)  
        -r, --reference <FILE> : reference genome file in FASTA, (i.e., genome.fasta)  
    * Optional fields  
        -g, --gene <FILE> : gene information file in BED indexed by Tabix  
        -v, --with-vcf <STRING> : query VCF [TRUE or FALSE], default is FALSE  
        -o, --folder <STRING> : output directory name, default is "OUTPUT"  
        -n, --cpu <NUMBER> : number of CPUs to use, default is 4  

 * Troubleshoot
   - BAM file should be indexed with SAMTools and please use the command below to build an index file:  
   ```$ samtools index Sample_001.bam```
   - To make sure that tview in samtools works properly, try the command below:  
   ```$ samtools tview -d T -p chr1:2409792-2409992 Sample_001.bam hg19.fa```
   - Make sure that the VCF is located in the same folder as the BAM file and that their filenames are identical.

  * For more information please visit: [CaReAl](http://www.snubi.org/~lootpiz/CaReAl/)  

  ![CaReAl](/imgs/logo.png) 


