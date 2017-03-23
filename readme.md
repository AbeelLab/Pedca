***PEDCA Tutorial (Ploidy Estimation by Dynamic Coverage Analysis) ***

You can find an illustrated pdf copy of this tutorial with useful images at [*https://github.com/AbeelLab/Pedca/tree/master/Documents*](https://github.com/AbeelLab/Pedca/tree/master/Documents)

PEDCA is a ploidy estimation algorithm that infers copy number of the
contigs submitted as input based on the read coverage that aligns to
them. It requires as an input only an alignment file in .bam or .sam
format of a library or set of libraries aligned to a reference file of
the contigs that will be estimated.

***Pre-processing the data (5 steps)***

We need to align the reads against a reference.

*Step 1. Index your reference.*

Example using bwa (all command in one single line):

&lt;*path\_to\_bwa\_aligner*&gt;/bwa index -a bwtsw
*&lt;path\_to\_reference\_file/your\_reference.fasta&gt;*

*Step 2. Align your reads to your reference*

Example using bwa and paired end reads (all command in one single line):

*&lt;bwa\_aligner\_path&gt;**/*****bwa mem**
*&lt;path\_reference\_file/your\_reference.fasta&gt;*
*&lt;readsPath/readsPairEnd1.fasta&gt; &lt; readsPath
/readsPairEnd2.fasta&gt;* **&gt;** *&lt;destination\_folder*
/*example.sam*&gt;

*Step 3. You might want to transform your .sam file into a .bam format *

Example using samTools (all command in one single line):

*&lt;samToolsPath&gt;* **/samtools view -Sb** *&lt;destination\_folder*
/*example.sam*&gt; **&gt;** *&lt;destination\_folder* /*example.bam*&gt;

PEDCA just accepts one input file. If you have several libraries you can
put all your bam files in a folder (or create a folder with symbolic
links to all files you want to merge) and then:

*&lt;samToolsPath&gt;***/samtools merge**
&lt;*bam***\_***destination\_folder*/finalBamFile.bam&gt; **\*.bam**

*Step 4. Sort the .bam/.sam file *

Example sorting a .bam file using samTools (all command in one single
line):

*&lt;samToolsPath&gt;* **/samtools sort -o** *&lt;destination\_folder*
/sorted\_*example.bam*&gt; **-O bam -T**
*&lt;temp\_folderPath*/*tempName*&gt; *&lt;destination\_folder*
/*example.bam*&gt;

*Step 5. Index the sorted .bam/.sam file *

Example indexing a sorted bam file using samTools (all command in one
single line):

*&lt;samToolsPath&gt;* **/samtools index** *&lt;destination\_folder*
/sorted\_*example.bam*&gt;

***Downloading PEDCA***

[*https://github.com/AbeelLab/Ploest*](https://github.com/AbeelLab/Ploest)

***Using PEDCA***

PEDCA has been designed to require minimal parameterization. It works by
running a sliding window over the genome and measuring the average depth
of coverage inside each bin. Most of its parameters are dependant of the
window length and have default values that allows PEDCA to function on
contigs &gt; 500 bp and &lt; 2.000 Kbp. Nevertheless, because each
genome has its own particular characteristics it is possible to tune in
the rest of the parameters. Here is a list of those options and how the
influence the output.

At any moment you can obtain the following guide using: **java –jar
PEDCA.jar -help**

(You can also use ‘**–h**’ or ‘**help**’)

+++++++++++++++++++++++++++++++++++++++++++++++++++++++

PEDCA -help:

USAGE: java –jar PEDCA.jar -p &lt;project name&gt; -i &lt; input sam/bam
File&gt; -o &lt;output Folder&gt; &lt;&lt;OPTIONAL PARAMETERS&gt;&gt;

REQUIRED PARAMETERS:

-p (project name) – (String) Prefix used to generate the results file.

-i (input file) - (Pathway to .bam/.sam file) Pathway (absolute if
necessary) to the input file containing the alignment file. Must be a
.bam or .sam file

-o (output Folder) - (String) Pathway (absolute if necessary) to the
auto-generated output folder that will contain the results('./' if the
folder points to the current directory)

OPTIONAL PARAMETERS:

-m (multi Run) - (no parameters) Runs a preselected set of default
window lengths {500,750,1000,2000,3000}

-w (windows length) - (int) Length of the sliding window, to measure the
coverage inside contig. Default 500 bp

-c (coverage rate) - (int) Rate factor for the coverage sampling in the
Read count distribution. Default 100. The smaller it is, the less bins
are sampled

-k (mode smoother window) - (int) Number of points over which the ploidy
estimation is smoothed. The mode over k numbers of windows is used to
average the values of the bin. Default=49

-s (significant min) - (double) Threshold to consider a cluster peak in
the read count to be significant. Default 0.1

-b (fitter bin factor) - (double) Affects the number of bins used to FIT
the read count distribution. Default 2.5; Recommended between min=2.0
and max=4.0

-v (allele frequencies) - (Pathway to .vcf file) Pathway (absolute if
necessary) to the file containing the variant calling. Must be a .vcf
file

+++++++++++++++++++++++++++++++++++++++++++++++++++++++

**REQUIRED PARAMETERS:**

The first three arguments are required for PEDCA to function:

-p &lt;project name&gt; -i &lt; input sam/bam File&gt; -o &lt;output
Folder&gt;

PEDCA creates a folder named by the concatenation of the project name
and the size of the window length at the output pathway indicated by the
user. The output has the following structure:

./&lt;OutputFolderPath&gt;

./BaseCall

. BaseCallHistogramCluster\_1.jpg

> . BaseCallHistogramCluster\_2.jpg
>
> . Matrix1stCluster.vcf
>
> . Matrix2ndCluster.vcf

./&lt;Project Name&lt;*wl&gt;&gt;*

*./Ploidy\_Estimation\_Charts*

*.PEDCA*&lt;Project Name&lt;*wl&gt;&gt;PloidyEstimation.txt*

*.PEDCA*&lt;Project
Name&lt;*wl&gt;&gt;PloidyEstimation\_2nd\_Round\_.txt*

*.* *readsDistribution.jpg*

*.readsDistributionFittedFINALRESULT.jpg*
**OPTIONAL PARAMETERS:**

-w &lt;window length&gt;

Is worth noticing that the window length (*wl*), despite being the main
parameter in PEDCA, is not a mandatory field. If no other preference is
indicated, PEDCA runs with the default value of *wl*=500 bp. Even if the
parameter is not required, the advantage of using PEDCA is to have a
customizable window length so it is advised to use it with different
values and compare the results.

A short *wl* provides more coverage data points to estimate the ploidy,
you might want to shorten the *wl* if your ploidy estimation plot is too
discontinuous or if it doesn’t have much coverage information to support
a reliable estimation . On the other hand, you might want a larger *wl*
if your coverage/estimation plot looks overcrowded with coverage data
with too much variation, which leads to a fragmented discontinuous copy
number estimation . The minimum size of the window is 16 bp.

The *wl* also affects the sampling in the read count distribution, if it
is too big, it will lead to a irregular sampling with unrecognizable
clusters and false cluster ratios . If it is too small the read count
distribution will have its clusters merged together with long and thick
tails that might hide undetected peaks .

-v &lt;Pathway to .vcf file&gt;

If the option is selected and a .vcf file submitted, a folder named
BaseCall is also created at the same address., containing the allele
frequencies plots and a matrix with the positions and frequencies o each
base (order A,C,G,T).


-k &lt;mode smoother window&gt;

The coverage data has a certain degree of variation that we don’t want
to see reflected in the copy number estimation. In order to avoid
undesired jumps in the ploidy plot, the points are averaged by the mode
value over a bin of length *k*. If *k* is too small it might lead to
fragmented ploidy estimation in regions with noisy coverage (). The
continuity is smoothed with the default *k* value of 50 bp (). The
correct length of *k* depends on the required precision, and can be
parameterized. If *k* is too big, it might lead to the non detection of
regions with different ploidies (i.e. large structural variations found
in hybrid genomes)

-c &lt;coverage rate&gt;

The coverage rate is the definition with which the read count
distribution is drawn. It affects the number of bins in the plot. The
default value is 100. In some cases, when the plot is too irregular and
sawed, it is convenient to reduce this rate to have a fit that doesn’t
identify false peaks . If the sampling rate is too low, the bins might
merge clusters, so it’s not recommended to go below a value of 50.

-s &lt;significant min&gt;

When the read count is fitted, only peaks detected above a certain
threshold are taken into account, otherwise isolated peaks could be
considered a cluster. The default values is preset to s= 0.1 % of the
normalized number of windows for a given read count. For some genomes
this value might be too big and real clusters would be missed with a
potential misinterpretation of the correct cluster ratio. In the other
hand, if the distribution has a long tail with isolated values that are
not considered clusters, it is important to raise the threshold to
ignore false peaks in that region that would also jeopardize finding the
appropriate cluster ratio. In the example in many micro peaks are
detected in the long tail of the distribution. With a 0.2 significant
minimum all peaks below the red line are discarded. If instead the
default value was used, an error message would be displayed because the
peaks would not make sense:

+++++++++++ bestScore.candidateUnit: No CN mixture was able to satisfy
the constraints. Result == null

-m

This parameter enables multirun mode. Instead of running PEDCA with one
single *wl* value, it automatically runs it five times with the preset
values {500, 750, 1000, 2000, 3000} and output the respective results to
the output folder. These values work well for contigs larger than 500 bp
and up to 2.000 Kbp.

-b &lt;fitter bin factor&gt;

Default 2.5; We fit the read count distribution with 25 bins, that is
2.5 x the maximum number of ploidies that PEDCA can detect. That number
is adapted to detect a few clusters that are not very spread over the x
axis of the read count distribution. If the clusters are far away from
each other a higher number might better fit the distribution. It is
recommended to remain between the values min=2.0 and max=4.0