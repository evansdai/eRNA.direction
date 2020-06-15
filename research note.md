# C1 Cage project 

## Feb 24
### Download alignments
1. I downloaded data from https://fantom.gsc.riken.jp/zenbu/gLyphs/#config=NMT9yTLnH59gIVssI9WRfD;loc=hg19::chr12:20966910..20983875+
   - Normalized and filtered data:
      - [C1 CAGE] TGF-β Timecourse libraries normalised, QC filtered and grouped by timepoint
      - osc files, bed12 blocks
      - C1 CAGE TGF-beta timecourse libraries CAGE_4, CAGE_5, CAGE_6 (respectively libraries 1, 2 and 3 in Kouno et al., 2018), This track contains all the cells that passed quality control ( "Keep" = TRUE). They are grouped by time point  (t00, t06 or t24). Expression values are normalised in tags per million.
      - fixed_bin_size: 1 to get accurate TTS rather than bed12 blocks

   - Non-normalized data (raw counts):
     - [C1 CAGE] TGF-β Timecourse libraries. Non-filtered, non-normalised.
     - osc files, bed12 blocks
     - C1 CAGE TGF-beta timecourse libraries CAGE_4, CAGE_5, CAGE_6 (respectively libraries 1, 2 and 3 in Kouno et al., 2018). This track contains all the cells, which can be filtered with metadata such as "Keep" (TRUE/FALSE), "Timepoint" (t00, t06 or t24), etc. The chambers in these C1 runs were imaged at multiple focal planes, to detect cell doublets.
     - fixed_bin_size: 1 to get accurate TTS rather than bed12 blocks


2. I transcribe osc into bedGraphs format
   - osc2bedgraphs.py

3. I transcribe normalized&filtered bedGraphs into BigWig format
   - my_bedgraph2bw.sh
     - remove mitochondrial loci
     - sort bedgrahp

4. I make initial design matrix for bw files
  - bw2table.py

5. I extract CTSS using CageFightr
   ```
   ### CTSS summary ###
   Number of samples: 151
   Number of CTSSs: 1.357 millions
   Sparsity: 98.95 %
   Final object size: 2.46 GB
   ```

6. further analysis:
   - TODO why tagcount pm becomes integar in CTSSs variable?
   - download count bed file


## Feb 25, 2020
### redownload count data
1. Download filtered tagcount data 
      - [C1 CAGE] TGF-β Timecourse libraries normalised, QC filtered and grouped by timepoint
      - osc files, bed12 blocks
      - tagcount
      - C1 CAGE TGF-beta timecourse libraries CAGE_4, CAGE_5, CAGE_6 (respectively libraries 1, 2 and 3 in Kouno et al., 2018), This track contains all the cells that passed quality control ( "Keep" = TRUE). They are grouped by time point  (t00, t06 or t24). Expression values are normalised in tags per million.
      - fixed_bin_size: 1 to get accurate TTS rather than bed12 blocks

2. Redo the previous format transcription 2,3,4

3. Use CageFightr
  - quantifyCTSSs
    ```
    ### CTSS summary ###
    Number of samples: 151
    Number of CTSSs: 1.357 millions
    Sparsity: 98.53 %
    Final object size: 2.46 GB
    ```

  - calclulate  TPM
    - Matrix::t(Matrix::t(assay(object, 
        inputAssay))/(colData(object)[, totalTags]/1e+06))

  - quickTSS, quickEnhancer:
    not working on single cell data because most of single cells are unidirectional

  - pool cells with the same timepoints&run to make 9 bulk data
    use bulk data to define enhancers and TSSs
    - editted bw2table.py 
    - editted osc2bedgraphs
    - generate new bw files

  - Annotation
    - Many antisense transcritps.
      - Many of them are classified as enhancers
        - To avoid this 
          ```
           RSE <- combineClusters(object1=TSSs, object2 = Enhancers, removeIfOverlapping="object2")
          ```
          no use

  - Promoters(highTSS)
    - enrichment of GC islands?

  - How to estimate enhancers?
  - How to use estimated enhancers to construct a new SummarizedExperiment?
    - That can tell 3 groups of enhancers from each other


## Feb 26
   - Debug to learn functions in enhancers.R
     Use insertSource to insert browser() in source code in CAGEfightR
     ```
      insertSource("C:/Users/evans/Documents/albin_lab/single cell cage/scripts/CAGEfightR_1.6.0.tar/CAGEfightR_1.6.0/CAGEfightR/R/enhancers.R", 
                package = "CAGEfightR", functions = c("findBidirRegions","coverageWindows"))

     ```

    - Notice that GRanges have data of every base
     ```
        - # GRanges object with 641709 ranges and 1 metadata column:
         #            seqnames    ranges strand |              score
         #               <Rle> <IRanges>  <Rle> |          <numeric>
         #        [1]     chr1     16923      + | 0.0646536118093701
         #        [2]     chr1     16924      + |   13.9662760258504
         #        [3]     chr1     16925      + |   1.13871625184785
         #        [4]     chr1     17902      + | 0.0649472783972882
         #        [5]     chr1     17938      + | 0.0324736391986441
     ```
     - Remove 0 count lines in bedGraph files 
       - edit osc2bedgraphs
         ```    
         def batch_save_bedgraph(self):
           for strand, df in [['+', self.df_plus],['-',self.df_minus]]:
               for i in self.expression_name:
                   df_temp = df.loc[:, ['eedb:chrom','eedb:start.0base','eedb:end',i]]
                   # add this line
                   df_temp2 = df_temp.loc[df_temp[i]>0,:]
                   print('{}.{}'.format(i, strand))
                   df_temp2.to_csv('{}\{}.{}.bedGraph'.format(self.workfolder,i,strand), header = False, index = False, sep = '\t')
         ```

        - 
         ```
         ### CTSS summary ###
         Number of samples: 151
         Number of CTSSs: 1.357 millions
         Sparsity: 98.53 %
         Final object size: 41.6 MB

         > rowRanges(CTSSs)
         UnstitchedGPos object with 1356573 positions and 0 metadata columns:
                     seqnames       pos strand
                        <Rle> <integer>  <Rle>
                 [1]     chr1     11473      +
                 [2]     chr1     12173      +
                 [3]     chr1     12174      +
                 [4]     chr1     13119      +
                 [5]     chr1     13621      +
                 ...      ...       ...    ...
           [1356569]     chrY  59354808      -
           [1356570]     chrY  59355813      -
           [1356571]     chrY  59355814      -
           [1356572]     chrY  59359626      -
           [1356573]     chrY  59360555      -
           -------
           seqinfo: 93 sequences (1 circular) from hg19 genome


         ```

## Feb 27
- Build a function that extract Enhancers from pooled data and calculate ratio of bidirectional TSS based on single cell data
   - Use this to extract sequences
       - Enhancer %>% rowRanges %>% swapRanges() %>% promoters(upstream=50, downstream=50) %>% getSeq(bsg, .)

    - TODO Use this to classify enhancers
      (number_plus - number_minus) / (number_plus + number_minus)

- what about extract Enhancer regions from other database? 
  - https://bioconductor.org/packages/release/bioc/html/annotatr.html

- devide enhancers into groups:
  - pure composition of plus+minus unidirectional TSS
  - mix of bidirectional and unidirectional TSS
  - pure bidirectional TSS


## Feb 29
- Extract -100 ~ +100 sequences 
  - PCA of enhancers bidirectional expression features 
    - PC1: expression intensity
    - PC2: bidirectional expression
    - PC3: direction of unidirectional expression
  - Use 25% quantile of PCs to subset enhancers
  - Draw logo of each enhancer

- GC content drop in the direction of eRNA

- total expression score drop in single-cell-bidirection expression

- Filter enhancers with 100% bidirectional expression: (enhancer_features.csv)
  - 2050 total
  - 36 bidirectional expression
  - 

## Mar 3
- Finished slices "C:\Users\evans\Documents\albin_lab\Notes\Mar1_sequence_analysis.pptx"
  
- Idea: bidirectional enhancers ~ expression reaction rate
  - I can define cell-specific expression rate by the degree of outpace of pseudo-time than real-time
  - I can define enhancer-specific reaction rate by ...? Annotations?

- VCM: Variable chromatin module(What are they?)

- Malte told me that countScoredOverlaps recalc_thick_in_TCs 

## Mar 5
- How to train a classifier of enhancer based on single cell data

## Mar 6
Talked with Albin

- High noise in enhancer loci
  - validate by FANTOM 5

- Decide positive strand starting site by max peaks  

- Low GC% in bidirectional enhancers is solid
  - Looking into the detail: heatmap (enhancer by position, sorted by distance)

- PCA on enhancer features
  - Not enough motivation to do this

- Slices details:
  - what is the position


Some thoughts:
  Could bidirectional enhancers on two different alleles

Fantom enhancers are not necessarily bidirectional (C:\Users\evans\Documents\albin_lab\Notes\tracks\loosed_cutoff_validated)

Bidirectional enhancers are not necessarily Fantom ones

To Confirm:
  - bulk v.s. fantom: what proportion of the enhancers are fantom validated? ---> hope fantom sucks
  - single cell v.s. bulk
  - fantom "wide" annotations
  - fantom5 cage tss

algorithm:
  see single cell counts as a sample from bulk data: calculate the probability


## Mar 7
If most cells are unidirectional,
  why define minus left plus right
  can minus right plus left?

## Mar 8
Define enhancers based on intersect (bulk data, fantom label)
C:\Users\evans\Documents\albin_lab\Notes\tracks\bulk_validated

Looks like that single cell signals are sample from bulk
  How to reject this?
    Single cell signal ~ time series
    Bulk signal ~ time series

Define enhancer mid point by:
  bulk data

Define bidirectionality by:
  num of single cell having bidirectional transcripts / total cell


## Mar 9
- I found that I made a mistake in extracting -100~+100 sequences

- rewrite calc_bidirectproportion_cell_for_each_enhancers
  - calculate mid peak and mid avearage peak based on bulk data

- 644 enhancers that:
  - has single cell counts
  - was identified in bulk
  - validated in fantom

- bell shape GC% of all enhancers

- plot heatmap of enhancer GC%
  - distance = maxpeak_plus - maxpeak_minus
  - normalize each enhancer
  - "C:\Users\evans\Documents\albin_lab\Notes\Mar9_gc_heatmap_ordered_by_distance_midpeak.png"



- need to plot trend of gc content
  - trend?
    - loess regression ---> plot the y axis
    - hmm model ---> plot high gc region
    - GC% in overlapping sliding window (!)
    - normalize by gc of the whole enhancer


## Mar 10
- define singlecell enhancers in probabilistic scheme?

- sequence repeat in enhancer region?

- Divergence score

- visualize gc by heatmap 
- visualize expression by heatmap


- plot heatmap of enhancer GC%
  - define distance by single cell data?
  - positive control : ordered by gc% of each line 


- results:
  - plot heatmap using gc calculated by chunk

## Mar 13

- In 600+ enhancers identified in bulk data that are validated in fantom 
  - 63 has bidirectional expression
  - 221 was pure plus
  - 234 was pure minus
  - 126 was pure combination of plus and minus

- Use a wide definition of enhancers
  - http://bioinfo.vanderbilt.edu/AE/HACER/
  - 644 validated in fantom --> 799 validated in hacer

- Enhancer loci by bulk data: 4346
  - with single cell signal: 1679
    - validated in fantom: 644
    - validated in hacer: 799
  - without single cell signal: 2667 (What are they?)


## Mar 14

- Draw strand specific gc plots

- draw enhancers identified in pooled data with bulk CTSSs
  - It is a bad idea to identify enhancers based on pooled data

- Is there any difference between enhancers with and without single cell signal?
    - No in  GC content

- Sinlge cell Bidirectional enhancers have an average lower GC content

- align all peak sequences weighted by TPM

- Write a function calculating PSSM by strand


## Mar 15

- Are there any difference between bidirectional eRNA from unidirectional eRNA?
  - Sinle cell reflection of bulk?
    - 
  
- Are single cell bidirectional eRNA just sample from highly expressed bulk eRNAs?
  - Positive correlation between single cell expression and bulk expression

- Are single cell bidirectional eRNA  

- Write a function to calculate balanceBC for single cell data

- If use Enhancer_bulk, the data is so sparse that bulk and single cell data are very different

- Change strategy: use all fantom enhancers
  
## Mar 18 to Apr 8
- Change strategy: identify time-specific enhancers
  - t0: enhancers 1330 
        promoter	751	24.4		
        proximal	873	28.3		
        fiveUTR	37	1.2		
        threeUTR	38	1.2		
        CDS	35	1.1		
        exon	17	0.6		
        intron	778	25.3		
        antisense	0	0.0		
        intergenic	552	17.9
    	
  - t6: enhancers 1740 
        promoter	743	21.2		
        proximal	866	24.8		
        fiveUTR	41	1.2		
        threeUTR	37	1.1		
        CDS	48	1.4		
        exon	23	0.7		
        intron	983	28.1		
        antisense	0	0.0		
        intergenic	757	21.6	

  - t24: enhancers 1615 
        promoter	740	22.5		
        proximal	772	23.5		
        fiveUTR	40	1.2		
        threeUTR	45	1.4		
        CDS	45	1.4		
        exon	25	0.8		
        intron	923	28.1		
        antisense	0	0.0		
        intergenic	692	21.1	
      
    bulk: enhancers 4346 (at least supported by 3/9 bulk samples)
        promoter	818	11.9		
        proximal	1255	18.3		
        fiveUTR	68	1.0		
        threeUTR	176	2.6		
        CDS	148	2.2		
        exon	52	0.8		
        intron	2586	37.7		
        antisense	0	0.0		
        intergenic	1760	25.6	

    bulk: enhancers 834 (supported by all 9/9 bulk samples)
        promoter	725	30.3		
        proximal	728	30.4		
        fiveUTR	33	1.4		
        threeUTR	24	1.0		
        CDS	32	1.3		
        exon	16	0.7		
        intron	480	20.1		
        antisense	0	0.0		
        intergenic	354	14.8	

- Subset by with v.s. without single cell data
  - t0: 356 v.s. 974 (26.77%) 
  - t6: 409 v.s. 1331 (23.51%)
  - t24: 596 v.s. 1019 (36.90%) -------> because of increasing cell number?
                              --------> venn plots(Apr9) very time specific
  
- Validation with fantom ----------> increasing
  - t0: 624 / 1330  (46.92%)
  - t6: 826 / 1740 (47.47%)
  - t24: 793 / 1615 (49.10%)

- Validation with HACER ----------> increasing
  - t0: 738 / 1330 (55.49%)
  - t6: 995 / 1740 (57.18%)
  - t24: 956 / 1615 (59.20%)


- GC% (plots Apr8)
  - t0 : higher GC% -------> house holding ?

- Subset by single cell bidirection  
  - bidirectional: percent TPM bidirectional expression > 90%
  - plus: 100% plus
  - minus: 100% minus
  - others: the rest
  - bidirectional, minus, plus, others
  - t0:
    - 10, 152, 132, 62
  - t6:
    - 11, 174, 149, 75
  - t24:
    - 25, 235, 206, 130
  - chi-square test not significant
  

## Apr 8
Idea: check whether transcription happens upstream or downstream of a unidirectional eRNA 
- On pooled level: 
- On sc level:

- Direction makes sense:


## Apr 12

Look for eRNA context on Single cell level
- TODO
[-] Extract both upstream and downstream nearest promoters
[ ] validate promoters by public database
[-] see if upstream:plus, downstream: minus is specific for active enhancers
    - compare single cell active enhancer and non active enhancers

## Apr 13
Observed TSSs having larger average scores near active enhancer
- is this an artifact caused by averaging ? few cell has large TSSs, divided by large cell numbers
  - poisson process
    - are TSSs conserved among groups?

TODO
[ ] add features into extract_... function
  [ ] enhancer cell number (support)
  [ ] promoter support number (support2)
      average expression of promoter should = score / support2, rather than score / support 1
  [ ] num of TSSs on two sides

## Apr 14
[-] Extract single cell feature of each enhancer

PCA
- PC2 seperates minus and plus
  - [ ] whether time point conserved
- PC1 sepeartes?

How to inpromve separation?
- More gradient in thr_expr?
  feature engineering?

## Apr 16
Talk with Albin
- Albin says that there is a correlation between transcription and GC% rich
  - we can test if number of cage tag ~ GC% rich region

- Identify enhancers with strong burstness
  - using mean v.s. median to identify strong burst enhancer on single cell level 
  - are strong burst enhancers unidirectional? bidirectional?
  - Albin said we should still see the spread of position even if enhancers are transcriptional bursting

- More sequence analysis on minus / plus
  - look at genome browser for things I dont understand

- Transcriptional context upstream and downstream 
  - collect cage tag cumulatively up and down ward
  - where to stop? experiment

## Apr 18
- sc bidirectional enhancers are likely to be AAAA rich (by naked eye, genome browser)
  
- Define CpG islands:
  - https://www.bioinformatics.org/sms2/cpg_islands.html

  - https://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=824733819_lDxUW7TA3KjeCQFMhsEnFVorYKB1&g=cpgIslandSuper

  -  Obs/Exp CpG = Number of CpG * N / (Number of C * Number of G)
  -  where N = length of sequence.

- plus is likely to have CpG island on the left side
- minus is likely to have CpG island on the right side

- to test:
  - for each position, 1 if -100~+100 is a GC island, 0 if not
  - calculate the mean weighted position: is it <0 or >0?
    - should observe <0 for plus and >0 for minus

- CpG island test:
  are they distributed on diff sides?

## Apr 19

- word freq analysis
  - then pca

- to let minus transcripts makes sense, generate reverse sequence

- to compare up and down stream, split at 0 and extract feature separately

- to test
  - are bidirect CT rich?

- draw logo of TSS for each cell

- For each word, calculate the distribution along genome and try to find a word that distinguish directionality



## Apr 21

- enhancer
  - Active mammalian enhancers can initiate transcrption at both boundaries [1]. However few cells contain stable eRNA CAGE tags of both strands [2].
  - enhancers are half opened?
  - site is over-reused
  - reuse site is precise (do not change)

- suggest first plus, second minus can also be enhancers

- By selecting bursting enhancer and bursting cell, what do we capture
  - a solid bound between TF and enhancer

- look at the minus-keep-minus TSSs, are they really bidirectional on bulk data?

- Albin says that there is a correlation between transcription and GC% rich
  - we can test if number of cage tag ~ GC% rich region
  - on single enhancer level 
  - GC% normalized by G% and C%?


- repeatative elements: how to quantify?
  - Number of occurance in the sequence

## Apr 26

- 3mers as features:
  - Minus and plus not separatable


- define: main strand by bulk count
  - number of different pos ~ main strand
    - is this true on single-cell level? (definition of burstness)

  - Sequence analysis on main strand
    - if main plus, main minus, minor plus, minor minus
    - should give different enhancer the same weight but within each enhancer, the most expressed the largest weight (TODO)
    - edit extract_CTSS_seqs function

    - if we filter out TPM<1 transcripts
      - we get patterns
      - we are saying that: highly expressed enhancers have such pattern:
        - Major minus, minor plus and major plus are GC rich
        - Minor minus is near poly T and TATA box

- define burstness:
  - negatively related with number of diff pos
  - positively related with TPM
  - can be calculated at both single-cell and bulk-cell level
    - are they correlated?


## Apr 27
- Making everything in a table
  - enhancer: expression, main strand direction, burstness...


- weight each enhancer similarly
  - do not see any sequencial differences between major and minor strand transcripts

- Main strand = plus 
  - count plus big
  - plus transcripts are located at plus boundary


## Apr 28
- CpG ratio higher in the direction of main strand


## Apr 29

- Burstness:
  - Fraction of tag position re-occurance
    - issue of mapping
  - Mean - Median 


- 50% of chromosome are covered by gene (intron)
  - are enhancers located at intron?

- Burstness, is it defined before?
  - burst size
    - Number of RNA molecule
  - burst frequency (h^-1)
    - how to infer in our data?
    - support?

- Is burstness associate with width?

- Are bursting enhancers more validated?

- Is bursting adaptive?

- Count number of bursting event in cells / enhancers
  - ~ frequency 

- Are there any cell full of burst events?

- burst on enhancer associated with burst on other? co-burst?

- examine sc direction preference: are bidirectional not possible in single cell
  - for each enhancer with at least 1 cell > 2 pos:(how many)
    - are they same direction or different direction?

## Apr 30
Pooled single cell signals (>= 2 pos) are also very unidirectional bursting (plus minus gini large)
  - should not extend this because bursting uses a non-stocastic model

Cell - specific signals are more unidirectional than pooled signals.
- How to test this?
- How to visualize this?


## May 1

Look at enhancers with low gini pooled but high gini single-cell
- chr22:30606240-30606705: high GC, intron, H3K27Ac signal, HACER validated
- chr17:7745936-7746377: high GC, intron, H3K27Ac signal, HACER validated
- chr17:8113067-8113546: high GC, intron, H3K27Ac signal, HACER validated
Look
- chr19:58912475-58912923: CpG island, 


Look at enhancers with largest bursting size
- bursting size may related to histone modification


Enhancers prefer one direction on single cell level
  - Not either plus or minus
    - but plus 


- Enhancers with a small Pooled Gini / Mean gini: 
  - bidirectional on bulk and unidirectional on single-cell
    - looks like they are conserved on sequence



TODO:
  what are the relationship between 
  - burst size
  - burst direction gini
  - burst freq
  - v.s.
  - GC content
  - intron / integenic (txType)
  - context
  - repeats
  - conservation

what is a bursting enhancer?
- enhancer that is 
  - active
  - near Pol II clustering

Why is it half open?
- it's half active? determined by sequence? no sequence feature
- determined by histone? then we should see transcriptional context
- TF occupy space

what is the use of half open?


## May 2

Write function to test GC difference between enhancers


- bursting events: 1, 1~10, 10~100, 100+
  - not any difference between -1000~+1000 logo of top 50 and bottom 50
  - top 10% mean burst size biased to plus, bottom 10% burst size biased to minus
    - why? (TODO)

- no gc% difference between groups

- no gc% difference 


repeats  enhance

## May 4
- Everything is bursting
  - burst size not related to burst frequency
  - burst freq ~ bulk express (TODO)
  - burst freq ~ GC 
    - filter(total_counts_max>100) more significant
  - burst freq ~ sequence (TODO)
  - burst freq ~ txType
  

- In developmental gene, if you really want to shut down one gene, you will have essentially no genes around it

## May 5
count word freq in -5 ~ +1 for different burst freq level

## May 6
when one-way bursting is happening, is the other way closed?

Is enhancer isolated in a sparse region? 

## May 7
Some enhancers are close to each other, what are they? 
  - Are they also simultaneously active on sc level?

## May 10
- correct the intron intergenic labelling
- are npos ~ ncell, burst size?
- is ncell~ mean burstsize?
- what are special with the bidirectional bursting enhancers?
  - do they have bidirectional context?
  - do they have more GC rich regions?

bidirectional enhancers: are they enhancers with GGGG and CCCC on two side?
  - write function to extract all pos with GGGG and CCCC.

## May 12 

- calculate enrichment score for motifs
  - enrich in the center of enhancer when averaging
  - enrich like a down-up-down-up-down pattern in each single enhancer

hmm model

- mistake: 
  - CTSSs_bulk_t0 had different # of ctss than others

- lose the cutoff of enhancer to include more BC clusters
  - what is bulk support?
    - are highly supported enhancers:
      - more single cells
      - highly expressed?
      - more fantom validated?

- one way to find divergent motifs: test for bimodal distribution
  - should be one modal in large scale (e.g. -1000 ~ +1000: enriched in enhancers)
  - should be bi-modal in small scale (e.g. -400 ~ +400)

## May 15
- Use MEME-suite to discover motifs associate with
  - burst freq
  - burst direction
  - ..

- MA plot
  - M: log(exp1 / exp2)
  - A: 0.5* log(exp1 + exp2)


## May 16
Do frequently bursting enhancers have more "certain motif" throughout the whole region?
- no 
- why? because the meme-chip only uses trimmed centered 100 bp of the whole sequence set
- and it only samples 600 sequences into the DREME, which means we should enter small datasets

Do they have more certain moitf on sc TSSs?
- yes
- are they expected?
  - for some reason they are more bursting, do they naturally will have more specific motifs?
  - it really depends on how TSSs pos/ TSSs size resemble each other across cells
    - if they do not resemble, we are actually choosing enhancers with more size

Do they have more certain moitf on bulk TSSs?


## May 17

- We should use TSSs centered sequences (100bp length) to discover motifs
  - or else we change the -ccut flag to 0
  - to what ever length we want
  
- The size of input fasta should < 600 sequences
  - -meme-searchsize 0
  - or else we change the -searchsize flat to 0
  - how to do it?
  - in extract

- Use high order marcov model 
  - -order 5
  - change markov_order to 5 to normalize specific pattern in non-coding region
    - Using a higher order background model can often result in more sensitive detection of motifs. This is because the background model more accurately models non-motif sequence, allowing MEME to discriminate against it and find the true motifs.
  - does not make any difference

- Possible (fast and accurate) solution:
  - extract motifs from TSSs and then scan enhancer whole region for motifs
  - use the scan result to say sth about the enhancer features

- should centered at -20 bp
  - extract -70 ~ 30 bp 
  - center cancel -centrimo-local


## May 18

- reverse trancriptase often adds extra G in the 5' end
  - quality control

- to prove direction preference is sequence defined:
  - turn every sequence to their preference direction
  - make groups: 1. half of the sequence on preference; 2. half of the sequence off preference
  - use meme-chip to find motif related

- positive control of meme-chip
  - low gc content enhancer and high gc content enhancers

- to prove enhancers do not have direction in cell: 
  - compare prefer-plus enhancers and prefer-minus enhancers


- Preference related motif
  - TGANTCA (t24)
  - G rich motif (t0)


- to make more sense:
  - classify enhancers to 3 instead of 2 classes: minus, bidirectional, plus
  - when extract CTSSs from enhancer: sample once instead of 3 times
  - when extract CTSSs from enhancer: sample from bulk data
  - compare bidirectional enhancers to unidirectional ones (combine plus and minus TSSs)
    - bidirectional ones should have more flat GC bell shape (less enrichment of GC in the center)


## May 20

- use heatmap to visualize distribution of motif
- also visualize the repeat mask

## May 21

further explore the repeat masker:
- does it work on 100bp long sequences?
- or I have to compare and extract from masked longer regions

My ideal report
- Talk about bidirectional enhancers
- Single-cell as a snap shot
  - How?
    - burst
    - unidirection
  - sequence base
    - treatment specific TF motifs
      - ~ burst freq
    - treatment intact Pol ii motifs
  - snap shot: transcriptional context
    - theoritical backup(?)
    - Cloudy with a Chance of Insights: Context Dependent Gene Regulation and Implications for Evolutionary Studies
    - The context of gene expression regulation


## May 22
-  draw Figure 1b of 1. Andersson, R. et al. An atlas of active enhancers across human cell types and tissues. Nature 507, 455–461 (2014).

- draw figure2d of 1. Andersson, R. et al. An atlas of active enhancers across human cell types and tissues. Nature 507, 455–461 (2014).
  - Are single cell more likely to be enhancers?
  - Are they more CGI located?
  - Are they more repeated?

- draw cumulative expr context like figure 2b
  - increase window size but fix one side

- plot the actual TSSs and 5'Splice motif, TATA box etc.
  - expect TATA box on plus and minus direction (line plot)
  - transcription terminate site

- see top Top enriched/depleted motifs

- why eRNA negatively correlated with context?
  - competition of Pol ii ?
  - rapid degradation

- Where does TSSs bind?
  - a better background model
  - order 5 marcov model
  - compare TSS -100 +10 with randomly selected 100 length regions on enhancers

  - a better positive control
    - strong promoter TSSs
    - 1000-1000 around strong promoters

  - what to search
    - JASPAR pol II
    - https://pdfs.semanticscholar.org/22f3/07237483498a07240d55d5257b2d67ab3822.pdf

- Are ubiquitous enhancers 
  - having more sc expr?
  - burst size?
  - burst freq?


- bulk and single cell correlation 
  - bulk: positive correlation: more chances open together
  - single cell: negative correlation: competition
  - CTCF positions https://rdrr.io/bioc/sevenC/man/motif.hg19.CTCF.html

- cell develoment association with the enhancers?
  - This dataset was intended to study ...,  so we identify over-represented FOXA motif
    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7025859/
    https://www.nature.com/articles/srep34962.pdf?origin=ppub
    
      master transcriptional factors
      pioneer transcriptional factors
        does it infer order of trancription? Mediator and
    Overexpression of FOXA protein
      https://www.sciencedirect.com/science/article/pii/S037811191831165X?casa_token=ip8HkJXzpnYAAAAA:kZ_CSTaK-35XSZmjzuX4oT1hrlRmLZjaSyIpB3ivSuMDhwW1PYxhy9KA32UMpltFUftABVRsQQ

    Enhancer reprogramming in tumor progression: a new route towards cancer cell plasticity
Cohesin
          https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3653129/
  - but we are more curious about ...
  - which is independent of ...

- how to predict enhancer - TSSs interaction?
  - No. We build our assumption on 40% enhancers regulate nearest TSSs
  - or CTCF positions

- Antisense transcription associate with genome structure
  - https://www.tandfonline.com/doi/full/10.1080/19490992.2015.1130779

- Are the direction preference associate with direction of nearest TSSs?
  - how to explain this?
    - gradient of histone?

- Are bidirectional enhancers more open?
  - e.g. bidirectional TSSs neighbors
    - more near?
    - more abundance?

- How to define neighbours?
  - TAD https://www.nature.com/articles/s41576-018-0060-8?WT.feed_name=subjects_chromosomes
  - hg 19 data http://promoter.bx.psu.edu/hi-c/publications.html


## May 23 

- Find retrotransposons repeat surrounding enhancers
  - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4914423/
  - HOW RETROTRANSPOSONS SHAPE GENOME REGULATION

  - to do next:
    - repeat category enrichment 
    - TSS relative location with repeat    not overlapping 
    - repeat direction and TSS direction correlation   cannot intepret
    - repeat classification   mostly SINE
    - distribution of repeat class ( fill in point from start to end) 
    - Pol ii related repeat class
    - intron intergenic


## May 24
- Build a linear model:
  - How is a burst event determined?
    - event direction ~ context direction + bulk preference + sequence reason...
    - how does context direction determined?
      - cell state?
      - sequence?


- Is eRNA expression minor compared with nearby promoter?

- Is eRNA pos wide spread compared with nearby promoter?

- how are those widespread differ from those narrow spread?

## May 28

- Use the count=1 TSSs as predictor of open chromatin
- - the figure looks good (May27_context_cumulative_TPM1_t24)
  - what is the backup support?

- Context direction:
  - bidirect enhancers v.s. unidirect enhancers
  - Compare with promoters

- Active enhancer v.s. non-active enhancer
  - Good results!
    - May28_context_cumulative_countpos_t24.png
    - May28_context_cumulative_TPM_t24.png

- sequence analysis of count 1 TSSs:
  - are they different from others


- PCA to discribe context

- high repeatability of promoter context ( and enhancer context)
  - how to quantify the repeatability


## May 30

lm(temp %>%
    lm(bias_d~direction+cid+eid,data=.))


## May 31

can context predict active enhancer?

very hard 


## Jun 1

Start to write my report:

1. The first article finding eRNA expression, suggest Enhancer DNA sequences are replete with clusters of binding sites for transcription factors
 Orom, U. A. & Shiekhattar, R. Long noncoding RNAs usher in a new era in the biology of enhancers. Cell 154, 1190 (2013).

2. What is the function of eRNA expression. Why is it wide spread?

3. We use the alignment of the NC paper, but find enhancers using cage fightr


4. 30% human 5′-capped transcripts use repetitive sequence- associated TSSs (27).
   1. Most bidirectional ... use repetitive sequence associated TSSs.
   2. compare with promoter

5. enhancers defined by open chromatin, H3K4me1 ...
   1. now we define by CAGE
   2. For example, one study found that eRNA served as a better marker of active enhancers than H3K27ac (Tyssowski et al., 2018). 
   3. eRNA transcription per se is able to identify in vivo enhancers that are “caught in the act” of regulating mRNA/gene expression and does so independently of the local chromatin state or transcription factor binding
  
6. bidirectionality is a marker as a degree of enhancer or promoter 
   https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5828394/pdf/42.pdf


7. recent studies indicate that key lineage factors FOXA1 in hepatoma cells and GATA1 in erythroid cells remain associated with select sites on mitotic chromosomes, some of which have enhancer markings (Caravaca et al., 2013; Kadauke et al., 2012). Moreover,

8. For enhancer looping as key, and likely required, for transcription activation.
Enhancer function: Mechanistic and genome-wide insights come together
new enhancer loops pre-cedes transcription activation in the b-globin locus

Negative correlation eRNA and mRNA suggest --> 
eRNA expression happens pre-cedes transcription activation

8. TAD
   1. Megabase-sized TADs (0.5–3 Mb)


9.  
Extract 10000bp near enhancers Find the boundary of repeat elements

10. Enhancer activity ~ DNA methelation
  A consistent theme appears to be loss of DNA methylation and transition to H3K4me1 and either H3K27me3 or H3K27ac to poise or activate enhancers, respectively (Gifford et al., 2013; Sheaffer et al., 2014).

11. Transcription of an eRNA upstream of the Arc gene requires the intact promoter and gene hinting that eRNA transcription and enhancer-promoter communication are linked (Kim et al., 2010).

12. eRNA stablized the promoter-enhancer contacts, and Mediator interactions


stabilization of loops between enhancers and target genes.

13. Take home message:
    1.  mix different proportion of something makes artifacts

14. Bidirectional transcription at enhancer sites generates comparatively shorter (0.5-2kb) and non-polyadenylated eRNAs. (https://en.wikipedia.org/wiki/Enhancer_RNA)

using cage, we can capture both eRNAs

"Diversity and Emerging Roles of Enhancer RNA in Regulation of Gene Expression and Cell Fate Preston"
15. 
traditionally function independently of their position and orientation with
respect to their target gene


16. actively transcribed (De Santa et al., 2010; Kim et al., 2010). 

17. the role 
    1.  byproduct
    2.  

18. it remains unclear whether eRNA which functions are generalizable.

19. nd hepatocytes among those with the highest abundance of cell-specific enhancers, a higher ratio of enhancers/gene, and high enhancer transcription (Andersson et al., 2014a). Indeed,

why we select lung cancer cell line

cell type specific ---> development 

20. occurs very early in the gene
transcription process, preceding mRNA expression from adjacent
cis loci

21. excellent marker for segregating active versus quiescent enhancers,
    1.  eRNAs were prominent at looped enhancers, but whether eRNA played a role in that looping process was less clear (Wang et al., 2011; Sanyal et al., 2012).

22. The details and mechanisms of enhancer transcription,
initiation, and elongation are similar to those found at promoter
sites, and the similarities between eRNA, lncRNA, PROMPT, and
mRNA transcription have been thoroughly reviewed

23. PAS which are bi-directionally found immediately downstream of enhancer TSSs

24. This biological coupling of eRNA initiation and termination suggests careful regulation of eRNA and implicates eRNA molecules as more than simple accidental transcriptional noise.
    1.  our inherit of shared enhancers also suggest not simple accidental noise

25. most eRNAs are bidirectionally transcribed, unspliced, and non-polyadenylated (Andersson et al., 2014a).
    1.  produce eRNA bidirectionally but with the predominant production coming from either the negative or positive strand
    1.  we study this kind of eRNA
    2.  but we discover that they are also uni 
            exclusively unidirectionally transcribed from either
            strand (Kouno et al., 2019).
        
    3.  We further discover that burst freq do not burst size 
    4.  the predominant production can be an artifact

26. The nuclear RNA exosome is known to degrade the 3? end of eRNAs (Andersson et al., 2014a,b; Pefanis et al., 2015; Imamura et al., 2018). Because the distance between TSS and PAS is inversely correlated with RNA exosome sensitivity
    1.  due to fast degrade
        1.  Our high abundance eRNA
            1.  is either special on 3' end sequence polyadenylated
                1.  not possible to test in CAGE data
            2.  or temporary high burst
            3.  integrate as a super enhancer
                1.  how to test?


27. some promoters may act as enhancers depending on the context (Mikhaylichenko et al., 2018; Rennie et al., 2018).
    1.  do not have TATA box


28. A central issue regarding our current understanding of eRNA is whether the transcription of eRNA or the eRNA itself is primarily responsible for observed fuctnions.

29. On a genomic scale, eRNA transcription and induction of mRNA transcription at neighboring genes are correlated
    1.  but on single cell scale, it is not correlated

30. which eRNAs are functional, and how eRNA function is linked to structure and localization.
    1.  use a regression model 
        1. TODO  bias_d = active * eid 

31. specificity of eRNA-promoter interactions

32. eRNA-protein interactions may play roles in protein recruitment, altering protein interactions, and providing scaffolding.

33. “transcription factor trapping” wherein nascent RNA from enhancers and promoters increase the affinity of otherwise weak DNA-TF interactions creating a kinetic sink that would “trap” escaped transcription factors (Sigova et al., 2015)

34. seldomly highly expressed
    1.  may be very functional important

35. knockdown of eRNA has been shown to decrease the accessibility of their respective enhancer
regions (as measured by DNase-seq or ATAC-seq),

36. what model does it support?
    1.  maybe enhancer - promoter interaction
        1. TODO  how to test?  direction ~ active * bias


37. This leads to eviction of NELF from chromatin and RNA elongation (Schaukowitch et al., 2014; Shii et al., 2017).
    1.  Integrator has been shown to interact with NELF and DSIF pausing machinery at the proximal promoter and to impact transcription initiation and RNAP II pausing (Gardini et al., 2014; Skaar et al., 2015).
    2.  Our data measures burst frequency on cell frequency
    3.  Integrator is uniquely positioned to promote productive eRNA interactions with transcriptional machinery (Lai et al., 2015).


38. draw enhancer - cell heatmap

39. What enables enhancers to transcribe bidirectionally?
    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4874897/
    - TATA?
    - CCAAT?
    - CpG island?
    -  


40. Cryptic transcripts?
  - TODO To remove these transcripts, we need to subset CTSS object

41. transcriptional burst review
    https://www.annualreviews.org/doi/pdf/10.1146/annurev-biochem-011520-105250

42. chromatin accessibility review
    https://www.nature.com/articles/s41576-018-0089-8.pdf

43. single cell ATAC
    https://www.nature.com/articles/s41467-018-07771-0.pdf
    https://www.nature.com/articles/nature14590.pdf

44. condensation droplet burst phase separation
      https://www.biorxiv.org/content/10.1101/737387v2.full

45. determinant of burst
    https://pubs.rsc.org/en/content/articlelanding/2017/mb/c7mb00154a/unauth#!divAbstract

46. Transcriptional Burst Initiation and Polymerase
Pause Release Are Key Control Points
of Transcriptional Regulation
    https://www.cell.com/molecular-cell/pdf/S1097-2765(18)30946-8.pdf

47. genomic encoding of burst
    https://www.nature.com/articles/s41586-018-0836-1.pdf

TODO longer range on the two sides of enhancers


1.  The fact that promoters can be unidirectional in some tissue types and bidirectional in others [22] supports the notion that core promoter sequences allow for bidirectional transcription but that this capacity is then regulated by secondary mechanisms such as cell-type specific TFs that promote either one or both transcripts in the pair in response to the different needs of the cell


TODO redraw the GC trend plot by direction class
analyze direction related motifs  

42. In the first, antisense transcription is a by-product of forward transcription and an NDR, meaning that the cell may need to utilize mechanisms to prevent the accumulation of detrimental antisense transcripts. In the second, there is a function for antisense transcripts with a number of possible effects including the regulation of sense transcripts. 

43. over 50% enhancer promoter interactions happen within 100kb, 1E6 is sufficient
    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4252644/


44.  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4553066/pdf/nihms709552.pdf
  At present, the possible
localization and function of paused Pol II at active enhancers is not clear.

A direct role of eRNA in the release of paused Pol II has been recently proposed [109]. 


## Jun 2

- Check if hg38 or hg19
  - hg19


- Write about GGGG?
  - maybe pause related
    - because
  - can hardly be artifact, cause the control  
  - what is that then?
    - strand invasion events
      - maybe biased to low expression TSSs like enhancers
      - https://www.bioconductor.org/packages/release/bioc/manuals/CAGEr/man/CAGEr.pdf
      - findstrandinvaders
        - IF artifact is not removed from the bed12 file,
          - many low support enhancers can have mapping errors,
          - maybe they are not capped.
          - we define our enhancers using bulk enhancers

- advantage of single cell cage: use single cell control 
- advantage of CAGE: at single molecule resolution

Experimental design
- multiplexing strategy based on colored stains


bulk enhancer  --> conserved
  - single cell signal: not frequent
  - how to interpret?
    - hetero...


how to remove pca replicates
  - so they are different molecules
  - collapses the pairs that have exactly the same alignment coordinates

how to regress out the batch effect
- hard to imagine the batch effect cause directional change
- but anyway, we can use the color batch


"We also suggest for the first time that while eRNAs appear lowly expressed in bulk data, they can be expressed at similar levels to gene promoters within single-cells, although they are expressed in a more restricted subset of cells—i.e., displaying transcriptional bursting"



Why unidirectional?
  "steric occupancy of RNA PolII on each strand at any given time,"
  "proximal to chromatin loop structures such as CTCF"


More support should come from multiple TSSs.


Introduction:

1. eRNA perfect indicator for activation of enhancer
   1. now we propose it can be potential indicator for sc topology

2. there is a historical debate: bidirectional or unidirectional (not for enhancer)
3. more evidence suggest directional 
4. but it is dangerous to make conclusion: at given time point, is enhancer able to transcribe both sides?
5. Why do direction matter?
   1. suggest CTCF
   2. suggest occupancy
  
6. what do eRNA mean? 


Results:
1. identification of bidirectional enhancers
2. unidirectional eRNAs
   1. how many enhancers are either minus
3. sc preference ~ bulk preference
4. sequence analysis


Discussion:

1. What is interacting with this enhancers?

2. How does enhancer function?
   1. A crazy idea: enhancer opens downstream chromatin
   2. 

3. How does eRNA function?

4. Why eRNA transcribed before mRNA ?

4. eRNA not only perfect indicator of enhancer activity
   1. the direction related with linear epigenome. histone distribution on chromatin  (turnover)
      1. TODO https://www.nature.com/articles/s41576-018-0089-8.pdf
   2. 3D topology? maybe, but we
   3. CTCF distribution
   4. indicator of chromosome-accessibility in real time
      1. co-analysis
      2. https://www.nature.com/articles/nature14590.pdf
   5. will be further useful in chromatin state description

5. Can we do sc ATAF and sc CAGE seq together?

6. Pol ii recruitment only occurs during transcriptional burst
   1. https://www.cell.com/molecular-cell/pdf/S1097-2765(18)30946-8.pdf
      1. Recently, some have suggested that bursting could be related to condensation or phase separation of transcriptional proteins, and others may wish to investigate this possibility (Hnisz et al., 2017; Cho et al., 2018). 
      2. CITE THIS https://www.biorxiv.org/content/10.1101/737387v2.full
   2. crazy: hence eRNA can recruit Polii ? 


7. TODO How is bidirectional enhancers special in sequence feature?

8. Transcription is regulated at several levels, including enhancer–
promoter interactions11, the formation of the transcription pre-initiation
complex, recruitment of RNA polymerase (Pol) II26, initiation of Pol II27
and control of Pol II elongation
https://www.nature.com/articles/s41586-018-0836-1.pdf

9. Another crazy model:
   1.  phage separation ~ enhancer direction
       1.  enhancer transport Pol ii to promoter
       2.  supported by increased frequency


10. eRNA- mRNA co-expression is infrequent
    1.  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5389544/
        1.   only a fraction of the cells respond to the stimulus at any given time. 
             1.   (i) that co-expressing alleles are rare, and (ii) that co-expression infrequently occurs in a closed enhancer–promoter configuration. One interpretation of this observation is that enhancer–promoter interactions are dynamic
             2.    feed-forward loop, whereby basal eRNA transcription facilitates the recruitment of TFs and co-regulators, which then further remodel chromatin and increase the frequency of eRNA and mRNA transcription initiation
             3.    their transcription is mostly repressed once looping interactions are established. Elucidating the dynamics of enhancer–promoter interactions in the presence or absence of eRNAs will require live cell analyses at high spatial and temporal resolution.
             4.    figure 5c
             5.    *** Our data supports this observation very much
                1.   and against a model that eRNA is necessary for maintain promoter-enhancer interaction during mRNA transcription
        2.   eRNAs are thought to facilitate the transition of paused RNAPII into productive elongation by mediating the release of the negative elongation factor NELF (24).
        3.   The low transcription frequency of eRNAs observed at the FOXC1 and P2RY2 enhancers, however, does not favor such a model
        4.    eRNAs in facilitating or stabilizing enhancer-promoter loops through interactions with cohesin and the mediator complex


11. previous studies (smFISH) focus on limited pair of enhancer-promoters
    1.  CAGE view the genome
    2.  We suggest, on the genome level, eRNA functions in a way to open the downstream chromatin
        1.  if downstream is open first, then we will see more TSSs
        2.  if downstream is open then, we predict more 



12. enhancers control of transcriptional bursting
    1.   https://www.cell.com/cell/pdfExtended/S0092-8674(16)30573-6
    2.   enhancers regulate the frequency of transcriptional bursts.
    3.   1 enhancer sychronize 2 promoters

13. r measured eRNA half-lives (7.5 min) (65, 69).
    1.  https://www.annualreviews.org/doi/pdf/10.1146/annurev-biochem-011520-105250
    2.  most eRNAs are unstable, HENCE NOT LIKELY to act in trans

14. eRNA transcription is correlated with developmental enhancer activity
    1.  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5828394/pdf/42.pdf
    2.  When transcription is initiating at an enhancer, the enhancer is very likely to be active; however, the converse does not always hold true
    3. Transcription cannot be detected at all active enhancers, suggesting that at least for this subset, eRNA is not mechanistically required for their activity
    4. Andersson et al. (2014a) found that 20%–33% of nontranscribed regulatory regions can activate transcription in an in vitro enhancer assay, suggesting that they may also function as enhancers without eRNA transcription.
    5. This indicates that promoter activity has intrinsic directionality and suggests the presence of directional sequence motifs within enhancer elements.
       1. crazy: this may related to bidirectional distribution of TF motifs


15.  the role of transcription factors in initiating accessibility remodelling
     1.   https://www.nature.com/articles/s41576-018-0089-8.pdf
     2.   n contrast to closed chromatin, permissive chromatin is sufficiently dynamic for transcription factors to initiate
sequence-specific accessibility remodelling and establish an open chromatin conformation

16. pioneer factor: opens chromatin
    1.  https://www.nature.com/articles/nbt.2798.pdf
    2.  fos jun as a pioneer
        1.  ELK1: ETS pioneer,Nondirectional,
        2.  PF0009.1 TGAYRTCA , pioneer
        3.  sp1 pioneer   ---> TODO unidirectional  Sp/KLF!


17. GC motif provides a good location to pioneer motifs


## Jun 4
- plot repeat elements around promoter

- regression out txType factors
  - how to do this

- Plot shows that upstream same direction 
  - active ~ bias towards center + eid 
    - not significant

- statistics shows preference to downstream
  - we have not found preference to 



- todo Introduce CAGE


## Jun 5

How explain cannot find sequence related with burst size?
  - maybe the statistical power of the burst size is not enough
  - maybe cell is highly dynamic

It is unidirectional, but not so much

Check directional bias of motif
  - TATA
  - PyPu
  - INR

Availability or interactivity (TAD)

## Jun 7

- thoughts during writing method
  - TPM: "However, such simple approaches assume that there are no systematic variations between samples (which are not controlled by the experimenter) that may cause the absolute tag-counts to vary across experiments."
  - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2728533/
  - TODO check if there is systematic variations

## Jun 9 

- TODO
  - provide time-specific figures in supplementary data, provide pooled figures in report
    - and mention that time specific data is like the same

  - mention that in whole time figure duplicated enhancers are collapsed

  - Our novalty: the first single cell CAGE data. We are able to analyze enhancer context.

## Jun 10

- TODO
  - write about improving validation rate by time specific  

  - bidirectional or unidirectional? 
    - bidirectional on biological sequence
    - but unidirectional on physics (discrete on single cell level)

  - drawback:
    - A549 is a cancer cell line with massive chromosome rearrangement and replication.
    - it will not affect our results
    - but it may cause artifact of "bidirectional" expression 

  - solution:
    - use other cell lines
    - especially exosome knock out cell lines


## Jun 11
Try to do:
- centrimo on burst freq 
  - how to extract TSS
    - 1. CTSS(instead of enhancer) of different burst freq
    - 2. extract from bulk (since bulk ~ sc)

- or abandon the burst freq.
  - focus on directional biased motif
  - do not mention GGGG

- or abandon motif enrichment
  - draw CG% plot on different burst size (!)


- Bulk two side:
  - without INR, without TATA

- do burst freq ~ bidirection?
  - if more frequent, can it becomes bidirectional?


- consider the average half life of eRNA
  - If observe bidirectional: is it a remain