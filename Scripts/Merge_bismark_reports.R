



Sequence pairs analysed in total:
Number of paired-end alignments with a unique best hit:
Mapping efficiency:
Sequence pairs with no alignments under any condition:
Sequence pairs did not map uniquely:
Sequence pairs which were discarded because genomic sequence could not be extracted:

Number of sequence pairs with unique best (first) alignment came from the bowtie output:
CT/GA/CT:       1803630 ((converted) top strand)
GA/CT/CT:       0       (complementary to (converted) top strand)
GA/CT/GA:       0       (complementary to (converted) bottom strand)
CT/GA/GA:       1804101 ((converted) bottom strand)

Number of alignments to (merely theoretical) complementary strands being rejected in total:

Total number of C's analysed:
Total methylated C's in CpG context:
Total methylated C's in CHG context:
Total methylated C's in CHH context:
Total methylated C's in Unknown context:


Total unmethylated C's in CpG context:
Total unmethylated C's in CHG context:
Total unmethylated C's in CHH context:
Total unmethylated C's in Unknown context:


C methylated in CpG context:    77.6%
C methylated in CHG context:    1.2%
C methylated in CHH context:    1.4%
C methylated in unknown context (CN or CHN):    20.3%



total.cpg.content = "Total methylated C's in CpG context:"
total.chg.content ="Total methylated C's in CHG context:"
total.chh.content ="Total methylated C's in CHH context:"
total.cuknown.content ="Total methylated C's in Unknown context:"


(6501750/(6501750+1872658) ) *100
$percent_meCpG = sprintf("%.1f",100*$counting{total_meCpG_count}/($counting{total_meCpG_count}+$counting{total_unmethylated_CpG_count}));



