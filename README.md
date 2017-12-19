# EcoliCallMutationType
This is a module call Ecoli mutation types, such as intergenetic,synonymousSNV,nonsynonymousSNV,and so on,to run it ,we need fasta file ,snpfile, and refGene file

* For snpfile, with header,eachline has 4 tab-delimited column:$Chrom,$Postion,$Ref,$Alt,ALT should be added indel imfomation such as "+ATC"or "-GAGC,which means insert or del any bases.For refGene file, without header,each line has 16 tab-delimited columns: $bin, $name, $chr, $dbstrand, $txstart, $txend, $cdsstart, $cdsend, $exoncount, $exonstart, $exonend, $id, $name2, $cdsstartstat, $cdsendstat, $exonframes. The only real important thing is $chr (chromosome), $dbstrand (strand of the transcript in reference genome), $txstart, $txend, because always Ecoli trscript reion is same as cds region and exon region, and the rest columns can be filled anything
* optional annfile: the annotationfile with header, eachline has eight tab-delimited; coloumn:$species,$start,$end,$strand,$type,$locus,$productId,$description. For example, Ecoli,38,1451,+,protein,Ecolc_0001,ACA75689.1, chromosomal replication initiator protein DnaA. Only important thing is $start,$end,$type,$locus, the type(including tRNA or rRNA) and locus both are selectively.
* default only consider mutation type without indel, because indel is not accuracy with sequencing error according to our research.

```
Usage: judgemuttype.py [options] arg
Options:
  -h, --help            show this help message and exit
  -a REFFA, --reffa=REFFA
                        read data from REFFA
  -f SNPFILE, --snpfile=SNPFILE
                        read data from SNPFILE
  -r REFGENE, --refGene=REFGENE
                        read data from REFGENE
  -o OUTPUTFILE, --outputfile=OUTPUTFILE
                        output data into OUTPUTFILE
  -n ANNOTATION, --annotation=ANNOTATION
                        annotation file including locus or RNA/tRNA(optional)
  -l, --addlocus        choose to add locus from annotationfile you
                        have(optional)
  -R, --RNAtype         choose to add RNA and tRNA type from annotationfile
                        you have(optional)
  -d, --description     choose to add description(optional)
  -i, --ignoreindel     selectively choose to ignore indel(optional,default is
                        True)
  -q, --quiet           no standard output
```

### Ex:
```./judgemuttype.py -a Ecoli_ATCC8739.fa -f test.snp.txt -r ecoli_refGene.txt -o MutType.txt -n Ecoli.anno -l```

### output
```
Chrom   Position        Ref     Alt     Muttype Locus   Strand  Ingenepos       Altpropos       Refpro  Altpro
Ecoli   39      G       A       intergenetic    .       .       .       .       .       .
Ecoli   100     C       T       intergenetic    EcolC_0001      .       .       .       .       .
Ecoli   1072    G       C       intergenetic    EcolC_0001      .       .       .       .       .
Ecoli   2891    G       A       nonsynonymousSNV        EcolC_0003      +       336/1074        111/358 M       I
```
