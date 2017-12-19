fafile="Ecoli_ATCC8739.fa"
refGene="ecoli_refGene.txt"
outputfile="retrieveseq.txt"
def retrieveseq(fafile,refGene,outputfile):
    f=open(refGene,"r")
    g=open(fafile,"r")
    h=open(outputfile,"w+")
    allseq=g.readlines()[1]     
    for eachline in f:
        lst=eachline.strip().split("\t")
        strand=lst[3]
        txstart=int(lst[4])
        txend=int(lst[5])
        if strand=="+":
            seq=allseq[txstart-1:txend]
        else:
            seq=revcom(allseq[txstart-1:txend])
        pro=transform(seq)
        h.write("%d\t%d\t%s\t%s\t%s\n" % (txstart,txend,strand,seq,pro))
def revcom(seq):
    revdic={"A":"T","G":"C","T":"A","C":"G"}
    i=0
    tmpseq=""
    while i< len(seq):
        tmpseq=tmpseq+revdic[seq[i]]
        i+=1
    return tmpseq[::-1]
  
def transform(seq):
    codedic={"TTT":"F", "TTC":"F", "TCT":"S", "TCC":"S", "TAT":"Y", "TAC":"Y", "TGT":"C", "TGC":"C", "TTA":"L", "TCA":"S", "TAA":"*", "TGA":"*", "TTG":"L", "TCG":"S", "TAG":"*", "TGG":"W", "CTT":"L", "CTC":"L", "CCT":"P", "CCC":"P", "CAT":"H", "CAC":"H", "CGT":"R", "CGC":"R", "CTA":"L", "CTG":"L", "CCA":"P", "CCG":"P", "CAA":"Q", "CAG":"Q", "CGA":"R", "CGG":"R", "ATT":"I", "ATC":"I", "ACT":"T", "ACC":"T", "AAT":"N", "AAC":"N", "AGT":"S", "AGC":"S", "ATA":"I", "ACA":"T", "AAA":"K", "AGA":"R", "ATG":"M", "ACG":"T", "AAG":"K", "AGG":"R", "GTT":"V", "GTC":"V", "GCT":"A", "GCC":"A", "GAT":"D", "GAC":"D", "GGT":"G", "GGC":"G", "GTA":"V", "GTG":"V", "GCA":"A", "GCG":"A", "GAA":"E", "GAG":"E", "GGA":"G", "GGG":"G"}
    pro="M"
    if len(seq)%3!=0:
        print "seq long is not 3 times"
    else:
        n=len(seq)/3
        i=1
        while i<n:
            codestart=i*3
            code=seq[codestart:codestart+3]
            pro=pro+codedic[code]
            i+=1
    return pro
if __name__ == "__main__":
    retrieveseq(fafile,refGene,outputfile)
