#!/usr/bin/env python

#This is a module call Ecoli mutation types, such as intergenetic,synonymousSNV,nonsynonymousSNV,amino acid changes, and so on,to run it ,we need fasta file ,snpfile, and refGene file;

#For snpfile, with header,eachline has 4 tab-delimited column:$Chrom,$Postion,$Ref,$Alt,ALT should be added indel imfomation such as "+ATC"or "-GAGC,which means insert or del any bases.For refGene file, without header,each line has 16 tab-delimited columns: $bin, $name, $chr, $dbstrand, $txstart, $txend, $cdsstart, $cdsend, $exoncount, $exonstart, $exonend, $id, $name2, $cdsstartstat, $cdsendstat, $exonframes. The only real important thing is $chr (chromosome), $dbstrand (strand of the transcript in reference genome), $txstart, $txend, because always Ecoli trscript reion is same as cds region and exon region, and the rest columns can be filled anything

#optional annfile: the annotationfile with header, eachline has eight tab-delimited; coloumn:$species,$start,$end,$strand,$type,$locus,$productId,$description. For example, Ecoli,38,1451,+,protein,Ecolc_0001,ACA75689.1, chromosomal replication initiator protein DnaA. Only important thing is $start,$end,$type,$locus, the type(including tRNA or rRNA) and locus both are selectively.

#default only consider mutation type without indel, because indel is not accuracy with sequencing error according to our research.

from optparse import OptionParser
import sys
import re
import os
from retrieveseq import retrieveseq,transform,revcom

usage = "usage: %prog [options] arg"
parser = OptionParser(usage)
parser.add_option("-a", "--reffa", dest="reffa",
                      help="read data from REFFA")
parser.add_option("-f", "--snpfile", dest="snpfile",
                      help="read data from SNPFILE")
parser.add_option("-r", "--refGene", dest="refGene",
                      help="read data from REFGENE")
parser.add_option("-o", "--outputfile", dest="outputfile",
                      help="output data into OUTPUTFILE")
parser.add_option("-n", "--annotation", dest="annotation",
                      help="annotation file including locus or RNA/tRNA(optional)")
parser.add_option("-l", "--addlocus",
                      action="store_false",dest="add",default=True,
                      help="choose to add locus from annotationfile you have(optional)")
parser.add_option("-R", "--RNAtype",
                      action="store_false",dest="RNAtype",default=True,
                      help="choose to add RNA and tRNA type from annotationfile you have(optional)")
parser.add_option("-d", "--description",
                      action="store_false",dest="descr",default=True,
                      help="choose to add description(optional)")
parser.add_option("-i", "--ignoreindel",
                      action="store_false",dest="ignoreindel",default=True,
                      help="selectively choose to ignore indel(optional,default is True)")

parser.add_option("-q", "--quiet",
                      action="store_false", dest="verbose",default=True,
                      help="no standard output")

(options, args) = parser.parse_args()


def readsnpfile(snpfile):
#read postion, ref,alt from snpfile ,then input them into judge fuction
    try:
        f=open(snpfile,"r")
        if options.verbose==True:
            print "reading %s and analyzing..." % options.snpfile
    except:
        print "snpfile is invalid"
        sys.exit(-1)
    f.next()
    for eachline in f:
        try:
            lst=eachline.strip().split("\t")
            pos=lst[1]
            ref=lst[2]
            alt=lst[3]
        except:
            print "snpfile format should be $Chrom,$Postition,$Ref,$Alt"
            sys.exit(-1)
        mut=pos,ref,alt,vsrefGene(options.refGene,pos)
        e.seek(0)
        result=judge(mut,e)
        if result==None:
            pass
        else:
            muttype=result[0]
            ingenepos=result[1]
            altpropos=result[2]
            refpro=result[3]
            altpro=result[4]
            strand=result[5]
            output(eachline,muttype,ingenepos,altpropos,refpro,altpro,strand)
    f.close()
def vsrefGene(refGene,pos):
#position mapped to region 
    try:
        g=open(options.refGene,"r")
        regionlst=g.readlines()
    except:
        print "refGene is invalid"
        sys.exit(-1)
    n=BinarySearch(regionlst,pos)
    return n[0],n[1]

def BinarySearch(regionlst,pos):
#get the pos region
    low = 0
    pos=int(pos)
    height = len(regionlst)-1
    while height-low>1:
        mid = (low+height)/2
        tmplst=regionlst[mid].strip().split("\t")
        midstart=int(tmplst[4])
        midend=int(tmplst[5])
        if midend < pos:
            low = mid
        elif midstart > pos:
            height = mid
        else: 
            return midstart,midend
    tmplst2=regionlst[height].strip().split("\t")
    heightstart=int(tmplst2[4])
    heightend=int(tmplst2[5])
    if heightstart<pos<heightend:
        return heightstart,heightend
    else:
        return "n","n"
def judge(mut,tmpfile):
    '''keys:
        mut = [pos,ref,alt,(anno.start,anno.end)]
        tmpfile: the retrieve seq ouputfile, is "tmp.txt"
    '''
#judge mutation type
    muttype=""
    #tmpfile.seek(0)
    if re.match('[ATGC]',mut[2]) is not None:
        switch=0
        if mut[3][0]=="n":
            muttype="intergenetic"
            ingenepos="."
            altpropos="."
            refpro="."
            altpro="."
            strand="." 
        else:
            ingenepos=int(mut[0])-mut[3][0]+1
            line = tmpfile.readline()
            if line=="":
                return
            tmplst=line.strip().split("\t")
            strand=tmplst[2]
            while switch==0:
                if mut[3][0]==int(tmplst[0]) and mut[3][1]==int(tmplst[1]):
                    mutseq=callmutseq(mut[0],mut[2],mut[3][0],mut[3][1])
                    seqlen=len(mutseq)
                    if tmplst[2]!="+":
                        mutseq=revcom(mutseq)
                        ingenepos=int(mut[0])-mut[3][0]+1
                    mutpro=transform(mutseq)
                    mutlst=compare(mutpro,tmplst[4])
                    muttype=mutlst[0]
                    if muttype=="synonymousSNV":
                        altpropos=(ingenepos-1)//3+1
                        altpro=mutpro[altpropos-1]
                        refpro=altpro
                    else:
                        altpropos=mutlst[1]
                        refpro=mutlst[2]
                        altpro=mutlst[3]
                    switch=1
                    tmpfile.seek(0)
                else:
                    line=tmpfile.readline()
                    if line=="":
                        tmpfile.seek(0)
                        return
                    tmplst=line.strip().split("\t")

    elif options.ignoreindel == True:
        switch=0
        var=mut[2][0]
        indelseq=mut[2][1::]
        num=len(indelseq)
        if mut[3][0]=="n": #if not between any start and end return n
            muttype="intergenetic"
            ingenepos="."
            altpropos="."
            refpro="."
            altpro="."
            strand="."
        else:
            line = tmpfile.readline()
            if line=="":
               tmp.seek(0)
               return
            tmplst=line.strip().split("\t")
            strand=tmplst[2]
            while switch==0:
                if mut[3][0]==int(tmplst[0]) and mut[3][1]==int(tmplst[1]):
                    mutseq=callmutseq(mut[0],mut[2],mut[3][0],mut[3][1])
                    seqlen=len(mutseq)
                    ingenepos=int(mut[0])-mut[3][0]+1
                    mutpro=transform(mutseq)
                    if var=="+":
                        if num%3==0:
                            muttype="nonframeshiftinsertion" 
                            if mutpro=="RNA":  #if the gene codes RNA
                                altpropos="."
                                refpro=mutpro
                                altpro="."
                                if tmplst[2]!="+": #consider if its direction ,which influents ingenepos and mutseq
                                    ingenepos=mut[3][1]-int(mut[0])+1
                                return muttype,ingenepos,altpropos,refpro,altpro,strand
                            if ingenepos % 3==1: #011100000,0 represent normal ,1 represent insertion
                                refseq=mutseq[ingenepos-1:ingenepos+2]
                                mutseq=mutseq[ingenepos-1]+indelseq+mutseq[ingenepos:ingenepos+2]
                                if tmplst[2]!="+": 
                                    ingenepos=mut[3][1]-int(mut[0])+1
                                    refseq=revcom(refseq)
                                    mutseq=revcom(mutseq)
                                if ingenepos!=1: 
                                    altpro=transform(mutseq,0)
                                    refpro=transform(refseq,0)
                                else:  #consider the situation of "ingenepos =1" ,whose referent base is M
                                    altpro=transform(mutseq,0)
                                    refpro=transform(refseq,1)
                            elif ingenepos %3 ==2: #001110000
                                refseq=mutseq[ingenepos-1:ingenepos+1]+mutseq[ingenepos+2]
                                mutseq=mutseq[ingenepos-1:ingenepos+1]+indelseq+mutseq[ingenepos+2]
                                if tmplst[2]!="+":
                                    ingenepos=mut[3][1]-int(mut[0])+1
                                if ingenepos!=2:
                                    altpro=transform(mutseq,0)
                                    refpro=transform(refseq,0)
                                else:
                                    altpro=transform(mutseq,0)
                                    refpro=transform(refseq,1)
                            else:#111000000
                                refpro="."
                                if tmplst[2]!="+":
                                    altpro=transform(revcom(indelseq),0)
                                else:
                                    altpro=transform(indelseq,0)
                    #compute at end because it derive defferent ingenepos,if its amount > 1 ,list it 
                            altpropos=(ingenepos-1)//3+1
                            t=len(altpro)
                            if t>1:
                                i=1
                                while i<t:
                                    altpropostmp=str(altpropos)+","+str(altpropos+1)
                                    i+=1
                                altpropos=altpropostmp
                        else:
                            muttype="frameshiftinsertion" 
                            if mutpro=="RNA":
                                altpropos="."
                                refpro=mutpro
                                altpro="."
                                if tmplst[2]!="+": 
                                    ingenepos=mut[3][1]-int(mut[0])+1
                                return muttype,ingenepos,altpropos,refpro,altpro,strand
                            altpropos="."
                            refpro="."
                            altpro="."
                    elif var=="-":
                        if int(mut[0])+num<int(mut[3][1])-2:
                            if num%3==0:
                                muttype="nonframeshiftdeletion"
                                if mutpro=="RNA":
                                    altpropos="."
                                    refpro=mutpro
                                    altpro="."
                                    if tmplst[2]!="+": 
                                        ingenepos=mut[3][1]-int(mut[0])+1
                                    return muttype,ingenepos,altpropos,refpro,altpro,strand
                                if ingenepos % 3==1:
                                    refseq=mutseq[ingenepos-1:ingenepos+num+2]
                                    mutseq=mutseq[ingenepos-1]+mutseq[ingenepos+num:ingenepos+num+2]
                                    if tmplst[2]!="+":
                                        ingenepos=mut[3][1]-int(mut[0])+1
                                        refseq=revcom(refseq)
                                        mutseq=revcom(mutseq)
                                    if ingenepos!=1:
                                        altpro=transform(mutseq,0)
                                        refpro=transform(refseq,0)
                                    else:
                                        altpro=transform(mutseq,0)
                                        refpro=transform(refseq,1)
                                elif ingenepos %3 ==2:
                                    refseq=mutseq[ingenepos-1:ingenepos+num+2]
                                    mutseq=mutseq[ingenepos-1:ingenepos+1]+mutseq[ingenepos+num+2]
                                    if tmplst[2]!="+":
                                        ingenepos=mut[3][1]-int(mut[0])+1
                                    if ingenepos!=1:
                                        altpro=transform(mutseq,0)
                                        refpro=transform(refseq,0)
                                    else:
                                        altpro=transform(mutseq,0)
                                        refpro=transform(refseq,1)
                                else:
                                    if tmplst[2]!="+":
                                        refpro=transform(revcom(indelseq),0)
                                    else:
                                        refpro=transform(indelseq,0)
                                    altpro="."
                                altpropos=(ingenepos-1)//3+1
                                t=len(altpro)
                                if t>1:
                                    i=1
                                    while i<t:
                                        altpropostmp=str(altpropos)+","+str(altpropos+1)
                                        i+=1
                                    altpropos=altpropostmp

                            else:
                                muttype="frameshiftdeletion"
                                if mutpro=="RNA":
                                    altpropos="."
                                    refpro=mutpro
                                    altpro="."
                                    if tmplst[2]!="+": 
                                        ingenepos=mut[3][1]-int(mut[0])+1
                                    return muttype,ingenepos,altpropos,refpro,altpro,strand
                                altpropos="."
                                refpro="."
                                altpro="."
                        else:
                            muttype="stoploss"
                            if mutpro=="RNA":
                                altpropos="."
                                refpro=mutpro
                                altpro="."
                                if tmplst[2]!="+": 
                                    ingenepos=mut[3][1]-int(mut[0])+1
                                return muttype,ingenepos,altpropos,refpro,altpro,strand
                            altpropos="."
                            refpro="*"
                            altpro="LOSS"
                    switch=1
                    tmpfile.seek(0)
                else:
                    line=tmpfile.readline()
                    if line=="":
                        tmpfile.seek(0)
                        return
                    tmplst=line.strip().split("\t")
                    strand=tmplst[2]
    else:
        return
#ingenepos=alterent base position / length of transcripted seq and altpropos=alterant protein position / length of translated protein
    if ingenepos !=".":
        ingenepos=str(ingenepos)+"/"+str(seqlen)
        if altpropos !=".":
            altpropos=str(altpropos)+"/"+str(seqlen/3)
    return muttype,ingenepos,altpropos,refpro,altpro,strand    
def callmutseq(pos,alt,txstart,txend,fafile=options.reffa):
#call the mutseq from snp and fafile
    g=open(fafile,'r')
    allseq=g.readlines()[1]
    #print type(txstart),type(pos)
    ipos=int(pos)
    if re.match('[ATGC]',alt):
        seq=allseq[txstart-1:ipos-1]+alt+allseq[ipos:txend]
    else:
        seq=allseq[txstart-1:txend]
    return seq
    g.close()

def compare(n1,n2):
#compare two protein sequence
    i=0
    m=len(n1)
    if n2=="RNA":
        muttype="nonsynonymousSNV"
        altpropos="."
        refpro="RNA"
        altpro="RNA" 
        return muttype,altpropos,refpro,altpro
    while i<m:
        if n1[i]!=n2[i]:
            if n1[i]=="*":
                muttype="stopgain"
                altpropos=i
                refpro=n2[i]
                altpro=n1[i]
            elif n2[i]=="*":
                muttype="stoploss"
                altpropos=i
                refpro=n2[i]
                altpro=n1[i]
            else:
                muttype="nonsynonymousSNV"
                altpropos=i
                refpro=n2[i]
                altpro=n1[i]
            i=m+1
        i=i+1
    if i==m:
        muttype="synonymousSNV"
        altpropos="."
        refpro="."
        altpro="."
    
    return muttype,altpropos,refpro,altpro
def outputheader(outfile):
#output the header and snpfile header is used
    f=open(outfile,'w+')
    line="Chrom"+"\t"+"Position"+"\t"+"Ref"+"\t"+"Alt"+"\t"+"Muttype"+"\t"+"Strand"+"\t"+"Ingenepos"+"\t"+"Altpropos"+"\t"+"Refpro"+"\t"+"Altpro"+"\n"
    f.write(line)
    f.close()
def output(snpline,muttype,ingenepos,altpropos,refpro,altpro,strand,outfile=options.outputfile):
    f=open(outfile,'a')
    snplst=snpline.strip().split("\t")
    if re.match('[+-][ATGC]+',snplst[3]): #if strand is "-", we should transform the indel seq into reverse complementation to read  more easily.
        typ=snplst[3][0]
        indel=snplst[3][1::]
        if strand=="-":
            snplst[3]=typ+revcom(indel)
    line=snplst[0]+"\t"+snplst[1]+"\t"+snplst[2]+"\t"+snplst[3]+"\t"+muttype+"\t"+strand+"\t"+str(ingenepos)+"\t"+str(altpropos)+"\t"+refpro+"\t"+altpro+"\n"
    f.write(line)

#option part: addlocus and tRNA,rRNA et.
def addlocus(rawoutputfile,annfile):
    '''selectively add locus or tRNA,rRNA types to the result.

    if you have annotationfile you can add locus where the alt gene falls.
   
    args:
        outputfile: the outputfile you have added locus
        rawoutputfile: the file you have finished previous work.
        annfile: the annotationfile with header, eachline has eight tab-delimited
            coloumn:$species,$start,$end,$strand,$type,$locus,$productId,
            $description. For example, Ecoli,38,1451,+,protein,Ecolc_0001,
            ACA75689.1, chromosomal replication initiator protein DnaA. Only 
            important thing is $start,$end,$type,$locus, the type(including 
            tRNA or rRNA) and locus both are selectively.
         
    '''
    os.rename(rawoutputfile,"tmp.txt2")
    try:
        a=open(annfile,'r')
        regionlst=a.readlines()
    except:
        print "annotation file is invalid"
        parser.print_help()
        sys.exit(-1)
    f=open("tmp.txt2",'r')
    line=f.readline().strip().split("\t")
    g=open(rawoutputfile,'w+')
    if options.add == False and options.RNAtype == False:
        line=line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]+"\t"+line[4]+"\t"+"ProductType"+"\t"+"Locus"+"\t"+line[5]+"\t"+line[6]+"\t"+line[7]+"\t"+line[8]+"\t"+line[9]+"\n"
    elif options.add == False and options.RNAtype == True:         
        line=line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]+"\t"+line[4]+"\t"+"Locus"+"\t"+line[5]+"\t"+line[6]+"\t"+line[7]+"\t"+line[8]+"\t"+line[9]+"\n"
    elif options.add == True and options.RNAtype == False:
        line=line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]+"\t"+line[4]+"\t"+"ProductType"+"\t"+line[5]+"\t"+line[6]+"\t"+line[7]+"\t"+line[8]+"\t"+line[9]+"\n"
    if options.descr == False:
        line=line.strip()+"\t"+"Description"+"\n"
    g.write(line)
    for eachline in f:
        line=eachline.strip()
        lst=eachline.strip().split("\t")
        pos=lst[1]
        mut=BinarySearch2(regionlst,pos) #start,end,producttype(selectivel,locus,proid.discr
        if options.add == False and options.RNAtype == True:
            line=lst[0]+"\t"+lst[1]+"\t"+lst[2]+"\t"+lst[3]+"\t"+lst[4]+"\t"+mut[3]+"\t"+lst[5]+"\t"+lst[6]+"\t"+lst[7]+"\t"+lst[8]+"\t"+lst[9]+"\n"
        elif options.RNAtype == False and options.add == False:
                if re.match("^.*RNA",mut[2]) is not None:
                    line=lst[0]+"\t"+lst[1]+"\t"+lst[2]+"\t"+lst[3]+"\t"+"nonsynonymousSNV"+"\t"+mut[2]+"\t"+mut[3]+"\t"+lst[5]+"\t"+lst[6]+"\t"+"."+"\t"+"."+"\t"+"."+"\n"
                else:
                    line=lst[0]+"\t"+lst[1]+"\t"+lst[2]+"\t"+lst[3]+"\t"+lst[4]+"\t"+mut[2]+"\t"+mut[3]+"\t"+lst[5]+"\t"+lst[6]+"\t"+lst[7]+"\t"+lst[8]+"\t"+lst[9]+"\n"
        elif options.RNAtype == False and options.add == True:
            if re.match("^.*RNA",mut[2]) is not None:
                line=lst[0]+"\t"+lst[1]+"\t"+lst[2]+"\t"+lst[3]+"\t"+"nonsynonymousSNV"+"\t"+mut[3]+"\t"+lst[5]+"\t"+lst[6]+"\t"+"."+"\t"+"."+"\t"+"."+"\n"
            else:
                line=lst[0]+"\t"+lst[1]+"\t"+lst[2]+"\t"+lst[3]+"\t"+lst[4]+"\t"+mut[3]+"\t"+lst[5]+"\t"+lst[6]+"\t"+lst[7]+"\t"+lst[8]+"\t"+lst[9]+"\n"
        if options.descr == False:
            line=line.strip()+"\t"+mut[-1]+"\n"
        g.write(line)
    f.close()
    os.remove("tmp.txt2")
def BinarySearch2(regionlst,pos):
    '''search pos in regionlst
    
    args:
        regionlst: a list of line in annotation file.
        pos: snp positon.

    returns:
        midstart:mid means result of binary search, a line number, whose number of start column and 
            end column fit the snp position, its start column number is midstart.
        midend:its end column number is midend.
        locus:its locus number
        proid:its proid column data
        discr:its discription column data
    '''
    low = 0
    pos=int(pos)
    height = len(regionlst)-1
    while height-low>1:
        mid = (low+height)/2
        tmplst=regionlst[mid].strip().split("\t")
        midstart=int(tmplst[1])
        midend=int(tmplst[2])
        producttype=tmplst[4]
        locus=tmplst[5]
        proid=tmplst[6]
        discr=tmplst[7]
        if midend < pos:
            low = mid
        elif midstart > pos:
            height = mid
        else:
            return midstart,midend,producttype,locus,proid,discr
    tmplst2=regionlst[height].strip().split("\t")
    heightstart=int(tmplst2[1])
    heightend=int(tmplst2[2])
    producttype=tmplst[4]
    locus=tmplst2[5]
    proid=tmplst2[6]
    discr=tmplst2[7]
    if heightstart<pos<heightend:
        return heightstart,heightend,producttype,locus,proid,discr
    else:
        return "n","n",".",".",".","."
#def addRNAtype(outputfile,annfile):
   
    
   
if __name__ == "__main__":
    if not(options.snpfile and options.refGene and options.outputfile and options.reffa):
        parser.print_help()
        sys.exit(0)
    if options.verbose==True:
        print "retrieve ref region seq from %s and %s into tmptxt ..." % (options.refGene,options.reffa)
    retrieveseq(options.reffa,options.refGene,"tmp.txt")#retrieve gene seq from ref.fa into tmp.txt
    e=open("tmp.txt","r")#open tmp.txt as a global variable
    outputheader(options.outputfile)#output header
    readsnpfile(options.snpfile)
    e.close()
    os.remove("tmp.txt")
    if options.add == False or options.annotation == False or options.descr == False:
        if options.add == False:
            print "adding locus..."
        if options.RNAtype ==  False:
            print "adding ProductType..."
        if options.descr == False:
            print "adding Description..."
        addlocus(options.outputfile,options.annotation)
    if options.verbose == True:
        print "call types of mutations successfully!"
