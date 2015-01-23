#copyright Jesse van Dam 

import vcf
import inspect
import re

accessionMap = {}
mappingfile = open("mappingfile2accession.csv",'r')
first = True
for line in mappingfile:
    if first:
        first = False
        continue
    temp = line.strip().split('\t')
    if len(temp) == 2:
        accessionMap[temp[0]] = temp[1]
mappingfile.close()

rdf = open("snps.n3","w")
vcf_reader = vcf.Reader(open('tomatoExample.vcf', 'r'))

prefix = "http://pbr.wur.nl/VCF/"
refVersion = "Tomatov1"
isType = "<http://www.w3.org/1999/02/22-rdf-syntax-ns#type>"

uniqId = {}
chromdone = {}

refId = "<" + prefix + refVersion + ">"
rdf.write(refId + " " + isType + " <" + prefix + "Reference> .\n")
rdf.write(refId + " <http://www.w3.org/1999/02/22-rdf-syntax-ns#comment> \"tomato reference assembly version 2.40(aka 2.31)\" .\n")

variant={}
subClassDef={}

firstRecord = True
for record in vcf_reader:
    SNPidRaw = "http://pbr.wur.nl/VCF/" + str(refVersion) + "/chrom" + str(record.CHROM) + "/" + str(record.POS) + "/"
    
    if not SNPidRaw in uniqId:
        uniqId[SNPidRaw] = 0;
    uniqId[SNPidRaw] = uniqId[SNPidRaw] + 1
    SNPidRaw = SNPidRaw + str(uniqId[SNPidRaw])
    SNPid = "<" + SNPidRaw + ">"
    rdf.write(SNPid + " " + isType + " " + "<" + prefix + "SNP> . \n")
    rdf.write(SNPid + " <" + prefix + "refSNP> \"" + record.REF + "\" . \n")
    skip = True
    for allel in record.alleles:
        if skip:
            skip = False
            continue
        rdf.write(SNPid + " <" + prefix + "allelySNP> \"" + str(allel) + "\" . \n")
    #double position object are automaticly merged
    posId = "<" + prefix + str(refVersion) + "/chrom" + str(record.CHROM) + "/" + str(record.POS) + ">"
    rdf.write(SNPid + " <" + prefix + "position> " + posId + " . \n")
    rdf.write(posId + " " + isType + " " + "<" + prefix + "Position> . \n")
    rdf.write(posId + " <" + prefix + "location> \"" + str(record.POS) + "\"^^<http://www.w3.org/2001/XMLSchema#integer> . \n")
    chromId = "<" + prefix + str(refVersion) + "/chrom" + str(record.CHROM) + ">"
    rdf.write(posId + " <" + prefix + "chromoson> " + chromId + " . \n")
    if not chromId in chromdone:
        chromdone[chromId] = 1
        rdf.write(chromId + " " + isType + " " + "<" + prefix + "Chromoson> . \n")
        scaffold = "<http://pbr.wur.nl/SCAFFOLD#SL2.31ch"
        if len(str(record.CHROM)) == 1:
            scaffold = scaffold + "0"
        #<http://pbr.wur.nl/SCAFFOLD#SL2.31ch00>
        scaffold = scaffold + str(record.CHROM) + ">" 
        rdf.write(chromId + " <" + prefix + "scaffold> " + scaffold + " . \n")
        rdf.write(chromId + " <" + prefix + "reference> " + refId + " . \n")
    rdf.write("")
    for sample in record.samples:
        sampleId = re.search(".*/(.*)",sample.sample).group(1)
        sampleIri = "<" + prefix + "sample/" + sampleId + ">"
        SampleSNPId = SNPidRaw + sampleId
        SampleSNPIri= "<" + SampleSNPId + ">"
        if firstRecord:
            rdf.write(sampleIri + " " + isType + " " + "<" + prefix + "Sample> . \n")
            rdf.write(sampleIri + " <" + prefix + "accession> <http://purl.org/cgngenis/accenumb/" + accessionMap[sampleId] +"> . \n")
        rdf.write(SNPid + " <" + prefix + "genotype> " + SampleSNPIri + " . \n")
        rdf.write(SampleSNPIri + " " + isType + " " + "<" + prefix + "SNPGenotype> . \n")
        rdf.write(SampleSNPIri + " <" + prefix + "sample> " + sampleIri + " . \n")
        rdf.write(SampleSNPIri + " <" + prefix + "phased> \"" + str(sample.phased).lower() + "\"^^<http://www.w3.org/2001/XMLSchema#boolean> . \n")
       
        #print(record)
        if sample.gt_type:
            count = 1;
            for allel in sample.gt_alleles:
                allelIri = "<" + SampleSNPId + "/allel_" + str(count) + ">"
                rdf.write(SampleSNPIri + " <" + prefix + "allel> " + allelIri + " . \n")
                rdf.write(allelIri + " " + isType + " " + "<" + prefix + "Allel> . \n")
                rdf.write(allelIri + " <" + prefix + "allelNumber> \"" + str(count) + "\"^^<http://www.w3.org/2001/XMLSchema#integer> . \n")
                rdf.write(allelIri + " <" + prefix + "variantAllel> \"" + str(record.alleles[int(allel)]) + "\" . \n")
                count = count + 1
      # rdf.write(SNPid + " <" + prefix + "genotype> " + SampleSNPId + " . \n")
    firstRecord = False
    count = 0
    for info in record.INFO['EFF']:
        match = re.search("(.*)\((.*\))",record.INFO['EFF'][0])
        variantType = match.group(1)
        data = match.group(2).split('|')
        #print(match.group(1))
        #if variantType == "missense_variant":
        
       # <http://pbr.wur.nl/GENE#
        variantClassType = "<" + prefix + variantType + "Annotation>"
        if not variantType in subClassDef:
            subClassDef[variantType] = False
            rdf.write(variantClassType + " <http://www.w3.org/2000/01/rdf-schema#subClassOf> <" + prefix + "Annotation> . \n")
            rdf.write(variantClassType + " " + isType + " <http://www.w3.org/2002/07/owl#Class> . \n")
        infoIRI = "<" + SNPidRaw + "/annot_" + str(count) + ">"
        rdf.write(infoIRI + " " + isType + " " + variantClassType +" . \n")
        #rdf.write(SampleSNPIri + " <" + prefix + "allel> " + allelIri + " . \n")
        rdf.write(SNPid + " <" + prefix + "annotation> " + infoIRI + " . \n")
        print(variantType + str(data) + data[3])
        geneName = data[8]
        if geneName != "":
            rdf.write(infoIRI + " <" + prefix + "gene> <http://pbr.wur.nl/GENE#" + geneName + "> . \n")
        hgvs_P = data[3]
        if hgvs_P != "":
            rdf.write(infoIRI + " <" + prefix + "hgvs_P> \"" + hgvs_P + "\" . \n")
        rdf.write(infoIRI + " <" + prefix + "importance> \"" + data[0] + "\" . \n")
        count = count + 1
       
  
    
    #print(dir(record.INFO['EFF'][0]))
    SNPidRaw
    rdf.write(SNPid + " <" + prefix + "refSNP> \"" + record.REF + "\" . \n")
   # break
    
rdf.close()

print("done")
#for record in vcf_reader:
#    print(record)
#    fields = dir(record)
#    for field,value in inspect.getmembers(record):
#        if field.startswith("__"):
#            continue;
#        print(field + " : " + str(value))
#    break

#ALT : [C]
#CHROM : 0
#FILTER : None
#FORMAT : GT:AD:DP:GQ:OG:PL
#ID : None
#INFO : {'NumGenotypesChanged': 0, 'BaseQRankSum': 0.874, 'MLEAC': [103], 'AN': 6, 'R2': 1.0, 'FS': 26.136, 'MLEAF': [0.725], 'HaplotypeScore': 2.4948, 'Dels': 0.0, 'MQRankSum': 2.235, 'MQ': 58.12, 'MQ0': 0, 'DP': 2520, 'AF': [0.667], 'EFF': ['intergenic_region(MODIFIER||||||||||1)'], 'ReadPosRankSum': -0.107, 'InbreedingCoeff': 0.9768, 'AC': [4], 'QD': 34.25}
#POS : 100944
#QUAL : 72819.21
#REF : T
#_sample_indexes : {'/ifshk5/PC_PA_EU/PMO/Tomato_reseq/01.BWA/SZAXPI009286-74': 2, '/ifshk5/PC_PA_EU/PMO/Tomato_reseq/01.BWA/SZAXPI009285-62': 1, '/ifshk5/PC_PA_EU/PMO/Tomato_reseq/01.BWA/SZAXPI009284-57': 0}
#aaf : [0.6666666666666666]
#add_filter : <bound method _Record.add_filter of <vcf.model._Record object at 0x7ffd3e401278>>
#add_format : <bound method _Record.add_format of <vcf.model._Record object at 0x7ffd3e401278>>
#add_info : <bound method _Record.add_info of <vcf.model._Record object at 0x7ffd3e401278>>
#alleles : ['T', C]
#call_rate : 1.0
#end : 100944
#genotype : <bound method _Record.genotype of <vcf.model._Record object at 0x7ffd3e401278>>
#get_hets : <bound method _Record.get_hets of <vcf.model._Record object at 0x7ffd3e401278>>
#get_hom_alts : <bound method _Record.get_hom_alts of <vcf.model._Record object at 0x7ffd3e401278>>
#get_hom_refs : <bound method _Record.get_hom_refs of <vcf.model._Record object at 0x7ffd3e401278>>
#get_unknowns : <bound method _Record.get_unknowns of <vcf.model._Record object at 0x7ffd3e401278>>
#heterozygosity : 0.4444444444444444
#is_deletion : False
##is_indel : False
#is_monomorphic : False
#is_snp : True
#is_sv : False
#is_sv_precise : False
#is_transition : True
#nucl_diversity : 0.5333333333333333
#num_called : 3
#num_het : 0
#num_hom_alt : 2
#num_hom_ref : 1
#num_unknown : 0
#samples : [Call(sample=/ifshk5/PC_PA_EU/PMO/Tomato_reseq/01.BWA/SZAXPI009284-57, CallData(GT=1|1, AD=[0, 32], DP=32, GQ=60, OG=None, PL=[1222, 96, 0])), Call(sample=/ifshk5/PC_PA_EU/PMO/Tomato_reseq/01.BWA/SZAXPI009285-62, CallData(GT=0|0, AD=[47, 0], DP=47, GQ=60, OG=None, PL=[0, 141, 1827])), Call(sample=/ifshk5/PC_PA_EU/PMO/Tomato_reseq/01.BWA/SZAXPI009286-74, CallData(GT=1|1, AD=[0, 31], DP=31, GQ=60, OG=None, PL=[1260, 93, 0]))]
#start : 100943
#sv_end : None
#var_subtype : ts
#var_type : snp
