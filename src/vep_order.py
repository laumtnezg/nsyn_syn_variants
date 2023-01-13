import re
path = "/home/lmartinezg/Documents/Laura/exons_dnds/"
exon = 'Primate_derived_upstream_exons_frame' 
file_name = path + exon + ".csv"

list_regions = []
dict_regions = {}
dict_chromosomes = {}
with open (file_name,'r') as fhandle:
    for line in fhandle:
        if line.startswith('Chromosome'):
            continue
        line = line.replace("\n", '')
        strand = line.split("\t")[3]
        chrom = line.split("\t")[0].split("chr")[1]
        frame = line.split("\t")[4]
        posi = chrom + "-" +str(line.split("\t")[1]) + "-" + str(line.split("\t")[2])
        list_regions.append(posi)
        dict_regions[posi] = [strand, frame]
        if chrom not in dict_chromosomes:
            dict_chromosomes[chrom] = []
        dict_chromosomes[chrom].append(posi)

out_name = path + "/results/" +  exon + "_wvep.gff"
out_no_name = path + "/results/" +  exon + "_notvep.gff"
out_no_annot = path + "/results/" +  exon + "_notannot.gff"
file_vep_name = "/home/lmartinezg/Documents/Laura/VEP/results/results_vep/all_results_filtered.txt"



with open(out_no_name, 'w') as outnohandle, open(file_vep_name) as vhandle, open(out_name, 'w') as outhandle:
    for line in vhandle:
        line = line.replace("\n","")
        if line.startswith("#"):
            continue
        pos = line.split("\t")[1].split(":")[1]
        chromosome = line.split("\t")[1].split(":")[0]
        if chromosome in dict_chromosomes:
            is_exon = False
            for exon in dict_chromosomes[chromosome]:
                start = exon.split("-")[1]
                end = exon.split("-")[2]
                if (start <= pos) and (pos <= end):
                    is_exon = True
                    outhandle.write(">" + line.split("\t")[3] + "\t" + exon + "\t" +  dict_regions[exon][0] +  "\t" + dict_regions[exon][1] + "\n")
                    if exon in list_regions:
                        list_regions.remove(exon)
                    outhandle.write(line + "\n")
            if is_exon == False:
                outnohandle.write(line + "\n")
        else:
            outnohandle.write(line + "\n")
print(len(list_regions))

with open(out_no_annot, 'w') as outno:
    for exon in list_regions:
        outno.write(exon +  "\t" + dict_regions[exon][0] +  "\t" + dict_regions[exon][1] + "\n")