path = "/home/lmartinezg/Documents/Laura/exons_dnds/"
exon = 'Primate_derived_upstream_exons_frame' 

file_name = path + "/results/" +  exon + "_wvep.gff"
out_name = path + "/results/" +  exon + "_wvep_merged.txt"
with open(file_name, 'r') as fhandle, open(out_name, 'w') as outhandle:
    exons = []
    for line in fhandle:
        line = line.replace("\n", "")
        if line.startswith(">"):
            exon = line.split("\t")[1]
            if exon in exons:
                continue
            else:
                outhandle.write(line + "\n")
                exons.append(exon)
        else:
            outhandle.write(line + "\n")