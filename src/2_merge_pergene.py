import sys
#path = sys.argv[1]
# file_name = sys.argv[2]
# out_name = sys.argv[3]
path = '/home/lmartinezg/Documents/Laura/exons_dnds/'
# load different files if youre doing a comparison analysis, inside each type theres usually principal and alt
NC_upstream = ['nc_upstream']
exons_list = [NC_upstream]
for exons_type in exons_list:
    print(exons_type)
    for exon_type in exons_type:
        file_name = path + "/results_/" +  exon_type + "_wvep.gff"
        out_name = path + "/results_/" +  exon_type + "_wvep_merged.txt"
        exons = []
        with open(file_name, 'r') as fhandle, open(out_name, 'w') as outhandle:
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
