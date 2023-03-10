import sys
#path = sys.argv[1]
#file_vep_name = sys.argv[2]
#out_name = sys.argv[3]
#out_no_name = sys.argv[4]
path = '/home/lmartinezg/Documents/Laura/exons_dnds/'
# load different files if youre doing a comparison analysis, inside each type theres usually principal and alt
NC_upstream = ['nc_upstream']
exons_list = [NC_upstream]
for exons_type in exons_list:
    print(exons_type)
    for exon_type in exons_type:
        print (exon_type)
        file_name = path + exon_type + ".tsv"
        fhandle = open(file_name)
        list_regions = []
        dict_regions = {}
        dict_chromosomes = {}
        for line in fhandle:
            line = line.replace("\n", '')
            try:
                strand = line.split("\t")[3]
            except:
                print(line)
                exit()
            chrom = line.split("\t")[0].split("chr")[1]
            posi = chrom + "-" +str(line.split("\t")[1]) + "-" + str(line.split("\t")[2])
            list_regions.append(posi)
            dict_regions[posi] = [strand]
            if chrom not in dict_chromosomes:
                dict_chromosomes[chrom] = []
            dict_chromosomes[chrom].append(posi)
        fhandle.close()
        #remove when in pipeline
        out_name = path + 'results_' + exon_type + 'wvep.gff'
        out_no_name = path + 'results_' + exon_type + 'novep.gff'
        outhandle = open(out_name, "w")
        file_vep_name = "/home/lmartinezg/Documents/Laura/VEP/results/results_vep/all_results_filtered.txt"
        with open(out_no_name, 'w') as outnohandle, open(file_vep_name) as vhandle:
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
                            outhandle.write(">" + line.split("\t")[3] + "\t" + exon + "\t" +  dict_regions[exon][0] + "\n")
                            outhandle.write(line + "\n")
                    if is_exon == False:
                        outnohandle.write(line + "\n")
                else:
                    outnohandle.write(line + "\n")
        outhandle.close()