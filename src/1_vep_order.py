import re
path = "/home/lmartinezg/Documents/Laura/appris_clinvar/ForRevision/"

'''exons_i = ["prin_appris_mane_longest_no_overlap", "alt_appris_mane_longest_no_overlap"]
exons_ii = ["prin_mane_unique_vs_l", "alt_appris_unique"]
exons_iii = ['prin_appris_unique', "alt_mane_unique_l"]
exons_ii_i =["prin_appris_unique", "alt_appris_unique"]
exons_iv = ["prin_longest_unique", "alt_longest_unique"]
exons_v = ['prin_appris_mane', "prin_appris_unique_vs_mane", "prin_mane_unique_vs_appris"] '''

exons_5 = ['set_5.mane_alt_no_overlap.len.undup', 'set_5.mane_no_overlap']
exons_3 = ['set_3.appris_alt_no_overlap.len.undup', 'set_3.principal_no_overlap']
exons_7 = ['set_7.longest_alt_no_overlap.len.undup', 'set_7.longest_no_overlap']
#exons_v = ['prin_appris_unique_vs_mane_FI', 'prin_appris_unique_vs_mane',
#'prin_mane_unique_vs_appris', 'prin_mane_unique_vs_appris_FI']

exons_list = [exons_3, exons_5, exons_7]
#exons_list = [exons_i, exons_ii, exons_iii, exons_iv, exons_v]
for exons_type in exons_list:
    print(exons_type)
    for exon_type in exons_type:
        print (exon_type)
        file_name = path + exon_type + ".txt"

        ''' if exons_type == exons_3:
            file_name = path + exon_type + "txt"
        elif exons_type == exons_5:
            file_name = path +  exon_type + ".txt"
        elif exons_type == exons_7:
            file_name = path + exon_type + ".txt" '''
        '''elif exons_type == exons_iv:
            file_name = path + exon_type + ".txt"
        elif exons_type == exons_v:
            file_name = path + exon_type + ".txt"
        elif exons_type == exons_ii_i:
            file_name = path + exon_type + ".txt" '''

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
        if exons_type == exons_3:
            out_name = path + "/results_3/" +  exon_type + "_wvep.gff"
        elif exons_type == exons_5:
            out_name = path + "/results_5/" +  exon_type + "_wvep.gff"
        elif exons_type == exons_7:
            out_name = path + "/results_7/" +  exon_type + "_wvep.gff"
        
        '''elif exons_type == exons_iv:
            out_name = path + "/results/iv/gencode.v37." +  exon_type + "_wvep.gff"
        elif exons_type == exons_v:
            out_name = path + "/results/v/gencode.v37." +  exon_type + "_wvep.gff"
        elif exons_type == exons_ii_i:
            out_name = path + "/results/ii_i/gencode.v37." +  exon_type + "_wvep.gff" '''

        outhandle = open(out_name, "w")
        if exons_type == exons_3:
            out_no_name = path + "/results_3/" +  exon_type + "_notvep.gff"
        elif exons_type == exons_5:
            out_no_name = path + "/results_5/" +  exon_type + "_notvep.gff"
        elif exons_type == exons_7:
            out_no_name = path + "/results_7/" +  exon_type + "_notvep.gff"
        
        ''' elif exons_type== exons_iv:
            out_no_name = path + "/results/iv/gencode.v37." +  exon_type + "_notvep.gff"
        elif exons_type == exons_v:
            out_no_name = path + "/results/v/gencode.v37." +  exon_type + "_notvep.gff"
        elif exons_type == exons_ii_i:
            out_no_name = path + "/results/ii_i/gencode.v37." +  exon_type + "_notvep.gff" '''

        #outnohandle = open(out_no_name, "w")
        file_vep_name = "/home/lmartinezg/Documents/Laura/VEP/results/results_vep/all_results_filtered.txt"
        #vhandle = open(file_vep_name)

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
        ''' if line.split("\t")[3] in dict_genes:
                is_exon = False
                for exon in dict_genes[line.split("\t")[3]]:
                    start = exon.split("-")[1]
                    end = exon.split("-")[2]
                    if (start <= pos) and (pos <= end):
                        is_exon = True
                        #outhandle.write('>' + exon + '\n')
                        outhandle.write(">" + dict_regions[exon][0] + "\t" + dict_regions[exon][1] + 
                                        "\t" + exon + "\t" + dict_regions[exon][2] + "\n")
                        outhandle.write(line + "\n")
                if is_exon == False:
                    outnohandle.write(line + "\n")
            else:
                outnohandle.write(line + "\n")
        vhandle.close()
        outhandle.close()
        outnohandle.close() '''