import re
path = "/home/lmartinezg/Documents/Laura/appris_clinvar/ForRevision/"

'''exons_i = ["prin_appris_mane_longest_no_overlap", "alt_appris_mane_longest_no_overlap"]
exons_ii = ["prin_mane_unique_vs_l", "alt_appris_unique"]
exons_ii_i =["prin_appris_unique", "alt_appris_unique"]
exons_iii = ['prin_appris_unique', "alt_mane_unique_l"]
exons_iv = ["prin_longest_unique", "alt_longest_unique"]
exons_v = ['prin_appris_mane', "prin_appris_unique_vs_mane", "prin_mane_unique_vs_appris"] '''
#exons_v = ['prin_appris_unique_vs_mane_FI', 'prin_appris_unique_vs_mane',
#'prin_mane_unique_vs_appris', 'prin_mane_unique_vs_appris_FI']

exons_5 = ['set_5.mane_alt_no_overlap.len.undup', 'set_5.mane_no_overlap']
exons_3 = ['set_3.appris_alt_no_overlap.len.undup', 'set_3.principal_no_overlap']
exons_7 = ['set_7.longest_alt_no_overlap.len.undup', 'set_7.longest_no_overlap']
exons_list = [exons_3, exons_5, exons_7]
#exons_list = [exons_i, exons_ii, exons_iii, exons_iv, exons_v]
for exons_type in exons_list:
    print(exons_type)
    for exon_type in exons_type:
        print (exon_type)
        if exons_type == exons_3:
            file_name = path + "/results_3/" +  exon_type + "_wvep.gff"
        elif exons_type == exons_5:
            file_name = path + "/results_5/" +  exon_type + "_wvep.gff"
        elif exons_type == exons_7:
            file_name = path + "/results_7/" +  exon_type + "_wvep.gff"

        '''elif exons_type == exons_iv:
            file_name = path + "/results/iv/gencode.v37." +  exon_type + "_wvep.txt"
        elif exons_type == exons_v:
            file_name = path + "/results/v/gencode.v37." +  exon_type + "_wvep.txt"
        elif exons_type == exons_ii_i:
            file_name = path + "/results/ii_i/gencode.v37." +  exon_type + "_wvep.gff" '''

        fhandle = open(file_name)
        if exons_type == exons_3:
            out_name = path + "/results_3/" +  exon_type + "_wvep_merged.txt"
        elif exons_type == exons_5:
            out_name = path + "/results_5/" +  exon_type + "_wvep_merged.txt"
        elif exons_type == exons_7:
            out_name = path + "/results_7/" +  exon_type + "_wvep_merged.txt"
        '''elif exons_type == exons_iv:
            out_name = path + "/results/iv/gencode.v37." + exon_type + "_wvep_merged.txt"
        elif exons_type == exons_v:
            out_name = path + "/results/v/gencode.v37." + exon_type + "_wvep_merged.txt"
        elif exons_type == exons_ii_i:
            out_name = path + "/results/ii_i/gencode.v37." + exon_type + "_wvep_merged.txt" '''

        outhandle = open(out_name, "w")
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
        outhandle.close()
        fhandle.close()