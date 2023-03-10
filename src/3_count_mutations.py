import re
import sys
# path = sys.argv[1]
# fname = sys.argv[2]
# outname = sys.argv[3]
# out_name = sys.argv[4]
path = '/home/lmartinezg/Documents/Laura/exons_dnds/'
# load different files if youre doing a comparison analysis, inside each type theres usually principal and alt
NC_upstream = ['nc_upstream']
exons_list = [NC_upstream]
for exons_type in exons_list:
    print(exons_type)
    for exon_type in exons_type:
        print(exon_type)
        fname = path + '/results/' + exon_type + "_wvep_merged.txt"
        outname = path + '/results/' + exon_type + "_wvep_filtered.gff"
        with open(fname) as khandle, open(outname, 'w') as ohandle:
            print(fname)
            print(outname)
            variations = []
            for line in khandle:
                line = line.replace('\n', '')
                if line.startswith('>'):
                    ohandle.write(line + '\n')
                else:
                    if line.split("\t")[0] not in variations:
                        ohandle.write(line + '\n')
                        variations.append(line.split("\t")[0])
                    else:
                        continue
    rare_syn = 0
    common_syn = 0
    rare_missense = 0
    common_missense = 0
    rare_high = 0
    common_high = 0
    frequency = None
    total = 0
    for exon_type in exons_type:
        rare_syn = 0
        common_syn = 0
        rare_missense = 0
        common_missense = 0
        rare_high = 0
        common_high = 0
        frequency = None
        total_common = 0
        total_rare = 0
        file_name = path + "/results/" + exon_type + "_wvep_filtered.gff"
        #file_name = sys.argv[3] MODIFY LATER  
        out_name = path + "/results/" + exon_type + "_wvep_results.gff"
        num_lines = sum(1 for line in open(file_name))
        with open(file_name, 'r') as fhandle, open(out_name, 'w') as outhandle: 
            outhandle.write("gene" + "\t" + "gene_name" + "\t" + "exon" + "\t" + "strand" + "\t" 
                            + "rare_synonymous" + "\t" + "common_synonymous" + "\t" + "rare_missense" + "\t" 
                            + "common_missense" + "\t" + "rare_high_impact" + "\t" + "common_high_impact" + "\t" + "total_rare" + "\t" 
                            + "total_common" + "\n")
            counter = 0
            for line in fhandle:
                counter += 1
                #print(counter)
                if line.startswith("\n") or line.startswith(">	"):
                    continue
                line = line.replace("\n","")
                if line.startswith(">"):
                    outhandle.write(str(rare_syn) + "\t" + str(common_syn) + "\t" + str(rare_missense) + "\t" +
                                    str(common_missense) + "\t" + str(rare_high) + "\t" + str(common_high) + "\t" + str(total_rare) +
                                    "\t" + str(total_common) + "\n")
                    outhandle.write(line + "\t")
                    rare_syn = 0
                    common_syn = 0
                    rare_missense = 0
                    common_missense = 0
                    rare_high = 0
                    common_high = 0
                    total_rare = 0
                    total_common = 0
                    continue
                search_frequency = re.search("AF=([0-9]{1}\.[0-9]*)", line)
                search_impact = re.search("IMPACT=HIGH;", line)
                if search_frequency:
                    if float(search_frequency.group(1)) < 0.05:
                        frequency = "rare"
                    else:
                        frequency = "common"
                if line.split("\t")[6] == "synonymous_variant" and frequency == "rare":
                    rare_syn += 1
                    total_rare += 1
                elif line.split("\t")[6] == "synonymous_variant" and frequency == "common":
                    common_syn += 1
                    total_common += 1
                elif line.split("\t")[6] != "synonymous_variant" \
                    and line.split("\t")[6] != "inframe_deletion" \
                    and line.split("\t")[6] != "frameshift_variant" \
                    and frequency == "rare":
                    rare_missense += 1
                    total_rare += 1
                elif line.split("\t")[6] != "synonymous_variant" \
                    and line.split("\t")[6] != "inframe_deletion" \
                    and line.split("\t")[6] != "frameshift_variant" \
                    and frequency == "common":
                    common_missense += 1
                    total_common += 1
                if search_impact and frequency == "rare" and line.split("\t")[6] != "inframe_deletion" and line.split("\t")[6] != "frameshift_variant":
                    rare_high += 1
                    total_rare += 1
                elif search_impact and frequency == "common" and line.split("\t")[6] != "inframe_deletion" and line.split("\t")[6] != "frameshift_variant":
                    common_high += 1
                    total_common += 1
                #total += 1
                if counter == num_lines:
                    outhandle.write(str(rare_syn) + "\t" + str(common_syn) + "\t" + str(rare_missense) + "\t" 
                    + str(common_missense) + "\t" + str(rare_high) + "\t" + str(common_high) + "\t" 
                    + str(total_rare) + "\t" + str(total_common) + "\n")