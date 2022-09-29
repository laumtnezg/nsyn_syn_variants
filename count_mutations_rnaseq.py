import re
exons_list = ['gencode.v37.prin_qsplice_unique', 'gencode.v37.alt_qsplice_unique']
for exon_type in exons_list:
    print(exon_type)
    fname = path + 'results/' + exon_type + '_wvep_merged.txt'
    outname = path + 'results/' + exon_type + '_wvep_filtered.txt'
    with open(fname, 'r') as fhandle, open(outname, 'w') as outhandle:
        variations = []
        for line in fhandle:
            line = line.replace('\n', '')
            if line.startswith('>'):
                outhandle.write(line + '\n')
            else:
                if line.split("\t")[0] not in variations:
                    outhandle.write((line + '\n'))
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
    rare_syn = 0
    frequency = None
    total_common = 0
    total_rare = 0
    out_new = path + 'results/' + exon_type + '_wvep_results.txt'
    num_lines = sum(1 for line in open(outname))
    with open(outname, 'r') as fhandle, open(out_new, 'w') as outhandle:
        outhandle.write("exon" + "\t" + "rare_synonymous" + "\t" + "common_synonymous" + "\t"
                         + "rare_missense" + "\t" + "common_missense" + "\t" + "rare_high_impact" + "\t"
                          + "common_high_impact" + "\t" + "total_rare" + "\t" + "total_common" + "\n")
        counter = 0
        for line in fhandle:
            counter += 1
            if line.startswith("\n") or line.startswith(">	"):
                continue
            line = line.replace("\n","")
            if line.startswith(">"):
                outhandle.write(str(rare_syn) + "\t" + str(common_syn) + "\t" + str(rare_missense) + "\t" +
                                str(common_missense) + "\t" + str(rare_high) + "\t" + str(common_high) + "\t" + str(total_rare) +
                                "\t" + str(total_common) + "\n")
                outhandle.write(line.split('>')[1] + "\t")
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
            if search_impact and frequency == "rare" \
                and line.split("\t")[6] != "inframe_deletion" \
                and line.split("\t")[6] != "frameshift_variant":
                rare_high += 1
                total_rare += 1
            elif search_impact and frequency == "common" \
                and line.split("\t")[6] != "inframe_deletion" \
                and line.split("\t")[6] != "frameshift_variant":
                common_high += 1
                total_common += 1
            if counter == num_lines:
                outhandle.write(str(rare_syn) + "\t" + str(common_syn) + "\t" + str(rare_missense) + "\t" 
                + str(common_missense) + "\t" + str(rare_high) + "\t" + str(common_high) + "\t" 
                + str(total_rare) + "\t" + str(total_common) + "\n")

