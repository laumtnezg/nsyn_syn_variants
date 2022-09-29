import re
exons_list = ['gencode.v37.prin_qsplice_unique', 'gencode.v37.alt_qsplice_unique']
vep_annot_file = path + '/data/all_results_filtered.txt'
for exon_type in exons_list:
    print(exon_type)
    chromosomes = {}
    fname = path + 'data/' + exon_type  + '.csv'
    with open(fname, 'r') as fhandle:
        for line in fhandle:
            region = line.split('\t')[0] + '-' + line.split('\t')[1] + '-' + line.split('\t')[2]# + '-' + line.split('\t')[3]
            positions = line.split('\t')[1] + '-' + line.split('\t')[2]
            if line.split('\t')[0].split('chr')[1] not in chromosomes:
                #print(line.split('\t')[0].split('chr')[1])
                chromosomes[line.split('\t')[0].split('chr')[1]] = [positions]
            else:
                chromosomes[line.split('\t')[0].split('chr')[1]].append(positions)
            #regions.append(region)
            outname = path + 'results/' + exon_type + '_wvep.txt'
            #out_no_name = path + 'results/' + exon_type + '_nowvep.txt'
        print(len(chromosomes['1']))
        #exit()
        with open(outname, 'w') as outhandle, open(vep_annot_file, 'r') as vep_file:#, open(out_no_name,'w') as out_no:
            for line in vep_file:
                line = line.replace('\n', '')
                if line.startswith('#'): continue
                #print(line)
                pos = line.split("\t")[1].split(":")[1]
                chrom = line.split("\t")[1].split(":")[0]
                for exon in chromosomes[chrom]:
                    start = exon.split('-')[0]
                    end = exon.split('-')[1]
                    if (start <= pos) and (pos <= end):
                        outhandle.write('>' + chrom + '-' + exon + '\n' + line + '\n')
        new_out = path + 'results/' + exon_type + '_wvep_merged.txt'
        with open(new_out,'w') as outhandle, open(outname, 'r') as fhandle:
            exons = []
            for line in fhandle:
                line = line.replace('\n', '')
                if line.startswith(">"):
                    exon = line.split('>')[1]
                    if exon in exons:
                        continue
                    else:
                        outhandle.write(line + "\n")
                        exons.append(exon)
                else:
                    outhandle.write(line + '\n')

            