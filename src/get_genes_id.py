########################
fhandle = open('/home/lmartinezg/Documents/Laura/appris_clinvar/exons/gencode.v37.alt_mane_appris.csv')
exons = []
for line in fhandle:
    exon = line.split('\t')[0] + '-' + line.split('\t')[1] + '-' + line.split('\t')[2]
    if exon not in exons:
        exons.append(exon)
    else:
        print(line)
fhandle.close()

khandle = open('/home/lmartinezg/Documents/Laura/appris_clinvar/Exons_old/gencode.v37.alternative_exons.csv')
outhandle = open('/home/lmartinezg/Documents/Laura/appris_clinvar/exons/gencode.v37.alt_mane_appris_genename.csv', 'w')
for line in khandle:
    exon = line.split('\t')[3] + '-' + line.split('\t')[4] + '-' + line.split('\t')[5]
    if exon in exons:
        outhandle.write(line)
        exons.remove(exon)
khandle.close()
outhandle.close()
print(str(len(exons)) + ' alt_mane_appris')

################################

phandle = open('/home/lmartinezg/Documents/Laura/appris_clinvar/exons/gencode.v37.alt_mane_no_overlap.csv')
exons_mane_no = []
for line in phandle:
    exon = line.split('\t')[0] + '-' + line.split('\t')[1] + '-' + line.split('\t')[2]
    if exon not in exons_mane_no:
        exons_mane_no.append(exon)
    else:
        print(line)
phandle.close()

mhandle = open('/home/lmartinezg/Documents/Laura/appris_clinvar/Exons_old/gencode.v37.mane_alt_exons.csv')
outhandle = open('/home/lmartinezg/Documents/Laura/appris_clinvar/exons/gencode.v37.alt_mane_no_overlap_genename.csv', 'w')
for line in mhandle:
    exon = line.split('\t')[3] + '-' + line.split('\t')[4] + '-' + line.split('\t')[5]
    if exon in exons_mane_no:
        outhandle.write(line)
        exons_mane_no.remove(exon)
mhandle.close()
outhandle.close()
print(str(len(exons_mane_no)) + ' mane_alt_exons')

#################################

dhandle = open('/home/lmartinezg/Documents/Laura/appris_clinvar/exons/gencode.v37.alt_mane_unique.csv')
exons_mane_unique = []
for line in dhandle:
    exon = line.split('\t')[0] + '-' + line.split('\t')[1] + '-' + line.split('\t')[2]
    if exon not in exons_mane_unique:
        exons_mane_unique.append(exon)
    else:
        print(line)
dhandle.close()

qhandle = open('/home/lmartinezg/Documents/Laura/appris_clinvar/Exons_old/gencode.v37.mane_alt_exons.csv')
outhandle = open('/home/lmartinezg/Documents/Laura/appris_clinvar/exons/gencode.v37.alt_mane_unique_genename.csv', 'w')
for line in qhandle:
    exon = line.split('\t')[3] + '-' + line.split('\t')[4] + '-' + line.split('\t')[5]
    if exon in exons_mane_unique:
        outhandle.write(line)
        exons_mane_unique.remove(exon)
qhandle.close()

print(str(len(exons_mane_unique)) + ' alt_mane_unique')

dhandle = open('/home/lmartinezg/Documents/Laura/appris_clinvar/Exons_old/gencode.v37.principal_exons.csv')
for line in dhandle:
    exon = line.split('\t')[3] + '-' + line.split('\t')[4] + '-' + line.split('\t')[5]
    if exon in exons_mane_unique:
        outhandle.write(line)
        exons_mane_unique.remove(exon)
print(str(len(exons_mane_unique)) + ' alt_mane_unique')
dhandle.close()
outhandle.close()

##############################

dhandle = open('/home/lmartinezg/Documents/Laura/appris_clinvar/exons/gencode.v37.alt_no_overlap.csv')
exons_alt_ov = []
for line in dhandle:
    exon = line.split('\t')[0] + '-' + line.split('\t')[1] + '-' + line.split('\t')[2]
    if exon not in exons_alt_ov:
        exons_alt_ov.append(exon)
    else:
        print(line)
dhandle.close()

qhandle = open('/home/lmartinezg/Documents/Laura/appris_clinvar/Exons_old/gencode.v37.alternative_exons.csv')
outhandle = open('/home/lmartinezg/Documents/Laura/appris_clinvar/exons/gencode.v37.alt_mane_unique_genename.csv', 'w')
for line in qhandle:
    exon = line.split('\t')[3] + '-' + line.split('\t')[4] + '-' + line.split('\t')[5]
    if exon in exons_alt_ov:
        outhandle.write(line)
        exons_mane_no.remove(exon)
qhandle.close()

print(str(len(exons_alt_ov)) + ' alt_mane_unique')