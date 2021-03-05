#!/usr/bin/env python

"""
Convert SNP file in phylip format to a SNAPP input file. SNPs are assumed
to be unlinked. Script will phase all ambiguity codes and identify all 
bi-allelic SNPs. Then all SNPs with at least one individual represented 
per OTU will be retained. SNPs will be converted to [0,1,2] and 
written to nexus file for analysis. Major (0) and minor (2) alleles are 
coded to estimate u and v in BEAUti.

author: J. Satler
date: 5 Mar 2021

usage: python SNPtoSNAPP.py infile.snps traits.file seqs_per_otu random|length
"""
import sys
import random

#dictionary of ambiguity codes
ambig = {"Y":["C", "T"], "R":["A", "G"], "W":["A", "T"],
         "S":["C", "G"], "K":["G","T"], "M":["A", "C"]}

#Accepted nucleotides
dna = ['A', 'C', 'G', 'T']

def pop_association(Traits):
    """Sets the individuals to their respective populations.
       Also returns the sample counts per population."""
    with open(Traits, 'r') as traits:
        Pops = {}
        Pop_counts = {}
        for line in traits:
            line = line.strip().split()
            Pops[line[0]] = line[1]
            if line[1] in Pop_counts:
                Pop_counts[line[1]] += 1
            else:
                Pop_counts[line[1]] = 1
        return Pops, Pop_counts

def snp(data, pops):
    """read in data and return list of SNPs for samples in traits file"""
    with open(data, 'r') as SNPs:
        return [line.strip().split() for line in SNPs
                if line.rstrip() if line.split()[0] in pops]

def sub(sp, pops, pop_counts, threshold, crit):
    """subsample individuals in species if over threshold"""
    filtered = []
    popAbove = {i:[] for i in set(pops.values()) if pop_counts[i] > threshold}

    for i in sp:
        #take all samples from species below threshold
        if pop_counts[pops[i[0]]] <= threshold:
            filtered.append(i)
        #subsample randomly or by seq length
        else:
            l = 0
            for j in i[1]:
                if j.upper() in ambig or j.upper() in dna:
                    l += 1
            popAbove[pops[i[0]]].append([i[0], i[1], l])

    for k, v in popAbove.items():
        #select individuals at random
        if crit.upper().startswith("R"):
            ch = random.sample(v, threshold)
            for sam in ch:
                filtered.append([sam[0], sam[1]])

        #select individuals by sequence length
        elif crit.upper().startswith("L"):
            v.sort(key=lambda x: x[2], reverse=True)
            for sam in range(threshold):
                filtered.append([v[sam][0], v[sam][1]])
        else:
            print "Do you want to subsample individuals randomly or by length?"
            sys.exit()
    return filtered

def biallelic(data):
    """check if SNP is biallelic, and retain those SNPs"""
    m = zip(*[i[1] for i in data])
    taxa = [i[0] for i in data]

    #Biallelic SNPs (starting with sample labels)
    m_bi = []
    for i in m:
        alleles = []
        for a in i:
            if a.upper() in dna:
                alleles.append(a.upper())
            elif a.upper() in ambig:
                alleles.extend(ambig[a.upper()])
        #if biallelic, add to new matrix
        if len(set(alleles)) == 2:
            m_bi.append(i)
    return m_bi, taxa

def threshold(bi_SNP_matrix, taxa, OTUs):
    """remove any SNP not represented in all species"""
    bi_thr_SNP = []
    for i in bi_SNP_matrix:
        counts = {i:0 for i in set(OTUs.values())}

        for j in range(len(i)):
            if i[j] != 'N' and i[j] != '-':
                counts[OTUs[taxa[j]]] += 1

        #add SNP if all pops meet threshold - 1 ind right now
        if all(val >= 1 for val in counts.values()) == True:
            bi_thr_SNP.append(i)
    return bi_thr_SNP

def all_code(Unlink_Matrix):
    """code alleles as [0,1,2] for SNAPP"""
    recoded_mat = []
    for snp in Unlink_Matrix:
        a1, a2, het = get_two_alleles(snp)

        recode_locus = []
        for i in snp:
            if i == a1:
                recode_locus.append('0')
            elif i == a2:
                recode_locus.append('2')
            elif i == het:
                recode_locus.append('1')
            else:
                recode_locus.append('-')
        recoded_mat.append(recode_locus)
    return recoded_mat

def get_two_alleles(locus):
    """get two alleles for a locus"""
    a = [allele for allele in set(locus) if not allele.upper() == 'N'
         and not allele == '-']
    ab = []
    het = ''
    for i in a:
        if i in dna:
            ab.append(i)
        else:
            het += i
            ab.extend(ambig[i])

    #the two alleles
    all = list(set(ab))
    a1 = all[0]
    a2 = all[1]
    if len(het) > 1:
        print "ERROR. Too many alleles in your locus due to an excess of het sites."
        sys.exit()

    #Get major and minor alleles
    aT = [ambig[allele] if allele in ambig else
          allele for allele in locus
          if not allele.upper() == 'N' and not allele == '-']
    #flattened list of all alleles
    aTf = [val for sublist in aT for val in sublist]

    #calculate and return major, minor, het
    if (aTf.count(a1) / float(len(aTf))) > 0.5:
        #print "Major: {0}, {1}\tMinor: {2}, {3}".format(a1, aTf.count(a1),
        #                                                a2, aTf.count(a2))
        return a1, a2, het
    else:
        #print "Major: {0}, {1}\tMinor: {2}, {3}".format(a2, aTf.count(a2),
        #                                                a1, aTf.count(a1))
        return a2, a1, het

def concat_sites(coded_m):
    """concat sites for each individual"""
    return zip(*coded_m)

def out_for_SNAPP(final_matrix, out_type, outfile, Pops, taxa):
    """write output nexus file"""
    preamble = """#NEXUS\n\nbegin data;\ndimensions ntax={0} nchar={1}; \
                  \nformat datatype={2} gap=-;\
                  \nmatrix\n""".format(len(taxa),
                                         len(final_matrix[0]),
                                         'integerdata symbols="012"'
                                         if out_type is "_snapp_"
                                         else 'dna missing=?')

    #turn matrix into dictionary so it can be sorted by OTU
    m = create_dict_matrix(final_matrix, taxa, Pops)

    with open(outfile[:-4] + out_type + sys.argv[3] + 'ind_' +
              sys.argv[4].upper()[0] + '.nex', 'w') as nex:
        nex.write(preamble)
        for k in sorted(Pops, key=Pops.get):
            if k in taxa:
                out = "{0}_{1}{2}{3}\n".format(Pops[k], k.replace("-", "_"),
                                               ' ' * (40 - len(Pops[k] + k)), m[k])
                nex.write(out)
        nex.write(";\nend;")

def create_dict_matrix(final_matrix, taxa, Pops):
    """create dictionary for final sort"""
    filtered = {}
    for p, ind in enumerate(final_matrix):
        seq = ''.join(ind)
        filtered[taxa[p]] = seq
    return filtered

def main():
    if len(sys.argv[1:]) != 4:
        print "python SNPtoSNAPP.py infile.snps traits.file seqs_per_otu random|length"
        sys.exit()

    Pops, Pop_counts = pop_association(sys.argv[2])
    d = snp(sys.argv[1], Pops)
    filtered = sub(d, Pops, Pop_counts, int(sys.argv[3]), sys.argv[4])
    bi_snps, taxa_sub = biallelic(filtered)
    bi_thr = threshold(bi_snps, taxa_sub, Pops)

    #recode matrix to [0,1,2]
    coded_mat = all_code(bi_thr)
    loci_recode = concat_sites(coded_mat)

    out_for_SNAPP(loci_recode, "_snapp_", sys.argv[1], Pops, taxa_sub)
    out_for_SNAPP(concat_sites(bi_thr), "_dna_", sys.argv[1], Pops, taxa_sub)

if __name__ == '__main__':
    main()
