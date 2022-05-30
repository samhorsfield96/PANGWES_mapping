import mappy as mp
import re
import pandas

def map_pairs(infile, outfile, *annotation_files):
    annotation_dict = {}
    for file in annotation_files:
        df = pandas.read_excel(file)
        for index, row in df.iterrows():
            annotation_dict[row[0]] = (row[2], row[3])

    a = mp.Aligner("SPARC_CDS_dna_sequences.mmi", preset="sr")

    unitig_pairs = []
    with open(infile, "r") as f:
        for line in f:
            split_line = line.strip()
            split_line = split_line.split(" ")
            pair = (split_line[-2], split_line[-1])
            unitig_pairs.append(pair)

    mapping_results = [None] * len(unitig_pairs)

    for index, pair in enumerate(unitig_pairs):
        mappings1 = []
        unitig1, unitig2 = pair
        for hit in a.map(unitig1):
            hit_ref = hit.ctg
            hit_ref = hit_ref.split("_")
            sample_accession = hit_ref[0] + "_" + hit_ref[1]
            COG_accession = hit_ref[3]
            cigar = hit.cigar_str
            matches = re.findall(r'(\d+)([A-Z]{1})', cigar)
            num_matches = 0
            for m in matches:
                if m[1] == 'M':
                    num_matches += int(m[0])
            prop_matches = num_matches / len(unitig1)
            mappings1.append((sample_accession, COG_accession, hit.r_st, hit.r_en, hit.cigar_str, prop_matches))
        mappings2 = []
        for hit in a.map(unitig2):
            hit_ref = hit.ctg
            hit_ref = hit_ref.split("_")
            sample_accession = hit_ref[0] + "_" + hit_ref[1]
            COG_accession = hit_ref[3]
            cigar = hit.cigar_str
            matches = re.findall(r'(\d+)([A-Z]{1})', cigar)
            num_matches = 0
            for m in matches:
                if m[1] == 'M':
                    num_matches += int(m[0])
            prop_matches = num_matches / len(unitig2)
            mappings2.append((sample_accession, COG_accession, hit.r_st, hit.r_en, hit.cigar_str, prop_matches))
        mapping_results[index] = (mappings1, mappings2)

    with open(outfile, "w") as o:
        o.write("Unitig_Pair_ID\tUnitig_seq1\tCOG_accession1\tGene_name1\tAnnotation1\tHit_start1\tHit_end1\tPerc_id1\tSource1\tUnitig_seq2\tCOG_accession2\tGene_name2\tAnnotation2\tHit_start2\tHit_end2\tPerc_id2\tSource2\n")
        for index, pair in enumerate(unitig_pairs):
            mappings1, mappings2 = mapping_results[index]
            mappings1 = sorted(mappings1, key=lambda i: i[-1], reverse=True)
            mappings2 = sorted(mappings2, key=lambda i: i[-1], reverse=True)

            # take top hit from each mappings
            map = [str(index), pair[0], "NA", "NA", "NA", "NA", "NA", "NA", "BLAST", pair[1], "NA", "NA", "NA", "NA", "NA", "NA", "BLAST"]
            if mappings1:
                mappings = mappings1[0]
                gene_name = "NA"
                annotation = "NA"
                COG_accession = mappings[1]
                if COG_accession in annotation_dict:
                    gene_name, annotation = annotation_dict[COG_accession]

                map[2] = COG_accession
                map[3] = gene_name
                map[4] = annotation
                map[5] = mappings[2]
                map[6] = mappings[3]
                map[7] = mappings[5]
                map[8] = "SPARC"
            if mappings2:
                mappings = mappings2[0]
                gene_name = "NA"
                annotation = "NA"
                COG_accession = mappings[1]
                if COG_accession in annotation_dict:
                    gene_name, annotation = annotation_dict[COG_accession]

                map[10] = COG_accession
                map[11] = gene_name
                map[12] = annotation
                map[13] = mappings[2]
                map[14] = mappings[3]
                map[15] = mappings[5]
                map[16] = "SPARC"

            for index, entry in enumerate(map):
                o.write(str(entry))
                if index < len(map) - 1:
                    o.write("\t")
            o.write("\n")

def parse_annotations(infile, outfile):
    in_df = pandas.read_excel(infile)
    d = {}
    columnsNamesList = list(in_df.columns.values)
    columnsNamesList = columnsNamesList[2:]
    for index, row in in_df.iterrows():
        CDS_name = row[0]
        gene_name = row[1]
        annotation_list = row[2:]
        non_zero = [index for index, entry in enumerate(annotation_list) if entry]
        annotation = ""
        for index in non_zero:
            annotation += columnsNamesList[index]
        d[CDS_name] = [gene_name, annotation]



    df = pandas.DataFrame(data=d)
    df = df.transpose()
    df.columns = ["Gene_name", "Annotation"]

    df = df.replace(to_replace="-",
               value="")
    df.to_excel(outfile)


if __name__ == "__main__":
    #parse_annotations("functional_annotation/pnas.1613937114.sd01.xlsx", "functional_annotation/pnas.1613937114.sd01.parsed.xlsx")
    map_pairs("unitigs/maela_k151.txt", "maela_k151_mapped.txt", "functional_annotation/NIHMS74007-supplement-Supplementary_Dataset_1.xls", "functional_annotation/NIHMS74007-supplement-Supplementary_Dataset_2.xls")
