import mappy as mp
import re

def map_pairs(infile, outfile):
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
            COG_accession = hit_ref[2] + "_" + hit_ref[3]
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
            COG_accession = hit_ref[2] + "_" + hit_ref[3]
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
        o.write("Unitig_Pair_ID\tUnitig_No\tUnitig_Seq\tReference_Sample_Accession\tReference_COG_ID\tAlignment_start\tAlignment_end\tPerc_Identity\n")
        for index, pair in enumerate(unitig_pairs):
            mappings1, mappings2 = mapping_results[index]
            for mapping in mappings1:
                o.write(str(index) + "\t" + "1" + "\t" + pair[0] + "\t" + str(mapping[0]) + "\t" + str(mapping[1]) + "\t" + str(mapping[2]) + "\t" + str(mapping[3]) + "\t" + str(mapping[5]) + "\n")
            if not mappings1:
                o.write(str(index) + "\t" + "1" + "\t" + pair[0] + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\n")
            for mapping in mappings2:
                o.write(str(index) + "\t" + "2" + "\t" + pair[1] + "\t" + str(mapping[0]) + "\t" + str(mapping[1]) + "\t" + str(mapping[2]) + "\t" + str(mapping[3]) + "\t" + str(mapping[5]) + "\n")
            if not mappings2:
                o.write(str(index) + "\t" + "2" + "\t" + pair[1] + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\n")

if __name__ == "__main__":
    map_pairs("unitigs/maela_k101.txt", "maela_k101_mapped.txt")
