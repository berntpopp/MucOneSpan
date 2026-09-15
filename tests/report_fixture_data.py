"""Small real reference/VCF fixtures shared by session and browser checks."""


def fixture_files(tmp_path):
    fasta = tmp_path / "reference.fa"
    fasta.write_text("".join(f">contig_{n}\n" + "A" * 2200 + "\n" for n in [1, 11, 21]))
    paths = {}
    for allele, contig in [("allele_1", "contig_11"), ("allele_2", "contig_21")]:
        vcf = tmp_path / f"{allele}.vcf"
        vcf.write_text(
            f"##fileformat=VCFv4.2\n##contig=<ID={contig},length=2200>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n{contig}\t700\t.\tA\tC\t60\tPASS\t.\n"
        )
        paths[allele] = vcf
    summary = {
        "alleles": {
            key: {"contig_name": f"contig_{n}", "reads": 100, "length": n + 9}
            for key, n in [("allele_1", 11), ("allele_2", 21)]
        }
    }
    return summary, fasta, paths
