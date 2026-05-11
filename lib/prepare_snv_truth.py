#!/usr/bin/env python3
import argparse
import sys
from collections import defaultdict


BASES = {"A", "C", "G", "T"}


def warn(message):
    print(f"[WARN] {message}", file=sys.stderr)


def fail(message):
    print(f"[ERROR] {message}", file=sys.stderr)
    sys.exit(1)


def read_regions(path):
    regions = defaultdict(list)
    active_chroms = []
    seen = set()
    with open(path, "r") as handle:
        for lineno, line in enumerate(handle, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 3:
                fail(f"{path}:{lineno}: regions BED needs at least 3 columns")
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            if end <= start:
                fail(f"{path}:{lineno}: region end must be greater than start")
            regions[chrom].append((start, end))
            if chrom not in seen:
                active_chroms.append(chrom)
                seen.add(chrom)
    if not regions:
        fail(f"{path}: no regions found")
    return regions, active_chroms


def overlaps_regions(chrom, start, end, regions):
    for region_start, region_end in regions.get(chrom, []):
        if start < region_end and end > region_start:
            return True
    return False


def read_snp_bed(path, label, regions):
    snps = {}
    skipped_outside = 0
    with open(path, "r") as handle:
        for lineno, line in enumerate(handle, 1):
            raw = line.strip()
            if not raw or raw.startswith("#"):
                continue
            fields = raw.split()
            if len(fields) < 5:
                fail(
                    f"{path}:{lineno}: SNP BED for {label} must have at least "
                    "5 columns: chrom start end SNP allele"
                )
            chrom = fields[0]
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError:
                fail(f"{path}:{lineno}: SNP coordinates must be integers")
            if end - start != 1:
                fail(f"{path}:{lineno}: only 1 bp SNP intervals are supported")
            if fields[3].upper() != "SNP":
                fail(f"{path}:{lineno}: expected variant type SNP, found {fields[3]}")
            allele = fields[4].upper()
            if allele not in BASES:
                fail(f"{path}:{lineno}: SNP allele must be A/C/G/T, found {fields[4]}")
            if not overlaps_regions(chrom, start, end, regions):
                skipped_outside += 1
                continue
            key = (chrom, start, end)
            if key in snps and snps[key] != allele:
                fail(
                    f"{path}:{lineno}: duplicate SNP {chrom}:{start}-{end} "
                    f"has conflicting alleles {snps[key]} and {allele}"
                )
            snps[key] = allele
    if skipped_outside:
        warn(f"{path}: skipped {skipped_outside} SNPs outside simulated regions")
    return snps


class FastaReader:
    def __init__(self, fasta_path, positions):
        self.fasta_path = fasta_path
        self.positions = positions
        self.fai = self._read_fai(f"{fasta_path}.fai")
        self.cache = {}
        self.handle = None
        if not self.fai:
            warn(f"{fasta_path}.fai not found; scanning FASTA to fetch SNV reference bases")
            self.cache = self._scan_fasta()
        else:
            self.handle = open(fasta_path, "rb")

    @staticmethod
    def _read_fai(path):
        fai = {}
        try:
            with open(path, "r") as handle:
                for line in handle:
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) < 5:
                        continue
                    fai[fields[0]] = {
                        "length": int(fields[1]),
                        "offset": int(fields[2]),
                        "line_bases": int(fields[3]),
                        "line_width": int(fields[4]),
                    }
        except FileNotFoundError:
            return {}
        return fai

    def _scan_fasta(self):
        wanted = defaultdict(set)
        for chrom, pos in self.positions:
            wanted[chrom].add(pos)
        found = defaultdict(dict)
        chrom = None
        offset = 0
        with open(self.fasta_path, "r") as handle:
            for line in handle:
                line = line.rstrip("\n")
                if line.startswith(">"):
                    chrom = line[1:].split()[0]
                    offset = 0
                    continue
                if chrom not in wanted:
                    offset += len(line)
                    continue
                for pos in list(wanted[chrom]):
                    idx = pos - offset - 1
                    if 0 <= idx < len(line):
                        found[chrom][pos] = line[idx].upper()
                        wanted[chrom].remove(pos)
                offset += len(line)
        missing = [(chrom, pos) for chrom, vals in wanted.items() for pos in vals]
        if missing:
            first = missing[0]
            fail(f"Reference FASTA does not contain requested position {first[0]}:{first[1]}")
        return found

    def base(self, chrom, pos):
        if self.fai:
            if chrom not in self.fai:
                fail(f"Reference FASTA index has no contig {chrom}")
            meta = self.fai[chrom]
            if pos < 1 or pos > meta["length"]:
                fail(f"Reference position out of bounds: {chrom}:{pos}")
            zero = pos - 1
            byte_offset = (
                meta["offset"]
                + (zero // meta["line_bases"]) * meta["line_width"]
                + (zero % meta["line_bases"])
            )
            self.handle.seek(byte_offset)
            return self.handle.read(1).decode("ascii").upper()
        return self.cache[chrom][pos]

    def close(self):
        if self.handle is not None:
            self.handle.close()


def write_bed(path, records):
    with open(path, "w") as handle:
        for chrom, start, end, allele in records:
            handle.write(f"{chrom}\t{start}\t{end}\tSNP\t{allele}\t0\n")


def vcf_header():
    return [
        "##fileformat=VCFv4.2\n",
        "##source=strand-creator\n",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tstrand_creator\n",
    ]


def write_vcfs(out_vcf, per_chrom_dir, records):
    header = vcf_header()
    with open(out_vcf, "w") as handle:
        handle.writelines(header)
        for rec in records:
            handle.write(rec["line"])

    by_chrom = defaultdict(list)
    for rec in records:
        by_chrom[rec["chrom"]].append(rec["line"])

    if per_chrom_dir:
        import os

        os.makedirs(per_chrom_dir, exist_ok=True)
        for chrom, lines in by_chrom.items():
            with open(os.path.join(per_chrom_dir, f"{chrom}.vcf"), "w") as handle:
                handle.writelines(header)
                handle.writelines(lines)


def main():
    parser = argparse.ArgumentParser(
        description="Prepare real heterozygous SNV truth from haplotype VISOR SNP BEDs."
    )
    parser.add_argument("--ref", required=True, help="Reference FASTA used by VISOR")
    parser.add_argument("--regions", required=True, help="Simulation regions BED")
    parser.add_argument("--h1", required=True, help="Haplotype 1 SNP BED")
    parser.add_argument("--h2", required=True, help="Haplotype 2 SNP BED")
    parser.add_argument("--out-h1", required=True, help="Filtered h1 VISOR SNP BED")
    parser.add_argument("--out-h2", required=True, help="Filtered h2 VISOR SNP BED")
    parser.add_argument("--out-vcf", required=True, help="Combined heterozygous SNV VCF")
    parser.add_argument("--per-chrom-dir", required=True, help="Directory for per-chrom VCFs")
    parser.add_argument("--positions", required=True, help="StrandPhaseR positions TSV")
    parser.add_argument("--summary", required=True, help="SNV summary TSV")
    parser.add_argument("--min-snps-per-chrom", type=int, default=5)
    args = parser.parse_args()

    regions, active_chroms = read_regions(args.regions)
    h1 = read_snp_bed(args.h1, "h1", regions)
    h2 = read_snp_bed(args.h2, "h2", regions)

    positions = [(chrom, end) for chrom, _start, end in sorted(set(h1) | set(h2))]
    if not positions:
        fail("No SNPs remain after filtering to simulated regions")

    fasta = FastaReader(args.ref, positions)
    h1_records = []
    h2_records = []
    vcf_records = []
    skipped = 0

    for idx, key in enumerate(sorted(set(h1) | set(h2)), 1):
        chrom, start, end = key
        pos = end
        ref = fasta.base(chrom, pos)
        if ref not in BASES:
            warn(f"Skipping {chrom}:{pos}; reference base is {ref}")
            skipped += 1
            continue

        alt1 = h1.get(key)
        alt2 = h2.get(key)
        hap1 = alt1 or ref
        hap2 = alt2 or ref

        if alt1 == ref:
            hap1 = ref
            alt1 = None
        if alt2 == ref:
            hap2 = ref
            alt2 = None

        if hap1 == hap2:
            skipped += 1
            continue

        if hap1 != ref and hap2 != ref:
            warn(f"Skipping {chrom}:{pos}; non-reference alleles on both haplotypes are unsupported")
            skipped += 1
            continue

        alt = hap1 if hap1 != ref else hap2
        gt = "1|0" if hap1 != ref else "0|1"
        snv_id = f"SC_SNV_{len(vcf_records) + 1}"
        vcf_records.append(
            {
                "chrom": chrom,
                "line": f"{chrom}\t{pos}\t{snv_id}\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{gt}\n",
            }
        )
        if hap1 != ref:
            h1_records.append((chrom, start, end, alt))
        if hap2 != ref:
            h2_records.append((chrom, start, end, alt))

    fasta.close()

    if not vcf_records:
        fail("No heterozygous SNVs remain after reference validation")

    write_bed(args.out_h1, h1_records)
    write_bed(args.out_h2, h2_records)
    write_vcfs(args.out_vcf, args.per_chrom_dir, vcf_records)

    counts = defaultdict(int)
    with open(args.positions, "w") as handle:
        for rec in vcf_records:
            fields = rec["line"].split("\t")
            counts[fields[0]] += 1
            handle.write(f"{fields[0]}\t{fields[1]}\n")

    with open(args.summary, "w") as handle:
        handle.write("chrom\theterozygous_snvs\tmin_snps_for_phasing\tstatus\n")
        for chrom in active_chroms:
            count = counts[chrom]
            status = "PASS" if count >= args.min_snps_per_chrom else "LOW_SNV_COUNT"
            handle.write(f"{chrom}\t{count}\t{args.min_snps_per_chrom}\t{status}\n")
            if status != "PASS":
                warn(
                    f"{chrom} has {count} heterozygous SNVs after filtering; "
                    f"MosaiCatcher default min_snps_for_phasing is {args.min_snps_per_chrom}"
                )

    if skipped:
        warn(f"Skipped {skipped} non-informative or unsupported SNP records")


if __name__ == "__main__":
    main()
