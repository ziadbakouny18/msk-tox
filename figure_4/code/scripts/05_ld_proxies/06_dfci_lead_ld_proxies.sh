#!/usr/bin/env bash
# 06_dfci_lead_ld_proxies.sh
#
# Build LD proxy sets for DFCI's genome-wide significant SNPs, mirroring the
# method Zeyun applied to our leads (results.docx, "Replication in the DFCI
# cohort"):
#
#   proxies = SNPs within +/-500 kb and r2 >= 0.8 of the lead, from
#   PLINK2 --r2-unphased --ld-snp-list --ld-window-kb 500 --ld-window-r2 0.8,
#   computed in two reference panels, merged, de-duplicated by genomic
#   position and annotated by source. GRCh37/hg19 throughout, matched by
#   chromosome and position.
#
# WHAT WE CAN AND CANNOT MATCH
#   Zeyun's panel (i) is in-sample PROFILE LD. We have no PROFILE genotypes, and
#   our own cohort is present here only as pancancer_gwas_cohort_hardcalls.bim
#   with no .bed/.fam, so NO in-sample panel is reproducible on this machine.
#   Panel (ii) is reproducible in full: 1000 Genomes Phase 3, EUR and EUR+AFR.
#   Our proxy sets are therefore a subset of the two-panel sets Zeyun used --
#   they can confirm a replication but cannot rule one out.
#
# Reuses the self-contained conda env and the streaming approach from
# ../../../generate_LD_1000g/run_ld_proxies.sh; 1000G windows are streamed over
# HTTPS from the public EBI mirror, so nothing local leaves this machine.
#
# Usage (from the `analysis` folder):
#   bash scripts/06_dfci_lead_ld_proxies.sh
set -euo pipefail

ANALYSIS="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
UPDATE="$(dirname "$ANALYSIS")"
RESULTS="$UPDATE/results"
ZIAD="$(dirname "$(dirname "$UPDATE")")"
LDTOOLS="$ZIAD/generate_LD_1000g"

export PATH="$LDTOOLS/env/bin:$PATH"
export HTS_USE_CURL=1

OUT="$ANALYSIS/output/ld_proxies"
WORK="$OUT/work"
mkdir -p "$OUT" "$WORK"

WINDOW_KB=500
R2MIN=0.8
BASE="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"
PANEL="$LDTOOLS/ref/integrated_call_samples_v3.20130502.ALL.panel"

command -v plink2   >/dev/null || { echo "plink2 not on PATH"; exit 1; }
command -v bcftools >/dev/null || { echo "bcftools not on PATH"; exit 1; }
[ -s "$PANEL" ] || { echo "1000G sample panel missing: $PANEL"; exit 1; }

# ---- 1) DFCI lead SNPs, de-duplicated -------------------------------------
# dfci_sig_snps.tsv has 15 rows but only 12 distinct variants (some SNPs are
# significant for more than one phenotype). LD depends only on the variant.
LEADS="$OUT/dfci_lead_snps.tsv"
awk -F'\t' 'NR==1 { for (i=1;i<=NF;i++) h[$i]=i; next }
            !seen[$h["RSID"]]++ {
              print $h["RSID"] "\t" $h["CHR"] "\t" $h["BP"] "\t" $h["REF"] "\t" $h["ALT"]
            }' "$RESULTS/dfci_sig_snps.tsv" \
  | sort -k2,2n -k3,3n > "$LEADS.body"
printf 'snp\tchr\tpos\ta1\ta2\n' > "$LEADS"
cat "$LEADS.body" >> "$LEADS"
rm -f "$LEADS.body"
echo "[*] $(( $(wc -l < "$LEADS") - 1 )) distinct DFCI lead SNPs"

# ---- 2) Fetch each 1000G window ONCE, in parallel --------------------------
# Streaming from EBI is the slow step by a wide margin, so it is separated from
# the LD computation: the window does not depend on the reference panel, and
# fetching serially (and twice, once per panel) is what made the first version
# of this script time out. Windows are cached, so re-running resumes.
FETCH_JOBS="${FETCH_JOBS:-4}"

# Only chr, pos and snp are passed as arguments. WORK, WINDOW_KB and BASE come
# through the environment because this path contains spaces, and xargs would
# otherwise split it into separate arguments.
fetch_window() {
  chr="$1"; pos="$2"; snp="$3"
  work="$FW_WORK"; window_kb="$FW_WINDOW_KB"; base="$FW_BASE"
  out="$work/${snp}_chr${chr}.vcf.gz"
  if [ -s "$out" ] && [ -s "$out.csi" ]; then
    echo "    cached  $snp chr$chr"
    return 0
  fi
  start=$(( pos - window_kb*1000 )); [ "$start" -lt 1 ] && start=1
  end=$(( pos + window_kb*1000 ))
  url="$base/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
  echo "    fetching $snp chr${chr}:${start}-${end}"
  if bcftools view -r "${chr}:${start}-${end}" "$url" -Oz -o "$out.part" 2>/dev/null \
     && bcftools index -f "$out.part" 2>/dev/null; then
    mv "$out.part" "$out"; mv "$out.part.csi" "$out.csi"
  else
    rm -f "$out.part" "$out.part.csi"
    echo "    !! fetch failed for $snp chr$chr"
  fi
}
export -f fetch_window
export FW_WORK="$WORK" FW_WINDOW_KB="$WINDOW_KB" FW_BASE="$BASE"

echo "[*] fetching 1000G windows (${FETCH_JOBS} in parallel) ..."
tail -n +2 "$LEADS" \
  | awk -F'\t' '{ print $2, $3, $1 }' \
  | xargs -P "$FETCH_JOBS" -n 3 bash -c 'fetch_window "$@"' _

n_win=$(find "$WORK" -name '*.vcf.gz' ! -name '*.part' | wc -l | tr -d ' ')
echo "[*] $n_win windows available"

# ---- 3) One proxy run per reference panel (local, fast) --------------------
for POPSPEC in EUR EUR+AFR; do
  POPTAG="${POPSPEC/+/_}"
  KEEP="$WORK/keep_${POPTAG}.txt"

  # Sample keep-list. POPSPEC may name more than one super-population, which is
  # how Zeyun's "EUR + AFR-ancestry reference samples" panel is built.
  awk -v pops="$POPSPEC" 'BEGIN { n=split(pops,a,"+"); for (i=1;i<=n;i++) want[a[i]]=1 }
                          NR>1 && ($3 in want) { print $1, $1 }' "$PANEL" > "$KEEP"
  echo "[*] panel $POPSPEC: $(wc -l < "$KEEP") samples"

  PROXY_TSV="$OUT/ld_proxies_1000G_${POPTAG}_r2_${R2MIN}.tsv"
  printf 'lead_snp\tlead_chr\tlead_pos\tlead_var_id\tpanel\tproxy_var_id\tproxy_chr\tproxy_pos\tr2\n' > "$PROXY_TSV"

  tail -n +2 "$LEADS" | while IFS=$'\t' read -r snp chr pos a1 a2; do
    [ -z "${snp:-}" ] && continue
    vcf="$WORK/${snp}_chr${chr}.vcf.gz"
    tag="${snp}_chr${chr}_${POPTAG}"
    if [ ! -s "$vcf" ]; then
      echo "    !! no window for $snp chr$chr - skipping"
      continue
    fi
    echo "==> $snp  chr${chr}:${pos}  [$POPSPEC]"

    # 3a) VCF -> PLINK bed, chosen panel, biallelic SNPs, ids chr:pos:ref:alt
    plink2 --vcf "$vcf" --double-id \
           --keep "$KEEP" \
           --snps-only --max-alleles 2 --min-alleles 2 \
           --set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len 200 \
           --rm-dup exclude-all \
           --make-bed --out "$WORK/${tag}" >/dev/null 2>&1

    # 3b) locate the lead by position (allele order in 1000G may differ)
    lead_id=$(awk -v p="$pos" '$4==p {print $2; exit}' "$WORK/${tag}.bim" || true)
    if [ -z "$lead_id" ]; then
      echo "    !! lead not found at chr${chr}:${pos} after filtering - skipping"
      continue
    fi

    # 3c) proxies, using the exact PLINK2 flags named in results.docx
    echo "$lead_id" > "$WORK/${tag}.leadlist"
    plink2 --bfile "$WORK/${tag}" \
           --r2-unphased --ld-snp-list "$WORK/${tag}.leadlist" \
           --ld-window-kb "$WINDOW_KB" --ld-window-r2 "$R2MIN" \
           --out "$WORK/${tag}_r2_${R2MIN}" >/dev/null 2>&1

    # 3d) append. The lead may sit in either column of the .vcor output.
    vcor="$WORK/${tag}_r2_${R2MIN}.vcor"
    if [ -s "$vcor" ]; then
      awk -v lead="$lead_id" -v snp="$snp" -v chr="$chr" -v pos="$pos" -v panel="$POPSPEC" '
        NR==1 { for (i=1;i<=NF;i++) h[$i]=i; next }
        {
          if ($h["ID_A"]==lead) { p=$h["ID_B"]; pc=$h["CHROM_B"]; pp=$h["POS_B"] }
          else                  { p=$h["ID_A"]; pc=$h["CHROM_A"]; pp=$h["POS_A"] }
          print snp"\t"chr"\t"pos"\t"lead"\t"panel"\t"p"\t"pc"\t"pp"\t"$h["UNPHASED_R2"]
        }' "$vcor" >> "$PROXY_TSV"
    fi
  done

  echo "[*] $POPSPEC: $(( $(wc -l < "$PROXY_TSV") - 1 )) proxy rows -> $PROXY_TSV"
done

echo "[*] DONE. Merge and lookup happens in 07_dfci_proxies_in_msk.R"
