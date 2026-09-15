#!/usr/bin/env bash
# Self-contained 1000G Phase 3 (hg19) LD-proxy finder.
# For each lead SNP in lead_snps.tsv: stream +/-500kb from 1000G over HTTPS,
# subset to a chosen super-population, compute proxies with r2 >= R2MIN in PLINK1.9.
# Everything stays inside this directory. No external project files are touched.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PATH="$HERE/env/bin:$PATH"
export HTS_USE_CURL=1   # let htslib/bcftools stream remote VCFs via libcurl

# ---- parameters -------------------------------------------------------------
POP="${POP:-EUR}"          # super-population: EUR/AFR/EAS/SAS/AMR
WINDOW_KB="${WINDOW_KB:-500}"
R2MIN="${R2MIN:-0.6}"
BASE="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"
PANEL_URL="$BASE/integrated_call_samples_v3.20130502.ALL.panel"
# -----------------------------------------------------------------------------

mkdir -p "$HERE/ref" "$HERE/work" "$HERE/out"
cd "$HERE"

# 1) Population panel + keep file ---------------------------------------------
if [ ! -s ref/integrated_call_samples_v3.20130502.ALL.panel ]; then
  echo "[*] downloading 1000G panel ..."
  curl -fsSL "$PANEL_URL" -o ref/integrated_call_samples_v3.20130502.ALL.panel
fi
awk -v pop="$POP" 'NR>1 && $3==pop {print $1, $1}' \
  ref/integrated_call_samples_v3.20130502.ALL.panel > "ref/keep_${POP}.txt"
echo "[*] $POP samples: $(wc -l < ref/keep_${POP}.txt)"

# 2) Per-SNP loop -------------------------------------------------------------
PROXY_TSV="out/ld_proxies_1000G_${POP}_r2_${R2MIN}.tsv"
echo -e "lead_snp\tlead_chr\tlead_pos\tlead_var_id\tpop\tproxy_var_id\tproxy_chr\tproxy_pos\tr2" > "$PROXY_TSV"

tail -n +2 lead_snps.tsv | while IFS=$'\t' read -r snp chr pos a1 a2; do
  [ -z "${snp:-}" ] && continue
  start=$(( pos - WINDOW_KB*1000 )); [ "$start" -lt 1 ] && start=1
  end=$(( pos + WINDOW_KB*1000 ))
  url="$BASE/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
  tag="${snp}_chr${chr}_${POP}"
  echo "==> $snp  chr${chr}:${start}-${end}"

  # 2a) stream just this window from the remote VCF
  bcftools view -r "${chr}:${start}-${end}" "$url" -Oz -o "work/${tag}.vcf.gz" 2>/dev/null
  bcftools index -f "work/${tag}.vcf.gz"

  # 2b) VCF -> PLINK bed, EUR only, biallelic SNPs, ids = chr:pos:ref:alt
  plink2 --vcf "work/${tag}.vcf.gz" --double-id \
         --keep "ref/keep_${POP}.txt" \
         --snps-only --max-alleles 2 --min-alleles 2 \
         --set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len 200 \
         --rm-dup exclude-all \
         --make-bed --out "work/${tag}" >/dev/null 2>&1

  # 2c) find the lead variant id by position (ref/alt order independent)
  lead_id=$(awk -v p="$pos" '$4==p {print $2; exit}' "work/${tag}.bim" || true)
  if [ -z "$lead_id" ]; then
    echo "    !! lead variant not found at chr${chr}:${pos} after filtering - skipping"
    continue
  fi
  echo "    lead variant id: $lead_id"

  # 2d) LD proxies in PLINK 1.9
  plink --bfile "work/${tag}" \
        --ld-snp "$lead_id" --r2 \
        --ld-window-kb "$WINDOW_KB" --ld-window 99999 --ld-window-r2 "$R2MIN" \
        --out "out/${tag}_r2_${R2MIN}" >/dev/null 2>&1

  # 2e) append proxies (lead can appear as SNP_A or SNP_B) to combined table
  if [ -s "out/${tag}_r2_${R2MIN}.ld" ]; then
    awk -v lead="$lead_id" -v snp="$snp" -v chr="$chr" -v pos="$pos" -v pop="$POP" '
      NR==1{ for(i=1;i<=NF;i++) h[$i]=i; next }
      {
        if ($h["SNP_A"]==lead){ p=$h["SNP_B"]; pc=$h["CHR_B"]; pp=$h["BP_B"] }
        else                  { p=$h["SNP_A"]; pc=$h["CHR_A"]; pp=$h["BP_A"] }
        print snp"\t"chr"\t"pos"\t"lead"\t"pop"\t"p"\t"pc"\t"pp"\t"$h["R2"]
      }' "out/${tag}_r2_${R2MIN}.ld" >> "$PROXY_TSV"
  fi
done

n=$(( $(wc -l < "$PROXY_TSV") - 1 ))
echo
echo "[*] DONE. $n proxy rows -> $PROXY_TSV"
