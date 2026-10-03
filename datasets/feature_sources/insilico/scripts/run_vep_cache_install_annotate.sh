#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd)"
VEP_RAW_DIR="$ROOT_DIR/datasets/feature_sources/insilico/raw/vep"
VEP_CACHE_DIR="$VEP_RAW_DIR/cache"
VEP_TARBALL="$VEP_RAW_DIR/downloads/homo_sapiens_vep_116_GRCh38.tar.gz"
SPLICEAI_VCF="$ROOT_DIR/datasets/feature_sources/insilico/raw/spliceai/downloads/spliceai_scores.raw.snv.ensembl_mane_v1.4.grch38.vcf.gz"
VEP_PLUGIN_DIR="$VEP_RAW_DIR/plugins"
VEP_FASTA_DIR="$VEP_RAW_DIR/fasta"
VEP_FASTA_GZ="$VEP_FASTA_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
VEP_FASTA="$VEP_FASTA_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
LOG_DIR="$VEP_RAW_DIR/logs"
INTERIM_DIR="$ROOT_DIR/datasets/feature_sources/insilico/interim/vep"
MATRIX="${VEP_MATRIX:-$ROOT_DIR/datasets/modeling/interim/final_modeling_table_local_features_gnomad_foldx.tsv}"

CONDA="$HOME/opt/miniforge3/bin/conda"
MAMBA="$HOME/opt/miniforge3/bin/mamba"
ENV_NAME="${VEP_ENV_NAME:-vep-116}"
CONDA_SUBDIR_OVERRIDE="${VEP_CONDA_SUBDIR:-}"
VEP_ENV_PREFIX="$HOME/opt/miniforge3/envs/$ENV_NAME"
ENABLE_HGVS="${VEP_ENABLE_HGVS:-0}"
ENABLE_SPLICEAI="${VEP_ENABLE_SPLICEAI:-0}"
RUN_SUFFIX="${VEP_RUN_SUFFIX:-vep}"
OUT_MATRIX="${VEP_OUT_MATRIX:-$ROOT_DIR/datasets/modeling/interim/final_modeling_table_local_features_gnomad_foldx_${RUN_SUFFIX}.tsv}"
RUN_PREFIX=()
if [[ "$CONDA_SUBDIR_OVERRIDE" == "osx-64" ]]; then
  RUN_PREFIX=(arch -x86_64)
fi

run_vep() {
  PATH="$VEP_ENV_PREFIX/bin:$PATH" "${RUN_PREFIX[@]}" "$VEP_ENV_PREFIX/bin/vep" "$@"
}

mkdir -p "$VEP_CACHE_DIR" "$LOG_DIR" "$INTERIM_DIR" "$VEP_FASTA_DIR"
exec > >(tee -a "$LOG_DIR/vep_cache_install_annotate.log") 2>&1

echo "[$(date)] Starting VEP setup + annotation"
echo "ROOT_DIR=$ROOT_DIR"

if [[ ! -s "$VEP_TARBALL" ]]; then
  echo "ERROR: VEP cache tarball missing: $VEP_TARBALL" >&2
  exit 1
fi

if [[ ! -d "$VEP_CACHE_DIR/homo_sapiens/116_GRCh38" ]]; then
  echo "[$(date)] Unpacking VEP cache into $VEP_CACHE_DIR"
  tar -xzf "$VEP_TARBALL" -C "$VEP_CACHE_DIR"
else
  echo "[$(date)] VEP cache already unpacked: $VEP_CACHE_DIR/homo_sapiens/116_GRCh38"
fi

if ! run_vep --help >/dev/null 2>&1; then
  echo "[$(date)] Installing VEP into conda env $ENV_NAME"
  if ! "$CONDA" env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
    if [[ -n "$CONDA_SUBDIR_OVERRIDE" ]]; then
      CONDA_SUBDIR="$CONDA_SUBDIR_OVERRIDE" "$MAMBA" create -y -n "$ENV_NAME" -c conda-forge -c bioconda ensembl-vep=116 perl-bio-db-hts tabix htslib
      "$CONDA" config --env --set subdir "$CONDA_SUBDIR_OVERRIDE" -n "$ENV_NAME" || true
    else
      "$MAMBA" create -y -n "$ENV_NAME" -c conda-forge -c bioconda ensembl-vep=116 perl-bio-db-hts tabix htslib
    fi
  else
    if [[ -n "$CONDA_SUBDIR_OVERRIDE" ]]; then
      CONDA_SUBDIR="$CONDA_SUBDIR_OVERRIDE" "$MAMBA" install -y -n "$ENV_NAME" -c conda-forge -c bioconda ensembl-vep=116 perl-bio-db-hts tabix htslib
    else
      "$MAMBA" install -y -n "$ENV_NAME" -c conda-forge -c bioconda ensembl-vep=116 perl-bio-db-hts tabix htslib
    fi
  fi
else
  echo "[$(date)] VEP already available in conda env $ENV_NAME"
fi

echo "[$(date)] VEP version"
run_vep --help | head -20 || true

echo "[$(date)] Building VCF input from current matrix"
python3 "$ROOT_DIR/datasets/feature_sources/insilico/scripts/make_vep_input_vcf.py" \
  --matrix "$MATRIX" \
  --out-vcf "$INTERIM_DIR/final_modeling_table_local_features_gnomad_foldx.vep_input.vcf"

cat > "$INTERIM_DIR/vep_smoke_input.vcf" <<'VCF'
##fileformat=VCFv4.2
##reference=GRCh38
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	74236138	1-74236138-T-C	T	C	.	.	.
VCF

COMMON_ARGS=(
  --offline
  --cache
  --dir_cache "$VEP_CACHE_DIR"
  --species homo_sapiens
  --assembly GRCh38
  --format vcf
  --tab
  --force_overwrite
  --no_stats
  --symbol
  --canonical
  --mane
  --protein
  --biotype
  --numbers
  --variant_class
)

FIELDS="Uploaded_variation,Location,Allele,Gene,Feature,Feature_type,Consequence,cDNA_position,CDS_position,Protein_position,Amino_acids,Codons,Existing_variation,IMPACT,SYMBOL,SYMBOL_SOURCE,HGNC_ID,BIOTYPE,CANONICAL,MANE_SELECT,MANE_PLUS_CLINICAL,VARIANT_CLASS"

if [[ "$ENABLE_HGVS" == "1" ]]; then
  if [[ ! -s "$VEP_FASTA" ]]; then
    if [[ ! -s "$VEP_FASTA_GZ" ]]; then
      echo "[$(date)] Downloading GRCh38 primary assembly FASTA for HGVS"
      (
        cd "$VEP_FASTA_DIR"
        curl --http1.1 --location --fail --continue-at - --remote-name \
          "https://ftp.ensembl.org/pub/release-116/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
      )
    fi
    echo "[$(date)] Decompressing FASTA for VEP HGVS"
    gzip -dkf "$VEP_FASTA_GZ"
  fi
  if [[ ! -s "$VEP_FASTA.fai" ]]; then
    echo "[$(date)] Indexing FASTA with samtools faidx"
    PATH="$VEP_ENV_PREFIX/bin:$PATH" "${RUN_PREFIX[@]}" "$VEP_ENV_PREFIX/bin/samtools" faidx "$VEP_FASTA"
  fi
  COMMON_ARGS+=(--hgvs)
  COMMON_ARGS+=(--fasta "$VEP_FASTA")
  FIELDS="$FIELDS,HGVSc,HGVSp"
else
  echo "[$(date)] HGVS disabled for this run; offline VEP requires a FASTA for --hgvs"
fi

PLUGIN_ARGS=()
if [[ "$ENABLE_SPLICEAI" == "1" && -s "$SPLICEAI_VCF" && -s "$SPLICEAI_VCF.tbi" && -s "$VEP_PLUGIN_DIR/SpliceAI.pm" ]]; then
  PLUGIN_ARGS=(--dir_plugins "$VEP_PLUGIN_DIR" --plugin "SpliceAI,snv=$SPLICEAI_VCF,split_output=1")
  FIELDS="$FIELDS,SpliceAI_pred_SYMBOL,SpliceAI_pred_DS_AG,SpliceAI_pred_DS_AL,SpliceAI_pred_DS_DG,SpliceAI_pred_DS_DL,SpliceAI_pred_DP_AG,SpliceAI_pred_DP_AL,SpliceAI_pred_DP_DG,SpliceAI_pred_DP_DL"
  echo "[$(date)] SpliceAI plugin enabled for SNVs"
else
  echo "[$(date)] SpliceAI plugin disabled for this run"
fi

echo "[$(date)] Running VEP smoke test"
if [[ ${#PLUGIN_ARGS[@]} -gt 0 ]]; then
  run_vep \
    "${COMMON_ARGS[@]}" \
    "${PLUGIN_ARGS[@]}" \
    --fields "$FIELDS" \
    --input_file "$INTERIM_DIR/vep_smoke_input.vcf" \
    --output_file "$INTERIM_DIR/vep_smoke_output.tsv"
else
  run_vep \
    "${COMMON_ARGS[@]}" \
    --fields "$FIELDS" \
    --input_file "$INTERIM_DIR/vep_smoke_input.vcf" \
    --output_file "$INTERIM_DIR/vep_smoke_output.tsv"
fi

echo "[$(date)] Running full VEP annotation"
if [[ ${#PLUGIN_ARGS[@]} -gt 0 ]]; then
  run_vep \
    "${COMMON_ARGS[@]}" \
    "${PLUGIN_ARGS[@]}" \
    --fields "$FIELDS" \
    --fork 4 \
    --input_file "$INTERIM_DIR/final_modeling_table_local_features_gnomad_foldx.vep_input.vcf" \
    --output_file "$INTERIM_DIR/final_modeling_table_local_features_gnomad_foldx.${RUN_SUFFIX}.tsv"
else
  run_vep \
    "${COMMON_ARGS[@]}" \
    --fields "$FIELDS" \
    --fork 4 \
    --input_file "$INTERIM_DIR/final_modeling_table_local_features_gnomad_foldx.vep_input.vcf" \
    --output_file "$INTERIM_DIR/final_modeling_table_local_features_gnomad_foldx.${RUN_SUFFIX}.tsv"
fi

echo "[$(date)] Parsing VEP output"
python3 "$ROOT_DIR/datasets/feature_sources/insilico/scripts/parse_vep_tab.py" \
  --matrix "$MATRIX" \
  --vep-tab "$INTERIM_DIR/final_modeling_table_local_features_gnomad_foldx.${RUN_SUFFIX}.tsv" \
  --out-features "$INTERIM_DIR/final_modeling_table_local_features_gnomad_foldx.${RUN_SUFFIX}_features.tsv" \
  --summary-json "$INTERIM_DIR/final_modeling_table_local_features_gnomad_foldx.${RUN_SUFFIX}_features.summary.json"

echo "[$(date)] Merging VEP features into modeling matrix"
python3 "$ROOT_DIR/datasets/feature_sources/insilico/scripts/merge_vep_features_into_matrix.py" \
  --matrix "$MATRIX" \
  --vep-features "$INTERIM_DIR/final_modeling_table_local_features_gnomad_foldx.${RUN_SUFFIX}_features.tsv" \
  --out-matrix "$OUT_MATRIX" \
  --summary-json "$ROOT_DIR/datasets/modeling/interim/final_modeling_table_local_features_gnomad_foldx_${RUN_SUFFIX}.summary.json"

echo "[$(date)] Done. Output matrix: $OUT_MATRIX"
