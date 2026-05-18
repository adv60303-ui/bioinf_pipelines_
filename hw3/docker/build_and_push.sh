#!/usr/bin/env bash

set -euo pipefail

REG="${DOCKER_USER:?Set DOCKER_USER=your-dockerhub-login}"
TAG="${DOCKER_TAG:-1.0.0}"
PUSH="${PUSH:-1}"
PLATFORM="${PLATFORM:-linux/amd64}"
SKIP_SRA="${SKIP_SRA:-0}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"

echo "Platform: ${PLATFORM}"

build() {
  local name="$1" dir="$2"
  local image="${REG}/hw3-${name}:${TAG}"
  echo "==> build ${image} (${PLATFORM})"
  docker build --platform "$PLATFORM" -t "$image" "${ROOT}/docker/${dir}"
  if [[ "$PUSH" == "1" ]]; then
    echo "==> push ${image}"
    docker push "$image"
  fi
}

if [[ "$SKIP_SRA" != "1" ]]; then
  build sra-tools sra-tools
else
  echo "==> skip hw3-sra-tools (SKIP_SRA=1)"
fi

if [[ -z "${ONLY:-}" ]] || [[ "$ONLY" == "fastqc" ]]; then
  build fastqc fastqc
fi
if [[ -n "${ONLY:-}" ]] && [[ "$ONLY" != "fastqc" ]]; then
  build "$ONLY" "$ONLY"
fi
if [[ -z "${ONLY:-}" ]]; then
  build trimmomatic  trimmomatic
  build megahit      megahit
  build bwa-samtools bwa-samtools
  build coverage     coverage
  build bcftools     bcftools
fi

echo ""
echo "Done. Run pipeline:"
echo "  nextflow run main.nf -profile container \\"
echo "    --container_registry ${REG} --container_tag ${TAG} \\"
echo "    --reads 'data/reads/DRR030302_{1,2}.fastq.gz' \\"
echo "    --reference data/ref/hiv_ref.fna --outdir results_container -resume"
