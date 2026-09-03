# MVA Hackathon 2026. Run `make help` for the list.
SHELL := /bin/bash
.DEFAULT_GOAL := help

help:  ## list targets
	@grep -hE '^[a-zA-Z0-9_-]+:.*?## ' $(MAKEFILE_LIST) \
	  | awk 'BEGIN{FS=":.*?## "}{printf "  \033[36m%-22s\033[0m %s\n", $$1, $$2}'

verify:  ## check the environment is what PROVENANCE.md claims
	@bash scripts/20_verify_environment.sh

phase0:  ## rerun Phase 0 recon end to end
	bash   scripts/00_inventory.sh
	bash   scripts/02_extract_clinical.sh
	python scripts/01_characterise.py
	python scripts/03_roh_proxy.py
	python scripts/04_verify_panel_tally.py
	python scripts/05_panel_depth_profile.py
	python scripts/06_gap_background_rate.py
	python scripts/07_write_provenance.py

resources:  ## rebuild panels and benchmark from public sources
	python scripts/10_harvest_clinvar_benchmark.py
	python scripts/11_build_mitotic_panel.py
	python scripts/12_join_constraint.py
	python scripts/13_harvest_splice_controls.py

downloads:  ## fetch reference files (large; resumable)
	mkdir -p refs benchmarks/background
	curl -L --retry 8 --retry-all-errors -C - -o refs/gnomad.v4.1.constraint_metrics.tsv \
	  https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv
	curl -L --retry 8 --retry-all-errors -C - -o refs/hgnc_complete_set.txt \
	  https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt
	curl -L --retry 8 --retry-all-errors -C - -o refs/Homo_sapiens.GRCh38.115.gtf.gz \
	  https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz
	curl -L --retry 8 --retry-all-errors -C - -o benchmarks/background/HG002_GRCh38_v4.2.1.vcf.gz \
	  https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
	@for f in refs/Homo_sapiens.GRCh38.115.gtf.gz benchmarks/background/HG002_GRCh38_v4.2.1.vcf.gz; do \
	  gzip -t $$f && echo "OK   $$f" || echo "TRUNCATED $$f - rerun 'make downloads' to resume"; \
	done

track2:  ## Track 2: direction audit, then the chemoprevention axis
	python scripts/14_track2_direction_audit.py
	python scripts/27_track2_chemoprevention.py

test:  ## run the test suite
	pytest -v

lint:  ## ruff and mypy
	ruff check src scripts tests
	mypy src

reproduce: verify phase0 resources test  ## the full path a judge should be able to run
	@echo "Reproduction complete. Compare results/ against the submitted artefact."

.PHONY: help verify phase0 resources downloads track2 test lint reproduce
