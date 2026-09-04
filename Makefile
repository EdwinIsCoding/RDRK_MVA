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
	curl -L --retry 8 --retry-all-errors -C - -o refs/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
	  https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
	@for f in refs/Homo_sapiens.GRCh38.115.gtf.gz benchmarks/background/HG002_GRCh38_v4.2.1.vcf.gz refs/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz; do \
	  gzip -t $$f && echo "OK   $$f" || echo "TRUNCATED $$f - rerun 'make downloads' to resume"; \
	done

downloads-track2:  ## the ONLY download Track 2 needs (HGNC, 16 MB)
	mkdir -p refs
	curl -L --retry 8 --retry-all-errors -C - -o refs/hgnc_complete_set.txt \
	  https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt
	@wc -l refs/hgnc_complete_set.txt

reproduce-track2: downloads-track2 track2 scalability structural-check resource test  ## Track 2 end to end, NO patient data required
	@echo
	@echo "Track 2 reproduced from public databases only."
	@echo "No file under data/ was read. Compare results/summaries/ against"
	@echo "the tables in submission/track2_nexusdwin_report.md."

track2:  ## Track 2: direction audit, then the chemoprevention axis
	python scripts/14_track2_direction_audit.py
	python scripts/27_track2_chemoprevention.py
	python scripts/28_track2_axis_availability.py

scalability:  ## run the same Track 2 pipeline on two comparator diseases
	python scripts/31_track2_scalability.py

resource:  ## publish the reusable directional availability table
	python scripts/30_publish_directional_availability.py

mito-axis:  ## follow through on the mitochondrial axis: drugs, screen, tissue gate
	python scripts/40_mitochondrial_axis_followthrough.py

mutect2:  ## Arm D re-run with a dedicated somatic caller (needs tools/, refs/)
	python scripts/39_mutect2_mosaic.py

signature:  ## LINCS signature reversal against a labelled BubR1-hypomorph proxy
	python scripts/41_signature_reversal.py

sv-screen:  ## calibrated breakpoint screen over the known MVA genes
	python scripts/36_sv_screen_panel.py

precedent:  ## curated MVA1 variants in the BubR1 kinase domain
	python scripts/37_kinase_domain_precedent.py

chemoprev-lit:  ## does chemoprevention evidence exist for this disease?
	python scripts/38_chemoprevention_literature.py

predictors:  ## in silico predictor panel for the two BUB1B alleles
	python scripts/34_missense_predictor_panel.py

dataset-revision:  ## recover the HuggingFace dataset revision for PROVENANCE
	python scripts/35_capture_dataset_revision.py

delete-plan:  ## dry run of the ETHICS 3b deletion obligation (deletes nothing)
	python scripts/33_delete_challenge_data.py

track2-drift:  ## has the live evidence moved since the Track 2 report was pinned?
	python scripts/29_track2_drift_check.py

pitch-wordcount:  ## spoken word count and runtime of the Track 2 pitch script
	@python3 -c "import pathlib; \
	t=pathlib.Path('submission/track2_nexusdwin_pitch.md').read_text(); \
	w=sum(len(l[1:].split()) for l in t.splitlines() if l.startswith('> ')); \
	print(f'{w} spoken words, {round(w/140*60)}s at 140 wpm, {round(w/160*60)}s at 160 wpm')"

test:  ## run the test suite
	pytest -v

lint:  ## ruff and mypy
	ruff check src scripts tests
	mypy src

reproduce: verify phase0 resources test  ## the full path a judge should be able to run
	@echo "Reproduction complete. Compare results/ against the submitted artefact."

.PHONY: help verify phase0 resources downloads downloads-track2 reproduce-track2 track2 scalability structural-check resource mito-axis mutect2 signature sv-screen precedent chemoprev-lit predictors dataset-revision delete-plan track2-drift pitch-wordcount test lint reproduce
