# MVA Hackathon 2026. Run `make help` for the list.
SHELL := /bin/bash
.DEFAULT_GOAL := help

help:  ## list targets
	@grep -hE '^[a-zA-Z_-]+:.*?## ' $(MAKEFILE_LIST) \
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

test:  ## run the test suite
	pytest -v

lint:  ## ruff and mypy
	ruff check src scripts tests
	mypy src

reproduce: verify phase0 resources test  ## the full path a judge should be able to run
	@echo "Reproduction complete. Compare results/ against the submitted artefact."

.PHONY: help verify phase0 resources test lint reproduce
