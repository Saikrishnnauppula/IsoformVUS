# IsoformVUS shortcuts. `test` needs the pinned deps (requirements.txt);
# use the repo .venv when it has them, else a venv of your own.
PYTHON ?= python3
PYTEST = PYTHONPATH=src $(PYTHON) -m pytest tests/ -q

.PHONY: test benchmark figures pipeline lint format

test:
	$(PYTEST)

# Held-out benchmark (stages 10-11). Logs every access to
# results/test_access.log; refuses a second run without --allow-rerun.
benchmark:
	$(PYTHON) scripts/10_11_benchmark.py --allow-rerun

figures:
	$(PYTHON) scripts/15_make_figures.py
	$(PYTHON) scripts/16_make_missing_figures.py

pipeline:
	bash scripts/run_all.sh

lint:
	ruff check scripts/ src/ tests/
	ruff format --check scripts/ src/ tests/
	black --check --exclude 'scripts/_archived_v1_root_18_exon_architecture_nmd.py' scripts/ src/ tests/

format:
	ruff check --fix scripts/ src/ tests/
	ruff format scripts/ src/ tests/
	black --exclude 'scripts/_archived_v1_root_18_exon_architecture_nmd.py' scripts/ src/ tests/
