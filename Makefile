.DEFAULT_GOAL := help
UV_RUN = uv run --locked --no-default-groups
UV_TEST = $(UV_RUN) --group test --extra report
UV_QUALITY = $(UV_RUN) --group quality --extra report
PYTHON_PATHS = src tests scripts
DOCKER_IMAGE ?= muconespan:local

.PHONY: help init install-uv install dev conda-setup test test-fast test-unit test-int lint lint-fix format format-check type-check file-size workflow-check quality check ci-check docs-check security-check build-check hooks clean generate-testdata lock sync docker-build docker-test docker-smoke

help:  ## Show available commands
	@awk 'BEGIN {FS = ":.*## "} /^[a-zA-Z_-]+:.*## / {printf "%-22s %s\n", $$1, $$2}' $(MAKEFILE_LIST)

init: dev hooks  ## Set up development dependencies and Git hooks (requires uv)

install-uv:  ## Show the uv installation documentation
	@echo "Install uv: https://docs.astral.sh/uv/getting-started/installation/"

install:  ## Install the locked runtime and report dependencies
	uv sync --locked --no-default-groups --extra report

dev:  ## Install all locked development, report, and documentation dependencies
	uv sync --locked --all-groups --all-extras

conda-setup:  ## Create bioinformatics tools environment
	conda env create -f conda/environment.yml

sync: dev  ## Synchronize the development environment without changing uv.lock

lock:  ## Resolve dependency changes deliberately
	uv lock

hooks: dev  ## Install pre-commit and pre-push hooks
	$(UV_RUN) --group dev pre-commit install --install-hooks

lint:  ## Lint package, tests, and helper scripts
	$(UV_QUALITY) ruff check $(PYTHON_PATHS)

lint-fix:  ## Apply Ruff lint fixes
	$(UV_QUALITY) ruff check --fix $(PYTHON_PATHS)

format:  ## Format package, tests, and helper scripts
	$(UV_QUALITY) ruff format $(PYTHON_PATHS)

format-check:  ## Verify formatting without editing files
	$(UV_QUALITY) ruff format --check $(PYTHON_PATHS)

type-check:  ## Check package and helper script types
	$(UV_QUALITY) mypy src/muc_one_span scripts

file-size:  ## Reject source/configuration files with 650 or more physical lines
	$(UV_RUN) python scripts/check_file_size.py

workflow-check:  ## Validate GitHub Actions syntax and expressions
	$(UV_QUALITY) actionlint -shellcheck=

quality: lint format-check type-check file-size workflow-check  ## Run the same static checks as CI

test:  ## Run all tests (external tests skip when prerequisites are absent)
	$(UV_TEST) pytest

test-fast:  ## Run unit tests without coverage
	$(UV_TEST) pytest tests/unit --no-cov -x

test-unit:  ## Run unit tests with the CI coverage gate
	$(UV_TEST) pytest tests/unit --cov-fail-under=80

test-int:  ## Run bioinformatics tool integration tests
	$(UV_TEST) pytest tests/integration -m integration --no-cov

ci-check: quality test-unit  ## Run CI static checks and unit coverage gate locally

check: ci-check docs-check security-check build-check  ## Run all portable release checks

docs-check:  ## Build documentation with warnings treated as errors
	$(UV_RUN) --extra docs mkdocs build --strict

security-check:  ## Audit every locked extra against published Python advisories
	uv export --locked --all-groups --all-extras --no-emit-project --format requirements-txt --output-file .audit-requirements.txt --quiet
	$(UV_RUN) --group security pip-audit --disable-pip --require-hashes -r .audit-requirements.txt

build-check:  ## Build wheel/sdist and verify bundled resources in an isolated environment
	uv build --clear
	$(UV_RUN) python scripts/check_distribution.py

docker-build:  ## Build the production runtime image locally
	docker build -f docker/Dockerfile -t $(DOCKER_IMAGE) .

docker-test:  ## Test the runtime in BuildKit without exporting/loading an image
	docker buildx build --target test --output type=cacheonly -f docker/Dockerfile .

docker-smoke:  ## Exercise an already-built image via Docker, including its entrypoint
	bash scripts/check_container.sh $(DOCKER_IMAGE)

generate-testdata:  ## Generate MucOneUp test data (requires external tools)
	$(UV_RUN) python scripts/generate_testdata.py

clean:  ## Remove build, coverage and documentation artifacts
	rm -rf build dist htmlcov site .coverage coverage.xml .audit-requirements.txt
