.PHONY: help check-uv fetch format format-check lint mypy python sync tests

ifeq ($(OS),Windows_NT)
  UV_VERSION := $(shell uv --version 2> NUL)
  PYTHON_VERSION := python
else
  UV_VERSION := $(shell command uv --version 2> /dev/null)
  PYTHON_VERSION := python3
endif

help: ## List all commands.
	@echo -e "\nThis makefile executes mostly uv commands. To view all uv commands available run 'uv help'."
	@echo -e "\n\033[1mAVAILABLE COMMANDS\033[0m"
	@awk 'BEGIN {FS = ":.*?## "} /^[a-zA-Z0-9 -]+:.*?## / {printf "  \033[36m%-20s\033[0m %s\n", $$1, $$2 | "sort"}' $(MAKEFILE_LIST)

check-uv: ## Check if uv is installed.
ifndef UV_VERSION
	@echo "Please install uv first. See https://docs.astral.sh/uv/ for instructions."
	@exit 1
else
	@echo "Found ${UV_VERSION}";
endif

fetch: check-uv sync ## Regenerate the dataset from USGS + EarthScope.
	uv run python fetch_events.py

format: check-uv ## Sort imports and format code.
	uv run ruff check --fix .
	uv run ruff format .

format-check: check-uv ## See what 'make format' would change.
	uv run ruff check --diff .
	uv run ruff format --diff .

lint: check-uv ## Run all linting checks.
	uv run ruff check .
	uv run ruff format --check .

mypy: check-uv ## Run the type checker.
	uv run mypy fetch_events.py tests

python: check-uv ## Start an interactive python shell in the project virtual environment.
	uv run python

sync: check-uv ## Install this project and its dependencies in a virtual environment.
	uv sync --locked

tests: check-uv mypy ## Run the dataset-integrity and unit tests.
	uv run pytest
