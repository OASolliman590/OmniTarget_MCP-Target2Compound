.PHONY: docs
docs:
	python scripts/generate_docs.py

.PHONY: lint
lint:
	@echo "Run your linter here (e.g., ruff/flake8)"

.PHONY: test
test:
	pytest -q

