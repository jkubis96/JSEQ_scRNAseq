.PHONY: install-deps format lint check all

install-deps:
	# R dependencies
	Rscript -e "if (!requireNamespace('styler', quietly = TRUE)) install.packages('styler')"
	# Python dependencies
	pip install --quiet black isort pylint
	# Shell formatter (shfmt)
	if ! command -v shfmt >/dev/null 2>&1; then \
		curl -sS https://webi.sh/shfmt | sh; \
	fi



format:
	Rscript -e "styler::style_dir('scripts')"
	isort scripts
	black scripts
	shfmt -w scripts





all: install-deps format



