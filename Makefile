.PHONY: build watch

OUT_MD_FILE := ./specification.md
OUT_HTML_FILE := ./specification.html

# builds the specification once
spec-build:
	cargo spec build --output-file $(OUT_MD_FILE)
	pandoc -f commonmark $(OUT_MD_FILE) --standalone  --output $(OUT_HTML_FILE)

# watches specification-related files and rebuilds them on the fly
spec-watch:
	cargo spec watch --output-file $(OUT_MD_FILE)
