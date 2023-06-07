.PHONY: spec-build spec-watch

OUT_MD_FILE := ./spec/specification.md
OUT_TEX_FILE := ./spec/specification.tex
OUT_PDF_FILE := ./spec/specification.pdf
OUT_SPEC_DIRECTORY := ./spec/

# builds the specification once
spec-build:
	cargo spec build --output-file $(OUT_MD_FILE)
	pandoc $(OUT_MD_FILE) --to=latex --standalone  --output $(OUT_TEX_FILE)
	pandoc $(OUT_TEX_FILE) --to=pdf --standalone  --output $(OUT_PDF_FILE)
# this gives lots of error and does not compile corretly
# pdftex -shell-escape -output-directory $(OUT_SPEC_DIRECTORY)  \\nonstopmode\\input specification.tex 

# watches specification-related files and rebuilds them on the fly
spec-watch:
	cargo spec watch --output-file $(OUT_MD_FILE)
