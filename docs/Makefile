SPHINXBUILD   = sphinx-build
SOURCEDIR     = .
BUILDDIR      = _build

.PHONY: help html latexpdf clean

help:
	@echo "  html      build HTML documentation"
	@echo "  latexpdf  build PDF via LaTeX"
	@echo "  clean     remove build artifacts"

html:
	$(SPHINXBUILD) -b html $(SOURCEDIR) $(BUILDDIR)/html

latexpdf:
	$(SPHINXBUILD) -b latex $(SOURCEDIR) $(BUILDDIR)/latex
	$(MAKE) -C $(BUILDDIR)/latex all-pdf

clean:
	rm -rf $(BUILDDIR)
