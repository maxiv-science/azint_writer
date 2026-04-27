# azint-writer

[![DOI](https://zenodo.org/badge/898945642.svg)](https://doi.org/10.5281/zenodo.16760895)

**azint-writer** is a Python package for writing HDF5 files using a new [NeXus](https://www.nexusformat.org/) extension tailored for azimuthal integration in diffraction experiments. It supports 1D/2D integration data, calibration info, and experiment metadata in a flexible, NeXus-compliant format.

- 📦 Install via: `conda install -c maxiv azint-writer`  
- 🔧 Features: HDF5 writing, metadata handling, integration results, and NeXus compatibility.  
- 💡 Contribute via pull requests (MIT licensed).

📚 Full docs: [maxiv-science.github.io/azint_writer](https://maxiv-science.github.io/azint_writer)

## Updating the documentation website

The website at [maxiv-science.github.io/azint_writer](https://maxiv-science.github.io/azint_writer) is built with [MkDocs](https://www.mkdocs.org/) from the `docs/` folder on `main`.

**To update the site:** edit `docs/index.md` (or other files in `docs/`) on `main` and push or merge a PR — the site will update automatically.

> **Important:** Do not edit files directly on the `gh-pages` branch it contains auto-generated HTML and any manual changes will be overwritten.
