# TEBrowser

Transposable Element Family Genome Browser

A web-based visualization tool for transposable element (TE) sequences that provides a genome browser-like interface for exploring TE family consensus sequences and their characteristics.

## Overview

TEBrowser is similar to a traditional genome browser, but instead of using a reference genome, it uses transposable element family consensus sequences or Stockholm alignment files as the reference. The tool displays multiple tracks showing various characteristics across the TE reference sequence, making it both a powerful reference visualization and a useful tool for curation of putative TE families.

## Features

- **Multi-track visualization** of TE characteristics including:
  - Self-similarity analysis
  - Instance alignment tracks
  - Protein similarity data
  - Alignments to other TE families
- **Flexible input formats**: Accepts both consensus sequences and Stockholm alignment files
- **Web-based interface** built on IGV.js with custom TE-specific tracks
- **Curation support** for refining and validating TE family classifications

## Installation

git submodule init
git submodule update
cd igv.js
npm install
npm run build

## Usage

Generate a TEBrowser visualization using the `genTEBrowser.pl` script:

### From a consensus sequence:
```bash
perl genTEBrowser.pl -consensus your_te_consensus.fasta [options]
```

### From a Stockholm alignment:
```bash
perl genTEBrowser.pl -stockholm your_te_alignment.sto [options]
```

## Input Formats

- **Consensus sequences**: FASTA format containing the TE family consensus sequence
- **Stockholm files**: Multiple sequence alignment format representing aligned instances from the TE family

## Output

The tool generates web-based visualization files that can be opened in a browser to explore the TE family characteristics interactively.

## Technology

TEBrowser is built on:
- **IGV.js** - Core genome browser functionality
- **Custom tracks** - Specialized tracks for TE-specific data visualization
- **Perl** - Backend processing and file generation

## Use Cases

- **Reference visualization**: Explore the structure and characteristics of TE families
- **Family curation**: Validate and refine putative TE family classifications
- **Comparative analysis**: Compare alignments and similarities across different TE families
- **Research**: Investigate transposable element biology and evolution

## Contributing

[Contributing guidelines to be added]

## License

[License information to be added]

## Support

[Contact/support information to be added]
