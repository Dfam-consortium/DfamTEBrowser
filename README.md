# DfamTEBrowser

Transposable Element Family Genome Browser

## Overview

A web-based visualization tool for transposable element (TE) families that treats the 
TE family consensus sequence as a reference genome within a genome browser framework. 
This approach builds on the concept pioneered by the UCSC Repeat Browser [1], which 
was released in 2020 and provides visualization of sequence homology with other TE families, 
protein alignments, and ChIP-seq datasets for KRAB Zinc Finger Proteins (KZNFs).

DfamTEBrowser extends this foundation by developing a standalone viewer with enhanced 
visualization tracks that provide clearer representation of inter- and intra-family relationships, 
seed alignments, and detailed protein alignment information.  In addition, DfamTEBrowser 
builds upon a extensible visualization framework (IGV.js [2]) that allows for the further addition 
of custom tracks and data types.

1. Fernandes, Jason D., et al. "The UCSC repeat browser allows discovery and visualization of evolutionary conflict across repeat families." Mobile DNA 11 (2020): 1-12. ( https://repeatbrowser.ucsc.edu )
2. James T Robinson, Helga Thorvaldsdottir, Douglass Turner, Jill P Mesirov, igv.js: an embeddable JavaScript implementation of the Integrative Genomics Viewer (IGV), Bioinformatics, Volume 39, Issue 1, January 2023, btac830, https://doi.org/10.1093/bioinformatics/btac830


## Features

- **Multi-track visualization** of TE characteristics including:
  - Self-similarity analysis
  - TE DNA similarity analysis
  - ULTRA Tandem Repeat analysis
  - Seed alignment track (showin details of instance coverage to the consensus)
  - Protein similarity analysis
- **Flexible input formats**: Accepts both consensus sequences and Stockholm alignment files
- **Web-based interface** built on IGV.js with custom TE-specific tracks

## Installation

The DfamTEBrowser has the following dependencies:

- Python3 + Numpy + Requests
    - `pip3 install numpy`
    - `pip3 install requests`
- RMBlast ( https://www.repeatmasker.org/rmblast/ )
    - download the binary release
- ULTRA (https://github.com/TravisWheelerLab/ULTRA)
    - download source release and follow instructions for building
- Dfam Reference Files ( Dfam-curated.fa, RepeatPeps.lib )

```bash
git clone git@github.com:Dfam-consortium/DfamTEBrowser.git
cd DfamTEBrowser
cd Libraries
# get the libraries with wget or curl
wget https://www.dfam.org/releases/current/families/Dfam-RepeatMasker.lib.gz -O Dfam-curated.fa.gz
wget https://www.dfam.org/releases/current/families/RepeatPeps.lib.gz 
# or
curl https://www.dfam.org/releases/current/families/Dfam-RepeatMasker.lib.gz -o Dfam-curated.fa.gz
curl https://www.dfam.org/releases/current/families/RepeatPeps.lib.gz -o RepeatPeps.lib.gz
# then
gunzip Dfam-curated.fa.gz
gunzip RepeatPeps.lib.gz
```

MacOS Instructions

MacOS includes Gatekeeper, a security feature that by default blocks applications that are not signed by an identified developer. If you've compiled an application by hand or downloaded one from an untrusted source, you may see errors like:
```
blastx cannot be opened because the developer cannot be verified.
```

To resolve this:
1. Attempt to open the application by double-clicking on it in the Finder.  Click "Cancel" when the error dialog appears.  This step is necessary to register the app with Gatekeeper's exception list.
2. Open System Settings -> Privacy & Security and scroll down the the security section.  You should see a message like "AppName was blocked from use because it is not from an identified developer.".  Click on the "Allow Anyway" button next to this message.
or you can grant the equivalent permissions via the command line as:

```bash
xattr -d com.apple.quarantine /path/to/blastx
chmod +x /path/to/blastx
```

You may need to do this with ultra, rmblastn, and blastx in order to fully run DfamTEBrowser.

[optional]
If you want to rebuild the IGV.js javascript component, you will need to install Node.js and npm.
The following commands will pull down the Dfam modified version of IGV.js and build the minified
files from scratch.

```bash
git submodule init
git submodule update
cd dfam_igv.js
npm install
npm run build
cp dist/igv.esm.min.js ../js/igv.esm.min.js
```

## Configuration

DfamTEBrowser can be configured one of three ways (in priority order):

  1. Using command line options to specify tool dependency locations and other parameters.
  2. Using environment variables to set dependency locations.
  3. Using a config file in INI format (specfied with -config <filename>)
  4. Using your PATH to automatically locate tools.

## Usage

Generate a DfamTEBrowser visualization using the `genTEBrowser.pl` script:

### From a consensus sequence:
```bash
python3 genTEBrowser.py your_te_consensus.fasta --output-dir output
```

### From a Stockholm alignment:
```bash
python3 genTEBrowser.py your_te_alignment.stk --output-dir output
```

### Given a Dfam accession:
Directly fetches the TE family data from the Dfam database using the provided accession number.
```bash
python3 genTEBrowser.py DF000000001 --output-dir output
python3 genTEBrowser.py DR002283232 --output-dir output
```

## Output

The tool generates web-based visualization files that can be opened in a browser to explore the TE family characteristics interactively.
These files include:

  - index.html -- The primary HTML file
  - ref.fa -- The reference sequence in FASTA format
  - seed.sam -- The seed alignment in SAM format (if input contained a seed alignment)
  - igv.esm.min.js -- The Dfam extended igv.js browser codebase

The main html file loads data additional files and therefore needs CORs functionality.  It is easiest to either write files to
a directory of an existing webserver, or simply start up one to exclusively serve the output directory like so:
```
cd DfamTEBrowser
python3 genTEBrowser.py DR002283232 --output-dir output
python3 -m http.server 8000 --directory output
```
Then open your browser to http://localhost:8000

## Technology

DfamTEBrowser is built on:
- **IGV.js** - Core genome browser functionality
- **Custom tracks** - Specialized tracks for TE-specific data visualization
- **Python** - Backend processing and file generation

## Use Cases

- **Reference visualization**: Explore the structure and characteristics of TE families
- **Family curation**: Validate and refine putative TE family classifications
- **Comparative analysis**: Compare alignments and similarities across different TE families
- **Research**: Investigate transposable element biology and evolution

## License

DfamTEBrowser: CC0 1.0 Universal
IGV.js:        The MIT License (MIT)

