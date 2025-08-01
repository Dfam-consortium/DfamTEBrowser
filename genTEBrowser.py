#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
Usage: ./genTEBrowser.py [--help]
                         [--config config.ini]
                         [--rmblast-dir /path/to/rmblast/]
                         [--ultra /path/to/ultra]
                         [--threads 8]
                         [--output-dir /path/to/output/]
                         [--base-url http://localhost:8000]
                         [--igv-js-url /path/to/igv.esm.min.js]
                         [--keep-temp]
                         family_file or Dfam accession

  Generate a web-based visualization tool for transposable element (TE)
  families that treats the TE family consensus sequence as a reference
  genome within a genome browser framework.

  This tool generates a set of HTML and data files that can be opened
  locally in a browser or served via a web server. The input may be
  a TE family consensus sequence (single sequence FASTA file), a
  TE family seed alignment (Stockholm format), or a valid Dfam accession.
  In the later case, the tool will attempt to download the seed alignment
  directly from Dfam.

  Examples:

    # Dfam accessions
    ./genTEBrowser.py DF000000001
    ./genTEBrowser.py DR002283232

    # FASTA file
    ./genTEBrowser.py myTE.fasta

    # Stockholm file
    ./genTEBrowser.py myTE.stk

  The output will be saved in the current working directory by default,
  although it is recommended that you pre-configure (e.g via a config.ini)
  or set the "--output-dir" to a directory that will contain the complete
  set of files generated.

    ./genTEBrowser.py myTE.stk --output-dir /path/to/output/


  Options:
    -h, --help            show this help message and exit
    --config CONFIG       Path to config.ini file for tool locations
    --rmblast-dir RMBLAST_DIR
                          Path to rmblast directory
    --ultra ULTRA         Path to ultra binary
    --threads THREADS     Maximum number of threads to use
    --output-dir OUTPUT_DIR
                          Output directory
    --base-url BASE_URL   Base URL for serving files (e.g., http://localhost:8000 or file:///path/to/output)
    --igv-js-url IGV_JS_URL
                          IGV JS file URL/path
    --keep-temp           Keep temporary directory for debugging

SEE ALSO: related_script.py
          Dfam: http://www.dfam.org

AUTHOR(S):
    Robert Hubley <rhubley@systemsbiology.org>

LICENSE:
    This code may be used in accordance with the Creative Commons
    Zero ("CC0") public domain dedication:
    https://creativecommons.org/publicdomain/zero/1.0/

DISCLAIMER:
  This software is provided ``AS IS'' and any express or implied
  warranties, including, but not limited to, the implied warranties of
  merchantability and fitness for a particular purpose, are disclaimed.
  In no event shall the authors or the Dfam consortium members be
  liable for any direct, indirect, incidental, special, exemplary, or
  consequential damages (including, but not limited to, procurement of
  substitute goods or services; loss of use, data, or profits; or
  business interruption) however caused and on any theory of liability,
  whether in contract, strict liability, or tort (including negligence
  or otherwise) arising in any way out of the use of this software, even
  if advised of the possibility of such damage.

"""
#
# Module imports
#
import argparse
import configparser
import os
import shutil
import re
import sys
import subprocess
import tempfile
import urllib.parse
from pathlib import Path
import requests

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)),'DfamLib'))
from MultAlign import MultAlign

from stkToSam import stk_to_sam

from analysis import (
    analyze_self,
    analyze_repeat,
    analyze_protein,
    analyze_ultra,
)
from output_html import generate_html

def _usage():
    """Print out docstring as program usage"""
    # Call to help/pydoc with scriptname ( sans path and file extension )
    help(os.path.splitext(os.path.basename(__file__))[0])
    sys.exit(0)

def _get_tool_path(tool_name, env_var, config, cli_arg=None, default=''):
    tool_path = cli_arg \
        or os.environ.get(env_var) \
        or (config.get('tools', tool_name, fallback=None) if config else None) \
        or shutil.which(tool_name) \
        or default
    if os.path.isfile(tool_path) and os.access(tool_path, os.X_OK):
        return tool_path
    raise ValueError(f"{tool_name} tool not found. Please use --{tool_name} command line option, set {env_var} environment variable or provide a config.ini file with the correct path.")

def _download_stockholm(accession, out_file):
    url = f"https://www.dfam.org/api/families/{accession}/seed?format=stockholm"
    print(f"  - Fetching {url} ...")
    response = requests.get(url)
    response.raise_for_status()
    with open(out_file, "w") as fh:
        fh.write(response.text)

def _fix_stk_reference(stk_path):
    print("\n\n*********FIXING RF*************\n\n")
    tmp_fixed = stk_path + ".fixed"
    msa = MultAlign(
               Stockholm=stk_path,
               checkForIllegalChars=True,
               internalGapChar='-',
               externalGapChar=' '
           )
    #cons = msa.consensus(with_inserts=True, only_ref_sites=True)
    cons = msa.consensus(with_inserts=True)
    cons = cons.replace('-', '.')
    prefix = "#=GC RF "
    with open(stk_path) as infile, open(tmp_fixed, "w") as outfile:
        for line in infile:
            if line.startswith(prefix):
                # Replace everything after prefix with the new RF data
                line = prefix + cons + "\n"
            outfile.write(line)
    shutil.move(tmp_fixed, stk_path)


def _parse_fasta(fasta_path):
    seq_id, seq_desc, seq = None, "", []
    with open(fasta_path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if seq_id is not None:
                    raise ValueError(f"FASTA file {fasta_path} contains more than one record.")
                header = line[1:]
                fields = header.split(maxsplit=1)
                seq_id = fields[0]
                seq_desc = fields[1] if len(fields) > 1 else ""
            elif line and not line.startswith(";"):
                seq.append(line)
    if seq_id is None:
        raise ValueError(f"FASTA file {fasta_path} does not contain a sequence record.")
    return seq_id, seq_desc, "".join(seq)

def _parse_stockholm(stk_path):
    rf_found = False
    rf = ""
    acc = None
    name = None
    with open(stk_path) as fh:
        for line in fh:
            if line.startswith("#=GC RF"):
                if rf_found:
                    raise ValueError(f"STK file {stk_path} contains more than one #=GC RF line.")
                rf = line.split(maxsplit=2)[2].strip().replace(".", "")
                rf_found = True
            elif line.startswith("#=GF AC"):
                parts = line.split(maxsplit=2)
                if len(parts) > 2:
                    acc = parts[2].strip()
                elif len(parts) > 1:
                    acc = parts[1].strip()
            elif line.startswith("#=GF ID"):
                parts = line.split(maxsplit=2)
                if len(parts) > 2:
                    name = parts[2].strip()
                elif len(parts) > 1:
                    name = parts[1].strip()
    if not rf_found:
        raise ValueError(f"STK file {stk_path} does not contain a #=GC RF line.")
    seq_id = acc if acc else name
    seq_desc = name if acc and name else ""
    return seq_id, seq_desc, rf

def _preprocess_input(input_arg, tmp_dir):
    accession_re = re.compile(r"^D[FR]\d{9}(\.\d+)?$")
    filetype, seq_id, seq_desc, consensus, fasta_path, stk_path = None, None, None, None, None, None

    if accession_re.match(input_arg):
        print("Looks like a Dfam accession (DR######### or DF#########).")
        stk_path = os.path.join(tmp_dir, "tmpAnnotSeqDfamSeed.stk")
        _download_stockholm(input_arg, stk_path)
        input_arg = stk_path
        filetype = "stockholm"
    elif os.path.isfile(input_arg):
        with open(input_arg) as fh:
            for line in fh:
                if line.startswith(">"):
                    filetype = "fasta"
                    break
                elif line.startswith("# STOCKHOLM"):
                    filetype = "stockholm"
                    break
        if not filetype:
            raise ValueError(f"Could not autodetect file type of {input_arg} (no FASTA or STOCKHOLM header).")
    else:
        raise ValueError(f"Input {input_arg} is neither a valid file nor a Dfam accession.  Dfam accessions\n" \
                          "       should contain 9-digits (e.g DF000000001 or DR002283232).  If you meant to use a\n" \
                          "       file, please provide a valid FASTA or STOCKHOLM file.")

    if filetype == "fasta":
        seq_id, seq_desc, consensus = _parse_fasta(input_arg)
        fasta_path = input_arg
        stk_path = None
    elif filetype == "stockholm":
        seq_id, seq_desc, rf = _parse_stockholm(input_arg)
        if rf and all(c in "xX" for c in rf):
            _fix_stk_reference(input_arg)
            seq_id, seq_desc, rf = _parse_stockholm(input_arg)
        consensus = rf
        fasta_path = os.path.join(tmp_dir, "tmpConsensus.fa")
        with open(fasta_path, "w") as fh:
            fh.write(f">{seq_id}\n{rf}\n")
        stk_path = input_arg

        stk_to_sam(stk_path, sam_path="tmpSeed.sam", ref_fa_path="tmpSamCons.fa")

        if not os.path.exists("tmpSeed.sam"):
            print("No SAM file (tmpSeed.sam) found for CRAM generation.")
            return False

        if os.path.exists("tmpSamCons.fa"):
            os.remove("tmpSamCons.fa")
    else:
        raise ValueError(f"Unknown file type for {input_arg}")

    return filetype, seq_id, seq_desc, consensus, fasta_path, stk_path

def _get_base_url(output_dir, base_url_arg):
    """
    Determine the base URL for serving files.
    Priority: CLI arg > file:// URL for absolute path > relative path assumption
    """
    if base_url_arg:
        return base_url_arg.rstrip('/')

    abs_output_dir = os.path.abspath(output_dir)
    return f"file://{abs_output_dir}"

def _move_output_files(tmp_dir, output_dir, seq_id, base_url):
    """
    Move generated files from temp directory to output directory and update track URLs.
    Returns updated file mappings for track URL updates.
    """
    os.makedirs(output_dir, exist_ok=True)

    file_mappings = {}

    # Define files to move with their output names
    files_to_move = [
        ("tmpBrowser.html", f"index.html"),
        ("tmpConsensus.fa", f"ref.fa"),
        ("tmpSeed.sam", f"seed.sam"),
    ]

    for temp_name, output_name in files_to_move:
        temp_path = os.path.join(tmp_dir, temp_name)
        output_path = os.path.join(output_dir, output_name)

        if os.path.exists(temp_path):
            shutil.move(temp_path, output_path)
            file_mappings[temp_name] = {
                'local_path': output_path,
                'url': f"{base_url}/{output_name}"
            }

    return file_mappings

def main():
    #
    # Options processing
    #
    #   There are two ways to document usage/command line
    #   arguments using this boilerplate.  The intended way
    #   is to document using docstrings at the top of the
    #   script.  This way the pydoc docs match the output
    #   produced by '-h' or '--help' using the argparse
    #   custom action class ( _CustomUsageAction ) defined
    #   below.  If you want the paired-down argparse default
    #   instead simply remove the "add_help=False" argument
    #   to the argparse constructor below and comment out
    #   the add_argment('-h', ...) line below.
    #
    class _CustomUsageAction(argparse.Action):
        def __init__(self, option_strings, dest, default=False, required=False, help=None):
            super(_CustomUsageAction, self).__init__(
                      option_strings=option_strings, dest=dest,
                      nargs=0, const=True, default=default,
                      required=required, help=help)
        def __call__(self, parser, args, values, option_string=None):
            _usage()

    parser = argparse.ArgumentParser( add_help=False )
    parser.add_argument('-h', '--help', action=_CustomUsageAction )
    parser.add_argument("input", help="Sequence file or Dfam accession (DF#########)")
    parser.add_argument("--config", help="Path to config.ini file for tool locations", default=None)
    parser.add_argument("--rmblast-dir", help="Path to rmblast directory", default=None)
    parser.add_argument("--ultra", help="Path to ultra binary", default=None)
    parser.add_argument("--threads", help="Maximum number of threads to use", default=8)
    parser.add_argument("--output-dir", help="Output directory", default=None)
    parser.add_argument("--base-url", help="Base URL for serving files (e.g., http://localhost:8000 or file:///path/to/output)", default=None)
    parser.add_argument("--igv-js-url", help="IGV JS file URL/path", default="./igv.esm.min.js")
    parser.add_argument("--keep-temp", action="store_true", help="Keep temporary directory for debugging")
    args = parser.parse_args()

    config = None
    if args.config:
        config = configparser.ConfigParser()
        config.read(args.config)

    install_dir = os.path.dirname(os.path.abspath(__file__))
    if not os.path.exists(os.path.join(install_dir, "Libraries", "Dfam-curated.fa")):
        print("Please finish installation of DfamTEBrowser by downloading the DfamLib/Dfam-curated.fa file (see README.md).")
        sys.exit(1)
    if not os.path.exists(os.path.join(install_dir, "Libraries", "RepeatPeps.lib")):
        print("Please finish installation of DfamTEBrowser by downloading the DfamLib/RepeatPeps.lib file (see README.md).")
        sys.exit(1)

    matrix_dir = os.path.join(install_dir, "Matrices/ncbi/nt")
    rmblast_dir =     os.environ.get('RMBLAST_DIR') \
                  or (config.get('tools', 'rmblast_dir', fallback=None) if config else None) \
                  or os.path.dirname(shutil.which('rmblastn')) if shutil.which('rmblastn') else None \
                  or "/usr/local/rmblast/bin"
    rmblastn = os.path.join(rmblast_dir,"rmblastn")
    makeblastdb = os.path.join(rmblast_dir, "makeblastdb")
    blastx = os.path.join(rmblast_dir, "blastx")
    if not os.path.exists(rmblastn) or not os.path.exists(makeblastdb) or not os.path.exists(blastx):
        print("RMBLAST tools not found. Please use -rmblast-dir command line option, set RMBLAST_DIR environment variable or provide a config.ini file with the correct paths.")
        sys.exit(1)
    ultra_prgm = _get_tool_path("ultra", "ULTRA", config, cli_arg=args.ultra, default="/usr/local/ultra/ultra")


    output_dir = None
    if not args.output_dir:
        if config and config.get('paths', 'output_dir', fallback=None):
            output_dir = config.get('paths', 'output_dir')
            if not os.path.exists(output_dir):
                os.makedirs(output_dir, exist_ok=True)
        else:
            users_home = os.path.expanduser("~")
            if os.path.exists(os.path.join(users_home, "public_html")):
                if not os.path.exists(os.path.join(users_home, "public_html", "DfamTEBrowser")):
                    os.makedirs(os.path.join(users_home, "public_html", "DfamTEBrowser"), exist_ok=True)
                output_dir = os.path.join(users_home, "public_html", "DfamTEBrowser")
            else:
                output_dir = "."
    else:
        if not os.path.exists(args.output_dir):
            os.makedirs(args.output_dir, exist_ok=True)
        output_dir = args.output_dir

    try:
        # Create a real temporary directory
        tmpdir = tempfile.mkdtemp(prefix="gentebrowser_")
        print("#\n# genTEBrowser :  Generate Dfam TE Visualization\n#")
        print(f"#   Input              : {args.input}")
        print(f"#   Temporary Directory: {tmpdir}")
        print(f"#   Output Directory   : {output_dir}")
        print(f"#   RMBlast            : {rmblast_dir}")
        print(f"#   Ultra              : {ultra_prgm}")

        input_val = args.input
        if os.path.isfile(input_val):
            input_val = os.path.abspath(input_val)

        # Change to temp directory for processing
        original_dir = os.getcwd()
        os.chdir(tmpdir)

        try:
            filetype, seq_id, seq_desc, consensus, fasta_path, stk_path = _preprocess_input(input_val, tmpdir)
            print("Filetype:", filetype)
            print("Sequence ID:", seq_id)
            #print("#   Description :", seq_desc)
            #print("#   Consensus:", consensus)
            #print("#   FASTA path:", fasta_path)
            #print("#   STK path:", stk_path)
        except Exception as e:
            print("Error:", e)
            sys.exit(1)

        tracks = []

        # If we have a seed alignment build alignment track
        sam_path = os.path.join(tmpdir, "tmpSeed.sam")
        if os.path.exists(sam_path):
            aln_track = {
                "name": "Seed Alignment",
                "type": "seedalign",
                "format": "sam",
                "sourceType": "sam",
                "url": "seed.sam",
                "displayMode": "SQUISHED",
                "fastaURL": "ref.fa",
                "autoHeight": True,
            }
            tracks.insert(0, aln_track)

        ultra_track = analyze_ultra(fasta_path, ultra_prgm, args.threads)
        if ultra_track and ultra_track.get("features"):
            if len(ultra_track["features"]) > 0:
                tracks.append(ultra_track)

        self_track = analyze_self(fasta_path, rmblastn, matrix_dir, args.threads)
        if self_track and self_track.get("features"):
            if len(self_track["features"]) > 0:
                tracks.append(self_track)

        dfam_curated_path  = os.path.join(install_dir, "Libraries", "Dfam-curated.fa")
        if not os.path.exists(dfam_curated_path + ".nsq"):
            cmd = f"{makeblastdb} -in {dfam_curated_path} -dbtype nucl"
            print(f"Building blastdb: {cmd}")
            subprocess.run(cmd, shell=True, check=True)

        repeat_track = analyze_repeat(
            fasta_path,
            dfam_curated_path,
            rmblastn,
            matrix_dir, args.threads
        )
        if repeat_track and repeat_track.get("features"):
            if len(repeat_track["features"]) > 0:
                tracks.append(repeat_track)


        te_proteins_path  = os.path.join(install_dir, "Libraries", "RepeatPeps.lib")
        if not os.path.exists(te_proteins_path + ".psq"):
            cmd = f"{makeblastdb} -in {te_proteins_path} -dbtype prot"
            print(f"Building blastdb: {cmd}")
            subprocess.run(cmd, shell=True, check=True)

        protein_track = analyze_protein(
            fasta_path,
            os.path.join(install_dir, "Libraries", "RepeatPeps.lib"),
            blastx, args.threads
        )
        if protein_track and protein_track.get("features"):
            if len(protein_track["features"]) > 0:
                tracks.append(protein_track)

        # Determine base URL for serving files
        base_url = _get_base_url(output_dir, args.base_url)

        print("Generating HTML output...")
        generate_html(
            seq_id=seq_id,
            title=seq_id + " Annotation",
            tracks=tracks,
            output_html_path="tmpBrowser.html",
            igv_js_url=args.igv_js_url,
        )

        # Return to original directory
        os.chdir(original_dir)

        # Move files to output directory
        file_mappings = _move_output_files(tmpdir, output_dir, seq_id, base_url)
        # Copy minified IGV JS file to output directory
        igv_js_src = os.path.join(install_dir, "js", "igv.esm.min.js")
        igv_js_dest = os.path.join(output_dir, "igv.esm.min.js")
        shutil.copy(igv_js_src, igv_js_dest)

        print(f"Output files saved to {os.path.abspath(output_dir)}\n")

        #if base_url.startswith("file://"):
        #    html_file = file_mappings.get('tmpBrowser.html', {}).get('local_path')
        #    if html_file:
        #        print(f"Open in browser: file://{os.path.abspath(html_file)}\n")
        print("If the output directory is not already served by an existing webserver")
        print("you can view your results by starting one in the output directory like so:")
        print("   python3 -m http.server 8000 --directory ", output_dir)
        print("Then open your browser and navigate to:")
        print(f"   http://localhost:8000\n\n")

    finally:
        # Clean up temporary directory unless keep-temp is specified
        if not args.keep_temp:
            shutil.rmtree(tmpdir, ignore_errors=True)
        else:
            print(f"Temporary directory preserved: {tmpdir}")

if __name__ == "__main__":
    main()
