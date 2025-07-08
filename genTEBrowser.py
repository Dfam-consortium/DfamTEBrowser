#!/usr/bin/python3
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
import pprint

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

def get_tool_path(tool_name, env_var, config, cli_arg=None, default=''):
    tool_path = cli_arg \
        or os.environ.get(env_var) \
        or (config.get('tools', tool_name, fallback=None) if config else None) \
        or shutil.which(tool_name) \
        or default
    if os.path.isfile(tool_path) and os.access(tool_path, os.X_OK):
        return tool_path
    raise ValueError(f"{tool_name} tool not found. Please use --{tool_name} command line option, set {env_var} environment variable or provide a config.ini file with the correct path.")

def generate_cram_from_stockholm(stk_path, fasta_path, cram_path, crai_path, samtools_path="samtools"):
    """
    Generate CRAM and CRAI files from Stockholm/SAM and reference FASTA.
    """

    stk_to_sam(stk_path, sam_path="tmpSeed.sam", ref_fa_path="tmpSeedCons.fa")

    if not os.path.exists("tmpSeed.sam"):
        print("No SAM file (tmpSeed.sam) found for CRAM generation.")
        return False

    if os.path.exists("tmpSeedCons.fa"):
        os.remove("tmpSeedCons.fa")

    try:
        # Create CRAM
        subprocess.run([
            samtools_path, "view", "-C", "-T", fasta_path, "-o", cram_path, "tmpSeed.sam"
        ], check=True)

        # Create CRAI index
        subprocess.run([
            samtools_path, "index", cram_path
        ], check=True)

        #os.remove("tmpSeed.sam")
        if os.path.exists(fasta_path + ".fai"):
            os.remove(fasta_path + ".fai")

        if os.path.exists(cram_path) and os.path.exists(crai_path):
            print(f"Created CRAM: {cram_path} and CRAI: {crai_path}")
            return True
    except Exception as e:
        print(f"Error generating CRAM: {e}")

    return False

def download_stockholm(accession, out_file):
    url = f"https://www.dfam.org/api/families/{accession}/seed?format=stockholm"
    print(f"  - Fetching {url} ...")
    response = requests.get(url)
    response.raise_for_status()
    with open(out_file, "w") as fh:
        fh.write(response.text)

def fix_stk_reference(stk_path):
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


def parse_fasta(fasta_path):
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

def parse_stockholm(stk_path):
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

def preprocess_input(input_arg, tmp_dir, samtools_path=None):
    accession_re = re.compile(r"^D[FR]\d{9}(\.\d+)?$")
    filetype, seq_id, seq_desc, consensus, fasta_path, stk_path = None, None, None, None, None, None

    if accession_re.match(input_arg):
        print("Looks like a Dfam accession (DR######### or DF#########).")
        stk_path = os.path.join(tmp_dir, "tmpAnnotSeqDfamSeed.stk")
        download_stockholm(input_arg, stk_path)
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
        raise ValueError(f"Input {input_arg} is neither a valid file nor a Dfam accession.")

    if filetype == "fasta":
        seq_id, seq_desc, consensus = parse_fasta(input_arg)
        fasta_path = input_arg
        stk_path = None
    elif filetype == "stockholm":
        seq_id, seq_desc, rf = parse_stockholm(input_arg)
        if rf and all(c in "xX" for c in rf):
            fix_stk_reference(input_arg)
            seq_id, seq_desc, rf = parse_stockholm(input_arg)
        consensus = rf
        fasta_path = os.path.join(tmp_dir, "tmpConsensus.fa")
        with open(fasta_path, "w") as fh:
            fh.write(f">{seq_id}\n{rf}\n")
        stk_path = input_arg

        # Generate CRAM files in temp directory
        cram_path = os.path.join(tmp_dir, "tmpSeed.cram")
        crai_path = os.path.join(tmp_dir, "tmpSeed.crai")
        generate_cram_from_stockholm(stk_path, fasta_path, cram_path, crai_path, samtools_path)
    else:
        raise ValueError(f"Unknown file type for {input_arg}")

    return filetype, seq_id, seq_desc, consensus, fasta_path, stk_path

def get_base_url(output_dir, base_url_arg):
    """
    Determine the base URL for serving files.
    Priority: CLI arg > file:// URL for absolute path > relative path assumption
    """
    if base_url_arg:
        return base_url_arg.rstrip('/')

    abs_output_dir = os.path.abspath(output_dir)
    return f"file://{abs_output_dir}"

def move_output_files(tmp_dir, output_dir, seq_id, base_url):
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
        ("tmpSeed.cram", f"seed.cram"),
        ("tmpSeed.cram.crai", f"seed.cram.crai"),
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

def update_track_urls(tracks, file_mappings, base_url):
    """
    Update track URLs to point to the final output location.
    """
    for track in tracks:
        if track.get("type") == "alignment":
            # Update CRAM track URLs
            if "tmpSeed.cram" in file_mappings:
                track["url"] = file_mappings["tmpSeed.cram"]["url"]
            if "tmpSeed.cram.crai" in file_mappings:
                track["indexURL"] = file_mappings["tmpSeed.cram.crai"]["url"]
            if "tmpConsensus.fa" in file_mappings:
                track["fastaURL"] = file_mappings["tmpConsensus.fa"]["url"]

def main():
    parser = argparse.ArgumentParser(
        description="GenTEBrowser Py: Generate Dfam TE Browser Input"
    )
    parser.add_argument("input", help="Sequence file or Dfam accession (DF#########)")
    parser.add_argument("--config", help="Path to config.ini file for tool locations", default=None)
    parser.add_argument("--rmblast-dir", help="Path to rmblast directory", default=None)
    parser.add_argument("--samtools", help="Path to samtools binary", default=None)
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
    print("INSTALL DIR ==", install_dir)
    if not os.path.exists(os.path.join(install_dir, "Libraries", "Dfam-curated.fa")):
        print("Please finish installation of DfamTEBrowser by downloading the DfamLib/Dfam-curated.fa file (see README.md).")
        sys.exit(1)
    if not os.path.exists(os.path.join(install_dir, "Libraries", "RepeatPeps.lib")):
        print("Please finish installation of DfamTEBrowser by downloading the DfamLib/RepeatPeps.lib file (see README.md).")
        sys.exit(1)

    matrix_dir = os.path.join(install_dir, "Matrices/ncbi/nt")
    rmblast_dir =     os.environ.get('RMBLAST_DIR') \
                  or (config.get('tools', 'RMBLAST_DIR', fallback=None) if config else None) \
                  or os.path.dirname(shutil.which('rmblastn')) if shutil.which('rmblastn') else None \
                  or "/usr/local/rmblast/bin"
    rmblastn = os.path.join(rmblast_dir,"rmblastn")
    makeblastdb = os.path.join(rmblast_dir, "makeblastdb")
    blastx = os.path.join(rmblast_dir, "blastx")
    if not os.path.exists(rmblastn) or not os.path.exists(makeblastdb) or not os.path.exists(blastx):
        print("RMBLAST tools not found. Please use -rmblast-dir command line option, set RMBLAST_DIR environment variable or provide a config.ini file with the correct paths.")
        sys.exit(1)
    samtools_path = get_tool_path("samtools", "SAMTOOLS", config, cli_arg=args.samtools, default="/usr/local/samtools/bin/samtools")
    ultra_prgm = get_tool_path("ultra", "ULTRA", config, cli_arg=args.ultra, default="/usr/local/ultra/ultra")


    output_dir = None
    if not args.output_dir:
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
        # Change to temp directory for processing
        original_dir = os.getcwd()
        os.chdir(tmpdir)

        print("#\n# genTEBrowser :  Generate Dfam TE Visualization\n#")
        print(f"#   Input              : {args.input}")
        print(f"#   Temporary Directory: {tmpdir}")
        print(f"#   Output Directory   : {output_dir}")
        print(f"#   RMBlast            : {rmblast_dir}")
        print(f"#   Samtools           : {samtools_path}")
        print(f"#   Ultra              : {ultra_prgm}")

        try:
            filetype, seq_id, seq_desc, consensus, fasta_path, stk_path = preprocess_input(args.input, tmpdir, samtools_path)
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

        # Check for CRAM files in temp directory
        cram_path = os.path.join(tmpdir, "tmpSeed.cram")
        crai_path = os.path.join(tmpdir, "tmpSeed.cram.crai")

        if os.path.exists(cram_path) and os.path.exists(crai_path):
            cram_track = {
                "name": "Seed Alignment",
                "type": "alignment",
                "format": "cram",
                "url": "seed.cram",
                "indexURL": "seed.cram.crai",
                "displayMode": "SQUISHED",
                "fastaURL": "ref.fa",
                "autoHeight": True,
            }
            tracks.insert(0, cram_track)

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
        #pprint.pprint(protein_track)

        # Determine base URL for serving files
        base_url = get_base_url(output_dir, args.base_url)

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
        file_mappings = move_output_files(tmpdir, output_dir, seq_id, base_url)
        # Copy minified IGV JS file to output directory
        igv_js_src = os.path.join(install_dir, "js", "igv.esm.min.js")
        igv_js_dest = os.path.join(output_dir, "igv.esm.min.js")
        shutil.copy(igv_js_src, igv_js_dest)

        print(f"Output files saved to {os.path.abspath(output_dir)}\n")

        if base_url.startswith("file://"):
            html_file = file_mappings.get('tmpBrowser.html', {}).get('local_path')
            if html_file:
                print(f"Open in browser: file://{os.path.abspath(html_file)}\n")

    finally:
        # Clean up temporary directory unless keep-temp is specified
        if not args.keep_temp:
            shutil.rmtree(tmpdir, ignore_errors=True)
        else:
            print(f"Temporary directory preserved: {tmpdir}")

if __name__ == "__main__":
    main()
