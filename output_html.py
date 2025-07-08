import os
import copy
import json
from string import Template
from colors import DISTINCT_COLORS

def html_escape(text):
    return (str(text).replace("&", "&amp;")
                     .replace("<", "&lt;")
                     .replace(">", "&gt;")
                     .replace('"', "&quot;")
                     .replace("'", "&#39;"))

HTML_TEMPLATE = Template("""
<html lang="en">
<head>
    <meta charset="utf-8">
    <title>$title</title>
</head>
<body>
    <h1>$title</h1>
    <div id="igvDiv" style="padding:10px; border:1px solid lightgray;"></div>

    <script type="module">
        import igv from "$igv_js_url";

        const options = {
            reference: { fastaURL: "$ref_fasta_url", indexed: false },
            locus: ["$seq_id"],
            tracks: [
$track_js
            ]
        };

        igv.createBrowser(document.getElementById("igvDiv"), options)
            .then(browser => console.log("IGV browser created."))
            .catch(error => console.error("Error creating IGV browser:", error));
    </script>
</body>
</html>
""")

def feature_to_js(feature, seq_id, track_type=None):
    if feature.get('type') == 'protein':
        f = feature
        subject_id = f['name']
        # Normalize coordinates and strand
        strand = '-' if (f.get('orient') == '-') else '+'

        start = (f['ref_start'] or 1) - 1
        end = f['ref_end'] or start + 1
        if start > end:
            start, end = end, start
        cdStart, cdEnd = start, end

        exon_json = {
            "start": start,
            "end": end,
            "cdStart": cdStart,
            "cdEnd": cdEnd
        }

        oChromStart = (f.get('cons_start') or 1) - 1
        oChromEnd = f.get('cons_end') or (oChromStart + 1)
        oStrand = f.get('orient') if f.get('orient') else "+"
        oChromSize = abs(oChromEnd - oChromStart)
        oChromStarts = f"{oChromStart},"
        oSequence = f.get('oseq', "")
        color = feature.get('color','#e6194b')

        return (
            "{ chr: \"" + html_escape(seq_id) + "\""
            + f", start: {start}"
            + f", end: {end}"
            + f", name: \"{html_escape(subject_id)}\""
            + f", score: {f.get('score', 0)}"
            + f", strand: \"{html_escape(strand)}\""
            + f", cdStart: {cdStart}"
            + f", cdEnd: {cdEnd}"
            + f", color: \"{html_escape(color)}\""
            + f", exons: [{json.dumps(exon_json)}]"
            + f", oChromStart: \"{html_escape(str(oChromStart))}\""
            + f", oChromEnd: \"{html_escape(str(oChromEnd))}\""
            + f", oStrand: \"{html_escape(oStrand)}\""
            + f", oChromSize: \"{html_escape(str(oChromSize))}\""
            + f", oChromStarts: \"{html_escape(oChromStarts)}\""
            + f", oSequence: \"{html_escape(oSequence)}\""
            + " }"
        )
    elif track_type == "selfpair":
        # Compute all fields in 0-based coordinates
        pstart = feature["ref_start"] - 1
        pend = feature["ref_end"]
        sstart = feature["cons_start"] - 1
        send = feature["cons_end"]
        start = min(pstart, sstart)
        end = max(pend, send)
        return ("{ chr: \"" + html_escape(seq_id) + "\""
                + f", start: {start}"
                + f", end: {end}"
                + f", pstart: {pstart}"
                + f", pend: {pend}"
                + f", sstart: {sstart}"
                + f", send: {send}"
                + f", color: \"{html_escape(feature.get('color','#e6194b'))}\""
                + f", strand: \"{html_escape(feature.get('orient','+'))}\""
                + " }")
    elif track_type == "chain":
        js = ("{ chr: \"" + html_escape(seq_id) + "\""
              + f", start: {feature['ref_start']-1}, end: {feature['ref_end']}"
              + f", name: \"{html_escape(feature.get('name',''))}\""
              + f", color: \"{html_escape(feature.get('color','gray'))}\""
              + f", strand: \"{html_escape(feature.get('orient','+'))}\""
              + (f", ostart: {feature['ostart']}" if 'ostart' in feature else "")
              + (f", oend: {feature['oend']}" if 'oend' in feature else "")
              + (f", osize: {feature['osize']}" if 'osize' in feature else ""))
        # Components
        if 'components' in feature and feature['components']:
            components = []
            for comp in feature['components']:
                comp0 = copy.deepcopy(comp)
                if 'start' in comp0:
                    comp0['start'] = comp0['start']-1
                if 'ostart' in comp0:
                    comp0['ostart'] = comp0['ostart']-1
                components.append(comp0)
            js += f", components: {json.dumps(components)}"
        js += " }"
        return js
    else:
        js = ("{ chr: \"" + html_escape(seq_id) + "\""
              + f", start: {feature['ref_start']-1}, end: {feature['ref_end']}"
              + f", name: \"{html_escape(feature.get('name',''))}\""
              + f", color: \"{html_escape(feature.get('color','gray'))}\""
              + (f", strand: \"{html_escape(feature.get('orient','+'))}\"" if 'orient' in feature else ""))
        # If this feature has 'components', emit them as JS with 0-based starts
        if 'components' in feature and feature['components']:
            # Adjust component start/ostart to 0-based
            components = []
            for comp in feature['components']:
                comp0 = copy.deepcopy(comp)
                if 'start' in comp0:
                    comp0['start'] = comp0['start'] - 1
                if 'ostart' in comp0:
                    comp0['ostart'] = comp0['ostart'] - 1
                components.append(comp0)
            js += f", components: {json.dumps(components)}"
        js += " }"
    return js

def cram_track_to_js(track):
    # Handles CRAM/CRAI ("alignment") track
    #    "format: 'cram'," +
    #    f"url: \"{html_escape(track['url'])}\"," +
    #    f"indexURL: \"{html_escape(track['indexURL'])}\"," +
    return (
        "{" +
        f"name: \"{html_escape(track.get('name','Seed Alignment'))}\"," +
        "type: 'alignment'," +
        "format: 'dfamsam'," +
        "url: 'seed.sam'," +
        "sourceType: 'dfamsam'," +
        f"fastaURL: \"{html_escape(track['fastaURL'])}\"," +
        "displayMode: 'SQUISHED',"
        "autoHeight: true" +
        "}"
    )

def build_track_js(tracks, seq_id):
    track_js_pieces = []
    for track in tracks:
        if track.get('type') == "alignment" and track.get('format') == "cram":
            # This is the CRAM track
            track_js_pieces.append(cram_track_to_js(track))
        else:
            features_js = ",\n".join(feature_to_js(f, seq_id, track.get('type')) for f in track['features'])
            js = f"""{{
                name: "{html_escape(track['name'])}",
                type: '{track.get("type", "annotation")}',
                displayMode: 'EXPANDED',
                autoHeight: true,
                features: [
{features_js}
                ]
            }}"""
            track_js_pieces.append(js)
    return ",\n".join(track_js_pieces)

def get_ref_fasta_url(tracks, default="ref.fa"):
    """
    If a CRAM track is present, use its fastaURL for reference.
    Otherwise, return the default.
    """
    for track in tracks:
        if track.get('type') == "alignment" and track.get('format') == "cram":
            return html_escape(track.get('fastaURL', default))
    return html_escape(default)

def generate_html(seq_id, title, tracks, output_html_path, igv_js_url="./igv.esm.min.js"):
    track_js = build_track_js(tracks, seq_id)
    ref_fasta_url = get_ref_fasta_url(tracks)
    html_content = HTML_TEMPLATE.substitute(
        title=html_escape(title),
        igv_js_url=igv_js_url,
        seq_id=html_escape(seq_id),
        track_js=track_js,
        ref_fasta_url=ref_fasta_url,
    )
    with open(output_html_path, "w") as out:
        out.write(html_content)
