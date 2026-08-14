"""
app.py — Shiny for Python app for genostruct.

User uploads a genome FASTA, a GFF3 annotation, and a folder of PDB structure
models (multi-file upload). The protein FASTA is derived automatically from
the genome + GFF3 (see pipeline/genostruct.py::derive_pep_fasta) — pure
Python, no external tools. The pipeline then runs server-side and the
resulting self-contained viewer HTML is rendered inline (with a download
button as a convenience).

Run locally:
    shiny run app.py

Deploy (shinyapps.io):
    rsconnect deploy shiny . --name <account> --title genostruct-viewer
"""
import os
import sys
import shutil
import tempfile
import traceback

from shiny import App, reactive, render, ui

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "pipeline"))
import genostruct
import build_viewer

ASSETS = os.path.join(os.path.dirname(__file__), "assets")
VIEWER_SRC = os.path.join(os.path.dirname(__file__), "viewer_src")

PDB_EXTS = (".pdb",)

# Shiny's input_file() has no folder-picker option, so the underlying <input
# type=file> for pdb_files gets the (non-standard but broadly supported:
# Chrome/Edge/Safari/Firefox) webkitdirectory/directory attributes bolted on
# after render. The browser then uploads every file in the chosen folder;
# non-.pdb files are dropped server-side (see _pdb_only below).
_FOLDER_PICKER_JS = """
document.addEventListener('DOMContentLoaded', function() {
  var el = document.getElementById('pdb_files');
  if (el) {
    el.setAttribute('webkitdirectory', '');
    el.setAttribute('directory', '');
    el.setAttribute('mozdirectory', '');
  }
});
"""

app_ui = ui.page_fluid(
    ui.head_content(ui.tags.script(_FOLDER_PICKER_JS)),
    ui.panel_title("Genostruct viewer builder"),
    ui.p(
        "Upload a genome FASTA, a GFF3 annotation, and a folder of PDB "
        "structure models. The protein FASTA is derived automatically from "
        "the genome + GFF3 automatically — no need to supply one."
    ),
    ui.layout_sidebar(
        ui.sidebar(
            ui.input_file("genome_fa", "Genome FASTA (.fa/.fasta)", accept=[".fa", ".fasta", ".fna"]),
            ui.input_file("gff", "GFF3 annotation (.gff/.gff3)", accept=[".gff", ".gff3"]),
            ui.input_file("pdb_files", "PDB structures folder", multiple=True),
            ui.input_radio_buttons(
                "match_mode", "Match PDBs to transcripts by",
                {
                    "sequence": "Sequence (recommended) — searches each model's own "
                                "residues against every candidate protein; slower but "
                                "doesn't depend on filenames meaning anything",
                    "filename": "Filename — fast, but only works if filenames reliably "
                                 "encode the transcript id",
                },
                selected="sequence",
            ),
            ui.panel_conditional(
                "input.match_mode === 'filename'",
                ui.input_text(
                    "id_regex", "Transcript-ID regex",
                    placeholder=r"e.g. (.+)_ranked_0",
                ),
                ui.p(
                    "Regex whose first capture group extracts the transcript id from "
                    "each PDB filename (without extension).",
                    class_="text-muted", style="font-size:0.85em;",
                ),
            ),
            ui.input_action_button("run", "Build viewer", class_="btn-primary"),
            ui.download_button("download", "Download viewer HTML"),
            width=350,
        ),
        ui.output_ui("status"),
        ui.output_ui("viewer_frame"),
    ),
)


def _pdb_only(uploads):
    """Keep only *.pdb files from a folder upload (drops subfolder junk,
    .DS_Store, README files, etc. that ride along with a directory picker)."""
    return [f for f in uploads if f["name"].lower().endswith(PDB_EXTS)]


def build_viewer_html(genome_path, gff_path, pdb_dir, work_dir, id_regex):
    out_dir = os.path.join(work_dir, "out")
    genostruct.build(
        genome_fa=genome_path,
        gff3=gff_path,
        pdb_dir=pdb_dir,
        out_dir=out_dir,
        pep_fa=None,  # derived automatically (pure Python, no external tools)
        id_regex=id_regex or None,
    )
    html_out = os.path.join(work_dir, "genostruct_viewer.html")
    build_viewer.build_html(
        out_dir=out_dir,
        molstar_js=os.path.join(ASSETS, "molstar.js"),
        molstar_css=os.path.join(ASSETS, "molstar.css"),
        pako_js=os.path.join(ASSETS, "pako_inflate.min.js"),
        template=os.path.join(VIEWER_SRC, "template.html"),
        app_js=os.path.join(VIEWER_SRC, "app.js"),
        app_css=os.path.join(VIEWER_SRC, "app.css"),
        html_out=html_out,
    )
    return html_out


def server(input, output, session):
    result_html = reactive.value(None)
    error_msg = reactive.value(None)

    @reactive.effect
    @reactive.event(input.run)
    def _run():
        result_html.set(None)
        error_msg.set(None)

        genome_file = input.genome_fa()
        gff_file = input.gff()
        pdb_uploads = _pdb_only(input.pdb_files() or [])
        if not genome_file or not gff_file or not pdb_uploads:
            error_msg.set(
                "Please upload a genome FASTA, a GFF3 file, and a folder "
                "containing at least one .pdb structure."
            )
            return

        id_regex = None
        if input.match_mode() == "filename":
            id_regex = (input.id_regex() or "").strip()
            if not id_regex:
                error_msg.set(
                    "Filename matching is selected but no transcript-ID regex was "
                    "given. Enter one, or switch to sequence matching."
                )
                return

        work_dir = tempfile.mkdtemp(prefix="genostruct_")
        try:
            pdb_dir = os.path.join(work_dir, "pdbs")
            os.makedirs(pdb_dir, exist_ok=True)
            for f in pdb_uploads:
                # folder uploads report name as a relative path (e.g. "pdbs/x.pdb");
                # flatten to the basename since we collect everything into pdb_dir.
                shutil.copy(f["datapath"], os.path.join(pdb_dir, os.path.basename(f["name"])))

            html_out = build_viewer_html(
                genome_file[0]["datapath"],
                gff_file[0]["datapath"],
                pdb_dir,
                work_dir,
                id_regex,
            )
            result_html.set(html_out)
        except Exception as ex:
            error_msg.set(f"{type(ex).__name__}: {ex}\n\n{traceback.format_exc()}")
        finally:
            pass  # work_dir cleanup left to OS temp GC; html is read into memory for rendering/download

    @render.ui
    def status():
        if error_msg():
            return ui.markdown(f"**Error:**\n```\n{error_msg()}\n```")
        if result_html():
            return ui.tags.p("Viewer built successfully.", class_="text-success")
        return ui.tags.p("Upload inputs and click 'Build viewer'.")

    @render.ui
    def viewer_frame():
        path = result_html()
        if not path:
            return None
        with open(path, encoding="utf-8") as fh:
            html = fh.read()
        return ui.tags.iframe(
            srcdoc=html,
            style="width:100%; height:85vh; border:1px solid #ccc;",
        )

    @render.download(filename="genostruct_viewer.html")
    def download():
        path = result_html()
        if not path:
            return
        with open(path, "rb") as fh:
            yield fh.read()


app = App(app_ui, server)
