#!/usr/bin/env python3
"""Build a browsable dashboard for combined-pipeline HTML outputs."""

from __future__ import annotations

import json
import re
import shutil
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = ROOT / "results"
EXCLUDED_TOP_LEVEL = {".git", "scripts", "motifs", "results"}
REGEX_DIRS = ("regex_hits", "isre_regex", "isre_hits_2kb")
SAN_REGEX_SUFFIX = "_Santhakumar_2018_isre_hits.tsv"
RUN_MARKERS = ("motif_map.html", "motif_map_uniq.html", "motif_map_all.html")

REPORT_ORDER = {
    "homerResults.html": 0,
    "knownResults.html": 1,
    "motif_map_uniq.html": 2,
    "motif_map.html": 3,
    "motif_map_heatmap_uniq_top5.html": 4,
    "motif_map_all.html": 5,
    "motif_map_heatmap_all_top5.html": 6,
    "motif_map_heatmap_top5.html": 7,
    "motif_map_heatmap_top10.html": 8,
    "geneOntology.html": 9,
}

REPORT_LABELS = {
    "motif_map.html": "Motif map",
    "motif_map_uniq.html": "Motif map unique",
    "motif_map_all.html": "Motif map all hits",
    "motif_map_heatmap_uniq_top5.html": "Unique heatmap top 5",
    "motif_map_heatmap_all_top5.html": "All hits heatmap top 5",
    "motif_map_heatmap_top5.html": "Heatmap top 5",
    "motif_map_heatmap_top10.html": "Heatmap top 10",
    "homerResults.html": "HOMER de novo motifs",
    "knownResults.html": "Known motif enrichment",
    "geneOntology.html": "Gene ontology",
}

TOKEN_LABELS = {
    "hs": "human",
    "ifn1": "IFN1",
    "orth": "orthologs",
    "ort": "orthologs",
    "2k": "2 kb",
    "2kb": "2 kb",
    "custombg": "custom bg",
    "builtinbg": "built-in bg",
    "cpg": "CpG",
    "upreg": "upregulated",
    "isg": "ISG",
    "gal": "chicken",
    "mus": "mouse",
    "uniq": "unique",
    "all": "all hits",
    "hits": "hits",
    "out": "out",
    "jaspar": "JASPAR",
    "paper": "paper",
    "compare": "compare",
}


def relpath(path: Path, start: Path) -> str:
    return path.relative_to(start).as_posix()


def copy_key(source_dir: Path) -> str:
    parts = list(source_dir.relative_to(ROOT).parts)
    species = parts[0]
    run_parts = [part.removeprefix("homer_") for part in parts[1:]]
    key = "_".join([species, *run_parts])
    return re.sub(r"[^A-Za-z0-9._-]+", "_", key).strip("_")


def titleize_slug(value: str) -> str:
    tokens = [token for token in re.split(r"[_\s-]+", value) if token]
    labels = [TOKEN_LABELS.get(token.lower(), token) for token in tokens]
    return " ".join(labels)


def run_priority(species: str, run: dict) -> tuple[int, str, str]:
    run_id = run["id"].lower()
    if species == "chicken":
        if "ifn1_builtinbg" in run_id:
            group = 0
        elif "builtinbg" in run_id:
            group = 1
        elif "custombg" in run_id:
            group = 2
        else:
            group = 3
    elif species == "human":
        if "builtinbg" in run_id:
            group = 0
        elif "custombg" in run_id:
            group = 1
        else:
            group = 2
    else:
        group = 0 if "builtinbg" in run_id else 1

    return (group, run["name"], run["id"])


def run_tags(run_name: str) -> list[str]:
    lower = run_name.lower()
    tags = []
    if "cpg" in lower:
        tags.append("CpG")
    if "orth" in lower or "ort" in lower:
        tags.append("orthologs")
    if "custombg" in lower:
        tags.append("custom bg")
    if "builtinbg" in lower:
        tags.append("built-in bg")
    if "ifn1" in lower:
        tags.append("IFN1")
    if "2kb" in lower or "2k" in lower:
        tags.append("2 kb")
    if "uniq" in lower:
        tags.append("unique")
    if "motif_map_all" in lower or lower.endswith("_all") or "/hits" in lower:
        tags.append("all hits")
    return tags


def report_label(report_path: Path, source_dir: Path) -> str:
    rel = report_path.relative_to(source_dir)
    filename = rel.name
    if filename in REPORT_LABELS:
        label = REPORT_LABELS[filename]
    else:
        stem = filename.removesuffix(".html")
        label = titleize_slug(stem)

    if len(rel.parts) > 1:
        prefix = titleize_slug("_".join(rel.parts[:-1]))
        if filename == "homerResults.html":
            label = "HOMER de novo motifs"
        return f"{prefix}: {label}"

    return label


def safe_filename(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("_")


def sort_reports(path: Path, source_dir: Path) -> tuple[int, int, str]:
    rel = path.relative_to(source_dir).as_posix()
    return (REPORT_ORDER.get(path.name, 99), len(path.relative_to(source_dir).parts), rel)


def is_report_html(path: Path, source_dir: Path) -> bool:
    rel = path.relative_to(source_dir)
    if len(rel.parts) == 1:
        return True
    if len(rel.parts) == 2:
        if rel.parts[0] in {"homerResults", "knownResults"}:
            return False
        return True
    return False


def discover_run_dirs() -> list[Path]:
    run_dirs = set()
    for top_level in ROOT.iterdir():
        if not top_level.is_dir() or top_level.name in EXCLUDED_TOP_LEVEL:
            continue
        for marker in RUN_MARKERS:
            for marker_path in top_level.rglob(marker):
                if any(part in EXCLUDED_TOP_LEVEL for part in marker_path.relative_to(ROOT).parts):
                    continue
                run_dirs.add(marker_path.parent)
    return sorted(run_dirs, key=lambda path: relpath(path, ROOT))


def discover_runs() -> list[dict]:
    runs = []
    for source_dir in discover_run_dirs():
        parts = source_dir.relative_to(ROOT).parts
        species = parts[0]
        if species in EXCLUDED_TOP_LEVEL:
            continue

        dest_name = copy_key(source_dir)
        html_reports = [
            report
            for report in source_dir.rglob("*.html")
            if is_report_html(report, source_dir)
        ]
        top_reports = sorted(html_reports, key=lambda path: sort_reports(path, source_dir))
        reports = [
            {
                "filename": report.relative_to(source_dir).as_posix(),
                "label": report_label(report, source_dir),
                "href": f"{dest_name}/{report.relative_to(source_dir).as_posix()}",
            }
            for report in top_reports
        ]

        source_tail = "_".join(part.removeprefix("homer_") for part in parts[1:])
        source_text = relpath(source_dir, ROOT)

        runs.append(
            {
                "id": dest_name,
                "species": species,
                "source_dir": source_text,
                "dest_dir": dest_name,
                "name": titleize_slug(source_tail),
                "raw_name": source_text,
                "tags": run_tags(source_text),
                "reports": reports,
            }
        )
    return sorted(runs, key=lambda run: (run["species"], run["dest_dir"]))


def parse_regex_summary(path: Path) -> dict | None:
    lines = [line.strip() for line in path.read_text(errors="replace").splitlines() if line.strip()]
    if not lines:
        return None

    last = lines[-1]
    fields = last.split("\t")
    ratio = fields[-1].strip()
    match = re.search(r"(\d+)\s*/\s*(\d+)", ratio)
    if not match:
        return None

    unique = int(match.group(1))
    total = int(match.group(2))
    percent = (unique / total * 100.0) if total else None
    source = path.name.removesuffix(SAN_REGEX_SUFFIX)
    if source == path.name:
        source = path.stem.replace("_Santhakumar_2018_isre_hits", "")
    if not source:
        source = path.parent.name

    dest_rel = Path("regex_hits") / safe_filename(relpath(path, ROOT).replace("/", "_"))

    return {
        "source": titleize_slug(source),
        "file": relpath(path, ROOT),
        "href": dest_rel.as_posix(),
        "unique": unique,
        "total": total,
        "ratio": f"{unique}/{total}",
        "percent": round(percent, 1) if percent is not None else None,
    }


def discover_regex_stats(species_names: set[str]) -> dict[str, list[dict]]:
    stats: dict[str, list[dict]] = {}
    for species in sorted(species_names):
        species_dir = ROOT / species
        if not species_dir.is_dir():
            continue

        files: list[Path] = []
        for regex_dir in REGEX_DIRS:
            candidate_dir = species_dir / regex_dir
            if candidate_dir.is_dir():
                files.extend(sorted(candidate_dir.glob(f"*{SAN_REGEX_SUFFIX}")))

        parsed = []
        seen = set()
        for path in files:
            if path in seen:
                continue
            seen.add(path)
            summary = parse_regex_summary(path)
            if summary is not None:
                parsed.append(summary)
        stats[species] = parsed
    return stats


def reset_results_dir() -> None:
    if RESULTS_DIR.exists():
        shutil.rmtree(RESULTS_DIR)
    RESULTS_DIR.mkdir()


def copy_run_reports(run: dict) -> None:
    source_dir = ROOT / run["source_dir"]
    dest_dir = RESULTS_DIR / run["dest_dir"]
    dest_dir.mkdir(parents=True, exist_ok=True)

    for report in run["reports"]:
        source_file = source_dir / report["filename"]
        dest_file = dest_dir / report["filename"]
        dest_file.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source_file, dest_file)


def copy_regex_summaries(regex_stats: dict[str, list[dict]]) -> None:
    for species_stats in regex_stats.values():
        for stat in species_stats:
            source_file = ROOT / stat["file"]
            dest_file = RESULTS_DIR / stat["href"]
            dest_file.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source_file, dest_file)


def build_species_payload(runs: list[dict], regex_stats: dict[str, list[dict]]) -> list[dict]:
    species_names = {run["species"] for run in runs} | set(regex_stats)
    species_payload = []
    for species in sorted(species_names):
        species_runs = [run for run in runs if run["species"] == species]
        species_runs.sort(key=lambda run: run_priority(species, run))
        species_payload.append(
            {
                "id": species,
                "name": species.capitalize(),
                "runs": species_runs,
                "regex": regex_stats.get(species, []),
            }
        )
    return species_payload


def build_html(species_payload: list[dict], total_reports: int) -> str:
    data_json = json.dumps(
        {
            "species": species_payload,
            "total_reports": total_reports,
        },
        indent=2,
        ensure_ascii=True,
    ).replace("</", "<\\/")

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>ISRE Analysis Results</title>
  <style>
    :root {{
      --bg: #f6f7f9;
      --panel: #ffffff;
      --text: #172026;
      --muted: #64707d;
      --border: #d9dee5;
      --accent: #0f766e;
      --accent-strong: #0b5f59;
      --accent-soft: #e3f4f1;
      --warm: #b45309;
      --warm-soft: #fff3df;
      --shadow: 0 1px 2px rgba(15, 23, 42, 0.08);
    }}

    * {{
      box-sizing: border-box;
    }}

    body {{
      margin: 0;
      min-height: 100vh;
      font-family: Arial, Helvetica, sans-serif;
      color: var(--text);
      background: var(--bg);
    }}

    .app {{
      display: grid;
      grid-template-columns: minmax(220px, 280px) minmax(260px, 360px) minmax(0, 1fr);
      min-height: 100vh;
    }}

    aside,
    .runs,
    main {{
      padding: 24px;
    }}

    aside {{
      border-right: 1px solid var(--border);
      background: #fbfcfd;
    }}

    .runs {{
      border-right: 1px solid var(--border);
      background: #ffffff;
    }}

    main {{
      overflow: auto;
    }}

    h1,
    h2,
    h3,
    p {{
      margin-top: 0;
    }}

    h1 {{
      font-size: 22px;
      line-height: 1.2;
      margin-bottom: 8px;
    }}

    h2 {{
      font-size: 20px;
      line-height: 1.25;
      margin-bottom: 14px;
    }}

    h3 {{
      font-size: 14px;
      margin-bottom: 10px;
      color: var(--muted);
      text-transform: uppercase;
      letter-spacing: 0;
    }}

    .summary {{
      color: var(--muted);
      font-size: 14px;
      line-height: 1.45;
      margin-bottom: 22px;
    }}

    .button-list {{
      display: grid;
      gap: 10px;
    }}

    button {{
      width: 100%;
      border: 1px solid var(--border);
      border-radius: 8px;
      background: var(--panel);
      color: var(--text);
      cursor: pointer;
      font: inherit;
      text-align: left;
      padding: 12px;
      box-shadow: var(--shadow);
    }}

    button:hover {{
      border-color: var(--accent);
    }}

    button.active {{
      border-color: var(--accent);
      background: var(--accent-soft);
    }}

    .button-title {{
      display: block;
      font-weight: 700;
      margin-bottom: 4px;
    }}

    .button-meta {{
      display: block;
      color: var(--muted);
      font-size: 13px;
      line-height: 1.35;
    }}

    .tags {{
      display: flex;
      flex-wrap: wrap;
      gap: 6px;
      margin-top: 8px;
    }}

    .tag {{
      display: inline-flex;
      align-items: center;
      min-height: 22px;
      padding: 3px 8px;
      border-radius: 999px;
      background: #edf2f7;
      color: #334155;
      font-size: 12px;
      line-height: 1;
    }}

    .stats {{
      display: grid;
      gap: 10px;
      margin-bottom: 24px;
    }}

    .stat-row,
    .report-link {{
      display: grid;
      gap: 6px;
      padding: 14px;
      border: 1px solid var(--border);
      border-radius: 8px;
      background: var(--panel);
      box-shadow: var(--shadow);
    }}

    .stat-row strong {{
      font-size: 18px;
    }}

    .stat-row a,
    .report-link a {{
      color: var(--accent-strong);
      font-weight: 700;
      text-decoration: none;
    }}

    .stat-row a:hover,
    .report-link a:hover {{
      text-decoration: underline;
    }}

    .stat-source,
    .report-source {{
      color: var(--muted);
      font-size: 13px;
      overflow-wrap: anywhere;
    }}

    .report-grid {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
      gap: 12px;
    }}

    .empty {{
      border: 1px dashed var(--border);
      border-radius: 8px;
      color: var(--muted);
      padding: 16px;
      background: rgba(255, 255, 255, 0.7);
    }}

    .source-line {{
      margin: 0 0 20px;
      color: var(--muted);
      font-size: 13px;
      overflow-wrap: anywhere;
    }}

    .run-heading {{
      display: flex;
      flex-wrap: wrap;
      gap: 8px;
      align-items: center;
      margin-bottom: 6px;
    }}

    .run-heading h2 {{
      margin: 0;
    }}

    .metric {{
      display: inline-flex;
      min-height: 28px;
      align-items: center;
      padding: 5px 10px;
      border-radius: 999px;
      background: var(--warm-soft);
      color: var(--warm);
      font-size: 13px;
      font-weight: 700;
    }}

    @media (max-width: 920px) {{
      .app {{
        grid-template-columns: 1fr;
      }}

      aside,
      .runs {{
        border-right: 0;
        border-bottom: 1px solid var(--border);
      }}
    }}
  </style>
</head>
<body>
  <div class="app">
    <aside>
      <h1>ISRE Analysis Results</h1>
      <p class="summary"><span id="speciesCount">0</span> species, <span id="reportCount">{total_reports}</span> HTML reports.</p>
      <h3>Species</h3>
      <div id="speciesList" class="button-list"></div>
    </aside>
    <section class="runs">
      <h3>Runs</h3>
      <div id="runList" class="button-list"></div>
    </section>
    <main>
      <section id="detail"></section>
    </main>
  </div>

  <script id="dashboard-data" type="application/json">{data_json}</script>
  <script>
    const DATA = JSON.parse(document.getElementById("dashboard-data").textContent);
    let activeSpeciesId = DATA.species[0] ? DATA.species[0].id : null;
    let activeRunId = null;

    const speciesList = document.getElementById("speciesList");
    const runList = document.getElementById("runList");
    const detail = document.getElementById("detail");
    document.getElementById("speciesCount").textContent = DATA.species.length;

    function escapeHtml(value) {{
      return String(value)
        .replaceAll("&", "&amp;")
        .replaceAll("<", "&lt;")
        .replaceAll(">", "&gt;")
        .replaceAll('"', "&quot;")
        .replaceAll("'", "&#039;");
    }}

    function currentSpecies() {{
      return DATA.species.find((species) => species.id === activeSpeciesId) || DATA.species[0];
    }}

    function currentRun(species) {{
      if (!species || !species.runs.length) {{
        return null;
      }}
      const existing = species.runs.find((run) => run.id === activeRunId);
      return existing || species.runs[0];
    }}

    function regexSummary(species) {{
      if (!species.regex.length) {{
        return "No regex summary";
      }}
      if (species.regex.length === 1) {{
        const stat = species.regex[0];
        return `${{stat.unique}}/${{stat.total}} unique genes`;
      }}
      return `${{species.regex.length}} regex summaries`;
    }}

    function renderSpeciesList() {{
      speciesList.innerHTML = "";
      DATA.species.forEach((species) => {{
        const button = document.createElement("button");
        button.className = species.id === activeSpeciesId ? "active" : "";
        button.type = "button";
        button.innerHTML = `
          <span class="button-title">${{escapeHtml(species.name)}}</span>
          <span class="button-meta">${{species.runs.length}} runs, ${{regexSummary(species)}}</span>
        `;
        button.addEventListener("click", () => {{
          activeSpeciesId = species.id;
          activeRunId = species.runs[0] ? species.runs[0].id : null;
          render();
        }});
        speciesList.appendChild(button);
      }});
    }}

    function renderRunList(species) {{
      runList.innerHTML = "";
      if (!species || !species.runs.length) {{
        runList.innerHTML = '<div class="empty">No HTML runs found for this species.</div>';
        return;
      }}

      species.runs.forEach((run) => {{
        const button = document.createElement("button");
        button.className = run.id === activeRunId ? "active" : "";
        button.type = "button";
        const tags = run.tags.map((tag) => `<span class="tag">${{escapeHtml(tag)}}</span>`).join("");
        button.innerHTML = `
          <span class="button-title">${{escapeHtml(run.name)}}</span>
          <span class="button-meta">${{run.reports.length}} reports</span>
          <span class="tags">${{tags}}</span>
        `;
        button.addEventListener("click", () => {{
          activeRunId = run.id;
          render();
        }});
        runList.appendChild(button);
      }});
    }}

    function renderRegexStats(species) {{
      if (!species.regex.length) {{
        return '<div class="empty">No Santhakumar_2018 regex summary TSV was found for this species.</div>';
      }}
      return `
        <div class="stats">
          ${{species.regex.map((stat) => `
            <div class="stat-row">
              <strong>${{stat.unique}} of ${{stat.total}} unique genes with hits</strong>
              <span>${{stat.percent === null ? "n/a" : stat.percent.toFixed(1) + "%"}} of genes in source set</span>
              <a href="${{escapeHtml(stat.href)}}">${{escapeHtml(stat.source)}}</a>
              <span class="stat-source">${{escapeHtml(stat.file)}}</span>
            </div>
          `).join("")}}
        </div>
      `;
    }}

    function renderReports(run) {{
      if (!run) {{
        return '<div class="empty">Select a species with HTML runs to show report links.</div>';
      }}
      return `
        <div class="run-heading">
          <h2>${{escapeHtml(run.name)}}</h2>
          <span class="metric">${{run.reports.length}} HTML reports</span>
        </div>
        <p class="source-line">Prepared from ${{escapeHtml(run.source_dir)}}</p>
        <div class="report-grid">
          ${{run.reports.map((report) => `
            <div class="report-link">
              <a href="${{escapeHtml(report.href)}}">${{escapeHtml(report.label)}}</a>
              <span class="report-source">${{escapeHtml(report.filename)}}</span>
            </div>
          `).join("")}}
        </div>
      `;
    }}

    function renderDetail(species, run) {{
      if (!species) {{
        detail.innerHTML = '<div class="empty">No results were found.</div>';
        return;
      }}
      detail.innerHTML = `
        <h2>${{escapeHtml(species.name)}}</h2>
        <h3>Santhakumar 2018 ISRE regex hits</h3>
        ${{renderRegexStats(species)}}
        <h3>Selected run</h3>
        ${{renderReports(run)}}
      `;
    }}

    function render() {{
      const species = currentSpecies();
      if (species && (!activeRunId || !species.runs.some((run) => run.id === activeRunId))) {{
        activeRunId = species.runs[0] ? species.runs[0].id : null;
      }}
      const run = currentRun(species);
      renderSpeciesList();
      renderRunList(species);
      renderDetail(species, run);
    }}

    render();
  </script>
</body>
</html>
"""


def main() -> None:
    runs = discover_runs()
    if not runs:
        raise SystemExit("No combined-pipeline run directories found.")

    species_names = {run["species"] for run in runs}
    top_level_species = {
        path.name
        for path in ROOT.iterdir()
        if path.is_dir() and path.name not in EXCLUDED_TOP_LEVEL
    }
    species_names |= top_level_species

    regex_stats = discover_regex_stats(species_names)
    species_payload = build_species_payload(runs, regex_stats)

    reset_results_dir()
    for run in runs:
        copy_run_reports(run)
    copy_regex_summaries(regex_stats)

    total_reports = sum(len(run["reports"]) for run in runs)
    (RESULTS_DIR / "index.html").write_text(
        build_html(species_payload, total_reports),
        encoding="utf-8",
    )

    print(f"Wrote {RESULTS_DIR / 'index.html'}")
    print(f"Copied {len(runs)} runs into {RESULTS_DIR}")
    print(f"Linked {total_reports} HTML reports")


if __name__ == "__main__":
    main()
