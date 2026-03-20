#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(biomaRt)
  library(httr)
  library(jsonlite)
  library(dplyr)
  library(readr)
  library(rtracklayer)
  library(stringr)
  library(tidyr)
  library(purrr)
})

default_go_roots <- c(
  "GO:0140888", # interferon-mediated signaling pathway
  "GO:0034340", # response to type I interferon
  "GO:0034341", # response to type II interferon
  "GO:0034342", # response to type III interferon
  "GO:0071357", # cellular response to type I interferon
  "GO:0071346", # cellular response to type II interferon
  "GO:0071358", # cellular response to type III interferon
  "GO:0032479", # regulation of type I interferon production
  "GO:0032649", # regulation of type II interferon production
  "GO:0034344"  # regulation of type III interferon production
)

# =========================
# HELPERS
# =========================
default_config <- function() {
  list(
    organism_dir = "/home/pospim/Desktop/isres/chicken/raw",
    gtf_file = NULL,
    gene2go_file = NULL,
    gene2ensembl_file = NULL,
    go_obo_file = "/home/pospim/Desktop/isres/scripts/misc/go-basic.obo",
    biomart_dataset = "ggallus_gene_ensembl",
    output_prefix = "chicken_homer_background",
    target_gene_file = NULL,
    go_root_file = NULL,
    discovered_fasta_files = character(),
    exclude_tables = list(
      list(
        path = "/home/pospim/Desktop/isres/chicken/raw/hs_ifn1_isg_galgal_orth.tsv",
        reason = "uploaded_human_isg_chicken_ortholog"
      ),
      list(
        path = "/home/pospim/Desktop/isres/chicken/raw/chicken_isgs_mapped.csv",
        reason = "uploaded_skinner_chicken_isg"
      )
    )
  )
}

normalize_id <- function(x) {
  x %>% as.character() %>% str_trim() %>% na_if("") %>% toupper()
}

usage_text <- function() {
  paste(
    "Usage:",
    "  Rscript scripts/build_chicken_homer_background.R [options]",
    "",
    "Primary workflow:",
    "  --organism-dir=DIR",
    "",
    "Required for non-default species:",
    "  --organism-dir=DIR or --gtf=PATH",
    "",
    "Optional:",
    "  --organism-dir=DIR",
    "  --output-prefix=PREFIX",
    "  --biomart-dataset=DATASET",
    "  --gene2go=PATH",
    "  --gene2ensembl=PATH",
    "  --go-obo=PATH",
    "  --go-root-file=PATH",
    "  --target-gene-file=PATH",
    "  '--exclude-table=PATH|REASON'   (repeatable; quote because of |)",
    "  --help",
    "",
    "Example:",
    "  Rscript scripts/build_chicken_homer_background.R \\",
    "    --organism-dir=mouse/raw \\",
    "    --output-prefix=mouse_homer_background \\",
    "    '--exclude-table=mouse/raw/ifn1_upreg.tsv|ifn1_upregulated'",
    sep = "\n"
  )
}

parse_args <- function(args) {
  cfg <- default_config()
  using_default_exclude_tables <- TRUE

  if (length(args) == 0) {
    return(cfg)
  }

  for (arg in args) {
    if (arg %in% c("--help", "-h")) {
      cat(usage_text(), "\n")
      quit(save = "no", status = 0)
    }

    if (!startsWith(arg, "--")) {
      stop("Unrecognized argument: ", arg, "\n\n", usage_text())
    }

    pieces <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    key <- pieces[1]
    value <- if (length(pieces) > 1) paste(pieces[-1], collapse = "=") else NULL

    if (is.null(value) || identical(value, "")) {
      stop("Argument requires a value: --", key)
    }

    if (key == "gtf") {
      cfg$gtf_file <- value
    } else if (key %in% c("organism-dir", "annotation-dir")) {
      cfg$organism_dir <- value
    } else if (key == "output-prefix") {
      cfg$output_prefix <- value
    } else if (key == "biomart-dataset") {
      cfg$biomart_dataset <- value
    } else if (key == "gene2go") {
      cfg$gene2go_file <- value
    } else if (key == "gene2ensembl") {
      cfg$gene2ensembl_file <- value
    } else if (key == "go-obo") {
      cfg$go_obo_file <- value
    } else if (key == "go-root-file") {
      cfg$go_root_file <- value
    } else if (key == "target-gene-file") {
      cfg$target_gene_file <- value
    } else if (key == "exclude-table") {
      if (using_default_exclude_tables) {
        cfg$exclude_tables <- list()
        using_default_exclude_tables <- FALSE
      }
      table_parts <- strsplit(value, "|", fixed = TRUE)[[1]]
      if (length(table_parts) != 2) {
        stop("--exclude-table must be in the form PATH|REASON")
      }
      cfg$exclude_tables[[length(cfg$exclude_tables) + 1]] <- list(
        path = table_parts[1],
        reason = table_parts[2]
      )
    } else {
      stop("Unknown option: --", key)
    }
  }

  cfg
}

discover_first_match <- function(files, pattern, label, required = FALSE) {
  matches <- files[str_detect(basename(files), regex(pattern, ignore_case = TRUE))]
  matches <- sort(matches)
  if (length(matches) == 0) {
    if (required) {
      stop("Could not auto-discover ", label, " from organism directory")
    }
    return(NULL)
  }
  if (length(matches) > 1) {
    warning("Multiple candidates found for ", label, ". Using: ", matches[[1]])
  }
  matches[[1]]
}

is_nonempty_file <- function(path) {
  !is.null(path) && nzchar(path) && file.exists(path) && isTRUE(file.info(path)$size > 0)
}

discover_species_files <- function(organism_dir) {
  files <- list.files(organism_dir, recursive = TRUE, full.names = TRUE, all.files = FALSE, no.. = TRUE)
  files <- files[file.info(files)$isdir %in% FALSE]

  list(
    gtf_file = discover_first_match(files, "\\.gtf(\\.gz)?$", "GTF"),
    gene2go_file = discover_first_match(files, "gene2go.*\\.tsv(\\.gz)?$", "gene2go"),
    gene2ensembl_file = discover_first_match(files, "gene2ensembl.*\\.tsv(\\.gz)?$", "gene2ensembl"),
    go_root_file = discover_first_match(files, "go.*exclude.*\\.(txt|tsv|csv)$", "GO root file"),
    fasta_files = sort(files[str_detect(basename(files), regex("\\.(fa|fasta)(\\.gz)?$", ignore_case = TRUE))])
  )
}

resolve_config <- function(cfg) {
  if (!is.null(cfg$organism_dir)) {
    if (!dir.exists(cfg$organism_dir)) {
      stop("Organism directory not found: ", cfg$organism_dir)
    }

    discovered <- discover_species_files(cfg$organism_dir)
    if (is.null(cfg$gtf_file)) {
      cfg$gtf_file <- discovered$gtf_file
    }
    if (is.null(cfg$gene2go_file)) {
      cfg$gene2go_file <- discovered$gene2go_file
    }
    if (is.null(cfg$gene2ensembl_file)) {
      cfg$gene2ensembl_file <- discovered$gene2ensembl_file
    }
    if (is.null(cfg$go_root_file)) {
      cfg$go_root_file <- discovered$go_root_file
    }
    cfg$discovered_fasta_files <- discovered$fasta_files
  }

  cfg
}

validate_args <- function(cfg) {
  if (is.null(cfg$gtf_file) || !nzchar(cfg$gtf_file)) {
    stop("No GTF file was provided or auto-discovered. Use --organism-dir or --gtf.")
  }

  required_files <- c(cfg$gtf_file)
  missing_required <- required_files[!file.exists(required_files)]
  if (length(missing_required) > 0) {
    stop("Required input file(s) missing: ", paste(missing_required, collapse = ", "))
  }

  for (tbl in cfg$exclude_tables) {
    if (!file.exists(tbl$path)) {
      stop("Exclude table not found: ", tbl$path)
    }
  }

  if (!is.null(cfg$target_gene_file) && !file.exists(cfg$target_gene_file)) {
    stop("Target gene file not found: ", cfg$target_gene_file)
  }

  if (!is.null(cfg$go_root_file) && !file.exists(cfg$go_root_file)) {
    stop("GO root file not found: ", cfg$go_root_file)
  }
}

slugify <- function(x) {
  x %>%
    tolower() %>%
    str_replace_all("[^a-z0-9]+", "_") %>%
    str_replace_all("^_+|_+$", "")
}

read_gene_table <- function(path) {
  ext <- tools::file_ext(path)
  if (tolower(ext) %in% c("tsv", "txt")) {
    suppressWarnings(readr::read_tsv(path, show_col_types = FALSE, progress = FALSE))
  } else {
    suppressWarnings(readr::read_csv(path, show_col_types = FALSE, progress = FALSE))
  }
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

extract_gtf_attr <- function(x, key) {
  out <- stringr::str_match(x %||% "", paste0(key, ' "([^"]+)"'))[, 2]
  out[is.na(out)] <- ""
  out
}

build_gene_table_from_gtf <- function(gtf_path) {
  gtf <- rtracklayer::import(gtf_path)
  genes <- gtf[gtf$type == "gene"]
  if (length(genes) == 0) {
    stop("No gene features were found in GTF: ", gtf_path)
  }

  attrs <- mcols(genes)
  gene_biotype <- attrs$gene_biotype %||% attrs$gene_type %||% extract_gtf_attr(attrs$gene_id, "gene_biotype")
  gene_name <- attrs$gene_name %||% extract_gtf_attr(attrs$gene_name, "gene_name")
  gene_id <- attrs$gene_id %||% extract_gtf_attr(attrs$gene_id, "gene_id")

  tibble(
    ensembl_gene_id = sub("\\..*$", "", as.character(gene_id)),
    external_gene_name = na_if(as.character(gene_name), ""),
    gene_biotype = as.character(gene_biotype),
    chromosome_name = as.character(seqnames(genes)),
    start_position = start(genes),
    end_position = end(genes),
    strand = case_when(
      as.character(strand(genes)) == "+" ~ 1L,
      as.character(strand(genes)) == "-" ~ -1L,
      TRUE ~ 0L
    )
  ) %>%
    filter(gene_biotype == "protein_coding", !is.na(ensembl_gene_id), ensembl_gene_id != "") %>%
    distinct()
}

extract_ids_from_uploaded_file <- function(df) {
  cols <- names(df)
  ensembl_cols <- cols[str_detect(tolower(cols), "ensmbl|ensembl")]
  symbol_cols  <- cols[str_detect(tolower(cols), "gene_id|gene name|gene_name|symbol")]

  ensembl_vals <- if (length(ensembl_cols) > 0) unlist(df[ensembl_cols], use.names = FALSE) else character(0)
  symbol_vals  <- if (length(symbol_cols)  > 0) unlist(df[symbol_cols],  use.names = FALSE) else character(0)

  tibble(
    input_value = c(as.character(ensembl_vals), as.character(symbol_vals)),
    input_type  = c(rep("ensembl_candidate", length(ensembl_vals)), rep("symbol_candidate", length(symbol_vals)))
  ) %>%
    mutate(input_value = str_trim(input_value)) %>%
    filter(!is.na(input_value), input_value != "") %>%
    distinct()
}

map_ids_local <- function(values, gene_table) {
  vals <- tibble(raw = unique(values)) %>% mutate(raw_norm = normalize_id(raw))
  if (nrow(vals) == 0) {
    return(tibble(input = character(), ensembl_gene_id = character(), external_gene_name = character(), mapping_source = character()))
  }

  gene_id_map <- gene_table %>%
    transmute(
      input = normalize_id(ensembl_gene_id),
      ensembl_gene_id,
      external_gene_name,
      mapping_source = "ensembl_gene_id"
    ) %>%
    distinct()

  symbol_map <- gene_table %>%
    filter(!is.na(external_gene_name), external_gene_name != "") %>%
    transmute(
      input = normalize_id(external_gene_name),
      ensembl_gene_id,
      external_gene_name,
      mapping_source = "external_gene_name"
    ) %>%
    distinct()

  matched_by_id <- vals %>%
    inner_join(gene_id_map, by = c("raw_norm" = "input")) %>%
    transmute(input = raw_norm, ensembl_gene_id, external_gene_name, mapping_source)

  unmatched_vals <- vals %>%
    anti_join(matched_by_id %>% distinct(input), by = c("raw_norm" = "input"))

  matched_by_symbol <- unmatched_vals %>%
    inner_join(symbol_map, by = c("raw_norm" = "input")) %>%
    transmute(input = raw_norm, ensembl_gene_id, external_gene_name, mapping_source)

  bind_rows(matched_by_id, matched_by_symbol) %>% distinct()
}

read_go_roots <- function(path = NULL) {
  if (is.null(path)) {
    return(default_go_roots)
  }

  lines <- readLines(path, warn = FALSE)
  roots <- stringr::str_extract(lines, "GO:[0-9]{7}")
  roots <- unique(roots[!is.na(roots)])
  if (length(roots) == 0) {
    stop("No GO IDs were found in go root file: ", path)
  }
  roots
}

safe_use_ensembl <- function(dataset) {
  mirrors <- c(NA, "useast", "asia")
  for (mirror in mirrors) {
    mart <- tryCatch(
      {
        if (is.na(mirror)) {
          useEnsembl(biomart = "genes", dataset = dataset)
        } else {
          useEnsembl(biomart = "genes", dataset = dataset, mirror = mirror)
        }
      },
      error = function(e) NULL
    )
    if (!is.null(mart)) {
      return(mart)
    }
  }
  NULL
}

safe_getBM <- function(attributes, mart, filters = NULL, values = NULL) {
  tryCatch(
    {
      getBM(attributes = attributes, filters = filters, values = values, mart = mart)
    },
    error = function(e) {
      warning("BioMart query failed: ", conditionMessage(e))
      NULL
    }
  )
}

parse_go_obo <- function(path) {
  lines <- readLines(path, warn = FALSE)
  records <- list()
  current <- NULL

  flush_current <- function(rec, out) {
    if (is.null(rec)) {
      return(out)
    }
    if (!identical(rec$namespace, "biological_process") || isTRUE(rec$is_obsolete) || is.null(rec$id)) {
      return(out)
    }
    out[[length(out) + 1]] <- tibble(
      go_id = rec$id,
      parent_id = rec$parents %||% character(0)
    )
    out
  }

  for (line in lines) {
    if (line == "[Term]") {
      records <- flush_current(current, records)
      current <- list(id = NULL, namespace = NULL, is_obsolete = FALSE, parents = character(0))
      next
    }
    if (is.null(current) || line == "" || startsWith(line, "[")) {
      next
    }
    if (startsWith(line, "id: GO:")) {
      current$id <- sub("^id: ", "", line)
    } else if (startsWith(line, "namespace: ")) {
      current$namespace <- sub("^namespace: ", "", line)
    } else if (startsWith(line, "is_a: GO:")) {
      current$parents <- c(current$parents, sub("^is_a: (GO:[0-9]+).*$", "\\1", line))
    } else if (startsWith(line, "relationship: part_of GO:")) {
      current$parents <- c(current$parents, sub("^relationship: part_of (GO:[0-9]+).*$", "\\1", line))
    } else if (startsWith(line, "is_obsolete: true")) {
      current$is_obsolete <- TRUE
    }
  }
  records <- flush_current(current, records)

  bind_rows(records) %>% distinct()
}

expand_go_descendants_local <- function(go_ids, obo_path) {
  go_edges <- parse_go_obo(obo_path)
  if (nrow(go_edges) == 0) {
    return(unique(go_ids))
  }

  children_by_parent <- split(go_edges$go_id, go_edges$parent_id)
  seen <- unique(go_ids)
  queue <- seen

  while (length(queue) > 0) {
    current <- queue[[1]]
    queue <- queue[-1]
    kids <- unique(children_by_parent[[current]] %||% character(0))
    new_kids <- setdiff(kids, seen)
    if (length(new_kids) > 0) {
      seen <- c(seen, new_kids)
      queue <- c(queue, new_kids)
    }
  }

  unique(seen)
}

build_go_excluded_local <- function(gene2go_path, gene2ensembl_path, obo_path, go_roots, gene_table) {
  if (!is_nonempty_file(gene2go_path) || !is_nonempty_file(gene2ensembl_path) || !is_nonempty_file(obo_path)) {
    return(NULL)
  }

  exclude_go_all <- expand_go_descendants_local(go_roots, obo_path)

  gene2go <- readr::read_tsv(
    gene2go_path,
    skip = 1,
    col_names = c("tax_id", "GeneID", "GO_ID", "Evidence", "Qualifier", "GO_term", "PubMed", "Category"),
    show_col_types = FALSE,
    progress = FALSE
  ) %>%
    transmute(
      GeneID = as.character(GeneID),
      go_id = GO_ID,
      category = Category
    ) %>%
    filter(category == "Process", go_id %in% exclude_go_all) %>%
    distinct()

  gene2ensembl <- readr::read_tsv(
    gene2ensembl_path,
    skip = 1,
    col_names = c(
      "tax_id",
      "GeneID",
      "Ensembl_gene_identifier",
      "RNA_nucleotide_accession.version",
      "Ensembl_rna_identifier",
      "protein_accession.version",
      "Ensembl_protein_identifier"
    ),
    show_col_types = FALSE,
    progress = FALSE
  ) %>%
    transmute(
      GeneID = as.character(GeneID),
      ensembl_gene_id = sub("\\..*$", "", Ensembl_gene_identifier)
    ) %>%
    filter(!is.na(ensembl_gene_id), ensembl_gene_id != "-") %>%
    distinct()

  gene2go %>%
    inner_join(gene2ensembl, by = "GeneID") %>%
    inner_join(gene_table %>% select(ensembl_gene_id, external_gene_name), by = "ensembl_gene_id") %>%
    select(ensembl_gene_id, external_gene_name, go_id) %>%
    distinct() %>%
    mutate(exclusion_reason = "ifn_related_go")
}

expand_quickgo_descendants <- function(go_ids) {
  base_url <- "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms"
  headers <- add_headers(Accept = "application/json")

  get_descendants_one <- function(go_id) {
    url <- paste0(base_url, "/", URLencode(go_id, reserved = TRUE), "/descendants")
    resp <- GET(url, headers)
    stop_for_status(resp)
    txt <- content(resp, as = "text", encoding = "UTF-8")
    js <- fromJSON(txt, simplifyDataFrame = TRUE)

    ids <- character(0)
    if (!is.null(js$results) && nrow(js$results) > 0) {
      if ("descendants" %in% names(js$results) && length(js$results$descendants[[1]]) > 0) {
        ids <- js$results$descendants[[1]]
      }
      if ("children" %in% names(js$results) && length(js$results$children[[1]]) > 0) {
        child_ids <- js$results$children[[1]]$id
        ids <- unique(c(ids, child_ids))
      }
    }
    unique(c(go_id, ids))
  }

  unique(unlist(lapply(go_ids, get_descendants_one)))
}

config <- parse_args(commandArgs(trailingOnly = TRUE))
config <- resolve_config(config)
validate_args(config)
exclude_go_roots <- read_go_roots(config$go_root_file)

message("Using inputs:")
message("  organism_dir: ", config$organism_dir %||% "<none>")
message("  gtf: ", config$gtf_file)
message("  gene2go: ", config$gene2go_file %||% "<none>")
message("  gene2ensembl: ", config$gene2ensembl_file %||% "<none>")
message("  go_obo: ", config$go_obo_file %||% "<none>")
message("  go_root_file: ", config$go_root_file %||% "<built-in defaults>")

# =========================
# BUILD GENE UNIVERSE
# =========================
gene_table <- build_gene_table_from_gtf(config$gtf_file)

# =========================
# OPTIONAL ENSEMBL CONNECTION FOR GO ANNOTATIONS
# =========================
mart <- if (!is.null(config$biomart_dataset) && nzchar(config$biomart_dataset)) {
  safe_use_ensembl(config$biomart_dataset)
} else {
  NULL
}

# =========================
# READ EXCLUSION TABLES
# =========================
exclude_maps <- list()

for (i in seq_along(config$exclude_tables)) {
  tbl_cfg <- config$exclude_tables[[i]]
  tbl_df <- read_gene_table(tbl_cfg$path)
  tbl_ids <- extract_ids_from_uploaded_file(tbl_df)
  tbl_map <- map_ids_local(tbl_ids$input_value, gene_table) %>%
    mutate(
      source_file = tbl_cfg$path,
      exclusion_reason = tbl_cfg$reason
    )
  exclude_maps[[i]] <- tbl_map
}

exclude_maps_combined <- if (length(exclude_maps) > 0) {
  bind_rows(exclude_maps)
} else {
  tibble(
    input = character(),
    ensembl_gene_id = character(),
    external_gene_name = character(),
    mapping_source = character(),
    source_file = character(),
    exclusion_reason = character()
  )
}

# =========================
# OPTIONAL TARGET FILE
# =========================
if (!is.null(config$target_gene_file) && file.exists(config$target_gene_file)) {
  target_vals <- read_lines(config$target_gene_file, progress = FALSE) %>% str_trim() %>% discard(~ .x == "")
  target_map <- map_ids_local(target_vals, gene_table) %>%
    mutate(exclusion_reason = "target_gene")
} else {
  target_map <- tibble(
    input = character(),
    ensembl_gene_id = character(),
    external_gene_name = character(),
    mapping_source = character(),
    exclusion_reason = character()
  )
}

# =========================
# GO-BASED EXCLUSIONS
# =========================
go_excluded <- tibble(
  ensembl_gene_id = character(),
  external_gene_name = character(),
  go_id = character(),
  exclusion_reason = character()
)

go_excluded_local <- build_go_excluded_local(
  gene2go_path = config$gene2go_file,
  gene2ensembl_path = config$gene2ensembl_file,
  obo_path = config$go_obo_file,
  go_roots = exclude_go_roots,
  gene_table = gene_table
)

if (!is.null(go_excluded_local)) {
  go_excluded <- go_excluded_local
} else if (!is.null(mart)) {
  exclude_go_all <- tryCatch(
    expand_quickgo_descendants(exclude_go_roots),
    error = function(e) {
      warning("QuickGO lookup failed: ", conditionMessage(e))
      exclude_go_roots
    }
  )

  universe_go <- safe_getBM(
    attributes = c("ensembl_gene_id", "external_gene_name", "go_id", "namespace_1003"),
    mart = mart
  )

  if (!is.null(universe_go)) {
    go_excluded <- universe_go %>%
      filter(namespace_1003 == "biological_process", go_id %in% exclude_go_all) %>%
      select(ensembl_gene_id, external_gene_name, go_id) %>%
      distinct() %>%
      mutate(exclusion_reason = "ifn_related_go")
  } else {
    warning("Skipping GO-based exclusions because the BioMart GO annotation query failed.")
  }
} else {
  warning("Skipping GO-based exclusions because no local GO cache files were found and BioMart is unavailable.")
}

# =========================
# COMBINE EXCLUSIONS
# =========================
all_excluded <- bind_rows(
  exclude_maps_combined %>% select(ensembl_gene_id, external_gene_name, exclusion_reason),
  target_map %>% select(ensembl_gene_id, external_gene_name, exclusion_reason),
  go_excluded %>% select(ensembl_gene_id, external_gene_name, exclusion_reason)
) %>%
  filter(!is.na(ensembl_gene_id), ensembl_gene_id != "") %>%
  distinct()

excluded_summary <- all_excluded %>%
  group_by(ensembl_gene_id, external_gene_name) %>%
  summarise(exclusion_reason = paste(sort(unique(exclusion_reason)), collapse = ";"), .groups = "drop")

background_candidates <- gene_table %>%
  anti_join(excluded_summary, by = c("ensembl_gene_id", "external_gene_name")) %>%
  distinct()

# =========================
# WRITE OUTPUTS
# =========================
excluded_summary_path <- paste0(config$output_prefix, "_excluded_genes_with_reasons.csv")
background_path <- paste0(config$output_prefix, "_background_candidate_genes.csv")
exclude_maps_path <- paste0(config$output_prefix, "_exclude_table_mappings.csv")
go_excluded_path <- paste0(config$output_prefix, "_go_excluded_genes.csv")

write_csv(excluded_summary, excluded_summary_path)
write_csv(background_candidates, background_path)
write_csv(exclude_maps_combined, exclude_maps_path)
write_csv(go_excluded, go_excluded_path)

if (length(exclude_maps) > 0) {
  for (i in seq_along(exclude_maps)) {
    tbl_cfg <- config$exclude_tables[[i]]
    tbl_path <- paste0(
      config$output_prefix,
      "_exclude_table_",
      sprintf("%02d", i),
      "_",
      slugify(tbl_cfg$reason),
      "_mapped.csv"
    )
    write_csv(exclude_maps[[i]], tbl_path)
  }
}

message("Done. Wrote:")
message("  ", excluded_summary_path)
message("  ", background_path)
message("  ", exclude_maps_path)
message("  ", go_excluded_path)
if (length(exclude_maps) > 0) {
  for (i in seq_along(exclude_maps)) {
    message(
      "  ",
      paste0(
        config$output_prefix,
        "_exclude_table_",
        sprintf("%02d", i),
        "_",
        slugify(config$exclude_tables[[i]]$reason),
        "_mapped.csv"
      )
    )
  }
}
if (is.null(config$target_gene_file)) {
  message("Target exclusion was skipped because no --target-gene-file was provided")
}
