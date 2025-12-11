#!/usr/bin/env Rscript

################################################################################
# Ibex Variant Classification Pipeline
#################################################################################
# This script:
#   1. Reads STAR ReadsPerGene.out.tab files and builds a gene x sample
#      count matrix.
#   2. Computes CPM using edgeR and flags "expressed" genes:
#        CPM >= cpm_threshold in >= min_samples_at_cutoff samples.
#   3. Parses a GTF file to obtain CDS (or gene) genomic intervals.
#   4. Filters intervals to those belonging to expressed genes.
#   5. For each *.pos file in positions_dir, outputs SNPs that fall in
#      expressed genes (CDS by default) as both TSV and BED.
#
# Example:
#   Rscript expressed_gene_SNP_overlap.R \
#       --star_root       path/to/star \
#       --gtf             CapIbe.FinalAnnott.gtf \
#       --column          4 \
#       --cpm_threshold   0.1 \
#       --min_samples     1 \
#       --region_type     CDS \
#       --positions_dir   path/to/position/files \
#       --pos_pattern     '\\.pos$' \
#       --out_dir         results_new_positions
##
# Author: Gözde Cilingir
# Last Updated: December 2025
################################################################################
# Environment Setup
################################################################################

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(edgeR)
  library(stringr)
})
################################################################################
# COMMAND-LINE OPTIONS
################################################################################
option_list <- list(
  make_option(
    c("--star_root"),
    type    = "character",
    default = ".",
    help    = "Root directory containing STAR ReadsPerGene.out.tab files (searched recursively) [default %default]"
  ),
  make_option(
    c("--gtf"),
    type    = "character",
    default = "CapIbe.FinalAnnott.gtf",
    help    = "GTF annotation file [default %default]"
  ),
  make_option(
    c("--column"),
    type    = "integer",
    default = 4,
    help    = "STAR column to use: 2 = unstranded, 3 = stranded sense, 4 = reverse-stranded [default %default]"
  ),
  make_option(
    c("--cpm_threshold"),
    type    = "double",
    default = 0.1,
    help    = "CPM threshold for calling a gene expressed [default %default]"
  ),
  make_option(
    c("--min_samples"),
    type    = "integer",
    default = 1,
    help    = "Minimum number of samples at or above CPM threshold [default %default]"
  ),
  make_option(
    c("--region_type"),
    type    = "character",
    default = "CDS",
    help    = "Region type to overlap with: 'CDS' or 'gene' [default %default]"
  ),
  make_option(
    c("--positions_dir"),
    type    = "character",
    default = "new_positions",
    help    = "Directory containing position files [default %default]"
  ),
  make_option(
    c("--pos_pattern"),
    type    = "character",
    default = "\\.pos$",
    help    = "Regex pattern for position files within positions_dir [default %default]"
  ),
  make_option(
    c("--out_dir"),
    type    = "character",
    default = "results_new_positions",
    help    = "Output directory [default %default]"
  )
)

opt_parser <- OptionParser(
  option_list = option_list,
  description = "Annotate SNP positions with expressed genes and write TSV/BED outputs."
)
opt <- parse_args(opt_parser)

# Map options to internal variables
star_root            <- opt$star_root
gtf_file             <- opt$gtf
use_star_column      <- opt$column
cpm_threshold        <- opt$cpm_threshold
min_samples_at_cutoff <- opt$min_samples
region_type          <- opt$region_type
positions_dir        <- opt$positions_dir
pos_pattern          <- opt$pos_pattern
out_dir              <- opt$out_dir

region_type <- toupper(region_type)

message("= PARAMETERS =")
message("  star_root        : ", star_root)
message("  gtf_file         : ", gtf_file)
message("  use_star_column  : ", use_star_column)
message("  cpm_threshold    : ", cpm_threshold)
message("  min_samples      : ", min_samples_at_cutoff)
message("  region_type      : ", region_type)
message("  positions_dir    : ", positions_dir)
message("  pos_pattern      : ", pos_pattern)
message("  out_dir          : ", out_dir)
message("========================================")

################################################################################
# HELPERS
################################################################################
# Find all STAR ReadsPerGene.out.tab files
find_star_files <- function(root) {
  list.files(root, pattern = "ReadsPerGene\\.out\\.tab$",
             recursive = TRUE, full.names = TRUE)
}

# Merge STAR count tables into a single matrix
merge_star_counts <- function(files, use_col) {
  mat <- NULL
  gene_ids <- NULL
  sample_names <- character()
  
  for (f in files) {
    sample <- basename(dirname(f))
    if (sample %in% c("", ".", basename(f))) {
      sample <- str_remove(basename(f), "\\.out\\.tab$")
    }
    
    dt <- fread(f, header = FALSE)
    setnames(dt, c("gene_id", "unstranded", "stranded", "reverse"))
    
    # drop summary rows
    dt <- dt[!startsWith(gene_id, "N_")]
    counts <- dt[[use_col]]
    
    if (is.null(gene_ids)) {
      gene_ids <- dt$gene_id
      mat <- as.matrix(counts)
    } else {
      if (!identical(gene_ids, dt$gene_id)) {
        stop("Gene ID order mismatch across STAR files: ", f)
      }
      mat <- cbind(mat, as.matrix(counts))
    }
    sample_names <- c(sample_names, sample)
  }
  colnames(mat) <- sample_names
  rownames(mat) <- gene_ids
  mat
}

# Hardened GTF parser (reads once)
parse_gtf_features <- function(gtf) {
  message("Reading GTF: ", gtf)
  con <- file(gtf, open = "r"); on.exit(close(con))
  
  gene_rows <- new.env(parent = emptyenv())
  exon_spans <- new.env(parent = emptyenv())
  cds_spans  <- new.env(parent = emptyenv())
  gene_spans <- new.env(parent = emptyenv())
  i <- 0L
  
  while (length(line <- readLines(con, n = 1L)) > 0) {
    if (startsWith(line, "#")) next
    fields <- strsplit(line, "\t")[[1]]
    if (length(fields) < 9) next
    
    chr <- fields[1]
    feat <- fields[3]
    s <- as.integer(fields[4])
    e <- as.integer(fields[5])
    attrs <- fields[9]
    
    gid <- stringr::str_match(attrs, 'gene_id "([^"]+)"')[, 2]
    if (is.na(gid)) next
    
    if (feat == "transcript") {
      gname <- stringr::str_match(attrs, 'gene_name "([^"]+)"')[, 2]
      gdesc <- stringr::str_match(attrs, 'gene_desc "([^"]+)"')[, 2]
      
      gene_rows[[gid]] <- list(
        gene_id   = gid,
        gene_name = ifelse(is.na(gname), "", gname),
        gene_desc = ifelse(is.na(gdesc), "", gdesc)
      )
      
      gs <- get0(gid, gene_spans, ifnotfound = NULL)
      cur <- list(chr = chr, s = s, e = e)
      if (is.null(gs)) {
        assign(gid, cur, envir = gene_spans)
      } else if (gs$chr == chr) {
        gs$s <- min(gs$s, s)
        gs$e <- max(gs$e, e)
        assign(gid, gs, envir = gene_spans)
      }
    }
    
    if (feat == "exon" || feat == "CDS") {
      dt <- data.table(chr = chr, start = s, end = e)
      env <- if (feat == "exon") exon_spans else cds_spans
      prev <- get0(gid, env, ifnotfound = NULL)
      if (is.null(prev)) {
        assign(gid, dt, envir = env)
      } else {
        assign(gid, rbind(prev, dt), envir = env)
      }
    }
    
    i <- i + 1L
    if (i %% 500000 == 0) {
      message("...processed ", i, " GTF lines")
    }
  }
  
  merge_intervals <- function(dt) {
    if (is.null(dt) || nrow(dt) == 0) return(NULL)
    dt <- as.data.table(dt)
    dt[, `:=`(
      start = as.integer(start),
      end   = as.integer(end),
      chr   = as.character(chr)
    )]
    setorder(dt, chr, start, end)
    
    res <- dt[, {
      cur_s <- start[1]
      cur_e <- end[1]
      out_s <- integer()
      out_e <- integer()
      
      if (.N > 1) {
        for (k in 2:.N) {
          if (start[k] <= cur_e + 1L) {
            cur_e <- max(cur_e, end[k])
          } else {
            out_s <- c(out_s, cur_s)
            out_e <- c(out_e, cur_e)
            cur_s <- start[k]
            cur_e <- end[k]
          }
        }
      }
      out_s <- c(out_s, cur_s)
      out_e <- c(out_e, cur_e)
      .(start = out_s, end = out_e)
    }, by = chr]
    
    res[]
  }
  
  gids <- ls(gene_rows)
  gene_dt <- data.table(
    gene_id   = gids,
    gene_name = vapply(gids, function(g) gene_rows[[g]]$gene_name, character(1)),
    gene_desc = vapply(gids, function(g) gene_rows[[g]]$gene_desc, character(1))
  )
  
  exon_len <- integer(length(gids))
  cds_len  <- integer(length(gids))
  has_cds  <- logical(length(gids))
  
  for (ii in seq_along(gids)) {
    g <- gids[ii]
    ex <- merge_intervals(get0(g, exon_spans, ifnotfound = NULL))
    cd <- merge_intervals(get0(g, cds_spans,  ifnotfound = NULL))
    
    exon_len[ii] <- if (!is.null(ex)) sum(ex$end - ex$start + 1L) else 0L
    cds_len[ii]  <- if (!is.null(cd)) sum(cd$end - cd$start + 1L) else 0L
    has_cds[ii]  <- cds_len[ii] > 0L
  }
  
  gene_dt[, `:=`(
    exon_length = exon_len,
    cds_length  = cds_len,
    has_CDS     = has_cds
  )]
  
  cds_tab <- rbindlist(lapply(gids, function(g) {
    cd <- merge_intervals(get0(g, cds_spans, ifnotfound = NULL))
    if (is.null(cd)) return(NULL)
    cd[, gene_id := g]
    cd[]
  }), fill = TRUE)
  
  if (!is.null(cds_tab)) {
    setcolorder(cds_tab, c("chr", "start", "end", "gene_id"))
  }
  
  gene_span_tab <- rbindlist(lapply(gids, function(g) {
    gs <- get0(g, gene_spans, ifnotfound = NULL)
    if (is.null(gs)) return(NULL)
    data.table(chr  = gs$chr,
               start = as.integer(gs$s),
               end   = as.integer(gs$e),
               gene_id = g)
  }), fill = TRUE)
  
  list(genes = gene_dt, cds = cds_tab, gene_spans = gene_span_tab)
}

# Read one positions file into a standard data.table(chr,start,end,pos)
read_positions <- function(file) {
  pos_raw <- fread(file, header = FALSE)
  # header detection: if 2nd col of row1 isn't an integer, assume header
  has_header <- suppressWarnings(is.na(as.integer(pos_raw[[2]][1])))
  
  if (has_header) {
    pos <- fread(file, header = TRUE)
    setnames(pos, old = names(pos)[1:2], new = c("chrom", "pos"))
  } else {
    pos <- pos_raw[, .(chrom = as.character(V1), pos = as.integer(V2))]
  }
  
  pos <- pos[!is.na(chrom) & !is.na(pos)]
  pos[, `:=`(start = as.integer(pos), end = as.integer(pos))]
  setnames(pos, c("chrom", "start", "end", "pos"))
  setcolorder(pos, c("chrom", "start", "end", "pos"))
  pos
}

# Overlap + write outputs for a single positions table
process_positions <- function(pos_file, feat_dt, flags_dt, out_dir, region_type) {
  prefix <- tools::file_path_sans_ext(basename(pos_file))
  snp_dt <- read_positions(pos_file)
  
  if (nrow(snp_dt) == 0) {
    message("[", prefix, "] No positions found; skipping.")
    return(invisible(NULL))
  }
  
  # keys for foverlaps
  snp_dt2 <- copy(snp_dt)
  setnames(snp_dt2, c("chrom", "start", "end"), c("chr", "start", "end"))
  setkey(snp_dt2, chr, start, end)
  
  hits <- foverlaps(snp_dt2, feat_dt, nomatch = 0L)
  if (nrow(hits) == 0) {
    message("[", prefix, "] 0 overlaps with expressed ", toupper(region_type), " features.")
  }
  
  hits <- merge(
    hits,
    flags_dt[, .(gene_id, gene_name, gene_desc)],
    by = "gene_id", all.x = TRUE
  )
  
  # outputs
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  out_tsv <- file.path(out_dir,
                       sprintf("%s.positions_expressed_%s.tsv",
                               prefix, tolower(region_type)))
  out_bed <- file.path(out_dir,
                       sprintf("%s.positions_expressed_%s.bed",
                               prefix, tolower(region_type)))
  
  fwrite(hits[, .(chrom = chr, pos = pos, gene_id, gene_name, gene_desc)],
         out_tsv, sep = "\t")
  
  bed <- hits[, .(
    chrom  = chr,
    start  = as.integer(pos) - 1L,
    end    = as.integer(pos),
    name   = fifelse(is.na(gene_name) | gene_name == "",
                     gene_id,
                     paste0(gene_id, ";", gene_name)),
    score  = 0L,
    strand = "."
  )]
  
  bed[start < 0L, start := 0L]
  setorder(bed, chrom, start, end)
  fwrite(bed, out_bed, sep = "\t", col.names = FALSE)
  
  message(sprintf("[%s] Wrote %d rows: %s ; %s",
                  prefix, nrow(hits), basename(out_tsv), basename(out_bed)))
  
  invisible(list(tsv = out_tsv, bed = out_bed, n = nrow(hits)))
}
################################################################################
# MAIN PIPELINE
################################################################################
# 1) STAR counts -> CPM
star_files <- find_star_files(star_root)
if (!length(star_files)) {
  stop("No ReadsPerGene.out.tab files found under: ", star_root)
}
message("Found ", length(star_files),
        " STAR count files. Using column ", use_star_column, ".")

counts <- merge_star_counts(star_files, use_star_column)

# 2) GTF
gtf <- parse_gtf_features(gtf_file)

# 3) Align annotation to counts order
ann <- gtf$genes[match(rownames(counts), gtf$genes$gene_id)]
stopifnot(nrow(ann) == nrow(counts))

# 4) CPM + expressed
dge <- DGEList(counts = counts)
cpm_mat <- cpm(dge)

expressed_flag <- apply(cpm_mat, 1, function(x) {
  sum(x >= cpm_threshold) >= min_samples_at_cutoff
})

flags_dt <- data.table(
  gene_id       = ann$gene_id,
  gene_name     = ann$gene_name,
  gene_desc     = ann$gene_desc,
  has_CDS       = ann$has_CDS,
  exon_length   = ann$exon_length,
  cds_length    = ann$cds_length,
  expressed_any = expressed_flag
)

# 5) Choose feature set once, keep only expressed genes
if (region_type == "CDS") {
  feat <- gtf$cds
  if (is.null(feat) || nrow(feat) == 0) {
    stop("No CDS intervals found in GTF.")
  }
} else if (region_type == "GENE") {
  feat <- gtf$gene_spans
  if (is.null(feat) || nrow(feat) == 0) {
    stop("No gene span intervals found in GTF.")
  }
} else {
  stop('region_type must be "CDS" or "gene"')
}

expressed_ids <- flags_dt[has_CDS == TRUE & expressed_any == TRUE, gene_id]
feat_dt <- as.data.table(feat)[gene_id %in% expressed_ids]
setnames(feat_dt, c("chr", "start", "end", "gene_id"))
setkey(feat_dt, chr, start, end)  # key once; reused for all files

# 6) Iterate all .pos files and write outputs with their prefixes
pos_files <- list.files(positions_dir, pattern = pos_pattern, full.names = TRUE)
if (!length(pos_files)) {
  stop("No position files found matching pattern in: ", positions_dir)
}

message("Processing ", length(pos_files), " position files in ", positions_dir, " ...")
invisible(lapply(pos_files, function(f) {
  try(process_positions(f, feat_dt, flags_dt, out_dir, region_type), silent = FALSE)
}))
message("All done.")

################################################################################
# SUMMARY: GENE EXPRESSION ONLY 
################################################################################
n_genes_all    <- nrow(flags_dt)
n_genes_cds    <- sum(flags_dt$has_CDS, na.rm = TRUE)
n_expr_any_all <- sum(flags_dt$expressed_any, na.rm = TRUE)
n_expr_any_cds <- sum(flags_dt$has_CDS & flags_dt$expressed_any, na.rm = TRUE)

cat(sprintf(
  "Genes: %d total (CDS-bearing: %d). Expressed (CPM ≥ %.3g in ≥%d sample[s]): %d total; %d with CDS.\n",
  n_genes_all, n_genes_cds,
  cpm_threshold, min_samples_at_cutoff,
  n_expr_any_all, n_expr_any_cds
))

