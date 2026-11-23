# Shared functions for Stage 3 (subtyping preparation) and the library builder.
suppressPackageStartupMessages({
  library(Biostrings)
  library(pwalign)
  library(dplyr)
  library(tibble)
})

# Translate each nucleotide sequence in the best of six frames, where "best" means the frame whose
# protein has the highest local alignment score (BLOSUM62, gap 10/1 as in the manual mafft step)
# to any reference. No motif is required, so a real KoRV lacking e.g. HLG is still translated.
# Returns a tibble with the oriented nucleotide sequence, protein ('*' for stops), strand, frame.
# n_frame_refs: the frame is identified against an evenly spread subset of the references
# (all subtypes are similar enough that any of them finds the right frame), which keeps this fast.
translate_best_frame <- function(nt, refs_aa, n_frame_refs = 4, cores = 1) {
  data("BLOSUM62", package = "pwalign", envir = environment())
  idx <- unique(round(seq(1, length(refs_aa), length.out = min(n_frame_refs, length(refs_aa)))))
  ref_list <- as.list(refs_aa[idx])
  rows <- parallel::mclapply(seq_along(nt), mc.cores = cores, FUN = function(i) {
    s <- nt[[i]]
    cands <- list()
    for (strand in c("+", "-")) {
      x <- if (strand == "+") s else reverseComplement(s)
      for (frame in 1:3) {
        y <- subseq(x, frame, frame + 3 * ((length(x) - frame + 1) %/% 3) - 1)
        cands[[length(cands) + 1]] <- list(strand = strand, frame = frame, nt = x,
                                           aa = translate(y, if.fuzzy.codon = "X"))
      }
    }
    aas <- AAStringSet(lapply(cands, `[[`, "aa"))
    best_per_ref <- vapply(ref_list, function(r)
      pairwiseAlignment(aas, r, type = "local", substitutionMatrix = BLOSUM62,
                        gapOpening = 10, gapExtension = 1, scoreOnly = TRUE), numeric(length(aas)))
    best_per_ref <- matrix(best_per_ref, nrow = length(aas))
    sc <- apply(best_per_ref, 1, max)
    k <- which.max(sc)
    tibble(id = names(nt)[i], strand = cands[[k]]$strand, frame = cands[[k]]$frame,
           nt_oriented = as.character(cands[[k]]$nt), protein = as.character(cands[[k]]$aa),
           frame_score = sc[k], runner_up_score = if (length(sc) > 1) sort(sc, TRUE)[2] else NA_real_)
  })
  bind_rows(rows)
}

# Number of stop codons that are not the final residue
n_internal_stops <- function(protein) {
  p <- sub("\\*$", "", protein)
  nchar(gsub("[^*]", "", p))
}

# Names safe for mafft / iqtree
safe_name <- function(x) gsub("[^A-Za-z0-9_.-]", "_", x)

# Subtype letter of an A-N reference name, e.g. "KoRV-A_[Hanger_2000]_(+3)" -> "A", "D35_(+3)" -> "D"
reference_subtype <- function(x) toupper(substr(sub("^KoRV-", "", x), 1, 1))

# Alignment columns of an anchored region: the residues between two anchor motifs in a given
# reference (used for the NQQHNN zone and for the variable region).
# R: reference rows of the alignment matrix, ref_safe: their (safe) names,
# vr: a rules entry with reference, left_anchor and right_anchor
anchored_region_columns <- function(R, ref_safe, vr) {
  k <- match(safe_name(vr$reference), ref_safe)
  if (is.na(k)) return(integer(0))
  vref <- R[k, ]; cols <- which(vref != "-")
  ung <- paste(vref[cols], collapse = "")
  l <- stringr::str_locate(ung, stringr::fixed(vr$left_anchor)); r <- stringr::str_locate(ung, stringr::fixed(vr$right_anchor))
  if (anyNA(l) || anyNA(r) || cols[r[1]] - cols[l[2]] <= 1) return(integer(0))
  (cols[l[2]] + 1):(cols[r[1]] - 1)
}

# Group used for colouring: all "possible new subtype (nearest X)" suggestions share one colour
group_of <- function(x) ifelse(grepl("^possible new", x), "possible new subtype", x)

# One colour per subtype group, shared by the figures and the interactive review page
subtype_palette <- function(groups) {
  base_pal <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#a6761d", "#1f78b4", "#00bcd4", "#33a02c", "#6a3d9a", "#e6ab02", "#8c564b", "#08306b", "#17becf", "#393b79", "#ff7f0e")
  pal <- setNames(rep_len(base_pal, length(groups)), groups)
  pal[intersect(c("Non-functional"), groups)] <- "#7f7f7f"
  pal[intersect(c("Not_KoRV"), groups)] <- "#bdbdbd"
  pal[intersect(c("possible new subtype"), groups)] <- "#e41a1c"
  pal
}

# ClustalX colouring, as described in the Jalview documentation of the Clustal scheme
# (https://www.jalview.org/help/html/colourSchemes/clustal.html). A residue is coloured only when
# its column meets one of the criteria for that residue; percentages are of the residues in the
# column, gaps excluded. Gaps, X (stop codons) and uncoloured residues stay white.
clustalx_fills <- function(aln) {
  col_hex <- c(blue = "#80a0f0", red = "#f01505", magenta = "#c048c0", green = "#15c015",
               orange = "#f09048", cyan = "#15a4a4", pink = "#f08080", yellow = "#c0c000", none = "#ffffff")
  hyd <- c("W", "L", "V", "I", "M", "A", "F", "C", "Y", "H", "P")
  aro <- c("W", "Y", "A", "C", "P", "Q", "F", "H", "I", "L", "M", "V")
  out <- matrix(col_hex[["none"]], nrow(aln), ncol(aln), dimnames = dimnames(aln))
  for (j in seq_len(ncol(aln))) {
    x <- aln[, j]; res <- x[x != "-"]
    if (length(res) == 0) next
    f <- function(set) mean(res %in% set)
    pick <- function(r) {
      if (r == "C" && f("C") > 0.85) return("pink")
      if (r %in% c("A", "C", "I", "L", "M", "F", "W", "V")) return(if (f(hyd) > 0.60) "blue" else "none")
      if (r %in% c("K", "R")) return(if (f(c("K", "R")) > 0.60 || f(c("K", "R", "Q")) > 0.85) "red" else "none")
      if (r == "E") return(if (f(c("K", "R")) > 0.60 || f(c("Q", "E")) > 0.50 || f(c("E", "D")) > 0.50 || f(c("E", "Q", "D")) > 0.85) "magenta" else "none")
      if (r == "D") return(if (f(c("K", "R")) > 0.60 || f(c("D", "E", "N")) > 0.85 || f(c("E", "D")) > 0.50) "magenta" else "none")
      if (r == "N") return(if (f("N") > 0.50 || f(c("N", "D")) > 0.85) "green" else "none")
      if (r == "Q") return(if (f(c("K", "R")) > 0.60 || f(c("Q", "E")) > 0.50 || f(c("Q", "T", "K", "R")) > 0.85) "green" else "none")
      if (r %in% c("S", "T")) return(if (f(hyd) > 0.60 || f(c("T", "S")) > 0.50 || f(c("S", "T")) > 0.85) "green" else "none")
      if (r == "G") return("orange")
      if (r == "P") return("yellow")
      if (r %in% c("H", "Y")) return(if (f(hyd) > 0.60 || f(aro) > 0.85) "cyan" else "none")
      "none"
    }
    for (r in unique(res)) out[x == r, j] <- col_hex[[pick(r)]]
  }
  out
}
