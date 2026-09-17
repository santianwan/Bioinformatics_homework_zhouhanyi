# =============================================================================
# Week 4, Question 4 — variant prioritization
# Completed from: Week 4/Homework/for_student/starter/q4_filter_starter.R
#
# Input : ../data/variants_q4.tsv   (synthetic teaching table, 12 variants)
# Output: ../results/q4_shortlist.tsv        surviving candidates, ranked
#         ../results/q4_filter_trace.tsv     per-variant reason for keep/drop
#
# Design note: every variant is scored on ALL criteria first and the reason for
# exclusion is recorded, rather than being silently dropped by a chained
# filter(). In a real triage the dropped rows are the ones you have to defend,
# so they must stay visible.
# =============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

variants_path <- file.path("..", "data", "variants_q4.tsv")
out_dir       <- file.path("..", "results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

variants <- read_tsv(variants_path, comment = "#", show_col_types = FALSE)
stopifnot(nrow(variants) == 12L,
          all(c("CHROM","POS","REF","ALT","FILTER","DP","GQ","AF",
                "GENE","CONSEQUENCE","CLINVAR_SIG","CLINVAR_ID") %in% names(variants)))

# --- Thresholds -------------------------------------------------------------
# Chosen before writing the code, and justified in the report:
min_dp       <- 20      # GATK best-practice territory for confident germline
                        # genotyping; below ~20x a het can be missed or invented.
min_gq       <- 30      # PHRED-scaled: P(wrong genotype) < 1e-3.
max_af       <- 0.001   # A penetrant dominant disease allele cannot be common.
                        # 0.1% is deliberately permissive vs the ~1e-4 used for
                        # rare dominant disease, so nothing borderline is lost.
require_pass <- TRUE    # FILTER != PASS means the caller itself distrusts it.

# Consequences that can plausibly change protein product or transcript level.
# Ordered by expected severity — used for ranking, not only for filtering.
impact_consequences <- c(
  "stop_gained",
  "frameshift_variant",
  "splice_acceptor_variant",
  "splice_donor_variant",
  "missense_variant"
)
severity_rank <- setNames(seq_along(impact_consequences), impact_consequences)

# ClinVar labels, as evidence weight. "Conflicting" is NOT treated as evidence
# for pathogenicity; it is treated as absence of consensus.
clinvar_weight <- c(
  "Pathogenic"                                   =  3,
  "Likely_pathogenic"                            =  2,
  "Uncertain_significance"                       =  1,
  "Conflicting_interpretations_of_pathogenicity" =  1,
  "Likely_benign"                                = -2,
  "Benign"                                       = -3
)

# --- Score every variant on every criterion (no silent drops) ---------------
trace <- variants %>%
  mutate(
    pass_filter      = FILTER == "PASS",
    pass_dp          = DP >= min_dp,
    pass_gq          = GQ >= min_gq,
    pass_af          = AF <= max_af,
    pass_consequence = CONSEQUENCE %in% impact_consequences,
    has_gene         = GENE != "." & !is.na(GENE),
    cv_weight        = ifelse(CLINVAR_SIG %in% names(clinvar_weight),
                              clinvar_weight[CLINVAR_SIG], 0),
    not_benign       = cv_weight > -2,
    kept             = pass_filter & pass_dp & pass_gq & pass_af &
                       pass_consequence & has_gene & not_benign,
    drop_reason = case_when(
      !pass_filter      ~ paste0("FILTER=", FILTER, " (caller distrusts the call)"),
      !pass_dp          ~ paste0("DP=", DP, " < ", min_dp, " (genotype not supportable)"),
      !pass_gq          ~ paste0("GQ=", GQ, " < ", min_gq),
      !pass_af          ~ paste0("AF=", AF, " > ", max_af, " (too common for a penetrant allele)"),
      !pass_consequence ~ paste0("CONSEQUENCE=", CONSEQUENCE, " (no expected protein/transcript effect)"),
      !has_gene         ~ "no gene symbol assigned (no interpretable target)",
      !not_benign       ~ paste0("ClinVar ", CLINVAR_SIG, " (positive evidence against)"),
      TRUE              ~ NA_character_
    ),
    severity = ifelse(CONSEQUENCE %in% names(severity_rank),
                      severity_rank[CONSEQUENCE], 99L)
  )

shortlist <- trace %>%
  filter(kept) %>%
  arrange(severity, desc(cv_weight), AF, desc(GQ)) %>%
  mutate(rank = row_number()) %>%
  select(rank, CHROM, POS, REF, ALT, GENE, CONSEQUENCE, CLINVAR_SIG, CLINVAR_ID,
         FILTER, DP, GQ, AF, severity, cv_weight)

# --- Report -----------------------------------------------------------------
cat("Input variants:            ", nrow(variants), "\n")
cat("Failed FILTER:             ", sum(!trace$pass_filter), "\n")
cat("Failed DP >= ", min_dp, ":         ", sum(!trace$pass_dp), "\n", sep = "")
cat("Failed GQ >= ", min_gq, ":         ", sum(!trace$pass_gq), "\n", sep = "")
cat("Failed AF <= ", max_af, ":     ", sum(!trace$pass_af), "\n", sep = "")
cat("Non-impactful consequence: ", sum(!trace$pass_consequence), "\n")
cat("No gene symbol:            ", sum(!trace$has_gene), "\n")
cat("ClinVar (Likely_)benign:   ", sum(!trace$not_benign), "\n")
cat("Surviving shortlist:       ", nrow(shortlist), "\n\n")

print(as.data.frame(shortlist))
cat("\n--- dropped, with reason ---\n")
print(as.data.frame(trace %>% filter(!kept) %>%
        select(CHROM, POS, GENE, CONSEQUENCE, drop_reason)))

write_tsv(shortlist, file.path(out_dir, "q4_shortlist.tsv"))
write_tsv(trace %>% select(CHROM, POS, REF, ALT, GENE, CONSEQUENCE, CLINVAR_SIG,
                           FILTER, DP, GQ, AF, kept, drop_reason),
          file.path(out_dir, "q4_filter_trace.tsv"))

cat("\nSession info:\n")
print(sessionInfo()$R.version$version.string)
print(vapply(c("readr","dplyr"), function(p) as.character(packageVersion(p)), ""))
