alpha <- read.csv("F:/subject/3_up/bioinformation/Bioinformatics_SUAT_2026_FALL-main/Week 6/Homework/unzipped_20260925-024344/tables/02_alpha_indices.csv",
                   check.names = FALSE)
meta <- read.csv("F:/EasyMultiProfiler-Web-9.0.4/EasyMultiProfiler-Web-9.0.4/tests/16S_mapping.csv", check.names = FALSE)

names(alpha)[1] <- "SampleID"
df <- merge(alpha, meta, by = "SampleID")
cat("merged rows:", nrow(df), "of", nrow(alpha), "alpha rows and", nrow(meta), "meta rows\n")

df$disease   <- sub("_.*", "", df$Group)
df$timepoint <- ifelse(grepl("after", df$Group), "after", "before")
df$response  <- sub(".*_", "", df$Group_sub)

cat("\n=== web tool's own Shannon, by group ===\n")
print(aggregate(shannon ~ disease + timepoint, data = df, FUN = function(x) round(mean(x), 3)))

cat("\nKruskal-Wallis (web tool shannon) ~ disease:  p =",
    format.pval(kruskal.test(shannon ~ disease, data = df)$p.value, digits = 3), "\n")
cat("Kruskal-Wallis (web tool shannon) ~ timepoint:p =",
    format.pval(kruskal.test(shannon ~ timepoint, data = df)$p.value, digits = 3), "\n")
cat("Kruskal-Wallis (web tool shannon) ~ response: p =",
    format.pval(kruskal.test(shannon ~ response, data = df)$p.value, digits = 3), "\n")

for (dz in c("IBS", "UC")) {
  sub <- df[df$disease == dz, ]
  kw <- kruskal.test(shannon ~ timepoint, data = sub)
  cat(sprintf("  web-tool Shannon ~ timepoint within disease=%s: p = %s (n=%d)\n",
              dz, format.pval(kw$p.value, digits = 3), nrow(sub)))
}

cat("\ncorrelation between my reference Shannon and the web tool's Shannon, same samples:\n")
ref <- read.csv("F:/subject/3_up/bioinformation/Bioinformatics_homework_zhouhanyi/Week6/outputs/week6_alpha_diversity.csv")
names(ref)[1] <- "SampleID"
cmp <- merge(df[, c("SampleID", "shannon")], ref[, c("SampleID", "shannon")], by = "SampleID",
             suffixes = c("_webtool", "_reference"))
cat("n compared:", nrow(cmp), "\n")
print(cor.test(cmp$shannon_webtool, cmp$shannon_reference, method = "spearman"))
