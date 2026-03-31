# Dane z tabeli
df <- data.frame(
  grSignature = c("globalUp", "globalDown", "NeuralCellsUp", "NeuralCellsDown",
                  "BloodCellsUp", "BloodCellsDown", "LungCellsUp", "LungCellsDown"),
  genesInSignature = c(208, 208, 147, 147, 124, 124, 208, 208),
  n_genes_p0.05 = c(0, 8, 25, 0, 30, 11, 9, 0),
  n_genes_p0.01 = c(0, 0, 12, 0, 16, 5, 0, 0)
)

# ---- Test dla progu p < 0.05 ----

# Suma globalnych i tkankowych
global_assoc <- sum(df$n_genes_p0.05[df$grSignature %in% c("globalUp", "globalDown")])
global_notassoc <- sum(df$genesInSignature[df$grSignature %in% c("globalUp", "globalDown")]) - global_assoc

tissue_assoc <- sum(df$n_genes_p0.05[!df$grSignature %in% c("globalUp", "globalDown")])
tissue_notassoc <- sum(df$genesInSignature[!df$grSignature %in% c("globalUp", "globalDown")]) - tissue_assoc

tbl <- matrix(c(tissue_assoc, tissue_notassoc,
                global_assoc, global_notassoc),
              nrow = 2, byrow = TRUE)

rownames(tbl) <- c("Tissue", "Global")
colnames(tbl) <- c("Associated", "Not_associated")

cat("=== Observed counts (p < 0.05) ===\n")
print(tbl)

# Test chi2
chi_result <- chisq.test(tbl, correct = FALSE)

cat("\n=== Expected counts ===\n")
print(round(chi_result$expected, 2))

cat("\n=== Chi-squared result ===\n")
print(chi_result)

# Iloraz szans (OR)
a <- tbl[1,1]; b <- tbl[1,2]; c <- tbl[2,1]; d <- tbl[2,2]
OR <- (a*d) / (b*c)
SE_log_OR <- sqrt(1/a + 1/b + 1/c + 1/d)
CI_low <- exp(log(OR) - 1.96 * SE_log_OR)
CI_high <- exp(log(OR) + 1.96 * SE_log_OR)

cat("\n=== Odds Ratio (OR) ===\n")
cat(sprintf("OR = %.3f (95%% CI: %.3f–%.3f)\n", OR, CI_low, CI_high))

# ---- Test dla progu p < 0.01 ----
cat("\n\n----------------------------\n")
cat("==== Now test for p < 0.01 ====\n")

global_assoc_01 <- sum(df$n_genes_p0.01[df$grSignature %in% c("globalUp", "globalDown")])
global_notassoc_01 <- sum(df$genesInSignature[df$grSignature %in% c("globalUp", "globalDown")]) - global_assoc_01

tissue_assoc_01 <- sum(df$n_genes_p0.01[!df$grSignature %in% c("globalUp", "globalDown")])
tissue_notassoc_01 <- sum(df$genesInSignature[!df$grSignature %in% c("globalUp", "globalDown")]) - tissue_assoc_01

tbl2 <- matrix(c(tissue_assoc_01, tissue_notassoc_01,
                 global_assoc_01, global_notassoc_01),
               nrow = 2, byrow = TRUE)

rownames(tbl2) <- c("Tissue", "Global")
colnames(tbl2) <- c("Associated", "Not_associated")

cat("\n=== Observed counts (p < 0.01) ===\n")
print(tbl2)

chi_result_01 <- chisq.test(tbl2, correct = FALSE)

cat("\n=== Expected counts ===\n")
print(round(chi_result_01$expected, 2))

cat("\n=== Chi-squared result ===\n")
print(chi_result_01)

# Iloraz szans (OR)
a <- tbl2[1,1]; b <- tbl2[1,2]; c <- tbl2[2,1]; d <- tbl2[2,2]
OR <- (a*d) / (b*c)
SE_log_OR <- sqrt(1/a + 1/b + 1/c + 1/d)
CI_low <- exp(log(OR) - 1.96 * SE_log_OR)
CI_high <- exp(log(OR) + 1.96 * SE_log_OR)

cat("\n=== Odds Ratio (OR) ===\n")
cat(sprintf("OR = %.3f (95%% CI: %.3f–%.3f)\n", OR, CI_low, CI_high))
