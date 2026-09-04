# Run from the repository root after preparing the input workbook.
library(readxl)
library(openxlsx)
source("R/calculate_mdtps.R", encoding = "UTF-8")

INPUT <- file.path("data", "mdtps", "gene_score_input.xlsx")
INPUT_SHEET <- 1  # Change to the actual worksheet name if needed.
OUTDIR <- file.path("results", "mdtps")
stopifnot(file.exists(INPUT))
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

dat <- read_excel(INPUT, sheet = INPUT_SHEET)
result <- calculate_mdtps(dat = dat)
main_cols <- c("gene", intersect(c("cell", "cancer"), names(result)),
               "n_supported_evidence", "B", "m_direction", "R", "C",
               "overall_direction", "D", "BaseScore", "DirectionScore",
               "DrugBonus", "MDTPS", "Rank")
write.xlsx(list(Full_results = result, Main_results = result[, main_cols]),
           file = file.path(OUTDIR, "MDTPS_results.xlsx"), rowNames = FALSE)
capture.output(sessionInfo(), file = file.path(OUTDIR, "sessionInfo.txt"))
message("Results saved in ", OUTDIR)
