test_that("MCBsites CpG names match those listed in MCBinformation", {
  data("demo_data")
  res <- IdentifyMCB(demo_data$realdata, verbose = FALSE)

  mcb_cpgs_from_info <- unique(unlist(strsplit(res$MCBinformation[, "CpGs"], " ")))

  # Every site in MCBsites must actually belong to an MCB block.
  # Before the fix, rownames(MethylationProfile) was indexed by correlation_res
  # positions, returning CpG names in the wrong (annotation) order.
  expect_true(
    all(res$MCBsites %in% mcb_cpgs_from_info),
    label = "all MCBsites appear in MCBinformation$CpGs"
  )

  # Counts must agree
  expect_equal(length(res$MCBsites), length(mcb_cpgs_from_info))
})

test_that("custom platform: no error and MCBsites are valid CpG names", {
  set.seed(6271)
  cpg_names <- paste0("cg", sprintf("%08d", 1:10))

  # Minimal annotation with required columns
  anno <- data.frame(
    chr = c(rep("chr1", 6), rep("chr2", 4)),
    pos = c(1000, 1200, 1400, 5000, 7000, 9000, 1000, 1100, 1200, 5000),
    UCSC_RefGene_Name = "",
    Relation_to_Island = "OpenSea",
    Islands_Name = "",
    row.names = cpg_names,
    stringsAsFactors = FALSE
  )

  n_samples <- 50
  meth <- matrix(
    runif(10 * n_samples), nrow = 10,
    dimnames = list(cpg_names, paste0("s", seq_len(n_samples)))
  )
  # Make first 3 CpGs highly correlated so an MCB forms
  base_signal <- runif(n_samples)
  meth[1, ] <- base_signal
  meth[2, ] <- base_signal + rnorm(n_samples, 0, 0.02)
  meth[3, ] <- base_signal + rnorm(n_samples, 0, 0.02)

  # Before the fix, passing a data.frame as platform caused:
  #   "the condition has length > 1" (from platform == "Illumina Methylation 450K")
  # followed by "object 'intersect_cpg' not found" in the else branch.
  expect_no_error(
    res <- IdentifyMCB(meth, platform = anno, verbose = FALSE)
  )

  # All returned CpG sites should be valid names from the custom annotation
  expect_true(all(res$MCBsites %in% cpg_names))

  # MCBsites should match CpGs listed in MCBinformation
  mcb_cpgs_from_info <- unique(unlist(strsplit(res$MCBinformation["CpGs"], " ")))
  expect_true(all(res$MCBsites %in% mcb_cpgs_from_info))
})

test_that("custom platform: CpGs absent from MethylationProfile are silently dropped", {
  set.seed(6271)
  cpg_names <- paste0("cg", sprintf("%08d", 1:10))

  anno <- data.frame(
    chr = c(rep("chr1", 6), rep("chr2", 4)),
    pos = c(1000, 1200, 1400, 5000, 7000, 9000, 1000, 1100, 1200, 5000),
    UCSC_RefGene_Name = "",
    Relation_to_Island = "OpenSea",
    Islands_Name = "",
    row.names = cpg_names,
    stringsAsFactors = FALSE
  )

  # Methylation profile is missing the last 3 CpGs from the annotation
  n_samples <- 50
  meth <- matrix(
    runif(7 * n_samples), nrow = 7,
    dimnames = list(cpg_names[1:7], paste0("s", seq_len(n_samples)))
  )
  base_signal <- runif(n_samples)
  meth[1, ] <- base_signal
  meth[2, ] <- base_signal + rnorm(n_samples, 0, 0.02)
  meth[3, ] <- base_signal + rnorm(n_samples, 0, 0.02)

  expect_no_error(
    res <- IdentifyMCB(meth, platform = anno, verbose = FALSE)
  )
  # Result should only contain CpGs present in the methylation profile
  expect_true(all(res$MCBsites %in% cpg_names[1:7]))
})

test_that("each MCB block carries independent metadata (no stale field bleed-over)", {
  # Two chromosomes, each with one MCB that has a distinct Islands_Name.
  # Without the MCB[] <- NA reset after recording each block, any field that
  # is conditionally not written in a future code change would silently retain
  # the previous block's value. This test guards against that regression.
  set.seed(3847)
  n_samples <- 60
  cpg2 <- paste0("cg", sprintf("%08d", 1:8))

  anno_two <- data.frame(
    chr               = c(rep("chr1", 4), rep("chr2", 4)),
    pos               = c(1000, 1200, 1400, 1600,
                          2000, 2200, 2400, 2600),
    UCSC_RefGene_Name = rep("", 8),
    Relation_to_Island = rep("OpenSea", 8),
    Islands_Name      = c(rep("IslandA", 4), rep("IslandB", 4)),
    row.names         = cpg2,
    stringsAsFactors  = FALSE
  )

  meth_two <- matrix(
    runif(8 * n_samples), nrow = 8,
    dimnames = list(cpg2, paste0("s", seq_len(n_samples)))
  )
  base1 <- runif(n_samples); base2 <- runif(n_samples)
  for (r in 1:4) meth_two[r, ] <- base1 + rnorm(n_samples, 0, 0.02)
  for (r in 5:8) meth_two[r, ] <- base2 + rnorm(n_samples, 0, 0.02)

  res <- IdentifyMCB(meth_two, platform = anno_two, verbose = FALSE)

  info <- res$MCBinformation
  expect_equal(nrow(info), 2)

  # Each block must carry the CGI_Coordinate from its own CpGs, not the
  # preceding block's value.
  block_chr1 <- info[info[, "chromosomes"] == "chr1", "CGI_Coordinate"]
  block_chr2 <- info[info[, "chromosomes"] == "chr2", "CGI_Coordinate"]
  expect_equal(block_chr1, "IslandA")
  expect_equal(block_chr2, "IslandB")

  # Chromosomes must also be independent
  expect_equal(block_chr1 != block_chr2, TRUE)
})

test_that("MCBsites are unchanged when input rows are in a different order", {
  data("demo_data")

  res_original <- IdentifyMCB(demo_data$realdata, verbose = FALSE)

  # Reverse the row order to force a different input ordering
  reversed <- demo_data$realdata[rev(seq_len(nrow(demo_data$realdata))), ]
  res_reversed <- IdentifyMCB(reversed, verbose = FALSE)

  expect_equal(sort(res_original$MCBsites), sort(res_reversed$MCBsites))
})
