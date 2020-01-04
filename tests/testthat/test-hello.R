context("test-hello")

test_that("multiplication works", {
  data("demo_data")
  res<-IdentifyMCB(demo_data$realdata)
  expect_equal(as.character(res$MCBinformation[,'location']), 
   c("chr2 : 74669349 - chr2 : 74669516","chr6 : 108444791 - chr6 : 108444920",
   "chr8 : 18541446 - chr8 : 18541627","chr8 : 61777711 - chr8 : 61778137"))
})
