test_that("plotPCA works with Eset", {
  library(ALL)
  data(ALL)

  #Default param
  p <- plotPCA(ALL)
  expect_silent(print(p))

  #color with column
  p <- plotPCA(ALL, color=age)
  expect_silent(print(p))

  #color with expression
  p <- plotPCA(ALL, color=age > 50)
  expect_silent(print(p))

  #multiple aesthetics
  p <- plotPCA(ALL, color=age > 50, shape=sex)
  expect_error(print(p), NA)

  #error when the column is not in the data
  p <- plotPCA(ALL, color=age2)
  expect_error(
    # We need to force eval here to throw the error
    print(p),
    "object 'age2' not found"
  )

  #will not error if the variable is present in the parent env
  age2 <- ALL$age
  p <- plotPCA(ALL, color=age2)
  expect_error(print(p), NA)
})

test_that("plotPCA row subset works as expected", {

  cu_se <- emtdata::cursons2018_se()
  gene_subset <- sample(rownames(cu_se), 50)

  expect_message( # not setting works but gives a message
    plt <- plotPCA(cu_se, assay = "logRPKM", colour = Subline),
    "Using the scater default of 500 features as input"
  )

  expect_silent(
    plt <- plotPCA(cu_se, assay = "logRPKM", subset_row = gene_subset, colour = Subline) # 50 stable genes
  )

  expect_error( # expect error for missing features
    plt <- plotPCA(cu_se, assay = "logRPKM", subset_row = c("bad", "names"), colour = Subline)
    )

})
