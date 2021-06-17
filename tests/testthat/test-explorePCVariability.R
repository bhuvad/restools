test_that("explorePCVariability works with Eset", {
  library(ALL)
  data(ALL)

  #Default param
  expect_error(
    explorePCVariability(ALL, c()),
    'No annotations present'
  )

  #assess with column
  expect_silent(explorePCVariability(ALL, age))
  expect_silent(explorePCVariability(ALL, c(age)))

  #multiple annotations
  expect_silent(explorePCVariability(ALL, c(age, sex)))


  #error when the column is not in the data
  expect_error(
    explorePCVariability(ALL, age2),
    "Column `age2` doesn't exist"
  )
})
