context("rawdata")

test_that("can round trip binary matrices",{
  
  ix <- matrix(sample(c(0L,1L),64,replace=T),nrow = 16,ncol = 4)
  ret_rows <- ceiling(nrow(ix)/8L)
  
  ret <-   snp2raw(ix)
  gret <- matrix(as.integer(rawToBits(ret)),nrow=16,ncol=4)
  expect_equal(gret,ix)
  
})




popcnt_fun <- function(x){
sum(  as.integer(rawToBits(x)))
  
}

test_that("can compute covariance for binary matrices",{
  
  ix <- matrix(sample(c(0L,1L),64,replace=T),nrow = 16,ncol = 4)
  tS <- cov(ix)

  


  ret <-   snp2raw(ix)
  
  bS <- covbin(ret)

  expect_equal(bS,tS)

  
})


test_that("can round trip binary matrices not multiples of 8",{
  
  ix <- matrix(sample(c(0L,1L),4*17,replace=T),nrow = 17,ncol = 4)
  ret_rows <- ceiling(nrow(ix)/8L)
  ret <-   snp2raw(ix)
  gret <- matrix(as.integer(rawToBits(ret)),nrow=8*3,ncol=4)
  expect_equal(gret[1:17,],ix)
  
  ix <- matrix(sample(c(0L,1L),4*19,replace=T),nrow = 19,ncol = 4)
  ret_rows <- ceiling(nrow(ix)/8L)
  ret <-   snp2raw(ix)
  gret <- matrix(as.integer(rawToBits(ret)),nrow=8*3,ncol=4)
  expect_equal(gret[1:19,],ix)
  
})

