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

test_that("can round trip binary matrices",{
  
  ix <- matrix(sample(c(0L,1L),64,replace=T),nrow = 16,ncol = 4)
  tS <- cov(ix)

  


  ret <-   snp2raw(ix)
  
  bS <- covbin(ret)
  
  pcv <- popcnt_v(ret)
  tbdf <- function(x,y){
    popcnt_fun(x&y)
  }
  
  cS <- matrix(0,4,4)
  for(i in 1:4){
    for(j in i:4){
      cS[i,j] <- tbdf(ret[,i],ret[,j])
      cS[j,i] <- cS[i,j]
    }
  }
  
  cS <- (16*cS-outer(pcv,pcv))/(16*15)
  expect_equal(tS,cS)  
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

