load("obsdata_s.Rdata")
cores <- 40
cl <- makeCluster(cores)
registerDoParallel(cl, cores=cores)
M.eval <- 600
success_count <- 0
start <- Sys.time()

result.diffAll <- foreach(i=1:M.eval, .combine='c', .packages=c("caret", "SuperLearner", "glmnet", "survSuperLearner", "CFsurvival")) %dopar% {
  if (success_count>=500) {
    return(NULL) 
  }
  dat.i <- dat.diffAll[[i]]
  tryCatch(
    {
      result <- list(FedCSA(data=dat.i, s=100*i))
      success_count <<- success_count+1
      result  
    },
    error=function(e) {
      NULL  
    }
  )
}
save(file="Res_diffAll_s.Rdata", result.diffAll)

end <- Sys.time()
end - start

stopCluster(cl)
