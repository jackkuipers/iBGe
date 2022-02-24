for (seedy in 101:200) {
  for (kk in 1:6) {
   filecommand<-paste("\nkk <- ",kk,"\nseed_number <- ",seedy,"\nsource(\"iBGe_perfect.R\")\n",sep="")
      write(filecommand, file = paste0("runeuler", kk, "seed", seedy, ".R"))
  }
}
