source("sim_ddpcr_bkm_orig.R")

#test 1 -------------------------

test1 <- sim_ddpcr_bkm_orig(0.5,n_exp=1,seed=1)
test1_res <- -log(1-test1[[1]][1]/test1[[1]][2])


#test 2 -------------------------
test2 <- sim_ddpcr_bkm_orig(0.5, n_exp = 1, seed = 2, sddropc = 500, mudropr = 0.7,
                            sddropr = 0.1, Pvar = TRUE, piperr = 0.02, dropsd = 0.2,
                            falpos = 0.001, falneg = 0.01)
str(test2)
test2_res <- -log(1-test2[[1]][1]/test2[[1]][2])

#test 3 -------------------------
test3 <- sim_ddpcr_bkm_orig(0.5, n_exp = 1, seed = 3, pos_sums = TRUE, fluo = TRUE,
                            sddropc = 500, mudropr = 0.7, sddropr = 0.1, Pvar = TRUE,
                            piperr = 0.02, dropsd = 0.2, falpos = 0.001, falneg = 0.01)
str(test3)
-log(1-sum(test3[[1]])/length(test3[[1]]))


test3b <- sim_ddpcr_bkm_orig(0.5, n_exp = 1, seed = 3, pos_sums = FALSE, fluo = TRUE,
                             sddropc = 500, mudropr = 0.7, sddropr = 0.1, Pvar = TRUE,
                             piperr = 0.02, dropsd = 0.2, falpos = 0.001, falneg = 0.01)
str(test3b)


# test 4 -------------------------------------
system.time(test4 <- sim_ddpcr_bkm_orig(0.5, n_exp = 1, seed = 4, pos_sums = TRUE, fluo = 10, sddropc = 500,
                                        mudropr = 0.7, sddropr = 0.1, Pvar = TRUE, piperr = 0.02,
                                        dropsd = 0.2, rain = 0.1))
test4_res <- -log(1-sum(test4[[1]])/length(test4[[1]]))


# test 5 -------------------------------------
system.time(test5 <- sim_ddpcr_bkm_orig(0.5, n_exp = 8, seed = 5, pos_sums = TRUE, fluo = TRUE,
                                        sddropc = 500, mudropr = 0.7, sddropr = 0.1, Pvar = TRUE,
                                        piperr = 0.02, dropsd = 0.2, falpos = 0.001, falneg = 0.01))
str(test5)

# test 6 -------------------------------------
system.time(test6 <- sim_ddpcr_bkm_orig(exp(-4:1), n_exp = 8, seed = 6, pos_sums = TRUE, fluo = TRUE,
                                        sddropc = 500, mudropr = 0.7, sddropr = 0.1, Pvar = TRUE,
                                        piperr = 0.02, dropsd = 0.2, falpos = 0.001, falneg = 0.01))
conc <- NULL
for(j in 1:6){
  conct <- NULL
  for(i in 1:8){
    conct <- c(conct, -log(1-sum(test6[[16*(j-1)+2*i-1]])/length(test6[[16*(j-1)+2*i-1]])))
  }
  conc <- cbind(conc, conct)
}
colnames(conc) <- round(exp(-4:1), 3)
conc
