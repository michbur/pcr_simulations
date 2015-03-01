source("sim_ddpcr_bkm_orig.R")
library(dpcR)

#test 1 -------------------------

test1 <- sim_ddpcr_bkm_orig(0.5,n_exp=1,seed=1)
test1_res <- -log(1-test1[[1]][1]/test1[[1]][2])

set.seed(1)
test1_2 <- sim_ddpcr_bkm(0.5, n_exp=1, type = "tnp")
test1_2_res <- as.vector(-log(1-test1_2/test1_2@n))

test1_res == test1_2_res


#test 2 -------------------------
test2 <- sim_ddpcr_bkm_orig(0.5, n_exp = 1, seed = 2, sddropc = 500, mudropr = 0.7,
                            sddropr = 0.1, Pvar = TRUE, piperr = 0.02, dropsd = 0.2,
                            falpos = 0.001, falneg = 0.01)
str(test2)
test2_res <- -log(1-test2[[1]][1]/test2[[1]][2])

set.seed(2)
test2_2 <- sim_ddpcr_bkm(0.5, n_exp=1, type = "tnp", sddropc = 500, mudropr = 0.7,
                         sddropr = 0.1, Pvar = TRUE, piperr = 0.02, dropsd = 0.2,
                         falpos = 0.001, falneg = 0.01)
test2_2_res <- as.vector(-log(1-test2_2/test2_2@n))
test2_res == test2_2_res
#difference between test2[[1]][1] (4476) and test2_2 (4463)


#test 3 -------------------------
test3 <- sim_ddpcr_bkm_orig(0.5, n_exp = 1, seed = 3, pos_sums = TRUE, fluo = TRUE,
                            sddropc = 500, mudropr = 0.7, sddropr = 0.1, Pvar = TRUE,
                            piperr = 0.02, dropsd = 0.2, falpos = 0.001, falneg = 0.01)
str(test3)
-log(1-sum(test3[[1]])/length(test3[[1]]))



set.seed(3)
test3_2 <- sim_ddpcr_bkm(0.5, n_exp = 1, type = "np",
                         sddropc = 500, mudropr = 0.7, sddropr = 0.1, Pvar = TRUE,
                         piperr = 0.02, dropsd = 0.2, falpos = 0.001, falneg = 0.01)
-log(1-sum(test3_2)/test3_2@n)

#difference between sum(test3[[1]]) (4476) and sum(test3_2) (6048)
plot(density(test3[[2]]))
#add way to convert fluo to np

test3b <- sim_ddpcr_bkm_orig(0.5, n_exp = 1, seed = 3, pos_sums = FALSE, fluo = TRUE,
                        sddropc = 500, mudropr = 0.7, sddropr = 0.1, Pvar = TRUE,
                        piperr = 0.02, dropsd = 0.2, falpos = 0.001, falneg = 0.01)
str(test3b)

# test 4 -------------------------------------
system.time(test4 <- sim_ddpcr_bkm_orig(0.5, n_exp = 1, seed = 4, pos_sums = TRUE, fluo = 10, sddropc = 500,
                                        mudropr = 0.7, sddropr = 0.1, Pvar = TRUE, piperr = 0.02,
                                        dropsd = 0.2, rain = 0.1))
test4_res <- -log(1-sum(test4[[1]])/length(test4[[1]]))


set.seed(4)
system.time(test4_2 <- sim_ddpcr_bkm(0.5, n_exp=1, type = "np", fluo_range = 10, sddropc = 500,
                                     mudropr = 0.7, sddropr = 0.1, Pvar = TRUE, piperr = 0.02,
                                     dropsd = 0.2, rain = 0.1))
#should be around 300 times faster
test4_2_res <- as.vector(-log(1-sum(test4_2)/test4_2@n))
#difference between sum(test4[[1]]) (4932) and sum(test4_2) (5042)