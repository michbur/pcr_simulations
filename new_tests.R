library(dpcR)


#test 1 -------------------------

set.seed(1)
test1_2 <- sim_ddpcr_bkm(0.5, n_exp=1, type = "tnp")
test1_2_res <- as.vector(-log(1-test1_2/test1_2@n))

test1_res == test1_2_res

#test 2 -------------------------

set.seed(2)
test2_2 <- sim_ddpcr_bkm(0.5, n_exp=1, type = "tnp", sddropc = 500, mudropr = 0.7,
                         sddropr = 0.1, Pvar = TRUE, piperr = 0.02, dropsd = 0.2,
                         falpos = 0.001, falneg = 0.01)
test2_2_res <- as.vector(-log(1-test2_2/test2_2@n))
test2_res == test2_2_res
#difference between test2[[1]][1] (4476) and test2_2 (4463)

#test 3 -------------------------

set.seed(3)
test3_2 <- sim_ddpcr_bkm(0.5, n_exp = 1, type = "np",
                         sddropc = 500, mudropr = 0.7, sddropr = 0.1, Pvar = TRUE,
                         piperr = 0.02, dropsd = 0.2, falpos = 0.001, falneg = 0.01)
-log(1-sum(test3_2)/test3_2@n)

#difference between sum(test3[[1]]) (4476) and sum(test3_2) (6048)
plot(density(test3[[2]]))
#add way to convert fluo to np

#test 4 -------------------------
set.seed(4)
system.time(test4_2 <- sim_ddpcr_bkm(0.5, n_exp=1, type = "np", fluo_range = 10, sddropc = 500,
                                     mudropr = 0.7, sddropr = 0.1, Pvar = TRUE, piperr = 0.02,
                                     dropsd = 0.2, rain = 0.1))
test4_2_res <- as.vector(-log(1-sum(test4_2)/test4_2@n))
#difference between sum(test4[[1]]) (4932) and sum(test4_2) (5042)


# test 5 -------------------------------------

set.seed(5)
system.time(test5_2 <- sim_ddpcr_bkm(0.5, n_exp = 8, type = "np",
                                     sddropc = 500, mudropr = 0.7, sddropr = 0.1, Pvar = TRUE,
                                     piperr = 0.02, dropsd = 0.2, falpos = 0.001, falneg = 0.01))

# test 6 -------------------------------------
set.seed(6)
system.time(test6_2 <- sim_ddpcr_bkm(exp(-4:1), n_exp = 8, type = "np",
                                     sddropc = 500, mudropr = 0.7, sddropr = 0.1, Pvar = TRUE,
                                     piperr = 0.02, dropsd = 0.2, falpos = 0.001, falneg = 0.01))

conc_2 <- summary(test6_2, print = FALSE)

conc_2 <- matrix(-log(1-conc_2$partitions$k/conc_2$partitions$n), nrow = 8) 

