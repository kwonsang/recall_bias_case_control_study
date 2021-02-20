##### 
# sim.data 
cor.cor.n800.beta00 = read.csv("sim1_n800_beta00.csv")
cor.cor.n800.beta05 = read.csv("sim1_n800_beta05.csv")
cor.cor.n800.beta10 = read.csv("sim1_n800_beta10.csv")
cor.cor.n2000.beta00 = read.csv("sim1_n2000_beta00.csv")
cor.cor.n2000.beta05 = read.csv("sim1_n2000_beta05.csv")
cor.cor.n2000.beta10 = read.csv("sim1_n2000_beta10.csv")

mis.cor.n800.beta00 = read.csv("sim2_n800_beta00.csv")
mis.cor.n800.beta05 = read.csv("sim2_n800_beta05.csv")
mis.cor.n800.beta10 = read.csv("sim2_n800_beta10.csv")
mis.cor.n2000.beta00 = read.csv("sim2_n2000_beta00.csv")
mis.cor.n2000.beta05 = read.csv("sim2_n2000_beta05.csv")
mis.cor.n2000.beta10 = read.csv("sim2_n2000_beta10.csv")

cor.mis.n800.beta00 = read.csv("sim3_n800_beta00.csv")
cor.mis.n800.beta05 = read.csv("sim3_n800_beta05.csv")
cor.mis.n800.beta10 = read.csv("sim3_n800_beta10.csv")
cor.mis.n2000.beta00 = read.csv("sim3_n2000_beta00.csv")
cor.mis.n2000.beta05 = read.csv("sim3_n2000_beta05.csv")
cor.mis.n2000.beta10 = read.csv("sim3_n2000_beta10.csv")

mis.mis.n800.beta00 = read.csv("sim4_n800_beta00.csv")
mis.mis.n800.beta05 = read.csv("sim4_n800_beta05.csv")
mis.mis.n800.beta10 = read.csv("sim4_n800_beta10.csv")
mis.mis.n2000.beta00 = read.csv("sim4_n2000_beta00.csv")
mis.mis.n2000.beta05 = read.csv("sim4_n2000_beta05.csv")
mis.mis.n2000.beta10 = read.csv("sim4_n2000_beta10.csv")

#############
log.mean.mat = rbind(apply(log(cor.cor.n800.beta00), 2, mean),
                     apply(log(cor.cor.n800.beta05), 2, mean),
                     apply(log(cor.cor.n800.beta10), 2, mean),
                     apply(log(cor.cor.n2000.beta00), 2, mean),
                     apply(log(cor.cor.n2000.beta05), 2, mean),
                     apply(log(cor.cor.n2000.beta10), 2, mean),
                     apply(log(mis.cor.n800.beta00), 2, mean),
                     apply(log(mis.cor.n800.beta05), 2, mean),
                     apply(log(mis.cor.n800.beta10), 2, mean),
                     apply(log(mis.cor.n2000.beta00), 2, mean),
                     apply(log(mis.cor.n2000.beta05), 2, mean),
                     apply(log(mis.cor.n2000.beta10), 2, mean),
                     apply(log(cor.mis.n800.beta00), 2, mean),
                     apply(log(cor.mis.n800.beta05), 2, mean),
                     apply(log(cor.mis.n800.beta10), 2, mean),
                     apply(log(cor.mis.n2000.beta00), 2, mean),
                     apply(log(cor.mis.n2000.beta05), 2, mean),
                     apply(log(cor.mis.n2000.beta10), 2, mean),
                     apply(log(mis.mis.n800.beta00), 2, mean),
                     apply(log(mis.mis.n800.beta05), 2, mean),
                     apply(log(mis.mis.n800.beta10), 2, mean),
                     apply(log(mis.mis.n2000.beta00), 2, mean),
                     apply(log(mis.mis.n2000.beta05), 2, mean),
                     apply(log(mis.mis.n2000.beta10), 2, mean))

colnames(log.mean.mat) = c("True", "Crude", "ML", "S_prop", "S_prog", "MH")
round(log.mean.mat, 3)




