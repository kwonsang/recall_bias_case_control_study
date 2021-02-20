##################################################################
## This Rcode is for replication of tables and figures
##################################################################

##################################################################
#### Figure 2
res.mat = read.csv("sensi_res.csv")

#jpeg("fig2.jpg", width = 8, height = 5, unit = "in", res = 300)

## M method: Estimates with CIs
plot(res.mat$zeta.seq, res.mat$ML, lty = 1, type ="l", ylim = c(0.8,3.3), xlab = expression(zeta[0] == zeta[1]), ylab = "Marginal COR")
points(res.mat$zeta.seq[c(1,11,21,31,41,51)], res.mat$ML[c(1,11,21,31,41,51)], pch=4)
lines(res.mat$zeta.seq, res.mat$ML.lci, lty = 2)
points(res.mat$zeta.seq[c(1,11,21,31,41,51)], res.mat$ML.lci[c(1,11,21,31,41,51)], pch=4)
lines(res.mat$zeta.seq, res.mat$ML.uci, lty = 2)
points(res.mat$zeta.seq[c(1,11,21,31,41,51)], res.mat$ML.uci[c(1,11,21,31,41,51)], pch=4)

## Stratification method: Estimates with CIs
lines(res.mat$zeta.seq, res.mat$S, lwd = 2, col = "blue")
points(res.mat$zeta.seq[c(1,11,21,31,41,51)], res.mat$S[c(1,11,21,31,41,51)], pch=20, col = "blue")
lines(res.mat$zeta.seq, res.mat$S.lci, lty = 2, col = "blue")
points(res.mat$zeta.seq[c(1,11,21,31,41,51)], res.mat$S.lci[c(1,11,21,31,41,51)], pch=20, col = "blue")
lines(res.mat$zeta.seq, res.mat$S.uci, lty = 2, col = "blue")
points(res.mat$zeta.seq[c(1,11,21,31,41,51)], res.mat$S.uci[c(1,11,21,31,41,51)], pch=20, col = "blue")

## Drawing the null line (i.e., Marginal COR = 1)
abline(h = 1, lty = 2, col = "red")

## Legend
legend(x=0, y=3.3, c("ML", "S"), bty = "n", lty=1, lwd=2, pch=c(4,20), col=c("black", "blue"))

# dev.off()

#### Table 3
ML.est = res.mat[1,2]
ML.sd = sqrt(res.mat[1,3])
ML.lci = res.mat[1,4]
ML.uci = res.mat[1,5]
round(c(ML.est, ML.sd, ML.lci, ML.uci), 2)

S.est = res.mat[1,6]
S.sd = sqrt(res.mat[1,7])
S.lci = res.mat[1,8]
S.uci = res.mat[1,9]
round(c(S.est, S.sd, S.lci, S.uci), 2)

#### Table 4
tab4 = cbind(res.mat[c(11,21,31,41,51), 1], round(res.mat[c(11,21,31,41,51),c(2,4,5,6,8,9)], 2))
colnames(tab4)[1] = "zeta"
tab4