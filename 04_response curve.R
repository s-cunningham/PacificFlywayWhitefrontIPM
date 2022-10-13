rm(list=ls())
library(boot)

load('data/ipm_data_noBBL.RData')
bpop <- dat$popcount
best <- dat$best
burn <- dat$burn2

load('data/ipmfit_noBBL_noCSE.RData')
nhs_mean1 <- c(fit$samples[[1]][,which(row.names(fit$summary)=='nhs_mean[1]')], 
               fit$samples[[2]][,which(row.names(fit$summary)=='nhs_mean[1]')], 
               fit$samples[[3]][,which(row.names(fit$summary)=='nhs_mean[1]')])
nhs_mean2 <- c(fit$samples[[1]][,which(row.names(fit$summary)=='nhs_mean[2]')], 
               fit$samples[[2]][,which(row.names(fit$summary)=='nhs_mean[2]')], 
               fit$samples[[3]][,which(row.names(fit$summary)=='nhs_mean[2]')])
logit_nhs_ddp <- c(fit$samples[[1]][,which(row.names(fit$summary)=='logit_nhs_ddp')], 
                   fit$samples[[2]][,which(row.names(fit$summary)=='logit_nhs_ddp')], 
                   fit$samples[[3]][,which(row.names(fit$summary)=='logit_nhs_ddp')])
logit_nhs_best <- c(fit$samples[[1]][,which(row.names(fit$summary)=='logit_nhs_best')], 
                    fit$samples[[2]][,which(row.names(fit$summary)=='logit_nhs_best')], 
                    fit$samples[[3]][,which(row.names(fit$summary)=='logit_nhs_best')])
ar_mean <- c(fit$samples[[1]][,which(row.names(fit$summary)=='ar_mean')], 
             fit$samples[[2]][,which(row.names(fit$summary)=='ar_mean')], 
             fit$samples[[3]][,which(row.names(fit$summary)=='ar_mean')])
log_ar_ddp <- c(fit$samples[[1]][,which(row.names(fit$summary)=='log_ar_ddp')], 
                fit$samples[[2]][,which(row.names(fit$summary)=='log_ar_ddp')], 
                fit$samples[[3]][,which(row.names(fit$summary)=='log_ar_ddp')])

npred <- 100

range(bpop[which(burn==1)], na.rm=T)
range(bpop[which(burn==2)], na.rm=T)
range(best[which(burn==1)], na.rm=T)
range(best[which(burn==2)], na.rm=T)

bpop_pred1 <- seq(80, 420, length.out=npred)
bpop_pred2 <- seq(280, 840, length.out=npred)
best_pred1 <- seq(-1.4, 3, length.out=npred)
best_pred2 <- seq(-3, 1.8, length.out=npred)
bpop_preda <- seq(80, 840, length.out=npred)

logn_mean <- mean(log(bpop),na.rm=T)
logn_sd <- sd(log(bpop),na.rm=T)

nhs_bpop1 <- nhs_bpop2 <- nhs_best1 <- nhs_best2 <- ar_bpop <- matrix(, length(ar_mean), npred)

for (i in 1:npred) {
  nhs_bpop1[,i] <- inv.logit(
    logit(nhs_mean1) + 
    logit_nhs_ddp * (log(bpop_pred1[i]) - logn_mean) / logn_sd)
  nhs_bpop2[,i] <- inv.logit(
    logit(nhs_mean2) + 
    logit_nhs_ddp * (log(bpop_pred2[i]) - logn_mean) / logn_sd)
  nhs_best1[,i] <- inv.logit(
    logit(nhs_mean1) + 
    logit_nhs_best * best_pred1[i])
  nhs_best2[,i] <- inv.logit(
    logit(nhs_mean2) + 
    logit_nhs_best * best_pred1[i])
  ar_bpop[,i] <- exp(
    log(ar_mean) + 
    log_ar_ddp * (log(bpop_preda[i]) - logn_mean) / logn_sd)
}

nhs_bpop_qt1 <- apply(nhs_bpop1, 2, quantile, probs=c(.5, .025, .975))
nhs_bpop_qt2 <- apply(nhs_bpop2, 2, quantile, probs=c(.5, .025, .975))
nhs_best_qt1 <- apply(nhs_best1, 2, quantile, probs=c(.5, .025, .975))
nhs_best_qt2 <- apply(nhs_best2, 2, quantile, probs=c(.5, .025, .975))
ar_bpop_qt <- apply(ar_bpop, 2, quantile, probs=c(.5, .025, .975))


par(mfrow=c(1,3))
par(mar=c(5,5,2,1))
par(oma=c(0,0,0,0))

# NHS x ddp
plot(nhs_bpop_qt1[1,] ~ bpop_pred1, type='n', 
     xlim=c(80,840), ylim=c(.4,1), axes=F, xlab='', ylab='')
for (t in 1:(length(bpop)-1)) {
  lines(x=c(bpop[t],bpop[t]), y=c(fit$q2.5$nhs[t],fit$q97.5$nhs[t]), 
        lwd=1, col='grey56')
}
points(fit$mean$nhs ~ bpop[-length(bpop)], cex=1.6, col='grey56', 
       pch=ifelse(burn==1, 4, 19))
lines(nhs_bpop_qt1[1,] ~ bpop_pred1, lwd=3, lty=1, col='grey16')
lines(nhs_bpop_qt1[2,] ~ bpop_pred1, lwd=2, lty=1, col='grey16')
lines(nhs_bpop_qt1[3,] ~ bpop_pred1, lwd=2, lty=1, col='grey16')
lines(nhs_bpop_qt2[1,] ~ bpop_pred2, lwd=3, lty=2, col='grey16')
lines(nhs_bpop_qt2[2,] ~ bpop_pred2, lwd=2, lty=2, col='grey16')
lines(nhs_bpop_qt2[3,] ~ bpop_pred2, lwd=2, lty=2, col='grey16')
axis(1, at=seq(100, 800, 100), cex.axis=1.4)
axis(2, at=seq(.4, 1, .2), las=2, cex.axis=1.4)
title(xlab='Population Size (x 1000)', cex.lab=1.6, line=2.8)
title(ylab='Non-harvest Survival', cex.lab=1.8, line=2.8)
mtext(text='(a)', side=3, outer=F, line=-1, adj=0, cex=1.8)
legend('bottomleft', col=c('grey56','grey16','grey56','grey16'), 
       pch=c(4,NA,19,NA), lty=c(NA,1,NA,2), bty='n', cex=1.4, lwd=1.8, 
       legend=c('IPM estimates, burn-allowed',
                'Fitted line, burn-allowed',
                'IPM estimates, burn-restricted',
                'Fitted line, burn-restricted'))

# NHS x BEST
plot(nhs_best_qt1[1,] ~ best_pred1, type='n', 
     xlim=c(-3,3), ylim=c(.4,1), axes=F, xlab='', ylab='')
for (t in 1:(length(bpop)-1)) {
  lines(x=c(best[t],best[t]), y=c(fit$q2.5$nhs[t],fit$q97.5$nhs[t]), 
        lwd=1, col='grey56')
}
points(fit$mean$nhs ~ best, cex=1.4, col='grey56', 
       pch=ifelse(burn==1, 4, 19))
lines(nhs_best_qt1[1,] ~ best_pred1, lwd=3, lty=1, col='grey16')
lines(nhs_best_qt1[2,] ~ best_pred1, lwd=2, lty=1, col='grey16')
lines(nhs_best_qt1[3,] ~ best_pred1, lwd=2, lty=1, col='grey16')
lines(nhs_best_qt2[1,] ~ best_pred2, lwd=3, lty=2, col='grey16')
lines(nhs_best_qt2[2,] ~ best_pred2, lwd=2, lty=2, col='grey16')
lines(nhs_best_qt2[3,] ~ best_pred2, lwd=2, lty=2, col='grey16')
axis(1, at=seq(-3, 3, 1), cex.axis=1.4)
axis(2, at=seq(.4, 1, .2), las=2, cex.axis=1.4)
title(xlab=expression(paste('El Ni', tilde(n), 'o-Southern Oscillation')), cex.lab=1.6, line=2.8)
title(ylab='Non-harvest Survival', cex.lab=1.8, line=2.8)
mtext(text='(b)', side=3, outer=F, line=-1, adj=0, cex=1.8)
legend('bottomleft', col=c('grey56','grey16','grey56','grey16'), 
       pch=c(4,NA,19,NA), lty=c(NA,1,NA,2), bty='n', cex=1.4, lwd=1.8, 
       legend=c('IPM estimates, burn-allowed',
                'Fitted line, burn-allowed',
                'IPM estimates, burn-restricted',
                'Fitted line, burn-restricted'))

# age ratio x ddp
plot(ar_bpop_qt[1,] ~ bpop_preda, type='n', 
     xlim=c(80,840), ylim=c(0,1), axes=F, xlab='', ylab='')
for (t in 1:(length(bpop)-1)) {
  lines(x=c(bpop[t],bpop[t]), y=c(fit$q2.5$ar[t],fit$q97.5$ar[t]), 
        lwd=1, col='grey56')
}
points(fit$mean$ar ~ bpop[-length(bpop)], cex=1.4, col='grey56', pch=15)
lines(ar_bpop_qt[1,] ~ bpop_preda, lwd=3, lty=1, col='grey16')
lines(ar_bpop_qt[2,] ~ bpop_preda, lwd=2, lty=1, col='grey16')
lines(ar_bpop_qt[3,] ~ bpop_preda, lwd=2, lty=1, col='grey16')
axis(1, at=seq(100, 800, 100), cex.axis=1.4)
axis(2, at=seq(0, 1, .2), las=2, cex.axis=1.4)
title(xlab='Population Size (x 1000)', cex.lab=1.6, line=2.8)
title(ylab='Age Ratio', cex.lab=1.8, line=2.8)
mtext(text='(c)', side=3, outer=F, line=-1, adj=0, cex=1.8)




