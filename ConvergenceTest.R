library(coda)

setwd("~/Desktop/Master_v2/conv")
posterior_a
post<-posterior_a$simul
post<-posterior_a

main_colour <- "navy"
main_colour2<- "deeppink3"
variable<-'UE_7'
country<-'USA'

mh.draws <- mcmc(post)
a<-geweke.diag(mh.draws)

for (i in a$z){
  n<-n+1
  ifelse(abs(i) > 1.96,print(a$z[n]),"FALSE")}

table(abs(a$z) > 1.96)
df<-as.data.frame(cbind(a$z))
png(file = paste0("Z_plot_",country,"_",variable,".png"), width = 10, height = 5, units ="in",res=300)
barplot(df[,1], main="Z-score for parameters", 
        xlab="parameter")
abline(h=c(-1.96,1.96), col=main_colour2)
dev.off()

n<-0
s<-0.5
colnames(post)
draw_conv<-list()
for (i in 1:ncol(post)){  
  mc <- mcmc(post[,i])
  if (abs(geweke.diag(mc)$z) > 1.96){ 
    print(paste(colnames(post)[i],geweke.diag(mc)$z))
    n<-n+1
    png(file = paste0("geweke_plot_",country,variable,"_",colnames(post)[i],".png"), width = 5, height = 3, units ="in",res=300)
    g<-geweke.plot(mc, frac1 = 0.1, frac2 = 0.5, nbins=40, pvalue=0.05,
                   cex.lab=s,cex.axis=s,cex=s,main="",
                   mar=c(2, 1, 1, 0) + 0.1,mgp=c(2, s, 0),
                   omi=c(0, 0, 0, 0))
    #c<-autocorr.plot(mc)
    #par(mfrow=c(5,5))
    #plot(g)
    #plot(c,add-TRUE)
    dev.off()
  }
}

geweke.diag(mcmc(posterior_a$chain))
summary(mcmc(posterior_a$chain))

mh.draws<-mcmc(posterior_a$chain)

png(file = paste0("trace_rho_",country,variable,".png"), width = 12, height = 5, units ="in",res=300)
plot( posterior_a$chain, type="s", xpd=NA, ylab="Parameter", xlab="Sample", las=1)
dev.off()
hist(posterior_a$chain, 50, freq=FALSE, main="", las=1,
     xlab="x", ylab="Probability density")

png(file = paste0("geweke_plot_rho_",country,variable,".png"), width = 6, height = 5, units ="in",res=300)
geweke.plot(mh.draws, frac1 = 0.1, frac2 = 0.5,
            nbins=40, pvalue=0.05)
dev.off()

png(file = paste0("autocorr_plot_rho_",country,variable,".png"), width = 6, height = 5, units ="in",res=300)
autocorr.plot(mh.draws)
dev.off()

conv_test<-summary(mh.draws)



n<-0
draw_geweke <- function(f="geweke.plot"){}
draw_conv<-list()
not_conv <-c()
for (i in 1:ncol(posterior_a)){  
  mc <- mcmc(posterior_a[,i])
  if (geweke.diag(mc)$z > 1.96){
    n<-n+1
    not_conv[n]<-i
    print(not_conv)
    print(geweke.diag(mc))}}
    
  temp_fun<-function(f,i)  {
    mc <- mcmc(posterior_a[,i])
    ifelse(f=="geweke.plot",
           geweke.plot(mc, frac1 = 0.1, frac2 = 0.5,nbins=40, pvalue=0.05),
          autocorr.plot(mc))
        }
    
    png(file = paste0(f,"_",country,variable,"_",i,".png"), width = 8.27, height = 11.69, units ="in",res=300)
    plots = lapply(not_conv, function(.x) temp_fun(f,.x))
    do.call(grid.arrange,c(plots, ncol=3))
    dev.off()
    
    #par(mfrow=c(5,5))
    #plot(g)
    #plot(c,add-TRUE)
    dev.off()
  }
}

conv_test<-summary(mh.draws)




library(ggmcmc)
S <- ggs(mh.draws)
print(ggs_geweke(S))
ggs_autocorrelation(S)
rows<-4
columns<-3
png(file = paste0("geweke_plot_",country,"_",variable,".png"), width = 8.27, height = 11.69, units ="in",res=300)
do.call(grid.arrange,c(draw_conv, ncol=3))
dev.off()

p<-marrangeGrob(draw_conv, nrow=rows, ncol=columns)
print(p)
dev.off()



