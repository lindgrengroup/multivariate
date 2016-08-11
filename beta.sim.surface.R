require(MASS)
require(foreach)
require(doRNG)
require(doMC)

registerDoMC(cores=16)

set.seed(14)

rho <- 0.5

n.sim <- 100
n.samples <- 100000

thresh <- 5e-8

#interval <- c(-0.03,-0.025,-0.02,-0.015,-0.01,-0.005,0,0.005,0.01,0.015,0.02,0.025,0.03)                                                                                                                                                                                                                                                                                 

interval <- seq(-0.03,0.03,by=0.001)

mat<- foreach(beta.A=interval,.combine=cbind) %do% {
  row <- foreach(beta.B=interval,.combine=c) %do% {

    name <- paste(beta.A,"_",beta.B,sep="")

    print(name)
    pleiotropy <- foreach(b.sim=1:n.sim) %dorng% {

      snp <- rbinom(n.samples,2,0.5)

      mean.A <- beta.A * snp
      mean.B <- beta.B * snp

      residual.variation <- mvrnorm(10000,c(0,0),Sigma=matrix(c(1,rho,rho,1),nrow=2))

      traits <- cbind((mean.A + residual.variation[,1]), (mean.B + residual.variation[,2]))

      traits <- cbind(traits,residuals(lm(traits[,1]~traits[,2])),residuals(lm(traits[,2]~traits[,1])))

      sample.names <- paste("ID",1:n.samples,sep="")
      sample.f <- cbind(sample.names,sample.names,rep(0,n.samples),traits)

      sample.f <- rbind(c(0, 0, 0, "P", "P", "P", "P"),sample.f)

      colnames(sample.f) <- c("ID_1", "ID_2", "missing", "TraitA", "TraitB", "TraitAadjTraitB", "TraitBadjTraitA")

      gen.f <- c(paste(name,"-SNP-",b.sim,sep=""), "rs1",format(100*b.sim,scientific=FALSE),"A","C",unlist(list(c(1,0,0),c(0,1,0),c(0,0,1))[snp+1]))

      write.table(matrix(gen.f,nrow=1),file=paste("in/in", name,"-",b.sim,".gen",sep=""),col.names=F, quote=F, row.names=F, sep=" ")

      write.table(sample.f,file=paste("in/in", name,"-",b.sim,".sample",sep=""),col.names=T, quote=F, row.names=F, sep=" ")

      #system(paste("/well/lindgren/alexd/UKBIOBANK/MULTIVAR/BIN/PLEIOTROPY --print_complex --print_covariance --betas --pheno_name TraitA --pheno_name TraitB -g in/in",name,"-", b.sim,".gen -s in/in",name,"-",b.sim,".sample -o out/out",name,"-",b.sim, sep=""),ignore.stdout=TRUE)                                                                                  
      system(paste("/apps/well/snptest/2.5.1/snptest_v2.5.1 -data in/in",name,"-", b.sim,".gen in/in",name,"-",b.sim,".sample -frequentist 1 -method expected -pheno TraitAadjTraitB -o out/out",name,"-",b.sim,".TraitAadjTraitB -use_raw_phenotypes",sep=""),ignore.stdout=TRUE)

      ##t <- read.table(paste("out/out",name,"-",b.sim,".result",sep=""), header=T)                                                                                                                                                                                                                                                                                       
      ##t <- read.table(paste("out/out",name,"-",b.sim,".result",sep=""), header=T)                                                                                                                                                                                                                                                                                       
      t <- read.table(paste("out/out",name,"-",b.sim,".TraitAadjTraitB",sep=""), header=T)

      file.remove(paste("in/in", name,"-",b.sim,".sample",sep=""))
      file.remove(paste("in/in", name,"-",b.sim,".gen",sep=""))
      ##file.remove(paste("out/out",name,"-",b.sim,".result",sep=""))                                                                                                                                                                                                                                                                                                     
      file.remove(paste("out/out",name,"-",b.sim,".TraitAadjTraitB",sep=""))

      ##t$"P.value"                                                                                                                                                                                                                                                                                                                                                       
      t$"frequentist_add_pvalue"
    }

    res <- sum(pleiotropy<thresh)
    print(res)
    res
  }
  row
}

save(mat,file="mat1000_rho05_AadjB_large.RData")

load(file="mat100_rho05_AadjB.RData")
mat.adj <- mat
x <- interval
cols <- colorRampPalette(c("orange","red"))(7)
png("~/temp/lines_both.png")
plot(x,mat.adj[,7],type="l",col=cols[1],ylim=c(0,100),xlab="betaB",ylab="Power",main="TraitAadjTraitB")
for(i in 8:13){lines(x,mat.adj[,i],col=cols[i-6])}
legend("topright",legend=x[7:13],col=cols,lty=1,title="betaA")
load(file="mat100_rho05.RData")
mat.scopa <- mat
cols <- colorRampPalette(c("lightblue","blue"))(7)
#png("~/temp/lines_05_scopa.png")                                                                                                                                                                                                                                                                                                                                         
lines(x,mat.scopa[,7],col=cols[1],ylim=c(0,100),xlab="betaB",ylab="Power",main="SCOPA TraitA+TraitB")
for(i in 8:13){lines(x,mat.scopa[,i],col=cols[i-6])}
legend("bottomleft",legend=x[7:13],col=cols,lty=1,title="betaA")
dev.off()

