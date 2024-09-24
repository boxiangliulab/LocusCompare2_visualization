## 1. sim eqtl

args = commandArgs(trailingOnly=TRUE)

phi = 0.6
neqtl = 2



infile = args[1]
gene = args[2]
d = read.table(infile)
X = d[,4:503]
X[is.na(X)] = 0
freq = apply(X,1, function(x) mean(x)/2)
index = which(freq > 0.05 & freq <0.95)
index = intersect(1:500, index)
eqtl_index = sample(index, neqtl, replace = F)
beta = rep(0,1500)
beta[eqtl_index] = rnorm(neqtl, sd=phi)

ye = t(X)%*%beta + rnorm(500)

r2 = summary(lm(ye~t(X[eqtl_index,])))$r.squared



pheno = c("pheno", gene, "expr", ye)
outd = rbind(pheno, as.matrix(d))
outt = cbind(rep(gene,2), d$V2[eqtl_index], beta[eqtl_index], rep(r2,2))
if(!dir.exists("sim_data/eqtl_sbams/")){
  dir.create("sim_data/eqtl_sbams/", showWarnings = FALSE)
}
if(!dir.exists("sim_data/eqtl_truth/")){
  dir.create("sim_data/eqtl_truth/", showWarnings = FALSE)
}
out_data_file = paste0("sim_data/eqtl_sbams/",gene,".expr.sbams.dat")
out_truth_file = paste0("sim_data/eqtl_truth/", gene, ".eqtl.truth")
write(file = out_data_file, t(outd), ncol=503)
write(file = out_truth_file, t(outt), ncol=4)


## 2. sim gwas

args = commandArgs(trailingOnly=TRUE)


infile = args[1]
eqtlfile = args[2]
gene = args[3]

d = read.table(infile)


a0 = -1.734601
a1 = 0

# additional gwas hits
ngwas = 1 

# effect size 
phi = 0.6


annot = rbinom(1, 1, 0.4)

pc = c( exp(a0)/(1+exp(a0)), exp(a0+a1)/(1+exp(a0+a1)) )

causal = rbinom(1,1, pc[annot+1])



# experssion
E = as.numeric(d[1,4:503])

Geno = d[2:1501,]
X = Geno[,4:503]
X[is.na(X)] = 0

freq = as.numeric(apply(X,1, function(x) mean(x)/2))
index = which(freq > 0.05 & freq <0.95)
index = intersect(1001:1500, index)

t = read.table(eqtlfile)
eqtl_snp = t$V2


ldxfile = paste0("sim_data/eqtl_truth/",gene,".gwas_ld_exclude.list")
ldx = try(read.table(ldxfile), silent=T);
if( inherits(ldx, "try-error") == FALSE){
    ld_snp = ldx$V1
    eqtl_snp = union(ld_snp, eqtl_snp);
}

eqtl_index = sapply(eqtl_snp, function(x) which(Geno$V2 == x))

gwas_index = sample(setdiff(index, eqtl_index),1)

beta = rep(0,1500)
beta[gwas_index] = rnorm(ngwas, sd=phi)
alpha = rnorm(1, sd=phi)
if(causal == 0){
    alpha = 0
}

yc = t(X)%*%beta + alpha*E + rnorm(500)

r2 = summary(lm(yc~t(X[c(gwas_index,eqtl_index),])))$r.squared



pheno = c("pheno", gene, "gwas", yc)
outd = rbind(pheno, as.matrix(Geno))

outt = c(gene, annot, alpha, causal, beta[gwas_index], Geno$V2[gwas_index], r2)

if(!dir.exists("sim_data/gwas_sbams/")){
  dir.create("sim_data/gwas_sbams/", showWarnings = FALSE)
}
if(!dir.exists("sim_data/gwas_truth/")){
  dir.create("sim_data/gwas_truth/", showWarnings = FALSE)
}
out_data_file = paste0("sim_data/gwas_sbams/",gene,".gwas.sbams.dat")
out_truth_file = paste0("sim_data/gwas_truth/", gene, ".gwas.truth")

write(file = out_data_file, t(outd), ncol=503)
write(file = out_truth_file, t(outt), ncol=7)


