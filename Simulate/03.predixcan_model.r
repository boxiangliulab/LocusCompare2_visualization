library(glmnet)
args = commandArgs(trailingOnly=T)

infile = args[1]
gene = args[2]

d = read.table(infile)

y = as.numeric(d[1,4:dim(d)[2]])
G = t(as.matrix(d[2:dim(d)[1], 4:dim(d)[2]]))

G[is.na(G)] = 0

SNP = d[2:dim(d)[1], 2]


enet.fit = cv.glmnet(G,y)
wts = as.numeric(coef(enet.fit))
wts = wts[2:length(wts)]

index = which(wts!=0)

outd = cbind(SNP[index], wts[index])


outfile = paste0("sim_rst/predixcan_wts/",gene,".enet.wts")

write(file=outfile, t(outd), ncol=2)


