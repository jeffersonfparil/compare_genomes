args = commandArgs(trailingOnly=TRUE)
# args = c("GlutathioneStransferase-OG0005237.aln-4.pw.kaks.tmp", "0.001")
fname_input = args[1]
alpha = as.numeric(args[2])
cex=1
fname_output_svg = gsub(".tmp", ".svg", fname_input)
fname_output_sig = gsub(".tmp", "-SIGNIFICANT_PEAKS.csv", fname_input)

dat = read.delim(fname_input, header=TRUE)
dat$Sequence = gsub("_rna-", "_rna_", dat$Sequence)
X = matrix(unlist(strsplit(as.character(dat$Sequence), "-")), ncol=12, byrow=TRUE)
Y = matrix(unlist(strsplit(X[,11], "\\(")), ncol=2, byrow=TRUE)
df = data.frame(
    Gene=X[,2],
    Orthogroup=X[,3],
    Species1=X[,1],
    Species_alignment1=X[,4],
    Species2=X[,8],
    Species_alignment2=Y[,1],
    Position_ini=as.numeric(Y[,2]),
    Position_fin=as.numeric(gsub("\\)", "", X[,12])),
    KaKs=dat$Ka.Ks,
    p=dat$P.Value.Fisher.
)
df$KaKs[is.na(df$KaKs)] = 1.0

alignments = unique(df$Species_alignment2)
n = ceiling(sqrt(length(alignments)))
m = ceiling(length(alignments)/n)

svg(fname_output_svg, height=5*n, width=8*m)
par(mfrow=c(n,m), cex=cex)
for (aln in alignments) {
    # aln = alignments[1]
    subdf = df[df$Species_alignment2 == aln, ]
    subdf = subdf[order(subdf$Position_ini, decreasing=FALSE), ]
    plot(subdf$Position_ini, subdf$KaKs, type="l",
        main=paste0(subdf$Species1[1], "::", subdf$Species_alignment1[1], " vs\n", subdf$Species2[1], "::", subdf$Species_alignment2[1]),
        sub=paste0(subdf$Gene[1], " (", subdf$Orthogroup[1], ")"),
        xlab="Position (bp)", ylab="Ka/Ks")
    grid()
    idx = c(1:nrow(subdf))[!is.na(subdf$KaKs) & (subdf$p <= alpha) & (subdf$KaKs > 1)]
    for (j in idx){
        text(x=subdf$Position_ini[j], y=subdf$KaKs[j], lab="*", col="red", cex=2)
    }
    abline(h=1, lty=2, col="red")
    if (!exists("outdf")){
        outdf = subdf[idx, ]
    } else {
        outdf = rbind(outdf, subdf[idx, ])
    }
}
dev.off()
write.table(outdf, file=fname_output_sig, sep=",", quote=FALSE, row.names=FALSE)
return(paste0("OUTPUTS: ", fname_input, ".{svg, sigout}"))
}
