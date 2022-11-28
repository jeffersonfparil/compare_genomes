args = commandArgs(trailingOnly=TRUE)
# args = c("EPSPS-OG0006242.aln-2.pw.kaks.tmp", "0.001")
# args = c("IRRELEVANT_INPUT_ARGUMENT", "0.001", "TRUE")
fname_input = args[1]
alpha = as.numeric(args[2])
final_plot = tryCatch(as.logical(args[3]), error=function(e){FALSE})
if (is.na(final_plot)){
    final_plot = FALSE
}

PLOT = function(dat, fname_input, alpha, cex=1){
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

    svg(paste0(fname_input, ".svg"), height=5*n, width=8*m)
    par(mfrow=c(n,m), cex=cex)
    for (aln in alignments) {
        # aln = unique(df$Species_alignment2)[1]
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
    write.table(outdf, file=paste0(fname_input, ".sigout"), quote=FALSE, row.names=FALSE)
    return(paste0("OUTPUTS: ", fname_input, ".{svg, sigout}"))
}

if (!final_plot) {
    dat = read.delim(fname_input, header=TRUE)
    PLOT(dat, fname_input)
} else {
    fname_EPSPS = "EPSPS-OG0006242.aln-2.pw.kaks.tmp"
    dat = read.delim(fname_EPSPS, header=TRUE)
    idx = list(Lolium_perenne =         grepl("Lolium_perenne", dat$Sequence),
               Oryza_sativa =           grepl("Oryza_sativa", dat$Sequence),
               Lolium_rigidum =         grepl("XP_047075249.1", dat$Sequence),
               Sorghum_bicolor =        grepl("Sorghum_bicolor", dat$Sequence),
               Arabidopsis_thaliana_1 = grepl("NP_973996.1", dat$Sequence),
               Arabidopsis_thaliana_2 = grepl("NP_001324851.1", dat$Sequence))
    for (i in 1:length(idx)){
        if (!exists("subdat")){
            subdat = dat[idx[[i]], ]
        } else (
            subdat = rbind(subdat, dat[idx[[i]], ])
        )
    }
    PLOT(subdat, "KaKs_ratio_EPSPS_gene", alpha, cex=1.5)
}