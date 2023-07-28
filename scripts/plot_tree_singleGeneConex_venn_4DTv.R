args = commandArgs(trailingOnly=TRUE)
# args = c("ORTHOGROUPS_SINGLE_GENE.NT.timetree.nex",
#          "CONTRACTION_EXPANSION.txt",
#          "PROTEOMES/orthogroups_summarised_gene_counts.csv",
#          "PROTEOMES/orthogroups_gene_counts_families_go.out",
#          "/data-weedomics-1/COMPARE_GENOMES_TEST",
#          ".4DTv",
#          "ORTHOGROUPS_SINGLE_GENE.NT.4DTv",
#          "/data-weedomics-1/compare_genomes/config/comparisons_4DTv.txt",
#          "/data-weedomics-1/compare_genomes/config/venn_species_max_5.txt",
#          "test.svg")
fname_tree = args[1]
fname_conex = args[2]
fname_gene_groups = args[3]
fname_gene_counts = args[4]
dir_name_4DTv = args[5]
extension_name_4DTv = args[6]
fname_4DTv_singlecopy = args[7]
fname_comparisons = args[8]
fname_venn_species_max_5 = args[9]
fname_svg_output = args[10]

library(ape)
library(VennDiagram)
library(png)
library(grid)

svg(fname_svg_output, width=12, height=9)

pushViewport(plotViewport(layout=grid.layout(2, 2), margins=c(0,2,5,5)))
layout(matrix(c(rep(1,times=6), rep(2,times=2), rep(3,times=4), rep(4,times=2),
                rep(5,times=7), rep(6,times=7)), byrow=TRUE, nrow=2))

### Tree: dendrogram
tree = read.nexus(fname_tree)
species_order = tree$tip.label
# tree = ape::rotate(tree, 9)
tree$edge[,1] = rev(tree$edge[,1])
tree$edge[,2] = rev(tree$edge[,2])
tree$edge.length = rev(tree$edge.length)
tree$node.label = rev(tree$node.label)

par(mar=c(0,0,0,0))
plot(x=c(0,1), y=c(0,1), type="n", xaxt="n", yaxt="n", xlab="", ylab="", bty="n")
text(x=0, y=1, lab="a", cex=2.5, font=2)

par(new=TRUE, mar=c(5,2,5,0))
plt = plot.phylo(tree, cex=1.5, direction="rightward")
x_axis = round(seq(0, max(tree$edge.length), length=10))
axis(side=1, line=1.5, at=max(x_axis)-x_axis, lab=x_axis)
mtext(text="Million years ago", side=1, line=4.5, at=median(x_axis))

### Expansion / Contraction: middle area text
par(mar=c(5,0,5,0))
plot(x=plt$x.lim, y=plt$y.lim, type="n", bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
adj_frac = 0.04
conex = read.table(fname_conex, header=TRUE)
conex = conex[order(conex$Species), ]
tips.order = unlist(lapply(conex$Species, FUN=function(x) {
                            which(x == species_order)}))
conex$order = tips.order
conex = conex[order(conex$order, decreasing=TRUE), ]
conex$Expansion = formatC(conex$Expansion, format="d", big.mark=",")
conex$Contraction = formatC(conex$Contraction, format="d", big.mark=",")
conex_lab = paste0(conex$Expansion, " : ", conex$Contraction)
# text(x=(plt$x.lim[2]-(adj_frac*plt$x.lim[2])), y=seq(plt$y.lim[1], plt$y.lim[2]), adj=0.5, lab=conex_lab, cex=1.2)
text(x=median(plt$x.lim), y=seq(plt$y.lim[1], plt$y.lim[2]), adj=0.5, lab=conex_lab, cex=1.5)
mtext(side=3, line=1, at=median(plt$x.lim), adj=0.449, text="Expansion : Contraction")

### Gene classifications: bar plot
gene_groups = read.csv(fname_gene_groups)
m = ncol(gene_groups)
gene_groups = gene_groups[order(gene_groups$Species), ]
rownames(gene_groups) = c(1:nrow(gene_groups))
gene_groups$order = tips.order
gene_groups = gene_groups[order(gene_groups$order, decreasing=TRUE), ]
X = t(as.matrix(gene_groups[, 3:m]))
rownames(X) = gsub("_", " ", colnames(gene_groups)[3:m])
colnames(X) = gene_groups$Species
colours = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c")
par(mar=c(3.5, 3.5, 3.5, 0.0))
barplot(X, col=colours, bord=NA, horiz=TRUE, yaxt="n", xaxt="n", xlim=c(0, signif(max(gene_groups$Total),0)))
x_axis = seq(0, signif(max(gene_groups$Total),0), length=5)
axis(side=1, at=x_axis, lab=formatC(x_axis, format="d", big.mark=","))
mtext(text="Gene counts", side=1, line=3, at=median(x_axis))
### Legend of the barplot
par(mar=c(0,0,0,0))
plot(0, 0, type="n",  xaxt="n", yaxt="n", bty="n")
legend("center", legend=rownames(X), fill=colours, bty="n")


### Venn diagram of shared gene families
venn_species_max_5 = gsub(" ", "_", read.delim(fname_venn_species_max_5, header=FALSE)$V1)
gene_counts = read.delim(fname_gene_counts, header=TRUE)
X = gene_counts[, 2:(ncol(gene_counts)-4)] > 0
Y = list()
colours = rainbow(length(species_order))
for (name in species_order){
    j = c(1:ncol(X))[colnames(X)==name]
    if (sum(species_order %in% name) == 0){
        next
    }
    name = gsub("_", " ", name)
    if (nchar(name) > 12) {
        vec_name = unlist(strsplit(name, " "))
        name = paste0(unlist(strsplit(vec_name[1], ""))[1], ". ", paste(tail(vec_name, -1), collapse=" "))
    }
    eval(parse(text=paste0("Y$`", name, "` = c(1:nrow(X))[X[, j]]")))
}

idx = c(1:length(species_order))[species_order %in% venn_species_max_5]
vp = venn.diagram(Y[idx], col=colours[idx], fill=colours[idx], alpha=0.3, filename=NULL)
par(mar=c(0,0,0,0))
plot(x=c(0,1), y=c(0,1), xaxt="n", yaxt="n", xlab="", ylab="", type="n", bty="n")
text(x=0.05, y=0.95, lab="b", cex=2.5, font=2)
pushViewport(plotViewport(layout.pos.col=1, layout.pos.row=2))
grid.draw(vp)
### Distribution of 4DTv (fraction of transverions among 4-fold degenerate codons - correlated with time from whole genome duplication using dual-copy paralogs and single-copy orthologs)
###@@@ Extract within genome 4DTv (pairwise paralogs)
FDTv_files = file.path(dir_name_4DTv, list.files(path=dir_name_4DTv, pattern=paste0(extension_name_4DTv, "$")))
FDTv_files = FDTv_files[grep(fname_4DTv_singlecopy, FDTv_files, invert=TRUE)]
id = c()
x = c()
y = c()
for (f in FDTv_files){
    # f = FDTv_files[1]
    df = read.delim(f, header=FALSE, na.string="NA")
    df[is.na(df)] = 0
    d = density(df$V4)
    id = c(id, gsub("_", " ", rep(gsub(".4DTv", "", basename(f)), length(d$x))))
    x = c(x, d$x)
    # y = c(y, (d$y - min(d$y))/diff(range(d$y)))
    y = c(y, d$y)
}
###@@@ Extract across genomes 4DTv (pairwise single-copy orthologs) and append to x and y
dat = read.delim(fname_4DTv_singlecopy, header=TRUE)
dat$SPECIES_1 = as.factor(dat$SPECIES_1)
dat$SPECIES_2 = as.factor(dat$SPECIES_2)
for (species1 in levels(dat$SPECIES_1)){
    # species1 = levels(dat$SPECIES_1)[1]
    for (species2 in levels(dat$SPECIES_2)){
        # species2 = levels(dat$SPECIES_2)[2]
        if (species1 == species2){
            next
        }
        df = droplevels(dat[(dat$SPECIES_1==species1) & (dat$SPECIES_2==species2), ])
        df[is.na(df)] = 0
        d = density(df$X4DTv)
        id = c(id, gsub("_", " ", rep(paste0(species1, " X ", species2), length(d$x))))
        x = c(x, d$x)
        # y = c(y, (d$y - min(d$y))/diff(range(d$y)))
        y = c(y, d$y)
    }
}
df = data.frame(id=as.factor(id), x=x, y=y)

species_labels = gsub("_", " ", species_order)
species_labels = expand.grid(x=species_labels, y=species_labels)
species_labels$x = as.character(species_labels$x)
species_labels$y = as.character(species_labels$y)
species_list = c()
for (i in 1:nrow(species_labels)){
    # i = 1
    row = species_labels[i, ]
    if (row[1] == row[2]) {
        species_list = c(species_list, row[1])
    } else {
        species_list = c(species_list, paste0(row[1], " X ", row[2]))
    }
}
df = droplevels(df[df$id %in% species_list, ])

### Keep only the comparisons you require
vec_comparisons = read.csv(fname_comparisons, header=FALSE)[,1]
vec_comparisons = gsub("_", " ", vec_comparisons)
df = droplevels(df[df$id %in% vec_comparisons, ])
if (nrow(df) == 0) {
    df = data.frame(id=as.factor(id), x=x, y=y)
    df = droplevels(df[df$id %in% species_list, ])
    df = droplevels(df[df$id %in% vec_comparisons, ])
}

par(mar=c(5, 5, 5, 2))
species_list = levels(df$id)
n = nlevels(df$id)
# colours = c(colours[c(1,2,4)], "#969696", "#66c2a5")
colours = rainbow(n)
plot(0, 0, xlim=range(df$x), ylim=range(df$y), xlab="4DTv", ylab="Density", type="n")
for (i in 1:n){
    # i = 1
    id = species_list[i]
    subdf = df[df$id==id, ]
    lines(x=subdf$x, y=subdf$y, col=colours[i], lty=i, lwd=2.5)
}
grid()
for (i in 1:length(species_list)){
    # i = 1
    entry = species_list[i]
    if (nchar(entry) > 12) {
        vec_species_names = unlist(strsplit(entry, " X "))
        for (j in 1:length(vec_species_names)) {
            # j = 1
            name = vec_species_names[j]
            vec_name = unlist(strsplit(name, " "))
            vec_species_names[j] = paste0(unlist(strsplit(vec_name[1], ""))[1], ". ", paste(tail(vec_name, -1), collapse=" "))
        }
        species_list[i] = paste(vec_species_names, collapse=" X ")
    }
}
legend("topright", legend=species_list, col=colours, lty=c(1:n), lwd=2.5, cex=1)

par(new=TRUE, mar=c(0,0,0,0))
plot(x=c(0,1), y=c(0,1), xaxt="n", yaxt="n", xlab="", ylab="", type="n", bty="n")
text(x=0.05, y=0.95, lab="c", cex=2.5, font=2)

dev.off()

### Clean-up
system("rm VennDiagram.*.log")
