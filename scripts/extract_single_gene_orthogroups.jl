### Extract rows 
fname = ARGS[1]
nspecies = parse(Int, ARGS[2])
# fname = "PROTEOMES/orthogroups_gene_counts_families_go.out"; nspecies = 2

f = open(fname, "r")
seekend(f); n = position(f); seekstart(f)
out = open("single_gene_list.grep", "a")
header = readline(f)
while !eof(f)
    line = split(readline(f), '\t')
    counts = line[2:(nspecies+1)]
    if sum(counts .== "1") == length(counts)
        write(out, string(line[1], '\n'))
    end
end
close(f)
close(out)
