fname = ARGS[1]
fname_out = ARGS[2]
using CSV, DataFrames, ProgressMeter
file = open(fname, "r")
df = CSV.read(file, DataFrames.DataFrame, header=true)
close(file)
counts = df[:, 2:(end-4)]
species_names = names(counts)

# Unassigned gene
function count_unassigned_genes(counts, species_names)
    println("Count unassigned genes")
    idx = df.Total .== 1
    unassigned = []
    for species in species_names
        # species = species_names[1]
        push!(unassigned, sum(counts[idx, species_names .== species][:,1]))
    end
    return(unassigned)
end
@time unassigned = count_unassigned_genes(counts, species_names)

# Unique paralogs
function count_unique_paralogs(df, counts, species_names)
    println("Count genes belonging to unique paralogs for each species")
    unique_paralogs = []
    for species in species_names
        # species = species_names[1]
        idx = species_names .== species
        push!(unique_paralogs, sum((counts[:, idx] .== df.Total)[:,1]))
    end
    return(unique_paralogs)
end
@time unique_paralogs = count_unique_paralogs(df, counts, species_names)

# Single-copy gene orthologs
function count_single_copy_gene_orthologs(count, species_names)
    println("Count single-copy gene orthologs")
    idx = repeat([true], nrow(counts))
    for j in 1:ncol(counts)
        idx = idx .& (counts[:, j] .== 1)
    end
    single_copy_orthologs = repeat([sum(idx)], length(species_names))
    return(single_copy_orthologs)
end
@time single_copy_orthologs = count_single_copy_gene_orthologs(count, species_names)

# Multiple-copy orthologs
function count_multicopy_orthologs(counts, species_names)
    println("Count multi-copy orthologs")
    total_per_species = []
    for species in species_names
        push!(total_per_species, sum(counts[:, species.==species_names][:,1]))
    end
    multiple_orthologs = total_per_species - (unassigned + unique_paralogs + single_copy_orthologs)
    return(total_per_species, multiple_orthologs)
end
@time total_per_species, multiple_orthologs = count_multicopy_orthologs(counts, species_names)

# Merge gene count classifications per species and save
out = DataFrames.DataFrame(Species=species_names,
                           Total=total_per_species,
                           Multiple_Orthologs=multiple_orthologs,
                           Single_Copy_Orthologs=single_copy_orthologs,
                           Unique_Paralogs=unique_paralogs,
                           Unassigned_genes=unassigned
                          )
file = open(fname_out, "w")
CSV.write(file, out)
close(file)
