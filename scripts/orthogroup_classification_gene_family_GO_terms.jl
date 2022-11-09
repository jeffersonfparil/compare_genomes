using CSV, DataFrames, ProgressMeter

fname_orthogroup_family_hits =  ARGS[1]
fname_family_GO =               ARGS[2]
fname_paralogs =                ARGS[3]
fname_orthogroup_gene_counts =  ARGS[4]
fname_unassigned_genes =        ARGS[5]
fname_output =                  ARGS[6]

# fname_orthogroup_family_hits = "ORTHOGROUPS/orthogroups.pthr"
# fname_family_GO = "PantherHMM_17.0/Panther17.0_HMM_familyIDs.txt"
# fname_paralogs = "all_orthogroups.tmp"
# fname_orthogroup_gene_counts = "ORTHOGROUPS/OrthoFinder/Results_May05/Orthogroups/Orthogroups.GeneCount.tsv"
# fname_unassigned_genes = "ORTHOGROUPS/OrthoFinder/Results_May05/Orthogroups/Orthogroups_UnassignedGenes.tsv"
# fname_output = "ORTHOGROUPS/orthogroups_gene_counts_families_go.out"

# Load all orthogroup IDs
file = open(fname_paralogs, "r")
all_orthogroups = readlines(file)
close(file)
rm(fname_paralogs) # clean-up

# Load orthogroup hits
file = open(fname_orthogroup_family_hits, "r")
seekend(file); n = position(file); seekstart(file)
orthogroup = []
gene_name = []
family_ID = []
evalue = []
pb = Progress(n)
while !eof(file)
    line = split(readline(file), " "[1])
    seqName = split(line[1], ":"[1])
    speciesAndGeneName = split(seqName[2], "|"[1])
    geneName = join(split(speciesAndGeneName[2], "-"[1])[2:end], "-"[1]) ### remove GeMoMa prefix (hyphen-delimited)
    if geneName == ""
        geneName = speciesAndGeneName[2]
    end
    geneName = replace(geneName, Regex("_R0\$")=>"") ### remove GeMoMa suffix (Regex("_R0\$"))
    push!(orthogroup, seqName[1])
    push!(gene_name, geneName)
    push!(family_ID, replace(line[2], ".orig.30.pir"=>""))
    push!(evalue, parse(Float64, line[3]))
    update!(pb, position(file))
end
close(file)

# Load PantherHMM family descriptions
file = open(fname_family_GO, "r")
seekend(file); n = position(file); seekstart(file)
PTHR_family_ID = []
PTHR_family_name = []
PTHR_GO_term = []
pb = Progress(n)
while !eof(file)
    line = split(readline(file), "\t"[1])
    push!(PTHR_family_ID, line[1])
    push!(PTHR_family_name, line[2])
    push!(PTHR_GO_term, line[3])
    update!(pb, position(file))
end
close(file)

# Load gene counts orthogroup per species
df_counts = CSV.read(open(fname_orthogroup_gene_counts), DataFrames.DataFrame)

# Load list of unassigned genes, i.e. orthogroups with a single gene specific to each species
df_unassigned = CSV.read(open(fname_unassigned_genes), DataFrames.DataFrame)

# For each gene, set the family ID as the one with the lowest E-value
function extract_best_matching_family_per_orthogroup(gene_name, orthogroup, family_ID, evalue)::Tuple{Vector{String}, Vector{String}}
    idx = sortperm(gene_name)
    orthogroup = orthogroup[idx]
    gene_name = gene_name[idx]
    family_ID = family_ID[idx]
    evalue = evalue[idx]
    all_gene = sort(unique(gene_name))
    all_family_ID = []
    all_orthogroup = []
    i = 1
    @showprogress for g in all_gene
        # g = all_gene[1]
        t = g == gene_name[i]
        while t == false
            i += 1
            t = g == gene_name[i]
        end
        f = []
        o = []
        e = []
        while t
            push!(f, family_ID[i])
            push!(o, orthogroup[i])
            push!(e, evalue[i])
            i += 1
            t = try
                    g == gene_name[i]
                catch
                    false
                end
        end
        f = f[e .== minimum(e)][1]
        o = o[e .== minimum(e)][1]
        push!(all_family_ID, f)
        push!(all_orthogroup, o)
    end
    return(all_family_ID, all_orthogroup)
end
_, _ = extract_best_matching_family_per_orthogroup(gene_name[1:5], orthogroup[1:5], family_ID[1:5], evalue[1:5])
all_family_ID, all_orthogroup = extract_best_matching_family_per_orthogroup(gene_name, orthogroup, family_ID, evalue)

# Identify the PantherHMM gene family names and GO terms
function identiy_gene_family_names_and_GO_terms(all_family_ID, PTHR_family_name, PTHR_GO_term)::Tuple{Vector{String}, Vector{String}}
    all_family = []
    all_GO = []
    @showprogress for ID in all_family_ID
        # ID = all_family_ID[1]
        idx = ID .== PTHR_family_ID
        if sum(idx) == 1
            push!(all_family, PTHR_family_name[idx][1])
            push!(all_GO, PTHR_GO_term[idx][1])
        else
            # for unclassified orthogroups
            push!(all_family, "UNKNOWN")
            push!(all_GO, "")
        end
    end
    return(all_family, all_GO)
end
_, _ = identiy_gene_family_names_and_GO_terms(all_family_ID[1:5], PTHR_family_name, PTHR_GO_term)
all_family, all_GO = identiy_gene_family_names_and_GO_terms(all_family_ID, PTHR_family_name, PTHR_GO_term)

# Summarise gene families per orthogroup
function summarise_gene_families_per_orthogroup(all_family_ID, all_orthogroup, all_family, all_GO)::DataFrames.DataFrame
    idx = sortperm(all_orthogroup)
    all_orthogroup = all_orthogroup[idx]
    all_family_ID = all_family_ID[idx]
    all_family = all_family[idx]
    all_GO = all_GO[idx]
    final_orthogroup = unique(all_orthogroup)
    sort!(final_orthogroup)
    final_family_ID = []
    final_family = []
    final_GO = []
    i = 1
    @showprogress for o in final_orthogroup
        t = o == all_orthogroup[i]
        while t == false
            i += 1
            t = o == all_orthogroup[i]
        end
        fid = []
        fam = []
        fgo = []
        while t
            push!(fid, all_family_ID[i])
            push!(fam, all_family[i])
            push!(fgo, all_GO[i])
            i += 1
            t = try
                    o == all_othogroup[i]
                catch
                    false
                end
        end
        push!(final_family_ID, join(unique(fid), ";"[1]))
        push!(final_family, join(unique(fam), ";"[1]))
        push!(final_GO, join(unique(fgo), ";"[1]))
    end

    df_ID = DataFrames.DataFrame(Orthogroup=final_orthogroup,
                            Family_ID=final_family_ID,
                            Family=final_family,
                            GO=final_GO)
    return(df_ID)
end
_ = summarise_gene_families_per_orthogroup(all_family_ID[1:5], all_orthogroup[1:5], all_family[1:5], all_GO[1:5])
df_ID = summarise_gene_families_per_orthogroup(all_family_ID, all_orthogroup, all_family, all_GO)

# Generate gene counts per species for the set of unassigned genes
df_append_unassigned = Int.(.!ismissing.(df_unassigned[:, 2:end]))
df_append_unassigned.Total = repeat([1], nrow(df_unassigned))
df_append_unassigned.Orthogroup = df_unassigned.Orthogroup

# Append unassigned orthogroups into the orthogroup gene counts dataframe
df_counts = vcat(df_counts, df_append_unassigned)

# Merge the orthogroup ID and orthogroup gene counts and save into a file
df = outerjoin(df_counts, df_ID, on=:Orthogroup)
CSV.write(open(fname_output, "w"), df, delim="\t"[1])
