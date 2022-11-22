### Calculate 4DTv: ratio of transversions in 4-fold degenerate sites

### 4-fold degenerate codons:
### (1) Ala - GCN
### (2) Arg - CGN (etc)
### (3) Gly - GGN
### (4) Leu - CTN (etc)
### (5) Pro - CCN
### (6) Ser - TCN (etc)
### (7) Thr - ACN
### (8) Val - GTN

## Transversions:
## A <-> C
## A <-> T
## G <-> C
## G <-> T

fname = ARGS[1]
fname_output =  try
                    ARGS[2]
                catch
                    string(fname, ".4DTv")
                end

struct Seq
    name::String
    seq::String
end

struct Cod
    name::String
    cod::Vector{String}
end

struct CodPair
    name::String
    cod1::Vector{String}
    cod2::Vector{String}
end

function extract_sequences(fname)::Vector{Seq}
    file = open(fname, "r")
    name = ""
    seq = ""
    vec_seq = []
    while !eof(file)
        line = readline(file)
        if (line[1] == '>') & (name == "")
            name = replace(line, ">"=>"")
        elseif (line[1] == '>') & (seq != "")
            push!(vec_seq, Seq(name, seq))
            name = replace(line, ">"=>"")
            seq = ""
        elseif (line[1] != '>')
            seq = string(seq, line)
        end
    end
    push!(vec_seq, Seq(name, seq)) ### Final sequence
    close(file)
    return(vec_seq)
end

function sequence_to_codons(seq::Seq)::Cod
    if length(seq.seq) % 3 != 0
        println("Input sequence is not a multiple of 3!")
        println("We cannot divide the sequence into codons!")
        println("Please check the input file.")
        exit()
    end
    cod = []
    for i in 1:3:length(seq.seq)
        push!(cod, seq.seq[i:(i+2)])
    end
   return(Cod(seq.name, cod))
end

function pair_up_sequences_and_convert_into_codons(vec_seq::Vector{Seq})::Vector{CodPair}
    vec_cod_pairs = []
    for i in 1:(length(vec_seq)-1)
        cod1 = sequence_to_codons(vec_seq[i])
        for j in (i+1):length(vec_seq)
            cod2 = sequence_to_codons(vec_seq[j])
            if length(cod1.cod) != length(cod2.cod)
                println("Sequences are not the same length!")
                println("Please make sure the sequences have been aligned.")
                exit()
            end
            name = string(cod1.name, ":", cod2.name)
            push!(vec_cod_pairs, CodPair(name, cod1.cod, cod2.cod))
        end
    end
    return(vec_cod_pairs)
end

function calculate_4DTv(vec_cod_pairs::Vector{CodPair})::Tuple{Vector{String}, Vector{Int64}, Vector{Int64}, Vector{Float64}}
    names = []
    n_4D = []
    n_4DTv = []
    for i in 1:length(vec_cod_pairs)
        # i = 1
        cod_pair =vec_cod_pairs[i]
        push!(names, cod_pair.name)
        push!(n_4D, 0)
        push!(n_4DTv, 0)
        for j in 1:length(cod_pair.cod1)
            # j = 1
            cod1 = cod_pair.cod1[j]
            cod2 = cod_pair.cod2[j]
            # Skip alignment gaps
            if (match(Regex("-"), cod1)!=nothing) | (match(Regex("-"), cod2)!=nothing)
                continue
            end
            # Identify 4-fold degenerate codons
            bool1 = ((match(Regex("GC."), cod1)!=nothing) & (match(Regex("GC."), cod2)!=nothing)) |
                    ((match(Regex("CG."), cod1)!=nothing) & (match(Regex("CG."), cod2)!=nothing)) |
                    ((match(Regex("GG."), cod1)!=nothing) & (match(Regex("GG."), cod2)!=nothing)) |
                    ((match(Regex("CT."), cod1)!=nothing) & (match(Regex("CT."), cod2)!=nothing)) |
                    ((match(Regex("CC."), cod1)!=nothing) & (match(Regex("CC."), cod2)!=nothing)) |
                    ((match(Regex("TC."), cod1)!=nothing) & (match(Regex("TC."), cod2)!=nothing)) |
                    ((match(Regex("AC."), cod1)!=nothing) & (match(Regex("AC."), cod2)!=nothing)) |
                    ((match(Regex("GT."), cod1)!=nothing) & (match(Regex("GT."), cod2)!=nothing))
            # Identify transversion
            if bool1
                n_4D[i] += 1
                bool2 = ((cod1[3] == 'A') & (cod2[3] == 'T')) |
                        ((cod1[3] == 'A') & (cod2[3] == 'C')) |
                        ((cod1[3] == 'G') & (cod2[3] == 'T')) |
                        ((cod1[3] == 'G') & (cod2[3] == 'C')) |
                        ((cod1[3] == 'T') & (cod2[3] == 'A')) |
                        ((cod1[3] == 'T') & (cod2[3] == 'G')) |
                        ((cod1[3] == 'C') & (cod2[3] == 'A')) |
                        ((cod1[3] == 'C') & (cod2[3] == 'G'))
                if bool2
                    n_4DTv[i] += 1
                end
            end
        end
    end
    return(names, n_4D, n_4DTv, n_4DTv ./ n_4D)
end

### Extract sequences
vec_seq = extract_sequences(fname)
### Transform into codon pairs
vec_cod_pairs = pair_up_sequences_and_convert_into_codons(vec_seq)
### Calculate 4DTv
names, n_4D, n_4DTv, out_4DTv = calculate_4DTv(vec_cod_pairs)
### Save 4DTv (Tab-delimited: Pair name (delimited by ":"), number of 4-fold degenerate codons, number of 4-fold degenerate codons with transversion, ratio of the previous 2 columns, i.e. 4DTv)
file = open(fname_output, "a")
for i in 1:length(out_4DTv)
    line = string(join([names[i],
                        n_4D[i],
                        n_4DTv[i],
                        out_4DTv[i]], "\t")
                 , "\n")
    write(file, line)
end
close(file)
