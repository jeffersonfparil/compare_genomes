# ARGS = ["SOD-OG0016003.aln.pw", "30", "6"]
fname_input = ARGS[1]
window_size = parse(Int, ARGS[2])
slide_size = parse(Int, ARGS[3])
fname_output = try
        ARGS[4]
    catch
    string(fname_input, ".window_slide.", window_size, "_", slide_size, ".pw")
    end

### Exit if window_size or slide_size are not a multiple of 3
if (window_size % 3 > 0) | (slide_size % 3 > 0)
    println("Window size and slide size should be a multiple of 3!")
    exit()
end

### Sequence alignment pair struct
struct Pair
    name::String
    seq1::String
    seq2::String
end

### Sequence pair splitting function
function split(pair::Pair, window_size::Int64, slide_size::Int64)::Vector{Pair}
    # pair = x[1]
    out = []
    n = length(pair.seq1)
    for i in 1:slide_size:n
        fin = i + (window_size-1)
        fin < n ? j = fin : j = n
        new_name=string(pair.name, "(", i, "-", j, ")")
        new_seq1 = pair.seq1[i:j]
        new_seq2 = pair.seq2[i:j]
        if length(new_seq1) < window_size
            break
        else
            push!(out, Pair(new_name, new_seq1, new_seq2))
        end
    end
    return(out)
end

### Save vector of sequence pairs into a file
function save(y::Vector{Pair}, fname_output::String)
    file = open(fname_output, "a")
    for p in y
        write(file, string(p.name, "\n"))
        write(file, string(p.seq1, "\n"))
        write(file, string(p.seq2, "\n"))
        write(file, string("\n"))
    end
    close(file)
end

### Extract pairs of sequence alignments
file = open(fname_input, "r")
x = []
while !eof(file)
    name = readline(file)
    seq1 = readline(file)
    seq2 = readline(file)
    _    = readline(file)
    push!(x, Pair(name, seq1, seq2))
end
close(file)
### Split and save windows of sequences
for p in x
    # p = x[1]
    y = split(p, window_size, slide_size)
    save(y, fname_output)
end
