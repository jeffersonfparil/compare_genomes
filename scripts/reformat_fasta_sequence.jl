using ProgressMeter

### BE SURE TO REPLACE "|" in the sequence_name_query input

fasta_input = ARGS[1]
c = parse(Int, ARGS[2])
fasta_output = try 
    ARGS[3]
catch
    string(join(split(fasta_input, ".")[1:(end-1)], "."), "-REFORMATTED.fasta")
end

function write_per_line(file_output, line, c)
    file_output = open(fasta_output, "a")
    x = collect(line)
    for i in 1:c:length(x)
        try 
            write(file_output, string(join(x[i:(i+c-1)], ""), "\n"))
        catch
            write(file_output, string(join(x[i:end], ""), "\n"))
        end
    end
    close(file_output)
    return(0)
end

file_input = open(fasta_input, "r")
seekend(file_input); n = position(file_input)
seekstart(file_input)
pb = Progress(n)
while !eof(file_input)
    line = readline(file_input)
    if line[1] == '>'
        file_output = open(fasta_output, "w")
        write(file_output, string(line, '\n'))
        close(file_output)
        line = readline(file_input)
        SEQ=line
        bool_test = line[1] != '>'
        while bool_test
            line = readline(file_input)
            SEQ = string(SEQ, line)
            bool_test = try
                line[1] != '>'
            catch
                false
            end
            update!(pb, position(file_input))
        end
        write_per_line(fasta_output, SEQ, c)
    end
    update!(pb, position(file_input))
end
close(file_input)
