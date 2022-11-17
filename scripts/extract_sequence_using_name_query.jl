using ProgressMeter

### BE SURE TO REPLACE "|" in the sequence_name_query input

fasta_input = ARGS[1]
sequence_name_query_vec = ARGS[2:(length(ARGS)-3)]
fasta_output = try 
                    ARGS[length(ARGS)-2]
                catch
                    ""
                end
new_sequence_name = try 
                    ARGS[length(ARGS)-1]
                catch
                    ""
                end
add_gene_coordinates = try 
                    parse(Bool, ARGS[length(ARGS)-0])
                catch
                    false
                end

if sequence_name_query_vec isa Vector
    if length(sequence_name_query_vec) > 1
        sequence_name_query = join(sequence_name_query_vec, " ")
    else
        sequence_name_query = sequence_name_query_vec[1]
    end
end

if fasta_output == ""
    fasta_output = string(join(split(fasta_input, ".")[1:(end-1)], "."), "-", sequence_name_query, ".out")
end

### Add escape characters in front of "."
if match(Regex("\\."), sequence_name_query) != nothing
    sequence_name_query = replace(sequence_name_query, "."=> "\\.")
end
### Remove return character "\r"
if match(Regex("\\r"), sequence_name_query) != nothing
    sequence_name_query = replace(sequence_name_query, "\r"=> "")
end

file_input = open(fasta_input, "r")
seekend(file_input); n = position(file_input)
seekstart(file_input)
pb = Progress(n)
while !eof(file_input)
    line = readline(file_input)
    if line[1] == '>'
        while match(Regex(sequence_name_query), replace(line, "|"=>":")) != nothing
            file_output = open(fasta_output, "a")
            vec_line = split(line, " ")
            if (new_sequence_name != "")
                line = string(">", new_sequence_name)
            end
            if add_gene_coordinates
                coordinates = try
                        vec_line[(match.(Regex("interval="), vec_line) .!= nothing) .| (match.(Regex("location="), vec_line) .!= nothing)][1]
                    catch
                        ""
                    end
                line = string(line, "(", coordinates, ")")
            end
            write(file_output, string(line, '\n'))
            line = readline(file_input)
            bool_test = line[1] != '>'
            while bool_test
                write(file_output, line)
                line = readline(file_input)
                bool_test = try
                    line[1] != '>'
                catch
                    false
                end
                update!(pb, position(file_input))
            end
            write(file_output, '\n')
            close(file_output)
        end
    end
    update!(pb, position(file_input))
end
close(file_input)
