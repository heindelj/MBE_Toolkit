function read_xyz(ifile::String)
    header = Array{String, 1}()
    atom_labels = Array{Array{String, 1}, 1}()
    geoms = Array{Array{Float64, 2}, 1}()
    open(ifile) do io
        for line in eachline(io)
            if isa(tryparse(Int, line), Int)
                # allocate the geometry for this frame
                N = parse(Int, line)
                labels = String[]
                geom = zeros((3, N))
                # store the header for this frame
                head = string(line, readline(io))
                push!(header, head)
                # loop through the geometry storing the vectors and atom labels as you go
                for i = 1:N
                    coords = split(readline(io))
                    push!(labels, coords[1])
                    geom[:,i] = parse.(Float64, coords[2:end])
                end
                push!(atom_labels, labels)
                push!(geoms, geom)
            end
        end
    end
return header, atom_labels, geoms
end