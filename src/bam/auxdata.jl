# BAM Auxiliary Data
# ==================
# Note: Because all this code is inherently type unstable, using multiple
# dispatch will result in less efficient code than spamming if/else statements.
# Therefore the code here is distinctly un-Julian.

struct AuxDataIterator
    start::Int
    stop::Int
    data::Vector{UInt8}
end

function AuxDataIterator(record::BAMRecord)
    return AuxDataIterator(auxdata_position(record), data_size(record), record.data)
end

function Base.iterate(aux::AuxDataIterator, pos=aux.start::Int)
    if pos > aux.stop
        return nothing
    end
    tag1, tag2, type1, type2, value, pos = get_next_auxpair(aux.data, pos)
    return Pair{String,Any}(String([tag1, tag2]), value), pos
end

function auxdata(record::BAMRecord)::Dict{String,Any}
    return Dict{String,Any}(AuxDataIterator(record))
end

function Base.getindex(record::BAMRecord, tag::AbstractString)
    pos = findauxtag(record, tag)
    if pos === nothing
        throw(KeyError(tag))
    end
    return load_auxpair(record.data, pos)[4]
end

function Base.get(record::BAMRecord, tag::AbstractString, default::Any)
    pos = findauxtag(record, tag)
    if pos === nothing
        return default
    end
    return load_auxpair(record.data, pos)[4]
end

function Base.haskey(record::BAMRecord, tag::AbstractString)::Bool
    return findauxtag(record, tag) !== nothing
end

function Base.pop!(record::BAMRecord, tag::AbstractString)
    pos = findauxtag(record, tag)
    if pos === nothing
        throw(KeyError(tag))
    end
    return pop_pos!(record, pos)
end

function Base.pop!(record::BAMRecord, tag::AbstractString, default)
    pos = findauxtag(record, tag)
    return pos === nothing ? default : pop_pos!(record, pos)
end

function Base.delete!(record::BAMRecord, tag::AbstractString)
    pos = findauxtag(record, tag)
    if pos === nothing
        throw(KeyError(tag))
    end
    pop_pos!(record, pos)
    return record
end

const TYPESTRINGS = Dict{DataType,String}(UInt8=>"C", Vector{UInt8}=>"BC",
                       Int8=>"c", Vector{Int8}=>"Bc",
                       UInt16=>"S", Vector{UInt16}=>"BS",
                       Int16=>"s", Vector{Int16}=>"Bs",
                       UInt32=>"I", Vector{UInt32}=>"BI",
                       Int32=>"i", Vector{Int32}=>"Bi",
                       Float32=>"f", Vector{Float32}=>"Bf",
                       String=>"Z", Char=>"A")

function Base.setindex!(record::BAMRecord, value, key::AbstractString)
    data = record.data
    tag = get(TYPESTRINGS, typeof(value), nothing)
    if tag === nothing
        throw(ArgumentError("Cannot encode data of type $(typeof(value)) in BAM aux fields."))
    end
    # If tag already in data, we can overwrite the data if it's a bitstype,
    # or else we must delete it, then re-add it.
    pos = findauxtag(record, key)
    if pos !== nothing
        if isbits(value) && load_aux_bitstype(data, pos + 2) === typeof(value)
            unsafe_store!(Ptr{typeof(value)}(pointer(data, pos+3)), value)
            return value
        else
            pop_pos!(record, pos)
        end
    end
    # Update size of data array to fit
    index = data_size(record)
    # Null byte for string, length int32 for vector, char is only 1 byte, not 4
    padding = value isa String ? 1 : value isa Vector ? 4 : value isa Char ? -3 : 0
    data_increase =  2 + sizeof(tag) + sizeof(value) + padding
    resize!(data, index + data_increase)
    # Add tag
    data[index + 1] = codeunit(key, 1)
    data[index + 2] = codeunit(key, 2)
    index += 3
    # Add type characters
    for char in tag
        data[index] = UInt8(char)
        index += 1
    end
    write_value(data, index, value)
    record.block_size += data_increase
    return value
end

function write_value(data::Vector{UInt8}, index::Int, value::String)
    unsafe_copyto!(pointer(data, index), pointer(value), sizeof(value))
    data[end] = 0x00 # null terminate
end

function write_value(data::Vector{UInt8}, index::Int, value::Vector)
    # Store first Int32 with length, then vector itself
    unsafe_store!(Ptr{Int32}(pointer(data, index)), length(value))
    unsafe_copyto!(Ptr{eltype(value)}(pointer(data, index+4)), pointer(value), length(value))
end

function write_value(data::Vector{UInt8}, index::Int, value::Char)
    @inbounds data[index] = UInt8(value)
end

function write_value(data::Vector{UInt8}, index::Int, value)
    unsafe_store!(Ptr{typeof(value)}(pointer(data, index)), value)
end

# Internals
# ---------
function get_next_auxpair(data::Vector{UInt8}, pos::Int)::Tuple{UInt8, UInt8, UInt8, UInt8, Any, Int}
    @inbounds tag1 = data[pos]
    @inbounds tag2 = data[pos+1]
    pos, type1, type2, value = load_auxpair(data, pos)
    return tag1, tag2, type1, type2, value, pos
end

function unsafe_write_to_buffer(buffer::IOBuffer, data::Vector{UInt8}, pos::Int)
    tag1, tag2, type1, type2, value, pos = get_next_auxpair(data, pos)
    write(buffer, tag1, tag2, UInt8(':'))
    if type1 == UInt8('C')
        write(buffer, UInt8('i'), UInt8(':'))
        write_int(buffer, value::UInt8)
    elseif type1 == UInt8('c')
        write(buffer, UInt8('i'), UInt8(':'))
        write_int(buffer, value::Int8)
    elseif type1 == UInt8('A')
        write(buffer, UInt8('i'), UInt8(':'))
        print(buffer, value::Char)
    elseif type1 == UInt8('S')
        write(buffer, UInt8('i'), UInt8(':'))
        write_int(buffer, value::UInt16)
    elseif type1 == UInt8('s')
        write(buffer, UInt8('i'), UInt8(':'))
        write_int(buffer, value::Int16)
    elseif type1 == UInt8('I')
        write(buffer, UInt8('i'), UInt8(':'))
        write_int(buffer, value::UInt32)
    elseif type1 == UInt8('i')
        write(buffer, UInt8('i'), UInt8(':'))
        write_int(buffer, value::Int32)
    elseif type1 == UInt8('f')
        write(buffer, type1, UInt8(':'))
        print(buffer, value::Float32)
    elseif type1 == UInt8('Z')
        write(buffer, type1, UInt8(':'))
        print(buffer, value::String)
    elseif type1 == UInt8('B')
        write(buffer, type1, UInt8(':'), type2)
        elt = eltype(value::Vector)
        @inbounds for i in 1:length(value::Vector)
            write(buffer, UInt8(','))
            print(buffer, value[i]::elt)
        end
    end
    return pos
end

function findauxtag(record::BAMRecord, tag::AbstractString)
    if sizeof(tag) != 2
        throw(ArgumentError("tag length must be 2"))
    end
    start = auxdata_position(record)
    stop = data_size(record)
    return findauxtag(record.data, start, stop, codeunit(tag, 1), codeunit(tag, 2))
end

function findauxtag(data::Vector{UInt8}, start::Int, stop::Int, t1::UInt8, t2::UInt8)
    pos = start
    while pos ≤ stop && !(data[pos] == t1 && data[pos+1] == t2)
        pos = next_tag_position(data, pos)
    end
    return pos > stop ? nothing : pos
end

function load_aux_bitstype(data::Vector{UInt8}, p::Int)
    b = data[p]
    return b == UInt8('A') ? Char  :
    b == UInt8('c') ? Int8  :
    b == UInt8('C') ? UInt8 :
    b == UInt8('s') ? Int16 :
    b == UInt8('S') ? UInt16 :
    b == UInt8('i') ? Int32 :
    b == UInt8('I') ? UInt32 :
    b == UInt8('f') ? Float32 :
    throw(ValueError("invalid type tag: '$(Char(b))'"))
end

function load_auxpair(data::Vector{UInt8}, p)::Tuple{Int, UInt8, UInt8, Any}
    p += 2
    @inbounds type1 = data[p]
    @inbounds type2 = data[p+1]
    p += 1
    if type1 == UInt8('A')
        return p + 1, Char(type2)
    elseif type1 == UInt8('C')
        return p + sizeof(UInt8), type1, type2, type2
    elseif type1 == UInt8('c')
        return p + sizeof(Int8), type1, type2, reinterpret(Int8, type2)
    elseif type1 == UInt8('S')
        return p + sizeof(UInt16), type1, type2, unsafe_load(Ptr{UInt16}(pointer(data, p)))
    elseif type1 == UInt8('s')
        return p + sizeof(Int16), type1, type2, unsafe_load(Ptr{Int16}(pointer(data, p)))
    elseif type1 == UInt8('I')
        return p + sizeof(UInt32), type1, type2, unsafe_load(Ptr{UInt32}(pointer(data, p)))
    elseif type1 == UInt8('i')
        return p + sizeof(Int32), type1, type2, unsafe_load(Ptr{Int32}(pointer(data, p)))
    elseif type1 == UInt8('f')
        return p + sizeof(Float32), type1, type2, unsafe_load(Ptr{Float32}(pointer(data, p)))
    elseif type1 == UInt8('Z')
        str = unsafe_string(pointer(data, p))
        return p + sizeof(str) + 1, type1, type2, str
    elseif type1 == UInt8('H')
        str = unsafe_string(pointer(data, p))
        return p + sizeof(str) + 1, type1, type2, hex2bytes(str)
    elseif type1 == UInt8('B')
        eltype = load_aux_bitstype(data, p)
        p += 1
        n = unsafe_load(Ptr{Int32}(pointer(data, p)))
        vector = Vector{eltype}(undef, n)
        unsafe_copyto!(pointer(vector), Ptr{eltype}(pointer(data, p)), n)
        return p + n * sizeof(eltype) + 4, type1, type2, vector
    else
        throw(ArgumentError("invalid type tag: '$(Char(t))'"))
    end
end

const TYPESIZES = let TYPESIZES = fill(Int8(0), 256)
    TYPESIZES[Int('A')] = 1
    TYPESIZES[Int('c')] = 1
    TYPESIZES[Int('C')] = 1
    TYPESIZES[Int('s')] = 2
    TYPESIZES[Int('S')] = 2
    TYPESIZES[Int('i')] = 4
    TYPESIZES[Int('I')] = 4
    TYPESIZES[Int('f')] = 4
    TYPESIZES[Int('B')] = -1
    TYPESIZES[Int('H')] = -2
    TYPESIZES[Int('Z')] = -2
    Tuple(TYPESIZES)
end

# Find next starting position of tag after the tag at p.
# (data[p], data[p+1]) is supposed to be a current tag.
function next_tag_position(data::Vector{UInt8}, p::Int)::Int
    @inbounds typ = data[p+2]
    @inbounds size = TYPESIZES[typ]
    p += 3
    if size == 0
        error("invalid type tag: '$(Char(typ))'")
    elseif size > 0
        p += size
    elseif size == -2 # NULL-terminalted string or null terminated hex
        @inbounds while data[p] != 0x00
            p += 1
        end
        p += 1
    else # Array
        @inbounds eltyp = data[p]
        @inbounds elsize = TYPESIZES[eltyp]
        if elsize < 1
            error("invalid type tag: '$(Char(eltyp))'")
        end
        p += 1
        n = unsafe_load(Ptr{Int32}(pointer(data, p)))
        p += 4 + elsize * n
    end
    return p
end

# Remove the tag from position. Unsafe (must be present.)
function pop_pos!(record::BAMRecord, pos::Int)
    nextpos, type1, type2, val = load_auxpair(record.data, pos)
    datasize = data_size(record)
    # If the tag was not the last, we shift the data leftward
    if nextpos <= datasize
        endsize = length(record.data) - nextpos
        record.data[pos:pos+endsize] = record.data[nextpos:nextpos+endsize]
    end
    # Update the block size to exclude the removed data
    record.block_size -= (nextpos - pos)
    return val
end
