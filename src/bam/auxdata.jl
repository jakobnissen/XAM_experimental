# BAM Auxiliary Data
# ==================

struct AuxDataIterator
    start::Int
    stop::Int
    data::Vector{UInt8}
end

function AuxDataIterator(record::BAMRecord)
    return AuxDataIterator(auxdata_position(record), data_size(record), record.data)
end

function Base.iterate(aux::AuxDataIterator, pos=aux.start::Int)
    data = aux.data
    if pos > aux.stop
        return nothing
    end
    @label doit
    t1 = data[pos]
    t2 = data[pos+1]
    pos, typ = loadauxtype(data, pos + 2)
    pos, value = loadauxvalue(data, pos, typ)
    if t1 == t2 == 0xff
        @goto doit
    end
    return Pair{String,Any}(String([t1, t2]), value), pos
end

function auxdata(record::BAMRecord)::Dict{String,Any}
    return Dict{String,Any}(AuxDataIterator(record))
end

function Base.getindex(record::BAMRecord, tag::AbstractString)
    pos = findauxtag(record, tag)
    if pos === nothing
        throw(KeyError(tag))
    end
    return loadauxvalue(record.data, pos)
end

function Base.get(record::BAMRecord, tag::AbstractString, default::Any)
    pos = findauxtag(record, tag)
    if pos === nothing
        return default
    end
    return loadauxvalue(record.data, pos)
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
                       String=>"Z")

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
        if isbits(value) && loadauxtype(data, pos + 2)[2] === typeof(value)
            unsafe_store!(Ptr{typeof(value)}(pointer(data, pos+3)), value)
            return value
        else
            pop_pos!(record, pos)
        end
    end
    # Update size of data array to fit
    index = data_size(record)
    padding = value isa String ? 1 : value isa Vector ? 4 : 0
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
    # Write value itself
    if value isa String
        unsafe_copyto!(pointer(data, index), pointer(value), sizeof(value))
        data[end] = 0x00 # null terminate
    elseif value isa Vector
        # Store first Int32 with length, then vector itself
        unsafe_store!(Ptr{Int32}(pointer(data, index)), length(value))
        unsafe_copyto!(Ptr{eltype(value)}(pointer(data, index+4)), pointer(value), length(value))
    else
        unsafe_store!(Ptr{typeof(value)}(pointer(data, index)), value)
    end
    record.block_size += data_increase
    return value
end

# Internals
# ---------
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

# Returns (position after tag, type of value)
function loadauxtype(data::Vector{UInt8}, p::Int)
    function auxtype(b)
        return (
            b == UInt8('A') ? Char  :
            b == UInt8('c') ? Int8  :
            b == UInt8('C') ? UInt8 :
            b == UInt8('s') ? Int16 :
            b == UInt8('S') ? UInt16 :
            b == UInt8('i') ? Int32 :
            b == UInt8('I') ? UInt32 :
            b == UInt8('f') ? Float32 :
            b == UInt8('Z') ? String :
            error("invalid type tag: '$(Char(b))'"))
    end
    t = data[p]
    if t == UInt8('B')
        return p + 2, Vector{auxtype(data[p+1])}
    else
        return p + 1, auxtype(t)
    end
end

function loadauxvalue(data::Vector{UInt8}, pos::Int)
    pos, T = loadauxtype(data, pos + 2)
    _, val = loadauxvalue(data, pos, T)
    return val
end

function loadauxvalue(data::Vector{UInt8}, p::Int, ::Type{T}) where T
    return p + sizeof(T), unsafe_load(Ptr{T}(pointer(data, p)))
end

function loadauxvalue(data::Vector{UInt8}, p::Int, ::Type{Char})
    return p + 1, Char(unsafe_load(pointer(data, p)))
end

function loadauxvalue(data::Vector{UInt8}, p::Int, ::Type{Vector{T}}) where T
    n = unsafe_load(Ptr{Int32}(pointer(data, p)))
    p += 4
    xs = Vector{T}(undef, n)
    unsafe_copyto!(pointer(xs), Ptr{T}(pointer(data, p)), n)
    return p + n * sizeof(T), xs
end

function loadauxvalue(data::Vector{UInt8}, p::Int, ::Type{String})
    dataptr = pointer(data, p)
    endptr = ccall(:memchr, Ptr{Cvoid}, (Ptr{Cvoid}, Cint, Csize_t),
                   dataptr, '\0', length(data) - p + 1)
    q::Int = p + (endptr - dataptr) - 1
    return q + 2, String(data[p:q])
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
    typ = data[p+2]
    @inbounds size = TYPESIZES[typ]
    p += 3
    if size == 0
        error("invalid type tag: '$(Char(typ))'")
    elseif size > 0
        p += size
    elseif size == -2 # NULL-terminalted string of hex
        while data[p] != 0x00
            p += 1
        end
        p += 1
    else # Array
        eltyp = data[p]
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
    val = loadauxvalue(record.data, pos)
    nextpos = next_tag_position(record.data, pos)
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
