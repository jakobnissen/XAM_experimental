# SAM Record
# ==========

const FIXED_FIELD_BYTES = 192

mutable struct SAMRecord <: XAMRecord
    # Indices
    qname::UnitRange{Int}
    flag::UnitRange{Int}
    rname::UnitRange{Int}
    pos::UnitRange{Int}
    mapq::UnitRange{Int}
    cigar::UnitRange{Int}
    rnext::UnitRange{Int}
    pnext::UnitRange{Int}
    tlen::UnitRange{Int}
    seq::UnitRange{Int}
    qual::UnitRange{Int}
    filled::UnitRange{Int}
    # data
    data::Vector{UInt8}
    fields::Vector{UnitRange{Int}}
end

"""
    SAM.SAMRecord()

Create an unfilled SAM record.
"""
function SAMRecord()
    return SAMRecord(
        # qname-mapq
        1:0, 1:0, 1:0, 1:0, 1:0,
        # cigar-seq
        1:0, 1:0, 1:0, 1:0, 1:0,
        # qual and fields
        1:0, 1:0,
        UInt8[],
        UnitRange{Int}[])
end

"""
    SAM.SAMRecord(data::Vector{UInt8})

Create a SAM record from `data`.
This function verifies the format and indexes fields for accessors.
Note that the ownership of `data` is transferred to a new record object.
"""
function SAMRecord(data::Vector{UInt8})
    return convert(SAMRecord, data)
end

function Base.convert(::Type{SAMRecord}, data::Vector{UInt8})
    return SAMRecord(
        # qname-mapq
        1:0, 1:0, 1:0, 1:0, 1:0,
        # cigar-seq
        1:0, 1:0, 1:0, 1:0, 1:0,
        # qual and fields
        1:0, 1:0,
        data,
        UnitRange{Int}[])
    index!(record)
    return record
end

"""
    SAM.SAMRecord(str::AbstractString)

Create a SAM record from `str`.
This function verifies the format and indexes fields for accessors.
"""
function SAMRecord(str::AbstractString)
    return convert(SAMRecord, str)
end

function Base.convert(::Type{SAMRecord}, str::AbstractString)
    return SAMRecord(Vector{UInt8}(str))
end

function Base.show(io::IO, record::SAMRecord)
    print(io, summary(record), ':')
    println(io)
    print(io,   "      template name: "); show(io, tempname(record))
    print(io, "\n               flag: "); show(io, flag(record))
    print(io, "\n          reference: "); show(io, refname(record))
    print(io, "\n           position: "); show(io, position(record))
    print(io, "\n    mapping quality: "); show(io, mappingquality(record))
    print(io, "\n              CIGAR: "); show(io, cigar(record))
    print(io, "\n     next reference: "); show(io, nextrefname(record))
    print(io, "\n      next position: "); show(io, nextposition(record))
    print(io, "\n    template length: "); show(io, templength(record))
    # Sequence and base quality
    LEFT_PADDING = 20
    width = displaysize()[2] - LEFT_PADDING
    seq = sequence(record)
    if seq === missing
        print(io, "\n           sequence: "); show(io,  seq)
        print(io, "\n       base quality: "); show(io, missing)
    elseif length(seq) <= width
        print(io, "\n           sequence: "); show(io,  seq)
        print(io, "\n       base quality: "); print(io, quality_string(quality(record)))
    else
        half = div(width, 2) - 1
        print(io, "\n           sequence: ", seq[1:half], '…', seq[end-half:end])
        qual = quality(record)
        print(io, "\n       base quality: ", quality_string(qual[1:half]), '…',
                  quality_string(qual[end-half:end]))
    end
    # Auxiliary fields - don't show if too long
    print(io, "\n     auxiliary data:")
    n_aux_bytes = 0
    for field in record.fields
        n_aux_bytes += length(field)
    end
    if n_aux_bytes > 500
        print(io, "<", n_aux_bytes, " bytes auxiliary data>")
    else
        for field in keys(auxdata(record))
            print(io, ' ', field, '='); show(record[field])
        end
    end
end

function Base.print(io::IO, record::SAMRecord)
    write(io, record)
    return nothing
end

function Base.write(io::IO, record::SAMRecord)
    return unsafe_write(io, pointer(record.data, first(record.filled)), length(record.filled))
end

function Base.copy(record::SAMRecord)
    return SAMRecord(
        record.qname,
        record.flag,
        record.rname,
        record.pos,
        record.mapq,
        record.cigar,
        record.rnext,
        record.pnext,
        record.tlen,
        record.seq,
        record.qual,
        record.filled,
        copy(record.data),
        copy(record.fields))
end

# ============== Accessor functions =====================---------
"""
    flag(record::SAMRecord)::UInt16

Get the bitwise flag of `record`.
"""
function flag(record::SAMRecord)::UInt16
    return isempty(record.flag) ? 0x004 : unsafe_parse_decimal(UInt16, record.data, record.flag)
end

# ================= Boolean functions ===================


"""
    ismapped(record::SAMRecord)::Bool

Test if `record` is mapped.
"""
function ismapped(record::SAMRecord)::Bool
    return flag(record) & FLAG_UNMAP == 0
end

"""
    isprimary(record::SAMRecord)::Bool

Test if `record` is a primary line of the read.

This is equivalent to `flag(record) & 0x900 == 0`.
"""
function isprimary(record::SAMRecord)::Bool
    return flag(record) & 0x900 == 0
end

"""
    refname(record::SAMRecord)::String

Get the reference sequence name of `record`.
"""
function refname(record::SAMRecord)
    if isempty(record.rname)
        return missing
    end
    name = record.data[record.rname]
    if isstar(name)
        return missing
    else
        return String(record.data[record.rname])
    end
end

"""
    position(record::SAMRecord)::Int

Get the 1-based leftmost mapping position of `record`.
"""
function position(record::SAMRecord)::Union{Int, Missing}
    if isempty(record.pos)
        return missing
    end
    pos = unsafe_parse_decimal(Int, record.data, record.pos)
    return pos == 0 ? missing : pos
end

"""
    rightposition(record::SAMRecord)::Int

Get the 1-based rightmost mapping position of `record`.
"""
function rightposition(record::SAMRecord)
    return position(record) + alignlength(record) - 1
end


"""
    isnextmapped(record::SAMRecord)::Bool

Test if the mate/next read of `record` is mapped.
"""
function isnextmapped(record::SAMRecord)::Bool
    return flag(record) & FLAG_MUNMAP == 0
end

"""
    nextrefname(record::SAMRecord)::String

Get the reference name of the mate/next read of `record`.
"""
function nextrefname(record::SAMRecord)::Union{String, Missing}
    if isempty(record.rname)
        return missing
    end
    name = record.data[record.rname]
    if isstar(name)
        return missing
    elseif @inbounds sizeof(name) == 1 && name[1] == UInt8('=')
        return refname(record)
    else
        return String(record.data[record.rname])
    end
end

"""
    nextposition(record::SAMRecord)::Int

Get the position of the mate/next read of `record`.
"""
function nextposition(record::SAMRecord)::Union{Int, Missing}
    if isempty(record.pnext)
        return missing
    end
    pos = unsafe_parse_decimal(Int, record.data, record.pnext)
    return pos == 0 ? missing : pos
end

"""
    mappingquality(record::SAMRecord)::UInt8

Get the mapping quality of `record`.
"""
function mappingquality(record::SAMRecord)::Union{UInt8, Missing}
    if isempty(record.mapq)
        return missing
    end
    qual = unsafe_parse_decimal(UInt8, record.data, record.mapq)
    return qual == 0xff ? missing : qual
end

"""
    cigar(record::SAMRecord)::String

Get the CIGAR string of `record`.
"""
function cigar(record::SAMRecord)::Union{String, Missing}
    if isempty(record.cigar)
        return missing
    end
    cig = record.data[record.cigar]
    if isstar(cig)
        return missing
    else
        return String(record.data[record.cigar])
    end
end

"""
    alignment(record::SAMRecord)::BioAlignments.Alignment

Get the alignment of `record`.
"""
function alignment(record::SAMRecord)::BioAlignments.Alignment
    if ismapped(record)
        return BioAlignments.Alignment(cigar(record), 1, position(record))
    else
        return BioAlignments.Alignment(BioAlignments.AlignmentAnchor[])
    end
end

"""
    alignlength(record::SAMRecord)::Int

Get the alignment length of `record`.
"""
function alignlength(record::SAMRecord)::Int
    if length(record.cigar) == 1 && record.data[first(record.cigar)] == UInt8('*')
        return 0
    end
    ret::Int = 0
    len = 0  # operation length
    for i in record.cigar
        c = record.data[i]
        if c ∈ UInt8('0'):UInt8('9')
            len = len * 10 + (c - UInt8('0'))
        else
            op = convert(BioAlignments.Operation, Char(c))
            if BioAlignments.ismatchop(op) || BioAlignments.isdeleteop(op)
                ret += len
                len = 0
            end
        end
    end
    return ret
end

"""
    tempname(record::SAMRecord)::String

Get the query template name of `record`.
"""
function tempname(record::SAMRecord)::Union{String, Missing}
    if isempty(record.qname)
        return missing
    end
    name = record.data[record.qname]
    if isstar(name)
        return missing
    else
        return String(name)
    end
end

"""
    templength(record::SAMRecord)::Int

Get the template length of `record`.
"""
function templength(record::SAMRecord)::Union{Int,Missing}
    if isempty(record.tlen)
        return missing
    end
    len = unsafe_parse_decimal(Int, record.data, record.tlen)
    return len == 0 ? missing : len
end

"""
    sequence(record::SAMRecord)::BioSequences.DNASequence

Get the segment sequence of `record`.
"""
function sequence(record::SAMRecord)::Union{BioSequences.DNASequence, Missing}
    if isempty(record.seq)
        return missing
    end
    seqdata = record.data[record.seq]
    if isstar(seqdata)
        return missing
    else
        seqlen = length(record.seq)
        ret = BioSequences.DNASequence(seqlen)
        BioSequences.encode_copy!(ret, 1, record.data, first(record.seq), seqlen)
        return ret
    end
end

"""
    sequence(::Type{String}, record::SAMRecord)::String

Get the segment sequence of `record` as `String`.
"""
function sequence(::Type{String}, record::SAMRecord)::String
    return String(record.data[record.seq])
end

"""
    seqlength(record::SAMRecord)::Int

Get the sequence length of `record`.
"""
function seqlength(record::SAMRecord)::Union{Int, Missing}
    if isempty(record.seq)
        return missing
    end
    return length(record.seq)
end

"""
    quality(record::SAMRecord)::Vector{UInt8}

Get the Phred-scaled base quality of `record`.
"""
function quality(record::SAMRecord)::Union{Vector{UInt8}, Missing}
    if isempty(record.qual)
        return missing
    end
    qual = record.data[record.qual]
    if isstar(qual)
        return missing
    else
        for i in 1:lastindex(qual)
            @inbounds qual[i] -= 33
        end
        return qual
    end
end

"""
    quality(::Type{String}, record::SAMRecord)::String

Get the ASCII-encoded base quality of `record`.
"""
function quality(::Type{String}, record::SAMRecord)::Union{String,Missing}
    if isempty(record.qual)
        return missing
    end
    return String(record.data[record.qual])
end

"""
    auxdata(record::SAMRecord)::Dict{String,Any}

Get the auxiliary data (optional fields) of `record`.
"""
function auxdata(record::SAMRecord)::Dict{String,Any}
    return Dict(k => record[k] for k in keys(record))
end

function Base.haskey(record::SAMRecord, tag::AbstractString)
    return findauxtag(record, tag) > 0
end

function Base.getindex(record::SAMRecord, tag::AbstractString)
    i = findauxtag(record, tag)
    if i == 0
        throw(KeyError(tag))
    end
    field = record.fields[i]
    # data type
    typ = record.data[first(field)+3]
    lo = first(field) + 5
    if i == lastindex(record.fields)
        hi = last(field)
    else
        hi = first(record.fields[i+1]) - 2
    end
    if typ == UInt8('A')
        @assert lo == hi
        return Char(record.data[lo])
    elseif typ == UInt8('i')
        return unsafe_parse_decimal(Int, record.data, lo:hi)
    elseif typ == UInt8('f')
        # TODO: Call a C function directly for speed?
        return parse(Float32, SubString(record.data[lo:hi]))
    elseif typ == UInt8('Z')
        return String(record.data[lo:hi])
    elseif typ == UInt8('H')
        return parse_hexarray(record.data, lo:hi)
    elseif typ == UInt8('B')
        return parse_typedarray(record.data, lo:hi)
    else
        throw(ArgumentError("type code '$(Char(typ))' is not defined"))
    end
end

function Base.keys(record::SAMRecord)
    return [String(record.data[first(f):first(f)+1]) for f in record.fields]
end

function Base.values(record::SAMRecord)
    return [record[k] for k in keys(record)]
end


# Bio Methods
# -----------
function BioCore.isfilled(record::SAMRecord)
    return !isempty(record.filled)
end

function BioCore.seqname(record::SAMRecord)
    return tempname(record)
end

function BioCore.hasseqname(record::SAMRecord)
    return hastempname(record)
end

function BioCore.sequence(record::SAMRecord)
    return sequence(record)
end

function BioCore.hassequence(record::SAMRecord)
    return hassequence(record)
end

function BioCore.rightposition(record::SAMRecord)
    return rightposition(record)
end

function BioCore.hasrightposition(record::SAMRecord)
    return hasrightposition(record)
end

function BioCore.leftposition(record::SAMRecord)
    return position(record)
end

function BioCore.hasleftposition(record::SAMRecord)
    return hasposition(record)
end


# Helper Functions
# ----------------

function initialize!(record::SAMRecord)
    record.filled = 1:0
    record.qname = 1:0
    record.flag = 1:0
    record.rname = 1:0
    record.pos = 1:0
    record.mapq = 1:0
    record.cigar = 1:0
    record.rnext = 1:0
    record.pnext = 1:0
    record.tlen = 1:0
    record.seq = 1:0
    record.qual = 1:0
    empty!(record.fields)
    return record
end

function findauxtag(record::SAMRecord, tag::AbstractString)
    if sizeof(tag) != 2
        return 0
    end
    t1, t2 = UInt8(tag[1]), UInt8(tag[2])
    for (i, field) in enumerate(record.fields)
        p = first(field)
        if record.data[p] == t1 && record.data[p+1] == t2
            return i
        end
    end
    return 0
end

function parse_hexarray(data::Vector{UInt8}, range::UnitRange{Int})
    @assert iseven(length(range))
    ret = Vector{UInt8}(length(range) >> 1)
    byte2hex(b) = b ∈ 0x30:0x39 ? (b - 0x30) : b ∈ 0x41:0x46 ? (b - 0x41 + 0x0A) : error("not in [0-9A-F]")
    j = 1
    for i in first(range):2:last(range)-1
        ret[j] = (byte2hex(data[range[i]]) << 4) | byte2hex(data[range[i+1]])
        j += 1
    end
    return ret
end

function parse_typedarray(data::Vector{UInt8}, range::UnitRange{Int})
    # format: [cCsSiIf](,[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)+
    t = data[first(range)]
    xs = split(String(data[first(range)+2:last(range)]))
    if t == UInt8('c')
        return [parse(Int8, x) for x in xs]
    elseif t == UInt8('C')
        return [parse(UInt8, x) for x in xs]
    elseif t == UInt8('s')
        return [parse(Int16, x) for x in xs]
    elseif t == UInt8('S')
        return [parse(UInt16, x) for x in xs]
    elseif t == UInt8('i')
        return [parse(Int32, x) for x in xs]
    elseif t == UInt8('I')
        return [parse(UInt32, x) for x in xs]
    elseif t == UInt8('f')
        return [parse(Float32, x) for x in xs]
    else
        throw(ArgumentError("type code '$(Char(t))' is not defined"))
    end
end

function isstar(x::Vector{UInt8})
    @inbounds return sizeof(x) == 1 && x[1] == UInt8('*')
end
