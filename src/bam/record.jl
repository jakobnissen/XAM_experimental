# BAM Record
# ==========

# the data size of fixed-length fields of the BAM object
const FIXED_FIELDS_BYTES = 36

"""
    BAM.Record()

Create an unfilled BAM record.
"""
mutable struct Record
    # fixed-length fields (see BAM specs for the details)
    block_size::Int32
    refid::Int32
    pos::Int32
    bin_mq_nl::UInt32
    flag_nc::UInt32
    l_seq::Int32
    next_refid::Int32
    next_pos::Int32
    tlen::Int32
    # variable length data
    data::Vector{UInt8}
    reader::Union{Reader, Nothing}

    function Record()
        # An empty record is a legitimate BAM record.
        block_size = FIXED_FIELDS_BYTES - sizeof(Int32) + 1
        flag_nc = UInt32(SAM.FLAG_UNMAP) << 16
        # Per specs, unknown MAPQ is 0xff, unknown length x01 (NULL terminator)
        bin_mq_nl = 0x0000ff01
        # Only the null terminator of the query name sequence
        data = UInt8[0x00]
        return new(block_size, -1, -1, bin_mq_nl, flag_nc, 0, -1, -1, 0, data, nothing)
    end
end

function Base.empty!(record::Record)
    record.block_size = FIXED_FIELDS_BYTES - sizeof(Int32) + 1
    record.refid      = -1
    record.pos        = -1
    record.bin_mq_nl  = 0x0000ff01
    record.flag_nc    = UInt32(SAM.FLAG_UNMAP) << 16
    record.l_seq      = 0
    record.next_refid = -1
    record.next_pos   = -1
    record.tlen       = 0
    record.data[1]    = 0x00
    record.reader     = nothing
    return record
end

function Record(data::Vector{UInt8})
    return convert(Record, data)
end

function Base.convert(::Type{Record}, data::Vector{UInt8})
    length(data) < FIXED_FIELDS_BYTES && throw(ArgumentError("data too short"))
    record = Record()
    dst_pointer = Ptr{UInt8}(pointer_from_objref(record))
    unsafe_copyto!(dst_pointer, pointer(data), FIXED_FIELDS_BYTES)
    dsize = data_size(record)
    length(data) < dsize + FIXED_FIELDS_BYTES && throw(ArgumentError("data too short"))
    resize!(record.data, dsize)
    unsafe_copyto!(record.data, 1, data, FIXED_FIELDS_BYTES + 1, dsize)
    return record
end

function Base.copy(record::Record)
    copy = Record()
    copy.block_size = record.block_size
    copy.refid      = record.refid
    copy.pos        = record.pos
    copy.bin_mq_nl  = record.bin_mq_nl
    copy.flag_nc    = record.flag_nc
    copy.l_seq      = record.l_seq
    copy.next_refid = record.next_refid
    copy.next_pos   = record.next_pos
    copy.tlen       = record.tlen
    copy.data       = record.data[1:data_size(record)]
    copy.reader     = record.reader
    return copy
end

function quality_string(quals::Vector{UInt8})
    characters = Vector{Char}(undef, length(quals))
    for i in eachindex(quals)
        @inbounds qual = quals[i]
        if qual < 10
            char = ' '
        elseif qual < 15
            char = '▁'
        elseif qual < 20
            char = '▂'
        elseif qual < 25
            char = '▃'
        elseif qual < 30
            char = '▄'
        elseif qual < 35
            char = '▆'
        elseif qual < 40
            char = '▇'
        elseif qual < 255
            char = '█'
        else
            char = '?'
        end
        @inbounds characters[i] = char
    end
    return join(characters)
end

function compact_string(sequence)
    LEFT_PADDING = 21
    width = displaysize()[2] - LEFT_PADDING
    if length(sequence) <= width
        return string(sequence)
    else
        half = div(width - 1, 2)
        return string(sequence[1:half]) * '…' * string(sequence[end-half:end])
    end
end

function Base.show(io::IO, record::Record)
    println(io, summary(record), ':')
    print(io,   "      template name: "); show(io, tempname(record))
    print(io, "\n               flag: "); show(io, flag(record))
    print(io, "\n       reference ID: "); show(io, refid(record))
    print(io, "\n           position: "); show(io, position(record))
    print(io, "\n    mapping quality: "); show(io, mappingquality(record))
    print(io, "\n              CIGAR: "); show(io, cigar(record))
    print(io, "\n  next reference ID: "); show(io, nextrefid(record))
    print(io, "\n      next position: "); show(io, nextposition(record))
    print(io, "\n    template length: "); show(io, templength(record))
    # Sequence and base quality
    LEFT_PADDING = 21
    width = displaysize()[2] - LEFT_PADDING
    seq = sequence(record)
    if seq === nothing
        print(io, "\n           sequence: "); show(io,  seq)
        print(io, "\n       base quality: "); show(io, nothing)
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
    n_aux_bytes = data_size(record) - auxdata_position(record) + 1
    if n_aux_bytes > 500
        print(io, "<", n_aux_bytes, " bytes auxiliary data>")
    else
        for field in keys(auxdata(record))
            print(io, ' ', field, '='); show(record[field])
        end
    end
end

function Base.read!(reader::Reader, record::Record)
    return _read!(reader, record)
end

# ============== Accessor functions =====================
"""
    flag(record::Record)::UInt16

Get the bitwise flag of `record`.
"""
function flag(record::Record)::UInt16
    return UInt16(record.flag_nc >> 16)
end

"""
    refid(record::Record)

Get the reference sequence ID of `record`.

The ID is 1-based (i.e. the first sequence is 1). Returns `nothing` if not defined.
Note that unmapped records can still have a reference.
"""
# According to BAM specs, unmapped records CAN have reference and position
function refid(record::Record)::Union{Int,Nothing}
    return record.refid == -1 ? nothing : Int(record.refid + 1)
end

# Return the length of the read name. NOT including the null byte terminator
function seqname_length(record::Record)::UInt8
    return UInt8(record.bin_mq_nl & 0xff) - 0x01
end

"""
    nextrefid(record::Record)

Get the next/mate reference sequence ID of `record`.

The ID is 1-based (i.e. the first sequence is 1) and is `nothing` for an unpaired
record. Note that unmapped records can still have a reference.
"""
function nextrefid(record::Record)::Union{Int,Nothing}
    ispaired = flag(record) & SAM.FLAG_PAIRED == SAM.FLAG_PAIRED
    if !ispaired || record.next_refid == -1
        return nothing
    end
    return  Int(record.next_refid + 1)
end

"""
    refname(record::Record)::String

Get the reference sequence name of `record`.

See also: `BAM.refid`.
"""
function refname(record::Record)::Union{String,Nothing}
    id = refid(record)
    if id === nothing || record.reader === nothing
        return nothing
    end
    return record.reader.refseqnames[id]
end

"""
    nextrefname(record::Record)::String

Get the reference name of the mate/next read of `record`.

See also: `BAM.nextrefid`
"""
function nextrefname(record::Record)::Union{String,Nothing}
    id = nextrefid(record)
    if id === nothing || record.reader === nothing
        return nothing
    end
    return record.reader.refseqnames[id]
end

"""
    reflen(record::Record)::Int

Get the length of the reference sequence this record applies to. Note that
unmapped records can still have a reference.
"""
function reflen(record::Record)::Union{Int,Nothing}
    id = refid(record)
    if id === nothing || record.reader === nothing
        return nothing
    end
    return record.reader.refseqlens[id]
end

"""
    position(record::Record)

Get the 1-based leftmost mapping position of `record`. Note that unmapped records
can still have a reference position.
"""
function position(record::Record)::Union{Int,Nothing}
    return record.pos == -1 ? nothing : Int(record.pos + 1)
end

"""
    nextposition(record::Record)::Int

Get the 1-based leftmost mapping position of the next/mate read of `record`.
Returns `nothing` if the read is unmapped.
"""
function nextposition(record::Record)::Union{Int,Nothing}
    ispaired = flag(record) & SAM.FLAG_PAIRED == SAM.FLAG_PAIRED
    if !ispaired || record.next_pos == -1
        return nothing
    end
    return Int(record.next_pos + 1)
end

"""
    rightposition(record::Record)::Int

Get the 1-based rightmost mapping position of `record`.
"""
function rightposition(record::Record)::Union{Int,Nothing}
    pos = position(record)
    return pos === nothing ? nothing : pos + alignlength(record) - 1
end

"""
    mappingquality(record::Record)::UInt8

Get the mapping quality of `record`. Returns nothing if the MAPQ field is set to
0xff, indicating unknown value.
"""
function mappingquality(record::Record)::Union{UInt8,Nothing}
    # Mapping qual of 0xff means unknown as per BAM specs
    as_uint8 = UInt8((record.bin_mq_nl >> 8) & 0xff)
    return as_uint8 == 0xff ? nothing : as_uint8
end

"""
    tempname(record::Record)::String

Get the query template name of `record`.
"""
function tempname(record::Record)::Union{String,Nothing}
    seqlen = seqname_length(record)
    return seqlen == 0 ? nothing : unsafe_string(pointer(record.data), seqlen)
end

"""
    seqlength(record::Record)::Int

Get the sequence length of `record`.
"""
function seqlength(record::Record)::Int
    return Int(record.l_seq)
end

"""
    templength(record::Record)

Get the template length of `record`. Returns `nothing` when unavailable.
"""
function templength(record::Record)::Union{Int,Nothing}
    return record.tlen == 0 ? nothing : Int(record.tlen)
end

# ================= Boolean functions ===================
"""
    ismapped(record::Record)::Bool

Test if `record` is mapped.
"""
function ismapped(record::Record)::Bool
    return flag(record) & SAM.FLAG_UNMAP == 0
end

"""
    isprimary(record::Record)::Bool

Test if `record` is a primary line of the read.

This is equivalent to `flag(record) & 0x900 == 0`.
"""
function isprimary(record::Record)::Bool
    return flag(record) & (SAM.FLAG_SECONDARY | SAM.FLAG_SUPPLEMENTARY) == 0
end

"""
    ispositivestrand(record::Record)::Bool

Test if `record` is aligned to the positive strand.

This is equivalent to `flag(record) & 0x10 == 0`.
"""
function ispositivestrand(record::Record)::Bool
    return flag(record) & SAM.FLAG_REVERSE == 0
end

"""
    isnextmapped(record::Record)::Bool

Test if the mate/next read of `record` is mapped.
"""
function isnextmapped(record::Record)::Bool
    return flag(record) & (SAM.FLAG_MUNMAP | SAM.FLAG_PAIRED) == 1
end

# ==============   SAM-compatibility functions ==========
function hastempname(record::Record)
    return seqname_length(record) > 0
end

# TODO: Update this - do we need it? What does the SAM equivalent do?
# Review this once I've rewritten SAM
function hasseqlength(record::Record)
    return true
end

function hassequence(record::Record)
    return seqlength(record) > 0
end

function hasrefid(record::Record)
    return record.refid > -1
end

function hasrefname(record::Record)
    return record.refid > -1
end

function hasposition(record::Record)
    return record.pos > -1
end

function hasrightposition(record::Record)
    return record.pos > -1
end

function hasmappingquality(record::Record)
    return mappingquality(x) !== nothing
end

function hasalignment(record::Record)
    return cigar_rle(record) !== nothing
end

function hasnextrefid(record::Record)
    return record.next_refid > -1
end

function hasnextrefname(record::Record)
    return record.next_refid > -1
end

function hasnextposition(record::Record)
    return record.next_pos > -1
end

function hastemplength(record::Record)
    return record.tlen != 0
end

function hasquality(record::Record)
    quals = quality(record)
    # undefined quals are 0xff per specs
    return quals !== nothing && any(i != 0xff for i in quality(record))
end

function hasauxdata(record::Record)
    return auxdata_position(record) < data_size(record)
end

# ================= BioCore functions =========================

function BioCore.isfilled(record::Record)
    return true
end

function BioCore.seqname(record::Record)
    return tempname(record)
end

function BioCore.hasseqname(record::Record)
    return hastempname(record)
end

function BioCore.sequence(record::Record)
    return sequence(record)
end

function BioCore.hassequence(record::Record)
    return hassequence(record)
end

function BioCore.leftposition(record::Record)
    return position(record)
end

function BioCore.hasleftposition(record::Record)
    return hasposition(record)
end

function BioCore.rightposition(record::Record)
    return rightposition(record)
end

function BioCore.hasrightposition(record::Record)
    return hasrightposition(record)
end

# =============== Helper functions ============================

# Return the size of the `.data` field.
function data_size(record::Record)
    return record.block_size - FIXED_FIELDS_BYTES + sizeof(record.block_size)
end

function auxdata_position(record::Record)
    seqlen = seqlength(record)
    #        template_name        null_term     n_cigar_op      per_cigar_op    seq         qual
    offset = seqname_length(record) + 1 + n_cigar_op(record, false) * 4 + cld(seqlen, 2) + seqlen
    return offset + 1
end


#######################################################



"""
    n_cigar_op(record::Record, checkCG::Bool = true)

Return the number of operations in the CIGAR string of `record`.

Note that in the BAM specification, the field called `cigar` typically stores
the cigar string of the record.
However, this is not always true, sometimes the true cigar is very long,
and due to  some constraints of the BAM format, the actual cigar string is
stored in an extra tag: `CG:B,I`, and the `cigar` field stores a pseudo-cigar
string.

Calling this method with `checkCG` set to `true` (default) this method will
always yield the number of operations in the true cigar string, because this is
probably what you want, the vast majority of the time.

If you have a record that stores the true cigar in a `CG:B,I` tag, but you still
want to get the number of operations in the `cigar` field of the BAM record,
then set `checkCG` to `false`.
"""
function n_cigar_op(record::Record, checkCG::Bool = true)
    return cigar_position(record, checkCG)[2]
end

"""
    cigar(record::Record)::String

Get the CIGAR string of `record`.

Note that in the BAM specification, the field called `cigar` typically stores
the cigar string of the record.
However, this is not always true, sometimes the true cigar is very long,
and due to  some constraints of the BAM format, the actual cigar string is
stored in an extra tag: `CG:B,I`, and the `cigar` field stores a pseudo-cigar
string.

Calling this method with `checkCG` set to `true` (default) this method will
always yield the true cigar string, because this is probably what you want
the vast majority of the time.

If you have a record that stores the true cigar in a `CG:B,I` tag, but you still
want to access the pseudo-cigar that is stored in the `cigar` field of the BAM
record, then you can set checkCG to `false`.

See also `BAM.cigar_rle`.
"""
function cigar(record::Record, checkCG::Bool = true)::Union{String,Nothing}
    rle = cigar_rle(record, checkCG)
    if rle === nothing
        return nothing
    else
        buf = IOBuffer()
        for (op, len) in zip(rle...)
            print(buf, len, convert(Char, op))
        end
        return String(take!(buf))
    end
end

"""
    cigar_rle(record::Record, checkCG::Bool = true)::Tuple{Vector{BioAlignments.Operation},Vector{Int}}

Get a run-length encoded tuple `(ops, lens)` of the CIGAR string in `record`.

Note that in the BAM specification, the field called `cigar` typically stores
the cigar string of the record.
However, this is not always true, sometimes the true cigar is very long,
and due to  some constraints of the BAM format, the actual cigar string is
stored in an extra tag: `CG:B,I`, and the `cigar` field stores a pseudo-cigar
string.

Calling this method with `checkCG` set to `true` (default) this method will
always yield the true cigar string, because this is probably what you want
the vast majority of the time.

If you have a record that stores the true cigar in a `CG:B,I` tag, but you still
want to access the pseudo-cigar that is stored in the `cigar` field of the BAM
record, then you can set checkCG to `false`.

See also `BAM.cigar`.
"""
function cigar_rle(record::Record, checkCG::Bool = true)::Union{Tuple{Vector{BioAlignments.Operation},Vector{Int}},Nothing}
    idx, nops = cigar_position(record, checkCG)
    if nops == 0
        return nothing
    else
        ops, lens = extract_cigar_rle(record.data, idx, nops)
        return ops, lens
    end
end

function extract_cigar_rle(data::Vector{UInt8}, offset, n)
    ops = Vector{BioAlignments.Operation}()
    lens = Vector{Int}()
    for i in offset:4:offset + (n - 1) * 4
        x = unsafe_load(Ptr{UInt32}(pointer(data, i)))
        op = BioAlignments.Operation(x & 0x0F)
        push!(ops, op)
        push!(lens, x >> 4)
    end
    return ops, lens
end

function cigar_position(record::Record, checkCG::Bool = true)::Tuple{Int, Int}
    cigaridx, nops = seqname_length(record) + 2, record.flag_nc & 0xFFFF
    if !checkCG
        return cigaridx, nops
    end
    if nops != 2
        return cigaridx, nops
    end
    x = unsafe_load(Ptr{UInt32}(pointer(record.data, cigaridx)))
    if x != UInt32(seqlength(record) << 4 | 4)
        return cigaridx, nops
    end
    start = auxdata_position(record)
    stop = data_size(record)
    tagidx = findauxtag(record.data, start, stop, UInt8('C'), UInt8('G'))
    if tagidx == 0
        return cigaridx, nops
    end
    # Tag exists, validate type is BI.
    typ = unsafe_load(Ptr{UInt16}(pointer(record.data, tagidx += 2)))
    if typ != (UInt16('I') << 8 | UInt16('B'))
        return cigaridx, nops
    end
    # If got this far, the CG tag is valid and contains the cigar.
    # Get the true n_cigar_ops, and return it and the idx of the first
    nops = UInt32(unsafe_load(Ptr{Int32}(pointer(record.data, tagidx += 2))))
    tagidx += 4
    return tagidx, nops
end

"""
    alignment(record::Record)::BioAlignments.Alignment

Get the alignment of `record`.
"""
function alignment(record::Record)::BioAlignments.Alignment
    if !ismapped(record)
        return BioAlignments.Alignment(BioAlignments.AlignmentAnchor[])
    end
    seqpos = 0
    refpos = position(record) - 1
    anchors = [BioAlignments.AlignmentAnchor(seqpos, refpos, BioAlignments.OP_START)]
    for (op, len) in zip(cigar_rle(record)...)
        if BioAlignments.ismatchop(op)
            seqpos += len
            refpos += len
        elseif BioAlignments.isinsertop(op)
            seqpos += len
        elseif BioAlignments.isdeleteop(op)
            refpos += len
        else
            error("operation $(op) is not supported")
        end
        push!(anchors, BioAlignments.AlignmentAnchor(seqpos, refpos, op))
    end
    return BioAlignments.Alignment(anchors)
end

"""
    alignlength(record::Record)::Int

Get the alignment length of `record`.
"""
function alignlength(record::Record)::Int
    offset = seqname_length(record) + 1
    length::Int = 0
    for i in offset + 1:4:offset + n_cigar_op(record, false) * 4
        x = unsafe_load(Ptr{UInt32}(pointer(record.data, i)))
        op = BioAlignments.Operation(x & 0x0F)
        if BioAlignments.ismatchop(op) || BioAlignments.isdeleteop(op)
            length += x >> 4
        end
    end
    return length
end

"""
    sequence(record::Record)::BioSequences.DNASequence

Get the segment sequence of `record`. Returns `nothing` if not available.
"""
function sequence(record::Record)::Union{BioSequences.DNASequence, Nothing}
    seqlen = seqlength(record)
    if seqlen == 0
        return nothing
    else
        data = Vector{UInt64}(undef, cld(seqlen, 16))
        src::Ptr{UInt64} = pointer(record.data, seqname_length(record) + n_cigar_op(record, false) * 4 + 2)
        for i in 1:lastindex(data)
            # copy data flipping high and low nybble
            x = unsafe_load(src, i)
            data[i] = (x & 0x0f0f0f0f0f0f0f0f) << 4 | (x & 0xf0f0f0f0f0f0f0f0) >> 4
        end
        return BioSequences.DNASequence(data, 1:seqlen, false)
    end
end

"""
    quality(record::Record)::Vector{UInt8}

Get the base quality of `record`, or `nothing` if not available.
"""
function quality(record::Record)::Union{Vector{UInt8}, Nothing}
    seqlen = seqlength(record)
    if seqlen == 0
        return nothing
    else
        offset = seqname_length(record) + 1 + n_cigar_op(record, false) * 4 + cld(seqlen, 2)
        return [reinterpret(UInt8, record.data[i+offset]) for i in 1:seqlen]
    end
end
