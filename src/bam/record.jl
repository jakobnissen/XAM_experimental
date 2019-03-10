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
    reader::Reader

    function Record()
        # An empty record is a legitimate BAM record.
        block_size = FIXED_FIELDS_BYTES - sizeof(Int32)
        flag_nc = UInt32(SAM.FLAG_UNMAP) << 16
        # Unknown mapping quality is set to 0xff as per specs
        bin_mq_nl = 0x0000ff00
        return new(block_size, -1, -1, bin_mq_nl, flag_nc, 0, -1, -1, 0, UInt8[])
    end
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
    resize!(record.data, dsize)
    length(data) < dsize + FIXED_FIELDS_BYTES && throw(ArgumentError("data too short"))
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
    if isdefined(record, :reader)
        copy.reader = record.reader
    end
    return copy
end

function Base.show(io::IO, record::Record)
    _refid = refid(record) === nothing ? "" : refid(record)
    _position = position(record) === nothing ? "" : position(record)
    _mappingquality = mappingquality(record) === nothing ? "" : mappingquality(record)
    _nextrefid = nextrefid(record) === nothing ? "" : nextrefid(record)
    _nextposition = nextposition(record) === nothing ? "" : nextposition(record)

    print(io, summary(record), ':')
    println(io)
    println(io, "      template name: ", tempname(record))
    println(io, "               flag: ", flag(record))
    println(io, "       reference ID: ", _refid)
    println(io, "           position: ", _position)
    println(io, "    mapping quality: ", _mappingquality)
    println(io, "              CIGAR: ", cigar(record))
    println(io, "  next reference ID: ", _nextrefid)
    println(io, "      next position: ", _nextposition)
    println(io, "    template length: ", templength(record))
    println(io, "           sequence: ", sequence(record))
    # TODO: pretty print base quality
    println(io, "       base quality: ", quality(record))
      print(io, "     auxiliary data:")
    for field in keys(auxdata(record))
        print(io, ' ', field, '=', record[field])
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
    refid(record::Record)::Int

Get the reference sequence ID of `record`.

The ID is 1-based (i.e. the first sequence is 1). Note that unmapped records can
still have a reference.
"""
# According to BAM specs, unmapped records CAN have reference and position
function refid(record::Record)::Int
    return Int(record.refid + 1)
end

# Return the length of the read name.
function seqname_length(record::Record)::UInt8
    return UInt8(record.bin_mq_nl & 0xff)
end

"""
    nextrefid(record::Record)

Get the next/mate reference sequence ID of `record`.

The ID is 1-based (i.e. the first sequence is 1) and is `nothing` for an unpaired
record. Note that unmapped records can still have a reference.
"""
function nextrefid(record::Record)::Union{Int,Nothing}
    ispaired = flag(record) & SAM.FLAG_PAIRED == SAM.FLAG_PAIRED
    return ispaired ? Int(record.next_refid + 1) : nothing
end

"""
    refname(record::Record)::String

Get the reference sequence name of `record`.

See also: `BAM.refid`.
"""
function refname(record::Record)::String
    id = refid(record)
    check_refid(id)
    return record.reader.refseqnames[id]
end

"""
    nextrefname(record::Record)::String

Get the reference name of the mate/next read of `record`.

See also: `BAM.nextrefid`
"""
function nextrefname(record::Record)::String
    check_paired(record)
    id = nextrefid(record)
    check_refid(id)
    return record.reader.refseqnames[id]
end

"""
    reflen(record::Record)::Int

Get the length of the reference sequence this record applies to. Note that
unmapped records can still have a reference.
"""
function reflen(record::Record)::Int
    id = refid(record)
    check_refid(id)
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
function nextposition(record::Record)
    ispaired = flag(record) & SAM.FLAG_PAIRED == SAM.FLAG_PAIRED
    return ispaired ? Int(record.next_pos + 1) : nothing
end

"""
    rightposition(record::Record)::Int

Get the 1-based rightmost mapping position of `record`.
"""
function rightposition(record::Record)::Int
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
function tempname(record::Record)::String
    return unsafe_string(pointer(record.data), max(seqname_length(record) - 1, 0))
end

"""
    seqlength(record::Record)::Int

Get the sequence length of `record`.
"""
function seqlength(record::Record)::Int
    return Int(record.l_seq)
end

"""
    templength(record::Record)::Int

Get the template length of `record`.
"""
function templength(record::Record)::Int
    return Int(record.tlen)
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
    return ismapped(record)
end

function hasrefname(record::Record)
    return ismapped(record)
end

function hasposition(record::Record)
    return ismapped(record)
end

function hasrightposition(record::Record)
    return ismapped(record)
end

function hasmappingquality(record::Record)
    return ismapped(record)
end

function hasalignment(record::Record)
    return ismapped(record)
end

function hasnextrefid(record::Record)
    return isnextmapped(record)
end

function hasnextrefname(record::Record)
    return isnextmapped(record)
end

function hasnextposition(record::Record)
    return isnextmapped(record)
end

function hastemplength(record::Record)
    return record.tlen != 0
end

function hasquality(record::Record)
    return any(quality(record) .!= 0xff) # undefined quals are 0xff per specs
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

function check_refid(refid)
    if refid < 1
        throw(ArgumentError("Invalid refid"))
    end
end

function check_paired(record)
    if flag(record) & SAM.FLAG_PAIRED == 0
        throw(ArgumentError("Unpaired record"))
    end
end

function check_mapped(record)
    if !ismapped(record)
        throw(ArgumentError("Record is not mapped"))
    end
end

# Return the size of the `.data` field.
function data_size(record::Record)
    return record.block_size - FIXED_FIELDS_BYTES + sizeof(record.block_size)
end

function auxdata_position(record::Record)
    seqlen = seqlength(record)
    return seqname_length(record) + n_cigar_op(record, false) * 4 + cld(seqlen, 2) + seqlen + 1
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
function cigar(record::Record, checkCG::Bool = true)::String
    buf = IOBuffer()
    for (op, len) in zip(cigar_rle(record, checkCG)...)
        print(buf, len, convert(Char, op))
    end
    return String(take!(buf))
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
function cigar_rle(record::Record, checkCG::Bool = true)::Tuple{Vector{BioAlignments.Operation},Vector{Int}}
    idx, nops = cigar_position(record, checkCG)
    ops, lens = extract_cigar_rle(record.data, idx, nops)
    return ops, lens
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
    cigaridx, nops = seqname_length(record) + 1, record.flag_nc & 0xFFFF
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
    offset = seqname_length(record)
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

Get the segment sequence of `record`.
"""
function sequence(record::Record)::BioSequences.DNASequence
    seqlen = seqlength(record)
    data = Vector{UInt64}(undef, cld(seqlen, 16))
    src::Ptr{UInt64} = pointer(record.data, seqname_length(record) + n_cigar_op(record, false) * 4 + 1)
    for i in 1:lastindex(data)
        # copy data flipping high and low nybble
        x = unsafe_load(src, i)
        data[i] = (x & 0x0f0f0f0f0f0f0f0f) << 4 | (x & 0xf0f0f0f0f0f0f0f0) >> 4
    end
    return BioSequences.DNASequence(data, 1:seqlen, false)
end

"""
    quality(record::Record)::Vector{UInt8}

Get the base quality of  `record`.
"""
function quality(record::Record)::Vector{UInt8}
    seqlen = seqlength(record)
    offset = seqname_length(record) + n_cigar_op(record, false) * 4 + cld(seqlen, 2)
    return [reinterpret(UInt8, record.data[i+offset]) for i in 1:seqlen]
end

"""
    auxdata(record::Record)::BAM.AuxData

Get the auxiliary data of `record`.
"""
function auxdata(record::Record)
    return AuxData(record.data[auxdata_position(record):data_size(record)])
end

function Base.getindex(record::Record, tag::AbstractString)
    checkauxtag(tag)
    start = auxdata_position(record)
    stop = data_size(record)
    return getauxvalue(record.data, start, stop, UInt8(tag[1]), UInt8(tag[2]))
end

function Base.haskey(record::Record, tag::AbstractString)
    checkauxtag(tag)
    start = auxdata_position(record)
    stop = data_size(record)
    return findauxtag(record.data, start, stop, UInt8(tag[1]), UInt8(tag[2])) > 0
end

function Base.keys(record::Record)
    return collect(keys(auxdata(record)))
end

function Base.values(record::Record)
    return [record[key] for key in keys(record)]
end
