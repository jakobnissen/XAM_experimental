# BAM Record
# ==========

# the data size of fixed-length fields of the BAM object
const FIXED_FIELDS_BYTES = 36

"""
    BAM.BAMRecord()

Create an unfilled BAM record.
"""
mutable struct BAMRecord <: XAMRecord
    # fixed-length fields (see BAM specs for the details)
    block_size::Int32
    refid::Int32
    pos::Int32
    l_read_name::UInt8
    mapq::UInt8
    bin::UInt16
    n_cigar_op::UInt16
    flag::UInt16
    l_seq::Int32
    next_refid::Int32
    next_pos::Int32
    tlen::Int32
    # variable length data
    data::Vector{UInt8}
    reader::Union{Reader, Nothing}

    function BAMRecord()
        # An empty record is a legitimate BAM record.
        block_size = FIXED_FIELDS_BYTES - sizeof(Int32) + 1
        # Only the null terminator of the query name sequence
        data = UInt8[0x00]
        return new(block_size, -1, -1, 0x01, 0xff, 0, 0, 0x004, 0, -1, -1, 0, data, nothing)
    end
end

function BAMRecord(data::Vector{UInt8})
    return convert(BAMRecord, data)
end

function Base.convert(::Type{BAMRecord}, data::Vector{UInt8})
    length(data) < FIXED_FIELDS_BYTES && throw(ArgumentError("data too short"))
    record = BAMRecord()
    dst_pointer = Ptr{UInt8}(pointer_from_objref(record))
    unsafe_copyto!(dst_pointer, pointer(data), FIXED_FIELDS_BYTES)
    dsize = data_size(record)
    length(data) < dsize + FIXED_FIELDS_BYTES && throw(ArgumentError("data too short"))
    resize!(record.data, dsize)
    unsafe_copyto!(record.data, 1, data, FIXED_FIELDS_BYTES + 1, dsize)
    return record
end

function Base.copy(record::BAMRecord)
    copy = BAMRecord()
    dst_pointer = Ptr{UInt8}(pointer_from_objref(copy))
    src_pointer = Ptr{UInt8}(pointer_from_objref(record))
    unsafe_copyto!(dst_pointer, src_pointer, FIXED_FIELDS_BYTES)
    copy.data       = record.data[1:data_size(record)]
    copy.reader     = record.reader
    return copy
end

function Base.show(io::IO, record::BAMRecord)
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
    n_aux_bytes = data_size(record) - auxdata_position(record) + 1
    if n_aux_bytes > 500
        print(io, "<", n_aux_bytes, " bytes auxiliary data>")
    else
        for field in keys(auxdata(record))
            print(io, ' ', field, '='); show(record[field])
        end
    end
end

function Base.read!(reader::Reader, record::BAMRecord)
    return _read!(reader, record)
end

const NUMBERARR = Vector{UInt8}(undef, 32)
const DNALETTERS = Tuple([UInt8(i) for i in "=ACMGRSVTWYHKDBN"])

function write_int(io::IO, x::Integer)
    uint = unsigned(abs(x))
    neg = x < 0
    len = neg + ndigits(x)
    i = len
    while i > neg
        @inbounds NUMBERARR[i] = 48+rem(uint,10)
        uint = oftype(uint,div(uint,10))
        i -= 1
    end
    if neg; @inbounds NUMBERARR[1]=0x2d; end
    unsafe_write(io, pointer(NUMBERARR), len)
end

function write_to_buffer(buffer::IO, record::BAM.BAMRecord)
    data = record.data

    # Template name
    len = BAM.seqname_length(record)
    if len == 0
        write(buffer, UInt8('*'))
    else
        unsafe_write(buffer, pointer(data), len)
    end
    write(buffer, UInt8('\t'))

    # Flag
    write_int(buffer, record.flag)
    write(buffer, UInt8('\t'))

    # Reference name
    rn = BAM.refname(record)
    if ismissing(rn)
        write(buffer, UInt8('*'))
    else
        write(buffer, rn)
    end
    write(buffer, UInt8('\t'))

    # Pos, mapq
    write_int(buffer, record.pos + 1)
    write(buffer, UInt8('\t'))
    write_int(buffer, record.mapq)
    write(buffer, UInt8('\t'))

    # Cigar
    idx, nops = BAM.cigar_position(record, true)
    for i in idx:4:idx + (nops - 1) * 4
        x = unsafe_load(Ptr{UInt32}(pointer(data, i)))
        op = BioAlignments.Operation(x & 0x0F)
        write_int(buffer, x >>> 4)
        write(buffer, UInt8(convert(Char, op)))
    end
    write(buffer, UInt8('\t'))

    # Next reference
    if record.next_refid == -1
        write(buffer, UInt8('*'))
    elseif record.next_refid == record.refid
        write(buffer, UInt8('='))
    else
        nextref = BAM.nextrefname(record)
        unsafe_write(buffer, pointer(data), len)
    end
    write(buffer, UInt8('\t'))

    # Next position, TLEN
    write_int(buffer, record.next_pos + 1)
    write(buffer, UInt8('\t'))
    write_int(buffer, record.tlen)
    write(buffer, UInt8('\t'))

    # Sequence + qual
    if record.l_seq == 0
        write(buffer, [UInt8('*'), UInt8('\t'), UInt8('*')])
    else
        seqbuf = Vector{UInt8}(undef, record.l_seq)
        start = BAM.seqname_length(record) + BAM.n_cigar_op(record, false) * 4 + 2
        @inbounds for i in 1:record.l_seq
            pos = start + ((i-1) >>> 1)
            code = (record.data[pos] >>> (4*isodd(i))) & 0x0f
            seqbuf[i] = DNALETTERS[code + 1]
        end
        write(buffer, seqbuf, UInt8('\t'))
        start = start+((record.l_seq+1) >>> 1)
        @inbounds for i in 1:record.l_seq
            seqbuf[i] = data[i+start-1] + UInt8(33)
        end
        write(buffer, seqbuf)
    end

    # Auxiliary data
    auxiter = BAM.AuxDataIterator(record)
    pos = auxiter.start
    while pos <= auxiter.stop
        write(buffer, '\t')
        pos = BAM.unsafe_write_to_buffer(buffer, data, pos)
    end
    return buffer
end

function Base.convert(::Type{Vector{UInt8}}, record::BAMRecord)
    buffer = IOBuffer(sizehint=1024)
    write_to_buffer(buffer, record)
    return take!(buffer)
end

function Base.convert(::Type{String}, record::BAMRecord)
    bytes = convert(Vector{UInt8}, record)
    return String(bytes)
end

function Base.convert(::Type{SAM.SAMRecord}, record::BAMRecord)
    bytes = convert(Vector{UInt8}, record)
    return SAM.SAMRecord(bytes)
end

# ============== Accessor functions =====================
"""
    flag(record::BAMRecord)::UInt16

Get the bitwise flag of `record`.
"""
function flag(record::BAMRecord)::UInt16
    return record.flag
end

"""
    refid(record::BAMRecord)

Get the reference sequence ID of `record`.

The ID is 1-based (i.e. the first sequence is 1). Returns `missing` if not defined.
Note that unmapped records can still have a reference.
"""
# According to BAM specs, unmapped records CAN have reference and position
function refid(record::BAMRecord)::Union{Int,Missing}
    return record.refid == -1 ? missing : Int(record.refid + 1)
end

# Return the length of the read name. NOT including the null byte terminator
function seqname_length(record::BAMRecord)::UInt8
    return record.l_read_name - 0x01
end

"""
    nextrefid(record::BAMRecord)

Get the next/mate reference sequence ID of `record`.

The ID is 1-based (i.e. the first sequence is 1) and is `missing` for an unpaired
record. Note that unmapped records can still have a reference.
"""
function nextrefid(record::BAMRecord)::Union{Int,Missing}
    ispaired = flag(record) & SAM.FLAG_PAIRED == SAM.FLAG_PAIRED
    if !ispaired || record.next_refid == -1
        return missing
    end
    return  Int(record.next_refid + 1)
end

"""
    refname(record::BAMRecord)::String

Get the reference sequence name of `record`.

See also: `BAM.refid`.
"""
function refname(record::BAMRecord)::Union{String,Missing}
    id = refid(record)
    if id === missing || record.reader === nothing
        return missing
    end
    return record.reader.refseqnames[id]
end

"""
    nextrefname(record::BAMRecord)::String

Get the reference name of the mate/next read of `record`.

See also: `BAM.nextrefid`
"""
function nextrefname(record::BAMRecord)::Union{String,Missing}
    id = nextrefid(record)
    if id === missing || record.reader === nothing
        return missing
    end
    return record.reader.refseqnames[id]
end

"""
    reflen(record::BAMRecord)::Int

Get the length of the reference sequence this record applies to. Note that
unmapped records can still have a reference.
"""
function reflen(record::BAMRecord)::Union{Int,Missing}
    id = refid(record)
    if id === missing || record.reader === nothing
        return missing
    end
    return record.reader.refseqlens[id]
end

"""
    position(record::BAMRecord)

Get the 1-based leftmost mapping position of `record`. Note that unmapped records
can still have a reference position.
"""
function position(record::BAMRecord)::Union{Int,Missing}
    return record.pos == -1 ? missing : Int(record.pos + 1)
end

"""
    nextposition(record::BAMRecord)::Int

Get the 1-based leftmost mapping position of the next/mate read of `record`.
Returns `missing` if the read is unmapped.
"""
function nextposition(record::BAMRecord)::Union{Int,Missing}
    ispaired = flag(record) & SAM.FLAG_PAIRED == SAM.FLAG_PAIRED
    if !ispaired || record.next_pos == -1
        return missing
    end
    return Int(record.next_pos + 1)
end

"""
    rightposition(record::BAMRecord)::Int

Get the 1-based rightmost mapping position of `record`.
"""
function rightposition(record::BAMRecord)::Union{Int,Missing}
    pos = position(record)
    return pos === missing ? missing : pos + alignlength(record) - 1
end

"""
    mappingquality(record::BAMRecord)::UInt8

Get the mapping quality of `record`. Returns missing if the MAPQ field is set to
0xff, indicating unknown value.
"""
function mappingquality(record::BAMRecord)::Union{UInt8,Missing}
    # Mapping qual of 0xff means unknown as per BAM specs
    return record.mapq == 0xff ? missing : record.mapq
end

"""
    tempname(record::BAMRecord)::String

Get the query template name of `record`. Returns `missing` if unavailable.
"""
function tempname(record::BAMRecord)::Union{String,Missing}
    seqlen = seqname_length(record)
    return seqlen == 0 ? missing : unsafe_string(pointer(record.data), seqlen)
end

"""
    seqlength(record::BAMRecord)::Int

Get the sequence length of `record`.
"""
function seqlength(record::BAMRecord)::Int
    return Int(record.l_seq)
end

"""
    templength(record::BAMRecord)

Get the template length of `record`. Returns `missing` when unavailable.
"""
function templength(record::BAMRecord)::Union{Int,Missing}
    return record.tlen == 0 ? missing : Int(record.tlen)
end

# ================= Boolean functions ===================
"""
    ismapped(record::BAMRecord)::Bool

Test if `record` is mapped.
"""
function ismapped(record::BAMRecord)::Bool
    return flag(record) & SAM.FLAG_UNMAP == 0
end

"""
    isprimary(record::BAMRecord)::Bool

Test if `record` is a primary line of the read.

This is equivalent to `flag(record) & 0x900 == 0`.
"""
function isprimary(record::BAMRecord)::Bool
    return flag(record) & (SAM.FLAG_SECONDARY | SAM.FLAG_SUPPLEMENTARY) == 0
end

"""
    ispositivestrand(record::BAMRecord)::Bool

Test if `record` is aligned to the positive strand.

This is equivalent to `flag(record) & 0x10 == 0`.
"""
function ispositivestrand(record::BAMRecord)::Bool
    return flag(record) & SAM.FLAG_REVERSE == 0
end

"""
    isnextmapped(record::BAMRecord)::Bool

Test if the mate/next read of `record` is mapped.
"""
function isnextmapped(record::BAMRecord)::Bool
    return flag(record) & (SAM.FLAG_MUNMAP | SAM.FLAG_PAIRED) == 1
end

# ==============   SAM-compatibility functions ==========
function hastempname(record::BAMRecord)
    return seqname_length(record) > 0
end

# TODO: Update this - do we need it? What does the SAM equivalent do?
# Review this once I've rewritten SAM
function hasseqlength(record::BAMRecord)
    return true
end

function hassequence(record::BAMRecord)
    return seqlength(record) > 0
end

function hasrefid(record::BAMRecord)
    return record.refid > -1
end

function hasrefname(record::BAMRecord)
    return record.refid > -1
end

function hasposition(record::BAMRecord)
    return record.pos > -1
end

function hasrightposition(record::BAMRecord)
    return record.pos > -1
end

function hasmappingquality(record::BAMRecord)
    return mappingquality(x) !== missing
end

function hasalignment(record::BAMRecord)
    return cigar_rle(record) !== missing
end

function hasnextrefid(record::BAMRecord)
    return record.next_refid > -1
end

function hasnextrefname(record::BAMRecord)
    return record.next_refid > -1
end

function hasnextposition(record::BAMRecord)
    return record.next_pos > -1
end

function hastemplength(record::BAMRecord)
    return record.tlen != 0
end

function hasquality(record::BAMRecord)
    quals = quality(record)
    # undefined quals are 0xff per specs
    return quals !== missing && any(i != 0xff for i in quality(record))
end

function hasauxdata(record::BAMRecord)
    return auxdata_position(record) < data_size(record)
end

# ================= BioCore functions =========================

function BioCore.isfilled(record::BAMRecord)
    return true
end

function BioCore.seqname(record::BAMRecord)
    return tempname(record)
end

function BioCore.hasseqname(record::BAMRecord)
    return hastempname(record)
end

function BioCore.sequence(record::BAMRecord)
    return sequence(record)
end

function BioCore.hassequence(record::BAMRecord)
    return hassequence(record)
end

function BioCore.leftposition(record::BAMRecord)
    return position(record)
end

function BioCore.hasleftposition(record::BAMRecord)
    return hasposition(record)
end

function BioCore.rightposition(record::BAMRecord)
    return rightposition(record)
end

function BioCore.hasrightposition(record::BAMRecord)
    return hasrightposition(record)
end

# =============== Helper functions ============================

# Return the size of the `.data` field.
function data_size(record::BAMRecord)
    return record.block_size - FIXED_FIELDS_BYTES + sizeof(record.block_size)
end

function auxdata_position(record::BAMRecord)
    seqlen = seqlength(record)
    #        template_name        null_term     n_cigar_op      per_cigar_op    seq         qual
    offset = seqname_length(record) + 1 + n_cigar_op(record, false) * 4 + cld(seqlen, 2) + seqlen
    return offset + 1
end


#######################################################



"""
    n_cigar_op(record::BAMRecord, checkCG::Bool = true)

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
function n_cigar_op(record::BAMRecord, checkCG::Bool = true)
    return cigar_position(record, checkCG)[2]
end

"""
    cigar(record::BAMRecord)::String

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
function cigar(record::BAMRecord, checkCG::Bool = true)::Union{String,Missing}
    rle = cigar_rle(record, checkCG)
    if rle === missing
        return missing
    else
        buf = IOBuffer()
        for (op, len) in zip(rle...)
            print(buf, len, convert(Char, op))
        end
        return String(take!(buf))
    end
end

"""
    cigar_rle(record::BAMRecord, checkCG::Bool = true)::Tuple{Vector{BioAlignments.Operation},Vector{Int}}

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
function cigar_rle(record::BAMRecord, checkCG::Bool = true)::Union{Tuple{Vector{BioAlignments.Operation},Vector{Int}},Missing}
    idx, nops = cigar_position(record, checkCG)
    if nops == 0
        return missing
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

function cigar_position(record::BAMRecord, checkCG::Bool = true)::Tuple{Int, Int}
    cigaridx, nops = seqname_length(record) + 2, record.n_cigar_op
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
    alignment(record::BAMRecord)::BioAlignments.Alignment

Get the alignment of `record`.
"""
function alignment(record::BAMRecord)::BioAlignments.Alignment
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
    alignlength(record::BAMRecord)::Int

Get the alignment length of `record`.
"""
function alignlength(record::BAMRecord)::Int
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
    sequence(record::BAMRecord)::BioSequences.DNASequence

Get the segment sequence of `record`. Returns `missing` if not available.
"""
function sequence(record::BAMRecord)::Union{BioSequences.DNASequence, Missing}
    seqlen = seqlength(record)
    if seqlen == 0
        return missing
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
    quality(record::BAMRecord)::Vector{UInt8}

Get the base quality of `record`, or `missing` if not available.
"""
function quality(record::BAMRecord)::Union{Vector{UInt8}, Missing}
    seqlen = seqlength(record)
    if seqlen == 0
        return missing
    else
        offset = seqname_length(record) + 1 + n_cigar_op(record, false) * 4 + cld(seqlen, 2)
        return [reinterpret(UInt8, record.data[i+offset]) for i in 1:seqlen]
    end
end
