@testset "BAM" begin
    bamdir = joinpath(get_bio_fmt_specimens("master", false), "BAM")

    @testset "AuxData" begin
        auxdata = BAM.AuxData(UInt8[])
        @test isempty(auxdata)

        buf = IOBuffer()
        write(buf, "NM", UInt8('s'), Int16(1))
        auxdata = BAM.AuxData(take!(buf))
        @test length(auxdata) == 1
        @test auxdata["NM"] === Int16(1)
        @test collect(auxdata) == ["NM" => Int16(1)]

        buf = IOBuffer()
        write(buf, "AS", UInt8('c'), Int8(-18))
        write(buf, "NM", UInt8('s'), Int16(1))
        write(buf, "XA", UInt8('f'), Float32(3.14))
        write(buf, "XB", UInt8('Z'), "some text\0")
        write(buf, "XC", UInt8('B'), UInt8('i'), Int32(3), Int32[10, -5, 8])
        auxdata = BAM.AuxData(take!(buf))
        @test length(auxdata) == 5
        @test auxdata["AS"] === Int8(-18)
        @test auxdata["NM"] === Int16(1)
        @test auxdata["XA"] === Float32(3.14)
        @test auxdata["XB"] == "some text"
        @test auxdata["XC"] == Int32[10, -5, 8]
        @test convert(Dict{String,Any}, auxdata) == Dict(
            "AS" => Int8(-18),
            "NM" => Int16(1),
            "XA" => Float32(3.14),
            "XB" => "some text",
            "XC" => Int32[10, -5, 8])
    end

    @testset "Record" begin
        record = BAM.Record()
        @test !isfilled(record)
        @test repr(record) == "BioAlignments.BAM.Record: <not filled>"
        @test_throws ArgumentError BAM.flag(record)
    end

    @testset "Reader" begin
        reader = open(BAM.Reader, joinpath(bamdir, "ce#1.bam"))
        @test isa(reader, BAM.Reader)
        @test eltype(reader) === BAM.Record
        @test startswith(repr(reader), "BioAlignments.BAM.Reader{IOStream}:")

        # header
        h = header(reader)
        @test isa(h, SAM.Header)

        # first record
        record = BAM.Record()
        read!(reader, record)
        @test BAM.ismapped(record)
        @test BAM.isprimary(record)
        @test ! BAM.ispositivestrand(record)
        @test BAM.refname(record) == "CHROMOSOME_I"
        @test BAM.refid(record) === 1
        @test BAM.hasnextrefid(record)
        @test BAM.nextrefid(record) === 0
        @test BAM.hasposition(record) === hasleftposition(record) === true
        @test BAM.position(record) === leftposition(record) === 2
        @test BAM.hasnextposition(record)
        @test BAM.nextposition(record) === 0
        @test rightposition(record) == 102
        @test BAM.hastempname(record) === hasseqname(record) === true
        @test BAM.tempname(record) == seqname(record) == "SRR065390.14978392"
        @test BAM.hassequence(record) === hassequence(record) === true
        @test BAM.sequence(record) == sequence(record) == dna"""
        CCTAGCCCTAACCCTAACCCTAACCCTAGCCTAAGCCTAAGCCTAAGCCT
        AAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAA
        """
        @test BAM.seqlength(record) === 100
        @test BAM.hasquality(record)
        @test eltype(BAM.quality(record)) == UInt8
        @test BAM.quality(record) == [Int(x) - 33 for x in "#############################@B?8B?BA@@DDBCDDCBC@CDCDCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"]
        @test BAM.flag(record) === UInt16(16)
        @test BAM.cigar(record) == "27M1D73M"
        @test BAM.alignment(record) == Alignment([
            AlignmentAnchor(  0,   1, OP_START),
            AlignmentAnchor( 27,  28, OP_MATCH),
            AlignmentAnchor( 27,  29, OP_DELETE),
            AlignmentAnchor(100, 102, OP_MATCH)])
        @test record["XG"] == 1
        @test record["XM"] == 5
        @test record["XN"] == 0
        @test record["XO"] == 1
        @test record["AS"] == -18
        @test record["XS"] == -18
        @test record["YT"] == "UU"
        @test keys(record) == ["XG","XM","XN","XO","AS","XS","YT"]
        @test values(record) == [1, 5, 0, 1, -18, -18, "UU"]
        @test eof(reader)
        close(reader)

        # Test conversion from byte array to record
        dsize = BAM.data_size(record)
        array = Vector{UInt8}(undef, BAM.FIXED_FIELDS_BYTES + dsize)
        GC.@preserve array record begin
            ptr = Ptr{UInt8}(pointer_from_objref(record))
            unsafe_copyto!(pointer(array), ptr, BAM.FIXED_FIELDS_BYTES)
            unsafe_copyto!(array, BAM.FIXED_FIELDS_BYTES + 1, record.data, 1, dsize)
        end
        new_record = convert(BAM.Record, array)
        @test record.bin_mq_nl == new_record.bin_mq_nl
        @test record.block_size == new_record.block_size
        @test record.flag_nc == new_record.flag_nc
        @test record.l_seq == new_record.l_seq
        @test record.next_refid == new_record.next_refid
        @test record.next_pos == new_record.next_pos
        @test record.refid == new_record.refid
        @test record.pos == new_record.pos
        @test record.tlen == new_record.tlen
        @test record.data == new_record.data

        # iterator
        @test length(collect(open(BAM.Reader, joinpath(bamdir, "ce#1.bam")))) == 1
        @test length(collect(open(BAM.Reader, joinpath(bamdir, "ce#2.bam")))) == 2

        # IOStream
        @test length(collect(BAM.Reader(open(joinpath(bamdir, "ce#1.bam"))))) == 1
        @test length(collect(BAM.Reader(open(joinpath(bamdir, "ce#2.bam"))))) == 2
    end

    @testset "Read long CIGARs" begin
        function get_cigar_lens(rec::BAM.Record)
            cigar_ops, cigar_n = BAM.cigar_rle(rec)
            field_ops, field_n = BAM.cigar_rle(rec, false)
            cigar_l = length(cigar_ops)
            field_l = length(field_ops)
            return cigar_l, field_l
        end

        function check_cigar_vs_field(rec::BAM.Record)
            cigar = BAM.cigar(rec)
            field = BAM.cigar(rec, false)
            cigar_l, field_l = get_cigar_lens(rec)
            return cigar != field && cigar_l != field_l
        end

        function check_cigar_lens(rec::BAM.Record, field_len, cigar_len)
            cigar_l, field_l = get_cigar_lens(rec)
            return cigar_l == cigar_len && field_l == field_len
        end

        reader = open(BAM.Reader, joinpath(bamdir, "cigar-64k.bam"))
        rec = BAM.Record()
        read!(reader, rec)
        @test !check_cigar_vs_field(rec)
        read!(reader, rec)
        @test check_cigar_vs_field(rec)
        @test check_cigar_lens(rec, 2, 72091)
    end

    function compare_records(xs, ys)
        if length(xs) != length(ys)
            return false
        end
        for (x, y) in zip(xs, ys)
            if !(
                x.block_size == y.block_size &&
                x.refid      == y.refid &&
                x.pos        == y.pos &&
                x.bin_mq_nl  == y.bin_mq_nl &&
                x.flag_nc    == y.flag_nc &&
                x.l_seq      == y.l_seq &&
                x.next_refid == y.next_refid &&
                x.next_pos   == y.next_pos &&
                x.tlen       == y.tlen &&
                x.data[1:BAM.data_size(x)] == y.data[1:BAM.data_size(y)])
                return false
            end
        end
        return true
    end

    @testset "Round trip" begin
        for specimen in YAML.load_file(joinpath(bamdir, "index.yml"))
            filepath = joinpath(bamdir, specimen["filename"])
            mktemp() do path, _
                # copy
                if occursin("bai", get(specimen, "tags", ""))
                    reader = open(BAM.Reader, filepath, index=filepath * ".bai")
                else
                    reader = open(BAM.Reader, filepath)
                end
                writer = BAM.Writer(
                    BGZFStream(path, "w"),
                    BAM.header(reader, fillSQ=isempty(findall(header(reader), "SQ"))))
                records = BAM.Record[]
                for record in reader
                    push!(records, record)
                    write(writer, record)
                end
                close(reader)
                close(writer)
                @test compare_records(open(collect, BAM.Reader, path), records)
            end
        end
    end

    @testset "Random access" begin
        filepath = joinpath(bamdir, "GSE25840_GSM424320_GM06985_gencode_spliced.head.bam")
        reader = open(BAM.Reader, filepath, index=filepath * ".bai")

        @test isa(eachoverlap(reader, "chr1", 1:100), BAM.OverlapIterator)
        @test isa(eachoverlap(reader, GenomicFeatures.Interval("chr1", 1, 100)), BAM.OverlapIterator)

        # expected values are counted using samtools
        for (refname, interval, expected) in [
                ("chr1", 1_000:10000,      21),
                ("chr1", 8_000:10000,      20),
                ("chr1", 766_000:800_000, 142),
                ("chr1", 786_000:800_000, 1),
                ("chr1", 796_000:800_000, 0)]
            intsect = eachoverlap(reader, refname, interval)
            @test eltype(intsect) == BAM.Record
            @test count(_ -> true, intsect) == expected
            # check that the intersection iterator is stateless
            @test count(_ -> true, intsect) == expected
        end

        # randomized tests
        for n in 1:50
            refindex = 1
            refname = "chr1"
            range = randrange(1:1_000_000)
            seekstart(reader)
            # linear scan
            expected = filter(collect(reader)) do record
                BAM.compare_intervals(record, (refindex, range)) == 0
            end
            # indexed scan
            actual = collect(eachoverlap(reader, refname, range))
            @test compare_records(actual, expected)
        end
        close(reader)

        filepath = joinpath(bamdir, "R_12h_D06.uniq.q40.bam")
        reader = open(BAM.Reader, filepath, index=filepath * ".bai")
        @test isempty(collect(eachoverlap(reader, "chr19", 5823708:5846478)))
        close(reader)
    end
end
