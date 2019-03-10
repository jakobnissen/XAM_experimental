using Test
import BioSequences: @dna_str, @aa_str
using BioAlignments
import BioCore.Exceptions: MissingFieldException
import BioCore.Testing.get_bio_fmt_specimens
import BGZFStreams: BGZFStream
import GenomicFeatures
import YAML

include("sam.jl")
include("bam.jl")
