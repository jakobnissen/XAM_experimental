# BAM File Format
# ===============

module BAM

import BGZFStreams
import BioAlignments
import XAM: XAM, SAM
import GenomicFeatures: GenomicFeatures, Interval
import BioSequences
import BioCore: BioCore, isfilled

include("../common.jl")
include("bai.jl")
include("reader.jl")
include("record.jl")
include("auxdata.jl")
include("writer.jl")
include("overlap.jl")

end
