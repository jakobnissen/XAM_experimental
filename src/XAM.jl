module XAM

export
    SAM,
    BAM

import BioCore: BioCore, distance, header, isfilled, seqname, hasseqname, sequence, hassequence, leftposition, rightposition, hasleftposition, hasrightposition
import BioSequences
import BioSymbols

include("sam/sam.jl")
include("bam/bam.jl")

end # module
