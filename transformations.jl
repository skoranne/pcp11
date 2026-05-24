################################################################################
# File   : transformations.jl
# Author : Sandeep Koranne (C) 2026.
# Purpose: coordinate transformations and confusion regarding conventions
#        : even though the DEF manual states its the llx, lly which is located
#        : at PLACED coordinates, we are sure its the (0,0) location of the macro
#        : in the case of FS, the (0,0) has to be compensated such that the physical
#        : location of the MACRO remains at same place.
################################################################################
const My2DBox = NTuple{2,NTuple{2,Float64}}

baremodule Orientation 
using Base: @enum
    @enum T N E S W FS FN FE FW
end

const _STR2ORI = Dict{String,Orientation.T}(
    "N"  => Orientation.N,
    "E"  => Orientation.E,
    "S"  => Orientation.S,
    "W"  => Orientation.W,
    "FS" => Orientation.FS,
    "FN" => Orientation.FN,
    "FE" => Orientation.FE,
    "FW" => Orientation.FW,
)

"""
    parse_orientation(s::AbstractString) -> Orientation

Convert a string such as `"N"` or `"FS"` into the corresponding `Orientation` enum.
Throws a `ArgumentError` if the string is not a known orientation.
The conversion is **case‑sensitive** (exact match); you can call `uppercase(s)` or
`lowercase(s)` beforehand if you want a case‑insensitive version.
"""
function parse_orientation(s::AbstractString)::Orientation.T
    ori = get(_STR2ORI, s, nothing)
    ori === nothing && throw(ArgumentError("$(s) is not a valid orientation"))
    return ori
end

function ComputeDEF2GDS(def_placement::Tuple{Float64,Float64},
                        orientation::String,
                        cellbox::My2DBox)
    X, Y = def_placement
    ((llx, lly), (urx, ury)) = cellbox

    if orientation == "N"
        return (X, Y)
    elseif orientation == "S"
        return (X + llx + urx, Y + lly + ury)
    elseif orientation == "W"
        return (X + llx + ury, Y + lly - llx)
    elseif orientation == "E"
        return (X + llx - lly, Y + lly + urx)
    elseif orientation == "FN"
        return (X + llx + urx, Y)
    elseif orientation == "FS"
        return (X, Y + lly + ury)
    elseif orientation == "FW"
        return (X + llx - lly, Y + lly - llx)
    elseif orientation == "FE"
        return (X + llx + ury, Y + lly + urx)
    else
        throw(ArgumentError("Unknown DEF orientation: $orientation"))
    end
end

#=
DEF Orientation,GDS Reflection,GDS Angle,Math Result
N (North),OFF,0°,"(x,y)"
S (South),OFF,180°,"(−x,−y)"
W (West),OFF,90°,"(−y,x)"
E (East),OFF,270°,"(y,−x)"
FS (Flip South),ON,0°,"(x,−y)"
FN (Flip North),ON,180°,"(−x,y)"
FE (Flip East),ON,90°,"(y,x)"
FW (Flip West),ON,270°,"(−y,−x)"
=#

function ComputeDEF2GDSUnitTest()
    #given the GDS (0,0) xmin,ymin we can compute the GDS location (2555.76,2350.08)
    println("Expected: 372.60 375.36")
    ibox = (372.6,372.64)
    cbox = My2DBox(((-0.19,-0.24),(0.65,2.96)))
    
    for ori in ["N","FS"]
        #ibox = (2555.76, 2352.8)
        #cbox = ((-0.19,-0.24),(0.65,2.96))
        tbox = ComputeDEF2GDS(ibox,ori,cbox)
        println("CBOX: ", cbox, " @ ", ibox, " ORI = ", ori, " tbox ", tbox)
    end
end
ComputeDEF2GDSUnitTest()

#=
I have a question about DEF PLACED transform vs GDS SREF coordinate. I have a cell whose bounding box
is ( -0.19 │  -0.24 │   0.65 │   2.96 ) and which is placed in DEF using PLACED (372600 , 372640)
using ORIENTATION "FS". The GDS location is showing as (372600   375360) and I cannot understand
how this translation is working. Can you please explain in detail first for "FS" and then for "N" how this works ?
=#
