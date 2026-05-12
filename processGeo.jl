################################################################################
# File    : processGeo.jl
# Author  : Sandeep Koranne (C) 2026 All rights reserved.
# Purpose : caplet gds2geo produces a simple text file, and we would like to
#         : convert this to GMsh format for analysis
################################################################################
#using Plots #for debug
using Printf
function split_numbers(line::AbstractString, ::Type{T}=Int) where {T}
    # remove possible spaces, split on commas, convert
    return parse.(T, split(strip(line), ','))
end

"""
The format is given below:
9 < layer count
0, 0, 119 < layer id, zmin, zmax
1, 119, 238
2, 326, 505
6 < via_count layer count continues
9, 505, 936, 2, 3 < zmin, zmax and connection point
10, 1036, 1376, 3, 4
0 < layer id
1 < Number of structure ids which follow
5 < vx count, first = last
-190, 1305
1570, 1305
1570, 2910
-190, 2910
-190, 1305

"""
##

"""
    TokenStream(io::IO)

Wraps an `IO` object and provides `peek` / `next` methods that return
the next *non‑empty* line (already stripped).  The stream remembers the
last line that was peeked so that a later `next!` will return the same line.
"""
mutable struct TokenStream
    io          :: IO
    buffer      :: Union{String,Nothing}   # holds a peeked line or `nothing`
    line_number :: Int                     # for nice error messages
end

TokenStream(io::IO) = TokenStream(io, nothing, 0)

"""
    peek!(ts::TokenStream) -> Union{String,Nothing}

Return the next non‑empty line **without** consuming it.  If the underlying
iterator is exhausted, return `nothing`.
"""
function peek!(ts::TokenStream)
    # If we already have something in the buffer, just return it
    isnothing(ts.buffer) || return ts.buffer

    # Otherwise read until we find a non‑empty line or reach EOF
    for raw in eachline(ts.io)
        ts.line_number += 1
        stripped = strip(raw)
        isempty(stripped) && continue
        ts.buffer = stripped
        return stripped
    end
    return nothing          # EOF
end

"""
    next!(ts::TokenStream) -> Union{String,Nothing}

Consume and return the next line (the same line that `peek!` would have
returned).  Returns `nothing` at EOF.
"""
function next!(ts::TokenStream)
    line = peek!(ts)
    ts.buffer = nothing     # discard the buffered line
    return line
end

# Helper that converts a comma‑separated line into a vector of numbers.
_split_numbers(line::AbstractString, ::Type{T}=Int) where {T} = parse.(T, split(line, ','))

# --------------------------------------------------------------
# 3.1  Helper to raise a nice error
# --------------------------------------------------------------
function parse_error(ts::TokenStream, msg::String)
    throw(ArgumentError("Parse error at line $(ts.line_number): $msg"))
end

# --------------------------------------------------------------
# 3.2  Parse a single integer line (e.g. the “9” that starts the
#       node block).  Returns the integer.
# --------------------------------------------------------------
function parse_int_line(ts::TokenStream)
    line = next!(ts)
    isnothing(line) && parse_error(ts, "unexpected end of file while expecting an integer")
    try
        return parse(Int, line)
    catch
        parse_error(ts, "expected an integer, got ‘$(line)’")
    end
end

function parse_integers(ts::TokenStream)
    line = next!(ts)
    nums = _split_numbers(line, Int)
    if( length(nums) == 3 )
        return (nums[1], nums[2], nums[3])
    elseif( length(nums) == 5 )
        return (nums[1], nums[2], nums[3], nums[4], nums[5])
    else
        parse_error(ts, "layer line must have 3 numbers, got $(length(nums))")
    end
end


# --------------------------------------------------------------
# 3.7  Parse a vertex line: “x, y” (floating point)
# --------------------------------------------------------------
function parse_vertex_line(ts::TokenStream)
end

# --------------------------------------------------------------
# 3.8  Parse a polygon block:
#        pid
#        extra
#        n_vertices
#        {vertex_line}ⁿᵛ
# --------------------------------------------------------------
function parse_polygon_block(ts::TokenStream)
    pid   = parse_int_line(ts)            # first line of the block
    extra = parse_int_line(ts)            # second line
    nvert = parse_int_line(ts)            # how many vertices follow
    verts = Tuple{Float64,Float64}[]
    for _ = 1:nvert
        push!(verts, parse_vertex_line(ts))
    end
    return Polygon(pid, extra, verts)
end



mutable struct LayerStack
    layers :: Vector{NTuple{3,Int}}
end

function parse_layers(ts::TokenStream)
    num_layers = parse_int_line(ts)
    println("Expecting ", num_layers, " layers.")
    layers = Vector{NTuple{3,Int}}()
    for l in 1:num_layers
        (lnum,zmin,zmax) = parse_integers(ts)
        push!( layers, (lnum,zmin,zmax))
        @show (lnum,zmin,zmax)
    end
    return layers
end
function find_layer(lnum, layers)
    @show layers
    retval_idx = findfirst( x->x[1] == lnum, layers)
    return retval_idx
    #return isnothing(retval_idx) ? nothing : layers[retval_idx][1]
end

function parse_vias(ts::TokenStream,layers)
    num_vias = parse_int_line(ts)
    println("Expecting ", num_vias, " layers.")
    vias = Vector{NTuple{5,Int}}()
    for l in 1:num_vias
        (lnum,zmin,zmax,lbot,ltop) = parse_integers(ts)
        push!( vias, (lnum,zmin,zmax,lbot,ltop))
        @assert layers[lbot][1] >= 0
        @assert layers[ltop][1] >= 0
        @show "Connecting VIA between: $(layers[lbot][1])-$(layers[ltop][1])'"
        @show (lnum,zmin,zmax,lbot,ltop)
    end
    return vias
end

################################################################################
"""
    horizontal_crossings(poly, ymid) -> Vector{Float64}

Return all x‑coordinates where the horizontal line `y = ymid` crosses the
polygon `poly`.  The returned vector is sorted increasingly.
"""
function horizontal_crossings(poly::Vector{Tuple{Float64,Float64}}, ymid::Float64)
    xs = Float64[]
    n = length(poly)

    for i in 1:n
        (x1, y1) = poly[i]
        (x2, y2) = poly[mod1(i+1, n)]   # next vertex (wrap around)

        # We only care about edges that straddle the horizontal line.
        # Horizontal edges are ignored (they either lie on the line or are
        # parallel to it – they do not produce a crossing).
        if (y1 < ymid && y2 > ymid) || (y2 < ymid && y1 > ymid)
            # Edge is vertical (since polygon is orthogonal) → x is constant.
            push!(xs, x1)   # x1 == x2 for a vertical edge
        end
    end

    sort!(xs)
    return xs
end
"""
    orthogonal_decompose(poly) -> Vector{NTuple{4,Float64}}

`poly` – vector of (x, y) vertices of a simple orthogonal polygon
          (counter‑clockwise, first vertex may be repeated or not).

Return a vector of rectangles.  Each rectangle is a 4‑tuple
`(xmin, ymin, xmax, ymax)`.
"""
function orthogonal_decompose(poly::Vector{Tuple{Float64,Float64}})
    # ------------------------------------------------------------------
    # 1) collect the distinct y‑coordinates that define slab boundaries
    # ------------------------------------------------------------------
    yset = Set{Float64}(y for (_, y) in poly)
    ylist = sort(collect(yset))          # e.g. [0.0, 1.5, 3.0, …]

    # If the polygon is closed by repeating the first vertex, the set will
    # already contain the closing y‑value; otherwise we still have all needed
    # coordinates because every edge contributes both its end‑points.

    rects = NTuple{4,Float64}[]          # result container

    # ------------------------------------------------------------------
    # 2) walk over every slab
    # ------------------------------------------------------------------
    for i in 1:(length(ylist)-1)
        y_low  = ylist[i]
        y_high = ylist[i+1]
        Δy = y_high - y_low
        # Skip degenerate zero‑height slabs (should not happen for a valid polygon)
        Δy == 0 && continue

        ymid = (y_low + y_high) / 2.0   # any value strictly inside the slab

        xs = horizontal_crossings(poly, ymid)

        # xs must appear in an even number of entries; otherwise the polygon
        # is not orthogonal or not simple.
        @assert iseven(length(xs)) "Invalid polygon – odd number of crossings at y=$ymid"

        # Pair them up → rectangles
        for j in 1:2:length(xs)
            x_left  = xs[j]
            x_right = xs[j+1]
            # Discard zero‑area intervals (should not occur)
            if x_right > x_left
                push!(rects, (x_left, y_low, x_right, y_high))
            end
        end
    end

    return rects
end

function parse_geometries(layers, vias, all_polygon_roots, ts::TokenStream)
    layer_id = parse_int_line(ts)
    num_polygons = parse_int_line(ts)
    polygon_roots = all_polygon_roots[layer_id+1]
    @assert length(polygon_roots) == 0
    if(layer_id >= length(layers) )
        println("Expecting ", num_polygons, " polygons on via layer: ", layer_id)        
    else
        println("Expecting ", num_polygons, " polygons on layer: ", layer_id)                
    end
    #layer_id can be nothing if this shape is a VIA

    boxes = Vector{NTuple{4,Float64}}()
    for p in 1:num_polygons
        polygon = Vector{Tuple{Float64,Float64}}()        
        num_vx = parse_int_line(ts)
        #@assert num_vx == 5
        minx = Inf;miny = Inf;maxx = -Inf;maxy = -Inf;
        for v in 1:num_vx
            line = next!(ts)
            nums = _split_numbers(line, Float64)
            length(nums) == 2 || parse_error(ts, "vertex line must have 2 numbers, got $(length(nums))")
            minx = min( minx, nums[1] )
            miny = min( miny, nums[2] )
            maxx = max( maxx, nums[1] )
            maxy = max( maxy, nums[2] )
            push!(polygon, (nums[1], nums[2]) )
        end
        if( num_vx > 5 )
            #@show polygon            
            pboxes = orthogonal_decompose(polygon)
            next_root = length(boxes)+1
            for p in 1:length(pboxes)
                push!(polygon_roots, next_root)
            end
            push!(boxes, pboxes...) # splat operation
            #plot_decomposition(polygon, pboxes)
        else
            #singleton
            push!(boxes, (minx,miny,maxx,maxy))
            push!(polygon_roots,length(boxes))
        end
    end
    return (layer_id, boxes)
end

function plot_decomposition(poly, rects)
    plt = plot(; aspect_ratio = :equal, legend = false, grid = false)
    # draw the original polygon
    if( length(poly) > 1 )
        xs = [p[1] for p in poly]
        ys = [p[2] for p in poly]
        plot!(xs, ys, lw = 2, linecolor = :black)
    end
    # draw each rectangle with a semi‑transparent fill
    for (i, (xmin, ymin, xmax, ymax)) in enumerate(rects)
        rect_x = [xmin, xmax, xmax, xmin, xmin]
        rect_y = [ymin, ymin, ymax, ymax, ymin]
        plot!(rect_x, rect_y, seriestype = :shape,
              fillcolor = RGBA(rand(), rand(), rand(), 0.4),
              linecolor = :gray)
    end
    display(plt)
end

function parse_file(ts::TokenStream)
    layers   = parse_layers(ts)
    vias     = parse_vias(ts,layers)
    total_layers = length(layers)+length(vias)
    all_geoms= Vector{ Vector{NTuple{4,Float64}}}(undef, total_layers)
    polygon_roots = Vector{Vector{Int}}()
    for i in 1:total_layers
        push!(polygon_roots, Vector{Int}() )
    end
    while true
        (lid, geoms)  = parse_geometries(layers, vias, polygon_roots, ts)
        all_geoms[lid+1] = geoms
        line = peek!(ts)
        isnothing(line) && break            # no more polygons
        # Otherwise a polygon block starts here.
    end
    return (layers, vias, all_geoms, polygon_roots)
end

"""
    read_custom(path::AbstractString) 

Open `path` with the `open(... ) do` idiom, create a `TokenStream`,
run the recursive‑descent parser, and return the resulting `CustomData`.
"""
function read_custom(path::AbstractString)
    open(path, "r") do io
        ts = TokenStream(io)
        return parse_file(ts)
    end
end

################################################################################
# Now that we have (layers, vias, geoms) we have to write it
# as a mesh file
################################################################################
"""
    write_gmsh_boxes(
        layer_vec      :: Vector{LayerSpec},
        boxes_by_layer :: Dict{Int, Vector{Box2D}},
        filename       :: AbstractString;
        physical_name  :: String = "boxes"
    ) -> Nothing

Write a Gmsh *.geo* file that contains a `Box` primitive for every 2‑D box,
extruded to the `z`‑range that belongs to its layer.

# Arguments
* `layer_vec` – a vector of `(layer, zmin, zmax)`.  The order in the vector
  does **not** matter; the function builds a dictionary internally.
* `boxes_by_layer` – a dictionary `layer_id => Vector{Box2D}`.  Every 2‑D box
  will be turned into a 3‑D box whose bottom is at `zmin` and top at `zmax`
  of the corresponding layer.
* `filename` – name of the *.geo* file that will be created (overwrites if it
  already exists).
* `physical_name` – optional name of the Gmsh *Physical Volume* group that
  will contain all generated boxes.

# Output
A text file in Gmsh *.geo* format.  The function returns `nothing`.
    """
const LayerSpec = Vector{Tuple{Int,Float64,Float64}};
const Box2D     = NTuple{4,Float64};
function empty_box()
    return Box2D((Inf,Inf,-Inf,-Inf))
end
function merge_box(a::Box2D, b::Box2D)
    return Box2D( (min(a[1],b[1]), min(a[2],b[2]), max(a[3],b[3]), max( a[4],b[4]) ) )
end
function merge_boxes(l::Vector{Box2D}, starting_box::Box2D)
    foldl( merge_box, l; init=starting_box)
end
function count_frequencies(array)
    # Create an empty dictionary that maps the array's element type to Integers
    freq_dict = Dict{eltype(array), Int}()
    
    for num in array
        # get(dict, key, default) returns the current count, or 0 if it doesn't exist yet
        freq_dict[num] = get(freq_dict, num, 0) + 1
    end
    
    return freq_dict
end
function write_gmsh_boxes( layers, vias, 
                           boxes_by_layer :: Vector{Vector{ Box2D }},
                           polygon_roots  :: Vector{Vector{Int}},
                           filename       :: AbstractString;
                           physical_name  :: String = "boxes"
    ) :: Nothing

    # --------------------------------------------------------------
    # 2.1  Build a lookup table `layer_id → (zmin, zmax)`
    # --------------------------------------------------------------
    layer_lookup = Dict{Int, Tuple{Float64,Float64}}()
    for (lid, zmin, zmax) in layers
        @assert zmax ≥ zmin "zmax must be ≥ zmin for layer $lid"
        layer_lookup[lid+1] = (zmin, zmax)
    end
    for (lid, zmin, zmax) in vias
        @assert zmax ≥ zmin "zmax must be ≥ zmin for layer $lid"
        layer_lookup[lid+1] = (zmin, zmax)
    end
    @show layer_lookup
    full_extent = empty_box()
    for i in 1:length(layers)
        full_extent = merge_boxes( boxes_by_layer[i], full_extent )
    end
    BLOAT = 500.0
    full_extent = (full_extent[1] - BLOAT, full_extent[2] - BLOAT, full_extent[3] + BLOAT, full_extent[4] + BLOAT)
    
    total_boxes = sum( [length(v) for v in boxes_by_layer] )
    @assert sum(map(length,polygon_roots)) == total_boxes
    @show total_boxes
    # --------------------------------------------------------------
    # 2.2  Open the file and write the header
    # --------------------------------------------------------------
    box_counter = 0
    open(filename, "w") do io
        println(io, "// --------------------------------------------------")
        println(io, "// Gmsh geometry generated by `write_gmsh_boxes`")
        println(io, "// --------------------------------------------------")
        println(io, "//   number of layers   : $(length(layers))")
        println(io, "//   total 2‑D boxes   : $(sum(length(v) for v in values(boxes_by_layer)))")
        println(io, "// --------------------------------------------------")
        println(io)   # blank line
        println(io,"SetFactory(\"OpenCASCADE\");")
        println(io,"lc = 0.5; // Characteristic length (mesh size)")
        (xmin,ymin,xmax,ymax) = full_extent
        dx = xmax - xmin
        dy = ymax - ymin
        @show full_extent
        # Now generate dielectric physical groups
        for layer_id in 1:length(layers)
            (zmin, zmax) = layer_lookup[layer_id]
            if( layer_id+1 <= length(layers) )
                zmax = layer_lookup[layer_id+1][2]
            end
            dz = zmax - zmin
            @printf(io, "Box(%d) = {% .12g, % .12g, % .12g, % .12g, % .12g, % .12g}; // PhysicalGroup=%d, layer=%d\n",
                    total_boxes+layer_id, xmin, ymin, zmin, dx, dy, dz, total_boxes+layer_id, layer_id)
        end
        println(io)
        # ----------------------------------------------------------
        # 2.3  Write every box as a Gmsh `Box` primitive
        # ----------------------------------------------------------
        for layer_id in 1:length(boxes_by_layer)
            boxlist = boxes_by_layer[layer_id]
            # retrieve the vertical limits for this layer
            (zmin, zmax) = layer_lookup[layer_id]

            dz = zmax - zmin
            @assert dz > 0 "Zero‑height layer $layer_id – cannot create a volume."

            for (xmin, ymin, xmax, ymax) in boxlist
                dx = xmax - xmin
                dy = ymax - ymin
                @assert dx > 0 && dy > 0 "Degenerate 2‑D box in layer $layer_id."

                box_counter += 1
                # Gmsh syntax:  Box(x0, y0, z0; dx, dy, dz);
                @printf(io, "Box(%d) = {% .12g, % .12g, % .12g, % .12g, % .12g, % .12g}; // id=%d, layer=%d\n",
                        box_counter, xmin, ymin, zmin, dx, dy, dz, box_counter, layer_id)
            end
        end

        println(io)   # blank line
        starting_point = 1
        for layer_id in 1:length(layers) #TODO: via layers are not subtracted yet
            """
            BooleanFragments{ Volume{34}; Delete; }{ Volume{1};  }
            BooleanFragments{ Volume{35}; Delete; }{ Volume{2,3};  }
            BooleanFragments{ Volume{36}; Delete; }{ Volume{4,5,6};  }
            """
            if( length(boxes_by_layer[layer_id]) == 0 )
                continue
            end
            
            vol_ids = join(starting_point : starting_point+length(boxes_by_layer[layer_id])-1, ",")
            starting_point += length(boxes_by_layer[layer_id])
            let
                num = total_boxes + layer_id
                println(io,"my_fragments",num,"() = BooleanFragments{ Volume{",num,"}; Delete; }{ Volume{$vol_ids}; Delete; };")
            end
        end
        #for i in 1:total_boxes
        #    @printf(io,"Physical Volume(\"Rect%d\",%d ) = {%d};\n", i, i, i)
        #end
        println(io)
        starting_point = 1
        sp_volume = total_boxes+100
        for layer_id in 1:length(layers)
            root_dict = count_frequencies( polygon_roots[layer_id] )
            for (root,fq) in pairs(root_dict)
                ending_point = starting_point+fq-1
                println(io,"Physical Volume(\"L_",layer_id,"_R_",root,"\",",sp_volume," ) = {$starting_point : $ending_point};")
                starting_point += fq
                sp_volume += 1
            end
        end
        
        for layer_id in 1:length(layers)
            if( length(boxes_by_layer[layer_id]) == 0 )
                continue
            end
            
            let
                num = total_boxes + layer_id
                #Physical Volume("Diel1", 34) = my_fragments34();
                @printf(io,"Physical Volume(\"Diel%d\", %d) = my_fragments%d();\n", layer_id, num,num)
            end
        end
        #@printf(io,"Physical Volume(\"Diel\", %d) = {%d:%d};\n", total_boxes+length(layers), total_boxes+1, total_boxes+length(layers))


        # ----------------------------------------------------------
        # 2.4  (Optional) put everything into a Physical Volume group
        # ----------------------------------------------------------
        if box_counter > 0
            # Physical Volume group – the numbers 1…box_counter are the volumes
            # that were just created.
            vol_ids = join(1:box_counter, ", ")
            println(io, "Physical Volume(\"$physical_name\") = {$vol_ids};")
        end
        println(io,"// Define the mesh size for all points")
        println(io,"Mesh.MeshSizeMin = lc;")
        #println(io,"Mesh.MeshSizeMax = lc;")
        println(io,"// Optional: Synchronize to ensure the CAD engine ")
        println(io,"// communicates with the Gmsh mesh generator")
        #println(io,"Coherence;")
        println(io)

    end

    @info "Wrote $box_counter boxes to `$filename`"
    return nothing
end
