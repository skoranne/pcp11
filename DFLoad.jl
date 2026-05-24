################################################################################
# File   : DFLoad.jl
# Author : Sandeep Koranne (C) 2026. All rights reserved.
# Purpose: Load a DuckDB into a DataFrame and do some analytics on it.
################################################################################

module DFLoad

using DataFrames
using DataFramesMeta
using OrderedCollections
using Printf
using DuckDB
using DBInterface
using SpatialIndexing
using Plots
using OrderedCollections
#using GLMakie #for PlotInstances
export LoadDesign
const MyDDB = NamedTuple{(:geoms,:instances,:texts,:bbox,:defcomponents,:parentchild),
                         NTuple{6,DataFrames.DataFrame}}
const My2DBox = NTuple{2,NTuple{2,Float64}}
struct TreeBox
    inum :: Int
end

mutable struct LoadedDesign
    mw         :: MyDDB
    name2num   :: OrderedDict{String,Int}
    bbox_cache :: Vector{My2DBox}
    rtree
end


function ReadGeomFromDB(db_path::String, phase::Int=1) :: MyDDB
    configuration = DuckDB.Config()
    DuckDB.set_config( configuration, "access_mode", "READ_ONLY" )
    conn = DuckDB.DB(db_path, configuration)
    df_geom = DataFrame( DBInterface.execute( conn, "SELECT * FROM geometries") )
    df_inst = DataFrame( DBInterface.execute( conn, "SELECT * FROM instances") )
    # we want to assign inum for each instance and then assign defcomponent name to that index
    @transform!(groupby(df_inst,:CellName), :inum = eachindex(:CellName))
    df_text = DataFrame( DBInterface.execute( conn, "SELECT * FROM textlabels") )
    df_bbox = DataFrame( DBInterface.execute( conn, "SELECT * FROM bbox") )
    # see notes on how to build the defcomponents table in the DB
    df_defcomponents = DataFrame( DBInterface.execute( conn, "SELECT * FROM defcomponents") )
    df_parentchild = DataFrame()
    if( phase > 0 )
        println("Loading parent-child relationship")
        df_parentchild = DataFrame( DBInterface.execute( conn, "SELECT * FROM parentchild") )
    end
    close( conn )
    MyDDB((df_geom, df_inst, df_text, df_bbox, df_defcomponents, df_parentchild))
end

# use the parentchild df
function ConstructName2NumDict(mw :: MyDDB)::OrderedDict{String,Int}
    df = mw.parentchild
    # --------------------------------------------------------------------
    # 1️⃣  Verify that the required columns exist
    # --------------------------------------------------------------------
    @assert ("CellName"  in names(df))  "DataFrame must contain a column named CellName"
    @assert ("CellNumber" in names(df)) "DataFrame must contain a column named CellNumber"

    # --------------------------------------------------------------------
    # 2️⃣  Pull the two columns out as vectors (this is zero‑copy)
    # --------------------------------------------------------------------
    names_col   = df[!, "CellName"]    # Vector{<:AbstractString}
    numbers_col = df[!, "CellNumber"]  # Vector{<:Integer or Real}

    # --------------------------------------------------------------------
    # 3️⃣  Allocate the dictionary with the (optional but fast)
    # --------------------------------------------------------------------
    n = length(names_col)
    d = OrderedDict{String,Int}()
    sizehint!(d,n) # pre‑allocate capacity

    # --------------------------------------------------------------------
    # 4️⃣  Fill it – skip rows that contain `missing`
    # --------------------------------------------------------------------
    for i in eachindex(names_col, numbers_col)
        name   = names_col[i]
        number = numbers_col[i]
        # If either entry is missing, ignore that row
        if name === missing || number === missing
            @warn "$name is missing or $number is missing"
            continue
        end
        # the exact types we want
        d[String(name)] = Int(number)
    end
    d
end


# Please note this bbox cache is in the CellNumber order
function ConstructBBOXCache( mw :: MyDDB, name2num :: OrderedDict{String,Int} ) ::Vector{My2DBox}    
    df = mw.bbox
    # --------------------------------------------------------------------
    # 1️⃣  Verify that the required columns exist
    # --------------------------------------------------------------------
    @assert ("CellName"  in names(df))  "DataFrame must contain a column named CellName"
    # --------------------------------------------------------------------
    # 2️⃣  Pull the two columns out as vectors (this is zero‑copy)
    # --------------------------------------------------------------------
    names_col   = df[!, "CellName"]    # Vector{<:AbstractString}
    #vec_of_tuples = NTuple{2, NTuple{2, Int}}.(zip(df.XMIN, df.YMIN, df.XMAX, df.YMAX)) do (a, b, c, d)
    #    ((a, b), (c, d))
    #end
    matrix_view = Matrix(df[:, [:XMIN, :YMIN, :XMAX, :YMAX]])
    vec_of_tuples = map(r -> ((matrix_view[r, 1], matrix_view[r, 2]), 
                              (matrix_view[r, 3], matrix_view[r, 4])), 
                        1:size(matrix_view, 1))
    p = map(str -> name2num[str], names_col)
    invpermute!( vec_of_tuples, p)
    vec_of_tuples
    #retval[p] = vec_of_tuples #inplace scatter
end

function TransformBox(sref_pt::Tuple{Float64, Float64},
                      cell_bbox::My2DBox,
                      angle::Float64,
                      reflected::Bool)
    (X, Y) = sref_pt
    ((xmin, ymin), (xmax, ymax)) = cell_bbox
    
    # 1. Map the 4 local corners of the bounding box
    corners = [
        (xmin, ymin),
        (xmax, ymin),
        (xmin, ymax),
        (xmax, ymax)
    ]
    
    # Convert angle to radians
    rad = deg2rad(angle)
    c, s = cos(rad), sin(rad)
    mirror = reflected ? -1.0 : 1.0
    
    # 2. Transform each corner into the top level space
    tx_corners = map(corners) do (x, y)
        # Apply reflection on Y first, then rotate, then shift by SREF origin
        y_mod = y * mirror
        x_top = X + (x * c - y_mod * s)
        y_top = Y + (x * s + y_mod * c)
        return (x_top, y_top)
    end
    
    # 3. Extract the true min and max bounds from the transformed points
    xs = map(p -> p[1], tx_corners)
    ys = map(p -> p[2], tx_corners)
    
    return ((minimum(xs), minimum(ys)), (maximum(xs), maximum(ys)))
end


function ConstructBBOXTree( mw::MyDDB, box_cache::Vector{My2DBox},
                            TopCellName::String)
    df   = mw.instances
    gdf = groupby( df, :CellName)
    sub_df = gdf[(CellName = TopCellName,)] #this is called the NamedTuple syntax
    rtree = SpatialIndexing.RTree{Float64, 2}(TreeBox)
    for (i,row) in enumerate( eachrow( sub_df ) )
        @assert row.RefCellNumber <= size( box_cache, 1 )
        box = box_cache[row.RefCellNumber]
        #tbox = ((box[1][1] + row.X, box[1][2] + row.Y), (box[2][1] + row.X, box[2][2] + row.Y))
        tbox = TransformBox( ( row.X, row.Y ), box, row.ANGLE, row.IsXReflect )
        insert!( rtree, SpatialIndexing.Rect{Float64,2}(tbox[1], tbox[2]), TreeBox(i) )
    end
    
    rtree
end

function LoadDesign( fileName :: String, topCellName :: String )
    mw = ReadGeomFromDB( fileName )
    name2num = ConstructName2NumDict( mw )
    box_cache = ConstructBBOXCache( mw, name2num )
    rtree = ConstructBBOXTree( mw, box_cache, topCellName )
    retval = LoadedDesign( mw, name2num, box_cache, rtree)
    retval
end

function PlotInstances(design, topCellName::String)
    instances = design.mw.instances
    box_cache = design.bbox_cache
    gdf = groupby( instances, :CellName )
    sub_df = gdf[(CellName = topCellName,)]
    rtree = design.rtree
    plt = plot(; aspect_ratio = :equal, legend = false, grid = false)    
    for (i,row) in enumerate( eachrow( sub_df ) )
        box = box_cache[row.RefCellNumber]
        b = TransformBox( ( row.X, row.Y ), box, row.ANGLE, row.IsXReflect )        
        (xmin,xmax,ymin,ymax)=(b[1][1],b[2][1],b[1][2],b[2][2])
        rect_x = [xmin, xmax, xmax, xmin, xmin]
        rect_y = [ymin, ymin, ymax, ymax, ymin]
        plot!(rect_x, rect_y, lw = 1, linecolor = :black)
        if( i > 100 )
            break
        end
        
    end
    display(plt)
end

end # of module DFLoad
