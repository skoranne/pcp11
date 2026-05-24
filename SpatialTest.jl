################################################################################
# File   : SpatialTest.jl
# Author : Sandeep Koranne (C) 2026. All rights reserved.
# Purpose: test for loading RTree data
################################################################################

#using CondaPkg

module ProcessGDSII

using PythonCall
using DataFrames
using Printf
using DuckDB
using DataFrames
using DataFramesMeta
using DBInterface
using OrderedCollections

gdstk = nothing
function SetupGDSTK()
    # before that do CondaPkg; and using CondaPkg
    # CondaPkg.add_pip("gdstk")
    global gdstk = pyimport("gdstk")
end


function LoadGDSFile(fileName::String)
    # Load the library using Python's gdstk.read_gds
    lib = gdstk.read_gds(fileName)
    
    # Print basic library info
    println("Successfully loaded: ", fileName)
    println("Library name: ", pyconvert(String, lib.name))
    
    # Loop through and print top-level cells
    # Note: top_level() returns a Python list of cells
    top_cells = lib.top_level()
    println("\nTop-level cells:")
    for cell in top_cells
        cell_name = pyconvert(String, cell.name)
        println("  - ", cell_name)
    end
    
    # Return the library object for further processing
    return lib
end

function Main(fileName)
    SetupGDSTK()
    gdslib = LoadGDSFile( fileName )
    # Load the file and get the top cell
    cell = gdslib.top_level()[0]  # Grab the first top-level cell
    # Extract all polygons in the cell (including paths converted to polygons)
    polygons = cell.get_polygons(include_paths=true)
    @printf("Cell '%s' contains %d polygons:", cell.name, length(polygons) )
    for (i, poly) in enumerate(polygons)
        layer = pyconvert(Int, poly.layer)
        datatype = pyconvert(Int, poly.datatype)
        
        # Convert the points to a native Julia Matrix or Vector of coordinates
        points = pyconvert(Matrix{Float64}, poly.points)
        
        println("\nPolygon #$i [Layer: $layer, DataType: $datatype]")
        println("  Vertices (X, Y):")
        # Display the coordinates
        display(points)
    end
end
const MyDB  = NamedTuple{(:geoms,:instances,:texts,:bbox), NTuple{4,DataFrames.DataFrame}}
const MyDDB = NamedTuple{(:geoms,:instances,:texts,:bbox,:defcomponents,:parentchild),
                         NTuple{6,DataFrames.DataFrame}}

function LoadGDSIntoDF(fileName) :: MyDB
    SetupGDSTK()
    lib = gdstk.read_gds(fileName)

    df_geom = DataFrame(
        CellName=String[], Layer=Int64[], DataType=Int64[], 
        ShapeType=String[],
        X_JSON=String[],
        Y_JSON=String[]
    )

    df_inst = DataFrame(
        CellName=String[], RefCellName=String[], X=Float64[], Y=Float64[], 
        MAG=Float64[], ANGLE=Float64[], IsXReflect=Bool[], IsArray=Bool[], 
        Columns=Int64[], Rows=Int64[], DX=Float64[], DY=Float64[]
    )
    df_text = DataFrame(
        CellName=String[], Layer=Int64[], TextType=Int64[], 
        X=Float64[], Y=Float64[], Label = String[]
    )
    df_bbox = DataFrame(
        CellName=String[], XMIN=Float64[], YMIN=Float64[], XMAX=Float64[], YMAX=Float64[]
    )
    # Iterate through all cells in the GDS file
    for cell in lib.cells
        c_name = pyconvert(String, cell.name)
        geom_cell = c_name
        bbox = cell.bounding_box()
        #println("Processing cell: ", geom_cell, " BBOX: ", bbox)
        # since we are using Python's structure we have to use 0 based indexing
        push!(df_bbox, (geom_cell,
                        pyconvert(Float64,bbox[0][0]),
                        pyconvert(Float64,bbox[0][1]),
                        pyconvert(Float64,bbox[1][0]),
                        pyconvert(Float64,bbox[1][1]) ) )
        geom_layer = Int64[]
        geom_datatype = Int64[]
        geom_shape = String[]
        x_geom_json = String[]
        y_geom_json = String[]        

        inst_parent = String[]
        inst_child = String[]
        inst_x = Float64[]
        inst_y = Float64[]
        inst_mag = Float64[]
        inst_angle = Float64[]
        inst_x_reflection = Bool[]
        inst_is_array = Bool[]
        inst_cols = Int64[]
        inst_rows = Int64[]
        inst_dx = Float64[]
        inst_dy = Float64[]

        text_layer = Int64[]
        text_texttype = Int64[]
        text_x = Float64[]
        text_y = Float64[]
        text_label = String[]
        
        # depth=nothing flattens all internal sub-cell references recursively
        py_polygons = cell.get_polygons(include_paths=true, depth=0)
        
        for poly in py_polygons
            lyr = pyconvert(Int64, poly.layer)
            dt = pyconvert(Int64, poly.datatype)
            
            # Convert Python points array to a Julia matrix Nx2
            pts_matrix = pyconvert(Matrix{Float64}, poly.points)
            
            # Convert matrix rows to a clean Vector of Tuples: [(x1, y1), (x2, y2), ...]
            x_pts = pts_matrix[:, 1]
            y_pts = pts_matrix[:, 1]            
            
            # Clean up 5-point rectangles where the last point repeats the first point
            if length(x_pts) == 5 && x_pts[1] == x_pts[5] && y_pts[1] == y_pts[5]
                pop!(y_pts) # Remove the redundant 5th point
                pop!(x_pts) # Remove the redundant 5th point
            end
            
            # Classify the geometry type
            sh_type = (length(x_pts) == 4) ? "Box" : "Polygon"
            # Convert numeric vector blocks directly into raw flat string arrays
            # Example: [1.0, 2.5] -> "[1.0, 2.5]"
            x_str = "[" * join(x_pts, ", ") * "]"
            y_str = "[" * join(y_pts, ", ") * "]"
              
            # Push data to columns
            #push!(cell_names, c_name)
            push!(geom_layer, lyr)
            push!(geom_datatype, dt)
            push!(geom_shape, sh_type)
            push!(x_geom_json, x_str)
            push!(y_geom_json, y_str)            
        end
        # 2. Capture structural References (SREF and AREF instances)
        # cell.references returns a list of gdstk.Reference objects
        for ref in cell.references
            ref_cell_name = pyconvert(String, ref.cell.name)
            
            # Origin point of insertion
            origin = pyconvert(Vector{Float64}, ref.origin)
            x_val, y_val = origin[1], origin[2]
            
            # Geometric transformations
            mag_val = pyconvert(Float64, ref.magnification === nothing ? 1.0 : ref.magnification)
            # Convert radians to degrees to stay native to GDSII notation
            angle_val = pyconvert(Float64, ref.rotation === nothing ? 0.0 : ref.rotation) * (180.0 / π)
            x_reflect_val = pyconvert(Bool, ref.x_reflection)
            cols_val, rows_val = 1, 1
            dx_val, dy_val = 0.0, 0.0
            is_arr = false
            rep = ref.repetition
            if rep !== nothing && pyconvert(Bool, rep != pybuiltins.None)
                if pyconvert(Bool, pyhasattr(rep, "columns")) && pyconvert(Bool, pyhasattr(rep, "rows"))
                    
                    # Safe check: Only convert if the fields are not Python None
                    if pyconvert(Bool, rep.columns != pybuiltins.None) && pyconvert(Bool, rep.rows != pybuiltins.None)
                        cols_val = pyconvert(Int64, rep.columns)
                        rows_val = pyconvert(Int64, rep.rows)
                        is_arr = (cols_val > 1 || rows_val > 1)
                    end
                    
                    if is_arr && pyconvert(Bool, pyhasattr(rep, "spacing")) && pyconvert(Bool, rep.spacing != pybuiltins.None)
                        spacing = pyconvert(Vector{Float64}, rep.spacing)
                        dx_val = spacing[1]
                        dy_val = spacing[2]
                    end
                end
            end

            push!(inst_parent, c_name)
            push!(inst_child, ref_cell_name)
            push!(inst_x, x_val)
            push!(inst_y, y_val)
            push!(inst_mag, mag_val)
            push!(inst_angle, angle_val)
            push!(inst_x_reflection, x_reflect_val)
            push!(inst_is_array, is_arr)
            push!(inst_cols, cols_val)
            push!(inst_rows, rows_val)
            push!(inst_dx, dx_val)
            push!(inst_dy, dy_val)
        end
        # Text processing
        text_data = []
        py_labels = cell.get_labels(depth=0)
        for label in py_labels
            x,y = label.origin
            push!(text_layer, pyconvert(Int64, label.layer))
            push!(text_texttype, pyconvert(Int64, label.texttype))
            push!(text_x, pyconvert(Float64, x))
            push!(text_y, pyconvert(Float64, y))
            push!(text_label, pyconvert(String, label.text))
        end
        if !isempty(geom_cell)
            cell_text_df = DataFrame( CellName=geom_cell, Layer=text_layer, TextType = text_texttype,
                                      X=text_x, Y=text_y, Label=text_label)
            append!(df_text, cell_text_df)
        end
        
        # 3. Create cell sub-groups and stream them into the master tables
        if !isempty(geom_cell)
            cell_geom_df = DataFrame(CellName=geom_cell, Layer=geom_layer, DataType=geom_datatype,
                                     ShapeType=geom_shape,
                                     X_JSON=x_geom_json,
                                     Y_JSON=y_geom_json)                                     
            append!(df_geom, cell_geom_df)
        else
            @error "How is this possible"
        end
        
        if !isempty(inst_parent)
            cell_inst_df = DataFrame(CellName=inst_parent, RefCellName=inst_child, X=inst_x, Y=inst_y,
                                     MAG=inst_mag, ANGLE=inst_angle, IsXReflect=inst_x_reflection,
                                     IsArray=inst_is_array, Columns=inst_cols, Rows=inst_rows, DX=inst_dx, DY=inst_dy)
            append!(df_inst, cell_inst_df)
        end
        
        # Free up memory reference loops explicitly for Julia garbage collector
        cell_geom_df = nothing
        cell_inst_df = nothing
        cell_text_df = nothing
        
    end # end of cell processing
    return MyDB((df_geom, df_inst, df_text, df_bbox))
end

#=
TODO:
1) Read DEF file to build correspondence table for instances
2) Convert the cell-names to unique index number and then build a reverse topological
   graph order so that 1..n is bottom to top. We have to see how we can add either VIEWS or
   new columns in existing DF and then perform grouping
3) Parse TEXT and PATH and assert on missing constructs
4) Parse meta data
5) Put all strings in the whole library in a Dict and index them
6) Compute transform matrices for SREF/AREF and then compute bounding box for cells/refs
7) Insert this into SpatialIndexing RTree
8) KLayout Python boolean OR cleanup
9) Strict mode cleanup and index table
=#

#=
from the DEF file
    - _348_ sky130_fd_sc_hd__buf_4 + PLACED ( 324300 46240 ) FS ;
    - _349_ sky130_fd_sc_hd__nand2_8 + PLACED ( 508300 65280 ) N ;
    - _350_ sky130_fd_sc_hd__buf_2 + PLACED ( 78660 59840 ) FN ;
    - _351_ sky130_fd_sc_hd__clkinv_2 + PLACED ( 67620 29920 ) FS ;
    - _352_ sky130_fd_sc_hd__clkbuf_4 + PLACED ( 525780 100640 ) S ;
    - _353_ sky130_fd_sc_hd__inv_2 + PLACED ( 134320 48960 ) FN ;
=#

function WriteGeomDFToDB(db_path::String, geoms::MyDB)
    df_geom = geoms.geoms
    df_inst = geoms.instances
    df_text = geoms.texts
    df_bbox = geoms.bbox
    # 1. Open a persistent DuckDB database file (creates it if it doesn't exist)
    db = DuckDB.DB(db_path)
    #The following disables RLE and blows up the db
    #DBInterface.execute(db, "SET force_compression = 'zstd'")    
    try
        # 2. Register the DataFrames as virtual views inside DuckDB.
        # This allows DuckDB's engine to read the DataFrame memory directly.
        DuckDB.register_data_frame(db, df_geom, "view_geom")
        DuckDB.register_data_frame(db, df_inst, "view_inst")
        DuckDB.register_data_frame(db, df_text, "view_text")
        DuckDB.register_data_frame(db, df_bbox, "view_bbox")                
        
        # 3. Use DBInterface to execute standard SQL commands.
        # This creates real, physical tables on disk populated by your views.
        DBInterface.execute(db, "CREATE TABLE geometries AS SELECT * FROM view_geom")
        DBInterface.execute(db, "CREATE TABLE instances AS SELECT * FROM view_inst")
        DBInterface.execute(db, "CREATE TABLE textlabels AS SELECT * FROM view_text")
        DBInterface.execute(db, "CREATE TABLE bbox AS SELECT * FROM view_bbox")                
        # ─── ADD THESE TWO LINES TO CLEAN UP DISK METADATA ───
        DBInterface.execute(db, "DROP VIEW IF EXISTS view_geom")
        DBInterface.execute(db, "DROP VIEW IF EXISTS view_inst")
        DBInterface.execute(db, "DROP VIEW IF EXISTS view_text")
        DBInterface.execute(db, "DROP VIEW IF EXISTS view_bbox")        
        #DBInterface.execute(db, """
        #                        COPY geometries TO 'geometries_layout.parquet' 
        #                        (FORMAT PARQUET, COMPRESSION 'ZSTD', COMPRESSION_LEVEL 7)
        #                        """)
    
        #DBInterface.execute(db, """
        #                    COPY instances TO 'instances_layout.parquet' 
        #                    (FORMAT PARQUET, COMPRESSION 'ZSTD')
        #                    """)
        
        println("Successfully wrote layout tables to database file: ", db_path)
        
    finally
        # 4. Always close the database connection to flush changes to disk safely
        DBInterface.close!(db)
    end
end

#=
DuckDB v1.5.2 (Variegata)
Enter ".help" for usage hints.
MW D show tables;
┌────────────┐
│    name    │
│  varchar   │
├────────────┤
│ bbox       │
│ geometries │
│ instances  │
│ textlabels │
└────────────┘
MW D select * from instances where Columns > 1 LIMIT 10;
┌──────────┬─────────────┬────────┬────────┬────────┬────────┬─────────┬─────────┬───────┬────────┬────────┐
│ CellName │ RefCellName │   X    │   Y    │  MAG   │ ANGLE  │ IsArray │ Columns │ Rows  │   DX   │   DY   │
│ varchar  │   varchar   │ double │ double │ double │ double │ boolean │  int64  │ int64 │ double │ double │
└──────────┴─────────────┴────────┴────────┴────────┴────────┴─────────┴─────────┴───────┴────────┴────────┘
=#


function UpdateDBForInstances(conn,df)
    DuckDB.register_data_frame(conn, df, "temp_view")
    DBInterface.execute(conn, "CREATE OR REPLACE TABLE instances AS SELECT * FROM temp_view")
    DuckDB.unregister_data_frame(conn, "temp_view")
end

#=
julia> @time mw=ProcessGDSII.ReadGeomFromDB("SMALL_ADDER.duck.db");
julia> include("/home/skoranne/GITHUB/pcp11/TopologicalOrder.jl")
julia> (relabel,odf,name2finalnum)=TopologicalProcessing.AnalyzeParentChild!(mw.instances,:CellName, :RefCellName);
julia> ProcessGDSII.AnalyzeForInstanceNames!(mw,name2finalnum,"SMALL_ADDER")
#ETL in Julia
transform!(df, :score => (s -> s .≤ 10 ? s .+ 5 : s) => :score)
=#

struct StringIndexer
    forward::Dict{String, Int}
    inverse::Vector{String}
end

# Constructor
StringIndexer() = StringIndexer(Dict{String, Int}(), String[])

# Insert and map tracking function
function get_or_assign_id!(indexer::StringIndexer, str::String)
    return get!(indexer.forward, str) do
        push!(indexer.inverse, str)
        length(indexer.inverse) # The new index ID
    end
end

# Usage:
#idx = StringIndexer()
#get_or_assign_id!(idx, "cherry") # returns 1

# Lookups:
# String -> Integer
#id = idx.forward["cherry"] # 1

# Integer -> String
#str = idx.inverse[1]       # "cherry"

# Now that we have ability to restore the data from disk, we should see how we can
# use this for checking something

function AnalysisTextLabel(df::MyDB)
    df_text = df.texts
    cellTexts = groupby( df_text, :CellName )
    idx = StringIndexer()
    for( cellName, textLabels ) in pairs( cellTexts )
        # we can convert the textLabels into a column table for easier/performant
        textMatrix = Tables.columntable( textLabels )
        # at this time textMatrix is accessible
        @inbounds for i in 1:nrow( textMatrix )
            # The two ways of accessing this matrix rows/cols, by name and col-id
            #println("::", textMatrix[1][i], " == ", textMatrix.CellName[i] )
            get_or_assign_id!( idx, textMatrix.Label[i] )
        end
    end
    idx
end

# Next we combine the placement coordinates with the DEF data
# use grep -n "END COMPO" and grep -A 5 -m 1 COMPONE to find and select the placements
#awk 'NF > 3 {printf $2 " " $3 " " $(NF-4) " " $(NF-3) " " $(NF-1) "\n"}' MW_COMP.def.txt
# Using SQL
#create table defcomponents (InstanceName Varchar, CellName Varchar, X INTEGER, Y INTEGER, Orientation Varchar);
#insert into defcomponents select * from read_csv('MW.COMPONENTS.txt', delim=' ', header = false );
#select * from defcomponents where InstanceName LIKE '%_wrapper_320';
#=
memory D select * from defcomponents where InstanceName LIKE '%_wrapper_320';
┌──────────────────────────┬─────────────────────────┬───────┬────────┬─────────────┐
│       InstanceName       │        CellName         │   X   │   Y    │ Orientation │
│         varchar          │         varchar         │ int32 │ int32  │   varchar   │
├──────────────────────────┼─────────────────────────┼───────┼────────┼─────────────┤
│ user_project_wrapper_320 │ sky130_fd_sc_hd__conb_1 │  6900 │ 239360 │ N           │
└──────────────────────────┴─────────────────────────┴───────┴────────┴─────────────┘
memory D ATTACH 'MW.duck.db' AS diskdb;
memory D COPY FROM DATABASE memory TO diskdb;
memory D DETACH diskdb;
MW D select * from defcomponents where InstanceName LIKE 'user_project_wrapper_320';
┌──────────────────────────┬─────────────────────────┬───────┬────────┬─────────────┐
│       InstanceName       │        CellName         │   X   │   Y    │ Orientation │
│         varchar          │         varchar         │ int32 │ int32  │   varchar   │
├──────────────────────────┼─────────────────────────┼───────┼────────┼─────────────┤
│ user_project_wrapper_320 │ sky130_fd_sc_hd__conb_1 │  6900 │ 239360 │ N           │
└──────────────────────────┴─────────────────────────┴───────┴────────┴─────────────┘
MW D select * from instances where RefCellName LIKE '%hd__conb_1' and X == 6.9 and Y == 239.360;
┌───────────┬─────────────────────────┬────────┬────────┬────────┬────────┬─────────┬─────────┬───────┬────────┬────────┐
│ CellName  │       RefCellName       │   X    │   Y    │  MAG   │ ANGLE  │ IsArray │ Columns │ Rows  │   DX   │   DY   │
│  varchar  │         varchar         │ double │ double │ double │ double │ boolean │  int64  │ int64 │ double │ double │
├───────────┼─────────────────────────┼────────┼────────┼────────┼────────┼─────────┼─────────┼───────┼────────┼────────┤
│ MICROWATT │ sky130_fd_sc_hd__conb_1 │    6.9 │ 239.36 │    1.0 │    0.0 │ false   │       1 │     1 │    0.0 │    0.0 │
└───────────┴─────────────────────────┴────────┴────────┴────────┴────────┴─────────┴─────────┴───────┴────────┴────────┘
MW D 
# We have to be careful matching non 'N' orientation; we have to use the BBOX and convert
MW D select * from defcomponents where InstanceName LIKE 'TAP_61613';
┌──────────────┬────────────────────────────────┬─────────┬─────────┬─────────────┐
│ InstanceName │            CellName            │    X    │    Y    │ Orientation │
│   varchar    │            varchar             │  int32  │  int32  │   varchar   │
├──────────────┼────────────────────────────────┼─────────┼─────────┼─────────────┤
│ TAP_61613    │ sky130_fd_sc_hd__tapvpwrvgnd_1 │ 2555760 │ 2358240 │ FS          │
└──────────────┴────────────────────────────────┴─────────┴─────────┴─────────────┘
MW D select * from instances where RefCellName LIKE 'sky130_fd_sc_hd__tapvpwrvgnd_1' and X == 2555.76 and Y > 2350 and Y < 2360;
┌───────────┬────────────────────────────────┬─────────┬─────────┬────────┬────────┬─────────┬─────────┬───────┬────────┬────────┐
│ CellName  │          RefCellName           │    X    │    Y    │  MAG   │ ANGLE  │ IsArray │ Columns │ Rows  │   DX   │   DY   │
│  varchar  │            varchar             │ double  │ double  │ double │ double │ boolean │  int64  │ int64 │ double │ double │
├───────────┼────────────────────────────────┼─────────┼─────────┼────────┼────────┼─────────┼─────────┼───────┼────────┼────────┤
│ MICROWATT │ sky130_fd_sc_hd__tapvpwrvgnd_1 │ 2555.76 │ 2350.08 │    1.0 │    0.0 │ false   │       1 │     1 │    0.0 │    0.0 │
│ MICROWATT │ sky130_fd_sc_hd__tapvpwrvgnd_1 │ 2555.76 │ 2355.52 │    1.0 │    0.0 │ false   │       1 │     1 │    0.0 │    0.0 │
└───────────┴────────────────────────────────┴─────────┴─────────┴────────┴────────┴─────────┴─────────┴───────┴────────┴────────┘
SDT_16 D select * from instances limit 10;
┌──────────┬──────────────────────────┬────┬────────┬────────┬────────┬────────────┬─────────┬─────────┬───────┬────────┬────────┐
│ CellName │       RefCellName        │  X │   Y    │  MAG   │ ANGLE  │ IsXReflect │ IsArray │ Columns │ Rows  │   DX   │   DY   │
│ varchar  │         varchar          │    │ double │ double │ double │  boolean   │ boolean │  int64  │ int64 │ double │ double │
├──────────┼──────────────────────────┼────┼────────┼────────┼────────┼────────────┼─────────┼─────────┼───────┼────────┼────────┤
│ SDT_16   │ sky130_fd_sc_hd__and2b_2 │7.36│  16.32 │    1.0 │    0.0 │ true       │ false   │       1 │     1 │    0.0 │    0.0 │
├──────────┼──────────────────────────┼────┼────────┼────────┼────────┼────────────┼─────────┼─────────┼───────┼────────┼────────┤
MW D select * from defcomponents where CellName LIKE 'sky130_fd_sc_hd__tapvpwrvgnd_1' and X == 2555760 and Y > 2350000 and Y < 2360000;
┌──────────────┬────────────────────────────────┬─────────┬─────────┬─────────────┐
│ InstanceName │            CellName            │    X    │    Y    │ Orientation │
│   varchar    │            varchar             │  int32  │  int32  │   varchar   │
├──────────────┼────────────────────────────────┼─────────┼─────────┼─────────────┤
│ TAP_61388    │ sky130_fd_sc_hd__tapvpwrvgnd_1 │ 2555760 │ 2352800 │ FS          │
│ TAP_61613    │ sky130_fd_sc_hd__tapvpwrvgnd_1 │ 2555760 │ 2358240 │ FS          │
└──────────────┴────────────────────────────────┴─────────┴─────────┴─────────────┘
MW D select * from bbox where CellName == 'sky130_fd_sc_hd__tapvpwrvgnd_1';
┌────────────────────────────────┬────────┬────────┬────────┬────────┐
│            CellName            │  XMIN  │  YMIN  │  XMAX  │  YMAX  │
│            varchar             │ double │ double │ double │ double │
├────────────────────────────────┼────────┼────────┼────────┼────────┤
│ sky130_fd_sc_hd__tapvpwrvgnd_1 │  -0.19 │  -0.24 │   0.65 │   2.96 │
└────────────────────────────────┴────────┴────────┴────────┴────────┘
# you can see 2.96 - 0.24 (location of Y origin) = 2.72, and FS flip south
# Using DF if you want to print some columns
# df[ df.CellName .== "RAM32_1RW1R", [:RefCellName,:inum]]
# sub_df = gdf[(CellName = "MICROWATT",)][:, [:CellName, :inum]]
# df_subset = @select(df, :x, :y, :cell_name)
df_result = @chain df begin
    @subset(:x .== 1458200)
    select(:x, :y, :cell_name)
end
julia> @chain mw.instances begin
       @subset (:CellName .== "MICROWATT")
       select(:X,:Y,:InstanceName)
       end

julia> stats = @chain mw.instances begin
           @subset(:CellName .== "MICROWATT")
           @combine(:min_x = minimum(:X), :max_x = maximum(:X))
       end
1×2 DataFrame
 Row │ min_x    max_x   
     │ Float64  Float64 
─────┼──────────────────
julia> stats = @chain mw.instances begin
           @subset(:CellName .== "MICROWATT")
           @combine(:min_i = minimum(:InstanceName), :max_i = maximum(:InstanceName))
       end
1×2 DataFrame
 Row │ min_i         max_i                    
     │ String        String                   
─────┼────────────────────────────────────────
   1 │ FILLER_0_100  user_project_wrapper_320


=#
const My2DBox = NTuple{2,NTuple{2,Float64}}
function EmptyBox() :: My2DBox
    return ((Inf,Inf),(-Inf,-Inf))
end

"""
    ComputeGDS2DEF(gds_coords::Tuple{Float64, Float64}, Box::Tuple{Float64, Float64, Float64, Float64}, orient::String; dbu::Int=1000)

Calculates the absolute DEF coordinates (in DBU) for a cell placement across all 8 standard 
DEF orientations (N, S, E, W, FN, FS, FE, FW) from its GDS coordinate reference point and local BBOX.
"""
function ComputeGDS2DEF(gds_coords::Tuple{Float64, Float64},
                        Box::Tuple{Float64, Float64, Float64, Float64},
                        orient::String; dbu::Int=1000)
    # 1. Convert everything to DBU integers
    X_gds = round(Int64, gds_coords[1] * dbu)
    Y_gds = round(Int64, gds_coords[2] * dbu)
    
    x_min = round(Int64, Box[1] * dbu)
    y_min = round(Int64, Box[2] * dbu)
    x_max = round(Int64, Box[3] * dbu)
    y_max = round(Int64, Box[4] * dbu)
    
    # 2. Define the transformation matrices [xx, xy, yx, yy] for each orientation
    # where: x' = xx*x + xy*y,  y' = yx*x + yy*y
    transformations = Dict(
        "N"  => [ 1,  0,  0,  1], # 0 deg
        "E"  => [ 0, -1,  1,  0], # 90 deg CW
        "S"  => [-1,  0,  0, -1], # 180 deg
        "W"  => [ 0,  1, -1,  0], # 270 deg CW
        "FN" => [-1,  0,  0,  1], # Flipped X (mirrored across Y axis)
        "FE" => [ 0,  1,  1,  0], # Flipped X, then 90 deg CW
        "FS" => [ 1,  0,  0, -1], # Flipped X, then 180 deg CW
        "FW" => [ 0, -1, -1,  0]  # Flipped X, then 270 deg CW
    )
    if !haskey(transformations, orient)
        throw(ArgumentError("Unknown DEF orientation: \$orient. Use N, S, E, W, FN, FS, FE, or FW."))
    end
    m = transformations[orient]
    xx, xy, yx, yy = m[1], m[2], m[3], m[4]
    # 3. Transform all 4 corners of the local cell bounding box
    local_corners = [
        (x_min, y_min),
        (x_max, y_min),
        (x_max, y_max),
        (x_min, y_max)
    ]
    
    tx_corners_x = Int64[]
    tx_corners_y = Int64[]
    
    for (cx, cy) in local_corners
        push!(tx_corners_x, xx * cx + xy * cy)
        push!(tx_corners_y, yx * cx + yy * cy)
    end
    
    # 4. Identify the lower-left corner of the newly transformed shape
    new_local_xmin = minimum(tx_corners_x)
    new_local_ymin = minimum(tx_corners_y)
    
    # 5. Extract core physical standard cell layout height (e.g., 2720 DBU for sky130)
    cell_height = y_max - abs(y_min)
    
    # 6. Resolve final toolchain tracking offset
    # X matches directly to local origin changes, Y normalizes tracking height 
    def_x = X_gds + new_local_xmin - x_min
    def_y = Y_gds + new_local_ymin - cell_height
    
    return (def_x, def_y)
end

"""
    def_to_gds_placement(def_placement::Tuple{Real, Real}, bbox::Tuple{Real, Real, Real, Real}, orientation::AbstractString)

Calculates the raw GDSII placement coordinates (X_gds, Y_gds) given the 
DEF placement point `def_placement` = (X_def, Y_def), the cell's local 
bounding box `bbox` = (llx, lly, urx, ury), and its DEF `orientation`.
    """
function ComputeDEF2GDS(def_placement,
                        orientation::String,
                        cellbox::My2DBox) :: NTuple{2,Int}
    X, Y = def_placement
    ((llx, lly), (urx, ury)) = cellbox
    llx = Int(round(llx*1000.0))
    lly = Int(round(lly*1000.0))
    urx = Int(round(urx*1000.0))
    ury = Int(round(ury*1000.0))
    
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
function ComputeDEF2GDSUnitTest()
    #given the GDS (0,0) xmin,ymin we can compute the GDS location (2555.76,2350.08) 
    ComputeDEF2GDS((2555.76, 2352.8),"FS",((-0.19,-0.24),(0.65,2.96)))
end

using Base: searchsortedfirst, searchsortedlast

"""
    FindNearestMatch(v, target_xy, target_s1) -> Union{Int, Nothing}
Searches `v` for a tuple whose first two numbers equal `target_xy`
and whose third element equals `target_s1`.

* Returns the **index** of the first matching element.
* Returns `nothing` if no such element exists.
"""
function FindNearestMatch(v::Vector{Tuple{Int,Int, String,String}},
                          target::Tuple{Int,Int},
                          cellName::String)
    isempty(v) && return nothing    
    lo = searchsortedfirst(v, (target[1],target[2]), by = t -> (t[1],t[2]))
    hi = searchsortedlast(v, (target[1],target[2]), by = t -> (t[1],t[2]))
    @assert lo <= length(v)
    @assert hi <= length(v)        
    # Scan the range as before …    
    rng = lo:hi #numeric_range(v, target_xy)          # O(log N) step
    for i in rng                               # O(k) – usually tiny
        if v[i][3] == cellName                 # exact string match
            return (index = i, element = v[i])
        end
    end
    return nothing
end

function AnalyzeForInstanceNames!(mw::MyDDB,
                                  Name2Num::OrderedDict{String,Int64},
                                  TopCellName::String)
    df   = mw.instances
    df.InstanceName  = [String("") for c in df[!,:RefCellName]]
    dc   = mw.defcomponents
    bbox = mw.bbox
    gdf = groupby( df, :CellName)
    gdc = groupby( dc, :CellName)    
    sub_df = gdf[(CellName = TopCellName,)] #this is called the NamedTuple syntax
    # Julia also allows gdf[(TopCellName,)] #both give a SubDataFrame for that cell
    box_cache = fill( EmptyBox(), size( bbox,1 ) )
    for i in eachrow( bbox )
        cellNumber = Name2Num[i.CellName]
        box_cache[cellNumber] = ((i.XMIN,i.YMIN),(i.XMAX,i.YMAX))
    end
    def_coord_cache = Vector{Tuple{Int,Int,String,String}}(undef, size( dc, 1 ) )
    for (i,row) in enumerate( eachrow( dc ) )
        cellNumber = Name2Num[row.CellName]
        @assert cellNumber <= size( box_cache, 1 )
        txcoord = ComputeDEF2GDS((Int(round(row.X)), Int(round(row.Y))), row.Orientation, box_cache[ cellNumber ])
        def_coord_cache[i] = (txcoord[1],txcoord[2],row.CellName, row.InstanceName)
    end
    sort!(def_coord_cache, by = x -> (x[1], x[2]))

    #sub-df has CellName, RefCellName, X,Y,MAG,ANGLE,IsArray,Columns,Rows,DX,DY,inum,CellNumber,RefCellNumber    
    for (i,row) in enumerate( eachrow( sub_df ) )
        @assert row.RefCellNumber <= size( box_cache, 1 )
        (idx,ele) = FindNearestMatch( def_coord_cache, (Int(round(row.X*1000.0)),Int(round(row.Y*1000.0))), row.RefCellName )
        @assert idx <= length( def_coord_cache )
        row.InstanceName = ele[4]
    end
    # Sanity check
    let k = extrema(sub_df.InstanceName)
        @assert k[1] != ""
    end
    stats = @chain mw.instances begin
        @subset(:CellName .== TopCellName)
        @combine(:min_i = minimum(:InstanceName), :max_i = maximum(:InstanceName))
    end
end

################################################################################
# The process to load the design into DB has become more complicated than we like
# Step 0: assemble the DEF components, run grep and awk
# Step 1: use mw = LoadGDSIntoDF(fileName)
# Step 2: ProcessGDSII.WriteGeomDFToDB("SMALL_ADDER.duck.db", mw);
# Step 4: now exit Julia and go to DuckDB
# This one is incomplete DB
# Then we have to import the defcomponents.txt into this db
#create table defcomponents (InstanceName Varchar, CellName Varchar, X INTEGER, Y INTEGER, Orientation Varchar);
#insert into defcomponents select * from read_csv('MW.COMPONENTS.txt', delim=' ', header = false );
# Step 5 : mw = DFLoad.ReadGeomFromDB("SMALL_ADDER.duck.db",0); #dont load parentchild
# NOStep 6 : name2num=DFLoad.ConstructName2NumDict(mw)
# Step 6 : (relabel,odf,name2finalnum)=TopologicalProcessing.AnalyzeParentChild!(mw.instances,:CellName, :RefCellName);
# ^^ This adds two extra columns into the mw.instances df which are required
# Step 7 : ProcessGDSII.AnalyzeForInstanceNames!(mw,name2num,"SMALL_ADDER")
# Step 8 : actually now mw.instances is good to go, but we want to persist it and db is READ_ONLY
# Step 9: CSV.write("FOO.csv",mw.instances);
# Step 3 : write odf into CSV.write("MW.odf.csv",odf);
# Step 10: Now load duck-db and drop instances; .import 'FOO.csv' instances
#.import 'SDT_16.odf.csv' parentchild # this is not needed until final load
# Step 11: gDesign=DFLoad.LoadDesign("SDT_16.duck.db","SDT_16")
# Step 12: DFLoad.PlotInstances(gDesign,"SDT_16");
# Now we are sure.

#=
SDT_16 D show tables;
┌───────────────┐
│     name      │
│    varchar    │
├───────────────┤
│ bbox          │
│ defcomponents │
│ geometries    │
│ instances     │
│ parentchild   │
│ textlabels    │
└───────────────┘
SDT_16 D select InstanceName from instances limit 10;
┌─────────────────────────┐
│      InstanceName       │
│         varchar         │
├─────────────────────────┤
│ _83876_                 │
│ _83875_                 │
│ PHY_EDGE_ROW_1_Left_387 │
│ PHY_EDGE_ROW_0_Left_386 │
│ FILLER_1_3              │
│ FILLER_0_3              │
│ _83872_                 │
│ _83360_                 │
│ input743                │
│ _84132_                 │
└─────────────────────────┘
=#

end #module end
