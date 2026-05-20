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
using DBInterface

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

function LoadGDSIntoDF(fileName)
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
        MAG=Float64[], ANGLE=Float64[], IsArray=Bool[], 
        Columns=Int64[], Rows=Int64[], DX=Float64[], DY=Float64[]
    )    
    # Iterate through all cells in the GDS file
    for cell in lib.cells
        c_name = pyconvert(String, cell.name)
        geom_cell = c_name
        println("Processing cell: ", geom_cell)
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
        inst_is_array = Bool[]
        inst_cols = Int64[]
        inst_rows = Int64[]
        inst_dx = Float64[]
        inst_dy = Float64[]        
        # depth=nothing flattens all internal sub-cell references recursively
        py_polygons = cell.get_polygons(include_paths=true, depth=nothing)
        
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
            push!(inst_is_array, is_arr)
            push!(inst_cols, cols_val)
            push!(inst_rows, rows_val)
            push!(inst_dx, dx_val)
            push!(inst_dy, dy_val)
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
            cell_inst_df = DataFrame(CellName=inst_parent, RefCellName=inst_child, X=inst_x, Y=inst_y, MAG=inst_mag, ANGLE=inst_angle, IsArray=inst_is_array, Columns=inst_cols, Rows=inst_rows, DX=inst_dx, DY=inst_dy)
            append!(df_inst, cell_inst_df)
        end
        
        # Free up memory reference loops explicitly for Julia garbage collector
        cell_geom_df = nothing
        cell_inst_df = nothing
        
    end # end of cell processing
    return (df_geom, df_inst)
    # Assemble the final unified DataFrame
    """
    Flat reader
    df = DataFrame(
        CellName = cell_names,
        Layer = layers,
        DataType = datatypes,
        ShapeType = shape_types,
        Coordinates = coordinates
    )
    return df
    """
end

function WriteGeomDFToDB(db_path::String, df_geom::DataFrame, df_inst::DataFrame)
    # 1. Open a persistent DuckDB database file (creates it if it doesn't exist)
    db = DuckDB.DB(db_path)
    #The following disables RLE and blows up the db
    #DBInterface.execute(db, "SET force_compression = 'zstd'")    
    try
        # 2. Register the DataFrames as virtual views inside DuckDB.
        # This allows DuckDB's engine to read the DataFrame memory directly.
        DuckDB.register_data_frame(db, df_geom, "view_geom")
        DuckDB.register_data_frame(db, df_inst, "view_inst")
        
        # 3. Use DBInterface to execute standard SQL commands.
        # This creates real, physical tables on disk populated by your views.
        DBInterface.execute(db, "CREATE TABLE geometries AS SELECT * FROM view_geom")
        DBInterface.execute(db, "CREATE TABLE instances AS SELECT * FROM view_inst")
        # ─── ADD THESE TWO LINES TO CLEAN UP DISK METADATA ───
        DBInterface.execute(db, "DROP VIEW IF EXISTS view_geom")
        DBInterface.execute(db, "DROP VIEW IF EXISTS view_inst")
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

end #module end

"""
DuckDB v1.5.2 (Variegata)
Enter ".help" for usage hints.
MW D show tables;
┌────────────┐
│    name    │
│  varchar   │
├────────────┤
│ geometries │
│ instances  │
└────────────┘
    """
