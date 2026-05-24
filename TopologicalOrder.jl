################################################################################
# File   : TopologicalOrder.jl
# Author : Sandeep Koranne (C) 2026. All rights reserved.
# Purpose: Construct relabel with topological order
################################################################################
#=
I am writing Julia program to process a dataframe which has the following schema: ParentName ChildName which shows that Child node is instantiated within parent, an example could be
Top Frame
Top Core
Frame Pad1
Frame Pad2
Core Row1
Core Row2
Core Row3
Row1 AND
Row1 INV
Row2 OR
Row2 INV
Row3 AND
Row3 INV
Write a Julia program to read this DataFrame, convert all strings to numbers, and construct a Graph with parent child as a Directed Acyclic Graph and perform topological ordering and then write the Dictionary of name/number such that numbers 1..n are in topological order. One example can be INV=1, AND=2, OR=3, Pad1=4, Pad2=5,Row1=6,Row2=7,Row3=8, Frame=9,Top=10. Please be thorough in your analysis and do not make mistakes.
=#
# ------------------------------------------------------------
#  Julia script: parent‑child → topological numbering
# ------------------------------------------------------------
module TopologicalProcessing
export AnalyzeParentChild!
using DataFrames
using DataFramesMeta
using CSV
using Graphs
using OrderedCollections   # for deterministic Dict order (optional)
#using Bijections

function UnitTest()
    # ------------------------------------------------------------------
    # 1️⃣  INPUT ---------------------------------------------------------
    # ------------------------------------------------------------------
    # Replace the block below with CSV.read("myfile.csv", DataFrame)
    raw = """
          ParentName,ChildName
          Top,Frame
          Top,Core
          Frame,Pad1
          Frame,Pad2
          Core,Row1
          Core,Row2
          Core,Row3
          Row1,AND
          Row1,INV
          Row2,OR
          Row2,INV
          Row2,INV
          Row3,AND
          Row3,INV
          """
    df = CSV.read(IOBuffer(raw), DataFrame)   # two‑column DataFrame
end

function AnalyzeParentChild!(df, ParentName::Symbol, ChildName::Symbol)
    # ------------------------------------------------------------------
    # 2️⃣  BUILD a mapping name → temporary integer id (1…n) ------------
    # ------------------------------------------------------------------
    all_names = unique(vcat(df[!,ParentName], df[!,ChildName]))
    sort!(all_names)                         # deterministic alphabetical order
    name2tempid = Dict(name => i for (i, name) in enumerate(all_names))
    tempid2name = Dict(i => name for (name, i) in name2tempid)

    # ------------------------------------------------------------------
    # 3️⃣  CONSTRUCT the DAG  (edge direction: child → parent) ---------
    # ------------------------------------------------------------------
    n = length(all_names)
    g = DiGraph(n)

    for row in eachrow(df)
        parent_id = name2tempid[row[ParentName]]
        child_id  = name2tempid[row[ChildName]]
        add_edge!(g, child_id, parent_id)    # child → parent
    end

    # ------------------------------------------------------------------
    # 4️⃣  CHECK that the graph is a DAG --------------------------------
    # ------------------------------------------------------------------
    if is_cyclic(g)
        error("Cycle detected in the parent‑child relationships. " *
            "Cannot produce a topological order.")
    end

    # ------------------------------------------------------------------
    # 5️⃣  TOPOLGICAL SORT  (children first, parents later) -------------
    # ------------------------------------------------------------------
    topo_temp_ids = topological_sort(g)      # vector of temporary ids

    # Build the final name → number dictionary (numbers 1…n in topological order)
    name2finalnum = OrderedDict{String,Int}()
    for (new_number, temp_id) in enumerate(topo_temp_ids)
        name = tempid2name[temp_id]
        name2finalnum[name] = new_number
    end

    # ------------------------------------------------------------------
    # 6️⃣  DISPLAY the result -------------------------------------------
    # ------------------------------------------------------------------
    #println("\n=== Topological numbering (child → parent) ===")
    #for (name, num) in name2finalnum
    #    println(rpad(name, 8), " => ", num)
    #end

    # ------------------------------------------------------------------
    # 7️⃣  OPTIONAL: write the mapping to a CSV file --------------------
    # ------------------------------------------------------------------
    out_df = DataFrame(Name = collect(keys(name2finalnum)),
                       Number = collect(values(name2finalnum)))
    #CSV.write("name_to_number.csv", out_df)
    #println("\nMapping saved to `name_to_number.csv`")

    # ------------------------------------------------------------------
    # 8️⃣  SANITY‑CHECK: every edge respects the numbering ---------------
    # ------------------------------------------------------------------
    function check_ordering(df, mapping)
        for row in eachrow(df)
            parent_num = mapping[row[ParentName]]
            child_num  = mapping[row[ChildName]]
            if child_num >= parent_num
                @warn "Ordering violation: child $(row.ChildName) (num=$child_num) " *
                    "should be < parent $(row.ParentName) (num=$parent_num)"
                return false
            end
        end
        return true
    end
    df.CellNumber = [name2finalnum[p] for p in df[!,ParentName]]
    df.RefCellNumber  = [name2finalnum[c] for c in df[!,ChildName]]

    #println("\nResulting DataFrame:")
    #println(df)
    all_nodes = unique(vcat(df[!,ParentName], df[!,ChildName]))                        
    master_df = DataFrame(CellName = all_nodes)

    parent_children_num = combine(groupby(df, ParentName),
                                  :RefCellNumber => (c -> [unique(collect(c))]) => :ChildrenNum)
    println( parent_children_num )
    # Join the grouped children onto the complete master list
    final_df = leftjoin(master_df, parent_children_num, on = :CellName => ParentName)
    final_df.CellNumber = [name2finalnum[p] for p in final_df[!,:CellName]]
    # Replace the `missing` arrays (the leaves) with actual empty arrays `[]`
    final_df.ChildrenNum = coalesce.(final_df.ChildrenNum, Ref(String[]))
    #println( final_df )
    println("\nSanity check → ", check_ordering(df, name2finalnum) ? "OK" : "FAILED")
    return (out_df,final_df,name2finalnum)
end

#=
Using this AnalyzeParentChild, we can get topological ordering as well as canonical number
julia> @subset relabel :Number .> 249
16×2 DataFrame
 Row │ Name                        Number 
     │ String                      Int64  
─────┼────────────────────────────────────
   1 │ sky130_fd_sc_hd__a211oi_1      250
   2 │ sky130_fd_sc_hd__a211o_4       251
   3 │ sky130_fd_sc_hd__a211o_2       252
   4 │ sky130_fd_sc_hd__a211o_1       253
   5 │ sky130_fd_sc_hd__a2111oi_4     254
   6 │ sky130_fd_sc_hd__a2111oi_2     255
   7 │ sky130_fd_sc_hd__a2111oi_1     256
   8 │ sky130_fd_sc_hd__a2111o_4      257
   9 │ sky130_fd_sc_hd__a2111o_2      258
  10 │ sky130_fd_sc_hd__a2111o_1      259
  11 │ sky130_ef_sc_hd__decap_12      260
  12 │ multiply_add_64x64             261
  13 │ RAM512                         262
  14 │ RAM32_1RW1R                    263
  15 │ Microwatt_FP_DFFRFile          264
  16 │ MICROWATT                      265
=#
end #end of module TopologicalProcessing

