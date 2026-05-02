################################################################################
# File    : GraphAlgo.jl
# Author  : Sandeep Koranne (C) All rights reserved.
# Purpose : Example use of Julia packages and functionality, which is interesting
#
################################################################################
using Krylov, KrylovPreconditioners, Preconditioners, IterativeSolvers, LinearAlgebra, SparseArrays, ILUZero, Graphs, UnicodePlots, LinearSolve, MatrixMarket, MUMPS, SimpleWeightedGraphs, InteractiveUtils, Metis

function TrickFunction(a,b,c)
    s=a+b+c
    a/(s-a)+b/(s-b)+c/(s-c)==4
end

function CheckValues()
    [(a,b,c) for a in -1:10 for b in -1:10 for c in -10:10 if f(a,b,c)]
end
#	a=4373612677928697257861252602371390152816537558161613618621437993378423467772036
#	b=36875131794129999827197811565225474825492979968971970996283137471637224634055579
#	c=154476802108746166441951315019919837485664325669565431700026634898253202035277999
# Use BigInt
const MyFloat = Float64

function FourthOrderDifferentialOperator(N)
    return spdiagm( 0  => 6*ones(MyFloat,N),
	            -1 => -4*ones(MyFloat,N-1),
	            1  => -4*ones(MyFloat,N-1),
	            -2 => ones(MyFloat,N-2),
	            2  => ones(MyFloat,N-2) )
end

function ResOp(N)
    return spdiagm( 0  => 6*ones(MyFloat,N),
	            -1 => 4*ones(MyFloat,N-1),
	            1  => 4*ones(MyFloat,N-1),
	            -2 => ones(MyFloat,N-2),
	            2  => ones(MyFloat,N-2) )
end


function Laplacian(size;dimension=2,order=2)
    In = spdiagm( 0 => ones(MyFloat,size) )
    D  = FourthOrderDifferentialOperator(size)
    return ( kron(D,In,In) + kron(In,D,In) + kron(In,In,D) )
end

function Res(N)
    In = spdiagm( 0 => ones(MyFloat,N) )
    D  = ResOp(N)
    return ( kron(D,In,In) + kron(In,D,In) + kron(In,In,D) )
end


function CheckSP(n)
    A = Res(n)
    display(spy(A, title="Original graph"))
    g = SimpleWeightedGraph(A)
    edges_at_1 = [Edge(1, n) for n in neighbors(g, 1)]
    for ei in edges_at_1
	u,v = src(ei), dst(ei)
	w = get_weight(g,u,v)
	println("Edge $u -> $v has weight: $w")
    end
    Base.summarysize(A)
    Base.summarysize(g)
    varinfo()
    @elapsed state = dijkstra_shortest_paths(g, 1)
    target = nv(g)
    dist = state.dists[target]
    
    println("Distance from 1 to $target: ", dist)
    # Count how many nodes have an infinite distance
    unreachable_count = count(isinf, state.dists)
    
    println("Number of unreachable nodes: ", unreachable_count)    
    #@elapsed state = bellman_ford_shortest_paths(g, 1)
    p = lineplot(1:length(state.dists), state.dists, 
                 title = "Shortest Distances from Source", 
                 xlabel = "Node ID", 
                 ylabel = "Distance")

    display(p)    
end

# How to partition a graph using Metis
function Partition(g)
    display(spy(g, title="Original graph"))    
    perm, iperm = Metis.permutation(g)
    reordered_g = g[perm, perm]
    display(spy(reordered_g, title="Reordered graph"))
end

