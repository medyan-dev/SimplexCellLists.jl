# Test Painter against Naive methods.

using StaticArrays
using LinearAlgebra
using SimplexCellLists
using Random
using Test

include("common.jl")



Random.seed!(1234)

@testset "just points" begin
    N = 100000
    points = rand(SVector{1,SVector{3,Float32}},N)
    naive = SimplexCellLists.Naive(1, 0, 0)
    painter = SimplexCellLists.Painter(1,0,0;
        grid_start= SA[0.1,0.1,0.1],
        grid_size= SA[10,10,10],
        voxel_length= 1/10,
        max_range= 0.1,
    )
    SimplexCellLists.setElements!(naive,[points],[],[])
    SimplexCellLists.setElements!(painter,[points],[],[])
    SimplexCellLists.addElement!(naive,1,points[1])
    SimplexCellLists.addElement!(painter,1,points[1])
    delid = rand(1:(N+1))
    SimplexCellLists.deactivate!(painter, 1, delid, SimplexCellLists.Point)
    SimplexCellLists.deactivate!(naive, 1, delid, SimplexCellLists.Point)
    cutoff = 0.09f0
    f(x,y,i,j,d2,output) = output+((√(d2)-cutoff)^2)
    naive_out = SimplexCellLists.mapSimplexElements(
        f,
        0.0,
        naive,
        points[1],
        1,
        SimplexCellLists.Point,
        cutoff,
    )
    painter_out = SimplexCellLists.mapSimplexElements(
        f,
        0.0,
        painter,
        points[1],
        1,
        SimplexCellLists.Point,
        cutoff,
    )
    @test naive_out ≈ painter_out atol=1E-6
end

@testset "just lines" begin
    for trial in 1:100
        #@show trial
        N = 1000
        lines = rand(SVector{2,SVector{3,Float32}},N)
        naive = SimplexCellLists.Naive(0, 1, 0)
        painter = SimplexCellLists.Painter(0, 1, 0;
            grid_start= SA[0.0,0.0,0.0],
            grid_size= SA[10,10,10],
            voxel_length= 1/10,
            max_range= 0.1,
        )
        SimplexCellLists.setElements!(naive,[],[lines],[])
        SimplexCellLists.setElements!(painter,[],[lines],[])
        SimplexCellLists.addElement!(naive,1,lines[1])
        SimplexCellLists.addElement!(painter,1,lines[1])
        push!(lines,lines[1])
        delid = rand(eachindex(lines))
        SimplexCellLists.deactivate!(painter, 1, delid, SimplexCellLists.Line)
        SimplexCellLists.deactivate!(naive, 1, delid, SimplexCellLists.Line)
        cutoff = 0.09f0
        function f!(x,y,i,j,d2,output)
            d2 = SimplexCellLists.distSqr(lines[1], lines[j])
            if √(d2) < cutoff-1E-5
                output[j] += 1
            end
            output
        end
        naive_out = SimplexCellLists.mapSimplexElements(
            f!,
            zeros(length(lines)),
            naive,
            lines[1],
            1,
            SimplexCellLists.Line,
            cutoff,
        )
        painter_out = SimplexCellLists.mapSimplexElements(
            f!,
            zeros(length(lines)),
            painter,
            lines[1],
            1,
            SimplexCellLists.Line,
            cutoff,
        )
        @test naive_out == painter_out
    end
end

@testset "points and lines" begin
    for trial in 1:100
        # @show trial
        N = 1000
        lines = rand(SVector{2,SVector{3,Float32}},N)
        points = rand(SVector{1,SVector{3,Float32}},N)
        naive = SimplexCellLists.Naive(1, 1, 0)
        painter = SimplexCellLists.Painter(1, 1, 0;
            grid_start= SA[0.0,0.0,0.0],
            grid_size= SA[10,10,10],
            voxel_length= 1/10,
            max_range= 0.1,
        )
        SimplexCellLists.setElements!(naive,[points],[lines],[])
        SimplexCellLists.setElements!(painter,[points],[lines],[])
        SimplexCellLists.addElement!(naive,1,lines[1])
        SimplexCellLists.addElement!(painter,1,lines[1])
        SimplexCellLists.addElement!(naive,1,points[1])
        SimplexCellLists.addElement!(painter,1,points[1])
        push!(lines,lines[1])
        push!(points,points[1])
        cutoff = 0.1f0
        function pointline_f!(x,y,i,j,d2,output)
            d2 = SimplexCellLists.distSqr(points[1], lines[j])
            if √(d2) < cutoff-1E-5
                output[j] += 1
            end
            output
        end
        function linepoint_f!(x,y,i,j,d2,output)
            d2 = SimplexCellLists.distSqr(lines[1], points[j])
            if √(d2) < cutoff-1E-5
                output[j] += 1
            end
            output
        end
        pointline_naive_out = SimplexCellLists.mapSimplexElements(
            pointline_f!,
            zeros(length(lines)),
            naive,
            points[1],
            1,
            SimplexCellLists.Line,
            cutoff,
        )
        linepoint_naive_out = SimplexCellLists.mapSimplexElements(
            linepoint_f!,
            zeros(length(points)),
            naive,
            lines[1],
            1,
            SimplexCellLists.Point,
            cutoff,
        )
        pointline_painter_out = SimplexCellLists.mapSimplexElements(
            pointline_f!,
            zeros(length(lines)),
            painter,
            points[1],
            1,
            SimplexCellLists.Line,
            cutoff,
        )
        linepoint_painter_out = SimplexCellLists.mapSimplexElements(
            linepoint_f!,
            zeros(length(points)),
            painter,
            lines[1],
            1,
            SimplexCellLists.Point,
            cutoff,
        )
        @test linepoint_naive_out == linepoint_painter_out
        @test pointline_naive_out == pointline_painter_out
    end
end

@testset "points and lines edge cases 1" begin
    cutoff = 1f-1
    painter = SimplexCellLists.Painter(1, 0, 0;
        grid_start= SA[-5.0,-5.0,-5.0],
        grid_size= SA[10,10,10],
        voxel_length= 1.0,
        max_range= SA[Float64[0.1],Float64[],Float64[]],
    )
    diag = SA_F32[1,1,1]
    diag_l = norm(diag)
    diag_hat = normalize(diag)
    point = 1//2*diag - (SimplexCellLists.MAX_SAMPLE_SPACING + cutoff)*diag_hat
    SimplexCellLists.addElement!(painter,1,SA[point])
    line_a = point + 0.9f0*cutoff*diag_hat
    line_b = line_a + 2*SimplexCellLists.MAX_SAMPLE_SPACING*diag_hat
    f(x,y,i,j,d2,output) = output + 1
    out = SimplexCellLists.mapSimplexElements(
        f,
        0,
        painter,
        SA[line_a, line_b],
        1,
        SimplexCellLists.Point,
        cutoff,
    )
    @test out == 1
end

@testset "points and lines edge cases 2" begin
    cutoff = 1f-1
    painter = SimplexCellLists.Painter(1, 0, 0;
        grid_start= SA[-5.0,-5.0,-5.0],
        grid_size= SA[10,10,10],
        voxel_length= 1.0,
        max_range= 0.1,
    )
    diag = SA_F32[1,1,1]
    diag_l = norm(diag)
    diag_hat = normalize(diag)
    point = 1//2*diag - (SimplexCellLists.MAX_SAMPLE_SPACING-0.001f0 + cutoff-0.001f0)*diag_hat
    SimplexCellLists.addElement!(painter,1,SA[point])
    line_a = 1//2*diag - (SimplexCellLists.MAX_SAMPLE_SPACING-0.001f0)*diag_hat
    line_b = line_a + 1.9999f0*SimplexCellLists.MAX_SAMPLE_SPACING*diag_hat
    f(x,y,i,j,d2,output) = output + 1
    out = SimplexCellLists.mapSimplexElements(
        f,
        0,
        painter,
        SA[line_a, line_b],
        1,
        SimplexCellLists.Point,
        cutoff,
    )
    @test out == 1
end

@testset "points and lines edge cases 3" begin
    cutoff = 1f-1
    painter = SimplexCellLists.Painter(1, 0, 0;
        grid_start= SA[-5.0,-5.0,-5.0],
        grid_size= SA[10,10,10],
        voxel_length= 1.0,
        max_range= 0.1,
    )
    diag = SA_F32[1,1,1]
    diag_l = norm(diag)
    diag_hat = normalize(diag)
    point = 1//2*diag - (SimplexCellLists.MAX_SAMPLE_SPACING-0.001f0 + cutoff-0.001f0)*diag_hat
    SimplexCellLists.addElement!(painter,1,SA[point])
    line_a = 1//2*diag - (SimplexCellLists.MAX_SAMPLE_SPACING-0.001f0)*diag_hat
    line_b = line_a + 2.1f0*SimplexCellLists.MAX_SAMPLE_SPACING*diag_hat
    f(x,y,i,j,d2,output) = output + 1
    out = SimplexCellLists.mapSimplexElements(
        f,
        0,
        painter,
        SA[line_a, line_b],
        1,
        SimplexCellLists.Point,
        cutoff,
    )
    @test out == 1
end

@testset "points and lines edge cases 4" begin
    cutoff = 1f-1
    painter = SimplexCellLists.Painter(1, 0, 0;
        grid_start= SA[-5.0,-5.0,-5.0],
        grid_size= SA[10,10,10],
        voxel_length= 1.0,
        max_range= 0.1,
    )
    diag = SA_F32[1,1,1]
    diag_l = norm(diag)
    diag_hat = normalize(diag)
    point = 1//2*diag - (SimplexCellLists.MAX_SAMPLE_SPACING-0.001f0 + cutoff-0.001f0)*diag_hat
    SimplexCellLists.addElement!(painter,1,SA[point])
    line_a = 1//2*diag - (SimplexCellLists.MAX_SAMPLE_SPACING-0.001f0)*diag_hat
    line_b = line_a - 2.1f0*SimplexCellLists.MAX_SAMPLE_SPACING*diag_hat
    f(x,y,i,j,d2,output) = output + 1
    out = SimplexCellLists.mapSimplexElements(
        f,
        0,
        painter,
        SA[line_a, line_b],
        1,
        SimplexCellLists.Point,
        cutoff,
    )
    @test out == 1
end

@testset "parallel lines" begin
    cutoff = 2.5f0
    painter = SimplexCellLists.Painter(0, 1, 0;
        grid_start= SA[-5.0,-5.0,-5.0],
        grid_size= SA[10,10,10],
        voxel_length= 1.0,
        max_range= 2.5f0,
    )
    x = SA_F32[1,0,0]
    y = SA_F32[0,1,0]
    z = SA_F32[0,0,1]
    SimplexCellLists.addElement!(painter,1,SA[-100x, 100x])
    f(x,y,i,j,d2,output) = output + 1
    out = SimplexCellLists.mapSimplexElements(
        f,
        0,
        painter,
        SA[-100x+y, 100x+y],
        1,
        SimplexCellLists.Line,
        cutoff,
    )
    @test out == 1
end

@testset "almost parallel lines that are closest outside the grid" begin
    cutoff = 2.5f0
    painter = SimplexCellLists.Painter(0, 1, 0;
        grid_start= SA[-5.0,-5.0,-5.0],
        grid_size= SA[10,10,10],
        voxel_length= 1.0,
        max_range= 2.5f0,
    )
    x = SA_F32[1,0,0]
    y = SA_F32[0,1,0]
    z = SA_F32[0,0,1]
    SimplexCellLists.addElement!(painter,1,SA[-100x, 100x])
    f(x,y,i,j,d2,output) = output + 1
    out = SimplexCellLists.mapSimplexElements(
        f,
        0,
        painter,
        SA[-100x+y, 100x],
        1,
        SimplexCellLists.Line,
        cutoff,
    )
    @test out == 1
end

@testset "points outside the grid" begin
    cutoff = 1f-1
    painter = SimplexCellLists.Painter(1, 0, 0;
        grid_start= SA[-5.0,-5.0,-5.0],
        grid_size= SA[10,10,10],
        voxel_length= 1.0,
        max_range= 0.1,
    )
    x = SA_F32[1,0,0]
    y = SA_F32[0,1,0]
    z = SA_F32[0,0,1]
    SimplexCellLists.addElement!(painter,1,SA[-100x])
    f123(x,y,i,j,d2,output) = output + 1
    out = SimplexCellLists.mapSimplexElements(
        f123,
        0,
        painter,
        SA[-100.05f0x],
        1,
        SimplexCellLists.Point,
        cutoff,
    )
    @test out == 1
end

@testset "map pairs of simplex $M in the same group" for M in (1:3)
    simplex_type = SimplexCellLists.Simplex{M}
    for trial in 1:20
        N = 1000
        scale = Float32(exp(rand()*4.0 - 2))
        elements = scale .* rand(simplex_type, N)
        naive = SimplexCellLists.Naive(1, 1, 1)
        painter = SimplexCellLists.Painter(1, 1, 1;
            min_point = SA[0.0,0.0,0.0],
            max_point = SA[scale,scale,scale],
            voxel_length= scale*1/10,
            max_range= scale*0.1,
        )
        total_elements = [[[]],[[]],[[]]]
        total_elements[M][1] = elements
        SimplexCellLists.setElements!(naive,total_elements...)
        SimplexCellLists.setElements!(painter,total_elements...)
        SimplexCellLists.addElement!(naive,1,elements[1])
        SimplexCellLists.addElement!(painter,1,elements[1])
        push!(elements,elements[1])
        delid = rand(eachindex(elements))
        SimplexCellLists.deactivate!(painter, 1, delid, simplex_type)
        SimplexCellLists.deactivate!(naive, 1, delid, simplex_type)
        cutoff = scale*0.099f0
        f!(x...) = test_f_mapPairElements!(x..., elements,cutoff)
        naive_out = SimplexCellLists.mapPairElements(
            f!,
            zeros(length(elements),length(elements)),
            naive,
            1,
            simplex_type,
            cutoff,
        )
        painter_out = SimplexCellLists.mapPairElements(
            f!,
            zeros(length(elements),length(elements)),
            painter,
            1,
            simplex_type,
            cutoff,
        )
        if naive_out != painter_out
            println("Simplex $M")
            dif_inds = findall(naive_out .!= painter_out)
            for ind in dif_inds
                @show ind
                @show naive_out[ind], painter_out[ind]
                d2 = SimplexCellLists.distSqr(lines[ind[1]], lines[ind[2]])
                @show d2
                @show √(d2)
                @show cutoff
                println()
            end
            println()
        end
        @test naive_out == painter_out
    end
end

@testset "map pairs in the same group" begin
    for trial in 1:20
        N = 1000
        triangles = trial .* rand(SVector{3,SVector{3,Float32}},N)
        lines = trial .* rand(SVector{2,SVector{3,Float32}},N)
        points = trial .* rand(SVector{1,SVector{3,Float32}},N)
        naive = SimplexCellLists.Naive(1, 1, 1)
        painter = SimplexCellLists.Painter(1, 1, 1;
            grid_start= SA[0.0,0.0,0.0],
            grid_size= SA[10,10,10],
            voxel_length= trial*1/10,
            max_range= trial*0.1,
        )
        SimplexCellLists.setElements!(naive,[points],[lines],[triangles])
        SimplexCellLists.setElements!(painter,[points],[lines],[triangles])
        SimplexCellLists.addElement!(naive,1,triangles[1])
        SimplexCellLists.addElement!(painter,1,triangles[1])
        SimplexCellLists.addElement!(naive,1,lines[1])
        SimplexCellLists.addElement!(painter,1,lines[1])
        SimplexCellLists.addElement!(naive,1,points[1])
        SimplexCellLists.addElement!(painter,1,points[1])
        push!(triangles,triangles[1])
        push!(lines,lines[1])
        push!(points,points[1])
        delid = rand(eachindex(lines))
        SimplexCellLists.deactivate!(painter, 1, delid, SimplexCellLists.Line)
        SimplexCellLists.deactivate!(naive, 1, delid, SimplexCellLists.Line)
        cutoff = trial*0.099f0
        #@show trial
        function LLf!(x,y,i,j,d2,output)
            d2 = SimplexCellLists.distSqr(lines[i], lines[j])
            if √(d2) < cutoff-1E-3
                output[i,j] += 1
                output[j,i] += 1
            end
            output
        end
        function LPf!(x,y,i,j,d2,output)
            d2 = SimplexCellLists.distSqr(lines[i], points[j])
            if √(d2) < cutoff-1E-3
                output[i,j] += 1
                output[j,i] += 1
            end
            output
        end
        function PPf!(x,y,i,j,d2,output)
            d2 = SimplexCellLists.distSqr(points[i], points[j])
            if √(d2) < cutoff-1E-3
                output[i,j] += 1
                output[j,i] += 1
            end
            output
        end
        naive_out_LL = SimplexCellLists.mapPairElements(
            LLf!,
            zeros(length(lines),length(lines)),
            naive,
            1,
            SimplexCellLists.Line,
            cutoff,
        )
        painter_out_LL = SimplexCellLists.mapPairElements(
            LLf!,
            zeros(length(lines),length(lines)),
            painter,
            1,
            SimplexCellLists.Line,
            cutoff,
        )

        naive_out_PP = SimplexCellLists.mapPairElements(
            PPf!,
            zeros(length(points),length(points)),
            naive,
            1,
            SimplexCellLists.Point,
            cutoff,
        )
        painter_out_PP = SimplexCellLists.mapPairElements(
            PPf!,
            zeros(length(points),length(points)),
            painter,
            1,
            SimplexCellLists.Point,
            cutoff,
        )

        naive_out_LP = SimplexCellLists.mapElementsElements(
            LPf!,
            zeros(length(lines),length(points)),
            naive,
            1,
            SimplexCellLists.Line,
            1,
            SimplexCellLists.Point,
            cutoff,
        )
        painter_out_LP = SimplexCellLists.mapElementsElements(
            LPf!,
            zeros(length(lines),length(points)),
            painter,
            1,
            SimplexCellLists.Line,
            1,
            SimplexCellLists.Point,
            cutoff,
        )
        for (painter_out, naive_out, name) in ((painter_out_LL, naive_out_LL, "Line Line"), (painter_out_LP, naive_out_LP, "Line Point"), (painter_out_PP, naive_out_PP, "Point Point"))
            if naive_out != painter_out
                println(name)
                dif_inds = findall(naive_out .!= painter_out)
                for ind in dif_inds
                    @show ind
                    @show naive_out[ind], painter_out[ind]
                    d2 = SimplexCellLists.distSqr(lines[ind[1]], lines[ind[2]])
                    @show d2
                    @show √(d2)
                    @show cutoff
                    println()
                end
                println()
            end
            @test naive_out == painter_out
        end
    end
end