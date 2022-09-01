# Test Painter

using StaticArrays
using LinearAlgebra
using MultiShapeCellLists
using Random
using Test


Random.seed!(1234)

@testset "just points" begin
    N = 100000
    points = rand(SVector{1,SVector{3,Float32}},N)
    naive = MultiShapeCellLists.Naive(1, 0)
    painter = MultiShapeCellLists.Painter(1,0;
        grid_start= SA[0.0,0.0,0.0],
        grid_size= SA[10,10,10],
        voxel_length= 1/10,
        max_range= SA[[0.1],Float64[]],
    )
    MultiShapeCellLists.setElements(naive,[points],[])
    MultiShapeCellLists.setElements(painter,[points],[])
    MultiShapeCellLists.addElement(naive,1,points[1])
    MultiShapeCellLists.addElement(painter,1,points[1])
    cutoff = 0.09f0
    f(x,y,i,j,d2,output) = output+((√(d2)-cutoff)^2)
    naive_out = MultiShapeCellLists.mapSimplexElements!(
        f,
        0.0,
        naive,
        1,
        points[1],
        MultiShapeCellLists.Point,
        cutoff,
    )
    painter_out = MultiShapeCellLists.mapSimplexElements!(
        f,
        0.0,
        painter,
        1,
        points[1],
        MultiShapeCellLists.Point,
        cutoff,
    )
    @test naive_out ≈ painter_out atol=1E-6
end


@testset "just lines" begin
    for trial in 1:100
        #@show trial
        N = 1000
        lines = rand(SVector{2,SVector{3,Float32}},N)
        naive = MultiShapeCellLists.Naive(0, 1)
        painter = MultiShapeCellLists.Painter(0, 1;
            grid_start= SA[0.0,0.0,0.0],
            grid_size= SA[10,10,10],
            voxel_length= 1/10,
            max_range= SA[Float64[],Float64[0.1]],
        )
        MultiShapeCellLists.setElements(naive,[],[lines])
        MultiShapeCellLists.setElements(painter,[],[lines])
        MultiShapeCellLists.addElement(naive,1,lines[1])
        MultiShapeCellLists.addElement(painter,1,lines[1])
        push!(lines,lines[1])
        cutoff = 0.09f0
        function f!(x,y,i,j,d2,output)
            d2 = MultiShapeCellLists.distSqr(lines[1], lines[j])
            if √(d2) < cutoff-1E-6
                output[j] += 1
            end
            output
        end
        naive_out = MultiShapeCellLists.mapSimplexElements!(
            f!,
            zeros(length(lines)),
            naive,
            1,
            lines[1],
            MultiShapeCellLists.Line,
            cutoff,
        )
        painter_out = MultiShapeCellLists.mapSimplexElements!(
            f!,
            zeros(length(lines)),
            painter,
            1,
            lines[1],
            MultiShapeCellLists.Line,
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
        naive = MultiShapeCellLists.Naive(1, 1)
        painter = MultiShapeCellLists.Painter(1, 1;
            grid_start= SA[0.0,0.0,0.0],
            grid_size= SA[10,10,10],
            voxel_length= 1/10,
            max_range= SA[Float64[0.1],Float64[0.1]],
        )
        MultiShapeCellLists.setElements(naive,[points],[lines])
        MultiShapeCellLists.setElements(painter,[points],[lines])
        MultiShapeCellLists.addElement(naive,1,lines[1])
        MultiShapeCellLists.addElement(painter,1,lines[1])
        MultiShapeCellLists.addElement(naive,1,points[1])
        MultiShapeCellLists.addElement(painter,1,points[1])
        push!(lines,lines[1])
        push!(points,points[1])
        cutoff = 0.1f0
        function pointline_f!(x,y,i,j,d2,output)
            d2 = MultiShapeCellLists.distSqr(points[1], lines[j])
            if √(d2) < cutoff-1E-6
                output[j] += 1
            end
            output
        end
        function linepoint_f!(x,y,i,j,d2,output)
            d2 = MultiShapeCellLists.distSqr(lines[1], points[j])
            if √(d2) < cutoff-1E-6
                output[j] += 1
            end
            output
        end
        pointline_naive_out = MultiShapeCellLists.mapSimplexElements!(
            pointline_f!,
            zeros(length(lines)),
            naive,
            1,
            points[1],
            MultiShapeCellLists.Line,
            cutoff,
        )
        linepoint_naive_out = MultiShapeCellLists.mapSimplexElements!(
            linepoint_f!,
            zeros(length(points)),
            naive,
            1,
            lines[1],
            MultiShapeCellLists.Point,
            cutoff,
        )
        pointline_painter_out = MultiShapeCellLists.mapSimplexElements!(
            pointline_f!,
            zeros(length(lines)),
            painter,
            1,
            points[1],
            MultiShapeCellLists.Line,
            cutoff,
        )
        linepoint_painter_out = MultiShapeCellLists.mapSimplexElements!(
            linepoint_f!,
            zeros(length(points)),
            painter,
            1,
            lines[1],
            MultiShapeCellLists.Point,
            cutoff,
        )
        @test linepoint_naive_out == linepoint_painter_out
        @test pointline_naive_out == pointline_painter_out
    end
end

@testset "points and lines edge cases 1" begin
    cutoff = 1f-1
    painter = MultiShapeCellLists.Painter(1, 0;
        grid_start= SA[-5.0,-5.0,-5.0],
        grid_size= SA[10,10,10],
        voxel_length= 1.0,
        max_range= SA[Float64[0.1],Float64[]],
    )
    diag = SA_F32[1,1,1]
    diag_l = norm(diag)
    diag_hat = normalize(diag)
    point = 1//2*diag - (MultiShapeCellLists.MAX_SAMPLE_SPACING + cutoff)*diag_hat
    MultiShapeCellLists.addElement(painter,1,SA[point])
    line_a = point + 0.9f0*cutoff*diag_hat
    line_b = line_a + 2*MultiShapeCellLists.MAX_SAMPLE_SPACING*diag_hat
    f(x,y,i,j,d2,output) = output + 1
    out = MultiShapeCellLists.mapSimplexElements!(
        f,
        0,
        painter,
        1,
        SA[line_a, line_b],
        MultiShapeCellLists.Point,
        cutoff,
    )
    @test out == 1
end

@testset "points and lines edge cases 2" begin
    cutoff = 1f-1
    painter = MultiShapeCellLists.Painter(1, 0;
        grid_start= SA[-5.0,-5.0,-5.0],
        grid_size= SA[10,10,10],
        voxel_length= 1.0,
        max_range= SA[Float64[0.1],Float64[]],
    )
    diag = SA_F32[1,1,1]
    diag_l = norm(diag)
    diag_hat = normalize(diag)
    point = 1//2*diag - (MultiShapeCellLists.MAX_SAMPLE_SPACING-0.001f0 + cutoff-0.001f0)*diag_hat
    MultiShapeCellLists.addElement(painter,1,SA[point])
    line_a = 1//2*diag - (MultiShapeCellLists.MAX_SAMPLE_SPACING-0.001f0)*diag_hat
    line_b = line_a + 1.9999f0*MultiShapeCellLists.MAX_SAMPLE_SPACING*diag_hat
    f(x,y,i,j,d2,output) = output + 1
    out = MultiShapeCellLists.mapSimplexElements!(
        f,
        0,
        painter,
        1,
        SA[line_a, line_b],
        MultiShapeCellLists.Point,
        cutoff,
    )
    @test out == 1
end

@testset "points and lines edge cases 3" begin
    cutoff = 1f-1
    painter = MultiShapeCellLists.Painter(1, 0;
        grid_start= SA[-5.0,-5.0,-5.0],
        grid_size= SA[10,10,10],
        voxel_length= 1.0,
        max_range= SA[Float64[0.1],Float64[]],
    )
    diag = SA_F32[1,1,1]
    diag_l = norm(diag)
    diag_hat = normalize(diag)
    point = 1//2*diag - (MultiShapeCellLists.MAX_SAMPLE_SPACING-0.001f0 + cutoff-0.001f0)*diag_hat
    MultiShapeCellLists.addElement(painter,1,SA[point])
    line_a = 1//2*diag - (MultiShapeCellLists.MAX_SAMPLE_SPACING-0.001f0)*diag_hat
    line_b = line_a + 2.1f0*MultiShapeCellLists.MAX_SAMPLE_SPACING*diag_hat
    f(x,y,i,j,d2,output) = output + 1
    out = MultiShapeCellLists.mapSimplexElements!(
        f,
        0,
        painter,
        1,
        SA[line_a, line_b],
        MultiShapeCellLists.Point,
        cutoff,
    )
    @test out == 1
end

@testset "points and lines edge cases 4" begin
    cutoff = 1f-1
    painter = MultiShapeCellLists.Painter(1, 0;
        grid_start= SA[-5.0,-5.0,-5.0],
        grid_size= SA[10,10,10],
        voxel_length= 1.0,
        max_range= SA[Float64[0.1],Float64[]],
    )
    diag = SA_F32[1,1,1]
    diag_l = norm(diag)
    diag_hat = normalize(diag)
    point = 1//2*diag - (MultiShapeCellLists.MAX_SAMPLE_SPACING-0.001f0 + cutoff-0.001f0)*diag_hat
    MultiShapeCellLists.addElement(painter,1,SA[point])
    line_a = 1//2*diag - (MultiShapeCellLists.MAX_SAMPLE_SPACING-0.001f0)*diag_hat
    line_b = line_a - 2.1f0*MultiShapeCellLists.MAX_SAMPLE_SPACING*diag_hat
    f(x,y,i,j,d2,output) = output + 1
    out = MultiShapeCellLists.mapSimplexElements!(
        f,
        0,
        painter,
        1,
        SA[line_a, line_b],
        MultiShapeCellLists.Point,
        cutoff,
    )
    @test out == 1
end

@testset "parallel lines" begin
    cutoff = 2.5f0
    painter = MultiShapeCellLists.Painter(0, 1;
        grid_start= SA[-5.0,-5.0,-5.0],
        grid_size= SA[10,10,10],
        voxel_length= 1.0,
        max_range= SA[Float64[],Float64[2.5f0]],
    )
    x = SA_F32[1,0,0]
    y = SA_F32[0,1,0]
    z = SA_F32[0,0,1]
    MultiShapeCellLists.addElement(painter,1,SA[-100x, 100x])
    f(x,y,i,j,d2,output) = output + 1
    out = MultiShapeCellLists.mapSimplexElements!(
        f,
        0,
        painter,
        1,
        SA[-100x+y, 100x+y],
        MultiShapeCellLists.Line,
        cutoff,
    )
    @test out == 1
end

@testset "almost parallel lines that are closest outside the grid" begin
    cutoff = 2.5f0
    painter = MultiShapeCellLists.Painter(0, 1;
        grid_start= SA[-5.0,-5.0,-5.0],
        grid_size= SA[10,10,10],
        voxel_length= 1.0,
        max_range= SA[Float64[],Float64[2.5f0]],
    )
    x = SA_F32[1,0,0]
    y = SA_F32[0,1,0]
    z = SA_F32[0,0,1]
    MultiShapeCellLists.addElement(painter,1,SA[-100x, 100x])
    f(x,y,i,j,d2,output) = output + 1
    out = MultiShapeCellLists.mapSimplexElements!(
        f,
        0,
        painter,
        1,
        SA[-100x+y, 100x],
        MultiShapeCellLists.Line,
        cutoff,
    )
    @test out == 1
end

@testset "points outside the grid" begin
    cutoff = 1f-1
    painter = MultiShapeCellLists.Painter(1, 0;
        grid_start= SA[-5.0,-5.0,-5.0],
        grid_size= SA[10,10,10],
        voxel_length= 1.0,
        max_range= SA[Float64[0.1],Float64[]],
    )
    x = SA_F32[1,0,0]
    y = SA_F32[0,1,0]
    z = SA_F32[0,0,1]
    MultiShapeCellLists.addElement(painter,1,SA[-100x])
    f123(x,y,i,j,d2,output) = output + 1
    out = MultiShapeCellLists.mapSimplexElements!(
        f123,
        0,
        painter,
        1,
        SA[-100.05f0x],
        MultiShapeCellLists.Point,
        cutoff,
    )
    @test out == 1
end