using SafeTestsets

@safetestset "test min distances" begin include("test_mindistance.jl") end

@safetestset "test Painter" begin include("test_painter.jl") end

