using TestItemRunner

@run_package_tests

@testitem "test min distances" begin include("test_mindistance.jl") end

@testitem "test Naive" begin include("test_naive.jl") end

@testitem "test Painter basic" begin include("test_painter-basic.jl") end

@testitem "test Painter" begin include("test_painter.jl") end

