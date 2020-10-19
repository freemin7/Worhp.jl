using SafeTestsets

@safetestset "RawBindingTest" begin
    include("RawBindingTest.jl")
end

@safetestset "BindingTest" begin
    include("BindingTest.jl")
end
