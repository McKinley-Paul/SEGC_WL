using Test
import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using segc_wl   # or the module name inside segc_wl.jl

@testset "Basic loading" begin
    @test true
end
