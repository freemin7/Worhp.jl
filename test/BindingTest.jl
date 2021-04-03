using WorhpOpt
using Test

function UserF(opt, wsp, par, cnt)
    x = unsafe_wrap(Array{Float64,1}, opt.X, 4; own=false)
    opt.F = wsp.ScaleObj * (x[1]*x[1] + 2.0*x[2]*x[2] - x[3])
end

function UserG(opt, wsp, par, cnt)
    x = unsafe_wrap(Array{Float64,1}, opt.X, 4; own=false)
    G = unsafe_wrap(Array{Float64,1}, opt.G, 3; own=false)

    G[1] = x[1]*x[1] + x[3]*x[3] + x[1] * x[3]
    G[2] = x[3] - x[4]
    G[3] = x[2] + x[4]
end

function UserDF(opt, wsp, par, cnt)
    x = unsafe_wrap(Array{Float64,1}, opt.X, 4; own=false)
    DF = unsafe_wrap(Array{Float64,1}, wsp.DF.val, 3; own=false)

    DF[1] = wsp.ScaleObj * 2. * x[1]
    DF[2] = wsp.ScaleObj * 4. * x[2]
    DF[3] = wsp.ScaleObj * -1.
end

function UserDG(opt, wsp, par, cnt)
    x = unsafe_wrap(Array{Float64,1}, opt.X, 4; own=false)
    DG = unsafe_wrap(Array{Float64,1}, wsp.DG.val, 6; own=false)

    DG[1] = 2.0 * x[1] + x[3]
    DG[2] = 1.
    DG[3] = x[1] + 2. * x[3]
    DG[4] =  1.
    DG[5] = -1.
    DG[6] =  1.

end

function UserHM(opt, wsp, par, cnt)
    x = unsafe_wrap(Array{Float64,1}, opt.X, 4; own=false)
    HM = unsafe_wrap(Array{Float64,1}, wsp.HM.val, 5; own=false)

    μ₀ = unsafe_load(opt.Mu,1)

    HM[1] = μ₀;
    HM[2] = wsp.ScaleObj * 2. + 2. * μ₀
    HM[3] = wsp.ScaleObj * 4.
    HM[4] = 2. * μ₀;
    HM[5] = 0.

end

println("Before WorhpProblem")

p = WorhpProblem(
    Int32(4),Int32(3),Int32(3),Int32(6),Int32(5),
    true, true, true,
    Vector{Int32}([1,2,3]),
    Vector{Int32}([1,3,1,2,2,3]), Vector{Int32}([1,2,3,3,4,4]),
    Vector{Int32}([3,1,2,3,4]), Vector{Int32}([1,1,2,3,4]),
    #Multiple dispatch can be used to make that nice, such as passing a matrix
    [2.0, 2.0, 1.0, 0.0],
    [0.0, 0.0, 0.0, 0.0],[0.0, 0.0, 0.0],
    [-0.5, -2.0, 0.0, -2.0],[Inf64, Inf64, 2.0, 2.0],
    [1.0, -Inf64, 2.5],[1.0, -1.0, 5.0],
    LibWorhp.WORHP_QUADRATIC,
    Vector{Int32}([LibWorhp.WORHP_QUADRATIC, LibWorhp.WORHP_LINEAR, LibWorhp.WORHP_LINEAR])

    )

    println("Before solveProblem")

solveProblem(p, UserF, UserG, UserDF, UserDG,
UserHM)

f = p.opt.F
sol = unsafe_wrap(Array{Float64,1}, p.opt.X,4)

@test abs(f + 0.5) < 1.0e-6
@test sum(abs, sol - [0, 0.5, 1., 2.]) < 1e-6
