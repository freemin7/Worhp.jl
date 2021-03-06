using WorhpOpt
using Test
opt = LibWorhp.OptVarStruct(undef);
wsp = LibWorhp.Workspace(undef);
par = LibWorhp.Params(undef);
cnt = LibWorhp.Control(undef);

@test 0 == LibWorhp.CheckWorhpVersion(LibWorhp.WORHP_MAJOR, LibWorhp.WORHP_MINOR, LibWorhp.WORHP_PATCH)

status = Ref{Cint}(123);

optR, wspR, parR, cntR = Ref(opt), Ref(wsp), Ref(par), Ref(cnt)

LibWorhp.WorhpPreInit(optR, wspR, parR, cntR)

LibWorhp.InitParams(status, parR)

par.NLPprint = 1

LibWorhp.ReadParamsNoInit(status, "worhp.xml", parR)

@test status[] != LibWorhp.DataError
@test status[] != LibWorhp.InitErr

opt.n = 4
opt.m = 3

wsp.DF.nnz = 3
wsp.DG.nnz = 6
wsp.HM.nnz = 1 + opt.n

LibWorhp.WorhpInit(optR,wspR,parR,cntR)

@test LibWorhp.FirstCall == cnt.status

unsafe_wrap(Array{Float64,1}, opt.X, 4; own=false) .= [2.0, 2.0, 1.0, 0.0]

unsafe_wrap(Array{Float64,1}, opt.Lambda, 4; own=false) .= [0.0, 0.0, 0.0, 0.0]

unsafe_wrap(Array{Float64,1}, opt.Mu, 3; own=false) .= [0.0, 0.0, 0.0]

unsafe_wrap(Array{Float64,1}, opt.XL, 4; own=false) .= [-0.5, -2.0, 0.0, -2.0]

unsafe_wrap(Array{Float64,1}, opt.XU, 4; own=false) .= [par.Infty, par.Infty, 2.0, 2.0]

unsafe_wrap(Array{Float64,1}, opt.GL, 3; own=false) .= [1.0, -par.Infty, 2.5]

unsafe_wrap(Array{Float64,1}, opt.GU, 3; own=false).= [1.0, -1.0, 5.0]

if wsp.DF.NeedStructure == true
    row = unsafe_wrap(Array{Int32,1}, wsp.DF.row, 3; own=false)
    row .= [1,2,3]
end
if wsp.DG.NeedStructure == true
    row = unsafe_wrap(Array{Int32,1}, wsp.DG.row, 6; own=false)
    col = unsafe_wrap(Array{Int32,1}, wsp.DG.col, 6; own=false)
    row .= [1,3,1,2,2,3]
    col .= [1,2,3,3,4,4]
end

if wsp.HM.NeedStructure == true
    row = unsafe_wrap(Array{Int32,1}, wsp.HM.row, wsp.HM.nnz; own=false)
    col = unsafe_wrap(Array{Int32,1}, wsp.HM.col, wsp.HM.nnz; own=false)

    row[1] = 3
    col[1] = 1
    for i in 2:5
        row[i] = i-1
        col[i] = i-1
    end
end

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


while cnt.status < LibWorhp.TerminateSuccess && cnt.status > LibWorhp.TerminateError

    if LibWorhp.GetUserAction(cntR, LibWorhp.callWorhp) == true
        LibWorhp.Worhp(optR,wspR,parR,cntR)
    end
    if LibWorhp.GetUserAction(cntR, LibWorhp.iterOutput) == true
        LibWorhp.IterationOutput(optR,wspR,parR,cntR)
        LibWorhp.DoneUserAction(cntR, LibWorhp.iterOutput)
    end
    if LibWorhp.GetUserAction(cntR, LibWorhp.evalF) == true
        UserF(opt, wsp, par, cnt)
        LibWorhp.DoneUserAction(cntR,LibWorhp.evalF)
    end
    if LibWorhp.GetUserAction(cntR, LibWorhp.evalG) == true
        UserG(opt, wsp, par, cnt)
        LibWorhp.DoneUserAction(cntR, LibWorhp.evalG)
    end
    if LibWorhp.GetUserAction(cntR, LibWorhp.evalDF) == true
        UserDF(opt, wsp, par, cnt)
        LibWorhp.DoneUserAction(cntR, LibWorhp.evalDF)
    end
    if LibWorhp.GetUserAction(cntR, LibWorhp.evalDG) == true
        UserDG(opt, wsp, par, cnt)
        LibWorhp.DoneUserAction(cntR, LibWorhp.evalDG)
    end
    if LibWorhp.GetUserAction(cntR, LibWorhp.evalHM) == true
        UserHM(opt, wsp, par, cnt)
        LibWorhp.DoneUserAction(cntR, LibWorhp.evalHM)
    end
    if LibWorhp.GetUserAction(cntR, LibWorhp.fidif) == true
        LibWorhp.WorhpFidif(optR,wspR,parR,cntR)
    end

end

LibWorhp.StatusMsg(optR,wspR,parR,cntR)

@test abs(opt.F + 0.5) < 10e-6 

@test sum(abs,unsafe_wrap(Array{Float64,1}, opt.X,4) - [0, 0.5, 1., 2.]) < 10e-6
LibWorhp.WorhpFree(optR, wspR, parR, cntR);
