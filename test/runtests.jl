using Worhp
using Test
opt = Worhp.LibWorhp.OptVarStruct(undef);
wsp = Worhp.LibWorhp.Workspace(undef);
par = Worhp.LibWorhp.Params(undef);
cnt = Worhp.LibWorhp.Control(undef);

GC.@preserve opt wsp par cnt begin
@test 0 == Worhp.LibWorhp.CheckWorhpVersion(Worhp.LibWorhp.WORHP_MAJOR, Worhp.LibWorhp.WORHP_MINOR, Worhp.LibWorhp.WORHP_PATCH)

status = Ref{Cint}(123);

optR, wspR, parR, cntR = pointer_from_objref(opt), pointer_from_objref(wsp), pointer_from_objref(par), pointer_from_objref(cnt)

Worhp.LibWorhp.WorhpPreInit(optR, wspR, parR, cntR)
Worhp.LibWorhp.InitParams(status, parR)

par.NLPprint = 1

Worhp.LibWorhp.ReadParamsNoInit(status, "worhp.xml", parR)

@test status[] != Worhp.LibWorhp.DataError
@test status[] != Worhp.LibWorhp.InitErr

opt.n = 4
opt.m = 3

wsp.DF.nnz = 3
wsp.DG.nnz = 6
wsp.HM.nnz = 1 + opt.n

Worhp.LibWorhp.WorhpInit(optR,wspR,parR,cntR)

@test Worhp.LibWorhp.FirstCall == cnt.status

unsafe_wrap(Array{Float64,1}, opt.X, 4; own=false) .= [2.0, 2.0, 1.0, 0.0]

unsafe_wrap(Array{Float64,1}, opt.Lambda, 4; own=false) .= [0.0, 0.0, 0.0, 0.0]

unsafe_wrap(Array{Float64,1}, opt.Mu, 3; own=false) .= [0.0, 0.0, 0.0]

unsafe_wrap(Array{Float64,1}, opt.XL, 4; own=false) .= [-0.5, -2.0, 0.0, -2.0]

unsafe_wrap(Array{Float64,1}, opt.XU, 4; own=false) .= [par.Infty, par.Infty, 2.0, 2.0]

unsafe_wrap(Array{Float64,1}, opt.GL, 3; own=false) .= [1.0, -par.Infty, 2.5]

unsafe_wrap(Array{Float64,1}, opt.GU, 3; own=false).= [1.0, -1.0, 5.0]

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
    G[3] = x[1] + x[4]
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


while cnt.status < Worhp.LibWorhp.TerminateSuccess && cnt.status > Worhp.LibWorhp.TerminateError

    if Worhp.LibWorhp.GetUserAction(cntR, Worhp.LibWorhp.callWorhp) == true
        Worhp.LibWorhp.Worhp(optR,wspR,parR,cntR)
    end
    if Worhp.LibWorhp.GetUserAction(cntR, Worhp.LibWorhp.iterOutput) == true
        Worhp.LibWorhp.IterationOutput(optR,wspR,parR,cntR)
        Worhp.LibWorhp.DoneUserAction(cntR,Worhp.LibWorhp.iterOutput)
    end
    if Worhp.LibWorhp.GetUserAction(cntR, Worhp.LibWorhp.evalF) == true
        UserF(opt, wsp, par, cnt)
        Worhp.LibWorhp.DoneUserAction(cntR,Worhp.LibWorhp.evalF)
    end
    if Worhp.LibWorhp.GetUserAction(cntR, Worhp.LibWorhp.evalG) == true
        UserG(opt, wsp, par, cnt)
        Worhp.LibWorhp.DoneUserAction(cntR,Worhp.LibWorhp.evalG)
    end
    if Worhp.LibWorhp.GetUserAction(cntR, Worhp.LibWorhp.evalDF) == true
        UserDF(opt, wsp, par, cnt)
        Worhp.LibWorhp.DoneUserAction(cntR,Worhp.LibWorhp.evalDF)
    end
    if Worhp.LibWorhp.GetUserAction(cntR, Worhp.LibWorhp.evalDG) == true
        UserDG(opt, wsp, par, cnt)
        Worhp.LibWorhp.DoneUserAction(cntR,Worhp.LibWorhp.evalDG)
    end
    if Worhp.LibWorhp.GetUserAction(cntR, Worhp.LibWorhp.evalHM) == true
        UserHM(opt, wsp, par, cnt)
        Worhp.LibWorhp.DoneUserAction(cntR,Worhp.LibWorhp.evalHM)
    end
    if Worhp.LibWorhp.GetUserAction(cntR, Worhp.LibWorhp.fidif) == true
        Worhp.LibWorhp.WorhpFidif(optR,wspR,parR,cntR)
    end

end

Worhp.LibWorhp.StatusMsg(optR,wspR,parR,cntR)
Worhp.LibWorhp.WorhpFree(optR, wspR, parR, cntR);

end
