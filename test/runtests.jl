using Worhp
using Test
opt = Worph.LibWorhp.OptVarStruct(undef);
wsp = Worph.LibWorhp.Workspace(undef);
par = Worph.LibWorhp.Params(undef);
cnt = Worph.LibWorhp.Control(undef);

@test 0 == Worph.LibWorhp.CheckWorhpVersion(Worph.LibWorhp.WORHP_MAJOR, Worph.LibWorhp.WORHP_MINOR, Worph.LibWorhp.WORHP_PATCH)

status = Ref{Cint}(123)

Worph.LibWorhp.WorhpPreInit(Ref(opt), Ref(wsp), Ref(par), Ref(cnt))
Worph.LibWorhp.InitParams(status, Ref(par))

par.NLPprint = 1

Worph.LibWorhp.ReadParamsNoInit(status, "worhp.xml", Ref(par))

@test status[] != Worph.LibWorhp.DataError
@test status[] != Worph.LibWorhp.InitErr

opt.n = 4
opt.m = 3

wsp.DF.nnz = 3
wsp.DG.nnz = 6
wsp.HM.nnz = 1 + opt.n

Worph.LibWorhp.WorhpInit(Ref(opt),Ref(wsp),Ref(par),Ref(cnt))

@test Worph.LibWorhp.FirstCall == cnt.status

X = unsafe_wrap(Array{Float64,1}, opt.X, 4; own=false)
X .= [2.0, 2.0, 1.0, 0.0]

Lambda = unsafe_wrap(Array{Float64,1}, opt.Lambda, 4; own=false)
Lambda .= [0.0, 0.0, 0.0, 0.0]

Mu = unsafe_wrap(Array{Float64,1}, opt.Mu, 3; own=false)
Mu .= [0.0, 0.0, 0.0]

XL = unsafe_wrap(Array{Float64,1}, opt.XL, 4; own=false)
XL .= [-0.5, -2.0, 0.0, -2.0]

XU = unsafe_wrap(Array{Float64,1}, opt.XU, 4; own=false)
XU .= [par.Infty, par.Infty, 2.0, 2.0]

GL = unsafe_wrap(Array{Float64,1}, opt.GL, 3; own=false)
GL .= [1.0, -par.Infty, 2.5]

GU = unsafe_wrap(Array{Float64,1}, opt.GU, 3; own=false)
GU .= [1.0, -1.0, 5.0]

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


while cnt.status < Worph.LibWorhp.TerminateSuccess && cnt.status > Worph.LibWorhp.TerminateError

    if Worph.LibWorhp.GetUserAction(Ref(cnt), Worph.LibWorhp.callWorhp) == true
        Worph.LibWorhp.Worhp(Ref(opt),Ref(wsp),Ref(par),Ref(cnt))
    end
    if Worph.LibWorhp.GetUserAction(Ref(cnt), Worph.LibWorhp.iterOutput) == true
        Worph.LibWorhp.IterationOutput(Ref(opt),Ref(wsp),Ref(par),Ref(cnt))
        Worph.LibWorhp.DoneUserAction(Ref(cnt),Worph.LibWorhp.iterOutput)
    end
    if Worph.LibWorhp.GetUserAction(Ref(cnt), Worph.LibWorhp.evalF) == true
        UserF(opt, wsp, par, cnt)
        Worph.LibWorhp.DoneUserAction(Ref(cnt),Worph.LibWorhp.evalF)
    end
    if Worph.LibWorhp.GetUserAction(Ref(cnt), Worph.LibWorhp.evalG) == true
        UserG(opt, wsp, par, cnt)
        Worph.LibWorhp.DoneUserAction(Ref(cnt),Worph.LibWorhp.evalG)
    end
    if Worph.LibWorhp.GetUserAction(Ref(cnt), Worph.LibWorhp.evalDF) == true
        UserDF(opt, wsp, par, cnt)
        Worph.LibWorhp.DoneUserAction(Ref(cnt),Worph.LibWorhp.evalDF)
    end
    if Worph.LibWorhp.GetUserAction(Ref(cnt), Worph.LibWorhp.evalDG) == true
        UserDG(opt, wsp, par, cnt)
        Worph.LibWorhp.DoneUserAction(Ref(cnt),Worph.LibWorhp.evalDG)
    end
    if Worph.LibWorhp.GetUserAction(Ref(cnt), Worph.LibWorhp.evalHM) == true
        UserHM(opt, wsp, par, cnt)
        Worph.LibWorhp.DoneUserAction(Ref(cnt),Worph.LibWorhp.evalHM)
    end
    if Worph.LibWorhp.GetUserAction(Ref(cnt), Worph.LibWorhp.fidif) == true
        Worph.LibWorhp.WorhpFidif(Ref(opt),Ref(wsp),Ref(par),Ref(cnt))
    end

end

Worph.LibWorhp.StatusMsg(Ref(opt),Ref(wsp),Ref(par),Ref(cnt))
Worph.LibWorhp.WorhpFree(Ref(opt), Ref(wsp), Ref(par), Ref(cnt));
