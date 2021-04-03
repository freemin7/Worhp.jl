module WorhpOpt


baremodule LibWorhp
    using CBinding: 𝐣𝐥

    const size_t = 𝐣𝐥.Csize_t

    𝐣𝐥.Base.include((𝐣𝐥.@__MODULE__), 𝐣𝐥.joinpath(𝐣𝐥.dirname(𝐣𝐥.@__DIR__), "deps", "bindings.jl"))  #

    export TerminateSuccess, TerminateError, GetUserAction
    export callWorhp, iterOutput, evalF, evalG, evalDF, evalDG, evalHM, fidif, WorhpFidif
    export DoneUserAction, Worhp, IterationOutput, StatusMsg
    export OptVarStruct, Workspace, Params, Control
end

using .LibWorhp

function versionCheck()
    if 0 != LibWorhp.CheckWorhpVersion(LibWorhp.WORHP_MAJOR, LibWorhp.WORHP_MINOR, LibWorhp.WORHP_PATCH)
        @error "Worhp version of headers doesn't match shared object. If you upgraded Worhp run ]build again"
    end
end


mutable struct WorhpProblem
    opt::LibWorhp.OptVarStruct
    wsp::LibWorhp.Workspace
    par::LibWorhp.Params
    cnt::LibWorhp.Control
    function WorhpProblem(
        n::Int32,m::Int32,DFnz::Int32,DGnz::Int32,HMnz::Int32, 
        has_df::Bool, has_dg::Bool, has_hm::Bool,
        dfrow::Vector{Int32},
        dgrow::Vector{Int32}, dgcol::Vector{Int32},
        hmrow::Vector{Int32}, hmcol::Vector{Int32},
        #Multiple dispatch can be used to make that nice, such as passing a matrix
        X0::Vector{Float64},
        Lambda::Vector{Float64},Mu::Vector{Float64},
        XL::Vector{Float64},XU::Vector{Float64},
        GL::Vector{Float64},GU::Vector{Float64},
        F_Type::Int32,
        G_Type::Vector{Int32};
        config_path::String="worhp.xml"
        )

        opt = LibWorhp.OptVarStruct(undef);
        wsp = LibWorhp.Workspace(undef);
        par = LibWorhp.Params(undef);
        cnt = LibWorhp.Control(undef);

        status = Ref{Cint}(123) #can't be undef, i tried.

        versionCheck();

        LibWorhp.WorhpPreInit(Ref(opt), Ref(wsp), Ref(par), Ref(cnt))
        LibWorhp.InitParams(status, Ref(par))

        LibWorhp.ReadParamsNoInit(status, config_path, Ref(par))

        if (status[] == LibWorhp.DataError || status[] == LibWorhp.InitErr)
            @error "Initalization failed with error $(status[])"
        end

        opt.n = n;
        opt.m = m;
        wsp.DF.nnz = DFnz;
        wsp.DG.nnz = DGnz;
        wsp.HM.nnz = HMnz;

        par.UserHMstructure = Int32(1);  # Worph is allowed to remove wrong hessian values

        par.UserDF = has_df
        par.UserDG = has_dg
        par.UserHM = has_hm

 

        LibWorhp.WorhpInit(Ref(opt), Ref(wsp), Ref(par), Ref(cnt))

        opt.FType = F_Type

        if m>0
        if length(G_Type) == m 
            unsafe_wrap(Vector{Int32},opt.GType,m) .= G_Type
        else
            @error "G_Type has wrong length"
        end
        end

        if cnt.status != LibWorhp.FirstCall
            @error "Main: Initialisation failed with status $(cnt.status)"
        end

        unsafe_wrap(Vector{Float64},opt.X,n) .= X0;
        unsafe_wrap(Vector{Float64},opt.XL,n) .= clamp.(XL,-par.Infty,par.Infty);
        unsafe_wrap(Vector{Float64},opt.XU,n) .= clamp.(XU,-par.Infty,par.Infty);;
        unsafe_wrap(Vector{Float64},opt.Lambda,n) .= Lambda;

        if m>0
            unsafe_wrap(Vector{Float64},opt.GL,m) .= clamp.(GL,-par.Infty,par.Infty);
            unsafe_wrap(Vector{Float64},opt.GU,m) .= clamp.(GU,-par.Infty,par.Infty);
            unsafe_wrap(Vector{Float64},opt.Mu,m) .= Mu;
        end

        if wsp.DF.NeedStructure == true
            #= if DFnz == LibWorhp.WorhpMatrix_Init_Dense or DFnz >= n
                @warn "Using dense first derivative"
                unsafe_wrap(Vector{Int32},wsp.DF.row,n) .= 1:n
            else =#
            if length(dfrow) != DFnz
                @error "Length of sparse gradient doesn't match number of non zeros given"
            else
                unsafe_wrap(Vector{Int32},wsp.DF.row,DFnz) .= dfrow
            end
        end

        if wsp.DG.NeedStructure == true
            if length(dgrow) != length(dgcol)
                @error "Length of dgrow and dgcol don't match"
            elseif length(dgrow) != DGnz
                @error "length of sparse Jacobi matrix of constraints in CS for format doesn't match number of non zeros given"
            #=
            elseif DGnz == LibWorhp.WorhpMatrix_Init_Dense or DGnz >= n*m
                local i = 0
                for c in 1:m
                    for r in  1:n
                        i = i+1;
                        unsafe_store!(wsp.DG.row,r,i)
                        unsafe_store!(wsp.DG.col,c,i)
                    end
                end
            =#
            elseif !issorted(dgcol)
                @error "Jacobi matrix is not column sorted"
            else
                unsafe_wrap(Vector{Int32},wsp.DG.row,DGnz) .= dgrow
                unsafe_wrap(Vector{Int32},wsp.DG.col,DGnz) .= dgcol
            end
        end


        if wsp.HM.NeedStructure == true
            if length(hmrow) != length(hmcol)
                @error "Length of hmrow and hmcol don't match"
            elseif length(hmrow) != HMnz
                @error "length of sparse Hessian  matrix in CS for format doesn't match number of zeros given"
            elseif !issorted(hmcol)
                @error "Hessian matrix is not column sorted"
            else
                unsafe_wrap(Vector{Int32},wsp.HM.row,HMnz) .= hmrow
                unsafe_wrap(Vector{Int32},wsp.HM.col,HMnz) .= hmcol
            end
        end

                println("HM")

        res = new(opt, wsp, par, cnt)

        return finalizer(ww->LibWorhp.WorhpFree(Ref(ww.opt), Ref(ww.wsp), Ref(ww.par), Ref(ww.cnt)), res)
    end
end


function solveProblem(prob::WorhpProblem, UserF, UserG, UserDF, UserDG,
    UserHM)
    opt, wsp, par, cnt = prob.opt, prob.wsp, prob.par, prob.cnt;
    optR, wspR, parR, cntR = Ref(prob.opt), Ref(prob.wsp), Ref(prob.par), Ref(prob.cnt)

    while prob.cnt.status < LibWorhp.TerminateSuccess && prob.cnt.status > LibWorhp.TerminateError
        if GetUserAction(cntR, callWorhp) == true
            Worhp(optR,wspR,parR,cntR)
        end
        if GetUserAction(cntR, iterOutput) == true
            IterationOutput(optR,wspR,parR,cntR)
            DoneUserAction(cntR, iterOutput)
        end
        if  GetUserAction(cntR, evalF) == true
            UserF(opt, wsp, par, cnt)
            DoneUserAction(cntR, evalF)
        end
        if  GetUserAction(cntR, evalG) == true
            UserG(opt, wsp, par, cnt)
            DoneUserAction(cntR, evalG)
        end
        if GetUserAction(cntR, evalDF) == true
            UserDF(opt, wsp, par, cnt)
            DoneUserAction(cntR, evalDF)
        end
        if GetUserAction(cntR, evalDG) == true
            UserDG(opt, wsp, par, cnt)
            DoneUserAction(cntR, evalDG)
        end
        if GetUserAction(cntR, evalHM) == true
            UserHM(opt, wsp, par, cnt)
            DoneUserAction(cntR, evalHM)
        end
        if GetUserAction(cntR, fidif) == true
            WorhpFidif(optR,wspR,parR,cntR)
        end

    end

    StatusMsg(optR,wspR,parR,cntR)
end

export solveProblem, WorhpProblem, LibWorhp

end # module
