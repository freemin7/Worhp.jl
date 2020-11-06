import MathOptInterface
const MOI = MathOptInterface
const MOIU = MathOptInterface.Utilities

mutable struct VariableInfo
    lower_bound::Float64  # May be -Inf even if has_lower_bound == true
    has_lower_bound::Bool # Implies lower_bound == Inf
    lower_bound_dual_start::Union{Nothing, Float64}
    upper_bound::Float64  # May be Inf even if has_upper_bound == true
    has_upper_bound::Bool # Implies upper_bound == Inf
    is_fixed::Bool        # Implies lower_bound == upper_bound and !has_lower_bound and !has_upper_bound.
    start::Union{Nothing, Float64}
    name::String
end

VariableInfo() = VariableInfo(-Inf, false, nothing, Inf, false, false, nothing, "")


mutable struct ConstraintInfo{F, S}
    func::F
    set::S
end

mutable struct Optimizer <: MOI.AbstractOptimizer
    inner::Union{WorhpProblem,Nothing}

    variable_info::Vector{VariableInfo}
    nlp_data::MOI.NLPBlockData
    sense::MOI.OptimizationSense
    objective::Union{MOI.SingleVariable,MOI.ScalarAffineFunction{Float64},MOI.ScalarQuadraticFunction{Float64},Nothing}
    linear_interval_constraints::Vector{ConstraintInfo{MOI.ScalarAffineFunction{Float64}, MOI.Interval{Float64}} }
    quadratic_interval_constrains::Vector{ConstraintInfo{MOI.ScalarQuadraticFunction{Float64}, MOI.Interval{Float64}}}


    # Parameters.
    silent::Bool
    options::Dict{String, Any}

    # Solution attributes.
    solve_time::Float64

end


struct EmptyNLPEvaluator <: MOI.AbstractNLPEvaluator end
MOI.features_available(::EmptyNLPEvaluator) = [:Grad, :Jac, :Hess]
MOI.initialize(::EmptyNLPEvaluator, features) = nothing
MOI.eval_objective(::EmptyNLPEvaluator, x) = NaN
function MOI.eval_constraint(::EmptyNLPEvaluator, g, x)
    @assert length(g) == 0
    return
end
function MOI.eval_objective_gradient(::EmptyNLPEvaluator, g, x)
    fill!(g, 0.0)
    return
end
MOI.jacobian_structure(::EmptyNLPEvaluator) = Tuple{Int64,Int64}[]
MOI.hessian_lagrangian_structure(::EmptyNLPEvaluator) = Tuple{Int64,Int64}[]
function MOI.eval_constraint_jacobian(::EmptyNLPEvaluator, J, x)
    @assert length(J) == 0
    return
end
function MOI.eval_hessian_lagrangian(::EmptyNLPEvaluator, H, x, σ, μ)
    @assert length(H) == 0
    return
end


empty_nlp_data() = MOI.NLPBlockData([], EmptyNLPEvaluator(), false)

function Optimizer(;options...)
    options_dict = Dict{String, Any}()
    # TODO: Setting options through the constructor could be deprecated in the
    # future.
    for (name, value) in options
        options_dict[string(name)] = value
    end
    return Optimizer(nothing, [], empty_nlp_data(), MOI.FEASIBILITY_SENSE,
                     nothing, [], [], false, options_dict, NaN)
end

MOI.supports(::Optimizer, ::MOI.NLPBlock) = true

function MOI.supports(::Optimizer,
                      ::MOI.ObjectiveFunction{MOI.SingleVariable})
    return true
end

function MOI.supports(::Optimizer,
    ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}})
    return true
end

function MOI.supports(::Optimizer,
    ::MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{Float64}})
    return true
end


function MOI.supports(::Optimizer, ::MOI.VariablePrimalStart,
                      ::Type{MOI.VariableIndex})
    return false
end

MOI.supports(::Optimizer, ::MOI.ObjectiveSense) = true

MOI.supports(::Optimizer, ::MOI.Silent) = true

const TIME_LIMIT = "max_cpu_time"
MOI.supports(::Optimizer, ::MOI.TimeLimitSec) = true
function MOI.set(model::Optimizer, ::MOI.TimeLimitSec, value::Real)
    MOI.set(model, MOI.RawParameter(TIME_LIMIT), Float64(value))
end
function MOI.set(model::Optimizer, attr::MOI.TimeLimitSec, ::Nothing)
    delete!(model.options, TIME_LIMIT)
end
function MOI.get(model::Optimizer, ::MOI.TimeLimitSec)
    return get(model.options, TIME_LIMIT, nothing)
end


function MOI.set(model::Optimizer, p::MOI.RawParameter, value)
    model.options[p.name] = value
    return
end

function MOI.get(model::Optimizer, p::MOI.RawParameter)
    if haskey(model.options, p.name)
        return model.options[p.name]
    end
    error("RawParameter with name $(p.name) is not set.")
end

MOI.get(model::Optimizer, ::MOI.SolveTime) = model.solve_time

MOI.supports_constraint(::Optimizer, ::Type{MOI.SingleVariable}, ::Type{MOI.Interval{Float64}}) = true
MOI.supports_constraint(::Optimizer, ::Type{MOI.ScalarAffineFunction{Float64}}, ::Type{MOI.Interval{Float64}}) = true
MOI.supports_constraint(::Optimizer, ::Type{MOI.ScalarQuadraticFunction{Float64}}, ::Type{MOI.Interval{Float64}}) = true


MOI.supports(::Optimizer, ::MOI.VariableName, ::Type{MOI.VariableIndex}) = true
function MOI.set(model::Optimizer, ::MOI.VariableName, vi::MOI.VariableIndex, name::String)
    model.variable_info[vi.value].name = name
    return
end
function MOI.get(model::Optimizer, ::MOI.VariableName, vi::MOI.VariableIndex)
    return model.variable_info[vi.value].name
end
function MOI.get(model::Optimizer, ::MOI.VariableIndex, vn::String)
    for i in eachindex(model.variable_info)
        if model.variable_info[i].name == vn
            return MOI.VariableIndex(i)
        end
    end
    return nothing
end


# supports_default_copy_to
MOIU.supports_default_copy_to(::Optimizer, copy_names::Bool) = !copy_names

function MOI.copy_to(model::Optimizer, src::MOI.ModelLike; copy_names = false)
    return MOIU.default_copy_to(model, src, copy_names)
end

MOI.get(::Optimizer, ::MOI.SolverName) = "Worhp"

MOI.get(model::Optimizer, ::MOI.ObjectiveFunctionType) = typeof(model.objective)

MOI.get(model::Optimizer, ::MOI.NumberOfVariables) = length(model.variable_info)


MOI.get(model::Optimizer, ::MOI.NumberOfConstraints{MOI.ScalarQuadraticFunction{Float64}, MOI.Interval{Float64}}) = length(model.quadratic_interval_constraints)
MOI.get(model::Optimizer, ::MOI.NumberOfConstraints{MOI.ScalarAffineFunction{Float64}, MOI.Interval{Float64}}) = length(model.linear_interval_constraints)
MOI.get(model::Optimizer, ::MOI.NumberOfConstraints{MOI.SingleVariable, MOI.Interval{Float64}}) = length(model.variable_info)


function MOI.get(model::Optimizer, ::MOI.ListOfVariableIndices)
    return [MOI.VariableIndex(i) for i in 1:length(model.variable_info)]
end

function MOI.get(model::Optimizer, ::MOI.ListOfConstraints)
    constraints = Set{Tuple{DataType, MOI.Interval{Float64}}}()
    for info in model.variable_info
        push!(constraints, (MOI.SingleVariable, MOI.Interval{Float64}))
    end

    if !isempty(model.linear_interval_constraints)
        push!(constraints, (MOI.ScalarAffineFunction{Float64}, MOI.Interval{Float64}))
    end
    if !isempty(model.quadratic_interval_constraints)
        push!(constraints, (MOI.ScalarQuadraticFunction{Float64}, MOI.Interval{Float64}))
    end

    return collect(constraints)
end


function MOI.get(
    model::Optimizer,
    ::MOI.ListOfConstraintIndices{MOI.ScalarAffineFunction{Float64}, MOI.Interval{Float64}}
)
    return MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}, MOI.Interval{Float64}}.(eachindex(model.linear_interval_constraints))
end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintFunction,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}, MOI.Interval{Float64}}
)
    return model.linear_interval_constraints[c.value].func
end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}, MOI.Interval{Float64}}
)
    return model.linear_interval_constraints[c.value].set
end

## No ScalarQuadraticFunction?

function MOI.get(model::Optimizer, ::MOI.ListOfConstraintIndices{MOI.SingleVariable, MOI.Interval{Float64}})
    return MOI.ConstraintIndex{MOI.SingleVariable, MOI.Interval{Float64}}.(1:length(model.variable_info))
end

# function MOI.get(model::Optimizer, ::MOI.ListOfConstraintIndices{MOI.SingleVariable, MOI.LessThan{Float64}})
#     dict = Dict(model.variable_info[i] => i for i in 1:length(model.variable_info))
#     filter!(info -> info.first.has_upper_bound, dict)
#     return MOI.ConstraintIndex{MOI.SingleVariable, MOI.LessThan{Float64}}.(values(dict))
# end
#
# function MOI.get(model::Optimizer, ::MOI.ListOfConstraintIndices{MOI.SingleVariable, MOI.EqualTo{Float64}})
#     dict = Dict(model.variable_info[i] => i for i in 1:length(model.variable_info))
#     filter!(info -> info.first.is_fixed, dict)
#     return MOI.ConstraintIndex{MOI.SingleVariable, MOI.EqualTo{Float64}}.(values(dict))
# end
#
# function MOI.get(model::Optimizer, ::MOI.ListOfConstraintIndices{MOI.SingleVariable, MOI.GreaterThan{Float64}})
#     dict = Dict(model.variable_info[i] => i for i in 1:length(model.variable_info))
#     filter!(info -> info.first.has_lower_bound, dict)
#     return MOI.ConstraintIndex{MOI.SingleVariable, MOI.GreaterThan{Float64}}.(values(dict))
# end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintFunction,
    c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.Interval{Float64}}
)
    return MOI.SingleVariable(MOI.VariableIndex(c.value))
end

# function MOI.get(
#     model::Optimizer,
#     ::MOI.ConstraintFunction,
#     c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.Interval{Float64}}
# )
#     return MOI.SingleVariable(MOI.VariableIndex(c.value))
# end
#
# function MOI.get(
#     model::Optimizer,
#     ::MOI.ConstraintFunction,
#     c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.EqualTo{Float64}}
# )
#     return MOI.SingleVariable(MOI.VariableIndex(c.value))
# end
#
# function MOI.get(
#     model::Optimizer,
#     ::MOI.ConstraintFunction,
#     c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.GreaterThan{Float64}}
# )
#     return MOI.SingleVariable(MOI.VariableIndex(c.value))
# end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.Interval{Float64}}
)
    return MOI.Interval{Float64}(model.variable_info[c.value].lower_bound, model.variable_info[c.value].upper_bound)
end

# function MOI.get(
#     model::Optimizer,
#     ::MOI.ConstraintSet,
#     c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.LessThan{Float64}}
# )
#     return MOI.LessThan{Float64}(model.variable_info[c.value].upper_bound)
# end
#
# function MOI.get(
#     model::Optimizer,
#     ::MOI.ConstraintSet,
#     c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.EqualTo{Float64}}
# )
#     return MOI.EqualTo{Float64}(model.variable_info[c.value].lower_bound)
# end
#
# function MOI.get(
#     model::Optimizer,
#     ::MOI.ConstraintSet,
#     c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.GreaterThan{Float64}}
# )
#     return MOI.GreaterThan{Float64}(model.variable_info[c.value].lower_bound)
# end


function MOI.get(
    model::Optimizer,
    ::MOI.ObjectiveFunction
)
    return model.objective
end


function MOI.set(model::Optimizer, ::MOI.ObjectiveSense,
                 sense::MOI.OptimizationSense)
    model.sense = sense
    return
end

MOI.get(model::Optimizer, ::MOI.ObjectiveSense) = model.sense

function MOI.set(model::Optimizer, ::MOI.Silent, value)
    model.silent = value
    return
end

MOI.get(model::Optimizer, ::MOI.Silent) = model.silent

# TODO Time limit

# TODO RawParameter

##

function MOI.empty!(model::Optimizer)
    model.inner = nothing
    empty!(model.variable_info)
    model.nlp_data = empty_nlp_data()
    model.sense = MOI.FEASIBILITY_SENSE
    model.objective = nothing
    empty!(model.linear_interval_constraints)
    empty!(model.quadratic_interval_constrains)
end

function MOI.is_empty(model::Optimizer)
    return isempty(model.variable_info) &&
           model.nlp_data.evaluator isa EmptyNLPEvaluator &&
           model.sense == MOI.FEASIBILITY_SENSE &&
           isempty(model.linear_interval_constraints) &&
           isempty(model.quadratic_interval_constrains)
end

function MOI.add_variable(model::Optimizer)
    push!(model.variable_info, VariableInfo())
    return MOI.VariableIndex(length(model.variable_info))
end
function MOI.add_variables(model::Optimizer, n::Int)
    return [MOI.add_variable(model) for i in 1:n]
end

MOI.is_valid(model::Optimizer, vi::MOI.VariableIndex) = vi.value in eachindex(model.variable_info)

function MOI.is_valid(model::Optimizer, ci::MOI.ConstraintIndex{MOI.SingleVariable, MOI.Interval{Float64}})
    vi = MOI.VariableIndex(ci.value)
    return MOI.is_valid(model, vi)
end

# function MOI.is_valid(model::Optimizer, ci::MOI.ConstraintIndex{MOI.SingleVariable, MOI.LessThan{Float64}})
#     vi = MOI.VariableIndex(ci.value)
#     return MOI.is_valid(model, vi) && has_upper_bound(model, vi)
# end
# function MOI.is_valid(model::Optimizer, ci::MOI.ConstraintIndex{MOI.SingleVariable, MOI.GreaterThan{Float64}})
#     vi = MOI.VariableIndex(ci.value)
#     return MOI.is_valid(model, vi) && has_lower_bound(model, vi)
# end
# function MOI.is_valid(model::Optimizer, ci::MOI.ConstraintIndex{MOI.SingleVariable, MOI.EqualTo{Float64}})
#     vi = MOI.VariableIndex(ci.value)
#     return MOI.is_valid(model, vi) && is_fixed(model, vi)
# end

function check_inbounds(model::Optimizer, var::MOI.SingleVariable)
    return MOI.throw_if_not_valid(model, var.variable)
end

function check_inbounds(model::Optimizer, aff::MOI.ScalarAffineFunction)
    for term in aff.terms
        MOI.throw_if_not_valid(model, term.variable_index)
    end
end

function check_inbounds(model::Optimizer, quad::MOI.ScalarQuadraticFunction)
    for term in quad.affine_terms
        MOI.throw_if_not_valid(model, term.variable_index)
    end
    for term in quad.quadratic_terms
        MOI.throw_if_not_valid(model, term.variable_index_1)
        MOI.throw_if_not_valid(model, term.variable_index_2)
    end
end

function has_upper_bound(model::Optimizer, vi::MOI.VariableIndex)
    return model.variable_info[vi.value].has_upper_bound
end

function has_lower_bound(model::Optimizer, vi::MOI.VariableIndex)
    return model.variable_info[vi.value].has_lower_bound
end

function is_fixed(model::Optimizer, vi::MOI.VariableIndex)
    return model.variable_info[vi.value].is_fixed
end

# function MOI.add_constraint(model::Optimizer, v::MOI.SingleVariable, lt::MOI.LessThan{Float64})
#     vi = v.variable
#     MOI.throw_if_not_valid(model, vi)
#     if isnan(lt.upper)
#         error("Invalid upper bound value $(lt.upper).")
#     end
#     if has_upper_bound(model, vi)
#         throw(MOI.UpperBoundAlreadySet{typeof(lt), typeof(lt)}(vi))
#     end
#     if is_fixed(model, vi)
#         throw(MOI.UpperBoundAlreadySet{MOI.EqualTo{Float64}, typeof(lt)}(vi))
#     end
#     model.variable_info[vi.value].upper_bound = lt.upper
#     model.variable_info[vi.value].has_upper_bound = true
#     return MOI.ConstraintIndex{MOI.SingleVariable, MOI.LessThan{Float64}}(vi.value)
# end
#
# function MOI.set(model::Optimizer, ::MOI.ConstraintSet,
#                  ci::MOI.ConstraintIndex{MOI.SingleVariable, MOI.LessThan{Float64}},
#                  set::MOI.LessThan{Float64})
#     MOI.throw_if_not_valid(model, ci)
#     model.variable_info[ci.value].upper_bound = set.upper
#     return
# end
#
# function MOI.delete(model::Optimizer, ci::MOI.ConstraintIndex{MOI.SingleVariable, MOI.LessThan{Float64}})
#     MOI.throw_if_not_valid(model, ci)
#     model.variable_info[ci.value].upper_bound = Inf
#     model.variable_info[ci.value].has_upper_bound = false
#     return
# end
#
# function MOI.add_constraint(model::Optimizer, v::MOI.SingleVariable, gt::MOI.GreaterThan{Float64})
#     vi = v.variable
#     MOI.throw_if_not_valid(model, vi)
#     if isnan(gt.lower)
#         error("Invalid lower bound value $(gt.lower).")
#     end
#     if has_lower_bound(model, vi)
#         throw(MOI.LowerBoundAlreadySet{typeof(gt), typeof(gt)}(vi))
#     end
#     if is_fixed(model, vi)
#         throw(MOI.LowerBoundAlreadySet{MOI.EqualTo{Float64}, typeof(gt)}(vi))
#     end
#     model.variable_info[vi.value].lower_bound = gt.lower
#     model.variable_info[vi.value].has_lower_bound = true
#     return MOI.ConstraintIndex{MOI.SingleVariable, MOI.GreaterThan{Float64}}(vi.value)
# end
#
# function MOI.set(model::Optimizer, ::MOI.ConstraintSet,
#                  ci::MOI.ConstraintIndex{MOI.SingleVariable, MOI.GreaterThan{Float64}},
#                  set::MOI.GreaterThan{Float64})
#     MOI.throw_if_not_valid(model, ci)
#     model.variable_info[ci.value].lower_bound = set.lower
#     return
# end
#
# function MOI.delete(model::Optimizer, ci::MOI.ConstraintIndex{MOI.SingleVariable, MOI.GreaterThan{Float64}})
#     MOI.throw_if_not_valid(model, ci)
#     model.variable_info[ci.value].lower_bound = -Inf
#     model.variable_info[ci.value].has_lower_bound = false
#     return
# end
#
# function MOI.add_constraint(model::Optimizer, v::MOI.SingleVariable, eq::MOI.EqualTo{Float64})
#     vi = v.variable
#     MOI.throw_if_not_valid(model, vi)
#     if isnan(eq.value)
#         error("Invalid fixed value $(eq.value).")
#     end
#     if has_lower_bound(model, vi)
#         throw(MOI.LowerBoundAlreadySet{MOI.GreaterThan{Float64}, typeof(eq)}(vi))
#     end
#     if has_upper_bound(model, vi)
#         throw(MOI.UpperBoundAlreadySet{MOI.LessThan{Float64}, typeof(eq)}(vi))
#     end
#     if is_fixed(model, vi)
#         throw(MOI.LowerBoundAlreadySet{typeof(eq), typeof(eq)}(vi))
#     end
#     model.variable_info[vi.value].lower_bound = eq.value
#     model.variable_info[vi.value].upper_bound = eq.value
#     model.variable_info[vi.value].is_fixed = true
#     return MOI.ConstraintIndex{MOI.SingleVariable, MOI.EqualTo{Float64}}(vi.value)
# end
#
# function MOI.add_constraint(model::Optimizer, v::MOI.SingleVariable, eq::MOI.EqualTo{Float64})
#     vi = v.variable
#     MOI.throw_if_not_valid(model, vi)
#     if isnan(eq.value)
#         error("Invalid fixed value $(gt.lower).")
#     end
#     if has_lower_bound(model, vi)
#         throw(MOI.LowerBoundAlreadySet{MOI.GreaterThan{Float64}, typeof(eq)}(vi))
#     end
#     if has_upper_bound(model, vi)
#         throw(MOI.UpperBoundAlreadySet{MOI.LessThan{Float64}, typeof(eq)}(vi))
#     end
#     if is_fixed(model, vi)
#         throw(MOI.LowerBoundAlreadySet{typeof(eq), typeof(eq)}(vi))
#     end
#     model.variable_info[vi.value].lower_bound = eq.value
#     model.variable_info[vi.value].upper_bound = eq.value
#     model.variable_info[vi.value].is_fixed = true
#     return MOI.ConstraintIndex{MOI.SingleVariable, MOI.EqualTo{Float64}}(vi.value)
# end
#
# function MOI.set(model::Optimizer, ::MOI.ConstraintSet,
#                  ci::MOI.ConstraintIndex{MOI.SingleVariable, MOI.EqualTo{Float64}},
#                  set::MOI.EqualTo{Float64})
#     MOI.throw_if_not_valid(model, ci)
#     model.variable_info[ci.value].lower_bound = set.value
#     model.variable_info[ci.value].upper_bound = set.value
#     return
# end
#
# function MOI.delete(model::Optimizer, ci::MOI.ConstraintIndex{MOI.SingleVariable, MOI.EqualTo{Float64}})
#     MOI.throw_if_not_valid(model, ci)
#     model.variable_info[ci.value].lower_bound = -Inf
#     model.variable_info[ci.value].upper_bound = Inf
#     model.variable_info[ci.value].is_fixed = false
#     return
# end

function MOI.add_constraint(model::Optimizer, v::MOI.SingleVariable, int::MOI.Interval{Float64})
    vi = v.variable
    MOI.throw_if_not_valid(model, vi)
    if isnan(int.upper) || isnan(int.lower) || int.lower >= int.upper
        error("Invalid fixed value ($(int.lower), $(int.upper)) .")
    end
    if has_lower_bound(model, vi)
        throw(MOI.LowerBoundAlreadySet{MOI.GreaterThan{Float64}, typeof(int)}(vi))
    end
    if has_upper_bound(model, vi)
        throw(MOI.UpperBoundAlreadySet{MOI.LessThan{Float64}, typeof(int)}(vi))
    end
    if is_fixed(model, vi)
        throw(MOI.LowerBoundAlreadySet{typeof(int), typeof(int)}(vi))
    end
    model.variable_info[vi.value].lower_bound = int.lower
    model.variable_info[vi.value].upper_bound = int.upper
    return MOI.ConstraintIndex{MOI.SingleVariable, MOI.EqualTo{Float64}}(vi.value)
end

function MOI.set(model::Optimizer, ::MOI.ConstraintSet,
                 ci::MOI.ConstraintIndex{MOI.SingleVariable, MOI.Interval{Float64}},
                 set::MOI.Interval{Float64})
    MOI.throw_if_not_valid(model, ci)
    model.variable_info[ci.value].lower_bound = set.lower
    model.variable_info[ci.value].upper_bound = set.upper
    return
end

function MOI.delete(model::Optimizer, ci::MOI.ConstraintIndex{MOI.SingleVariable, MOI.Interval{Float64}})
    MOI.throw_if_not_valid(model, ci)
    model.variable_info[ci.value].lower_bound = -Inf
    model.variable_info[ci.value].upper_bound = Inf

    model.variable_info[ci.value].has_lower_bound = false
    model.variable_info[ci.value].has_upper_bound = false
    model.variable_info[ci.value].is_fixed = false
    return
end

macro define_add_constraint(function_type, set_type, prefix)
    array_name = Symbol(string(prefix) * "_constraints")
    quote
        function MOI.add_constraint(model::Optimizer, func::$function_type, set::$set_type)
            check_inbounds(model, func)
            push!(model.$(array_name), ConstraintInfo(func, set))
            return MOI.ConstraintIndex{$function_type, $set_type}(length(model.$(array_name)))
        end
    end
end


@define_add_constraint(MOI.ScalarAffineFunction{Float64}, MOI.Interval{Float64},
                       linear_interval)
@define_add_constraint(MOI.ScalarQuadraticFunction{Float64},
                       MOI.Interval{Float64}, quadratic_interval)

# no dual constraints i think

function MOI.supports(::Optimizer, ::MOI.NLPBlockDualStart)
   return true
end
function MOI.set(model::Optimizer, ::MOI.NLPBlockDualStart, values)
   model.nlp_dual_start = -values
   return
end

function MOI.set(model::Optimizer, ::MOI.NLPBlock, nlp_data::MOI.NLPBlockData)
   model.nlp_data = nlp_data
   return
end

function MOI.set(model::Optimizer, ::MOI.ObjectiveFunction,
                func::Union{MOI.SingleVariable, MOI.ScalarAffineFunction,
                            MOI.ScalarQuadraticFunction})
   check_inbounds(model, func)
   model.objective = func
   return
end
