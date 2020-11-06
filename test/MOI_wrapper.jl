using Test

using WorhpOpt

using MathOptInterface
const MOI = MathOptInterface
const MOIT = MOI.Test
const MOIU = MOI.Utilities
const MOIB = MOI.Bridges

const optimizer = WorhpOpt.Optimizer()



const config_no_duals = MOIT.TestConfig(atol=1e-4, rtol=1e-4, duals=false,
                               optimal_status=MOI.LOCALLY_SOLVED)



@testset "SolverName" begin
    @test MOI.get(optimizer, MOI.SolverName()) == "Worhp"
end

@testset "supports_default_copy_to" begin
    @test MOIU.supports_default_copy_to(optimizer, false)
    # Use `@test !...` if names are not supported
    @test !MOIU.supports_default_copy_to(optimizer, true)
end

const config_no_duals = MOIT.TestConfig(atol=1e-4, rtol=1e-4, duals=false,
                                        optimal_status=MOI.LOCALLY_SOLVED)
const BRIDGED = MOIB.full_bridge_optimizer(WorhpOpt.Optimizer(), Float64)

@testset "Unit" begin
    exclude = ["delete_variable", # Deleting not supported.
               "delete_variables", # Deleting not supported.
               "getvariable", # Variable names not supported.
               "solve_zero_one_with_bounds_1", # Variable names not supported.
               "solve_zero_one_with_bounds_2", # Variable names not supported.
               "solve_zero_one_with_bounds_3", # Variable names not supported.
               "getconstraint", # Constraint names not suported.
               "variablenames", # Variable names not supported.
               "solve_with_upperbound", # loadfromstring!
               "solve_with_lowerbound", # loadfromstring!
               "solve_integer_edge_cases", # loadfromstring!
               "solve_affine_lessthan", # loadfromstring!
               "solve_affine_greaterthan", # loadfromstring!
               "solve_affine_equalto", # loadfromstring!
               "solve_affine_interval", # loadfromstring!
               "get_objective_function", # Function getters not supported.
               "solve_constant_obj",  # loadfromstring!
               "solve_blank_obj", # loadfromstring!
               "solve_singlevariable_obj", # loadfromstring!
               "solve_objbound_edge_cases", # ObjectiveBound not supported.
               "solve_affine_deletion_edge_cases", # Deleting not supported.
               "solve_unbounded_model", # `NORM_LIMIT`
               "number_threads", # NumberOfThreads not supported
               "delete_nonnegative_variables", # get ConstraintFunction n/a.
               "update_dimension_nonnegative_variables", # get ConstraintFunction n/a.
               "delete_soc_variables", # VectorOfVar. in SOC not supported
               "solve_result_index", # DualObjectiveValue not supported
               ]
    MOIT.unittest(
        BRIDGED,
        config_no_duals,
        exclude
        )
end

@testset "Modification" begin
    MOIT.modificationtest(BRIDGED, config_no_duals)
end

@testset "Continuous Linear" begin
    MOIT.contlineartest(BRIDGED, config_no_duals)
end

@testset "Continuous Conic" begin
    MOIT.contlineartest(BRIDGED, config_no_duals)
end

@testset "Integer Conic" begin
    MOIT.intconictest(BRIDGED, config_no_duals)
end


MOI.empty!(optimizer)

@testset "MOI QP/QCQP tests" begin
    qp_optimizer = MOIU.CachingOptimizer(MOIU.Model{Float64}(), optimizer)
    MOIT.qptest(qp_optimizer, config)
    exclude = []
    MOIT.qcptest(qp_optimizer, config_no_duals,
        #exclude
        )
end

MOI.empty!(optimizer)

@testset "MOI NLP tests" begin
    MOIT.nlptest(optimizer, config_no_duals)
end
