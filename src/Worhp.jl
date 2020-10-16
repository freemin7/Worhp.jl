module Worph


baremodule LibWorhp
    using CBinding: 𝐣𝐥

    const size_t = 𝐣𝐥.Csize_t
    𝐣𝐥.Base.include((𝐣𝐥.@__MODULE__), 𝐣𝐥.joinpath(𝐣𝐥.dirname(𝐣𝐥.@__DIR__), "deps", "bindings.jl"))  #
end

end # module
