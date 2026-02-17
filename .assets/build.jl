using Literate

# Config for the execution of notebooks
cfg = Dict(
    "credit" => false,
    "execute" => true,
    "continue_on_error" => true,
)

@info "Building to $(pwd())"
Literate.markdown("travail.jl", pwd(); config=cfg, flavor=Literate.CommonMarkFlavor())
