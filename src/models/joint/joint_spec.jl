function joint_spec_sha256(path::String)
    isfile(path) || throw(ArgumentError("joint variable specification does not exist: $path"))
    return bytes2hex(open(sha256, path))
end

function load_joint_covariate_spec(path::String, products::Vector{String})
    isfile(path) || throw(ArgumentError("joint variable specification does not exist: $path"))
    table = CSV.read(path, DataFrame)
    required = [
        :product, :family, :variable_group, :predictor_columns,
        :independent_selected, :joint_vif_retained, :role, :final_included,
    ]
    missing_columns = setdiff(required, propertynames(table))
    isempty(missing_columns) || throw(ArgumentError(
        "joint variable specification is missing columns: $(join(string.(missing_columns), ", "))",
    ))
    allunique(zip(String.(table.product), String.(table.variable_group))) ||
        throw(ArgumentError("joint variable specification contains duplicate product/group rows"))
    Set(String.(table.product)) == Set(products) ||
        throw(ArgumentError("joint variable specification products differ from benchmark products"))
    all(in(Set(JOINT_GROUP_ORDER)), String.(table.variable_group)) ||
        throw(ArgumentError("joint variable specification contains an unsupported variable group"))
    included = filter(:final_included => identity, table)
    for row in eachrow(included)
        Bool(row.independent_selected) && Bool(row.joint_vif_retained) ||
            throw(ArgumentError("included variable bypassed independent screening or joint VIF"))
        String(row.role) in ("local", "global") ||
            throw(ArgumentError("included variable has unresolved role"))
        expected_family = GROUP_FAMILY[String(row.variable_group)]
        String(row.family) == expected_family ||
            throw(ArgumentError("variable family differs for $(row.variable_group)"))
    end
    role_maps = Dict{String,Dict{String,String}}()
    for product in products
        rows = filter([:product, :final_included] =>
            (p, included) -> String(p) == product && included, table)
        role_maps[product] = Dict(String(row.variable_group) => String(row.role)
            for row in eachrow(rows))
    end
    return (; table, included, role_maps)
end
