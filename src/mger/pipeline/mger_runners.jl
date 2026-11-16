function run_spatial_kfold_pipeline(
	cfg::MGERConfig; k::Int=5, rng::AbstractRNG=Random.GLOBAL_RNG, fold_scheme::Symbol=:random,
)
	mkpath(cfg.outdir)

	st = load_station_meta(cfg.station_meta_path;
		station_id_col=cfg.station_id_col, lon_col=cfg.lon_col, lat_col=cfg.lat_col)
	products, common_product_ids, product_data = load_global_common_product_data(cfg)

	common_lonlat = build_X_lonlat(st, common_product_ids; station_id_col=cfg.station_id_col, lon_col=cfg.lon_col, lat_col=cfg.lat_col)
	folds = if fold_scheme == :random
		split_stations_kfold(common_product_ids; k=k, rng=rng)
	elseif fold_scheme == :spatial_block
		split_stations_spatial_block_kfold(common_product_ids, common_lonlat; k=k)
	else
		throw(ArgumentError("fold_scheme must be :random or :spatial_block"))
	end
	cv_scheme = fold_scheme == :spatial_block ? "station_spatial_block_$(k)fold" : "station_$(k)fold"
	write_validation_scope(
		cfg.outdir;
		cv_scheme="$(cv_scheme)_pooled",
		use_loocv_eval=cfg.use_loocv_eval,
		fold_scheme,
		slope_ridge_candidates=cfg.slope_ridge_candidates,
		min_scan_coverage=cfg.min_scan_coverage,
		kernel_candidates=cfg.kernels,
	)

	common_split = DataFrame(station_id=String[], fold=Int[], fold_scheme=String[])
	for (fold_idx, fold_ids) in enumerate(folds)
		for station_id in fold_ids
			push!(common_split, (station_id=station_id, fold=fold_idx, fold_scheme=string(fold_scheme)))
		end
	end
	CSV.write(joinpath(cfg.outdir, "split_common_spatial5fold.csv"), common_split)

	fold_summary_rows = NamedTuple[]
	pooled_summary_rows = NamedTuple[]
	for product in products
		data = product_data[product]
		times = data.times
		id_map_all = Dict(id => i for (i, id) in enumerate(data.ids))
		use_idx = [id_map_all[id] for id in common_product_ids]
		ids = common_product_ids
		Y_obs = data.Y_obs[use_idx, :]
		Y_sat = data.Y_sat[use_idx, :]
		lonlat = build_X_lonlat(st, ids; station_id_col=cfg.station_id_col, lon_col=cfg.lon_col, lat_col=cfg.lat_col)
		Y_corr_pooled = fill(NaN, size(Y_obs))
		selected_slope_ridges = Float64[]

		product_split = DataFrame(station_id=String[], fold=Int[], fold_scheme=String[])
		for (fold_idx, fold_ids) in enumerate(folds)
			for station_id in fold_ids
				push!(product_split, (station_id=station_id, fold=fold_idx, fold_scheme=string(fold_scheme)))
			end
		end
		CSV.write(joinpath(cfg.outdir, "split_$(product)_spatial5fold.csv"), product_split)

		for fold_idx in 1:k
			val_ids = folds[fold_idx]
			train_ids = String[]
			for other_fold in 1:k
				other_fold == fold_idx && continue
				append!(train_ids, folds[other_fold])
			end

			fold_dir = joinpath(cfg.outdir, "fold_$(fold_idx)")
			mkpath(fold_dir)
			split_df = DataFrame(
				station_id = vcat(train_ids, val_ids),
				split = vcat(fill("train", length(train_ids)), fill("validation", length(val_ids))),
				fold = fill(fold_idx, length(train_ids) + length(val_ids)),
				fold_scheme = fill(string(fold_scheme), length(train_ids) + length(val_ids)),
			)
			CSV.write(joinpath(fold_dir, "split_$(product)_spatialcv.csv"), split_df)

			println("[$product] $cv_scheme fold $fold_idx/$k: $(length(train_ids)) train sites, $(length(val_ids)) val sites")
			result = evaluate_spatial_holdout(
				product, ids, times, Y_obs, Y_sat, lonlat, train_ids, val_ids, cfg;
				scan_path=joinpath(fold_dir, "scan_$(product)_spatialcv.csv"),
				fail_path=joinpath(fold_dir, "scan_$(product)_spatialcv_failures.csv"),
			)
			Y_corr_pooled[result.val_idx, :] = result.Y_corr_val
			push!(selected_slope_ridges, result.slope_ridge)
			write_wide(joinpath(fold_dir, "corr_$(product)_spatialcv_val.csv"), times, val_ids, result.Y_corr_val)

			push!(fold_summary_rows, (
				product=product,
				fold=fold_idx,
				spatial_cv=true,
				cv_scheme=cv_scheme,
				n_train=result.n_train,
				n_val=result.n_val,
				n=result.n,
				pre_n=result.pre_n,
				post_n=result.post_n,
				common_n=result.common_n,
				total_n=result.total_n,
				coverage=result.coverage,
				scan_n=result.scan_n,
				scan_coverage=result.scan_coverage,
				kernel=result.kernel,
				adaptive=result.adaptive,
				bw=result.bw,
				slope_ridge=result.slope_ridge,
				RMSE_pre=result.RMSE_pre, RMSE_post=result.RMSE_post,
				MAE_pre=result.MAE_pre, MAE_post=result.MAE_post,
				Bias_pre=result.Bias_pre, Bias_post=result.Bias_post,
				r_pre=result.r_pre, r_post=result.r_post,
				POD_pre=result.POD_pre, POD_post=result.POD_post,
				FAR_pre=result.FAR_pre, FAR_post=result.FAR_post,
				CSI_pre=result.CSI_pre, CSI_post=result.CSI_post,
				negative_n=result.negative_n,
				negative_fraction=result.negative_fraction,
				min_corrected=result.min_corrected,
			))
		end

		write_wide(joinpath(cfg.outdir, "corr_$(product)_spatial5fold_val.csv"), times, ids, Y_corr_pooled)
		pre_mask = common_valid_mask(Y_obs, Y_sat)
		post_mask = common_valid_mask(Y_obs, Y_corr_pooled)
		eval_mask = common_valid_mask(Y_obs, Y_sat, Y_corr_pooled)
		@assert any(eval_mask) "[$product] five-fold validation has no common valid samples; cannot compare before/after correction"

		pre_c = metric_continuous(Y_obs, Y_sat; mask=eval_mask)
		pre_e = metric_event(Y_obs, Y_sat; thr=cfg.rain_threshold, mask=eval_mask)
		post_c = metric_continuous(Y_obs, Y_corr_pooled; mask=eval_mask)
		post_e = metric_event(Y_obs, Y_corr_pooled; thr=cfg.rain_threshold, mask=eval_mask)
		negative = negative_output_stats(Y_corr_pooled)
		push!(pooled_summary_rows, (
			product=product,
			spatial_cv=true,
			cv_scheme="$(cv_scheme)_pooled",
			k=k,
			n_station=length(ids),
			n=post_c.n,
			pre_n=count(pre_mask),
			post_n=count(post_mask),
			common_n=count(eval_mask),
			total_n=length(eval_mask),
			coverage=count(eval_mask) / length(eval_mask),
			slope_ridge_by_fold=join(selected_slope_ridges, ","),
			RMSE_pre=pre_c.RMSE, RMSE_post=post_c.RMSE,
			MAE_pre=pre_c.MAE, MAE_post=post_c.MAE,
			Bias_pre=pre_c.Bias, Bias_post=post_c.Bias,
			r_pre=pre_c.r, r_post=post_c.r,
			POD_pre=pre_e.POD, POD_post=post_e.POD,
			FAR_pre=pre_e.FAR, FAR_post=post_e.FAR,
			CSI_pre=pre_e.CSI, CSI_post=post_e.CSI,
			negative...,
		))
	end

	fold_summary_df = DataFrame(fold_summary_rows)
	pooled_summary_df = DataFrame(pooled_summary_rows)
	fold_stats_df = summarize_fold_metrics(fold_summary_df)
	CSV.write(joinpath(cfg.outdir, "summary_three_products_folds.csv"), fold_summary_df)
	CSV.write(joinpath(cfg.outdir, "summary_three_products_fold_stats.csv"), fold_stats_df)
	CSV.write(joinpath(cfg.outdir, "summary_three_products_pooled.csv"), pooled_summary_df)
	CSV.write(joinpath(cfg.outdir, "summary_three_products.csv"), pooled_summary_df)
	return pooled_summary_df
end


function add_kernel_metadata!(df::DataFrame, kernel::Int, name::AbstractString)
	if :kernel in propertynames(df)
		df[!, :kernel] .= kernel
	else
		insertcols!(df, 2, :kernel => fill(kernel, nrow(df)))
	end
	if :kernel_name in propertynames(df)
		df[!, :kernel_name] .= String(name)
	else
		insertcols!(df, 3, :kernel_name => fill(String(name), nrow(df)))
	end
	return df
end


"""
    run_multikernel_spatial_kfold_pipeline(
        cfg::MGERConfig; k::Int=5, seed::Integer, fold_scheme::Symbol=:random,
    )

Run one independent spatial K-fold experiment per configured kernel. Each kernel
uses the same station split and selects its bandwidth from the configured grids.
"""
function run_multikernel_spatial_kfold_pipeline(
	cfg::MGERConfig; k::Int=5, seed::Integer, fold_scheme::Symbol=:random,
)
	isempty(cfg.kernels) && throw(ArgumentError("cfg.kernels must not be empty"))
	length(unique(cfg.kernels)) == length(cfg.kernels) ||
		throw(ArgumentError("cfg.kernels contains duplicate kernel values"))
	isempty(cfg.bw_adaptive) && isempty(cfg.bw_fixed_km) &&
		throw(ArgumentError("at least one adaptive or fixed bandwidth is required"))
	isempty(cfg.slope_ridge_candidates) &&
		throw(ArgumentError("cfg.slope_ridge_candidates must not be empty"))
	all(isfinite(ridge) && ridge >= 0 for ridge in cfg.slope_ridge_candidates) ||
		throw(ArgumentError("cfg.slope_ridge_candidates must be finite and non-negative"))
	isfinite(cfg.min_scan_coverage) && 0.0 <= cfg.min_scan_coverage <= 1.0 ||
		throw(ArgumentError("cfg.min_scan_coverage must be between 0 and 1"))
	for kernel in cfg.kernels
		kernel_name(kernel) # Validate before creating any output.
	end

	mkpath(cfg.outdir)
	status_df = DataFrame(
		kernel=Int[], kernel_name=String[], status=String[],
		error=String[], outdir=String[],
	)
	pooled_tables = DataFrame[]
	fold_tables = DataFrame[]
	fold_stats_tables = DataFrame[]

	for kernel in cfg.kernels
		name = kernel_name(kernel)
		kernel_outdir = joinpath(cfg.outdir, name)
		kernel_cfg = config_for_kernel(cfg, kernel, kernel_outdir)
		try
			pooled = run_spatial_kfold_pipeline(
				kernel_cfg;
				k=k,
				rng=MersenneTwister(seed),
				fold_scheme=fold_scheme,
			)
			add_kernel_metadata!(pooled, kernel, name)
			folds = CSV.read(joinpath(kernel_outdir, "summary_three_products_folds.csv"), DataFrame)
			fold_stats = CSV.read(joinpath(kernel_outdir, "summary_three_products_fold_stats.csv"), DataFrame)
			add_kernel_metadata!(folds, kernel, name)
			add_kernel_metadata!(fold_stats, kernel, name)
			push!(pooled_tables, pooled)
			push!(fold_tables, folds)
			push!(fold_stats_tables, fold_stats)
			push!(status_df, (kernel, name, "success", "", kernel_outdir))
		catch e
			message = sprint(showerror, e)
			push!(status_df, (kernel, name, "failed", message, kernel_outdir))
			@error "Kernel experiment failed" kernel name exception=(e, catch_backtrace())
		end
	end

	CSV.write(joinpath(cfg.outdir, "kernel_run_status.csv"), status_df)
	if !isempty(pooled_tables)
		CSV.write(
			joinpath(cfg.outdir, "summary_five_kernels_pooled.csv"),
			vcat(pooled_tables...; cols=:union),
		)
		CSV.write(
			joinpath(cfg.outdir, "summary_five_kernels_folds.csv"),
			vcat(fold_tables...; cols=:union),
		)
		CSV.write(
			joinpath(cfg.outdir, "summary_five_kernels_fold_stats.csv"),
			vcat(fold_stats_tables...; cols=:union),
		)
	end

	failed = status_df[status_df.status .== "failed", :]
	if nrow(failed) > 0
		error("$(nrow(failed)) kernel experiment(s) failed; see $(joinpath(cfg.outdir, "kernel_run_status.csv"))")
	end
	return vcat(pooled_tables...; cols=:union)
end

"""
Per-product count of how many of the k folds picked each kernel, from a
`summary_three_products_folds.csv`-shaped table (must carry `product`, `fold`, `kernel` columns).
A kernel winning every fold means the choice is stable; a split vote means it is fold-dependent.
"""
function _kernel_selection_stability(folds::DataFrame)
	nrow(folds) == 0 && return DataFrame()
	rows = NamedTuple[]
	for group in groupby(folds, [:product, :kernel])
		push!(rows, (;
			product=String(group.product[1]), kernel=Int(group.kernel[1]),
			kernel_name=kernel_name(Int(group.kernel[1])), win_count=nrow(group),
		))
	end
	stability = DataFrame(rows)
	fold_counts = combine(groupby(stability, :product), :win_count => sum => :fold_count)
	stability = leftjoin(stability, fold_counts; on=:product)
	stability.selection_frequency = stability.win_count ./ stability.fold_count
	return sort!(stability, [:product, order(:win_count; rev=true)])
end

"""
    run_nested_kernel_spatial_kfold_pipeline(
        cfg::MGERConfig; k::Int=5, seed::Integer, fold_scheme::Symbol=:random,
    )

Select the kernel via nested cross-validation: every training fold picks its own best
`(kernel, bw, adaptive, slope_ridge)` by that fold's own LOOCV, using only training-fold data,
and the held-out fold stations validate the winning combination — reusing `scan_params`'s
existing joint scan over `cfg.kernels` unrestricted, rather than forcing one kernel per run.

Contrast with `run_multikernel_spatial_kfold_pipeline`, which forces one kernel per run (via
`config_for_kernel`) to give each kernel its own paired head-to-head comparison on the same
fold split; its pooled output answers "how does each kernel perform," not "which kernel should
be used," and reading it as a selection recommendation is exactly the leak this function closes.
"""
function run_nested_kernel_spatial_kfold_pipeline(
	cfg::MGERConfig; k::Int=5, seed::Integer, fold_scheme::Symbol=:random,
)
	length(cfg.kernels) > 1 || throw(ArgumentError(
		"run_nested_kernel_spatial_kfold_pipeline needs more than one kernel in cfg.kernels; " *
		"for a single fixed kernel use run_spatial_kfold_pipeline directly",
	))
	pooled = run_spatial_kfold_pipeline(cfg; k, rng=MersenneTwister(seed), fold_scheme)
	folds = CSV.read(joinpath(cfg.outdir, "summary_three_products_folds.csv"), DataFrame)
	stability = _kernel_selection_stability(folds)
	CSV.write(joinpath(cfg.outdir, "kernel_selection_stability.csv"), stability)
	return pooled
end
