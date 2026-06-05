########
# fit outcome models
########
# step 1: compute Eq: E(q|enroll at j), same format as "shares"
aa,bb,cc = getPosteriorQ(zeros(7))
Eq = Dict(:eqbm=>aa, )
Eshock = Dict(:eqbm=>bb, )
Ev = Dict(:eqbm=>cc, )
shares_eqbm = hcat( 
    getPrEnroll(zeros(J),1:J,CC,CC.ttt,CC.ttt)...
)
shares = OrderedDict(:eqbm => shares_eqbm,)

include_v = false 

# step 2a: UTA regressions
#covariates in auxiliary model
ww_UTA =  hcat( [[1.0; CC.zi[ii][7,J+1:end]; ttt[ii]*1.0; longhorn[ii]; century[ii]; ziRaw[ii,:classrank].^2; ziRaw[ii,:classrank].^3] for ii=1:II]...) |> transpose |> collect 
ww_TAMU = ww_UTA

function obj_outcome(θoutcome,
    coefs_data=coefs_UTA,
    zq = zq_UTA,
    ww = ww_UTA,
    wts = shares_eqbm[7,:],
    )
y = zq*θoutcome
ep = y .- ww*coefs_data
G_foc = ww'*Diagonal(wts)*ep 
sum( G_foc.^2) 
end
function obj_outcome_func(θoutcome,
    fun::F,
    coefs_data=coefs_UTA,
    zq = zq_UTA,
    ww = ww_UTA,
    wts = shares_eqbm[7,:],
    ) where {F}
y = zq*θoutcome
ep = fun.(y) .- ww*coefs_data
G_foc = ww'*Diagonal(wts)*ep 
sum( G_foc.^2) 
end

function runOutcomeModels(CC=CC,include_v=false)
    #covariates in "true" model
    if include_v
        nθoutcome = nz + 4
        zq_UTA = hcat( [ [1.0; CC.zi[ii][7,J+1:end]; longhorn[ii]; Eq[:eqbm][7,ii]; Ev[:eqbm][7,ii]] for ii=1:II]...) |> transpose |> collect
        zq_TAMU = hcat( [ [1.0; CC.zi[ii][6,J+1:end]; century[ii]; Eq[:eqbm][6,ii]; Ev[:eqbm][6,ii]] for ii=1:II]...) |> transpose |> collect
    else
        nθoutcome = nz + 3
        zq_UTA = hcat( [ [1.0; CC.zi[ii][7,J+1:end]; longhorn[ii]; Eq[:eqbm][7,ii]] for ii=1:II]...) |> transpose |> collect
        zq_TAMU = hcat( [ [1.0; CC.zi[ii][6,J+1:end]; century[ii]; Eq[:eqbm][6,ii]] for ii=1:II]...) |> transpose |> collect
    end

    res_GPA_UTA = Optim.optimize(θ->obj_outcome(θ,coefs_UTA,zq_UTA),zeros(nθoutcome),LBFGS(),autodiff=:forward)
    θGPA_UTA = res_GPA_UTA.minimizer
    gpa_hat_UTA = zq_UTA*θGPA_UTA

    res_persist_UTA = Optim.optimize(θ->obj_outcome(θ,coefs_UTA_persist,zq_UTA),zeros(nθoutcome),LBFGS(),autodiff=:forward)
    θpersist_UTA = res_persist_UTA.minimizer
    persist_hat_UTA = zq_UTA*θpersist_UTA

    res_stem_UTA = Optim.optimize(θ->obj_outcome(θ,coefs_UTA_stem,zq_UTA,ww_UTA,shares_eqbm[7,:]),zeros(nθoutcome),LBFGS(),autodiff=:forward)
    θstem_UTA = res_stem_UTA.minimizer
    stem_hat_UTA = zq_UTA*θstem_UTA

    res_GPA_TAMU = Optim.optimize(θ->obj_outcome(θ,coefs_TAMU,zq_TAMU,ww_TAMU,shares_eqbm[6,:]),zeros(nθoutcome),LBFGS(),autodiff=:forward)
    θGPA_TAMU = res_GPA_TAMU.minimizer
    gpa_hat_TAMU = zq_TAMU*θGPA_TAMU

    res_persist_TAMU = Optim.optimize(θ->obj_outcome(θ,coefs_TAMU_persist,zq_TAMU,ww_TAMU,shares_eqbm[6,:]),zeros(nθoutcome),LBFGS(),autodiff=:forward)
    θpersist_TAMU = res_persist_TAMU.minimizer
    persist_hat_TAMU = zq_TAMU*θpersist_TAMU

    res_stem_TAMU = Optim.optimize(θ->obj_outcome(θ,coefs_TAMU_stem,zq_TAMU,ww_TAMU,shares_eqbm[6,:]),zeros(nθoutcome),LBFGS(),autodiff=:forward)
    θstem_TAMU = res_stem_TAMU.minimizer
    stem_hat_TAMU = zq_TAMU*θstem_TAMU

    ff(x) = cdf(Normal(),x)
    res_persist_UTA_probit = Optim.optimize(θ->obj_outcome_func(θ,ff,coefs_UTA_persist,zq_UTA),zeros(nθoutcome),LBFGS(),autodiff=:forward)
    res_persist_TAMU_probit = Optim.optimize(θ->obj_outcome_func(θ,ff,coefs_TAMU_persist,zq_TAMU,ww_TAMU,shares_eqbm[6,:]),zeros(nθoutcome),LBFGS(),autodiff=:forward)
    res_stem_UTA_probit = Optim.optimize(θ->obj_outcome_func(θ,ff,coefs_UTA_stem,zq_UTA,ww_UTA,shares_eqbm[7,:]),zeros(nθoutcome),LBFGS(),autodiff=:forward)
    res_stem_TAMU_probit = Optim.optimize(θ->obj_outcome_func(θ,ff,coefs_TAMU_stem,zq_TAMU,ww_TAMU,shares_eqbm[6,:]),zeros(nθoutcome),LBFGS(),autodiff=:forward)

    dfoutcomes = DataFrame(
        gpa_uta = θGPA_UTA,
        persist_uta=θpersist_UTA,
        persist_uta_probit=res_persist_UTA_probit.minimizer,
        stem_uta = θstem_UTA,
        stem_uta_probit = res_stem_UTA_probit.minimizer,
        gpa_tamu = θGPA_TAMU,
        persist_tamu=θpersist_TAMU,
        persist_tamu_probit=res_persist_TAMU_probit.minimizer,
        stem_tamu = θstem_TAMU,
        stem_tamu_probit = res_stem_TAMU_probit.minimizer,
    )
    return zq_UTA, zq_TAMU, dfoutcomes
end

zq_UTA, zq_TAMU, dfoutcomes = runOutcomeModels(CC,false)
# WtWinv_Wt_UTA = inv(ww_UTA'*Diagonal(shares_eqbm[7,:])*ww_UTA)*ww_UTA'*Diagonal(shares_eqbm[7,:])
# function obj_outcome_alt(θoutcome,
#         zq = zq_UTA,
#         WtWinv_Wt = WtWinv_Wt_UTA,
#         coefs_data=coefs_UTA,
#         var_y_data = var(uta_rdsample.gpa)
#         )
#     y = zq*θoutcome
#     coefs_model = WtWinv_Wt*y
#     sum( (coefs_model .- coefs_data).^2) + (var(y) - var_y_data)^2
# end
