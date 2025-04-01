{smcl}
{* *! version 1.0.0 M. Park 28Mar2025}
{viewerjumpto "Wald_options" "roy_specification##Wald_options"}{...}
{viewerjumpto "LR_options" "roy_specification##LR_options"}{...}
{viewerjumpto "SM_options" "roy_specification##SM_options"}{...}
{viewerjumpto "Roy_reference" "roy_specification##Roy_reference"}{...}
{title:Title}

{phang}
{bf:roy_specification} {hline 1} Specification tests for Roy models

{title:Syntax}

{pstd}
Wald test
{p_end}

{phang2}
{cmd:roy_specification} {depvar} {indepvars} {ifin}, {opt method}({it:wald}) {opth sel:ect(varname)} [{cmd:}{help roy_specification##Wald_options:Wald_options}]
{p_end}

{pstd}
Likelihood ratio test
{p_end}

{phang2}
{cmd:roy_specification} {depvar} {indepvars} {ifin}, {opt method}({it:lr}) {opth sel:ect(varname)} [{cmd:}{help roy_specification##LR_options:LR_options}]
{pstd}
Stochastic monotonicity test
{p_end}

{phang2}
{cmd:roy_specification} {depvar} {ifin}, {opt method}({it:sm}) [{cmd:}{help roy_specification##SM_options:SM_options}]{synoptset 15 tabbed}{...}
{synopt:{space 2}{it:depvar}} outcome variable{p_end}
{synopt:{space 2}{it:indepvars}} covariates{p_end}
{synopt:{space 2}{it:select}} selection variable{p_end}

{title:Description}

{pstd}
{bf:roy_specification} implements specification tests for Roy models. The null model is the original Roy selection rule. Non-rejection of the null model indicates that the outcome variable is the only determinant of agents' choice. Three methods, Wald, LR, and SM are supported, and technical details are available in {help roy_specification##Roy_reference:Reference}.


{title:Options}

{marker Wald_options}{...}
{dlgtab: Wald_options}

{pmore}
{opt vce}({it:vcetype}) specifies the type of the covariance matrix used in the Wald statistic. The standard covariance matrix ({it:mle}), which is the default type, and the sandwich ({it:robust)} covariance matrix are available.
{p_end}

{pmore}
{opt absoff} specifies that the test statistic does not use the absolute values. If you do not use the absolute value, it is possible to not reject the null even though the selection rule is outcome minimization behavior.
{p_end}
{pmore}
{opt adj}({it:string}) makes a small-sample adjustment. The default is to not use the adjustment ({it:none}), and the adjustment for degrees of freedom ({it:df}) is available. The adjustment multiplies the test statistic by (K - k) / N, where k is the number of covariates including the intercept term.
{p_end}

{marker LR_options}{...}
{dlgtab: LR_options}

{pmore}
{opt exclude}({it:varlist}) specifies the exclusion restrictions. For the {it:varlist}, the maximum likelihood estimations assume that the coefficients in the selection equation are zeros. If the exclusion restrictions are misspecified, the MLE may not converge or the estimates are extreme.
{p_end}

{pmore}
{opt vce}({it:vcetype}) specifies the type of the covariance matrix used in the Wald statistic. The standard covariance matrix ({it:mle}), which is the default type, and the sandwich ({it:robust)} covariance matrix are available.
{p_end}
{pmore}
{opt tol}({it:real}) sets the tolerance parameters. See help mata moptimize_init_conv_ptol, moptimize_init_conv_vtol, and moptimize_init_conv_nrtol.
{p_end}

{marker SM_options}{...}
{dlgtab: SM_options}

{phang2}
{opt smiv}({it:varname}) sets the stochastically monotone instrumental variable (SMIV). Currently, only one SMIV is supported.
{p_end}

{phang2}
{opt cube}({it:string}) specifies the hypercubes in the generalized regression monotonicity test, and "countable" and "continuum" are available. For more details, see {help roy_specification##Roy_reference:Hsu et al 2019}.
{p_end}

{phang2}
{opt boot}({it:string}) sets the number of bootstrap for the critical value. Default is 500.
{p_end}

{pmore}
{bf:Note: }other user-chosen parameters of the generalized regression monotonicity test and the generalized moment selection method are set following the replication package of {help roy_specification##Roy_reference:Mourifie et al 2020} and suggested values. The parameters can be specified using {opt qn}({it:int}), {opt kn}({it:int}), and {opt bn}({it:int}). Users interested in these options can refer to the hidden options in the ado and mata source codes.
{p_end}

{title:Stored results}

{dlgtab: Wald}

{pmore}
r(kx): Number of covariates.
{p_end}

{pmore}
r(reject_*): Indicators for the rejection at the nominal level * / 100.
{p_end}

{pmore}
r(p): p-value.
{p_end}

{pmore}
r(df): Degrees of freedom of the null distribution. It is k if r(kx) > 1 and N - 7 if r(kx) = 1.
{p_end}

{pmore}
r(stat): Value of the test statistic.
{p_end}

{pmore}
r(N): Number of observations (used in the test).
{p_end}

{pmore}
r(adj): Type of the small-sample adjustment.
{p_end}

{pmore}
r(vce): Type of the covariance matrix.
{p_end}

{pmore}
r(V): Covariance matrix of the MLE.
{p_end}

{pmore}
r(theta): Maximum likelihood estimates.
{p_end}

{dlgtab: LR}

{pmore}
r(llr_*): Values of maximum likelihood functions. For *, u and r stands for unrestricted MLE and restricted MLE, respectively.
{p_end}

{pmore}
r(reject_*): Indicators for the rejection at the nominal level * / 100.
{p_end}

{pmore}
r(p): p-value.
{p_end}

{pmore}
r(stat): Value of the test statistic.
{p_end}

{pmore}
r(df): Degrees of freedom of the null distribution. It is k if r(k) > 1 and N - 7 if r(k) = 1.
{p_end}

{pmore}
r(kx): Number of covariates.
{p_end}

{pmore}
r(N): Number of observations (used in the test).
{p_end}

{pmore}
r(vce): Type of the covariance matrix.
{p_end}

{pmore}
r(V_*): Covariance matrix of the MLE. For *, u and r stands for unrestricted MLE and restricted MLE, respectively.
{p_end}

{pmore}
r(theta_*): Maximum likelihood estimates. For *, u and r stands for unrestricted MLE and restricted MLE, respectively.
{p_end}

{dlgtab: SM}

{pmore}
r(reject_*): Indicators for the rejection at the nominal level * / 100.
{p_end}

{pmore}
r(cv_*): Generalized moment selection critical values defined as the (1 - * / 100 + eta)-th quantiles of the null distribution plus eta. The default value of eta is 10^-16, and it can be adjusted using a hidden option {opt eta}({it:real}). Reject the null hypothesis at the nominal level * / 100, if the value of the test statistic is larger than cv_*.
{p_end}

{pmore}
r(Tn): Value of the test statistic.
{p_end}

{pmore}
r(N): Number of observations (used in the test).
{p_end}

{title:Examples}

{phang}
Wald test
{p_end}

	{cmd:webuse nlswork, clear}
	{cmd:gen agesq = age^2}
	
	1) Basic usage

	{cmd:roy_specification ln_wage age agesq ib3.race collgrad ttl_exp hours, method(wald) sel(msp)}
	{cmd:roy_specification ln_wage age agesq ib3.race collgrad ttl_exp hours, method(wald) sel(msp) vce(mle)}

	2) Robust covariance matrix

	{cmd:roy_specification ln_wage age agesq ib3.race collgrad ttl_exp hours, method(wald) sel(msp) vce(robust)}

	3) Adjustment for degrees of freedom

	{cmd:roy_specification ln_wage age agesq ib3.race collgrad ttl_exp hours, method(wald) sel(msp) adj(df)}

	4) No absolute value
	
	{cmd:roy_specification ln_wage age agesq ib3.race collgrad ttl_exp hours, method(wald) sel(msp) absoff}
	
        5) if and in

	{cmd:roy_specification ln_wage age agesq ib3.race collgrad ttl_exp hours if year == 68, method(wald) sel(msp)}
	{cmd:roy_specification ln_wage age agesq ib3.race collgrad ttl_exp hours in 1/20000, method(wald) sel(msp)}
	6) By group

	{cmd:keep if year == 68 | year == 69 | year == 70}
	{cmd:bysort year: roy_specification ln_wage age agesq ib3.race collgrad ttl_exp hours, method(wald) sel(msp)}

{phang}
Likelihood ratio test (note that the LR test can fail to converge if the identification assumptions are not satisfied)
{p_end}

	{cmd:clear all}
	{cmd:set seed 1}
	{cmd:quietly set obs 3000}

	{cmd:matrix muu = (0, 0)}
	{cmd:matrix sigmau = (1, 0.5 \ 0.5, 1)}
	{cmd:drawnorm u0 u1, mean(muu) cov(sigmau)}

	{cmd:gen X1 = runiform(0, 1)}
	{cmd:gen X2 = rnormal(0, 1)}
	{cmd:gen X3 = rnormal(0, 1)}

	{cmd:gen Y0 = 1 + X2 + X3 + u0}
	{cmd:gen Y1 = -1 + X1 - X2 - X3 + u1}
	{cmd:gen V  = -X1 // X2 and X3 are excluded} // X2 and X3 are excluded

	{cmd:gen D = Y1 - Y0 + V >= 0}
	{cmd:gen Y = D * Y1 + (1 - D) * Y0}

	1) Basic usage
	{cmd:roy_specification Y X1 X2 X3, method(lr) sel(D) exc(X2 X3)}
	{cmd:roy_specification Y X1 X2 X3, method(lr) sel(D) exc(X2)}

	2) Robust covariance matrix

	{cmd:roy_specification Y X1 X2 X3, method(lr) sel(D) exc(X2 X3) vce(robust)}	3) Wrong exclusion assumption
	{cmd:roy_specification Y X1 X2 X3, method(lr) sel(D) exc(X1)}

{phang}
Stochastic monotonicity test
{p_end}

	{cmd:clear all}
	{cmd:quietly set obs 500}
	{cmd:set seed 1}

	{cmd:gen U = rnormal(0, 0.1)}
	{cmd:gen Z = runiform(0, 1)}
	{cmd:gen Y = .}
	{cmd:gen G = rbinomial(1, 0.5)}
	{cmd:replace Y = U if G == 0} // i.e. Y is stochastically monotone w.r.t. Z if G == 0
	{cmd:replace Y = Z * (1 - Z) + U if G == 1} // i.e. Y is not stochastically monotone w.r.t. Z if G == 1
	
	1) Basic usage

	{cmd:roy_specification Y if G == 0, method(sm) smiv(Z)}
	{cmd:roy_specification Y if G == 1, method(sm) smiv(Z)}

	2) Continuum hypercubes

	{cmd:roy_specification Y if G == 0, method(sm) smiv(Z) cube(continuum)}
	{cmd:roy_specification Y if G == 1, method(sm) smiv(Z) cube(continuum)}
	3) Number of bootstrap

	{cmd:roy_specification Y if G == 0, method(sm) smiv(Z) boot(200)}
	{cmd:roy_specification Y if G == 1, method(sm) smiv(Z) boot(200)}	4) By group

	{cmd:bysort G: roy_specification Y, method(sm) smiv(Z) cube(continuum) boot(200)}

{marker Roy_reference}{...}
{title:Reference}

{pstd}
Björklund, A., & Moffitt, R. (1987). The estimation of wage gains and welfare gains in self-selection models. The Review of Economics and Statistics, 42-49.
{p_end}

{pstd}
Han, C., & Park, M. (2025). A simple specification test for Roy models. Working paper.
{p_end}

{pstd}
Hsu, Y. C., Liu, C. A., & Shi, X. (2019). Testing generalized regression monotonicity. Econometric Theory, 35(6), 1146-1200.
{p_end}

{pstd}
Mourifie, I., Henry, M., & Meango, R. (2020). Sharp bounds and testability of a Roy model of STEM major choices. Journal of Political Economy, 128(8), 3220-3283.
{p_end}

{title:Author}
{phang}
Minchul Park (Korea University Economics), mcpark1352@gmail.com
{p_end}