*! version 1.0.2 M. Park 4Apr2025

program define roy_specification, rclass byable(recall)
	version 18
	
	syntax varlist(min = 1 numeric fv) [if] [in], method(string) [SELect(varname) ///
	vce(string) absoff adj(string) /// Wald option
	 EXClude(varlist) tol(real 999) /// LR option
	 SMIV(varname) cube(string) boot(integer 500) /// SM option
	 tn(integer 50) qn(integer 999) eta(real 0.000001) eps(real 0.000001) /// sm hidden options
	 kn(integer 999) bn(integer 999) alpha(real 0.05)]

// ------------------------------------------------------------------ sample ---
	 marksample touse
	 
// -------------------------------------------------------------------- wald ---
	if ("`method'" == "wald") {
		if ("`select'" == "") {
			di as err "option {bf:select()} required"
			error 198
		}
		
		gettoken Y X : varlist
		gettoken D   : select
		gettoken V   : vce
		gettoken ADJ : adj
		
		markout `touse' `Y' `X' `D'
	
		_rmcoll `X', expand
		
		mata: instance = roy_wald()
		mata: instance.reportResult()
		mata: instance.saveResult()
		mata: mata drop instance // prevent crash

		return scalar N         = N
		return scalar stat      = stat
		return scalar df        = df
		return scalar p         = p
		return scalar reject_10 = reject_10
		return scalar reject_5  = reject_5
		return scalar reject_1  = reject_1
		return scalar kx        = kx
	
		return matrix theta     = theta
		return matrix V         = V
	
		return local  vce       = "`V'"
		return local  adj       = "`ADJ'"
	}

// ---------------------------------------------------------------------- lr ---
	if ("`method'" == "lr") {
		if ("`select'" == "") {
			di as err "option {bf:select()} required"
			error 198
		}
		
		gettoken Y X : varlist
		gettoken D   : select
		gettoken V   : vce
		gettoken EXC : exclude
		
		markout `touse' `Y' `X' `D'
	
		_rmcoll `X', expand

		mata: instance = roy_lr()
		mata: instance.report()
		mata: instance.returns()
		mata: mata drop instance // prevent crash

		return scalar N         = N
		return scalar kx        = kx
		return scalar df        = df
	
		return scalar stat      = stat
		return scalar p         = p
		return scalar reject_10 = reject_10
		return scalar reject_5  = reject_5
		return scalar reject_1  = reject_1
		return scalar llk_u     = llk_u
		return scalar llk_r     = llk_r
	
		return matrix theta_u   = theta_u
		return matrix theta_r   = theta_r
		return matrix V_u       = V_u
		return matrix V_r       = V_r
	
		return local  vce       = "`V'"
	}

// ---------------------------------------------------------------------- sm ---
	if ("`method'" == "sm") {
		if ("`smiv'" == "") {
			di as err "option {bf:smiv()} required"
			error 198
		}
		
		gettoken Y : varlist
		gettoken Z : smiv
		gettoken G : cube
	
		markout `touse' `Y' `Z'
	
		if ("`G'" == "") local G "countable"
	
		mata: instance = roy_sm()
		mata: instance.returns()
		mata: instance.report()
		mata: mata drop instance // prevent crash

		return scalar N         = N
		return scalar Tn        = Tn
		return scalar cv_10     = cv_10
		return scalar cv_5      = cv_5
		return scalar cv_1      = cv_1
		return scalar reject_10 = reject_10
		return scalar reject_5  = reject_5
		return scalar reject_1  = reject_1
	}
end

// -------------------------------------------------------------------- mata ---
include roy_specification.mata, adopath
