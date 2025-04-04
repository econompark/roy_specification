*! version 1.0.2 M. Park 4Apr2025

version 18
set matastrict on

// TODO - use better names

// -------------------------------------------------------------------- wald ---
mata:
void llk(transmorphic M, real rowvector b, real colvector lnf) {
	real colvector p1, p2, p3, p4, p5, p6, p7, y, d

	p1 = moptimize_util_xb(M, b, 1) // gamma
	p2 = moptimize_util_xb(M, b, 2) // beta0
	p3 = moptimize_util_xb(M, b, 3) // beta1
	p4 = moptimize_util_xb(M, b, 4) // sig0
	p5 = moptimize_util_xb(M, b, 5) // sig1
	p6 = moptimize_util_xb(M, b, 6) // rho0
	p7 = moptimize_util_xb(M, b, 7) // rho1

	y  = moptimize_util_depvar(M, 1)
	d  = moptimize_util_depvar(M, 2)

	lnf = (1 :- d) :* (-ln(p4) :- .5 :* ((y :- p2) :/ p4):^2 :+ ///
	lnnormal(-(p1 :+ (p6 :/ p4) :* (y :- p2)) :/ sqrt(1 :- p6:^2))) :+ ///
	d :* (-ln(p5) :- .5 :* ((y :- p3) :/ p5):^2 :+ ///
	lnnormal((p1 :+ (p7 :/ p5) :* (y :- p3)) :/ sqrt(1 :- p7:^2)))
}

class roy_wald {
	public:
	// initialize
	void new()
	
	// setup
	real scalar    N, kx, df
	real rowvector alp
	real matrix    X
	
	// estimates
	real scalar    sig0, sig1, rho0, rho1, sig
	real colvector gam, bet0, bet1
	real matrix    Om
	transmorphic   M
	void           mleu()
	
	// wald statistic
	real scalar    computeTestStatistic()
	real scalar    computePValue()
	real rowvector checkReject()
	real colvector computeCoefs()
	real colvector h()
	real matrix    Dh()

	// report, returns
	string matrix  saveCoefsName()
	void           reportResult()
	void           saveResult()
}

void roy_wald::new() {
	N = rows(st_data(., st_local("Y"), st_local("touse")))
	if (st_local("X") == "") {
		X  = J(N, 1, 1)
		kx = 1
	} else if (st_local("X") != "") {
		X  = st_data(., st_local("X"), st_local("touse"))
		X  = select(X, colsum(X :== 0) :!= N)
		kx = cols(X) + 1
	}
	if (kx == 1) {
		df = N - 7
	} else if (kx > 1) {
		df = kx
	}
	
	alp = (.1, .05, .01)
	M   = moptimize_init()
	
	mleu()

	gam  = moptimize_result_coefs(M)[1 .. kx]'
	bet0 = moptimize_result_coefs(M)[(kx + 1) .. (2 * kx)]'
	bet1 = moptimize_result_coefs(M)[(2 * kx + 1) .. (3 * kx)]'
	
	sig0 = moptimize_result_coefs(M)[(3 * kx + 1)]
	sig1 = moptimize_result_coefs(M)[(3 * kx + 2)]
	rho0 = moptimize_result_coefs(M)[(3 * kx + 3)]
	rho1 = moptimize_result_coefs(M)[(3 * kx + 4)]
	
	if (st_local("absoff") == "") {
		sig = abs(sig1 * rho1 - sig0 * rho0)
	} else if (st_local("absoff") != "") {
		sig = sig1 * rho1 - sig0 * rho0
	}
	
	if (st_local("V") == "" | st_local("V") == "mle") {
		Om = moptimize_result_V(M)
	} else if (st_local("V") == "robust") {
		Om = moptimize_result_V_robust(M)
	}
}

void roy_wald::mleu() {
	moptimize_init_evaluator(M, &llk())
	
	moptimize_init_depvar(M, 1, st_data(., st_local("Y"), st_local("touse")))
	moptimize_init_depvar(M, 2, st_data(., st_local("D"), st_local("touse")))

	if (kx == 1) {
		moptimize_init_eq_indepvars(M, 1, "")
		moptimize_init_eq_indepvars(M, 2, "")
		moptimize_init_eq_indepvars(M, 3, "")
	} else if (kx > 1) {
		moptimize_init_eq_indepvars(M, 1, X)
		moptimize_init_eq_indepvars(M, 2, X)
		moptimize_init_eq_indepvars(M, 3, X)
	}
	
	moptimize_init_eq_indepvars(M, 4, "")
	moptimize_init_eq_indepvars(M, 5, "")
	moptimize_init_eq_indepvars(M, 6, "")
	moptimize_init_eq_indepvars(M, 7, "")
	
	moptimize_init_tracelevel(M, "none")
	moptimize_init_valueid(M, "Log likelihood")
	
	moptimize_init_search(M, "on")
	moptimize_init_search_random(M, "off")
	
	moptimize(M)
}

real colvector roy_wald::computeCoefs() {
	return((gam \ bet0 \ bet1 \ sig0 \ sig1 \ rho0 \ rho1))
}

real colvector roy_wald::h() {
	return(gam :* sig :- (bet1 :- bet0))
}

real matrix roy_wald::Dh() {
	return((I(kx) :* sig, I(kx), -I(kx), ///
	-rho0 :* gam, rho1 :* gam, -sig0 :* gam, sig1 :* gam))
}

real scalar roy_wald::computeTestStatistic() {
	if (kx == 1) {
		// the Wald t-test statistic
		return(h() / sqrt(Dh() * Om * Dh()'))
	}
	if (st_local("ADJ") == "" | st_local("ADJ") == "none") {
		return((h()' * invsym(Dh() * Om * Dh()') * h()))
	} else if (st_local("ADJ") == "df") {
		return(((N - length(computeCoefs())) / N) * (h()' * invsym(Dh() * Om * Dh()') * h()))
	}
}

real scalar roy_wald::computePValue() {
	if (kx == 1) {
		return(2 * (1 - t(df, abs(computeTestStatistic()))))
	} else if (kx > 1) {
		return(1 - chi2(df, computeTestStatistic()))
	}
}

real rowvector roy_wald::checkReject() {
	return(computePValue() :< alp)
}

string matrix roy_wald::saveCoefsName() {
	string colvector colnames, eqs

	if (kx == 1) {
		colnames = ("gam" \ "bet0" \ "bet1" \ "sig0" \ "sig1" \ "rho0" \ "rho1")
	} else if (kx > 1) {
		colnames = colshape(J(1, kx, ("gam" \ "bet0" \ "bet1")) :+ ///
		("_" :+ strofreal((1 .. (kx - 1), 0))), 1) \ "sig0" \ "sig1" \ "rho0" \ "rho1"
	}
	eqs = J(length(colnames), 1, "")
	return((eqs, colnames))
}

void roy_wald::reportResult() {
	printf("\n")
	printf("{txt}%s \n", "{space 3}Wald test for Roy models")
	printf("{txt}%s \n", "{space 3}Ho: selection rule is the original Roy selection mechanism")
	printf("\n")
	if (kx == 1) {
		printf("{txt}%s {res}%g {txt}%s {res}%5.4f \n", "{space 3}t(", df, ")    = ", computeTestStatistic())
	} else if (kx > 1) {
		printf("{txt}%s {res}%g {txt}%s {res}%5.4f \n", "{space 3}chi2(", df, ")    = ", computeTestStatistic())
	}
	if (kx == 1) {
		printf("{txt}%s {res}%5.4f \n", "{space 3}Prob > |t|  = ", computePValue())
	} else if (kx > 1) {
		printf("{txt}%s {res}%5.4f \n", "{space 3}Prob > chi2  = ", computePValue())
	}
	if (checkReject()[2] == 1) {
		printf("{txt}%s \n", "{space 3}Ho is rejected at the 5% level")
	} else if (checkReject()[2] == 0) {
		printf("{txt}%s \n", "{space 3}Ho is not rejected at the 5% level")
	}
}

void roy_wald::saveResult() {
	st_numscalar("N", N)
	st_numscalar("stat", computeTestStatistic())
	st_numscalar("df", df)
	st_numscalar("p", computePValue())
	st_numscalar("reject_10", checkReject()[1])
	st_numscalar("reject_5", checkReject()[2])
	st_numscalar("reject_1", checkReject()[3])
	st_numscalar("kx", kx)
	
	st_matrix("theta", computeCoefs()')
	st_matrix("V", Om)
	st_matrixcolstripe("theta", saveCoefsName())
	st_matrixrowstripe("theta", J(1, 2, " "))
	st_matrixcolstripe("V", saveCoefsName())
	st_matrixrowstripe("V", J(3 * kx + 4, 2, " "))

	if (st_local("V") == "") st_local("V", "mle")
	if (st_local("ADJ") == "") st_local("ADJ", "none")
}
end

// ---------------------------------------------------------------------- lr ---
mata:
void llk_lr(transmorphic M, real rowvector b, real colvector lnf) {
	real colvector p1, p2, p3, p4, p5, p6, p7, p8, y, d
	real colvector z1, z2, z3, r1, r2
	p1 = moptimize_util_xb(M, b, 1) // alp
	p2 = moptimize_util_xb(M, b, 2) // beta0
	p3 = moptimize_util_xb(M, b, 3) // beta
	p4 = moptimize_util_xb(M, b, 4) // sig0^2
	p5 = moptimize_util_xb(M, b, 5) // sig1^2
	p6 = moptimize_util_xb(M, b, 6) // sig0e
	p7 = moptimize_util_xb(M, b, 7) // sig1e
	p8 = moptimize_util_xb(M, b, 8) // sige^2

	y  = moptimize_util_depvar(M, 1)
	d  = moptimize_util_depvar(M, 2)

	z1 = (y :- p2 :- p3) :/ sqrt(p5)
	z2 = -(p3 :+ p1) :/ sqrt(p8)
	z3 = (y :- p2) :/ sqrt(p4)

	r1 = (z2 :- (p7 :/ (sqrt(p5) :* sqrt(p8))) :* z1) :/ sqrt(1 :- (p7 :/ (sqrt(p5) :* sqrt(p8))):^2)
	r2 = (z2 :- (p6 :/ (sqrt(p4) :* sqrt(p8))) :* z3) :/ sqrt(1 :- (p6 :/ (sqrt(p4) :* sqrt(p8))):^2)

	lnf = d :* (lnnormal(-r1) :+ lnnormalden(z1) :- ln(sqrt(p5))) ///
	:+ (1 :- d) :* (lnnormal(r2) :+ lnnormalden(z3) :- ln(sqrt(p4)))
}

class roy_lr {
	public:
	// initialize
	void new()
	
	// setup
	real scalar    N, kx, df
	real rowvector a
	real matrix    X

	// initial value
	real scalar    excnbr
	real rowvector unexcidx
	real colvector bet0_init, bet1_init
	
	real matrix    Xc_init
	real colvector Xhat_init, eps0sq_init, eps1sq_init
	transmorphic   vns_init
	
	real scalar    sig_eta_init, alp1_init, alp2_init, alp0_init, sig0sq_init, sig1sq_init
	real colvector alp_init
	
	void mle()
	
	// estimates unrestricted
	real scalar    sig0sq_u, sig1sq_u, sig0e_u, sig1e_u, sigesq_u, llk_u
	real colvector alp_u, bet0_u, bet_u
	real matrix    Om_u
	transmorphic   M_u

	// estimates restricted
	real scalar    sig0sq_r, sig1sq_r, sig0e_r, sig1e_r, sigesq_r, llk_r
	real colvector alp_r, bet0_r, bet_r
	real matrix    Om_r
	transmorphic   M_r
	
	// lr statistic
	real scalar    stat()
	real scalar    pval()
	real rowvector reject()
	
	// coefficients
	real colvector coef_u()
	real colvector coef_r()

	// report, returns
	string matrix  coefsname()
	void           report()
	void           returns()
}

void roy_lr::new() {
	N        = rows(st_data(., st_local("Y"), st_local("touse")))
	X        = st_data(., st_local("X"), st_local("touse"))
	X        = select(X, colsum(X :== 0) :!= N)
	kx       = cols(X) + 1
	df       = (kx - 1) + 2
	a        = (.1, .05, .01)

	M_u      = moptimize_init()
	M_r      = moptimize_init()

	// exclusion restriction
	xstr     = st_local("X")
	excstr   = st_local("EXC")
	
	xtkn     = tokens(xstr)
	exctkn   = tokens(excstr)
	
	xnbr     = length(xtkn)
	excnbr   = length(exctkn)
	
	unexcidx = colsum(J(excnbr, 1, xtkn) :!= J(1, xnbr, exctkn')) :== excnbr
	unexcstr = invtokens(select(xtkn, unexcidx))
	unexcnbr = sum(unexcidx)

	// initial value
	// step 1
	stata("qui probit " + st_local("D") + st_local("X") + " " + st_local("if"))

	// step 2
	stata("qui predict xg, xb")
	stata("gen lam0 = normalden(xg) / (1 - normal(xg))")
	stata("gen lam1 = -normalden(xg) / normal(xg)")
	stata("qui regress " + st_local("Y") + st_local("X") + " lam0 " + " if" + " " + st_local("D") + " == 0")
	bet0_init    = st_matrix("e(b)")[(1 .. (kx - 1), kx + 1)]'
	stata("qui regress " + st_local("Y") + st_local("X") + " lam1 " + " if" + " " + st_local("D") + " == 1")
	bet1_init    = st_matrix("e(b)")[(1 .. (kx - 1), kx + 1)]'
	
	// step 3
	Xc_init      = X, J(rows(X), 1, 1)
	Xhat_init    = Xc_init * bet1_init - Xc_init * bet0_init
	eps0sq_init  = (st_data(., st_local("Y"), st_local("touse")) - Xc_init * bet0_init):^2
	eps1sq_init  = (st_data(., st_local("Y"), st_local("touse")) - Xc_init * bet1_init):^2
	vns_init     = st_addvar("double", ("Xhat_init", "eps0sq_init", "eps1sq_init"))
	st_store(., vns_init, st_local("touse"), (Xhat_init, eps0sq_init, eps1sq_init))
	
	stata("qui probit " + st_local("D") + " " + unexcstr + " Xhat_init" + " " + st_local("if"))
	xhatidx                             = unexcnbr + 1
	sig_eta_init                        = 1 / st_matrix("e(b)")[xhatidx]
	unexcidxcon                         = (unexcidx, 1)
	alp_init                            = J(1, kx, .)
	alp_init[selectindex(!unexcidxcon)] = J(1, sum(!unexcidxcon), 0)
	alp_init[selectindex(unexcidxcon)]  = st_matrix("e(b)")[((1 .. unexcnbr), unexcnbr + 2)]
	alp_init                            = alp_init'
	
	// step 4
	stata("qui gen xgMills1 = -xg * (normalden(xg) / normal(xg))")
	stata("qui gen xgMills0 = xg * (normalden(xg) / (1 - normal(xg)))")
	stata("qui regress eps0sq xgMills0 " + " if" + " " + st_local("D") + " == 0")
	sig0sq_init = st_matrix("e(b)")[2]
	stata("qui regress eps1sq xgMills1 " + " if" + " " + st_local("D") + " == 1")
	sig1sq_init = st_matrix("e(b)")[2]

	stata("drop xg lam0 lam1 Xhat_init eps0sq_init eps1sq_init xgMills1 xgMills0")
	
	mle()
	
	// unrestricted mle
	alp_u     = moptimize_result_coefs(M_u)[1 .. kx]'
	bet0_u    = moptimize_result_coefs(M_u)[(kx + 1) .. (2 * kx)]'
	bet_u     = moptimize_result_coefs(M_u)[(2 * kx + 1) .. (3 * kx)]'
	
	sig0sq_u  = moptimize_result_coefs(M_u)[(3 * kx + 1)]
	sig1sq_u  = moptimize_result_coefs(M_u)[(3 * kx + 2)]
	sig0e_u   = moptimize_result_coefs(M_u)[(3 * kx + 3)]
	sig1e_u   = moptimize_result_coefs(M_u)[(3 * kx + 4)]
	sigesq_u  = moptimize_result_coefs(M_u)[(3 * kx + 5)]

	llk_u     = moptimize_result_value(M_u)
	
	// restricted mle
	alp_r     = moptimize_result_coefs(M_r)[1 .. kx]'
	bet0_r    = moptimize_result_coefs(M_r)[(kx + 1) .. (2 * kx)]'
	bet_r     = moptimize_result_coefs(M_r)[(2 * kx + 1) .. (3 * kx)]'
	
	sig0sq_r  = moptimize_result_coefs(M_r)[(3 * kx + 1)]
	sig1sq_r  = moptimize_result_coefs(M_r)[(3 * kx + 2)]
	sig0e_r   = moptimize_result_coefs(M_r)[(3 * kx + 3)]
	sig1e_r   = moptimize_result_coefs(M_r)[(3 * kx + 4)]
	sigesq_r  = moptimize_result_coefs(M_r)[(3 * kx + 5)]

	llk_r     = moptimize_result_value(M_r)

	if (st_local("V") == "" | st_local("V") == "mle") {
		Om_u  = moptimize_result_V(M_u)
		Om_r  = moptimize_result_V(M_r)
	} else if (st_local("V") == "robust") {
		Om_u  = moptimize_result_V_robust(M_u)
		Om_r  = moptimize_result_V_robust(M_r)
	}
}

void roy_lr::mle() {
	real matrix    C1
	real colvector c1
	real rowvector idC, C2
	real scalar    idc, c2, tol
	
	tol = strtoreal(st_local("tol"))
	
	// unrestricted mle
	moptimize_init_evaluator(M_u, &llk_lr())
	
	moptimize_init_depvar(M_u, 1, st_data(., st_local("Y"), st_local("touse")))
	moptimize_init_depvar(M_u, 2, st_data(., st_local("D"), st_local("touse")))

	moptimize_init_eq_indepvars(M_u, 1, X)
	moptimize_init_eq_indepvars(M_u, 2, X)
	moptimize_init_eq_indepvars(M_u, 3, X)
	moptimize_init_eq_indepvars(M_u, 4, "")
	moptimize_init_eq_indepvars(M_u, 5, "")
	moptimize_init_eq_indepvars(M_u, 6, "")
	moptimize_init_eq_indepvars(M_u, 7, "")
	moptimize_init_eq_indepvars(M_u, 8, "")

	moptimize_init_tracelevel(M_u, "none")
	moptimize_init_valueid(M_u, "unrestricted log-likelihood")

	moptimize_init_search(M_u, "on")
	moptimize_init_search_random(M_u, "on")

	// identification constraint
	idx = (1 .. kx)[selectindex(!unexcidx)]
	idC = e(idx[1], kx)
	for (i = 2; i <= excnbr; i++) {
		idC = idC \ e(idx[i], kx)
	}
	idC = idC, J(excnbr, 2 * kx + 5, 0)
	idc = J(excnbr, 1, 0)
	moptimize_init_constraints(M_u, (idC, idc))
	
	// initial values
	moptimize_init_eq_coefs(M_u, 1, alp_init')
	moptimize_init_eq_coefs(M_u, 2, bet0_init')
	moptimize_init_eq_coefs(M_u, 3, (bet1_init :- bet0_init)')
	moptimize_init_eq_coefs(M_u, 4, sig0sq_init)
	moptimize_init_eq_coefs(M_u, 5, sig1sq_init)
	moptimize_init_eq_coefs(M_u, 8, sig_eta_init^2)

	if (tol != 999) {
		moptimize_init_conv_ptol(M_u, tol)
		moptimize_init_conv_vtol(M_u, tol)
		moptimize_init_conv_nrtol(M_u, tol)
	}

	moptimize(M_u)

	// restricted mle
	moptimize_init_evaluator(M_r, &llk_lr())

	moptimize_init_depvar(M_r, 1, st_data(., st_local("Y"), st_local("touse")))
	moptimize_init_depvar(M_r, 2, st_data(., st_local("D"), st_local("touse")))

	moptimize_init_eq_indepvars(M_r, 1, X)
	moptimize_init_eq_indepvars(M_r, 2, X)
	moptimize_init_eq_indepvars(M_r, 3, X)
	moptimize_init_eq_indepvars(M_r, 4, "")
	moptimize_init_eq_indepvars(M_r, 5, "")
	moptimize_init_eq_indepvars(M_r, 6, "")
	moptimize_init_eq_indepvars(M_r, 7, "")
	moptimize_init_eq_indepvars(M_r, 8, "")

	moptimize_init_tracelevel(M_r, "none")
	moptimize_init_valueid(M_r, "restricted log-likelihood")

	// constraints
	C1 = I(kx), J(kx, 2 * kx + 5, 0)
	c1 = J(kx, 1, 0)
	C2 = J(1, 3 * kx + 2, 0), 1, -1, 1
	c2 = 0
	C3 = J(1, 3 * kx, 0), -1, 1, -1, -1, 0
	c3 = 0

	moptimize_init_search(M_r, "on")
	moptimize_init_search_random(M_r, "on")

	moptimize_init_constraints(M_r, ((C1 \ C2 \ C3), (c1 \ c2 \c3)))
	
	// initial values
	moptimize_init_eq_coefs(M_r, 1, alp_init')
	moptimize_init_eq_coefs(M_r, 2, bet0_init')
	moptimize_init_eq_coefs(M_r, 3, (bet1_init :- bet0_init)')
	moptimize_init_eq_coefs(M_r, 4, sig0sq_init)
	moptimize_init_eq_coefs(M_r, 5, sig1sq_init)
	moptimize_init_eq_coefs(M_r, 8, sig_eta_init^2)
	
	if (tol != 999) {
		moptimize_init_conv_ptol(M_r, tol)
		moptimize_init_conv_vtol(M_r, tol)
		moptimize_init_conv_nrtol(M_r, tol)
	}
	
	moptimize(M_r)
}

real colvector roy_lr::coef_u() {
	return((alp_u \ bet0_u \ bet_u \ sig0sq_u \ sig1sq_u \ sig0e_u \ sig1e_u \ sigesq_u))
}

real colvector roy_lr::coef_r() {
	return((alp_r \ bet0_r \ bet_r \ sig0sq_r \ sig1sq_r \ sig0e_r \ sig1e_r \ sigesq_r))
}

real scalar roy_lr::stat() {
	return(-2 * (llk_r - llk_u))
}

real scalar roy_lr::pval() {
	return(1 - chi2(df, stat()))
}

real rowvector roy_lr::reject() {
	return(pval() :< a)
}

string matrix roy_lr::coefsname() {
	string colvector colnames, eqs

	colnames = colshape(J(1, kx, ("alp" \ "bet0" \ "bet")) :+ ///
	("_" :+ strofreal((1 .. (kx - 1), 0))), 1) \ "sig0sq" \ "sig1sq" \ "sig0e" \ "sig1e" \ "sigesq"
	eqs = J(length(colnames), 1, "")
	return((eqs, colnames))
}

void roy_lr::report() {
	printf("\n")
	printf("{txt}%s \n", "{space 3}LR test for Roy models")
	printf("{txt}%s \n", "{space 3}Ho: selection rule is the original Roy selection mechanism")
	printf("\n")
	printf("{txt}%s {res}%g {txt}%s {res}%5.4f \n", "{space 3}chi2(", df, ")    = ", stat())
	printf("{txt}%s {res}%5.4f \n", "{space 3}Prob > chi2  = ", pval())
	if (reject()[2] == 1) {
	printf("{txt}%s", "{space 3}Ho is rejected at the 5% level")
	} else if (reject()[2] == 0) {
		printf("{txt}%s", "{space 3}Ho is not rejected at the 5% level")
	}
}

void roy_lr::returns() {
	st_numscalar("N", N)
	st_numscalar("kx", kx)
	st_numscalar("df", df)
	
	st_numscalar("stat", stat())
	st_numscalar("p", pval())
	st_numscalar("reject_10", reject()[1])
	st_numscalar("reject_5", reject()[2])
	st_numscalar("reject_1", reject()[3])
	st_numscalar("llk_u", llk_u)
	st_numscalar("llk_r", llk_r)
	
	st_matrix("theta_u", coef_u()')
	st_matrix("theta_r", coef_r()')
	st_matrix("V_u", Om_u)
	st_matrix("V_r", Om_r)
	st_matrixcolstripe("theta_u", coefsname())
	st_matrixrowstripe("theta_u", J(1, 2, " "))
	st_matrixcolstripe("theta_r", coefsname())
	st_matrixrowstripe("theta_r", J(1, 2, " "))
	st_matrixcolstripe("V_u", coefsname())
	st_matrixrowstripe("V_u", J(3 * kx + 5, 2, " "))
	st_matrixcolstripe("V_r", coefsname())
	st_matrixrowstripe("V_r", J(3 * kx + 5, 2, " "))

	if (st_local("V") == "") st_local("V", "mle")
}
end

// ---------------------------------------------------------------------- sm ---
mata:
class roy_sm {
	public:
	// initialize
	void new()
	
	// setup
	string scalar G
	real   scalar N, tn, q1, eta, eps, B, kn, bn, alp
	
	// variable
	real colvector Y, Z
	
	// support transformation
	real colvector zstar
	void           zstar()
	
	// tau
	real rowvector tau
	void           tau()
	
	// ell
	real scalar    ln
	void           ell()
	
	// Q(tau, ell)
	real colvector ind_q
	void           Q()
	
	// f(W, tau)
	real matrix    f
	void           f()
	
	// test statistic, critical value
	real scalar    Tn, cv, cv10, cv5, cv1
	void           Tncv()
	
	// result
	void           report()
	void           returns()
}

void roy_sm::new() {
	// variable
	Y   = st_data(., st_local("Y"), st_local("touse"))
	Z   = st_data(., st_local("Z"), st_local("touse"))
	
	// setup
	N   = rows(Y)
	
	G   = st_local("G")
	tn  = strtoreal(st_local("tn"))
	eta = strtoreal(st_local("eta"))
	eps = strtoreal(st_local("eps"))
	B   = strtoreal(st_local("boot"))
	alp = strtoreal(st_local("alpha"))
	
	q1  = floor((N / 15))
	kn  = 0.15 * log(N)
	bn  = 0.85 * (log(N) / (log(log(N))))
	
	if (strtoreal(st_local("qn")) != 999) q1 = strtoreal(st_local("qn"))
	if (strtoreal(st_local("kn")) != 999) kn = strtoreal(st_local("kn"))
	if (strtoreal(st_local("bn")) != 999) bn = strtoreal(st_local("bn"))
	
	// support transformation
	zstar = .z
	
	// tau
	tau   = .z
	
	// ell
	ln    = .z
	
	// Q(tau, ell)
	ind_q = .z
	
	// f(W, tau)
	f     = .z
	
	// test statistic
	Tn    = .z
	
	// critical value
	cv   = .z
	cv10 = .z
	cv5  = .z
	cv1  = .z
}

void roy_sm::zstar() {
	real scalar zmax, zmin
	
	if (zstar == .z) {
		zmax  = max(Z)
		zmin  = min(Z)
		zstar = (Z :- zmin) :/ (zmax - zmin)
	}
}

void roy_sm::tau() {
	real scalar ymax, ymin, inc
	
	if (tau  == .z) {
		ymax = max(Y)
		ymin = min(Y)
		inc  = (ymax - ymin) / (tn + 1)
		tau  = ymin :+ inc :* (1 .. tn)
	}
}

void roy_sm::ell() {
	real scalar l, q
	
	if (ln == .z) {
		l = 0
		for (q = 2; q <= q1; q++) {
			l = l + q * (q - 1) / 2
		}
		ln = l
	}
}

void roy_sm::Q() {
	real scalar k, q, i, j
	real rowvector ind
	real matrix ind_ln, ind_l
	
	if (ind_q == .z) {
	if (G == "countable") {
		k = 1
		ind_ln = J(ln, 4, .)
		
		for (q = 2; q <= q1; q++) {
			for (i = 1; i <= q; i++) {
				for (j = (i + 1); j <= q; j++) {
					ind = (q, q * (q - 1) / 2, j, i)
					ind_ln[k, ] = ind
					k = k + 1
				}
			}
		}
		ind_l = (ind_ln[, 1]:^(-2)) :/ ind_ln[, 2]
	} else if (G == "continuum") {
		ind_l = (1 / ln) :* J(ln, 1, 1)
	}
		ind_q = J(tn, 1, 1 / tn) # ind_l
	}
}

void roy_sm::f() {
	if (f == .z) f = -(J(1, tn, Y) :<= J(N, 1, tau))
}

void roy_sm::Tncv() {
	real matrix wnmn_lntn, wnmn_ln, zstarq1, zstarq2, g0, mi0
	real colvector vn_lntn, Tn_lntn, sig_lntn, phi_lntn
	real colvector vn_ln, Tn_ln, sig_ln, phi_ln
	real colvector ft, mi2, wi2, mi1, wi1, wnmn
	real colvector cv_u, zero_lntn, ui, phiu_lntn, Tn_u1, Tn_u2
	real rowvector qi1, qi2, mn0, wn0
	real scalar t, k, q, i, j, r, qi0, qn, mn2, wn2, mn1, wn1, vn, sig2, sig, Tn1, Tn2, phi, t1, t2, Tn_u
	
	vn_lntn   = J(ln * tn, 1, .)
	wnmn_lntn = J(N, ln * tn, .)
	Tn_lntn   = J(ln * tn, 1, .)
	sig_lntn  = J(ln * tn, 1, .)
	phi_lntn  = J(ln * tn, 1, .)

	vn_ln     = J(ln, 1, .)
	wnmn_ln   = J(N, ln, .)
	Tn_ln     = J(ln, 1, .)
	sig_ln    = J(ln, 1, .)
	phi_ln    = J(ln, 1, .)

	for (t = 1; t <= tn; t++) {
		ft = f[, t]
		k  = 1
	
		for (q = 1; q <= (q1 - 1); q++) {
			if (G == "countable") {
				r   = 1 / (q + 1)
				qi0 = r
			} else if (G == "continuum") {
				r   = q / q1
				qi0 = 1 / q1
			}
			qi1     = (0 .. q) * r
			qi2     = qi1 :+ r
			qn      = length(qi1)
			zstarq1 = J(1, qn, zstar) :- J(N, 1, qi1)
			zstarq2 = J(1, qn, zstar) :- J(N, 1, qi2)
			g0      = (zstarq1 :>= 0) :* (zstarq2 :<= 0)
			
			mi0     = J(1, qn, ft) :* g0
			mn0     = mean(mi0)
			wn0     = mean(g0)

			for (i = 1; i <= qn; i++) {
				mi2 = mi0[, i]
				wi2 = g0[, i]
				mn2 = mn0[i]
				wn2 = wn0[i]
			
				for (j = (i + 1); j <= qn; j++) {
					mi1  = mi0[, j]
					wi1  = g0[, j]
					mn1  = mn0[j]
					wn1  = wn0[j]

					vn   = mn2 * wn1 - mn1 * wn2
					wnmn = (wn1 :* (mi2 :- mn2) + mn2 :* (wi1 :- wn1) - wn2 :* (mi1 :- mn1) - mn1 :* (wi2 :- wn2))
					sig2 = mean(wnmn:^2)
					sig2 = max((sig2, eps))
					sig  = sqrt(sig2)
					Tn1  = sqrt(N) * vn / sig
					Tn2  = max((Tn1, 0))^2
					phi  = -bn * (Tn1 < -kn)
					Tn_ln[k, 1]  = Tn2
					vn_ln[k, 1]  = vn
					phi_ln[k, 1] = phi
					sig_ln[k, 1] = sig
					wnmn_ln[, k] = wnmn
					k = k + 1
				}
			}
		}
		t1 = 1 + (t - 1) * ln
		t2 = t * ln
		Tn_lntn[t1::t2, 1]    = Tn_ln
		vn_lntn[t1::t2, 1]    = vn_ln
		phi_lntn[t1::t2, 1]   = phi_ln
		sig_lntn[t1::t2, 1]   = sig_ln
		wnmn_lntn[, t1 .. t2] = wnmn_ln
	}
	Tn = cross(Tn_lntn, ind_q)

	// critical value
	cv_u = J(B, 1, .)
	zero_lntn = J(ln * tn, 1, 0)
	
	for (b = 1; b <= B; b++) {
		ui = 2 * sqrt(3) :* runiform(N, 1, 0, 1) :- sqrt(3)
		phiu_lntn = (1 / sqrt(N)) :* cross(wnmn_lntn, ui)
		Tn_u1 = (phiu_lntn :/ sig_lntn) + phi_lntn
		Tn_u2 = rowmax((Tn_u1, zero_lntn)):^2
		Tn_u = cross(Tn_u2, ind_q)
		cv_u[b] = Tn_u
	}
	cv   = sort(cv_u, 1)[floor((B + 1) * (1 - alp + eta))] + eta
	cv10 = sort(cv_u, 1)[floor((B + 1) * (1 - 0.1 + eta))] + eta
	cv5  = sort(cv_u, 1)[floor((B + 1) * (1 - 0.05 + eta))] + eta
	cv1  = sort(cv_u, 1)[floor((B + 1) * (1 - 0.01 + eta))] + eta
}

void roy_sm::report() {
	printf("\n")
	printf("{txt}%s \n", "{space 3}SM test for Roy models")
	printf("{txt}%s \n", "{space 3}Ho: selection rule is the original Roy selection mechanism")
	printf("{txt}%s {res}%5.4f \n", "{space 3}Tn = ", Tn)
	printf("{txt}%s {res}%5.4f \n", "{space 3}cv = ", cv)
	if (Tn > cv) {
		printf("{txt}%s {res}%1.0f {txt}%s \n", "{space 3}Ho is rejected at the", alp * 100, "% level")
	} else if (Tn <= cv) {
		printf("{txt}%s {res}%1.0f {txt}%s \n", "{space 3}Ho is not rejected at the", alp * 100, "% level")
	}
}

void roy_sm::returns() {
	zstar()
	tau()
	ell()
	Q()
	f()
	Tncv()
	
	st_numscalar("N", N)
	st_numscalar("Tn", Tn)
	st_numscalar("cv_10", cv10)
	st_numscalar("cv_5", cv5)
	st_numscalar("cv_1", cv1)
	st_numscalar("reject_10", (Tn > cv10))
	st_numscalar("reject_5", (Tn > cv5))
	st_numscalar("reject_1", (Tn > cv1))
}
end
