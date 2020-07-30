clear all
 
*** Replicate Romer and Romer (2017) 
*** Folder & data
cd "/Users/Alex/Dropbox/CS/RRFinReplicationFiles"
use data

* Keep only US
keep if id == 24
 
*** Country and time tariables
gen time = _n

tsset time
 
 *** Dependent variable
gen y = ln(gdp)
local lhs "y"   /* LeftHandSide variables: variables that IRFs are created for */
local H = 11            /* Impulse response horizon */
local lag = 4       /* Lag length */  
 
 *** Crisis variable
gen x = crisis
gen z = ln(gdp)

foreach v in `lhs'  {
quietly:    tsrevar F(0/`H').`v'
quietly:    rename (`r(varlist)') `v'_#, addnumber
quietly:    gen `v' = `v'_1
}

* 3.7 EA, 7 US
foreach v in `lhs'  {
forvalues i=1/`H' {
ivreg2 `v'_`i'  L(0/`lag').x L(1/`lag').y, robust bw(auto)
quietly: est sto irf_`v'_`i'
* Rescale by 7, i.e. effect of a moderate financial crisis-minus
quietly: generate beta_`v'_`i' = _b[x]*7
quietly: generate se_`v'_`i' = _se[x]*7
}
}
 
*preserve
keep beta* se*
drop if _n > 1
gen i=_n
reshape long beta_y_ se_y_, i(i) j(time)
gen t = time -1

foreach v in `lhs'  {
gen ub_`v' = beta_`v'_ + 1.645*se_`v'_ /* 90% CI: beta +/- 1.645*SE */
gen lb_`v' = beta_`v'_ - 1.645*se_`v'_ 
twoway (line beta_`v'_ t ) (line ub_`v' t ) (line lb_`v' t), name(symplot_`v')  legend(off) ytitle(`v')
}
