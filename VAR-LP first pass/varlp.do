**************************************************************************************
** Barnichon, Matthes, Ziegenbein: Are the effects of financial shocks big or small? *
** Stata code to run the VAR-LP robustness check  ************************************
**************************************************************************************

clear all

****************************
*** Folder & data **********
****************************
cd "/Users/Alex/Dropbox/Creditpaper/"
use "data_lp.dta"

****************************
*** Variables *************
****************************

* gdp: annual GDP growth
* cpi: annual cpi inflation
* ebp: excess bond premium
* ffr: effective federal funds rate
* gz: Gilchrist & Zakrajsek (2012) credit spread
* baa: Baa spread

****** Generate time series of integer numbers to describe time
gen time = _n 
tsset time

* Only use pre-2007 data
*drop if date > 2007

* Run recursive VAR and obtain shocks
quietly: reg ebp L(1/6).ebp L(1/6).ffr L(0/6).gdp L(0/6).cpi
predict shock, r

****** Generate interaction terms: bad = contractionary shock to EBP
gen thr = 0
gen bad = .
replace bad = 0 if shock <= 0
replace bad = 1 if shock > 0


****************************
*** Regressions ************
****************************

***** LP settings
gen y = gdp			/* choose LHS variable */
local H = 60        /* Impulse response horizon */
local lag = 6 		/* Lag length */  
local rhs "y x ebp" /* LHS and RHS variables */
gen x = shock     /* Impulse of interest */  


foreach v in `rhs' {
*foreach l=0/`lag' {
quietly: gen bad_`v' = `v'*bad
quietly: gen good_`v' = `v'*(1-bad)
*}
}

***** Run regressions and save coefficients
forvalues i=0/`H' {
quietly: ivreg2 F`i'.y  good_x bad_x  L(1/`lag').y L(1/`lag').good_x, bw(auto) r 
quietly: generate betag_y_`i' = _b[good_x]
quietly: generate seg_y_`i' = _se[good_x]
quietly: generate betab_y_`i' = _b[bad_x]
quietly: generate seb_y_`i' = _se[bad_x]
}

***** Create data file that saves results -> just so we can make nicer figures in matlab later
keep beta* se* 
drop if _n > 1
gen i=_n
reshape long betag_y_ seg_y_ betab_y_ seb_y_, i(i) j(time)
gen t = time -1

***** Figures
gen ubg_y = betag_y_ + seg_y_
gen lbg_y = betag_y_ - seg_y_
gen ubb_y = betab_y_ + seb_y_
gen lbb_y = betab_y_ - seb_y_
set graphics off
twoway (line betag_y_ t ) (line ubg_y t ) (line lbg_y t), name(plotg_y)  legend(off) ytitle(good times)
twoway (line betab_y_ t ) (line ubb_y t ) (line lbb_y t), name(plotb_y)  legend(off) ytitle(bad times)


set graphics on
graph combine plotg_y plotb_y

outsheet using "varlpirf.xls", replace comma
** Important: IRFs have to be cumulated and rescaled by the effect of the shock on EBP.
