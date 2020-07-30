clear all
* Code to replicate Gilcrhsit and Zakrajsek (2012) 
cd "/Users/Alex/Dropbox/CS/RRFinReplicationFiles"
use gzdata

local repetitions = "100"
local spec        = "ebp_oa_dd_spec2"
gen time = _n
tsset time

** Get shock
quietly: reg ebp L(0/2).dlcons L(0/2).dlinvest L(0/2).dlgdp L(0/2).dlpgdp L(1/2).ebp L(1/2).vwx L(1/2).treas L(1/2).ffr
predict resid, r

* dlcons: log-difference of real consumption
* dlinvest: log-difference of real investment
* dlpgdp: log-difference of real gdp
* dlpgdp: log-difference of price deflator
* ebp_oa_avg: quarterly average of ebp
* vwx: VIX volatility index
* treas: 10 years treasury yield
* ffr: federal funds rate

local varlist     = "dlcons dlinvest dlgdp dlpgdp ebp_oa_avg vwx treas ffr "
local impulse1    = "ebp_oa_avg"

var `varlist', lags(1/2) level(95) dfk
varstable
irf create `spec', step (21) set(`spec') bsp reps(`repetitions') replace

*Save two files for two different sets
*(Impulse 1)
drop _all
use `spec'.irf
drop irf stdirf sirf stdsirf cirf stdcirf dm stddm cdm stdcdm mse sfevd stdsfevd irfname
order step response impulse oirf stdoirf fevd stdfevd coirf stdcoirf
sort impulse response step
keep if impulse == "`impulse1'"
save impulse_`impulse1'_`spec', replace
outsheet using "impulse_`impulse1'_`spec'.csv", replace comma
