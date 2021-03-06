clear all
set seed 12345678
local N_runs = 100000
local scenario=0
local scenariocount=1

/****** SIMULATION ******/
quietly { 
foreach true_tx_effect in  1   {	
	foreach true_interaction_effect in  0  1   {
		foreach n_trial in   200   500  1000  {
			foreach n_total in 2000   5000   10000  {
				foreach selection_effect in   0    1  {

tempname saved_results 

postfile `saved_results' double(	scenario true_tx_effect ///
									true_interaction_effect ///
									selection_effect ///
									n_trial n_total ///
									popS0_d trial_d ///
										IPW1_d IPW2_d OM_d DR1_d DR2_d DR3_d ///
									popS0_1 trial_1 ///
										IPW1_1 IPW2_1 OM_1 DR1_1 DR2_1 DR3_1 ///
									popS0_0 trial_0 ///
										IPW1_0 IPW2_0 OM_0 DR1_0 DR2_0 DR3_0) ///
								using linear_norm_composite_`scenariocount'.dta, replace 

local ++scenario
local ++scenariocount

forvalues i = 1/`N_runs' {
di in red " scenario=`scenario', simulation = `i', true_tx_effect=`true_tx_effect'"
clear
/**** Data generation ****/
local n_total = `n_total'
set obs `n_total'
generate X1 = rnormal()
generate X2 = rnormal() 
generate X3 = rnormal()
if `selection_effect' == 0 {
if `n_trial' == 200 & `n_total' == 2000 {
		local b_0 =  -2.880859375
}
if `n_trial' == 500 & `n_total' == 2000 {
		local b_0 =   -1.491699219 
}
if `n_trial' == 1000 & `n_total' == 2000 {
		local b_0 = 0
}
if `n_trial' == 200 & `n_total' == 5000 {
		local b_0 = -4.023437500
}
if `n_trial' == 500 & `n_total' == 5000 {
		local b_0 =   -2.880859375
}
if `n_trial' == 1000 & `n_total' == 5000 {
		local b_0 = -1.866760254   
}
if `n_trial' == 200 & `n_total' == 10000 {
		local b_0 =  -4.799804688
}
if `n_trial' == 500 & `n_total' == 10000 {
		local b_0 =  -3.754882812 
}
if `n_trial' == 1000 & `n_total' == 10000 {
		local b_0 =  -2.880859375  
}
}
if `selection_effect' == 1 {
if `n_trial' == 200 & `n_total' == 2000 {
		local b_0 =      -3.159180  
}
if `n_trial' == 500 & `n_total' == 2000 {
		local b_0 =   -1.645508
}
if `n_trial' == 1000 & `n_total' == 2000 {
		local b_0 = 0
}
if `n_trial' == 200 & `n_total' == 5000 {
		local b_0 =  -4.384766
}
if `n_trial' == 500 & `n_total' == 5000 {
		local b_0 = -3.159180 
}
if `n_trial' == 1000 & `n_total' == 5000 {
		local b_0 = -2.055359 
}
if `n_trial' == 200 & `n_total' == 10000 {
		local b_0 = -5.214844 
}
if `n_trial' == 500 & `n_total' == 10000 {
		local b_0 = -4.101562
}
if `n_trial' == 1000 & `n_total' == 10000 {
		local b_0 =  -3.159180
}
}
generate S = runiform() < invlogit(`b_0' + `selection_effect' * X1 + X2 +X3)
*assign tx, outcome for trial participants
generate A = rbinomial(1, 0.5) if S==1
replace A = 0 if A == .
/*potential outcomes*/
/* potential outcome for a=1*/ 
generate Y_1 = X1 +  X2 + X3  + `true_tx_effect' + `true_interaction_effect' * X1 + rnormal() 
/* potential outcome for a=0*/ 
generate Y_0 =  X1 +  X2 + X3  + rnormal()
generate Y_1_mean = X1 +  X2 + X3  + `true_tx_effect' + `true_interaction_effect' * X1 generate Y_0_mean =  X1 +  X2 + X3 
*whole target population
gen d_pop = Y_1 - Y_0  /*difference of treatment effect in population*/
gen d_pop_mean = Y_1_mean - Y_0_mean 
*target outside trial 
summ d_pop_mean if S==0 
	local popS0_d = r(mean)
summ Y_1_mean if S==0 
	local popS0_1 = r(mean)
summ Y_0_mean if S==0 
	local popS0_0 = r(mean) 
* observed outcome
generate Y = Y_1 * A + Y_0 * (1 - A) if S==1 
/**** Estimation ****/
/* estimate for trial */
regress  Y A if S == 1 
matrix estimates = e(b)
	local trial_d = estimates[1,1] 
	local trial_1 = estimates[1,1] + estimates[1,2]
	local trial_0 = estimates[1,2]
/* estimate the weights */
logistic S X1 X2 X3
predict ps 
logistic A X1 X2 X3 if S == 1 
predict pa if S == 1
generate w_corr2 = (A*S * (1-ps) ) /(ps*pa) + ( (1-A)*S * (1-ps) ) /(ps*(1-pa))
replace w_corr2 = 0 if w_corr2 == .
replace A = 0 if A == .
replace Y = 0 if Y == .
/* IP weighting, IPW1*/
summ S
	local one_over = 1 / ( 1 - r(mean) )
generate summand1 = A * S * w_corr2 * Y
summ summand1
local IPW1_1 = `one_over' * r(mean)
generate summand0 = (1 - A) * S * w_corr2 * Y
summ summand0
local IPW1_0 = `one_over' * r(mean)
local IPW1_d = `IPW1_1' - `IPW1_0'
/* IP weighting, IPW2*/
regress Y A [pw=w_corr2] if S==1 
matrix estimates = e(b)
local IPW2_d = estimates[1,1]
local IPW2_1 = estimates[1,1] + estimates[1,2]
local IPW2_0 = estimates[1,2]
/* Outcome regression (OR) */
regress Y X1 X2 X3  if A== 1 & S == 1
	predict p_1_12 
    summ p_1_12 if S == 0
    local OM_1 = r(mean)
regress  Y X1 X2 X3 if A == 0 & S == 1
    predict p_0_12 
	summ p_0_12 if S == 0
    local OM_0 = r(mean)
local OM_d = `OM_1' - `OM_0' /*estimate effect for A (treatment)*/

/* DR1 */
summ S
	local one_over = 1 / ( 1 - r(mean) )
generate first_term_a1 = A * S * w_corr2 * Y
generate second_term_a1 = (1 - S) * p_1_12
generate third_term_a1 = - A * S * w_corr2 * p_1_12
generate sum_a1 = first_term_a1 + second_term_a1 + third_term_a1
summ sum_a1
local DR1_1 = `one_over' * r(mean)
generate first_term_a0 = (1 - A) * S * w_corr2 * Y
generate second_term_a0 = (1 - S) * p_0_12
generate third_term_a0 = - (1 - A) * S * w_corr2 * p_0_12
generate sum_a0 = first_term_a0 + second_term_a0 + third_term_a0
summ sum_a0
local DR1_0 = `one_over' * r(mean)
local DR1_d = `DR1_1'-`DR1_0' 

/* DR2 */
replace sum_a1 = first_term_a1 + third_term_a1
summ sum_a1
	local sum_a1 = r(mean)
generate normalization1 = A * S * w_corr2
	sum normalization1
	local norm1 = 1 / r(mean) 
sum second_term_a1
	local sum_b1 = r(mean)
local DR2_1 = `one_over' *  `sum_b1'   + `norm1' * `sum_a1'
replace sum_a0 = first_term_a0  + third_term_a0
summ sum_a0
	local sum_a0 = r(mean)
generate normalization0 = (1 - A) * S * w_corr2
	sum normalization0
	local norm0 = 1 / r(mean) 
summ second_term_a0
	local sum_b0 = r(mean)
local DR2_0 = `one_over' * `sum_b0'  + `norm0' *  `sum_a0'
local DR2_d = `DR2_1'-`DR2_0' 

/* DR3 */
regress Y X1 X2 X3 [pw=w_corr2]  if A== 1 & S==1
	predict Y_11 
    summ Y_11 if S == 0
    local DR3_1 = r(mean)
regress Y X1 X2 X3 [pw=w_corr2]  if A== 0 & S==1
    predict Y_01 /*if R==0*/
    summ Y_01 if S == 0
    local DR3_0 = r(mean)
local DR3_d = `DR3_1'-`DR3_0' 

post `saved_results' (`scenario') (`true_tx_effect') (`true_interaction_effect') (`selection_effect') ///
					 (`n_trial') (`n_total') ///
					 (`popS0_d') (`trial_d') ///
						(`IPW1_d') (`IPW2_d') (`OM_d') (`DR1_d') (`DR2_d') (`DR3_d') ///
					 (`popS0_1') (`trial_1') ///
						(`IPW1_1') (`IPW2_1') (`OM_1') (`DR1_1') (`DR2_1') (`DR3_1') ///
					 (`popS0_0') (`trial_0') ///
						(`IPW1_0') (`IPW2_0') (`OM_0') (`DR1_0') (`DR2_0') (`DR3_0')
}
postclose `saved_results'
}
}	
}	
}
}
}
//
