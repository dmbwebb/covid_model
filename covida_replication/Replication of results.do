
//Do-file that replicates the figures / tables from the CoVIDA data source.

//Keep in mind that the stratum 2 includes individuals in stratums 1 & 2 and the stratum 5 includes individuals in the stratums 5 & 6

ssc install tabout
set more off

// In order to run the do file, you ned to change these globals with the location of the replication file in your computer

global processed "C:\Users\rlaaj\Dropbox (Uniandes)\PROJECTS\COVID-19\Research COVID-19\Inequalitiy in COVID\Replication"
global results "C:\Users\rlaaj\Dropbox (Uniandes)\PROJECTS\COVID-19\Research COVID-19\Inequalitiy in COVID\Replication/Results/"

****FIGURE S9: EFFECT OF HOUSEHOLD SIZE CONDITIONAL ON STRATUM 

use "${processed}/Datos_Salesforce_treated_mar03_newvars.dta", clear 

*Panel B
areg positive i.stratum#c.hhsize i.stratum, ro a(test_month)
outreg2 using "${results}/Figure S9 Panel B.xls", replace bracket  noni less(1) nor2 noaster 
//Figure S9 Panel B is constructed from the regression, refer to the excel "Additional figures.xls", sheet "Figure S9 Panel B"

****TABLE S1: REGRESSIONS ABOUT HOW CAREFUL THEY ARE (Self-declared prevention practices)
foreach x in  trip_protec lavadomanos_protec tiempolavadomanos_protec usatapabocasensudíaadía_protec usagelantibacterial_protec  {
qui areg `x' i.stratum , ro a(test_month)
local p_val_Ftest = Ftail(e(df_m), e(df_r), e(F))
outreg2 using "${results}/Table S1.xls", append bracket  noni less(1) nor2 addstat("p_val_Ftest", `p_val_Ftest') 
}


****TABLE S2: NUMBER OF DAYS WORKED OUTSIDE HOME AS A FUNCTION OF SICKNESS AND SES 
* We need to exclude days of holidays when a high proportion was not working 
reg work_outside_home i.stratum i.stratum#c.fraction_sick [aw=weight_ocup] if symptom == 1 & noholiday ==1, ro // excludes holidays after Dec 22nd until Jan 21s
outreg2 using "${results}/Table S2.xls", replace bracket  noni less(1) nor2 


****TABLE S3: STATISTICAL PREDICTIONS OF THE NUMBER OF DAYS WORKED USING RESULTS FROM TABLE S2
//Refer to the excel "Additional figures.xls", sheet "Table S3"


****TABLE S6 PANEL B: OLS CORRELATION BETWEEN AGE AND POSITIVITY RATE 
reg positive age , ro
outreg2 using "${results}/Table S6 Panel B.xls", replace bracket  less(1) 


****Table S7: Potential Determinants of Infection Estimated by SES

////Parameter: Infections outside the home

*Days working outside of home (in previous 10 working days)
tabform work_outside_home using "${results}\work_outside_tab.xls" if symptom ==0 & contact == 0 & noholiday ==1 [aw=weight_ocup], by(stratum) bdec(3)  se ci level(95) sdbracket cibrace vertical

	*Regressions to show F-stat 
qui reg work_outside_home i.stratum if symptom ==0 & contact == 0 & noholiday ==1 [aw=weight_ocup], ro
local p_val_Ftest = Ftail(e(df_m), e(df_r) , e(F))
outreg2 using "${results}\work_outside_reg.xls" , replace bracket  noni less(1) nor2 noaster addstat("p_val_Ftest", `p_val_Ftest') 

use "${processed}/TRACING_mar03.dta", clear

*Number of non-work contacts outside of home (time frame)
preserve
collapse (sum) outside_nonwork, by(index_cedula index_stratum)
tabform outside_nonwork using "${results}\contacts_nonwork_tabform.xls", by(index_stratum) se ci level(95) bdec(3) sdbracket cibrace vertical mtbdec(2)

	*Regressions to show F-stat 
qui reg outside_nonwork i.index_stratum , ro
local p_val_Ftest = Ftail(e(df_m), e(df_r), e(F))
outreg2 using "${results}\contacts_nonwork_reg.xls", replace bracket  noni less(1) nor2 noaster addstat("p_val_Ftest", `p_val_Ftest') 
restore

use "${processed}/TRACING_mar03.dta", clear

*Secondary attack rate (outside of home)
tabform out_positive using "${results}\SAR_out_tabform.xls", by(index_stratum) se ci level(95) bdec(3) sdbracket cibrace vertical mtbdec(2)

	*Regressions to show F-stat 
qui reg out_positive i.index_stratum , ro
local p_val_Ftest = Ftail(e(df_m), e(df_r), e(F))
outreg2 using "${results}\SAR_out_reg.xls", replace bracket  noni less(1) nor2 noaster addstat("p_val_Ftest", `p_val_Ftest') 

*Contact matrix structure
tabout index_stratum contact_stratum if household == 0 using "${results}\contact_matrix.xls" , replace  


////parameter: Infections inside home

*Secondary attack rate (inside home)
tabform hh_positive using "${results}\SAR_hh_tabform.xls", by(index_stratum) se ci level(95) bdec(3) sdbracket cibrace vertical mtbdec(2)

	*Regressions to show F-stat 
qui reg hh_positive i.index_stratum, ro
local p_val_Ftest = Ftail(e(df_m), e(df_r), e(F))
outreg2 using "${results}\SAR_out_reg.xls.xls", replace bracket  noni less(1) nor2 noaster addstat("p_val_Ftest", `p_val_Ftest') 

////Parameter: Isolation behavior

use "${processed}/Seguimiento.dta", clear 

*Isolation after positive test result (%self declared)
tabform isolates using "${results}\tested_isolates_tab.xls", by(stratum) se ci level(95) bdec(2) sdbracket cibrace vertical mtbdec(2)
	
	*Regressions to show F-stat 
qui reg isolates i.stratum  , ro
local p_val_Ftest = Ftail(e(df_m), e(df_r), e(F))
outreg2 using "${results}\tested_isolates_reg.xls", replace bracket  noni less(1) nor2 noaster addstat("p_val_Ftest", `p_val_Ftest')

*# of days worked when has symptoms
//Data from Table S3

*# days worked when knowing about positive contact 
//Result comes from the possitive detected cases dataset.

use "${processed}/Datos_Salesforce_treated_mar03_newvars.dta", clear 

*# days worked when someone is tested positive in same household
tabform work_outside_home using "${results}\tabform_work_covid_hogar.xls" if symptom ==0 & covid_hogar == 1 & noholiday ==1 [aw=weight_ocup], by(stratum) bdec(3)  se ci level(95) sdbracket cibrace vertical

	*Regressions to show F-stat (not significant)
qui reg work_outside_home i.stratum if symptom ==0 & covid_hogar == 1 & noholiday ==1 [aw=weight_ocup], ro
local p_val_Ftest = Ftail(e(df_m), e(df_r), e(F))
outreg2 using "${results}\work_outside_reg.xls" , append bracket  noni less(1) nor2 noaster addstat("p_val_Ftest", `p_val_Ftest') 


 

