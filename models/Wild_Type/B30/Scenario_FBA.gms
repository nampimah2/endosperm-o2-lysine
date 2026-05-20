*************************************************************
* Scenario Testing FBA for Lysine Enhancement
* Tests the impact of gene knockouts, overexpression, and
* o2-bound mimicry on lysine production in WT models
*************************************************************
* Run this from each Wild_Type/B{DAP}/ directory
*************************************************************
$INLINECOM /*  */

OPTIONS
        decimals = 8
        lp = cplex
        limrow = 0
        limcol = 0
        solprint = off
        sysout = off
;

*********Defining Sets**************************************
SETS
        i                       set of metabolites
$include "metabolites.txt"

        j                       set of reactions
$include "reactions.txt"

        sc                      scenario index /s1*s30/
;
*************************************************************

***********Defining Parameters*******************************
PARAMETERS
        S(i,j)                  stoichiometric matrix
$include "sij.txt"

        v_max(j)                maximum flux of v(j)
$include "v_max.txt"

        v_min(j)                minimum flux of v(j)
$include "v_min.txt"

* Save original bounds for restoration between scenarios
        v_max_orig(j)           original maximum flux
        v_min_orig(j)           original minimum flux

* Results storage
        res_biomass(sc)         biomass for each scenario
        res_lysine(sc)          R00451 flux for each scenario
        res_dhdps(sc)           R02292 flux for each scenario
        res_lysExch(sc)         Exchange_C00047 flux for each scenario
        res_status(sc)          solver status for each scenario
;

* Store originals
v_max_orig(j) = v_max(j);
v_min_orig(j) = v_min(j);

**************************************************************

*********Defining Equations***********************************
EQUATIONS
        objective               objective function
        mass_balance1(i)        steady state mass balance
        lower_bound(j)          lower bounds on reactions
        upper_bound(j)          upper bounds on reactions
;
**************************************************************

*********Defining Variables***********************************
FREE VARIABLES
        v(j)                    reaction flux
        obj                     objective value
;
****************************************************************

***************Defining Model***********************************
objective..             obj =e= v('Seed_Biomass[K]');
mass_balance1(i)..      sum(j, S(i,j) * v(j)) =e= 0;
lower_bound(j)..        v_min(j) =l= v(j);
upper_bound(j)..        v(j) =l= v_max(j);

Model endosperm_fba /all/;
endosperm_fba.optfile = 1;
endosperm_fba.holdfixed = 1;
******************************************************************

*************************************************************
* MACRO: Restore all bounds to original values
*************************************************************
$macro RESTORE_BOUNDS v_max(j) = v_max_orig(j); v_min(j) = v_min_orig(j);

*************************************************************
* MACRO: Record results for scenario sc_idx
*************************************************************
$macro RECORD_RESULTS(sc_idx) \
  res_biomass(sc_idx) = obj.l; \
  res_status(sc_idx) = endosperm_fba.modelstat; \
  res_lysine(sc_idx)  = v.l('R00451[K,p]'); \
  res_dhdps(sc_idx)   = v.l('R02292[K,p]'); \
  res_lysExch(sc_idx) = v.l('Exchange_C00047[K,L]');


*************************************************************
* === SCENARIO 1: BASELINE (no modifications) ===
*************************************************************
Solve endosperm_fba using lp maximizing obj;
RECORD_RESULTS('s1')

SCALAR baseline_lysine;
baseline_lysine = v.l('R00451[K,p]');

*************************************************************
* === SCENARIO 2: KO R00282[K,c] - Acetyl-CoA Carboxylase ===
* Blocks entry into fatty acid biosynthesis
*************************************************************
RESTORE_BOUNDS
v_max('R00282[K,c]') = 0;
v_min('R00282[K,c]') = 0;
Solve endosperm_fba using lp maximizing obj;
RECORD_RESULTS('s2')

*************************************************************
* === SCENARIO 3: KO R00216[K,p] + MR00216[K,p] ===
* Pyruvate dehydrogenase (plastidial) - both reversible & irreversible
*************************************************************
RESTORE_BOUNDS
v_max('R00216[K,p]') = 0;
v_min('R00216[K,p]') = 0;
v_max('MR00216[K,p]') = 0;
v_min('MR00216[K,p]') = 0;
Solve endosperm_fba using lp maximizing obj;
RECORD_RESULTS('s3')

*************************************************************
* === SCENARIO 4: KO Exe2[K] - Lipid export ===
*************************************************************
RESTORE_BOUNDS
v_max('Exe2[K]') = 0;
v_min('Exe2[K]') = 0;
Solve endosperm_fba using lp maximizing obj;
RECORD_RESULTS('s4')

*************************************************************
* === SCENARIO 5: KO Trans2[K] - Lipid transport ===
*************************************************************
RESTORE_BOUNDS
v_max('Trans2[K]') = 0;
v_min('Trans2[K]') = 0;
Solve endosperm_fba using lp maximizing obj;
RECORD_RESULTS('s5')

*************************************************************
* === SCENARIO 6: KO cpTransport_C00007[K] - O2 transport ===
*************************************************************
RESTORE_BOUNDS
v_max('cpTransport_C00007[K]') = 0;
v_min('cpTransport_C00007[K]') = 0;
Solve endosperm_fba using lp maximizing obj;
RECORD_RESULTS('s6')

*************************************************************
* === SCENARIO 7: OE R00883[K,c] - Aspartate kinase (2x vmax) ===
*************************************************************
RESTORE_BOUNDS
v_max('R00883[K,c]') = 2 * v_max_orig('R00883[K,c]');
Solve endosperm_fba using lp maximizing obj;
RECORD_RESULTS('s7')

*************************************************************
* === SCENARIO 8: OE R00512[K,c] - Aspartate semialdehyde DH (2x) ===
*************************************************************
RESTORE_BOUNDS
v_max('R00512[K,c]') = 2 * v_max_orig('R00512[K,c]');
Solve endosperm_fba using lp maximizing obj;
RECORD_RESULTS('s8')

*************************************************************
* === SCENARIO 9: OE R01015[K,c] - Dihydrodipicolinate reductase (2x) ===
*************************************************************
RESTORE_BOUNDS
v_max('R01015[K,c]') = 2 * v_max_orig('R01015[K,c]');
Solve endosperm_fba using lp maximizing obj;
RECORD_RESULTS('s9')

*************************************************************
* === SCENARIO 10: OE R00667[K,m] - Homoserine DH (2x) ===
*************************************************************
RESTORE_BOUNDS
v_max('R00667[K,m]') = 2 * v_max_orig('R00667[K,m]');
Solve endosperm_fba using lp maximizing obj;
RECORD_RESULTS('s10')

*************************************************************
* === SCENARIO 11: OE R01280[K,p] - Diaminopimelate epimerase (2x) ===
*************************************************************
RESTORE_BOUNDS
v_max('R01280[K,p]') = 2 * v_max_orig('R01280[K,p]');
Solve endosperm_fba using lp maximizing obj;
RECORD_RESULTS('s11')

*************************************************************
* === SCENARIO 12: O2-Mimicry MR00422[K,p] ===
* Replace WT vmax with O2 mutant vmax
* (O2 values loaded via parameter below)
*************************************************************
RESTORE_BOUNDS
* O2 mimicry value will be set by the shell script via --define
* Default: reduce by 15% as representative O2 effect
$if not set O2_MR00422_VMAX $set O2_MR00422_VMAX 0
SCALAR o2_mr00422_val / %O2_MR00422_VMAX% /;
$if "%O2_MR00422_VMAX%" == "0" o2_mr00422_val = 0.85 * v_max_orig('MR00422[K,p]');
v_max('MR00422[K,p]') = o2_mr00422_val;
Solve endosperm_fba using lp maximizing obj;
RECORD_RESULTS('s12')

*************************************************************
* === SCENARIO 13: O2-Mimicry MR00732[K,m] ===
*************************************************************
RESTORE_BOUNDS
$if not set O2_MR00732_VMAX $set O2_MR00732_VMAX 0
SCALAR o2_mr00732_val / %O2_MR00732_VMAX% /;
$if "%O2_MR00732_VMAX%" == "0" o2_mr00732_val = 1.1 * v_max_orig('MR00732[K,m]');
v_max('MR00732[K,m]') = o2_mr00732_val;
Solve endosperm_fba using lp maximizing obj;
RECORD_RESULTS('s13')

*************************************************************
* === SCENARIO 14: O2-Mimicry R00342[K,p] ===
* Aspartate aminotransferase - increase to O2 level
*************************************************************
RESTORE_BOUNDS
$if not set O2_R00342_VMAX $set O2_R00342_VMAX 0
SCALAR o2_r00342_vmax_val / %O2_R00342_VMAX% /;
$if "%O2_R00342_VMAX%" == "0" o2_r00342_vmax_val = 1.5 * v_max_orig('R00342[K,p]');
v_max('R00342[K,p]') = o2_r00342_vmax_val;
$if not set O2_R00342_VMIN $set O2_R00342_VMIN 0
SCALAR o2_r00342_vmin_val / %O2_R00342_VMIN% /;
$if "%O2_R00342_VMIN%" == "0" o2_r00342_vmin_val = -1.5 * v_max_orig('R00342[K,p]');
v_min('R00342[K,p]') = o2_r00342_vmin_val;
Solve endosperm_fba using lp maximizing obj;
RECORD_RESULTS('s14')

*************************************************************
* === SCENARIO 15: COMBO - KO R00282 + OE R00883 ===
* Block FA entry + boost aspartate kinase
*************************************************************
RESTORE_BOUNDS
v_max('R00282[K,c]') = 0;
v_min('R00282[K,c]') = 0;
v_max('R00883[K,c]') = 2 * v_max_orig('R00883[K,c]');
Solve endosperm_fba using lp maximizing obj;
RECORD_RESULTS('s15')

*************************************************************
* === SCENARIO 16: COMBO - KO R00282 + OE R00512 + OE R01015 ===
* Block FA entry + boost lysine pathway
*************************************************************
RESTORE_BOUNDS
v_max('R00282[K,c]') = 0;
v_min('R00282[K,c]') = 0;
v_max('R00512[K,c]') = 2 * v_max_orig('R00512[K,c]');
v_max('R01015[K,c]') = 2 * v_max_orig('R01015[K,c]');
Solve endosperm_fba using lp maximizing obj;
RECORD_RESULTS('s16')

*************************************************************
* === SCENARIO 17: COMBO - KO R00216[K,p] + MR00216[K,p] + OE R00883 ===
* Block plastidial pyruvate DH + boost aspartate kinase
*************************************************************
RESTORE_BOUNDS
v_max('R00216[K,p]') = 0;
v_min('R00216[K,p]') = 0;
v_max('MR00216[K,p]') = 0;
v_min('MR00216[K,p]') = 0;
v_max('R00883[K,c]') = 2 * v_max_orig('R00883[K,c]');
Solve endosperm_fba using lp maximizing obj;
RECORD_RESULTS('s17')

*************************************************************
* === SCENARIO 18: COMBO - KO Exe2 + KO Trans2 (block lipid export/transport) ===
*************************************************************
RESTORE_BOUNDS
v_max('Exe2[K]') = 0;
v_min('Exe2[K]') = 0;
v_max('Trans2[K]') = 0;
v_min('Trans2[K]') = 0;
Solve endosperm_fba using lp maximizing obj;
RECORD_RESULTS('s18')

*************************************************************
* === SCENARIO 19: COMBO - All 3 mimicry targets ===
* MR00422 + MR00732 + R00342 set to O2 levels simultaneously
*************************************************************
RESTORE_BOUNDS
v_max('MR00422[K,p]') = o2_mr00422_val;
v_max('MR00732[K,m]') = o2_mr00732_val;
v_max('R00342[K,p]') = o2_r00342_vmax_val;
v_min('R00342[K,p]') = o2_r00342_vmin_val;
Solve endosperm_fba using lp maximizing obj;
RECORD_RESULTS('s19')

*************************************************************
* === SCENARIO 20: COMBO - KO R00282 + all mimicry ===
* Block FA entry + mimic O2 bound changes
*************************************************************
RESTORE_BOUNDS
v_max('R00282[K,c]') = 0;
v_min('R00282[K,c]') = 0;
v_max('MR00422[K,p]') = o2_mr00422_val;
v_max('MR00732[K,m]') = o2_mr00732_val;
v_max('R00342[K,p]') = o2_r00342_vmax_val;
v_min('R00342[K,p]') = o2_r00342_vmin_val;
Solve endosperm_fba using lp maximizing obj;
RECORD_RESULTS('s20')

*************************************************************
* === SCENARIO 21: COMBO - All KOs (R00282 + Exe2 + Trans2) ===
* Complete FA pathway blockade
*************************************************************
RESTORE_BOUNDS
v_max('R00282[K,c]') = 0;
v_min('R00282[K,c]') = 0;
v_max('Exe2[K]') = 0;
v_min('Exe2[K]') = 0;
v_max('Trans2[K]') = 0;
v_min('Trans2[K]') = 0;
Solve endosperm_fba using lp maximizing obj;
RECORD_RESULTS('s21')

*************************************************************
* === SCENARIO 22: COMBO - All KOs + All OE ===
* Complete FA blockade + all lysine pathway OE
*************************************************************
RESTORE_BOUNDS
v_max('R00282[K,c]') = 0;
v_min('R00282[K,c]') = 0;
v_max('Exe2[K]') = 0;
v_min('Exe2[K]') = 0;
v_max('Trans2[K]') = 0;
v_min('Trans2[K]') = 0;
v_max('R00883[K,c]') = 2 * v_max_orig('R00883[K,c]');
v_max('R00512[K,c]') = 2 * v_max_orig('R00512[K,c]');
v_max('R01015[K,c]') = 2 * v_max_orig('R01015[K,c]');
v_max('R00667[K,m]') = 2 * v_max_orig('R00667[K,m]');
v_max('R01280[K,p]') = 2 * v_max_orig('R01280[K,p]');
Solve endosperm_fba using lp maximizing obj;
RECORD_RESULTS('s22')

*************************************************************
* === SCENARIO 23: OE All lysine pathway (R00883+R00512+R01015+R01280) ===
*************************************************************
RESTORE_BOUNDS
v_max('R00883[K,c]') = 2 * v_max_orig('R00883[K,c]');
v_max('R00512[K,c]') = 2 * v_max_orig('R00512[K,c]');
v_max('R01015[K,c]') = 2 * v_max_orig('R01015[K,c]');
v_max('R01280[K,p]') = 2 * v_max_orig('R01280[K,p]');
Solve endosperm_fba using lp maximizing obj;
RECORD_RESULTS('s23')

*************************************************************
* === SCENARIO 24: COMBO - KO R00282 + OE all lysine pathway ===
*************************************************************
RESTORE_BOUNDS
v_max('R00282[K,c]') = 0;
v_min('R00282[K,c]') = 0;
v_max('R00883[K,c]') = 2 * v_max_orig('R00883[K,c]');
v_max('R00512[K,c]') = 2 * v_max_orig('R00512[K,c]');
v_max('R01015[K,c]') = 2 * v_max_orig('R01015[K,c]');
v_max('R01280[K,p]') = 2 * v_max_orig('R01280[K,p]');
Solve endosperm_fba using lp maximizing obj;
RECORD_RESULTS('s24')

*************************************************************
* === SCENARIO 25: GRAND COMBO - All KOs + All OE + All mimicry ===
*************************************************************
RESTORE_BOUNDS
v_max('R00282[K,c]') = 0;
v_min('R00282[K,c]') = 0;
v_max('Exe2[K]') = 0;
v_min('Exe2[K]') = 0;
v_max('Trans2[K]') = 0;
v_min('Trans2[K]') = 0;
v_max('R00883[K,c]') = 2 * v_max_orig('R00883[K,c]');
v_max('R00512[K,c]') = 2 * v_max_orig('R00512[K,c]');
v_max('R01015[K,c]') = 2 * v_max_orig('R01015[K,c]');
v_max('R00667[K,m]') = 2 * v_max_orig('R00667[K,m]');
v_max('R01280[K,p]') = 2 * v_max_orig('R01280[K,p]');
v_max('MR00422[K,p]') = o2_mr00422_val;
v_max('MR00732[K,m]') = o2_mr00732_val;
v_max('R00342[K,p]') = o2_r00342_vmax_val;
v_min('R00342[K,p]') = o2_r00342_vmin_val;
Solve endosperm_fba using lp maximizing obj;
RECORD_RESULTS('s25')

*************************************************************
* === SCENARIO 26: KO R00282 + OE R00883 + Mimicry R00342 ===
* FA block + boost aspartate kinase + boost aspartate aminotransferase
*************************************************************
RESTORE_BOUNDS
v_max('R00282[K,c]') = 0;
v_min('R00282[K,c]') = 0;
v_max('R00883[K,c]') = 2 * v_max_orig('R00883[K,c]');
v_max('R00342[K,p]') = o2_r00342_vmax_val;
v_min('R00342[K,p]') = o2_r00342_vmin_val;
Solve endosperm_fba using lp maximizing obj;
RECORD_RESULTS('s26')

*************************************************************
* === SCENARIO 27: OE R00883 + OE R00667 (aspartate kinase + homoserine DH) ===
*************************************************************
RESTORE_BOUNDS
v_max('R00883[K,c]') = 2 * v_max_orig('R00883[K,c]');
v_max('R00667[K,m]') = 2 * v_max_orig('R00667[K,m]');
Solve endosperm_fba using lp maximizing obj;
RECORD_RESULTS('s27')

*************************************************************
* === SCENARIO 28: KO R00216+MR00216 + KO R00282 (all plastidial FA entry) ===
*************************************************************
RESTORE_BOUNDS
v_max('R00216[K,p]') = 0;
v_min('R00216[K,p]') = 0;
v_max('MR00216[K,p]') = 0;
v_min('MR00216[K,p]') = 0;
v_max('R00282[K,c]') = 0;
v_min('R00282[K,c]') = 0;
Solve endosperm_fba using lp maximizing obj;
RECORD_RESULTS('s28')

*************************************************************
* OUTPUT RESULTS TO CSV FILE
*************************************************************
FILE RESULTS /scenario_results.csv/;
PUT RESULTS;

* Header
PUT "Scenario,Description,Biomass,R00451_Lysine,R02292_DHDPS,Exchange_C00047,ModelStat,Lys_FoldChange" /;

* Scenario names
PUT "s1,Baseline,"
    res_biomass('s1'):0:8 ","
    res_lysine('s1'):0:8 ","
    res_dhdps('s1'):0:8 ","
    res_lysExch('s1'):0:8 ","
    res_status('s1'):0:0 ","
    "1.00000000" /;

PUT "s2,KO_R00282_ACC,"
    res_biomass('s2'):0:8 ","
    res_lysine('s2'):0:8 ","
    res_dhdps('s2'):0:8 ","
    res_lysExch('s2'):0:8 ","
    res_status('s2'):0:0 ",";
if(baseline_lysine > 0,
    PUT (res_lysine('s2')/baseline_lysine):0:8 /;
else
    PUT "Inf" /;
);

PUT "s3,KO_R00216_PDH_plastid,"
    res_biomass('s3'):0:8 ","
    res_lysine('s3'):0:8 ","
    res_dhdps('s3'):0:8 ","
    res_lysExch('s3'):0:8 ","
    res_status('s3'):0:0 ",";
if(baseline_lysine > 0,
    PUT (res_lysine('s3')/baseline_lysine):0:8 /;
else
    PUT "Inf" /;
);

PUT "s4,KO_Exe2_lipid_export,"
    res_biomass('s4'):0:8 ","
    res_lysine('s4'):0:8 ","
    res_dhdps('s4'):0:8 ","
    res_lysExch('s4'):0:8 ","
    res_status('s4'):0:0 ",";
if(baseline_lysine > 0,
    PUT (res_lysine('s4')/baseline_lysine):0:8 /;
else
    PUT "Inf" /;
);

PUT "s5,KO_Trans2_lipid_transport,"
    res_biomass('s5'):0:8 ","
    res_lysine('s5'):0:8 ","
    res_dhdps('s5'):0:8 ","
    res_lysExch('s5'):0:8 ","
    res_status('s5'):0:0 ",";
if(baseline_lysine > 0,
    PUT (res_lysine('s5')/baseline_lysine):0:8 /;
else
    PUT "Inf" /;
);

PUT "s6,KO_cpTransO2_transport,"
    res_biomass('s6'):0:8 ","
    res_lysine('s6'):0:8 ","
    res_dhdps('s6'):0:8 ","
    res_lysExch('s6'):0:8 ","
    res_status('s6'):0:0 ",";
if(baseline_lysine > 0,
    PUT (res_lysine('s6')/baseline_lysine):0:8 /;
else
    PUT "Inf" /;
);

PUT "s7,OE_R00883_AspKinase,"
    res_biomass('s7'):0:8 ","
    res_lysine('s7'):0:8 ","
    res_dhdps('s7'):0:8 ","
    res_lysExch('s7'):0:8 ","
    res_status('s7'):0:0 ",";
if(baseline_lysine > 0,
    PUT (res_lysine('s7')/baseline_lysine):0:8 /;
else
    PUT "Inf" /;
);

PUT "s8,OE_R00512_AspSemiDH,"
    res_biomass('s8'):0:8 ","
    res_lysine('s8'):0:8 ","
    res_dhdps('s8'):0:8 ","
    res_lysExch('s8'):0:8 ","
    res_status('s8'):0:0 ",";
if(baseline_lysine > 0,
    PUT (res_lysine('s8')/baseline_lysine):0:8 /;
else
    PUT "Inf" /;
);

PUT "s9,OE_R01015_DHDPR,"
    res_biomass('s9'):0:8 ","
    res_lysine('s9'):0:8 ","
    res_dhdps('s9'):0:8 ","
    res_lysExch('s9'):0:8 ","
    res_status('s9'):0:0 ",";
if(baseline_lysine > 0,
    PUT (res_lysine('s9')/baseline_lysine):0:8 /;
else
    PUT "Inf" /;
);

PUT "s10,OE_R00667_HomoserDH,"
    res_biomass('s10'):0:8 ","
    res_lysine('s10'):0:8 ","
    res_dhdps('s10'):0:8 ","
    res_lysExch('s10'):0:8 ","
    res_status('s10'):0:0 ",";
if(baseline_lysine > 0,
    PUT (res_lysine('s10')/baseline_lysine):0:8 /;
else
    PUT "Inf" /;
);

PUT "s11,OE_R01280_DAPepimerase,"
    res_biomass('s11'):0:8 ","
    res_lysine('s11'):0:8 ","
    res_dhdps('s11'):0:8 ","
    res_lysExch('s11'):0:8 ","
    res_status('s11'):0:0 ",";
if(baseline_lysine > 0,
    PUT (res_lysine('s11')/baseline_lysine):0:8 /;
else
    PUT "Inf" /;
);

PUT "s12,Mimicry_MR00422,"
    res_biomass('s12'):0:8 ","
    res_lysine('s12'):0:8 ","
    res_dhdps('s12'):0:8 ","
    res_lysExch('s12'):0:8 ","
    res_status('s12'):0:0 ",";
if(baseline_lysine > 0,
    PUT (res_lysine('s12')/baseline_lysine):0:8 /;
else
    PUT "Inf" /;
);

PUT "s13,Mimicry_MR00732,"
    res_biomass('s13'):0:8 ","
    res_lysine('s13'):0:8 ","
    res_dhdps('s13'):0:8 ","
    res_lysExch('s13'):0:8 ","
    res_status('s13'):0:0 ",";
if(baseline_lysine > 0,
    PUT (res_lysine('s13')/baseline_lysine):0:8 /;
else
    PUT "Inf" /;
);

PUT "s14,Mimicry_R00342,"
    res_biomass('s14'):0:8 ","
    res_lysine('s14'):0:8 ","
    res_dhdps('s14'):0:8 ","
    res_lysExch('s14'):0:8 ","
    res_status('s14'):0:0 ",";
if(baseline_lysine > 0,
    PUT (res_lysine('s14')/baseline_lysine):0:8 /;
else
    PUT "Inf" /;
);

PUT "s15,COMBO_KO_R00282+OE_R00883,"
    res_biomass('s15'):0:8 ","
    res_lysine('s15'):0:8 ","
    res_dhdps('s15'):0:8 ","
    res_lysExch('s15'):0:8 ","
    res_status('s15'):0:0 ",";
if(baseline_lysine > 0,
    PUT (res_lysine('s15')/baseline_lysine):0:8 /;
else
    PUT "Inf" /;
);

PUT "s16,COMBO_KO_R00282+OE_R00512+R01015,"
    res_biomass('s16'):0:8 ","
    res_lysine('s16'):0:8 ","
    res_dhdps('s16'):0:8 ","
    res_lysExch('s16'):0:8 ","
    res_status('s16'):0:0 ",";
if(baseline_lysine > 0,
    PUT (res_lysine('s16')/baseline_lysine):0:8 /;
else
    PUT "Inf" /;
);

PUT "s17,COMBO_KO_PDH+OE_R00883,"
    res_biomass('s17'):0:8 ","
    res_lysine('s17'):0:8 ","
    res_dhdps('s17'):0:8 ","
    res_lysExch('s17'):0:8 ","
    res_status('s17'):0:0 ",";
if(baseline_lysine > 0,
    PUT (res_lysine('s17')/baseline_lysine):0:8 /;
else
    PUT "Inf" /;
);

PUT "s18,COMBO_KO_Exe2+Trans2,"
    res_biomass('s18'):0:8 ","
    res_lysine('s18'):0:8 ","
    res_dhdps('s18'):0:8 ","
    res_lysExch('s18'):0:8 ","
    res_status('s18'):0:0 ",";
if(baseline_lysine > 0,
    PUT (res_lysine('s18')/baseline_lysine):0:8 /;
else
    PUT "Inf" /;
);

PUT "s19,COMBO_All3_mimicry,"
    res_biomass('s19'):0:8 ","
    res_lysine('s19'):0:8 ","
    res_dhdps('s19'):0:8 ","
    res_lysExch('s19'):0:8 ","
    res_status('s19'):0:0 ",";
if(baseline_lysine > 0,
    PUT (res_lysine('s19')/baseline_lysine):0:8 /;
else
    PUT "Inf" /;
);

PUT "s20,COMBO_KO_R00282+All_mimicry,"
    res_biomass('s20'):0:8 ","
    res_lysine('s20'):0:8 ","
    res_dhdps('s20'):0:8 ","
    res_lysExch('s20'):0:8 ","
    res_status('s20'):0:0 ",";
if(baseline_lysine > 0,
    PUT (res_lysine('s20')/baseline_lysine):0:8 /;
else
    PUT "Inf" /;
);

PUT "s21,COMBO_AllKOs_FA_block,"
    res_biomass('s21'):0:8 ","
    res_lysine('s21'):0:8 ","
    res_dhdps('s21'):0:8 ","
    res_lysExch('s21'):0:8 ","
    res_status('s21'):0:0 ",";
if(baseline_lysine > 0,
    PUT (res_lysine('s21')/baseline_lysine):0:8 /;
else
    PUT "Inf" /;
);

PUT "s22,COMBO_AllKOs+AllOE,"
    res_biomass('s22'):0:8 ","
    res_lysine('s22'):0:8 ","
    res_dhdps('s22'):0:8 ","
    res_lysExch('s22'):0:8 ","
    res_status('s22'):0:0 ",";
if(baseline_lysine > 0,
    PUT (res_lysine('s22')/baseline_lysine):0:8 /;
else
    PUT "Inf" /;
);

PUT "s23,OE_All_lysine_pathway,"
    res_biomass('s23'):0:8 ","
    res_lysine('s23'):0:8 ","
    res_dhdps('s23'):0:8 ","
    res_lysExch('s23'):0:8 ","
    res_status('s23'):0:0 ",";
if(baseline_lysine > 0,
    PUT (res_lysine('s23')/baseline_lysine):0:8 /;
else
    PUT "Inf" /;
);

PUT "s24,COMBO_KO_R00282+OE_AllLys,"
    res_biomass('s24'):0:8 ","
    res_lysine('s24'):0:8 ","
    res_dhdps('s24'):0:8 ","
    res_lysExch('s24'):0:8 ","
    res_status('s24'):0:0 ",";
if(baseline_lysine > 0,
    PUT (res_lysine('s24')/baseline_lysine):0:8 /;
else
    PUT "Inf" /;
);

PUT "s25,GRAND_COMBO_AllKO+AllOE+AllMim,"
    res_biomass('s25'):0:8 ","
    res_lysine('s25'):0:8 ","
    res_dhdps('s25'):0:8 ","
    res_lysExch('s25'):0:8 ","
    res_status('s25'):0:0 ",";
if(baseline_lysine > 0,
    PUT (res_lysine('s25')/baseline_lysine):0:8 /;
else
    PUT "Inf" /;
);

PUT "s26,COMBO_KO_R00282+OE_R00883+Mim_R00342,"
    res_biomass('s26'):0:8 ","
    res_lysine('s26'):0:8 ","
    res_dhdps('s26'):0:8 ","
    res_lysExch('s26'):0:8 ","
    res_status('s26'):0:0 ",";
if(baseline_lysine > 0,
    PUT (res_lysine('s26')/baseline_lysine):0:8 /;
else
    PUT "Inf" /;
);

PUT "s27,OE_R00883+R00667,"
    res_biomass('s27'):0:8 ","
    res_lysine('s27'):0:8 ","
    res_dhdps('s27'):0:8 ","
    res_lysExch('s27'):0:8 ","
    res_status('s27'):0:0 ",";
if(baseline_lysine > 0,
    PUT (res_lysine('s27')/baseline_lysine):0:8 /;
else
    PUT "Inf" /;
);

PUT "s28,KO_PDH+ACC_all_FA_entry,"
    res_biomass('s28'):0:8 ","
    res_lysine('s28'):0:8 ","
    res_dhdps('s28'):0:8 ","
    res_lysExch('s28'):0:8 ","
    res_status('s28'):0:0 ",";
if(baseline_lysine > 0,
    PUT (res_lysine('s28')/baseline_lysine):0:8 /;
else
    PUT "Inf" /;
);

PUTCLOSE;

DISPLAY res_biomass, res_lysine, res_dhdps, res_lysExch, res_status;
