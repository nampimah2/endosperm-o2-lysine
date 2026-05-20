*************************************************************
* Batch Overexpression Screen for Lysine Target ID
* Factor = 2.0x v_max
*************************************************************
$INLINECOM /*  */

OPTIONS
    decimals = 8
    lp = cplex
    solprint = off
    limrow = 0
    limcol = 0
    solvelink = 5
;

SETS
    i   set of metabolites
$include "metabolites.txt"
    j   set of reactions
$include "reactions.txt"
;

PARAMETERS
    S(i,j)   stoichiometric matrix
$include "sij.txt"
    v_max_orig(j) maximum flux
$include "v_max.txt"
    v_min_orig(j) minimum flux
$include "v_min.txt"

    v_max_mod(j) modified max bounds
    v_min_mod(j) modified min bounds
    result_biomass(j)    biomass when j is overexpressed
    result_r00451(j)     R00451 flux when j is overexpressed
    base_biomass         baseline biomass
    base_r00451          baseline R00451
;

EQUATIONS
    objective
    mass_balance1(i)
    lower_bound(j)
    upper_bound(j)
;

FREE VARIABLES
    v(j)
    obj
;

objective..         obj =e= v('Seed_Biomass[K]');
mass_balance1(i)..  sum(j, S(i,j) * v(j)) =e= 0;
lower_bound(j)..    v_min_mod(j) =l= v(j);
upper_bound(j)..    v(j) =l= v_max_mod(j);

Model fba_screen /all/;
fba_screen.optfile = 1;
fba_screen.holdfixed = 1;

* -------- Baseline solve --------
v_max_mod(j) = v_max_orig(j);
v_min_mod(j) = v_min_orig(j);
Solve fba_screen using lp maximizing obj;
base_biomass = v.l('Seed_Biomass[K]');
base_r00451  = v.l('R00451[K,p]');

* -------- Overexpression loop --------
ALIAS(j, jj);

loop(jj$(v_max_orig(jj) > 1e-12 or v_min_orig(jj) < -1e-12),
    v_max_mod(j) = v_max_orig(j);
    v_min_mod(j) = v_min_orig(j);
    
    v_max_mod(jj) = v_max_orig(jj) * 2.0;
    if(v_min_orig(jj) < 0,
        v_min_mod(jj) = v_min_orig(jj) * 2.0;
    );
    
    Solve fba_screen using lp maximizing obj;
    
    if(fba_screen.modelstat = 1 or fba_screen.modelstat = 2,
        result_biomass(jj) = v.l('Seed_Biomass[K]');
        result_r00451(jj)  = v.l('R00451[K,p]');
    else
        result_biomass(jj) = -999;
        result_r00451(jj)  = -999;
    );
);

* -------- Output Results --------
FILE RESULTS /'results_overexpression.txt'/;
PUT RESULTS;
PUT "BASELINE_BIOMASS ", base_biomass:0:10 /;
PUT "BASELINE_R00451 ", base_r00451:0:10 /;
PUT "REACTION BIOMASS R00451" /;
loop(jj$(result_r00451(jj) ne 0),
    put jj.tl:0:40, " ", result_biomass(jj):0:10, " ", result_r00451(jj):0:10 /;
);
PUTCLOSE;
