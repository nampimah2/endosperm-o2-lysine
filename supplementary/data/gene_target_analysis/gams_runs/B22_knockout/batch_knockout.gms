*************************************************************
* Batch Knockout Screen for Lysine Target ID
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
    result_biomass(j)
    result_r00451(j)
    base_biomass
    base_r00451
    base_flux(j)
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

Model fba_ko /all/;
fba_ko.optfile = 1;
fba_ko.holdfixed = 1;

* -------- Baseline solve --------
v_max_mod(j) = v_max_orig(j);
v_min_mod(j) = v_min_orig(j);
Solve fba_ko using lp maximizing obj;
base_biomass = v.l('Seed_Biomass[K]');
base_r00451  = v.l('R00451[K,p]');
base_flux(j) = v.l(j);

* -------- Knockout loop (only active reactions) --------
ALIAS(j, jj);

loop(jj$(abs(base_flux(jj)) > 1e-12
         and not sameas(jj, 'Seed_Biomass[K]')
         and not sameas(jj, 'R00451[K,p]')),
    v_max_mod(j) = v_max_orig(j);
    v_min_mod(j) = v_min_orig(j);
    
    v_max_mod(jj) = 0;
    v_min_mod(jj) = 0;
    
    Solve fba_ko using lp maximizing obj;
    
    if(fba_ko.modelstat = 1 or fba_ko.modelstat = 2,
        result_biomass(jj) = v.l('Seed_Biomass[K]');
        result_r00451(jj)  = v.l('R00451[K,p]');
    else
        result_biomass(jj) = -999;
        result_r00451(jj)  = -999;
    );
);

* -------- Output Results --------
FILE RESULTS /'results_knockout.txt'/;
PUT RESULTS;
PUT "BASELINE_BIOMASS ", base_biomass:0:10 /;
PUT "BASELINE_R00451 ", base_r00451:0:10 /;
PUT "REACTION BIOMASS R00451 ORIG_FLUX" /;
loop(jj$(result_r00451(jj) ne 0),
    put jj.tl:0:40, " ", result_biomass(jj):0:10, " ", result_r00451(jj):0:10, " ", base_flux(jj):0:10 /;
);
PUTCLOSE;
