*************************************************************
* FSEOF: Flux Scanning based on Enforced Objective Flux
* Target: R00451[K,p] (lysine biosynthesis)
* Steps: 10, Biomass fraction: 0.9
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
    s   steps /s0*s10/
;

PARAMETERS
    S(i,j)   stoichiometric matrix
$include "sij.txt"
    v_max_orig(j) maximum flux
$include "v_max.txt"
    v_min_orig(j) minimum flux
$include "v_min.txt"

    biomass_floor
    max_r00451
    r00451_target(s)
    flux_profile(s,j)
;

EQUATIONS
    objective
    mass_balance1(i)
    lower_bound(j)
    upper_bound(j)
    biomass_constraint
    r00451_enforce
;

FREE VARIABLES
    v(j)
    obj
;

SCALAR r00451_min /0/;

objective..         obj =e= v('Seed_Biomass[K]');
mass_balance1(i)..  sum(j, S(i,j) * v(j)) =e= 0;
lower_bound(j)..    v_min_orig(j) =l= v(j);
upper_bound(j)..    v(j) =l= v_max_orig(j);
biomass_constraint.. v('Seed_Biomass[K]') =g= biomass_floor;
r00451_enforce..    v('R00451[K,p]') =g= r00451_min;

* Step 1: Solve for max biomass
Model fba_base /objective, mass_balance1, lower_bound, upper_bound/;
fba_base.optfile = 1;
Solve fba_base using lp maximizing obj;

SCALAR base_biomass;
SCALAR base_r00451;
base_biomass = v.l('Seed_Biomass[K]');
base_r00451  = v.l('R00451[K,p]');
biomass_floor = 0.9 * base_biomass;

* Step 2: Max R00451 subject to biomass floor
Model fba_max_r00451 /mass_balance1, lower_bound, upper_bound, biomass_constraint/;
fba_max_r00451.optfile = 1;

VARIABLE obj2;
EQUATION obj_r00451;
obj_r00451.. obj2 =e= v('R00451[K,p]');
Model fba_max_lys /obj_r00451, mass_balance1, lower_bound, upper_bound, biomass_constraint/;
fba_max_lys.optfile = 1;
Solve fba_max_lys using lp maximizing obj2;
max_r00451 = v.l('R00451[K,p]');

* Step 3: Scan
SCALAR step_size;
step_size = (max_r00451 - base_r00451) / 10;

Model fba_fseof /objective, mass_balance1, lower_bound, upper_bound, r00451_enforce/;
fba_fseof.optfile = 1;

loop(s,
    r00451_min = base_r00451 + ord(s) * step_size;
    r00451_target(s) = r00451_min;
    
    Solve fba_fseof using lp maximizing obj;
    
    if(fba_fseof.modelstat = 1 or fba_fseof.modelstat = 2,
        flux_profile(s,j) = v.l(j);
    );
);

* Output
FILE RESULTS /'results_fseof.txt'/;
PUT RESULTS;
PUT "BASE_BIOMASS ", base_biomass:0:10 /;
PUT "BASE_R00451 ", base_r00451:0:10 /;
PUT "MAX_R00451 ", max_r00451:0:10 /;
PUT "STEP_SIZE ", step_size:0:10 /;
loop(s,
    PUT "STEP ", ord(s):0:0, " TARGET ", r00451_target(s):0:10 /;
    loop(j$(abs(flux_profile(s,j)) > 1e-10),
        PUT "  ", j.tl:0:40, " ", flux_profile(s,j):0:10 /;
    );
);
PUTCLOSE;
