$INLINECOM /*  */
OPTIONS
    decimals = 8
    lp = cplex
    limrow = 0
    limcol = 0
    solprint = off
    sysout = off
;

SETS
    i   set of metabolites
$include "metabolites.txt"
    j   set of reactions
$include "reactions.txt"
;

ALIAS(j, jj);

PARAMETERS
    S(i,j)   stoichiometric matrix
$include "sij.txt"
    v_max(j)   maximum flux
$include "v_max.txt"
    v_min(j)   minimum flux
$include "v_min.txt"
;

FREE VARIABLES
    v(j)    reaction flux
    obj     objective value
;

EQUATIONS
    objective
    mass_balance(i)
    lower_bound(j)
    upper_bound(j)
;

objective..              obj =e= v('R00451[K,p]');
mass_balance(i)..        sum(j, S(i,j) * v(j)) =e= 0;
lower_bound(j)..         v_min(j) =l= v(j);
upper_bound(j)..         v(j) =l= v_max(j);

Model fba /all/;
fba.optfile = 1;
fba.holdfixed = 1;

Solve fba using lp maximizing obj;

SCALAR target_max_flux;
target_max_flux = obj.l;

FILE RESULTS /fva_results.txt/;
PUT RESULTS;
PUT "TARGET_MAX_FLUX " target_max_flux:0:10 /;

v_min('R00451[K,p]') = target_max_flux * 0.01;

SCALAR target_lb;
target_lb = target_max_flux * 0.01;
PUT "TARGET_LB " target_lb:0:10 /;

PARAMETER obj_coeff(j);
obj_coeff(j) = 0;

EQUATIONS
    fva_objective
;
fva_objective..   obj =e= sum(j, obj_coeff(j) * v(j));

Model fva_model /mass_balance, lower_bound, upper_bound, fva_objective/;
fva_model.optfile = 1;
fva_model.holdfixed = 1;

PARAMETERS
    fva_max(j)
    fva_min(j)
;

loop(jj,
    obj_coeff(jj) = 1;

    Solve fva_model using lp maximizing obj;
    if(fva_model.modelstat le 2,
        fva_max(jj) = obj.l;
    else
        fva_max(jj) = 0;
    );

    Solve fva_model using lp minimizing obj;
    if(fva_model.modelstat le 2,
        fva_min(jj) = obj.l;
    else
        fva_min(jj) = 0;
    );

    obj_coeff(jj) = 0;
);

loop(j,
    PUT j.tl:0:50, " ", fva_min(j):0:10, " ", fva_max(j):0:10 /;
);

PUTCLOSE;
