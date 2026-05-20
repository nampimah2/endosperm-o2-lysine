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
    i set of metabolites
$include "metabolites.txt"
    j set of reactions
$include "reactions.txt"
    jcheck(j) /
        'R00355[K,p]'
        'R00480[K,p]'
        'R02291[K,p]'
        'R02292[K,p]'
        'R02735[K,p]'
        'Exe4[K]'
        'R04198[K,p]'
        'cpTransport_C00011[K]'
        'cpTransport_C00712[K]'
        'Trans4[K]'
    /
;

ALIAS(jcheck, jc);

PARAMETERS
    S(i,j) stoichiometric matrix
$include "sij.txt"
    v_max(j) maximum flux
$include "v_max.txt"
    v_min(j) minimum flux
$include "v_min.txt"
    v_min_orig(j)
;

v_min_orig(j) = v_min(j);

FREE VARIABLES
    v(j) reaction flux
    obj  objective value
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

PARAMETERS
    rev_min(j)
    rev_max(j)
;

FILE RESULTS /reverse_coupling_results.txt/;
PUT RESULTS;

loop(jc,
    v_min(j) = v_min_orig(j);
    v_min(jc) = max(0.001, v_min_orig(jc));

    Solve fba using lp maximizing obj;
    if(fba.modelstat le 2,
        rev_max(jc) = obj.l;
    else
        rev_max(jc) = -999;
    );

    Solve fba using lp minimizing obj;
    if(fba.modelstat le 2,
        rev_min(jc) = obj.l;
    else
        rev_min(jc) = -999;
    );

    v_min(jc) = v_min_orig(jc);

    PUT jc.tl:0:50, " ", rev_min(jc):0:10, " ", rev_max(jc):0:10 /;
);

PUTCLOSE;
