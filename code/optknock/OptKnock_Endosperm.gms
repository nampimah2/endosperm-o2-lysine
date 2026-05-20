*************************************************************
* OptKnock: Bilevel Optimization for Lysine Overproduction
* in Maize Endosperm Metabolic Model
*
* Based on: Burgard AP, Pharkya P, Maranas CD (2003)
*   "OptKnock: A bilevel programming framework for
*    identifying gene knockout strategies for microbial
*    strain optimization"
*   Biotechnology and Bioengineering 84(6):647-657
*
* Model: Maize endosperm FBA model (kernel development)
* Outer objective: Maximize lysine biosynthesis (R00451[K,p])
* Inner objective: Maximize seed biomass (Seed_Biomass[K])
*
* The bilevel program is converted to a single-level MILP
* using strong duality of the inner LP and McCormick
* linearization of bilinear terms.
*************************************************************

$INLINECOM /* */

* =========================================================
* Command-line parameters with defaults
* =========================================================
$if not set DATADIR $set DATADIR "."
$if not set MAX_KO  $set MAX_KO "5"
$if not set MIN_BIO $set MIN_BIO "0.01"
$if not set BIGM    $set BIGM "1000"
$if not set TIMELIM $set TIMELIM "3600"
$if not set DAPNAME $set DAPNAME "unknown"

OPTIONS
        decimals = 8
        optcr = 0.05
        reslim = %TIMELIM%
        mip = cplex
        lp = cplex
        limrow = 0
        limcol = 0
        solprint = off
;

* =========================================================
* PART 1: DATA INPUT
* =========================================================

SETS
        i       metabolites
$include "%DATADIR%/metabolites.txt"

        j       reactions
$include "%DATADIR%/reactions.txt"
;

ALIAS(i, ii);
ALIAS(j, jj);

PARAMETERS
        S(i,j)          stoichiometric matrix
$include "%DATADIR%/sij.txt"

        v_max(j)        maximum flux bound
$include "%DATADIR%/v_max.txt"

        v_min(j)        minimum flux bound
$include "%DATADIR%/v_min.txt"
;

* =========================================================
* PART 2: DEFINE REACTION SETS
* =========================================================

* Key reactions
SET biomass_rxn(j)      biomass reaction
        / 'Seed_Biomass[K]' /;

SET target_rxn(j)       lysine target reaction (DAP decarboxylase)
        / 'R00451[K,p]' /;

* Protected reactions: NOT available for knockout
SET protected(j)        protected reactions
$include "protected_reactions.txt"
;

* Candidate reactions: available for knockout
SET candidates(j)       knockout candidate reactions;
candidates(j) = yes$(not protected(j));

SCALAR  num_candidates  number of candidate reactions;
num_candidates = card(candidates);
DISPLAY num_candidates;

* =========================================================
* PART 3: BASELINE FBA (Wild-Type Reference)
* =========================================================

FREE VARIABLES
        v_fba(j)        FBA reaction flux
        z_fba           FBA objective
;

EQUATIONS
        eq_fba_obj      FBA objective
        eq_fba_mb(i)    FBA mass balance
        eq_fba_lb(j)    FBA lower bound
        eq_fba_ub(j)    FBA upper bound
;

eq_fba_obj..            z_fba =e= v_fba('Seed_Biomass[K]');
eq_fba_mb(i)..          sum(j, S(i,j) * v_fba(j)) =e= 0;
eq_fba_lb(j)..          v_fba(j) =g= v_min(j);
eq_fba_ub(j)..          v_fba(j) =l= v_max(j);

MODEL fba_wt / eq_fba_obj, eq_fba_mb, eq_fba_lb, eq_fba_ub /;
SOLVE fba_wt USING lp MAXIMIZING z_fba;

SCALAR wt_biomass       wild-type maximum biomass;
SCALAR wt_lysine        wild-type lysine flux at max biomass;
SCALAR min_biomass      minimum allowed biomass in OptKnock;

wt_biomass = z_fba.l;
wt_lysine  = v_fba.l('R00451[K,p]');
min_biomass = %MIN_BIO% * wt_biomass;

PARAMETER wt_flux(j)    wild-type optimal flux distribution;
wt_flux(j) = v_fba.l(j);

DISPLAY wt_biomass, wt_lysine, min_biomass;

* =========================================================
* PART 4: OptKnock BILEVEL MILP FORMULATION
* =========================================================
*
* Single-level reformulation via strong duality:
*
* max  v_target  (R00451[K,p] - lysine production)
* s.t.:
*   [Primal feasibility]
*   sum_j S(i,j)*v(j) = 0                          for all i
*   v_min(j)*y(j) <= v(j) <= v_max(j)*y(j)         for all j
*
*   [Dual feasibility]
*   sum_i S(i,j)*lambda(i) - mu_lo(j) + mu_up(j) = c(j)  for all j
*   mu_lo(j) >= 0, mu_up(j) >= 0
*   (c(j) = 1 for biomass, 0 otherwise)
*
*   [Strong duality - linearized]
*   v_biomass = sum_j (v_max(j)*w_up(j) - v_min(j)*w_lo(j))
*   where w_up(j) = y(j)*mu_up(j), w_lo(j) = y(j)*mu_lo(j)
*
*   [Linearization: w = y*mu via McCormick/Glover]
*   0 <= w_up(j) <= M*y(j)
*   w_up(j) <= mu_up(j)
*   w_up(j) >= mu_up(j) - M*(1-y(j))
*   (analogous for w_lo)
*
*   [Knockout limit]
*   sum_{j in candidates} (1 - y(j)) <= K
*
*   [Minimum biomass]
*   v_biomass >= min_biomass
*
*   y(j) in {0,1} for candidates, y(j)=1 for protected
* =========================================================

* --- Variables ---
FREE VARIABLES
        vv(j)           primal reaction flux
        lambda(i)       dual of mass balance (free)
        z_ok            OptKnock objective (lysine flux)
;

POSITIVE VARIABLES
        mu_lo(j)        dual of lower bound constraint (>= 0)
        mu_up(j)        dual of upper bound constraint (>= 0)
        w_lo(j)         linearization variable: y(j) * mu_lo(j)
        w_up(j)         linearization variable: y(j) * mu_up(j)
;

BINARY VARIABLES
        y(j)            knockout indicator (1=active 0=knocked out)
;

* Fix y=1 for all protected (non-candidate) reactions
y.fx(j)$(protected(j)) = 1;

SCALAR M_val    big-M constant / %BIGM% /;
SCALAR K_val    current max knockouts / %MAX_KO% /;

* --- Equations ---
EQUATIONS
        eq_ok_obj               OptKnock outer objective
        eq_ok_mb(i)             primal mass balance
        eq_ok_lb(j)             primal lower bound with knockout
        eq_ok_ub(j)             primal upper bound with knockout
        eq_ok_dual(j)           dual feasibility
        eq_ok_sd                strong duality (linearized)
        eq_ok_ko_limit          knockout count limit
        eq_ok_min_bio           minimum biomass

* McCormick linearization for w_up = y * mu_up
        eq_lin_wup_ub1(j)       w_up <= M * y
        eq_lin_wup_ub2(j)       w_up <= mu_up
        eq_lin_wup_lb(j)        w_up >= mu_up - M*(1-y)

* McCormick linearization for w_lo = y * mu_lo
        eq_lin_wlo_ub1(j)       w_lo <= M * y
        eq_lin_wlo_ub2(j)       w_lo <= mu_lo
        eq_lin_wlo_lb(j)        w_lo >= mu_lo - M*(1-y)
;

* Outer objective: maximize lysine biosynthesis
eq_ok_obj..
        z_ok =e= vv('R00451[K,p]');

* Primal: mass balance
eq_ok_mb(i)..
        sum(j, S(i,j) * vv(j)) =e= 0;

* Primal: flux bounds (when y=0, both bounds become 0 => vv=0)
eq_ok_lb(j)..
        vv(j) =g= v_min(j) * y(j);

eq_ok_ub(j)..
        vv(j) =l= v_max(j) * y(j);

* Dual feasibility: S^T lambda - mu_lo + mu_up = c
*   where c(j)=1 for biomass, 0 otherwise
eq_ok_dual(j)..
        sum(i, S(i,j) * lambda(i)) - mu_lo(j) + mu_up(j)
        =e= 1$biomass_rxn(j);

* Strong duality (linearized bilinear terms):
*   v_biomass = sum_j (v_max(j)*w_up(j) - v_min(j)*w_lo(j))
eq_ok_sd..
        vv('Seed_Biomass[K]') =e=
        sum(j, v_max(j) * w_up(j) - v_min(j) * w_lo(j));

* Knockout limit: at most K_val knockouts among candidates
eq_ok_ko_limit..
        sum(j$candidates(j), 1 - y(j)) =l= K_val;

* Minimum biomass: ensure organism viability
eq_ok_min_bio..
        vv('Seed_Biomass[K]') =g= min_biomass;

* McCormick linearization: w_up(j) = y(j) * mu_up(j)
eq_lin_wup_ub1(j)..    w_up(j) =l= M_val * y(j);
eq_lin_wup_ub2(j)..    w_up(j) =l= mu_up(j);
eq_lin_wup_lb(j)..     w_up(j) =g= mu_up(j) - M_val * (1 - y(j));

* McCormick linearization: w_lo(j) = y(j) * mu_lo(j)
eq_lin_wlo_ub1(j)..    w_lo(j) =l= M_val * y(j);
eq_lin_wlo_ub2(j)..    w_lo(j) =l= mu_lo(j);
eq_lin_wlo_lb(j)..     w_lo(j) =g= mu_lo(j) - M_val * (1 - y(j));

* --- Define the OptKnock model ---
MODEL optknock / eq_ok_obj, eq_ok_mb, eq_ok_lb, eq_ok_ub,
                 eq_ok_dual, eq_ok_sd,
                 eq_ok_ko_limit, eq_ok_min_bio,
                 eq_lin_wup_ub1, eq_lin_wup_ub2, eq_lin_wup_lb,
                 eq_lin_wlo_ub1, eq_lin_wlo_ub2, eq_lin_wlo_lb /;

optknock.optfile = 1;
optknock.optcr = 0.05;
optknock.reslim = %TIMELIM%;

* =========================================================
* PART 5: SOLVE OptKnock FOR MULTIPLE K VALUES
* =========================================================

SET k_idx           knockout count index / k1*k5 /;
PARAMETER k_values(k_idx) / k1 1, k2 2, k3 3, k4 4, k5 5 /;

* Storage for results across K values
PARAMETER
        ok_lysine(k_idx)        OptKnock lysine flux
        ok_biomass(k_idx)       OptKnock biomass
        ok_modelstat(k_idx)     model status
        ok_solvestat(k_idx)     solver status
        ok_gap(k_idx)           optimality gap
        ok_knocked(k_idx,j)     1 if reaction j knocked out at K
        ok_flux(k_idx,j)        reaction fluxes at each K
;

* Verification FBA storage
PARAMETER
        ver_biomass(k_idx)      verified biomass after KO
        ver_lysine(k_idx)       verified lysine after KO
;

* Pessimistic/Optimistic verification variables and models
* (must be declared outside loop)
FREE VARIABLE z_pess   pessimistic objective;
FREE VARIABLE z_opti   optimistic objective;
SCALAR ver_bio_opt     verified max biomass;
PARAMETER pess_lysine(k_idx)  pessimistic (min) lysine at max biomass;
PARAMETER opti_lysine(k_idx)  optimistic (max) lysine at max biomass;

EQUATIONS
        eq_pess_obj            pessimistic objective
        eq_pess_fix_bio        fix biomass for pessimistic check
        eq_opti_obj            optimistic objective
        eq_opti_fix_bio        fix biomass for optimistic check
;
eq_pess_obj..     z_pess =e= v_fba('R00451[K,p]');
eq_pess_fix_bio.. v_fba('Seed_Biomass[K]') =g= ver_bio_opt - 1e-8;
eq_opti_obj..     z_opti =e= v_fba('R00451[K,p]');
eq_opti_fix_bio.. v_fba('Seed_Biomass[K]') =g= ver_bio_opt - 1e-8;

MODEL fba_pessimistic / eq_pess_obj, eq_fba_mb, eq_fba_lb, eq_fba_ub, eq_pess_fix_bio /;
MODEL fba_optimistic  / eq_opti_obj, eq_fba_mb, eq_fba_lb, eq_fba_ub, eq_opti_fix_bio /;

* Only solve up to the requested maximum K
SCALAR max_k / %MAX_KO% /;

FILE RESULTS / results_optknock.txt /;
PUT RESULTS;
PUT "================================================================" /;
PUT "  OptKnock Results: Lysine Overproduction in Maize Endosperm" /;
PUT "  DAP: %DAPNAME%" /;
PUT "  Data: %DATADIR%" /;
PUT "  Date: March 2026" /;
PUT "================================================================" /;
PUT "  WT Biomass (Seed_Biomass[K]): ", wt_biomass:0:8 /;
PUT "  WT Lysine  (R00451[K,p]):     ", wt_lysine:0:10 /;
PUT "  Min Biomass Fraction:          %MIN_BIO%" /;
PUT "  Min Biomass Value:             ", min_biomass:0:8 /;
PUT "  Big-M:                         %BIGM%" /;
PUT "  Candidate reactions:           ", num_candidates:0:0 /;
PUT "================================================================" /;
PUT /;

LOOP(k_idx$(k_values(k_idx) le max_k),

        K_val = k_values(k_idx);

        PUT "-------------------------------------------------------" /;
        PUT "  Solving OptKnock with K = ", K_val:0:0, " knockouts" /;
        PUT "-------------------------------------------------------" /;

* Reset binary variables for candidates
        y.l(j)$candidates(j) = 1;

* Solve OptKnock MILP
        SOLVE optknock USING mip MAXIMIZING z_ok;

* Store results
        ok_lysine(k_idx)    = z_ok.l;
        ok_biomass(k_idx)   = vv.l('Seed_Biomass[K]');
        ok_modelstat(k_idx) = optknock.modelstat;
        ok_solvestat(k_idx) = optknock.solvestat;
        ok_gap(k_idx)       = abs(optknock.objest - z_ok.l)
                              / max(abs(z_ok.l), 1e-10);

        loop(j$candidates(j),
                ok_knocked(k_idx,j) = 1$(y.l(j) < 0.5);
        );
        ok_flux(k_idx,j) = vv.l(j);

* Report
        PUT "  Model Status:  ", ok_modelstat(k_idx):0:0;
        if(ok_modelstat(k_idx) = 1, PUT " (OPTIMAL)";
        elseif ok_modelstat(k_idx) = 2, PUT " (locally optimal)";
        elseif ok_modelstat(k_idx) = 8, PUT " (integer solution)";
        else PUT " (other)";
        );
        PUT /;
        PUT "  Solver Status: ", ok_solvestat(k_idx):0:0 /;
        PUT "  Gap:           ", ok_gap(k_idx):0:6 /;
        PUT "  Lysine Flux:   ", ok_lysine(k_idx):0:10 /;
        PUT "  Biomass:       ", ok_biomass(k_idx):0:8 /;

        if(wt_lysine > 1e-12,
                PUT "  Lysine Fold-Change vs WT: ",
                    (ok_lysine(k_idx)/wt_lysine):0:4 /;
        else
                PUT "  Lysine Fold-Change vs WT: INF (WT lysine ~ 0)" /;
        );
        PUT "  Biomass % of WT: ",
            (ok_biomass(k_idx)/wt_biomass*100):0:2, "%" /;
        PUT /;

* List knocked-out reactions
        PUT "  Knocked-out reactions (K=", K_val:0:0, "):" /;
        loop(j$candidates(j),
                if(ok_knocked(k_idx,j) = 1,
                        PUT "    KO: ", j.tl:0:60,
                            "  (WT flux: ", wt_flux(j):0:8, ")" /;
                );
        );
        PUT /;

* -----------------------------------------------------------
* VERIFICATION: Standard FBA + Pessimistic + Optimistic checks
* -----------------------------------------------------------

* --- STEP A: Maximize biomass with knockouts ---
        PUT "  Verification:" /;

* Set bounds for verification
        v_fba.lo(j) = v_min(j);
        v_fba.up(j) = v_max(j);

* Apply knockouts
        loop(jj$candidates(jj),
                if(ok_knocked(k_idx,jj) = 1,
                        v_fba.fx(jj) = 0;
                );
        );

        SOLVE fba_wt USING lp MAXIMIZING z_fba;

        ver_biomass(k_idx) = z_fba.l;
        ver_lysine(k_idx)  = v_fba.l('R00451[K,p]');

        PUT "    A. Max biomass FBA:" /;
        PUT "       Biomass: ", ver_biomass(k_idx):0:8 /;
        PUT "       Lysine:  ", ver_lysine(k_idx):0:10 /;

* --- STEP B: Fix biomass at max, MINIMIZE lysine (pessimistic) ---
*     This tells us: the MINIMUM lysine guaranteed at max biomass
ver_bio_opt = z_fba.l;

* Keep knockout bounds active
        v_fba.lo(j) = v_min(j);
        v_fba.up(j) = v_max(j);
        loop(jj$candidates(jj),
                if(ok_knocked(k_idx,jj) = 1,
                        v_fba.fx(jj) = 0;
                );
        );

        SOLVE fba_pessimistic USING lp MINIMIZING z_pess;

pess_lysine(k_idx) = z_pess.l;
        PUT "    B. Pessimistic (min lysine at max biomass):" /;
        PUT "       Min Lysine: ", pess_lysine(k_idx):0:10 /;
        if(wt_lysine > 1e-12,
                PUT "       Min Lysine Fold vs WT: ",
                    (pess_lysine(k_idx)/wt_lysine):0:4 /;
        );

* --- STEP C: Fix biomass at max, MAXIMIZE lysine (optimistic) ---

        v_fba.lo(j) = v_min(j);
        v_fba.up(j) = v_max(j);
        loop(jj$candidates(jj),
                if(ok_knocked(k_idx,jj) = 1,
                        v_fba.fx(jj) = 0;
                );
        );

        SOLVE fba_optimistic USING lp MAXIMIZING z_opti;

opti_lysine(k_idx) = z_opti.l;
        PUT "    C. Optimistic (max lysine at max biomass):" /;
        PUT "       Max Lysine: ", opti_lysine(k_idx):0:10 /;
        if(wt_lysine > 1e-12,
                PUT "       Max Lysine Fold vs WT: ",
                    (opti_lysine(k_idx)/wt_lysine):0:4 /;
        );

        PUT "    Summary: At max biomass lysine range = [",
            pess_lysine(k_idx):0:10, " , ", opti_lysine(k_idx):0:10, "]" /;
        if(pess_lysine(k_idx) > wt_lysine * 1.01,
                PUT "    ** REAL GROWTH COUPLING: min lysine > WT lysine **" /;
        else
                PUT "    Note: No guaranteed growth coupling (min lysine ~ WT)" /;
        );
        PUT /;

* Reset FBA bounds for next iteration
        v_fba.lo(j) = v_min(j);
        v_fba.up(j) = v_max(j);

);

* =========================================================
* PART 6: SUMMARY TABLE
* =========================================================

PUT /;
PUT "================================================================" /;
PUT "  SUMMARY TABLE" /;
PUT "================================================================" /;
PUT "K", @5, "OptKnock_Lys", @20, "Biomass", @35, "PessLys", @50,
    "OptiLys", @65, "Bio%WT", @77, "Coupled?", @87, "Status" /;
PUT "---", @5, "------", @20, "-------", @35, "-------", @50,
    "-------", @65, "------", @77, "--------", @87, "------" /;

LOOP(k_idx$(k_values(k_idx) le max_k),
        PUT k_values(k_idx):0:0, @5,
            ok_lysine(k_idx):0:10, @20,
            ok_biomass(k_idx):0:8, @35,
            pess_lysine(k_idx):0:10, @50,
            opti_lysine(k_idx):0:10, @65,
            (ok_biomass(k_idx)/wt_biomass*100):0:2, @77;
        if(pess_lysine(k_idx) > wt_lysine * 1.01,
                PUT "YES", @87;
        else
                PUT "no", @87;
        );
        PUT ok_modelstat(k_idx):0:0 /;
);

PUT /;
PUT "================================================================" /;
PUT "  KNOCKOUT DETAILS BY K" /;
PUT "================================================================" /;

LOOP(k_idx$(k_values(k_idx) le max_k),
        PUT /;
        PUT "K = ", k_values(k_idx):0:0, ":" /;
        loop(j$(ok_knocked(k_idx,j) = 1),
                PUT "  ", j.tl:0:60, @65,
                    "  WT_flux=", wt_flux(j):0:8,
                    "  OK_flux=", ok_flux(k_idx,j):0:8 /;
        );
);

* =========================================================
* PART 7: FULL FLUX COMPARISON (for K = max_k)
* =========================================================

PUT /;
PUT "================================================================" /;
PUT "  FLUX COMPARISON: WT vs OptKnock (K=%MAX_KO%)" /;
PUT "================================================================" /;
PUT "Reaction", @55, "WT_flux", @70, "OK_flux", @85, "y" /;

SCALAR k_max_idx;
k_max_idx = max_k;

* Find the k_idx corresponding to max_k
LOOP(k_idx$(k_values(k_idx) = max_k),
        loop(j$(abs(wt_flux(j)) > 1e-8 or abs(ok_flux(k_idx,j)) > 1e-8),
                PUT j.tl:0:52, @55,
                    wt_flux(j):0:8, @70,
                    ok_flux(k_idx,j):0:8, @85;
                if(ok_knocked(k_idx,j) = 1,
                        PUT "KO";
                else
                        PUT "1";
                );
                PUT /;
        );
);

PUTCLOSE;

* =========================================================
* PART 8: CSV OUTPUT FOR PROGRAMMATIC ANALYSIS
* =========================================================

FILE CSV_SUMMARY / results_optknock_summary.csv /;
PUT CSV_SUMMARY;
PUT "DAP,K,OptKnock_Lysine,Biomass,Pessimistic_Lysine,Optimistic_Lysine,";
PUT "Biomass_PctWT,Growth_Coupled,ModelStatus,SolverStatus,Gap" /;

LOOP(k_idx$(k_values(k_idx) le max_k),
        PUT "%DAPNAME%,",
            k_values(k_idx):0:0, ",",
            ok_lysine(k_idx):0:10, ",",
            ok_biomass(k_idx):0:8, ",",
            pess_lysine(k_idx):0:10, ",",
            opti_lysine(k_idx):0:10, ",",
            (ok_biomass(k_idx)/wt_biomass*100):0:4, ",";
        if(pess_lysine(k_idx) > wt_lysine * 1.01,
                PUT "YES,";
        else
                PUT "no,";
        );
        PUT ok_modelstat(k_idx):0:0, ",",
            ok_solvestat(k_idx):0:0, ",",
            ok_gap(k_idx):0:6 /;
);
PUTCLOSE;

FILE CSV_KO / results_optknock_knockouts.csv /;
PUT CSV_KO;
PUT "DAP,K,Reaction,WT_Flux,OK_Flux" /;

LOOP(k_idx$(k_values(k_idx) le max_k),
        loop(j$(ok_knocked(k_idx,j) = 1),
                PUT "%DAPNAME%,",
                    k_values(k_idx):0:0, ',"',
                    j.tl:0:60, '",',
                    wt_flux(j):0:8, ",",
                    ok_flux(k_idx,j):0:8 /;
        );
);
PUTCLOSE;

* =========================================================
* DONE
* =========================================================
DISPLAY "OptKnock analysis complete.";
DISPLAY ok_lysine, ok_biomass, ok_modelstat;
