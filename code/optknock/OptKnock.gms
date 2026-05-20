*************************************************************
* OptKnock: Bilevel Optimization for Target Metabolite
*           Overproduction via Gene Knockouts
*
* Reference: Burgard AP, Pharkya P, Maranas CD (2003)
*   "OptKnock: A bilevel programming framework for
*    identifying gene knockout strategies for microbial
*    strain optimization"
*   Biotechnology and Bioengineering 84(6):647-657
*
* The bilevel program is converted to a single-level MILP
* using strong duality of the inner LP and McCormick
* linearization of bilinear terms (y * mu).
*
* ==========================================================
* REQUIRED INPUT FILES (in DATADIR):
*   reactions.txt    - GAMS set of reaction identifiers
*   metabolites.txt  - GAMS set of metabolite identifiers
*   sij.txt          - GAMS parameter: stoichiometric matrix
*   upper_bound.txt  - GAMS parameter: upper flux bounds
*   lower_bound.txt  - GAMS parameter: lower flux bounds
*   (bound filenames are configurable via UBFILE / LBFILE)
*
* OPTIONAL INPUT FILE (in working directory):
*   protected_reactions.txt - GAMS set of reactions NOT
*       allowed to be knocked out. If absent, only the
*       biomass and target reactions are protected.
*
* CONFIGURATION (optknock_config.inc in working directory):
*   Because GAMS CLI cannot handle brackets or commas in
*   parameter values, set BIOMASS and TARGET in this file:
*     $setglobal BIOMASS  YourBiomassReactionName
*     $setglobal TARGET   YourTargetReactionName
*     $setglobal UBFILE   upper_bound.txt
*     $setglobal LBFILE   lower_bound.txt
*   The run script generates this file automatically.
*
* CONFIGURATION:
*   Create optknock_config.inc in the working directory with:
*     $setglobal BIOMASS YourBiomassReactionName
*     $setglobal TARGET  YourTargetReactionName
*   (The run script generates this file automatically.)
*
* COMMAND-LINE PARAMETERS:
*   --DATADIR   Path to folder with input files   [default: .]
*   --MAX_KO    Maximum number of simultaneous KOs [default: 5]
*   --MIN_BIO   Min biomass as fraction of WT      [default: 0.01]
*   --BIGM      Big-M constant for McCormick       [default: 1000]
*   --TIMELIM   Solver time limit in seconds       [default: 3600]
*   --DAPNAME   Label for this run (appears in CSV)[default: run]
*
* EXAMPLE USAGE:
*   echo '$setglobal BIOMASS Biomass'  >  optknock_config.inc
*   echo '$setglobal TARGET EX_lys_e' >> optknock_config.inc
*   gams OptKnock.gms --DATADIR=./data --MAX_KO=3
*        --MIN_BIO=0.10 --DAPNAME=WT lo=3
*************************************************************

$INLINECOM /* */

* =========================================================
* 0. COMMAND-LINE PARAMETERS WITH DEFAULTS
* =========================================================
$if not set DATADIR $set DATADIR "."
$if not set MAX_KO  $set MAX_KO  "5"
$if not set MIN_BIO $set MIN_BIO "0.01"
$if not set BIGM    $set BIGM    "1000"
$if not set TIMELIM $set TIMELIM "3600"
$if not set DAPNAME $set DAPNAME "run"

* BIOMASS and TARGET reaction names.
*
* Because GAMS command-line parsing does not handle brackets,
* commas, or other special characters in parameter values,
* these are best set in an include file (optknock_config.inc).
*
* Option A - include file (recommended for complex names):
*   Create optknock_config.inc containing:
*     $setglobal BIOMASS Seed_Biomass[K]
*     $setglobal TARGET  R00451[K,p]
*   The run script generates this automatically.
*
* Option B - direct set (works when names have no commas):
*   Pass --BIOMASS=Biomass --TARGET=EX_lys on the command line.
*
$if exist "optknock_config.inc" $include "optknock_config.inc"

* Fallback defaults if neither config file nor CLI was used
$if not set BIOMASS $set BIOMASS "__BIOMASS_NOT_SET__"
$if not set TARGET  $set TARGET  "__TARGET_NOT_SET__"
$if not set UBFILE  $set UBFILE  "upper_bound.txt"
$if not set LBFILE  $set LBFILE  "lower_bound.txt"

OPTIONS
        decimals = 8
        optcr    = 0.05
        reslim   = %TIMELIM%
        mip      = cplex
        lp       = cplex
        limrow   = 0
        limcol   = 0
        solprint = off
;

* =========================================================
* 1. DATA INPUT
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

        ub(j)           upper flux bound
$include "%DATADIR%/%UBFILE%"

        lb(j)           lower flux bound
$include "%DATADIR%/%LBFILE%"
;

* =========================================================
* 2. REACTION CLASSIFICATION
* =========================================================

* Key reactions (user-defined via command line)
SET biomass_rxn(j)      biomass reaction
        / '%BIOMASS%' /;

SET target_rxn(j)       target reaction to maximise
        / '%TARGET%' /;

* --- Validate that biomass and target exist in model ---
ABORT$( card(biomass_rxn) = 0 )
  "BIOMASS reaction not found in reactions.txt. "
  "Check the --BIOMASS parameter.";

ABORT$( card(target_rxn) = 0 )
  "TARGET reaction not found in reactions.txt. "
  "Check the --TARGET parameter.";

* --- Protected reactions (not available for knockout) ---
*   If protected_reactions.txt exists in the working directory,
*   load it; otherwise, we protect only biomass + target.
SET protected(j)        protected reactions
$if     exist "protected_reactions.txt"  $include "protected_reactions.txt"
$if not exist "protected_reactions.txt"  / /
;

* If no protected file was provided, seed with biomass + target
$if not exist "protected_reactions.txt"  protected(j)$(biomass_rxn(j)) = yes;
$if not exist "protected_reactions.txt"  protected(j)$(target_rxn(j))  = yes;

* Always ensure biomass and target are protected, even if
* the user's protected_reactions.txt omitted them
protected(j)$(biomass_rxn(j)) = yes;
protected(j)$(target_rxn(j))  = yes;

* Candidate reactions: everything that is NOT protected
SET candidates(j)       knockout candidate reactions;
candidates(j) = yes$(not protected(j));

SCALAR num_candidates  number of candidate reactions;
num_candidates = card(candidates);
DISPLAY num_candidates;

* =========================================================
* 3. BASELINE FBA (WILD-TYPE REFERENCE)
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

eq_fba_obj..            z_fba =e= v_fba('%BIOMASS%');
eq_fba_mb(i)..          sum(j, S(i,j) * v_fba(j)) =e= 0;
eq_fba_lb(j)..          v_fba(j) =g= lb(j);
eq_fba_ub(j)..          v_fba(j) =l= ub(j);

MODEL fba_wt / eq_fba_obj, eq_fba_mb, eq_fba_lb, eq_fba_ub /;
SOLVE fba_wt USING lp MAXIMIZING z_fba;

SCALAR wt_biomass       wild-type maximum biomass;
SCALAR wt_target        wild-type target flux at max biomass;
SCALAR min_biomass      minimum allowed biomass in OptKnock;

wt_biomass  = z_fba.l;
wt_target   = v_fba.l('%TARGET%');
min_biomass = %MIN_BIO% * wt_biomass;

PARAMETER wt_flux(j)    wild-type optimal flux distribution;
wt_flux(j) = v_fba.l(j);

DISPLAY wt_biomass, wt_target, min_biomass;

* =========================================================
* 4. OptKnock BILEVEL MILP FORMULATION
* =========================================================
*
* Single-level reformulation via strong duality:
*
* max  v(TARGET)
* s.t.:
*   [Primal feasibility]
*   sum_j S(i,j)*v(j) = 0                     for all i
*   lb(j)*y(j) <= v(j) <= ub(j)*y(j)          for all j
*
*   [Dual feasibility]
*   sum_i S(i,j)*lambda(i) - mu_lo(j) + mu_up(j) = c(j)
*   mu_lo(j) >= 0,  mu_up(j) >= 0
*   (c(j) = 1 for biomass, 0 otherwise)
*
*   [Strong duality -- linearised]
*   v(BIOMASS) = sum_j ( ub(j)*w_up(j) - lb(j)*w_lo(j) )
*   where w_up = y * mu_up,  w_lo = y * mu_lo
*
*   [McCormick / Glover linearisation of w = y * mu]
*   0 <= w_up(j) <= M * y(j)
*   w_up(j) <= mu_up(j)
*   w_up(j) >= mu_up(j) - M*(1-y(j))
*   (analogous for w_lo)
*
*   [Knockout budget]
*   sum_{j in candidates} (1 - y(j)) <= K
*
*   [Minimum biomass]
*   v(BIOMASS) >= min_biomass
*
*   y(j) in {0,1} for candidates,  y(j) = 1 for protected
*
* =========================================================

* --- Variables ---
FREE VARIABLES
        vv(j)           primal flux in bilevel model
        lambda(i)       dual of mass balance (free)
        z_ok            OptKnock objective (target flux)
;

POSITIVE VARIABLES
        mu_lo(j)        dual of lower bound ( >= 0 )
        mu_up(j)        dual of upper bound ( >= 0 )
        w_lo(j)         linearisation: y(j) * mu_lo(j)
        w_up(j)         linearisation: y(j) * mu_up(j)
;

BINARY VARIABLES
        y(j)            knockout indicator -- 1 active or 0 knocked out
;

* Fix y = 1 for all protected (non-candidate) reactions
y.fx(j)$(protected(j)) = 1;

SCALAR M_val    big-M constant / %BIGM% /;
SCALAR K_val    current max knockouts / %MAX_KO% /;

* --- Equations ---
EQUATIONS
        eq_ok_obj               outer objective: max target flux
        eq_ok_mb(i)             primal mass balance
        eq_ok_lb(j)             primal lower bound with knockout
        eq_ok_ub(j)             primal upper bound with knockout
        eq_ok_dual(j)           dual feasibility
        eq_ok_sd                strong duality (linearised)
        eq_ok_ko_limit          knockout budget
        eq_ok_min_bio           minimum biomass

* McCormick for w_up = y * mu_up
        eq_lin_wup_ub1(j)       w_up <= M * y
        eq_lin_wup_ub2(j)       w_up <= mu_up
        eq_lin_wup_lb(j)        w_up >= mu_up - M*(1-y)

* McCormick for w_lo = y * mu_lo
        eq_lin_wlo_ub1(j)       w_lo <= M * y
        eq_lin_wlo_ub2(j)       w_lo <= mu_lo
        eq_lin_wlo_lb(j)        w_lo >= mu_lo - M*(1-y)
;

* Outer objective: maximise target metabolite production
eq_ok_obj..
        z_ok =e= vv('%TARGET%');

* Primal: mass balance
eq_ok_mb(i)..
        sum(j, S(i,j) * vv(j)) =e= 0;

* Primal: flux bounds (y = 0 forces vv = 0 => reaction knocked out)
eq_ok_lb(j)..
        vv(j) =g= lb(j) * y(j);

eq_ok_ub(j)..
        vv(j) =l= ub(j) * y(j);

* Dual feasibility: S^T lambda - mu_lo + mu_up = c
*   c(j) = 1 for biomass, 0 otherwise
eq_ok_dual(j)..
        sum(i, S(i,j) * lambda(i)) - mu_lo(j) + mu_up(j)
        =e= 1$biomass_rxn(j);

* Strong duality (linearised bilinear terms)
eq_ok_sd..
        vv('%BIOMASS%') =e=
        sum(j, ub(j) * w_up(j) - lb(j) * w_lo(j));

* Knockout budget: at most K_val knockouts among candidates
eq_ok_ko_limit..
        sum(j$candidates(j), 1 - y(j)) =l= K_val;

* Minimum biomass: ensure viability
eq_ok_min_bio..
        vv('%BIOMASS%') =g= min_biomass;

* McCormick linearisation: w_up(j) = y(j) * mu_up(j)
eq_lin_wup_ub1(j)..    w_up(j) =l= M_val * y(j);
eq_lin_wup_ub2(j)..    w_up(j) =l= mu_up(j);
eq_lin_wup_lb(j)..     w_up(j) =g= mu_up(j) - M_val * (1 - y(j));

* McCormick linearisation: w_lo(j) = y(j) * mu_lo(j)
eq_lin_wlo_ub1(j)..    w_lo(j) =l= M_val * y(j);
eq_lin_wlo_ub2(j)..    w_lo(j) =l= mu_lo(j);
eq_lin_wlo_lb(j)..     w_lo(j) =g= mu_lo(j) - M_val * (1 - y(j));

MODEL optknock / eq_ok_obj, eq_ok_mb, eq_ok_lb, eq_ok_ub,
                 eq_ok_dual, eq_ok_sd,
                 eq_ok_ko_limit, eq_ok_min_bio,
                 eq_lin_wup_ub1, eq_lin_wup_ub2, eq_lin_wup_lb,
                 eq_lin_wlo_ub1, eq_lin_wlo_ub2, eq_lin_wlo_lb /;

optknock.optfile = 1;
optknock.optcr   = 0.05;
optknock.reslim  = %TIMELIM%;

* =========================================================
* 5. SOLVE OptKnock FOR K = 1 .. MAX_KO
* =========================================================

SET k_idx           knockout count index / k1*k5 /;
PARAMETER k_values(k_idx) / k1 1, k2 2, k3 3, k4 4, k5 5 /;

* Results storage
PARAMETER
        ok_target(k_idx)        OptKnock target flux
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
        ver_target(k_idx)       verified target flux after KO
;

* Pessimistic / optimistic verification variables and models
* (declared outside loop per GAMS rules)
FREE VARIABLE z_pess   pessimistic objective;
FREE VARIABLE z_opti   optimistic objective;
SCALAR ver_bio_opt     verified max biomass;
PARAMETER pess_target(k_idx)  pessimistic (min) target at max biomass;
PARAMETER opti_target(k_idx)  optimistic (max) target at max biomass;

EQUATIONS
        eq_pess_obj            pessimistic: min target at max biomass
        eq_pess_fix_bio        fix biomass for pessimistic check
        eq_opti_obj            optimistic:  max target at max biomass
        eq_opti_fix_bio        fix biomass for optimistic check
;
eq_pess_obj..     z_pess =e= v_fba('%TARGET%');
eq_pess_fix_bio.. v_fba('%BIOMASS%') =g= ver_bio_opt - 1e-8;
eq_opti_obj..     z_opti =e= v_fba('%TARGET%');
eq_opti_fix_bio.. v_fba('%BIOMASS%') =g= ver_bio_opt - 1e-8;

MODEL fba_pessimistic / eq_pess_obj, eq_fba_mb, eq_fba_lb, eq_fba_ub,
                         eq_pess_fix_bio /;
MODEL fba_optimistic  / eq_opti_obj, eq_fba_mb, eq_fba_lb, eq_fba_ub,
                         eq_opti_fix_bio /;

SCALAR max_k / %MAX_KO% /;

* =========================================================
* OUTPUT FILES
* =========================================================
FILE RESULTS / results_optknock.txt /;
PUT RESULTS;
PUT "================================================================" /;
PUT "  OptKnock Results" /;
PUT "  Run label:  %DAPNAME%" /;
PUT "  Data dir:   %DATADIR%" /;
PUT "  Biomass:    %BIOMASS%" /;
PUT "  Target:     %TARGET%" /;
PUT "================================================================" /;
PUT "  WT Biomass: ", wt_biomass:0:8 /;
PUT "  WT Target:  ", wt_target:0:10 /;
PUT "  Min Biomass Fraction: %MIN_BIO%" /;
PUT "  Min Biomass Value:    ", min_biomass:0:8 /;
PUT "  Big-M:                %BIGM%" /;
PUT "  Candidates:           ", num_candidates:0:0 /;
PUT "================================================================" /;
PUT /;

* =========================================================
* MAIN LOOP
* =========================================================
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
        ok_target(k_idx)    = z_ok.l;
        ok_biomass(k_idx)   = vv.l('%BIOMASS%');
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
        PUT "  Target Flux:   ", ok_target(k_idx):0:10 /;
        PUT "  Biomass:       ", ok_biomass(k_idx):0:8 /;

        if(wt_target > 1e-12,
                PUT "  Target Fold-Change vs WT: ",
                    (ok_target(k_idx)/wt_target):0:4 /;
        else
                PUT "  Target Fold-Change vs WT: INF (WT target ~ 0)" /;
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
* VERIFICATION: FBA + Pessimistic + Optimistic checks
* -----------------------------------------------------------

* --- STEP A: Maximise biomass with knockouts applied ---
        PUT "  Verification:" /;

        v_fba.lo(j) = lb(j);
        v_fba.up(j) = ub(j);

        loop(jj$candidates(jj),
                if(ok_knocked(k_idx,jj) = 1,
                        v_fba.fx(jj) = 0;
                );
        );

        SOLVE fba_wt USING lp MAXIMIZING z_fba;

        ver_biomass(k_idx) = z_fba.l;
        ver_target(k_idx)  = v_fba.l('%TARGET%');

        PUT "    A. Max biomass FBA:" /;
        PUT "       Biomass: ", ver_biomass(k_idx):0:8 /;
        PUT "       Target:  ", ver_target(k_idx):0:10 /;

* --- STEP B: Fix biomass at max, MINIMISE target (pessimistic) ---
ver_bio_opt = z_fba.l;

        v_fba.lo(j) = lb(j);
        v_fba.up(j) = ub(j);
        loop(jj$candidates(jj),
                if(ok_knocked(k_idx,jj) = 1,
                        v_fba.fx(jj) = 0;
                );
        );

        SOLVE fba_pessimistic USING lp MINIMIZING z_pess;

pess_target(k_idx) = z_pess.l;
        PUT "    B. Pessimistic (min target at max biomass):" /;
        PUT "       Min Target: ", pess_target(k_idx):0:10 /;
        if(wt_target > 1e-12,
                PUT "       Min Target Fold vs WT: ",
                    (pess_target(k_idx)/wt_target):0:4 /;
        );

* --- STEP C: Fix biomass at max, MAXIMISE target (optimistic) ---

        v_fba.lo(j) = lb(j);
        v_fba.up(j) = ub(j);
        loop(jj$candidates(jj),
                if(ok_knocked(k_idx,jj) = 1,
                        v_fba.fx(jj) = 0;
                );
        );

        SOLVE fba_optimistic USING lp MAXIMIZING z_opti;

opti_target(k_idx) = z_opti.l;
        PUT "    C. Optimistic (max target at max biomass):" /;
        PUT "       Max Target: ", opti_target(k_idx):0:10 /;
        if(wt_target > 1e-12,
                PUT "       Max Target Fold vs WT: ",
                    (opti_target(k_idx)/wt_target):0:4 /;
        );

        PUT "    Summary: At max biomass target range = [",
            pess_target(k_idx):0:10, " , ", opti_target(k_idx):0:10, "]" /;
        if(pess_target(k_idx) > wt_target * 1.01,
                PUT "    ** REAL GROWTH COUPLING: min target > WT target **" /;
        else
                PUT "    Note: No guaranteed growth coupling (min target ~ WT)" /;
        );
        PUT /;

* Reset FBA bounds for next iteration
        v_fba.lo(j) = lb(j);
        v_fba.up(j) = ub(j);

);

* =========================================================
* 6. SUMMARY TABLE
* =========================================================

PUT /;
PUT "================================================================" /;
PUT "  SUMMARY TABLE" /;
PUT "================================================================" /;
PUT "K", @5, "OptKnock_Tgt", @20, "Biomass", @35, "PessTgt", @50,
    "OptiTgt", @65, "Bio%WT", @77, "Coupled?", @87, "Status" /;
PUT "---", @5, "------", @20, "-------", @35, "-------", @50,
    "-------", @65, "------", @77, "--------", @87, "------" /;

LOOP(k_idx$(k_values(k_idx) le max_k),
        PUT k_values(k_idx):0:0, @5,
            ok_target(k_idx):0:10, @20,
            ok_biomass(k_idx):0:8, @35,
            pess_target(k_idx):0:10, @50,
            opti_target(k_idx):0:10, @65,
            (ok_biomass(k_idx)/wt_biomass*100):0:2, @77;
        if(pess_target(k_idx) > wt_target * 1.01,
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
* 7. FULL FLUX COMPARISON (K = max_k)
* =========================================================

PUT /;
PUT "================================================================" /;
PUT "  FLUX COMPARISON: WT vs OptKnock (K=%MAX_KO%)" /;
PUT "================================================================" /;
PUT "Reaction", @55, "WT_flux", @70, "OK_flux", @85, "y" /;

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
* 8. CSV OUTPUT FOR PROGRAMMATIC ANALYSIS
* =========================================================

FILE CSV_SUMMARY / results_optknock_summary.csv /;
PUT CSV_SUMMARY;
PUT "DAP,K,OptKnock_Target,Biomass,Pessimistic_Target,Optimistic_Target,";
PUT "Biomass_PctWT,Growth_Coupled,ModelStatus,SolverStatus,Gap" /;

LOOP(k_idx$(k_values(k_idx) le max_k),
        PUT "%DAPNAME%,",
            k_values(k_idx):0:0, ",",
            ok_target(k_idx):0:10, ",",
            ok_biomass(k_idx):0:8, ",",
            pess_target(k_idx):0:10, ",",
            opti_target(k_idx):0:10, ",",
            (ok_biomass(k_idx)/wt_biomass*100):0:4, ",";
        if(pess_target(k_idx) > wt_target * 1.01,
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
DISPLAY ok_target, ok_biomass, ok_modelstat;
