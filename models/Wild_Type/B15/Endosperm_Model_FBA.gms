*************************************************************
*Basic FBA model for Arabidopsis root
*************************************************************
*****************Lohani Beer************************
*************************************************************
$INLINECOM /*  */

OPTIONS

        decimals = 8
	lp = cplex
;

*********Defining Sets**************************************
SETS

	i					set of metabolites

$include "metabolites.txt"	

	j					set of reactions

$include "reactions.txt"


;
*************************************************************

***********Defining Parameters*******************************
PARAMETERS

	S(i,j)					stoichiometric matrix

$include "sij.txt"

	v_max(j)				maximum flux of v(j)
	
$include "v_max.txt"

	v_min(j)				minimum flux of v(j)

$include "v_min.txt"




;
**************************************************************

*********Defining Equations***********************************
EQUATIONS

	objective				objective function
	mass_balance1(i)			steady state mass balance
	
	lower_bound(j)				lower bounds on reactions
	upper_bound(j)				upper bounds on reactions
;
**************************************************************

*********Defining Variables***********************************
FREE VARIABLES

	v(j)					reaction flux
	obj					objective value
;

****************************************************************

***************Defining Model***********************************
objective..			obj =e= v('Seed_Biomass[K]');

mass_balance1(i)..		sum(j, S(i,j) * v(j)) =e= 0;


lower_bound(j)..		v_min(j) =l= v(j);

upper_bound(j)..		v(j) =l= v_max(j);

***************Set bounds***********************************


Model arabidopsis_root /all/;

******************************************************************

**********Solving Model*********************

Solve arabidopsis_root using lp maximizing obj;

arabidopsis_root.optfile = 1;

arabidopsis_root.holdfixed = 1;


********************************************
****************Output File*****************

FILE RESULTS /results_FBA.txt/;

PUT RESULTS;

PUT "The max Biomass value is : " obj.l:0:8//;

loop(j,
	put j.tl:0:30, "    ", v.l(j):0:8/;
);

PUTCLOSE;
**********************************************