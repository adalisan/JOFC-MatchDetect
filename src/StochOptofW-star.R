# TODO: Add comment
# Author: Sancar
# Stochastic optimization of w^*
###############################################################################

eval.beta.at.w.est(w.est){
  
run.result <-run.mc.replicate("gaussian",p, r, q, c.val,  ##try(
  			d           = d,
				  pprime1=pprime1,   pprime2=pprime2,
				Wchoice     = Wchoice, 
				pre.scaling = pre.scaling,
				oos         = oos,
				alpha       = alpha,
				n = n, m = m,
			
				old.gauss.model.param=old.gauss.model.param,
				separability.entries.w=separability.entries.w,
				compare.pom.cca=FALSE,
				oos.use.imputed=oos.use.imputed,
				level.mcnemar=level.mcnemar,
				def.w=def.w,
				rival.w=rival.w,
				proc.dilation=proc.dilation,
				assume.matched.for.oos=assume.matched.for.oos,
				w.vals=w.est,
				wt.equalize=wt.equalize,
				verbose=verbose,
				power.comparison.test=power.comparison.test)


}