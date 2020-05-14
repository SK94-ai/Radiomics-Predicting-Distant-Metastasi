# Based on Harrell's val.prob function
# - scaled Brier score by relating to max for average calibrated Null model
# - risk distribution according to outcome
# - 0 and 1 to indicate outcome label; set with d1lab="..", d0lab=".."
# - labels: y axis: "Observed Frequency"; Triangle: "Grouped patients"
# - confidence intervals around triangles
# - a cut-off can be plotted; set x coordinate

# work done by Yvonne Vergouwe & Ewout Steyerberg

val.prob.ci<-
function(p, y, logit, group, weights = rep(1, length(y)), normwt = F, pl = T, 
	smooth = T, logistic.cal = F, xlab = "Predicted probability", ylab = 
	"Observed frequency", xlim = c(-0.02, 1),ylim = c(-0.15,1), m, g, cuts, emax.lim = c(0, 1), 
	legendloc =  c(0.55 , 0.27), statloc = c(0,.85),dostats=c(3,11,15,2,12,13),
	riskdist = "predicted", cex = 0.75, mkh = 0.02, connect.group = 
	F, connect.smooth = T, g.group = 4, evaluate = 100, nmin = 0, d0lab="0", d1lab="1", cex.d01=0.7,
  dist.label=0.04, cutoff, cex.lab=1, las=1, length.seg=1)
{
	if(missing(p))
		p <- 1/(1 + exp( - logit))
	else logit <- log(p/(1 - p))
	if(length(p) != length(y))
		stop("lengths of p or logit and y do not agree")
	names(p) <- names(y) <- names(logit) <- NULL
	if(!missing(group)) {
		if(length(group) == 1 && is.logical(group) && group)
			group <- rep("", length(y))
		if(!is.factor(group))
			group <- if(is.logical(group) || is.character(group)) 
				  as.factor(group) else cut2(group, g = 
				  g.group)
		names(group) <- NULL
		nma <- !(is.na(p + y + weights) | is.na(group))
		ng <- length(levels(group))
	}
	else {
		nma <- !is.na(p + y + weights)
		ng <- 0
	}
	logit <- logit[nma]
	y <- y[nma]
	p <- p[nma]
	if(ng > 0) {
		group <- group[nma]
		weights <- weights[nma]
		return(val.probg(p, y, group, evaluate, weights, normwt, nmin)
			)
	}

	if(length(unique(p)) == 1) {
#22Sep94
		P <- mean(y)
		Intc <- log(P/(1 - P))
		n <- length(y)
		D <- -1/n
		L01 <- -2 * sum(y * logit - log(1 + exp(logit)), na.rm = T)
		L.cal <- -2 * sum(y * Intc - log(1 + exp(Intc)), na.rm = T)
		U.chisq <- L01 - L.cal
		U.p <- 1 - pchisq(U.chisq, 1)
		U <- (U.chisq - 1)/n
		Q <- D - U

		stats <- c(0, 0.5, 0, D, 0, 1, U, U.chisq, U.p, Q, mean((y - p[
			1])^2), Intc, 0, rep(abs(p[1] - P), 2))
		names(stats) <- c("Dxy", "C (ROC)", "R2", "D", "D:Chi-sq", 
			"D:p", "U", "U:Chi-sq", "U:p", "Q", "Brier", 
			"Intercept", "Slope", "Emax", "Eavg")
		return(stats)
	}
	
	i <- !is.infinite(logit)
	nm <- sum(!i)
	if(nm > 0)
		warning(paste(nm, 
			"observations deleted from logistic calibration due to probs. of 0 or 1"
			))
	f <- lrm.fit(logit[i], y[i])
	f2<-	lrm.fit(offset=logit[i], y=y[i])
	stats <- f$stats
	n <- stats["Obs"]
	predprob <- seq(emax.lim[1], emax.lim[2], by = 0.0005)
	lt <- f$coef[1] + f$coef[2] * log(predprob/(1 - predprob))
	calp <- 1/(1 + exp( - lt))
	emax <- max(abs(predprob - calp))
    if (pl) {
        plot(0.5, 0.5, xlim = xlim, ylim = ylim, type = "n", xlab = xlab, 
            ylab = ylab)
        abline(0, 1, lty = 2)
        lt <- 2
        leg <- "Ideal"
        marks <- -1
        if (logistic.cal) {
            lt <- c(lt, 1)
            leg <- c(leg, "Logistic calibration")
            marks <- c(marks, -1)
        }
        if (smooth) {
            Sm <- lowess(p, y, iter = 0)
            if (connect.smooth) {
                lines(Sm, lty = 3)
                lt <- c(lt, 3)
                marks <- c(marks, -1)
            }
            else {
                points(Sm)
                lt <- c(lt, 0)
                marks <- c(marks, 1)
            }
            leg <- c(leg, "Nonparametric")
            cal.smooth <- approx(Sm, xout = p)$y
            eavg <- mean(abs(p - cal.smooth))
        }
      if(!missing(m) | !missing(g) | !missing(cuts)) {
			if(!missing(m))
				q <- cut2(p, m = m, levels.mean = T, digits = 7)
			else if(!missing(g))
				q <- cut2(p, g = g, levels.mean = T, digits = 7)
			else if(!missing(cuts))
			q <- cut2(p, cuts = cuts, levels.mean = T, digits = 7)
			means <- as.single(levels(q))
			prop <- tapply(y, q, function(x)mean(x, na.rm = T))
			points(means, prop, pch = 2, cex=cex)
#18.11.02: CI triangles			
			ng	<-tapply(y, q, length)
			og	<-tapply(y, q, sum)
			ob	<-og/ng
			se.ob	<-sqrt(ob*(1-ob)/ng)
			g		<- length(as.single(levels(q)))
		
			#for (i in 1:g) lines(c(means[i], means[i]), c(prop[i],min(1,prop[i]+1.96*se.ob[i])), type="l")
			#for (i in 1:g) lines(c(means[i], means[i]), c(prop[i],max(0,prop[i]-1.96*se.ob[i])), type="l")

			if(connect.group) {
				lines(means, prop)
				lt <- c(lt, 1)
			}
			else lt <- c(lt, 0)
			leg <- c(leg, "Grouped patients")
			marks <- c(marks, 2)
		}
	}
	lr <- stats["Model L.R."]
  p.lr <- stats["P"]
  D <- (lr - 1)/n
  L01 <- -2 * sum(y * logit - logb(1 + exp(logit)), na.rm = TRUE)
  U.chisq <- L01 - f$deviance[2]
  p.U <- 1 - pchisq(U.chisq, 2)
  U <- (U.chisq - 2)/n
  Q <- D - U
  Dxy <- stats["Dxy"]
  C <- stats["C"]
  R2 <- stats["R2"]
  B <- sum((p - y)^2)/n
# ES 15dec08 add Brier scaled
  Bmax  <- mean(y) * (1-mean(y))^2 + (1-mean(y)) * mean(y)^2
  Bscaled <- 1 - B/Bmax
  stats <- c(Dxy, C, R2, D, lr, p.lr, U, U.chisq, p.U, Q, B, 
        f$coef, emax, Bscaled)
    names(stats) <- c("Dxy", "C (ROC)", "R2", "D", "D:Chi-sq", 
        "D:p", "U", "U:Chi-sq", "U:p", "Q", "Brier", "Intercept", 
        "Slope", "Emax", "Brier scaled")
 	if(smooth)
		stats <- c(stats, c(Eavg = eavg))
	
	# Cut off definition	
	if(!missing(cutoff)) {
    	arrows(x0=cutoff,y0=.1,x1=cutoff,y1=-0.025,length=.15)
		}	
	if(pl) {
		logit <- seq(-7, 7, length = 200)
		prob <- 1/(1 + exp( - logit))
		pred.prob <- f$coef[1] + f$coef[2] * logit
		pred.prob <- 1/(1 + exp( - pred.prob))
		if(logistic.cal) lines(prob, pred.prob, lty = 1)	
	#	pc <- rep(" ", length(lt))
#	pc[lt==0] <- "."
		lp <- legendloc
  if (!is.logical(lp)) {
         if (!is.list(lp)) 
                lp <- list(x = lp[1], y = lp[2])
                legend(lp, leg, lty = lt, pch = marks, cex = cex, bty = "n")
    }
		if(!is.logical(statloc)) {
			dostats <- dostats
			leg <- format(names(stats)[dostats])	#constant length
			leg <- paste(leg, ":", format(stats[dostats]), sep = 
				"")
			if(!is.list(statloc))
				statloc <- list(x = statloc[1], y = statloc[2])
			text(statloc, paste(format(names(stats[dostats])), 
				collapse = "\n"), adj = 0, cex = cex)
			text(statloc$x + 0.4 , statloc$y, paste(
				format(round(stats[dostats], 3)), collapse = 
				"\n"), adj = 1, cex = cex)	
		}
		if(is.character(riskdist)) {
			if(riskdist == "calibrated") {
				x <- f$coef[1] + f$coef[2] * log(p/(1 - p))
				x <- 1/(1 + exp( - x))
				x[p == 0] <- 0
				x[p == 1] <- 1
			}
			else x <- p
			bins <- seq(0, 1, length = 101)
			x <- x[x >= 0 & x <= 1]
#08.04.01,yvon: distribution of predicted prob according to outcome
			f0	<-table(cut(x[y==0],bins))
			f1	<-table(cut(x[y==1],bins))
			j0	<-f0 > 0
			j1	<-f1 > 0
			bins0 <-(bins[-101])[j0]
			bins1 <-(bins[-101])[j1]
			f0	<-f0[j0]
			f1	<-f1[j1]
			maxf <-max(f0,f1)
			f0	<-(0.1*f0)/maxf
			f1	<-(0.1*f1)/maxf
			segments(bins1,-0.05,bins1,length.seg*f1-0.05)
			segments(bins0,-0.05,bins0,length.seg*-f0-0.05)
			lines(c(min(bins0,bins1)-0.01,max(bins0,bins1)+0.01),c(-0.05,-0.05))
			text(max(bins0,bins1)+dist.label,-0.01,d1lab,cex=cex.d01)
			text(max(bins0,bins1)+dist.label,-0.08,d0lab,cex=cex.d01)
					}
			}
	stats
	print(means)
	print(prop)
}
