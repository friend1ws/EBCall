library(VGAM)

lambda <- 0.5;

marginalLikelihood <- function(params, data) {
		
	alpha <- params[1];
	beta <- params[2];
	vec <-data;
	dlen <- length(vec);
	
	As <- vec[seq(1, dlen, 2)];
	Bs <- vec[seq(2, dlen, 2)];
	
	# print(As);
	# print(Bs);
	
	ML <- 0;
	for (i in 1:length(As)) {
		ML <- ML + ( lgamma(As[i] + Bs[i] + 1) - lgamma(As[i] + 1) - lgamma(Bs[i] + 1));
		ML <- ML - ( lgamma(alpha + beta + As[i] + Bs[i]) - lgamma(alpha + As[i]) - lgamma(beta + Bs[i]) );
		ML <- ML + ( lgamma(alpha + beta) - lgamma(alpha) - lgamma(beta) );
	}
	
	ML <- ML - lambda * log(alpha + beta);

	return(-ML);	
	
}




my_trans <- function(x = 0) {
    return( max(0, -log10(x)));
}


bdata <- read.table(commandArgs()[5], sep="\t",header=F);

infoData <- bdata[,1:4];
tumData <- bdata[,5];
norData <- bdata[,6];
refData <- bdata[,7:ncol(bdata)];

positiveList <- rep(0, nrow(bdata));
resultList <- matrix(0, nrow(bdata), 16);

for (i in 1:nrow(infoData)) {

    # print(bdata[i,]);
	
	tumBases <- as.integer(unlist(strsplit(as.character(tumData[i]), ",")));
	norBases <- as.integer(unlist(strsplit(as.character(norData[i]), ",")));
	
	if (pbinom(norBases[2] + norBases[4], norBases[1] + norBases[3], 0.5) < 0.01) {
		
		derror_p <- rep(0, 2 * ncol(refData));
		derror_n <- rep(0, 2 * ncol(refData));
		
		for (j in 1:ncol(refData)) {
			
			refBases <- as.integer(unlist(strsplit(as.character(refData[i, j]), ",")));
			
			derror_p[2 * j - 1] <- refBases[2];
			derror_p[2 * j] <- refBases[1] - refBases[2];
			
			derror_n[2 * j - 1] <- refBases[4];
			derror_n[2 * j] <- refBases[3] - refBases[4];
			
						
		}

        # print(derror_p);
        # print(derror_n);


        alpha_p <- 0.1;
        beta_p <- 1;
        if (sum(derror_p) > 0) {		
		    res_p <- constrOptim(c(20, 20), marginalLikelihood, grad=NULL, ui=matrix(c(1, 0, 1, 0, 1, 1), 3, 2), ci=c(0.1, 1, 1), data=derror_p);
		    alpha_p <- res_p$par[1];
            beta_p <- res_p$par[2];
        }

        alpha_n <- 0.1;
        beta_n <- 1;
        if (sum(derror_n) > 0) {
            res_n <- constrOptim(c(20, 20), marginalLikelihood, grad=NULL, ui=matrix(c(1, 0, 1, 0, 1, 1), 3, 2), ci=c(0.1, 1, 1), data=derror_n);
		    alpha_n <- res_n$par[1];
		    beta_n <- res_n$par[2];
	    }	
		
		lpv_p <- 0;
		lpv_n <- 0;
		
        if (tumBases[2] > 0) {
            pv_p <- pbetabinom.ab(tumBases[2] - 1, tumBases[1], alpha_p, beta_p);
            if (1 - pv_p < 1e-60) {
                lpv_p <- 60;
            } else {
                lpv_p <- -log10(1 - pv_p);
            }
        }
        if (tumBases[4] > 0) {  
            pv_n <- pbetabinom.ab(tumBases[4] - 1, tumBases[3], alpha_n, beta_n);
            if (1 - pv_n < 1e-60) {
                lpv_n <- 60;
            } else {
                lpv_n <- -log10(1 - pv_n);
            }
        }


		misRatio_tum <- (tumBases[2] + tumBases[4]) / (tumBases[1] + tumBases[3]);
		misRatio_nor <- (norBases[2] + norBases[4]) / (norBases[1] + norBases[3]);

        ref_tum <- (tumBases[1] + tumBases[3]) - (tumBases[2] + tumBases[4]);
        alt_tum <- tumBases[2] + tumBases[4];
        ref_nor <- (norBases[1] + norBases[3]) - (norBases[2] + norBases[4]);
        alt_nor <- norBases[2] + norBases[4]; 
   
        Fmatrix <- matrix(as.integer(c(alt_tum, alt_nor, ref_tum, ref_nor)), 2, 2);
        fisPV <- max(0, -log10(fisher.test(Fmatrix)$p.value));
		
		strandRatio_tum <- tumBases[2] / (tumBases[2] + tumBases[4]);
		strandRatio_nor <- "---"
		if (norBases[2] + norBases[4] > 0) {
			strandRatio_nor <- norBases[2] / (norBases[2] + norBases[4]);
		}
				
		fsPv <- - log10(1 - pchisq(2 * lpv_p / log10(exp(1)) + 2 * lpv_n / log10(exp(1)), 4));

        if (fsPv == Inf) {
            fsPv <- 100;
        }

		if (fsPv > -log10(0.01)) {
			positiveList[i] <- 1
		}
		
		
		# print(c(misRatio_tum, misRatio_nor, strandRatio_tum, strandRatio_nor, fsPv, lpv_p, lpv_n, fisPV, alpha_p, beta_p, alpha_n, beta_n));
		# resultList[i,] <- c(misRatio_tum, misRatio_nor, strandRatio_tum, strandRatio_nor, fsPv, lpv_p, lpv_n, fisPV, alpha_p, beta_p, alpha_n, beta_n);
        # print(c(misRatio_tum, strandRatio_tum, alt_tum, misRatio_nor, strandRatio_nor, alt_nor, fsPv, lpv_p, lpv_n, fisPV, alpha_p, beta_p, alpha_n, beta_n));
        # print(c(misRatio_tum, strandRatio_tum, ref_tum + alt_tum, alt_tum, misRatio_nor, strandRatio_nor, ref_nor + alt_nor, alt_nor, fsPv, lpv_p, lpv_n, fisPV, alpha_p, beta_p, alpha_n, beta_n));
        resultList[i,] <- c(misRatio_tum, strandRatio_tum, ref_tum + alt_tum, alt_tum, misRatio_nor, strandRatio_nor, ref_nor + alt_nor, alt_nor, fsPv, lpv_p, lpv_n, fisPV, alpha_p, beta_p, alpha_n, beta_n);

	}
	
}



write.table(cbind(infoData[positiveList == 1,,drop=FALSE], tumData[positiveList == 1,drop=FALSE], norData[positiveList == 1,drop=FALSE], resultList[positiveList == 1,,drop=FALSE]), commandArgs()[6], sep="\t", row.names=F, col.names=F);






