normalization <-
function(count,method="median"){

    METHODS = c('total','median','quantile','TMM');
    
    if (is.na(pmatch(method,METHODS))){
          stop("invalid normalization method", paste("", method));
    }
    else if (method=='median'){
          mediani = apply(count,2,median);
          deltai_1 = mediani/exp(mean(log(mediani)));
    }
    else if (method=='total'){  
          sumi = apply(count,2,sum);
          deltai_1 = sumi/mean(sumi);
    }
    else if (method=='quantile'){
          quantilei = apply(count,2,quantile)[4,];
          deltai_1 = quantilei/mean(quantilei);
    }
    else if (method=='TMM'){
          deltai_1 = edgeR::calcNormFactors(count);
    }
    count = round(sweep(count,2,deltai_1,"/"));
    return(list(count=count,sizeFactor=deltai_1));
 }
