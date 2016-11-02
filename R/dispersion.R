dispersion <-
function(M,condition,matched,method="edgeR"){

    METHODS = c("edgeR","DESeq2");
    
    count = M$count;
    deltai_1 = M$sizeFactor;
    
    if (is.na(pmatch(method,METHODS))){
        stop("invalid normalization method", paste("", method));
    }
    else if (method=='edgeR'){
        if (!matched){
            # dispersion parameter estimate          
            condition = as.factor(condition);
            y = edgeR::DGEList(count=count, group=condition);
            y = edgeR::calcNormFactors(y);
            design = model.matrix(~condition);
            y = edgeR::estimateGLMCommonDisp(y, design);
            y = edgeR::estimateGLMTagwiseDisp(y, design);
            phyg_0 = y$common.dispersion;
            phyg_1 = y$tagwise.dispersion;
        }else{
            # dispersion parameter estimate
            R = length(condition[!duplicated(condition)]);
            I = dim(count)[2]/R;
            phyg_1 = phyg_0 = 0;
            condition0 = as.factor(rep(1,I));
            for (d in 1:R){
                count0 = count[,1:I+(d-1)*I];
                y = edgeR::DGEList(count=count0, group=condition0);
                y = edgeR::calcNormFactors(y);
                y = edgeR::estimateGLMCommonDisp(y);
                y = edgeR::estimateGLMTagwiseDisp(y);
                phyg_1 = phyg_1+y$tagwise.dispersion;
                phyg_0 = phyg_0+y$common.dispersion;
            } 
            phyg_1 = phyg_1/R;
            phyg_0 = phyg_0/R;
        }
    }else if (method=='DESeq2'){  
        if (!matched){
            condition = as.factor(condition);    
            colData = S4Vectors::DataFrame(condition=condition);
            d = DESeq2::DESeqDataSetFromMatrix(count, colData, design=~condition);
    
            d = DESeq2::estimateSizeFactors(d);
            d = DESeq2::estimateDispersions(d);  
        }else{
            R = length(condition[!duplicated(condition)]);
            I = dim(count)[2]/R;
            condition = as.factor(condition);
            sample = as.factor(rep(1:I,R));
    
            colData = S4Vectors::DataFrame(condition=condition,sample=sample)
            d = DESeq2::DESeqDataSetFromMatrix(count, colData, formula(~ condition+sample))
    
            d = DESeq2::estimateSizeFactors(d);
            d = DESeq2::estimateDispersions(d);
        }
        phyg_1 = DESeq2::dispersions(d);
        phyg_0 = mean(phyg_1);
    }
    return(list(count=count,sizeFactor=deltai_1,condition=condition,commonDispersion=phyg_0,tagwiseDispersion=phyg_1,matched=matched));
}
