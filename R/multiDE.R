multiDE <-
function(M,n.top=nrow(M$count)){
  
  condition = M$condition;
  R = length(condition[!duplicated(condition)]);
  is.matched = M$matched;  
  count = M$count;
  deltai_1 = M$sizeFactor;
    
  if(!is.matched){    
    count = round(count);
    G = nrow(count);
    
    phyg_1=M$tagwiseDispersion;
 
    # deal with zero count
    nd = as.numeric(table(condition[duplicated(condition)])+1);  
    for (g in 1:G){
      for (d in 1:R){
        if (mean(as.numeric(count[g,(1+sum(nd[1:d])-nd[d]):sum(nd[1:d])]))==0 & var(as.numeric(count[g,(1+sum(nd[1:d])-nd[d]):sum(nd[1:d])]))==0){
          count[g,]=count[g,]+1;
        }
      }
    }
    
    # model parameter estimate
    muldg_1 = matrix(NA,R,G);
    for (d in 1:R){
      muldg_1[d,] = count%*%c(rep(0,sum(nd[1:d])-nd[d]),rep(1/nd[d],nd[d]),rep(0,sum(nd[1:R])-sum(nd[1:d])));
    }                                                                             
    sigmadg_1 = sqrt(sweep(muldg_1^2,2,phyg_1,"*"));
    eitadg_1 = log(muldg_1);
    mul_1 = sum(sweep(eitadg_1,1,nd,"*"))/(sum(nd)*G);
    alphad_1 = rowSums(eitadg_1)/G-mul_1;
    betag_1 = (t(eitadg_1)%*%nd)/(sum(nd))-mul_1;
    gammadg_1 = sweep(sweep(eitadg_1-mul_1,1,alphad_1,"-"),2,betag_1,"-");
                
    # ud and vg iteration estimate
    
    condition = as.factor(condition);
    y = edgeR::DGEList(count=count, group=condition);
    y = edgeR::calcNormFactors(y);
    design = model.matrix(~condition);
    y = edgeR::estimateGLMCommonDisp(y, design);
    y = edgeR::estimateGLMTagwiseDisp(y, design);
    fit = edgeR::glmFit(y, design);
    result = edgeR::glmLRT(fit, coef=2:R);
    p.value = result$table$PValue;
    top.id = order(p.value)[1:n.top];
    
    ud_1 = c(1,rep(-1/(R-1),R-1));
    vg_1 = rep(0,G);
    k = 0;
    repeat{
      vg_2 = vg_1;
      vg_1 = (t(gammadg_1)%*%(nd*ud_1))/sum(nd*ud_1^2);
      ud_2 = ud_1;
      ud_1 = (gammadg_1[,top.id]%*%vg_1[top.id])/sum(vg_1[top.id]^2);
      ud_1 = ud_1/ud_1[1];
      k = k+1;
      if (sum(abs(vg_1-vg_2))<0.001 & sum(abs(ud_1-ud_2))<0.001)break
    }
    
    # paramater variance estimate
    varvg = as.vector(t(nd*ud_1^2)%*%sweep(muldg_1^(-1),2,phyg_1,"+"))/sum(nd*ud_1^2)^2;
    
    # perform wald test and generate pvalue
    tg = vg_1/sqrt(varvg);
    pvalue = rep(NA,G);
    pvalue = 2*pnorm(-abs(tg));
    log2fold.change = outer(as.numeric(vg_1),as.numeric(ud_1-ud_1[1]))/log(2);
    sd.log2foldchange = outer(as.numeric(sqrt(varvg)),abs(as.numeric(ud_1-ud_1[1])))/log(2);
      
    return(list(p.value=pvalue,log2FoldChange=log2fold.change,u=ud_1,sd.log2FoldChange=sd.log2foldchange));
  }else{
    if(ncol(count)%%R!=0){
      warning('The total sample number should be the number of condition times a positive integer.')
    }else{
      dm = dim(count);
      I = dm[2]/R;
      count = round(count);
      G = nrow(count);
      
      phyg_0=M$commonDispersion; phyg_1=M$tagwiseDispersion;
      
      # deal with zero count
      for (g in 1:G){
        for (d in 1:R){
          if (mean(as.numeric(count[g,(1+I*(R-d)):(I*(R-d+1))]))==0 & var(as.numeric(count[g,(1+I*(R-d)):(I*(R-d+1))]))==0){
            count[g,]=count[g,]+1;
          }
        }
      }
      
      # model parameter estimate
      muldg_1 = matrix(NA,R,G);
      for (d in 1:R){
        muldg_1[d,] = count%*%c(rep(0,I*(d-1)),rep(1/I,I),rep(0,I*(R-d)));
      }
      sigmadg_0 = sqrt(sweep(muldg_1^2,2,phyg_0,"*"));
      eitadg_1 = log(muldg_1);
      mul_1 = sum(I*eitadg_1)/(I*R*G);
      alphad_1 = rowSums(eitadg_1)/G-mul_1;
      betag_1 = colSums(I*eitadg_1)/(I*R)-mul_1;
      gammadg_1 = sweep(sweep(eitadg_1-mul_1,1,alphad_1,"-"),2,betag_1,"-");
      
      # ud and vg iteration estimate
      
      condition = as.factor(condition);
      sample = as.factor(rep(1:I,R));
      y = edgeR::DGEList(count=count, group=condition);
      y = edgeR::calcNormFactors(y);
      design = model.matrix(~condition+sample);
      y = edgeR::estimateGLMCommonDisp(y, design);
      y = edgeR::estimateGLMTagwiseDisp(y,design);
      fit = edgeR::glmFit(y, design);
      result = edgeR::glmLRT(fit, coef=2:R);
      p.value = result$table$PValue;
      top.id = order(p.value)[1:n.top];
      
      ud_1 = c(1,rep(-1/(R-1),R-1));
      vg_1 = rep(0,G);
      k = 0;
      repeat {
        vg_2 = vg_1;
        vg_1 = (t(gammadg_1)%*%(I*ud_1))/sum(I*ud_1^2);
        ud_2 = ud_1;
        ud_1 = (gammadg_1[,top.id]%*%vg_1[top.id])/sum(vg_1[top.id]^2);
        ud_1 = ud_1/ud_1[1];
        k = k+1;
        if (sum(abs(vg_1-vg_2))<0.001 & sum(abs(ud_1-ud_2))<0.001)break
      }
    
      # correlations estimate
      dat = array(NA,c(G,I,R));
      for (d in 1:R){
        dat[,,d] = count[,1:I+(d-1)*I];
      }
      correYg = array(NA,c(G,R,R));
      for (g in 1:G){
        x = dat[g,,];
        v = apply(x,2,var);
        if (any(v==0)){
          correYg[g,,] = diag(R);
        }else{
          correYg[g,,] = cor(x);
        }
      }
      n = I;
      tmp = exp(lgamma((n-1)/2)-lgamma(1/2)-lgamma((n-2)/2));
      
      cor.update = function(r){
        n.r = length(r);
        res = rep(1/tmp,n.r);
        for(i in 1:n.r){
          r0 = r[i];
          if(!r0 %in% c(0,1,-1)){
            f.integral <- function(t){
              return(t^{-1/2}*(1-t)^{(n-2)/2-1}/(1-t*(1-r0^2))^{1/2})
            }
            res[i] = integrate(f.integral,0,1)$value;
          }
        }
        res = tmp * res * r;
        return(res);
      }

      for (d1 in 1:R){
        for (d2 in d1:R){
          correYg[,d1,d2] = cor.update(correYg[,d1,d2]);
          correYg[,d2,d1] = correYg[,d1,d2]
        }
      }

      correYg_1 = apply(correYg,2:3,mean);
      correZ_1 = matrix(1,R,R);
      for (d1 in 1:(R-1)){
        for (d2 in (d1+1):R){
          czg = correYg_1[d1,d2]*sqrt(muldg_1[d1,]+sigmadg_0[d1,]^2)*sqrt(muldg_1[d2,]+sigmadg_0[d2,]^2)/sigmadg_0[d1,]/sigmadg_0[d2,];
          correZ_1[d1,d2] = mean(czg);
          correZ_1[d2,d1] = correZ_1[d1,d2];
        }
      }

      # paramater variance estimate
      covvg = 0;
      for (d1 in 1:(R-1)){
        for (d2 in (d1+1):R){
          covvg = covvg+ud_1[d1]*ud_1[d2]*phyg_1*correZ_1[d1,d2];   
        }
      }
      varvg = (as.vector(t(ud_1^2)%*%sweep(muldg_1^(-1),2,phyg_1,"+"))+2*covvg)/I/sum(ud_1^2)^2;
      # perform wald test and generate pvalue
      tg = vg_1/sqrt(varvg);
      pvalue = 2*pnorm(-abs(tg));
      
      log2fold.change = outer(as.numeric(vg_1),as.numeric(ud_1-ud_1[1]))/log(2);
      sd.log2foldchange = outer(as.numeric(sqrt(varvg)),abs(as.numeric(ud_1-ud_1[1])))/log(2);
      
      return(list(p.value=pvalue,log2FoldChange=log2fold.change,u=ud_1,sd.log2FoldChange=sd.log2foldchange));
    }
  }
}
