pMimHTML = function(output,dir = 'mirpathways',cutoff = 0.05,filename = 'pMimResults.html',GeneSymbol = NULL,outputHTML=FALSE){
  require(googleVis)
  results = output$Results
  results = results[results[,'Score'] <= cutoff,]
  
  results = data.frame(results, miRNA.DE.pvalue = 2*pnorm(-abs(output$Zmi[results[,'miRNA']])),genes = rep(NA,dim(results)[1]))
  results = results[,c(1,2,6,3,5,4)]
  colnames(results) = c('microRNA','Pathway','Genes','Direction','microRNA DE p-value','Intergrative Score')
  
  dir.create(dir)
  gene.files = NULL
  for(i in 1:dim(results)[1]){
    miRNA = as.character(results[i,1])
    path = as.character(results[i,2])
    genes = intersect(output$targets[[miRNA]],output$pathways[[path]])
    Correlation = output$cor[miRNA,genes]
    corP = pnorm(output$corTransform[miRNA,genes]*sign(output$Zmi[miRNA]))
    DE = pnorm(output$Zg[genes]*sign(output$Zmi[miRNA]))
    if(!is.null(GeneSymbol))genes = GeneSymbol[genes]
    res = data.frame(Genes = genes,'One Sided DE P-value' = DE,Correlation = Correlation,'One Sided Correlation P-value' = corP)
    form = as.list(rep('#.###',4))
    names(form) = colnames(res)
    file = paste(miRNA,'-',path,'.html',sep = '')
    file = gsub('-','_',file)
    file = paste(dir,file,sep= '/')
    gene.files[i] = file
    res = res[order(res[,4]),]
    print(gvisTable(res,format = form,options = list(width = '100em')),type = 'html',file = file)
  }
  
  
  results[,'Genes'] = paste('<a href=','"',gene.files,'" target="_blank">','genes','</a>',sep="")
  
  library(googleVis)
  form = as.list(rep('#.###',6))
  names(form) = colnames(results)
  results = as.data.frame(results)
  print(gvisTable(results,format = form,options = list(width = '100em')),type = 'html',file = filename)
  
  if(outputHTML == TRUE)gvisTable(as.data.frame(results),format = form,options = list(width = '100em'))
}
