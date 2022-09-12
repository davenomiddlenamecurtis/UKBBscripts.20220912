#!/share/apps/R-3.6.1/bin/Rscript
# find variants in sao file which might be individually significant

gene="HUWE1"

for (gg in 1:length(genes)){
  
  for (cc in 1:length(cohort_names)){
    gene=genes[gg]
    nValid=0
    nContDis=0
    nCaseDis=0
    nContDam=0
    nCaseDam=0
    prefix = sprintf(prefixSpec, cohort_names[cc])
    saoFile=sprintf("%s/%s.%s.sao",resultsFolder,prefix,gene)
    
    
    lines=readLines(saoFile)
    
    lResults=data.frame(matrix(NA, ncol = 7+length(AFs), nrow = length(lines))) #make sure to cut this down to nValid later on 
    colnames(lResults)=c(c("gene","chr","pos","DNAChange","aaChange","control","case"),AFs)
    locRow=0
    pResults=data.frame(matrix(NA, ncol = 8, nrow = length(lines)*maxTranscripts),stringsAsFactors=FALSE)
    colnames(pResults)=c("gene","transcript","amino_acid_position","amino_acid_change","control","case","chr","pos")
    
    r=1
    
    for (ll in 1:length(lines)) #add nValid = nValid +1 at the bottom of this loop
    {   #use ll=4 when testing as all conditions work 
      words=strsplit(lines[ll],"\\s+")[[1]] 
      if (length(words)!=17) next
      weight=as.numeric(words[16])
      if (weight < weight_threshold) next
      print(weight)
      control = as.numeric(words[4])
      case = as.numeric(words[10])
      if (case+control > maxHetsToUse) next
      
      current_line = lines[ll]
      split_string = unlist(strsplit(current_line, split=":")) 
      chromosome = split_string[1]
      
      
      if (grepl(".", split_string[2], fixed=TRUE)){
        pos = unlist(strsplit(split_string[2], "."))
        chromosome_position = pos[1]
      } else {
        chromosome_position = split_string[2]
      }
      locRow=locRow+1
      lResults$gene[locRow]=gene
      lResults$chr[locRow]=chromosome
      lResults$pos[locRow]=chromosome_position
      lResults$control[locRow]=control
      lResults$case[locRow]=case
      gnomadLine=getVCFLine(gnomADFile,sprintf("%s/%s",outputFolder,"tempGnomAD.vcf"),as.numeric(chromosome),as.numeric(chromosome_position))
      for (f in 1:length(AFs)) 
      {	
        if (gnomadLine=="") {
          lResults[locRow,7+f]=0
        } else {
          lResults[locRow,7+f]=sum(as.numeric(unlist(strsplit(getVCFEntry(gnomadLine,AFs[f]),",")))) # "0.2,0.15" for different alleles
        }
      }
      aWords=strsplit(current_line,"\\|")[[1]]
      t=-entriesPerTranscript
      while (t<length(aWords)-(18+entriesPerTranscript))
      {
        t=t+entriesPerTranscript
        aGene=aWords[t+4]
        if (gene!=aGene) next
        transcript=(aWords[t+7])
        amino_acid_position = (aWords[t+15]) # can be e.g. 1338-1339
        amino_acid_position = strsplit(amino_acid_position,"-")[[1]][1]
        amino_acid_change = aWords[t+16]
        DNA_change=aWords[t+17]
        
        
        if (amino_acid_change=="") next # e.g. upstream in this transcript
        print(gene)
        print(transcript)
        print(amino_acid_position)
        print(amino_acid_change)
        if (is.na(lResults$aaChange[locRow])) {
          lResults$aaChange[locRow]=amino_acid_change
          lResults$DNAChange[locRow]=DNA_change
        }
        pResults[r,1] = gene
        pResults[r,2] = transcript
        pResults[r,3] = amino_acid_position
        pResults[r,4] = amino_acid_change
        pResults[r,5] = control
        pResults[r,6] = case
        pResults[r,7] = as.numeric(chromosome)
        pResults[r,8] = as.numeric(chromosome_position)
        r = r+1
      }
      if (grepl("X",lResults$aaChange[locRow],fixed=TRUE) || grepl("*",lResults$aaChange[locRow],fixed=TRUE)) {
        nContDis=nContDis+control
        nCaseDis=nCaseDis+case
        r=r-1# do not use this row
      }else {
        nContDam=nContDam+control
        nCaseDam=nCaseDam+case
      }
      nValid = nValid + 1
    }
