#last update 20160831

###	1. process_subVCF_May2016_function													  40
###	2. simple.coverageARRAY.creator_function20150719										 200
###	3. Additional_elements_for_Cross_Validation_function__June2016						 400
###	4. gamma_cross_validation_procedure__function__May2016								 470
###	5. New.ACGTZ.array.creator__with.gamma.validation__function__May2016					 800
###	6. COLLAPSE_ACGTZ_ARRAYS__GAMMA__function_May2016									1100
###	7. Conservative_method__COLLAPSE_ACGTZ_ARRAYS__GAMMA__function__May2016				1500		
###	8. Site_elimiantion_method__COLLAPSE_ACGTZ_ARRAYS__GAMMA__function__May2016			1770
###	9. Fixed_Thresholds_method__COLLAPSE_ACGTZ_ARRAYS__GAMMA__function__May2016			2000
###	10. Improved_GAMMA_Create_ARRAY_collapsed_short_info_ACGTZ_function20150725			2240
###	11. gamma_validation_v_two__function_20150723										2320
###	12. GAMMA_collapse_acgtz_arrays_function20150725										2670
###	13. New_INFO_for_Collapse_acgtz_arrays__function_20160601							2770
###	14. Collapse_array__Conservative_method												2880
###	15. Collapse_array__Site_elimination_method											2960	
###	16. Collapse_array__Fixed_Thresholds_method											3080	
###	17. COLLAPSE_InDel_MatrixNA__function__May2016										3300	
###	18. Run.GapJumper.five.methods  				 [all prarametes]						3420
###	19. GapJumper.plus.all.methods														4650 
###	20. Run.Probabilistic.method															4780 
###	21. Run.Probabilistic.beta.method													4930
###	22. GapJumper.plus.cons.method														5060
###	23. GapJumper.plus.SiteElim.method													5220 
###	24. GapJumper.plus.FixedTr.method													5350
###	25. Calculate.Noise.threshold														5500	
###	26. Test.Noise.threshold															5570
###	27. Remove.Noise.Global				- old function . no use, but I kept it in case	5630
###	28. Remove.Noise.using.Noise.Threshold												5750
###	29. Remove.Noise.using.Acceptance.Threshold											6000
###	30. fuse.consensus																	6200
###	31. Compare.samples.using.multiSNV.files												6400
###	32. Calculate.Acc.using.Ns															6800





###	---------------------------------------------------------------------------------------------------
###	1. process_subVCF_May2016_function		40
###	---------------------------------------------------------------------------------------------------
process_subVCF_May2016_function<-function(
							In.directory	,
							Res.directory,	
							subVCF.list,
							RADloci.custom.names,
							RADloci.to.use,	
							RADloci.to.remove,
							Set.name	){

	#  1. info:
	#	is not here anymore
		
	#  2. load one file and prepare rad names:		
	setwd(In.directory)
	any.subVCF		<-read.table(subVCF.list[1], sep="\t",as.is=T,strip.white=T,colClasses="character") 
 	any.subVCF		<-as.matrix(any.subVCF)  				# it has loci, but not headers	
	rad.nr			<-1:nrow(any.subVCF)
	
	#  3. new matrices
	cluster_matrix	<-matrix(ncol=length(subVCF.list),nrow=nrow(any.subVCF)); colnames(cluster_matrix)<-subVCF.list
	repeat_matrix	<-matrix(ncol=length(subVCF.list),nrow=nrow(any.subVCF)); colnames(repeat_matrix)<-subVCF.list
	coding_matrix	<-matrix(ncol=length(subVCF.list),nrow=nrow(any.subVCF)); colnames(coding_matrix)<-subVCF.list
	snp_type_matrix	<-matrix(ncol=length(subVCF.list),nrow=nrow(any.subVCF)); colnames(snp_type_matrix)<-subVCF.list
	snp_matrix		<-matrix(ncol=length(subVCF.list),nrow=nrow(any.subVCF)); colnames(snp_matrix)<-subVCF.list
	frequency_matrix<-matrix(ncol=length(subVCF.list),nrow=nrow(any.subVCF)); colnames(frequency_matrix)<-subVCF.list
	coverage_matrix	<-matrix(ncol=length(subVCF.list),nrow=nrow(any.subVCF)); colnames(coverage_matrix)<-subVCF.list
	indel_matrix		<-matrix(ncol=length(subVCF.list),nrow=nrow(any.subVCF)); colnames(indel_matrix)<-subVCF.list
	###
	sum(as.numeric(RADloci.to.remove))
	
	#  4. columns to take (standard.col.nr.for.subVCF.v6):
	snp.type.col.nr		<-10
	snp.col.nr			<-11
	freq.col.nr			<-12
	gen.cov.col.nr		<-13
	
	#  5. info:
	cat("---","        ","\n")
	cat("*--","        combining data from separate subVCF files","\n")
	cat("---","        according to order and names in a table with sample names (1st column),","\n")
	cat("*--","        ",length(subVCF.list),"samples to combine","\n")
	cat("-->","        ")
	
	#  6. filling in:
	for (i in 1:length(subVCF.list)){
	  #i<-1	
	  subVCF_file			<- read.table(subVCF.list[i], sep="\t",as.is=T,strip.white=T,colClasses="character")
	  subVCF_file			<- as.matrix(subVCF_file)
	  snp_type_matrix[,i]	<-subVCF_file[,snp.type.col.nr]
	  snp_matrix[,i]			<-subVCF_file[,snp.col.nr] 
	  frequency_matrix[,i]	<-subVCF_file[,freq.col.nr]
	  coverage_matrix[,i]	<-subVCF_file[,gen.cov.col.nr]
	  cat(i,".  ")
	 };cat(" - ","end","\n");gc()

	#  7. indel matrix (indel==1):
	cat("---","        ","\n")
	cat("*--","        localizing indels:","\n")
	cat("-->","        ")
	for (col_nr in 1:ncol(snp_type_matrix)){
	  processing_vector		<- snp_type_matrix[,col_nr]
	  find_I					<- as.numeric(grepl("I",processing_vector ))
	  find_D					<- as.numeric(grepl("D",processing_vector ))
	  find_M					<- as.numeric(grepl("M",processing_vector ))
	  indel_matrix[,col_nr]	<- as.numeric((find_I + find_D + find_M)>0)
	cat(col_nr,".  ")};cat(" - ","end","\n")
	
	#  8. noindel matrix (snp==1):
	noindel_matrix<-indel_matrix; noindel_matrix[,]	<-0
	noindel_matrix[which(indel_matrix[,]==0)]		<-1

	#  9. final elements:
	ProcSubVCF_matrices <-c("cluster_matrix", "repeat_matrix", "coding_matrix", "snp_type_matrix", "snp_matrix", "frequency_matrix", "coverage_matrix", "indel_matrix")
	 cluster_matrixNA			<- cluster_matrix
	 repeat_matrixNA				<- repeat_matrix
	 coding_matrixNA				<- coding_matrix
	 snp_type_matrixNA			<- snp_type_matrix
	 snp_matrixNA				<- snp_matrix
	 frequency_matrixNA			<- frequency_matrix
	 coverage_matrixNA			<- coverage_matrix
	 indel_matrixNA				<- indel_matrix
	 noindel_matrixNA			<- noindel_matrix	
	rm(list=ProcSubVCF_matrices);gc()
	
	# 10. NA problem:
	cat("---","        ","\n")
	cat("*--","        removing potential problems with na:","\n")
 	cat("---","        NA instead of na","\n")
 	cluster_matrixNA[cluster_matrixNA=="na"]			<-NA
 	repeat_matrixNA[repeat_matrixNA=="na"]			<-NA
 	coding_matrixNA[coding_matrixNA=="na"]			<-NA
 	snp_type_matrixNA[snp_type_matrixNA=="na"]		<-NA
 	snp_matrixNA[snp_matrixNA=="na"]					<-NA
 	frequency_matrixNA[frequency_matrixNA=="na"]		<-NA
 	coverage_matrixNA[coverage_matrixNA=="na"]		<-NA						
	
	# 11. rownames:
	rad.name							<-RADloci.custom.names
 	rownames(cluster_matrixNA)		<-rad.name
 	rownames(repeat_matrixNA)		<-rad.name
 	rownames(coding_matrixNA)		<-rad.name
 	rownames(snp_type_matrixNA)		<-rad.name
 	rownames(snp_matrixNA)			<-rad.name
 	rownames(frequency_matrixNA)		<-rad.name
 	rownames(coverage_matrixNA)		<-rad.name
 	rownames(indel_matrixNA)		<-rad.name
 	rownames(noindel_matrixNA)		<-rad.name

	# 12. list:
	matrixNA_LIST	<- list(cluster_matrixNA, repeat_matrixNA, coding_matrixNA, snp_type_matrixNA, snp_matrixNA, frequency_matrixNA, coverage_matrixNA, indel_matrixNA, noindel_matrixNA)
	matrix_LIST_names	<-c("cluster_matrixNA", "repeat_matrixNA", "coding_matrixNA", "snp_type_matrixNA", "snp_matrixNA", "frequency_matrixNA", "coverage_matrixNA", "indel_matrixNA", "noindel_matrixNA")	
	names(matrixNA_LIST)	<- matrix_LIST_names		
	
	# 13. cut big matrix into smaller pieces
	NotSelected_matrixNA_LIST		<-matrixNA_LIST
	take.loci						<-as.numeric(RADloci.to.use)
	take.loci[which(as.numeric(RADloci.to.remove)[]>0)]<-0
	#
	RAD.loci.to.select			<-take.loci
	selected_matrixNA_LIST		<-NotSelected_matrixNA_LIST
	l.nr<-1
	for(l.nr in 1:length(NotSelected_matrixNA_LIST)){
			selected_matrixNA_LIST[[l.nr]]<-NotSelected_matrixNA_LIST[[l.nr]][which(RAD.loci.to.select[]==1),,drop=FALSE]}
			dim(selected_matrixNA_LIST[[l.nr]])
	matrixNA_LIST				<-selected_matrixNA_LIST
	
	# 13. save:
	setwd(Res.directory)	
	save.name					<-paste(Set.name	,"__matrixNA_LIST.RData",sep=""); save.name
	save(matrixNA_LIST,file=save.name)
	cat("---","        ","\n")	
	cat("*--","        matrixNA_LIST.RData  - saved","\n")
	cat("---","        with:","\n")
	cat("---","        ",nrow(matrixNA_LIST[[1]]),"positions","\n")
	cat("---","        ",ncol(matrixNA_LIST[[1]]),"samples","\n")
	cat("---","        ","\n")
	
	# 14. info:
	cat("***","        Done","\n")
	cat("***","        Now waith for more fun","\n")
	gc()
}















###	---------------------------------------------------------------------------------------------------
###	2. simple.coverageARRAY.creator_function20150719	(May 2016)
###	---------------------------------------------------------------------------------------------------

simple.coverageARRAY.creator_function20150719<-function(
	data_list,
	In.directory,
	Res.directory){

cat(".........................................................","\n")
cat("\n");cat(length(data_list)," - elements to process:","\n");cat("\n");
cat(".........................................................","\n")
DATA_LIST_NUMBER <-1
for(DATA_LIST_NUMBER in 1:length(data_list)){
		
		#...load
		setwd(In.directory); 
		load(data_list[DATA_LIST_NUMBER])		
		
		#...prefix.name	
	 	temp_Data_set_name		<-data_list[DATA_LIST_NUMBER]
	 	temp_Data_set_name		<-strsplit(temp_Data_set_name,split="_matrixNA_LIST.RData")
	 	temp_Data_set_name		<-temp_Data_set_name[[1]]
	 	temp_Data_set_name		<-temp_Data_set_name[1]	
		
		#...array.data	
		snp				<-matrixNA_LIST[[5]]
		snp.type			<-matrixNA_LIST[[4]]
		snp.freq			<-matrixNA_LIST[[6]]
		ind.pos			<-matrixNA_LIST[[8]]
		snp.cov			<-matrixNA_LIST[[7]]	
		snp.row.clust	<-as.numeric(((matrixNA_LIST[[1]])[,1])>0) #every ONE must be eliminated	 
	
		#...display:
		cat("\n");
		cat(".........................................................","\n")				
		cat("[",DATA_LIST_NUMBER,"] - ",temp_Data_set_name,"\n");cat("[...] - ",date(),"\n")	

				
		#...arrays and total frequency
		 array.names		<-c("A", "C", "G", "T","snp.filtered.cov","gen.cov","ind.cov","poly")		  
		 Empty_matrix		<-snp
		 Empty_matrix[,]	<-0
		 mode(Empty_matrix)	<-"numeric"
		 RADloci_names		<-rownames(snp)
		 Sample_names		<-colnames(snp)	 	 
		 array_matrices		<-c(rep(Empty_matrix,length(array.names)))
		 array_dim			<-c(dim(snp)[1],dim(snp)[2],length(array.names))
		 array_dim_names	<-list(rownames(snp),colnames(snp),array.names)	 	
		 ARRAY_CovACGT		<-array(array_matrices,dim=array_dim)
		 dimnames(ARRAY_CovACGT)<-array_dim_names	
		 mode(ARRAY_CovACGT)<-"numeric"	 		

		#sum frequency - only for indels and polyallelic positions - all the rest has zero value
		 total.freq	<-Empty_matrix
		 for(i in 1:ncol(snp.freq)){
		  #i<-1
		  row <-snp.freq[,i]	
		  row <-strsplit(row,split=",")	 
		  s.row <-1:length(row)
		  for(ii in 1:length(row)){
		   #ii<-1	
		   one.el <-row[[ii]]	
		   mode(one.el)<-"numeric"
		   one.el<-sum(one.el)
		   s.row[ii]<-one.el	
		   }#s.row
		  total.freq[,i]<-s.row	   	
		  }#end for. # aa<-total.freq; colnames(aa)<-1:ncol(aa); rownames(aa)<-1:nrow(aa); aa
		#######################################snp.total.freq.for.perc.calc
		snp.rel.cov<-snp.cov
		snp.rel.cov[which(total.freq[,]>0)]<-total.freq[which(total.freq[,]>0)]
		mode(snp.rel.cov)<-"numeric"
		# aa<-snp.rel.cov; colnames(aa)<-1:ncol(aa); rownames(aa)<-1:nrow(aa); aa
		
################################################# fast elements
#i - col.nr; ii- row.nr; ll- allele number in each matrixNA[ii,i]
ARRAY_CovACGT[,,5] <-snp.rel.cov	#dimnames(ARRAY_CovACGT)
ARRAY_CovACGT[,,6] <-snp.cov		#dimnames(ARRAY_CovACGT)
snp.cov<-snp.rel.cov; rm(snp.rel.cov); mode(snp.cov)<-"numeric"

#################################################i<-1; ii<-1; ll<-1
for(i in 1:ncol(snp)){
cat("matr,")
#take.one.row(i)-only 4 variables (snp, snp.type, snp.freq, snp.cov-newly done)	  
snp.row			<-snp[,i]
snp.type.row	<-snp.type[,i]
snp.freq.row 	<-snp.freq[,i]
snp.cov.row		<-as.numeric(snp.cov[,i])
#eliminate.NA
NA.snp.row<-is.na(snp.row)
snp.row[which(NA.snp.row[]==TRUE)]<-"empty"   	   	   
snp.type.row[which(snp.type.row[]=="")]<-"empty"
snp.freq.row[which(snp.freq.row[]=="")]<-0 	 
#split.allele.information
snp.row 		<-strsplit(snp.row,split=",")
snp.type.row	<-strsplit(snp.type.row,split=",")
snp.freq.row	<-strsplit(snp.freq.row,split=",")
	  	  	      
for(ii in 1:length(snp.type.row)){
#.el ==one.locus wth all snp's in it.	
snp.el			<-snp.row[[ii]]
snp.type.el		<-snp.type.row[[ii]]
snp.freq.el 	<-as.numeric(snp.freq.row[[ii]])	   	
#el.from.vector
snp.cov.el		<-snp.cov.row[ii]							
#poly.or.not-poly ->ARRAY_CovACGT[ii,i,8]
if(length(snp.type.el)>1){ 
ARRAY_CovACGT[ii,i,8]<-1
}else{
ARRAY_CovACGT[ii,i,8]<-0
}		 		 		
#snp and indels filtered coverage; dimnames(ARRAY_CovACGT) #ll<-1
for(ll in 1:length(snp.type.el)){	
		snp.all			<-snp.el[ll]
		snp.type.all	<-snp.type.el[ll]
		snp.freq.all 	<-as.numeric(snp.freq.el[ll])	   	 
		snp.cov.all		<-as.numeric(snp.cov.el)	  
	  	#snp.freq.all==snp.cov.all
		if(grepl("empty",snp.type.all)){		 			 	
			#when empty - there is no freq and no total.freq #snp.type.all<-"empty"			 	
		 	snp.el.freq<-100		 			 			 	
		 	if(grepl("A",snp.all)){
		 		ARRAY_CovACGT[ii,i,1]<-snp.cov.all				 				 				 	
		 	}else if (grepl("C",snp.all)){
		 		ARRAY_CovACGT[ii,i,2]<-snp.cov.all				 				 			 		
		 	}else if (grepl("G",snp.all)){
		 		ARRAY_CovACGT[ii,i,3]<-snp.cov.all				 				 					
		 	}else if (grepl("T",snp.all)){
		 		ARRAY_CovACGT[ii,i,4]<-snp.cov.all				 				 				 				 			
			}else{ 
				ARRAY_CovACGT[ii,i,5]<-snp.cov.all
			} #ARRAY_CovACGT[ii,i,] 
		#here snp.freq.all is individual		
		}else if(grepl("S",snp.type.all)){		 			 		 	
		 	if(grepl("A",snp.all)){
		 		ARRAY_CovACGT[ii,i,1]<-snp.freq.all				 				 				 	
		 	}else if (grepl("C",snp.all)){
		 		ARRAY_CovACGT[ii,i,2]<-snp.freq.all				 				 			 		
		 	}else if (grepl("G",snp.all)){
		 		ARRAY_CovACGT[ii,i,3]<-snp.freq.all				 				 					
		 	}else{
		 		ARRAY_CovACGT[ii,i,4]<-snp.freq.all	
		 	} #ARRAY_CovACGT[ii,i,]
		}else if(grepl("R",snp.type.all)){		 			 			 	
		 	if(grepl("A",snp.all)){
		 		ARRAY_CovACGT[ii,i,1]<-snp.freq.all				 				 				 	
		 	}else if (grepl("C",snp.all)){
		 		ARRAY_CovACGT[ii,i,2]<-snp.freq.all				 				 			 		
		 	}else if (grepl("G",snp.all)){
		 		ARRAY_CovACGT[ii,i,3]<-snp.freq.all				 				 					
		 	}else{
		 		ARRAY_CovACGT[ii,i,4]<-snp.freq.all
		 	} #ARRAY_CovACGT[ii,i,]; dimnames(ARRAY_CovACGT)
		}else if(grepl("D",snp.type.all)){
			temp<-as.numeric(ARRAY_CovACGT[ii,i,7])+as.numeric(snp.freq.all); ARRAY_CovACGT[ii,i,7]<-temp
		}else if(grepl("I",snp.type.all)){
			temp<-as.numeric(ARRAY_CovACGT[ii,i,7])+as.numeric(snp.freq.all); ARRAY_CovACGT[ii,i,7]<-temp
		}else{
			temp<-as.numeric(ARRAY_CovACGT[ii,i,7])+as.numeric(snp.freq.all); ARRAY_CovACGT[ii,i,7]<-temp
		}#end if-else					
}#end(ll)		
}#end(ii)
}#end(i)

mode(ARRAY_CovACGT)<-"numeric"
#end all i and ii
#aa<-ARRAY_CovACGT[1:100,1:10,];mode(aa)<-"numeric" ;colnames(aa)<-1:10; rownames(aa)<-1:100; aa[10,,]
#locus<-82;sample<-8;snp.cov<-matrixNA_LIST[[7]];snp[locus,sample];snp.type[locus,sample];snp.freq[locus,sample];snp.cov[locus,sample]
	  	
setwd(Res.directory); dir()
ARRAY_CovACGT_name<-paste(temp_Data_set_name,"_ARRAY_CovACGT",".RData",sep=""); ARRAY_CovACGT_name
save(ARRAY_CovACGT,file=ARRAY_CovACGT_name)

#display:
cat("\n");cat(ARRAY_CovACGT_name, ">>> has been saved","\n")
cat("[...] - --------------------------------------------------------------------------","\n")
}#DATA_LIST_NUMBER
}
























###	---------------------------------------------------------------------------------------------------
###	3. Additional_elements_for_Cross_Validation_function__June2016	(May 2016)
###	---------------------------------------------------------------------------------------------------
Additional_elements_for_Cross_Validation_function__June2016<-function(In.directory,Res.directory,data.list,z.classifier){		
		set.nr<-1	
		for(set.nr in 1:length(data.list)){			
			cat("***",set.nr,"-->",data.list[set.nr],"\n" )	#set.nr<-42		
			setwd(In.directory); 
			load(data.list[set.nr]); 
			set.name		<-(strsplit(data.list[set.nr],split=".RData")[[1]])[1];set.name
			save.name	<-paste(set.name,"_additional_elements_for_CV.RData",sep="");save.name			
			# exctract information			
			is.list		<-colnames(ARRAY_CovACGT)
			loci.list	<-rownames(ARRAY_CovACGT)								
						rm(ARRAY_CovACGT);	 gc()			
			# 	- 1 - isolate.names.to.use = long list with only RG group names
			  	isolate.names.to.use<-rep(0,length(is.list))
			  	for(is.n in 1:length(is.list)){isolate.names.to.use[is.n]<-paste((strsplit(is.list[is.n],split="__")[[1]])[1],"__",sep="")}
			#	 - 2 - isolates.to.test	 = short list with unique RG group names as table result
			  	isolates.to.test	 <-table(isolate.names.to.use)
			
			# 	- 3 - not.repeated.loci	 = short list with unique RG group names as table result
			#	it is an old function adapted to my particular experiment.
			#	find.repeated.or.not.repeated.loci.in.matrix_function20150116<-function(loci.names,vector.with.index.numbers){
			#		#any.array<-ar.s.bin; vector.with.index.numbers<-1
			#		main_categories_table<-c("__code10","__code00","__code11","__code01")
			#		result.vector<-rep(0,length(loci.names))
			#		for(i in 1:length(vector.with.index.numbers)){
			#		temp<-as.numeric(grepl(main_categories_table[vector.with.index.numbers[i]],loci.names)); result.vector<-result.vector+temp}
			#		return(result.vector)}
			#	not.repeated.index.nr	<-c(1,2)
			#	not.repeated.loci		<-as.numeric(find.repeated.or.not.repeated.loci.in.matrix_function20150116(loci.list, not.repeated.index.nr)[]>0)
			#	sum(not.repeated.loci)
			#	sum(not.repeated.loci)
							
			# 	- 3 - bis 20160601
				not.repeated.loci		<-as.numeric(as.numeric(z.classifier)[]>0)
				
			# 	list
				additional_elements_for_CV<-list(isolate.names.to.use,isolates.to.test,not.repeated.loci)	
				names(additional_elements_for_CV)<-c("isolate.names.to.use","isolates.to.test","not.repeated.loci")					
			#	###											
			setwd(Res.directory); dir()	
			save(additional_elements_for_CV,file=save.name)	
			cat("---- ","---"," --> SAVED: ",save.name,"\n" );
			cat("----","with","-->",sum(not.repeated.loci),"not repeated loci were found","\n")		
			cat("----","and ","-->",sum(not.repeated.loci[]==0),"repeated loci were found","\n")	
			cat("----","-!- ","-->","cross validation is done separately for each of them","\n")
			cat("----","-!- ","-->","less than 100 loci in each group can give inaccurate results","\n")
			cat("\n");cat("\n")
		}#for(set.nr ...
		}














###	---------------------------------------------------------------------------------------------------
###	4. gamma_cross_validation_procedure__function__May2016
###	---------------------------------------------------------------------------------------------------
gamma_cross_validation_procedure__function__May2016<-function(
		In.directory	,
		In.validation.directory,
		Res.directory,
		data.list,
		validation.list,
		min.coverage.accepted,
		fast.protocol,
		unique.count.protocol.nr,
		additional.name){
	### function		####################################################################	
		
	DATA_LIST_NUMBER <-1				
	for(DATA_LIST_NUMBER in 1:length(data.list)){	
	# 	--------------------------------------------------------------------------------------------------
	#.. constant:
		categories					<-c("All.pos","Rep.pos","NotRep.pos","VarRep.pos","VarNotRep","Stable.pos")
		snp.names					<-c("A","C","G","T","z")
		snp.values					<-c(1,2,4,8,16,32)
		names(snp.values)			<-c(snp.names,"n")	
		validation.procedure.name	<-"Gamma"	
		snp.names					<-snp.names
				
	# 	--------------------------------------------------------------------------------------------------
	#..	set.name
 		temp_Data_set_name		<-data.list[DATA_LIST_NUMBER]
 		temp_Data_set_name		<-strsplit(temp_Data_set_name,split="_ARRAY_CovACGT")
 		temp_Data_set_name		<-temp_Data_set_name[[1]]
 		temp_Data_set_name		<-temp_Data_set_name[1]	 					
	#..	display:		
		cat("[...] . ","\n"); cat("[...]",".......................................................","\n")
		cat("[",DATA_LIST_NUMBER,"] - ",temp_Data_set_name,"\n");	 
		cat("[...] - ",date(),"\n"); 					
		cat("[...] - ","done with:","\n"); 		
	#.. load covARRAY
		setwd(In.directory); 
		load(data.list[DATA_LIST_NUMBER])		#ARRAY_CovACGT	;dim(ARRAY_CovACGT)		
		#extract data (shorter.name:)
		ar<-ARRAY_CovACGT[,,,drop=FALSE]
		rm(ARRAY_CovACGT)	
	#.. load other elements:
		setwd(In.validation.directory); 
		load(validation.list[grepl(temp_Data_set_name,validation.list)])
		# and exctract.data:
		isolate.names.to.use			<-additional_elements_for_CV[[1]]
		isolates.to.test				<-additional_elements_for_CV[[2]]
		not.repeated.loci			<-additional_elements_for_CV[[3]]
		# names(additional_elements_for_CV)
	#.. info
		cat("[...] . ---> ",data.list[DATA_LIST_NUMBER],"\n")	
		cat("[...] . ---> ",validation.list[grepl(temp_Data_set_name,validation.list)],"\n")	
		cat("[...] . ","\n")		
	
	# 	--------------------------------------------------------------------------------------------------		
	# 	Cross.validation list of results for each isolate
	# 	--------------------------------------------------------------------------------------------------	
	#.. list for results: 
	#.. Samples.Cross.Refference				
		Samples.Cross.Refference									<-as.list(names(isolates.to.test))	#list to fill later on
		names(Samples.Cross.Refference)							<-names(isolates.to.test)		
		Samples.Cross.Refference.essencials.for.substitution		<-Samples.Cross.Refference
	
	# 	--------------------------------------------------------------------------------------------------
	#.. THE LOOP
	# 	--------------------------------------------------------------------------------------------------
		isolate.to.subset.nr<-1
		for(isolate.to.subset.nr in 1:length(isolates.to.test)){			
		##	info:	
			cat("\n")
			cat("     -----------------------------------------------------------------------","\n")
			cat("[pr+] ->",isolate.to.subset.nr," - ",
			names(isolates.to.test)[isolate.to.subset.nr]," -> repl.nr = ",
			isolates.to.test[isolate.to.subset.nr],"...at:",date(),"\n");								
		## 	get.out.necessary:
			isolate.to.subset.name		<-names(isolates.to.test)[isolate.to.subset.nr]
			find.isolate.in.ar			<-as.numeric(grepl(isolate.to.subset.name,isolate.names.to.use))
			ar.s							<-ar[,which(find.isolate.in.ar[]==1),,drop=FALSE]; dim(ar.s)		
		##	replicate.info:
			replicate.nr					<-dim(ar.s)[[2]]	; 	replicates.nr<-replicate.nr
			original.replicates.names	<-colnames(ar)[grepl(isolate.to.subset.name,isolate.names.to.use)]
			INFO.list					<-list(isolate.to.subset.name,replicate.nr,original.replicates.names)
			names(INFO.list)			<-c("isolate.to.subset.name","replicates.nr","original.replicates.names")
												
		##	make.binary.array:
			# here I introduced chnages for checking z sites with general coverage, not with filtered coverage !
			# 2015.07.17
			ar.s.bin							<-ar.s[,,1:6,drop=FALSE]; dim(ar.s.bin)
			dimnames(ar.s.bin)[[3]]			<-c(snp.names,"n")
			for(i in 1:4){
				temp							<-ar.s.bin[,,i,drop=FALSE]; 
				temp[which(ar.s.bin[,,i]>0)]<-1; 
				ar.s.bin[,,i]				<-temp}
			temp								<-ar.s[,,6,drop=FALSE]; dim(temp)
			temp[which(ar.s[,,6]==0)]		<-1; 
			temp[which(ar.s[,,6]>0)]			<-0; 
			ar.s.bin[,,5]					<-temp; 
			rm(temp)						
		##  find.n's (ghost.zone.more.less)	
			n.zone.from		<-0
			n.zone.till		<-as.numeric(min.coverage.accepted)-1		#this is why I give 10 coverage not 9 
			mat.s.n.zone		<-ar.s[,,6,drop=FALSE]; dim(mat.s.n.zone)
			mat.s.cov.temp	<-mat.s.n.zone
			mat.s.n.zone[,,]<-0	
			mat.s.n.zone[which(mat.s.cov.temp[,,]>=n.zone.from & mat.s.cov.temp[,,]<=n.zone.till)]<-1
			rm(mat.s.cov.temp)
			sum(mat.s.n.zone)
			ar.s.bin[,,6]	<-mat.s.n.zone		
		##  put.one.by.one.into.list:
			n.zone.list.for.each.sample			<-as.list(original.replicates.names)
			names(n.zone.list.for.each.sample)	<-paste(isolate.to.subset.name,"__",original.replicates.names,"-n.zone",sep="")
			for(is.nr in 1:ncol(mat.s.n.zone)){n.zone.list.for.each.sample[[is.nr]]<-as.numeric(mat.s.n.zone[,is.nr,])}#is.nr<-1
			#rm(mat.s.n.zone)		
		##	ar.s.nr
			ar.s.nr			<-ar.s.bin	
			for(i.dim in 1:dim(ar.s.nr)[3]){
				temp							<-ar.s.nr[,,i.dim,drop=FALSE]	
				temp[which(ar.s.bin[,,i.dim]>0)]<-snp.values[i.dim]
				ar.s.nr[,,i.dim]<-temp
				}			
			#s.t	<-ar.s.nr
			#colnames(s.t)<-1:ncol(s.t)
			#rownames(s.t)<-1:nrow(s.t)
			#cbind(s.t[1:50,1:5,1],s.t[1:50,1:5,2],s.t[1:50,1:5,3],s.t[1:50,1:5,4],s.t[1:50,1:5,5],s.t[1:50,1:5,6])			
		## constant.vs.variable.loci
			nb.of.samples				<-replicate.nr		
			if(nb.of.samples[]>1){								
				ar.s.nr.sum				<-ar.s.nr[,,1,drop=FALSE]
				ar.s.nr.sum[,,]			<-0
				i.dim<-1
				for(i.dim in 1:dim(ar.s.bin)[3]){
					temp							<-ar.s.nr[,,i.dim,drop=FALSE]
					temp[which(temp[,,]>0)]		<-snp.values[i.dim]
					ar.s.nr.sum					<-ar.s.nr.sum+temp}
				r.nr<-1
				constant.loci					<-rep(0,nrow(ar.s.nr.sum))
				ar.s.nr.sum						<-ar.s.nr.sum[,,1,drop=TRUE]; dim(ar.s.nr.sum)
				constant.loci					<-apply(ar.s.nr.sum,1,function(x){as.numeric(length(table(x))[]==1)})
				variable.loci					<-as.numeric(constant.loci[]==0)
			}else{
				constant.loci		<-as.numeric(rep(1,nrow(ar.s.bin)))
				variable.loci		<-as.numeric(rep(0,nrow(ar.s.bin)))
			}#if(nb.of.samples[]>1)		
		##	other.loci categories for cross-validation:
			All.loci						<-as.numeric(rep(1,length(constant.loci)))
			repeated.loci				<-as.numeric(not.repeated.loci[]==0)
			var.repeated.loci			<-as.numeric((repeated.loci+variable.loci)[]==2)
			var.not.repeated.loci		<-as.numeric((not.repeated.loci+variable.loci)[]==2)
		##  Vectors.with.categories.positions
			Vectors.with.categories.positions<-list(
												All.loci,		
												repeated.loci,
												not.repeated.loci,
												var.repeated.loci,
												var.not.repeated.loci,
												constant.loci)
			names(Vectors.with.categories.positions)<-paste(categories,"-vector.radloci",sep="")
		## 	little info
			tot<-c(ncol(ar.s.bin)*nrow(ar.s.bin))
			cat("-----","\n")
			cat("[+++] ->",round(c(ncol(ar.s.bin)*nrow(ar.s.bin)),digits=2)*100, "of loci in total among",nb.of.samples,"samples","\n")				
			cat("[+++] ->",round(sum(constant.loci)*nb.of.samples/tot,digits=2)*100,"% of loci with exactly same allele compositions within all samples have been found","\n")	
			cat("----- ->","these are your best data and validation procedure will give them full confidence (i.e val = 1)","\n")	
			cat("-----","\n")
			cat("[+++] ->",round(sum(variable.loci)*nb.of.samples/tot,digits=2)*100,"% of loci with exactly variable allele compositions have been found","\n")	
			cat("[+++] ->",round(c(sum(ar.s.bin[,,6])-sum(ar.s.bin[,,5]))/tot,digits=2)*100,
			"% of loci with cov >0 and was excluded from allele analysis due to prefiltration procedure","\n")	
			cat("[+++] ->",round(sum(ar.s.bin[,,5])/tot,digits=2)*100,"% of loci has coverage 0","\n")
			cat("-----","validating categories 1 to 6","\n")
		##	Loci.in.each.category
			Loci.nr.in.each.category				<-matrix(ncol=length(Vectors.with.categories.positions),nrow=1)
			colnames(Loci.nr.in.each.category)	<-categories
			rownames(Loci.nr.in.each.category)	<-"loci.nr"	
			for(ii in 1:length(Vectors.with.categories.positions)){Loci.nr.in.each.category[,ii]<-sum(Vectors.with.categories.positions[[ii]])}
		
		## 	Cross validation for each isolate - with changes done at 2015.07.20						
			#.1.list:				
				Samples.Cross.validation.Six.Categories			<-as.list(c(categories,c("Loci.nr.in.each.category","Info")))
				names(Samples.Cross.validation.Six.Categories)	<-c(paste(categories,".Cross.validation",sep=""),"Loci.nr.in.each.category","Info")	
				Samples.Cross.validation.Six.Categories[[length(Samples.Cross.validation.Six.Categories)]]		<-INFO.list
				Samples.Cross.validation.Six.Categories[[length(Samples.Cross.validation.Six.Categories)-1]]	<-Loci.nr.in.each.category			
			#.2.	validate samples:
				#unique.count.protocol.nr					<-1	# 	higher values, done only with each snp data
				#unique.count.protocol.nr					<-2 #	smaler values, done with all known positions 
			cat("----- ")
			
			if(fast.protocol[]=="no"){
			for(validation.category in 1:6){
				cat("   ",validation.category,";")
				positions.to.validate							<-as.numeric(Vectors.with.categories.positions[[validation.category]])			
				if(sum(positions.to.validate)[]>0){
					# TO SEE WHETHER WE HAVE ANYTHING TO VALIDATE !	
					n.ar										<-ar.s.bin[which(positions.to.validate[]==1),,1:6,drop=FALSE]; dim(n.ar)
					Samples.Cross.validation.Six.Categories[[validation.category]]<-gamma_validation_v_two__function_20150723(n.ar,unique.count.protocol.nr,positions.to.validate)	
				}else{Samples.Cross.validation.Six.Categories[[validation.category]]<-"empty"}
				};cat(" done ","\n")			
			#tes	t:
			#Samples.Cross.validation.Six.Categories[[4]][19:20]
			#Samples.Cross.validation.Six.Categories[[5]][19:20]
			#Samples.Cross.validation.Six.Categories[[6]][19:20]
			}else if(fast.protocol[]=="yes"){
			cat("----- fast protocl for 1:3; ")
			for(validation.category in 4:6){
				cat("   ",validation.category,";")
				positions.to.validate							<-as.numeric(Vectors.with.categories.positions[[validation.category]])			
				if(sum(positions.to.validate)[]>0){
					# TO SEE WHETHER WE HAVE ANYTHING TO VALIDATE !	
					n.ar										<-ar.s.bin[which(positions.to.validate[]==1),,1:6,drop=FALSE]; dim(n.ar)
					Samples.Cross.validation.Six.Categories[[validation.category]]<-gamma_validation_v_two__function_20150723(n.ar,unique.count.protocol.nr,positions.to.validate)	
				}else{Samples.Cross.validation.Six.Categories[[validation.category]]<-"empty"}
				};cat(" done ","\n")			
				Samples.Cross.validation.Six.Categories[1:3]<-"empty"
				}else{}
		#.(c).........Samples.Cross.Refference.essencials.for.substitution:.................................................	
		#..						
		#.c1.....isolate.position:						
		isolate.position.in.CovACGT.array		<-as.numeric(grepl(isolate.to.subset.name,isolate.names.to.use))
		#..
		#.c2....take.out.only.trust.levels:	ii<-1
		Trust.level.in.each.category				<-as.list(categories)
		names(Trust.level.in.each.category)		<-paste(categories,"-Trust.level",sep="")
		#..
		NA.error.level.in.each.category			<-as.list(categories)
		names(NA.error.level.in.each.category)	<-paste(categories,"-NA.error.level",sep="")		
		###
		for(ii in 1:6){
			temp		<-Samples.Cross.validation.Six.Categories[[ii]]			
			if((Samples.Cross.validation.Six.Categories[[ii]])[1]=="empty"){
					temp		<-temp
			}else{
					temp		<-temp[[19]]	# so It takes data only when I have them
			}	#names(temp)
			Trust.level.in.each.category[[ii]]<-temp}
		###
		for(ii in 1:6){
			temp		<-Samples.Cross.validation.Six.Categories[[ii]]			
			if((Samples.Cross.validation.Six.Categories[[ii]])[1]=="empty"){
					temp		<-temp
			}else{
					temp		<-temp[[20]]	# so It takes data only when I have them
			}	#names(temp)
			NA.error.level.in.each.category	[[ii]]<-temp}		
		###
			# remove.stable.pos.from.n.zone.list.for.each.sample - sometimes it was making problems
			i<-1
			stable.pos.bin	<-as.numeric(Vectors.with.categories.positions[[6]])
			for(i in 1:length(n.zone.list.for.each.sample)){
				temp		<-as.numeric(n.zone.list.for.each.sample[[i]])
				temp[which(stable.pos.bin[]==1)]<-0
				n.zone.list.for.each.sample[[i]]<-temp}		
		#..
		#.c3....Put this.all.into.one.list:				
		Basic.info.to.put.Cross.Comparison	<-list(
													replicate.nr,
													isolate.position.in.CovACGT.array,
													Vectors.with.categories.positions,
													Trust.level.in.each.category,
													INFO.list,
													n.zone.list.for.each.sample,
													NA.error.level.in.each.category,
													constant.loci
													)
		names(Basic.info.to.put.Cross.Comparison)<-c(
													"replicate.nr",
													"isolate.position.in.CovACGT.array",
													"Vectors.with.categories.positions",
													"Trust.level.in.each.category",
													"INFO.list",
													"n.zone.list.for.each.sample",
													"NA.error.level.in.each.category",
													"constant.loci")
		#.d.5....Finaly...>>>essencials	
		Samples.Cross.Refference.essencials.for.substitution[[isolate.to.subset.nr]]<-Basic.info.to.put.Cross.Comparison			
	}# Samples.Cross.Refference
	# 	--------------------------------------------------------------------------------------------------	
		#	test:
		#	Samples.Cross.Refference.essencials.for.substitution[[1]][c(4,7)]
		#	Samples.Cross.Refference.essencials.for.substitution[[2]][c(4,7)]
		#	Samples.Cross.Refference.essencials.for.substitution[[3]][c(4,7)]
		#	Samples.Cross.Refference.essencials.for.substitution[[4]][c(4,7)]
		
	#.[3].		-		Save:
	#......................................................................................			
	setwd(Res.directory);# getwd(); dir()	
		
		#.a) 	Samples.Cross.Refference
				save.name<-paste(temp_Data_set_name,additional.name,"_",validation.procedure.name,"ValidationProcedure_Samples.Cross.Refference.RData",sep="")
				save.name
				
				save(Samples.Cross.Refference,file=save.name)
		
		#.b) 	Samples.Cross.Refference.essencials.for.substitution
				save.name<-paste(temp_Data_set_name,additional.name,"_",validation.procedure.name,"ValidationProcedure_Samples.Cross.Refference.essencials.for.substitution.RData",sep="")
				save(Samples.Cross.Refference.essencials.for.substitution,file=save.name)
		gc(); dir()	
	#...
	cat("[...] . ","\n"); 
	cat("[...] - Samples.Cross.Refference - SAVED","\n");
	cat("[...] - Samples.Cross.Refference.essencials.for.substitution - SAVED","\n");
	cat("[...] - in:","\n"); 
	cat("[...] -",Res.directory,"\n");
	cat("[...] - ",date(),"\n"); cat("[...] . ","\n") 
	cat("[...]",".......................................................","\n")
	cat("[...] . ","\n")
}#DATA_LIST_NR
}






















###	---------------------------------------------------------------------------------------------------
###	5. New.ACGTZ.array.creator__with.gamma.validation__function__May2016
###	---------------------------------------------------------------------------------------------------
New.ACGTZ.array.creator__with.gamma.validation__function__May2016<-function(
		array.data.list, 
		cross.validation.essentials.for.array.list, 
		parameters.for.pbinom,
		Noise.value.for.zero.cov,
		In.directory	,
		In.Valid.essentials.dir,
		Res.directory,
		additional.name,
		coverage.for.pbinom){
		

cat("********************************************************","\n")
cat("\n"); cat(length(array.data.list),"-> acgtz array set's to create:","\n");cat("   -> started at:",date(),"\n"); cat("\n");
cat("********************************************************","\n")				
DATA_LIST_NUMBER <-1


for(DATA_LIST_NUMBER in 1:length(array.data.list)){			

# - 1 - to.start
 	# - -------------------------------------------------------------------------------------------
 	temp_Data_set_name		<-array.data.list[DATA_LIST_NUMBER]
 	temp_Data_set_name		<-strsplit(temp_Data_set_name,split="_ARRAY_CovACGT")
 	temp_Data_set_name		<-temp_Data_set_name[[1]]
 	temp_Data_set_name		<-temp_Data_set_name[1]; 
 							temp_Data_set_name	 				
	setwd(In.directory)	;
	load(array.data.list[DATA_LIST_NUMBER])		#ARRAY_CovACGT	;dim(ARRAY_CovACGT)
	setwd(In.Valid.essentials.dir);
	load(cross.validation.essentials.for.array.list[grepl(temp_Data_set_name,cross.validation.essentials.for.array.list)])	
							# Samples.Cross.Refference.essencials.for.substitution
							# names(Samples.Cross.Refference.essencials.for.substitution)					
	cat("[...] . ","\n"); cat("[...]",".......................................................","\n")
	cat("[",DATA_LIST_NUMBER,"] - ",temp_Data_set_name,"\n"); cat("[...] -  dim's: ",dim(ARRAY_CovACGT),"\n"); 
	cat("[...] -  colapsed.isolates.to.process:",length(Samples.Cross.Refference.essencials.for.substitution),"groups \n");
	cat("[...] - ",date(),"\n"); cat("[...] . ","\n")		

# - 2 - create.elements:
	# - -------------------------------------------------------------------------------------------
	ARRAY_ValidationACGTZ				<-ARRAY_CovACGT[,,1:5]; 
	ARRAY_ValidationACGTZ[,,]			<-0
	dimnames(ARRAY_ValidationACGTZ)[[3]]<-c("A","C","G","T","z")
	ARRAY_CovACGTZ							<-ARRAY_ValidationACGTZ
	ARRAY_BinaryACGTZ						<-ARRAY_ValidationACGTZ
	ARRAY_GenCovACGTZ						<-ARRAY_ValidationACGTZ
	ARRAY_FiltCovACGTZ						<-ARRAY_ValidationACGTZ
	dim(ARRAY_ValidationACGTZ)
	(dimnames(ARRAY_ValidationACGTZ)[[1]])[1:5]
	(dimnames(ARRAY_ValidationACGTZ)[[2]])	
	# a) - stable pos - trust and na.error are equal = 1 
	ARRAY_ValidationACGTZ[,,c(1:5)]		<-1
	
# - 3 - fill.in.cross.validation.array:	isolate.to.subset.nr<-6
	# - -------------------------------------------------------------------------------------------
	isolate.to.subset.nr<-1
	for(isolate.to.subset.nr in 1:length(Samples.Cross.Refference.essencials.for.substitution)){			
							   names(Samples.Cross.Refference.essencials.for.substitution)	
	# 1. ONE.ISOLATE.Valid.Ess:
		ONE.ISOLATE.Valid.Ess	<-Samples.Cross.Refference.essencials.for.substitution[[isolate.to.subset.nr]]	
		names(ONE.ISOLATE.Valid.Ess)
							Replicate.NUMBER					<-ONE.ISOLATE.Valid.Ess[[1]]
							Replicate.position.in.ColArray	<-ONE.ISOLATE.Valid.Ess[[2]]
							Category.position.in.RowArray	<-ONE.ISOLATE.Valid.Ess[[3]]	#names(Category.position.in.RowArray)[5]
							Trust.level.six.categories		<-ONE.ISOLATE.Valid.Ess[[4]]	#names(Trust.level.six.categories)
							NA.unique.level.six.categories	<-ONE.ISOLATE.Valid.Ess[[7]]
							sample.names						<-(ONE.ISOLATE.Valid.Ess[[5]])[[3]]
							n.zone.list.for.each.sample		<-ONE.ISOLATE.Valid.Ess[[6]]
														
		#reverse the order so it is compatible with the old script
		Category.position.in.RowArray.LIST						<-list(Category.position.in.RowArray[[6]],
												  				Category.position.in.RowArray[[5]],
												  				Category.position.in.RowArray[[4]])											  
		names(Category.position.in.RowArray.LIST)				<-c("Stable.pos","VarNotRep.pos","VarRep.pos")	
		#trust.levels:
		Trust.level.categories.LIST								<-list(Trust.level.six.categories[[6]],
										 						 Trust.level.six.categories[[5]],
										 				 		Trust.level.six.categories[[4]])
		names(Trust.level.categories.LIST)						<-c("Stable.pos","VarNotRep.pos","VarRep.pos")
		#na.error.level:
		na.error.level.categories.LIST							<-list(NA.unique.level.six.categories[[6]],
										 						NA.unique.level.six.categories[[5]],
										 				 		NA.unique.level.six.categories[[4]])
		names(na.error.level.categories.LIST)					<-c("Stable.pos","VarNotRep.pos","VarRep.pos")
		#clean:
		#rm(ONE.ISOLATE.Valid.Ess)
		#rm(Trust.level.six.categories)
		#rm(Category.position.in.RowArray)
		#gc()		
			
	# 1b. CONTROL STEP:
		#check.whether.each.category.is.lokated.only.once:		
		ttt<-Category.position.in.RowArray.LIST
		test.for.categories<-sum(as.numeric(ttt[[1]]>0))+sum(as.numeric(ttt[[2]]>0))+sum(as.numeric(ttt[[3]]>0))
		only.ones.at.each.locus<-(sum(test.for.categories)/nrow(ARRAY_CovACGT)); rm(ttt)
		#display:
		cat("[pr+] - ",names(Samples.Cross.Refference.essencials.for.substitution)[isolate.to.subset.nr]," with repl.nr:",Replicate.NUMBER,"\n"); 
		if(only.ones.at.each.locus==1){cat("[...] -  +--->  all.positions.in.array are covered.once:)",date(),"\n")
		}else{cat("[...] -  ERROR !!!","\n"); cat("[...] -  ERROR !!!","\n"); cat("[...] -  Test failed and val is == ",only.ones.at.each.locus,"\n")
		cat("[...] -  some positions are not matched (val < 0)","\n");cat("[...] -  some positions are not matched more than once (val > 0)","\n")}	
		
	# 2. put gamma validation values:
		# first put 1 at all positions - as it would be all filled with stable positions, 
		# so this category doens't have to be filled		
		# b) - trust level at first
			if(Replicate.NUMBER[]==1){
				#do nothing, trust is always 1	
			}else{	
				CAT.NR<-2
				for(CAT.NR in 2:length(Trust.level.categories.LIST)){
					test<-Trust.level.categories.LIST[[CAT.NR]]; test
					#to test whether the category is empty, than do nothing:)
					if(as.vector(test)[1]=="empty"){	
						# again do nothing because these was no positions to validate in this category
					}else{			
						one.cat.trust.level						<-(Trust.level.categories.LIST[[CAT.NR]])[,1:4]	
						one.cat.positions						<-as.numeric(Category.position.in.RowArray.LIST[[CAT.NR]])									
						one.cat.na.error.level					<-(na.error.level.categories.LIST[[CAT.NR]])[,1:4]		
						###
						SNP.NR<-1
						for(SNP.NR in 1:4){	
							one.snp.one.cat.trust.level			<-as.numeric(one.cat.trust.level[,SNP.NR])
							one.snp.one.cat.na.error.level		<-as.numeric(one.cat.na.error.level[,SNP.NR])
							one.snp.pos.to.substitute			<-which(one.cat.positions[]==1)
							s.pos.in.ar							<-which(Replicate.position.in.ColArray[]==1)		
							###
							r.nr									<-1
							for(r.nr	 in 1:length(s.pos.in.ar)){
								trust.one.snp.one.s				<-as.numeric(one.snp.one.cat.trust.level[r.nr])
								na.error.one.snp.one.s			<-as.numeric(one.snp.one.cat.na.error.level[r.nr])
								## - trust
									one.s.pos					<-s.pos.in.ar[r.nr]
									ARRAY_ValidationACGTZ[one.snp.pos.to.substitute,one.s.pos,SNP.NR]<-trust.one.snp.one.s
								## - n.zone.pos
									one.s.n.zone				<-as.numeric(n.zone.list.for.each.sample[[r.nr]])# same for all snp's
									one.s.n.zone.in.cat			<-which(as.numeric(one.s.n.zone+one.cat.positions)[]==2)				
									#
									ARRAY_ValidationACGTZ[one.s.n.zone.in.cat,one.s.pos,SNP.NR]<-na.error.one.snp.one.s
								}#end.r.nr
							}#end SNP.NR	
						}#end if else 	as.vector(test)[1]=="empty"
					}#end CAT.NR
				}#END if(Replicate.NUMBER[]==1)			
			}#end isolate.to.subset.nr
							#################################################
							#test
							#################################################
							#	aa<-ARRAY_ValidationACGTZ
							#	rownames(aa)<-1:nrow(aa)
							#	colnames(aa)<-1:ncol(aa)
							#	aa[one.snp.pos.to.substitute[1:10],1:12,1]
							#	aa[one.snp.pos.to.substitute,7:12,5]	
							#	aa[1:10,7:12,1]
							#################################################
	# 3. remove zero values .- they can be very dangerous for some future scripts
		ARRAY_ValidationACGTZ[which(ARRAY_ValidationACGTZ[,,]==0)]<-Noise.value.for.zero.cov
		#length(which(ARRAY_ValidationACGTZ[,,]==0))
		#length(which(ARRAY_ValidationACGTZ[,,]==Noise.value.for.zero.cov))	
	# 3. save:
		setwd(Res.directory); dir()	
		ARRAY_ValidationACGTZ_name	<-paste(temp_Data_set_name,"_ARRAY_ValidationACGTZ_RG",".RData",sep="");ARRAY_ValidationACGTZ_name
									save(ARRAY_ValidationACGTZ,file=ARRAY_ValidationACGTZ_name)	
		cat("[...] - ","ARRAY_ValidationACGTZ - saved","\n"); cat("[...] . ","\n")							
	
# - 4 - DATA_LIST_NUMBER.Part.2: ARRAZ_CovACGTZ;ARRAT_BinaryACGTZ;ARRAZ_GenCovACGTZ;ARRAZ_FiltCovACGTZ
	# - -------------------------------------------------------------------------------------------	

	
	#ARRAY_BinaryACGTZ			
	for(DIM.NR in 1:4){temp<-ARRAY_CovACGT[,,DIM.NR];temp[which(ARRAY_CovACGT[,,DIM.NR]>0)]<-1; ARRAY_BinaryACGTZ[,,DIM.NR]<-temp} 
	temp<-ARRAY_CovACGT[,,5];temp[,]<-0; temp[which(ARRAY_CovACGT[,,5]==0)]<-1; ARRAY_BinaryACGTZ[,,5]<-temp; #ARRAY_BinaryACGTZ[,4,]

	#ARRAY_CovACGTZ
	ARRAY_CovACGTZ[,,1:4]	<-ARRAY_CovACGT[,,1:4]	
	ARRAY_CovACGTZ[,,5]		<-ARRAY_BinaryACGTZ[,,5]; #ARRAY_CovACGTZ[,4,]
	
	#ARRAY_GenCovACGTZ	- for pbinom. because I have indels in my system!			
	temp<-ARRAY_CovACGT[,,6]; #dimnames(ARRAY_CovACGT)[[3]]
	for(DIM.NR in 1:5){ARRAY_GenCovACGTZ[,,DIM.NR]<-temp}

	#ARRAY_FiltCovACGTZ	- for weight because I dont want to have indels in my system!
	temp<-ARRAY_CovACGT[,,5]; #dimnames(ARRAY_CovACGT)[[3]]
	for(DIM.NR in 1:5){ARRAY_FiltCovACGTZ[,,DIM.NR]<-temp}	
	
	cat("[set] -  +--->  easy ACGTZ arrays done","\n")
	
	
#DATA_LIST_NUMBER.Part.3: pbinom calculations:	
################################################################################	
cat("[...] - pbinom calculations:","\n")
ARRAY_pbinom_ACGTZ_LIST<-as.list(parameters.for.pbinom)
names(ARRAY_pbinom_ACGTZ_LIST)<-paste("ARRAY_pbinom",parameters.for.pbinom,"ACGTZ",sep="_")
for(P.VALUE in 1:length(parameters.for.pbinom)){
	cat("[pr+] - p.succces:",parameters.for.pbinom[P.VALUE]," > ")	
	#P.VALUE<-1; SNP.NR<-1; SAMPLE.NR<-1; LOCUS.NR<-1	;#is.allele.present<-1
	p.succes					<-parameters.for.pbinom[P.VALUE]; 	
	ARRAY_pbinom_ACGTZ		<-ARRAY_BinaryACGTZ; 
	ARRAY_pbinom_ACGTZ[which(ARRAY_BinaryACGTZ[,,5]==0)]<-Noise.value.for.zero.cov	
	#pbinom.for.each.snp			
	for(SNP.NR in 1:4){
		cat(" ",SNP.NR,"snp; ")
		SNP.ARRAY_CovACGTZ					<-ARRAY_CovACGTZ[,,SNP.NR];	
		SNP.ARRAY_BinaryACGTZ				<-ARRAY_BinaryACGTZ[,,SNP.NR]
			if(coverage.for.pbinom[]=="general"){
			SNP.ARRAY_GenCovACGTZ			<-ARRAY_GenCovACGTZ[,,SNP.NR]		
			}else if(coverage.for.pbinom[]=="filtered"){
			SNP.ARRAY_GenCovACGTZ			<-ARRAY_FiltCovACGTZ[,,SNP.NR]
			}else{}
		for(SAMPLE.NR in 1:ncol(ARRAY_BinaryACGTZ)){
			Sample.SNP.ARRAY_CovACGTZ		<-as.numeric(SNP.ARRAY_CovACGTZ[,SAMPLE.NR])	
			Sample.SNP.ARRAY_BinaryACGTZ	<-as.numeric(SNP.ARRAY_BinaryACGTZ[,SAMPLE.NR])
			Sample.SNP.ARRAY_GenCovACGTZ	<-as.numeric(SNP.ARRAY_GenCovACGTZ[,SAMPLE.NR])		
			for(LOCUS.NR in 1:length(Sample.SNP.ARRAY_BinaryACGTZ)){
				is.allele.present			<-Sample.SNP.ARRAY_BinaryACGTZ[LOCUS.NR]
				if(is.allele.present[]==1){
						p.succes			<-p.succes
						k.success			<-Sample.SNP.ARRAY_CovACGTZ[LOCUS.NR]
						n.total				<-Sample.SNP.ARRAY_GenCovACGTZ[LOCUS.NR]
						Probability			<-pbinom(k.success,n.total,p.succes)
				}else{	Probability			<-Noise.value.for.zero.cov}
				#finally:
				ARRAY_pbinom_ACGTZ[LOCUS.NR,SAMPLE.NR,SNP.NR]<-Probability
			}#locus
		}#sample
	}#snp
	cat(" and 5th snp (i.e zero cov.) is done!")
	cat("\n")
	ARRAY_pbinom_ACGTZ_LIST[[P.VALUE]]<-ARRAY_pbinom_ACGTZ		
}#P.VALUE
#s.nr<-12	
#aa<-cbind(round((ARRAY_pbinom_ACGTZ_LIST[[1]])[,s.nr,],digits=2),round((ARRAY_pbinom_ACGTZ_LIST[[5]])[,s.nr,],digits=2))	
#rownames(aa)<-1:nrow(aa); aa	
cat("[...]\n")
cat("[...] -  +--->  ARRAY_pbinom_ACGTZ_LIST done at:",date(),"\n")

#DATA_LIST_NUMBER.Part.4: save:	
################################################################################
setwd(Res.directory)
ARRAY_pbinom_ACGTZ_LIST_name<-paste(temp_Data_set_name,additional.name,"_ARRAY_pbinom_ACGTZ_LIST_RG",".RData",sep="")	
ARRAY_ValidationACGTZ_name	<-paste(temp_Data_set_name,additional.name,"_ARRAY_ValidationACGTZ_RG",".RData",sep="")	
ARRAY_CovACGTZ_name			<-paste(temp_Data_set_name,additional.name,"_ARRAY_CovACGTZ_RG",".RData",sep="")	
ARRAY_BinaryACGTZ_name		<-paste(temp_Data_set_name,additional.name,"_ARRAY_BinaryACGTZ_RG",".RData",sep="")	
ARRAY_GenCovACGTZ_name		<-paste(temp_Data_set_name,additional.name,"_ARRAY_GenCovACGTZ_RG",".RData",sep="")	
ARRAY_FiltCovACGTZ_name		<-paste(temp_Data_set_name,additional.name,"_ARRAY_FiltCovACGTZ_RG",".RData",sep="")	
							save(ARRAY_pbinom_ACGTZ_LIST,	file=ARRAY_pbinom_ACGTZ_LIST_name)		
							save(ARRAY_CovACGTZ,				file=ARRAY_CovACGTZ_name)			
							save(ARRAY_BinaryACGTZ,			file=ARRAY_BinaryACGTZ_name)			
							save(ARRAY_GenCovACGTZ,			file=ARRAY_GenCovACGTZ_name)			
							save(ARRAY_FiltCovACGTZ,			file=ARRAY_FiltCovACGTZ_name)
															rm(ARRAY_pbinom_ACGTZ_LIST)		
															rm(ARRAY_CovACGTZ)		
															rm(ARRAY_BinaryACGTZ)		
															rm(ARRAY_GenCovACGTZ)		
															rm(ARRAY_FiltCovACGTZ)
															gc()
cat("[...] . ","\n"); 
cat("[.s.] - ARRAY_pbinom_ACGTZ_LIST and other - SAVED","\n"); 
cat("[.s.] - in:","\n"); 
cat("[.s.] -",Res.directory,"\n");
cat("[.s.] - ",date(),"\n"); cat("[.s.] . ","\n") 
cat("[.s.]",".......................................................","\n")
cat("[...] . ","\n")
}#DATA_LIST_NUMBER
}

































###	---------------------------------------------------------------------------------------------------
###	6. COLLAPSE_ACGTZ_ARRAYS__GAMMA__function_May2016		
###	---------------------------------------------------------------------------------------------------
COLLAPSE_ACGTZ_ARRAYS__GAMMA__function_May2016<-function(
	#	collapse.info
		Valid.dir,
		collapse.info,	
	#	arrays and results:
		acgtz.ar.dir,
		Res.directory,
	#	array lists:		
		pbinom.ACGTZ,
		valid.ACGTZ,
		cov.ACGTZ,
		gen.cov.ACGTZ,
	#	parameters:
		additional.name,
		pbinom.to.take,
		min.coverage.accepted,
		when.weight.is.zero.put,
		use.full.cov.spectrum.at.na,
	#	set.name
		run.part.without.beta.correction,
		run.part.with.beta.correction,
		prepare.info.arrays,
	#   add to name - my small custom made info	
		Name_add.method){



	DATA_LIST_NUMBER<-1
# ---------------------------------------------------------------------------------------------------------
	for(DATA_LIST_NUMBER in 1:length(pbinom.ACGTZ)){			
# ---------------------------------------------------------------------------------------------------------
 
# 	1. load data:
# 	-------------------------------------------------------------------------------------------------------
	setwd(acgtz.ar.dir);
 	temp_Data_set_name	<-pbinom.ACGTZ[DATA_LIST_NUMBER]; temp_Data_set_name
 	temp_Data_set_name	<-strsplit(temp_Data_set_name,split=".RData")
 	temp_Data_set_name	<-temp_Data_set_name[[1]]
 	temp_Data_set_name	<-temp_Data_set_name[1]; temp_Data_set_name		 	

	# - a - pbinom:		
		setwd(acgtz.ar.dir);load(pbinom.ACGTZ[DATA_LIST_NUMBER])
		ARRAY_pbinom_ACGTZ<-ARRAY_pbinom_ACGTZ_LIST[[pbinom.to.take]]
		I.took.pbinom<-names(ARRAY_pbinom_ACGTZ_LIST)[pbinom.to.take]
		mode(ARRAY_pbinom_ACGTZ)		<-"numeric"
	# - b - cov, validation, :	
		setwd(acgtz.ar.dir);load(valid.ACGTZ[DATA_LIST_NUMBER]); 	dim(ARRAY_ValidationACGTZ)
		setwd(acgtz.ar.dir);load(cov.ACGTZ[DATA_LIST_NUMBER]);		dim(ARRAY_CovACGTZ)
		setwd(acgtz.ar.dir);load(gen.cov.ACGTZ[DATA_LIST_NUMBER]);	dim(ARRAY_GenCovACGTZ)
		#setwd(acgtz.ar.dir);load(filt.cov[DATA_LIST_NUMBER]);		dim(ARRAY_FiltCovACGTZ)
		###
		mode(ARRAY_ValidationACGTZ)		<-"numeric"
		mode(ARRAY_CovACGTZ) 			<-"numeric"
		mode(ARRAY_GenCovACGTZ)			<-"numeric"
		#mode(ARRAY_FiltCovACGTZ)			<-"numeric"
	# - c - collapse info:
		setwd(Valid.dir);load(collapse.info[DATA_LIST_NUMBER]);	
		names(Samples.Cross.Refference.essencials.for.substitution)
		reff.list					<-Samples.Cross.Refference.essencials.for.substitution; 	
		collapse.group.names		<-as.list(names(reff.list)); names(collapse.group.names)<-names(reff.list)
		collapse.group.pos			<-collapse.group.names
		collapse.group.r.nr			<-collapse.group.names
		for(i in 1:length(reff.list)){collapse.group.r.nr[[i]]<-(reff.list[[i]])[[1]]; collapse.group.pos[[i]]<-(reff.list[[i]])[[2]]}
		names(reff.list[[1]])

#   2. info	
# 	-------------------------------------------------------------------------------------------------------
	cat("      -   ........................................................","\n")
	cat("[ ",DATA_LIST_NUMBER," ] - ",temp_Data_set_name,"\n")
	cat("[ ",DATA_LIST_NUMBER," ] -  with:",I.took.pbinom,"\n")
	cat("[ ",DATA_LIST_NUMBER," ] - ",date(),"\n")
	cat("[ ... ] -  elements are:","\n")
	cat("[ ... ] -  ",valid.ACGTZ[DATA_LIST_NUMBER],"\n")
	cat("[ ... ] -  ",cov.ACGTZ[DATA_LIST_NUMBER],"\n")
	cat("[ ... ] -  ",gen.cov.ACGTZ[DATA_LIST_NUMBER],"\n")
	cat("[ ... ] -  ",collapse.info[DATA_LIST_NUMBER],"\n")	
	cat("[ ... ] \n")	
	cat("[ .?. ] -  ","Info arrays 		-",prepare.info.arrays,"\n")	
	cat("[ .?. ] -  ","NoCorr array		-",run.part.without.beta.correction,"\n")	
	cat("[ .?. ] -  ","BetaCorr array 	-",run.part.with.beta.correction,"\n")	
	cat("[...] -   ........................................................","\n")


#   3. Zero.cov.MATRIX (Binary 1 == 0 coverage) from general coverage matrix ! - it is important
# 	-------------------------------------------------------------------------------------------------------
	ARRAY_GenCovACGTZ_original<-ARRAY_GenCovACGTZ[,,1,drop=FALSE]					
	#ARRAY_ValidationACGTZ[which(ARRAY_ValidationACGTZ[,,1:5]==0)]<-as.numeric(when.weight.is.zero.put)	
					#summary(ARRAY_ValidationACGTZ)		
	Zero.cov.MATRIX	<-ARRAY_GenCovACGTZ[,,5,drop=FALSE]
					Zero.cov.MATRIX[which(ARRAY_GenCovACGTZ[,,5]==0)]<-1
					Zero.cov.MATRIX[which(ARRAY_GenCovACGTZ[,,5]>0)]<-0	
																
#   4. Log(coverage,10) in log_ARRAY_GenCovACGTZ
# 	-------------------------------------------------------------------------------------------------------
	log_ARRAY_GenCovACGTZ								<-ARRAY_GenCovACGTZ; dim(log_ARRAY_GenCovACGTZ)
	log_ARRAY_GenCovACGTZ[,,]							<-0
	if(use.full.cov.spectrum.at.na[]=="yes"){
		n.col<-1
		for(n.col in 1:ncol(ARRAY_GenCovACGTZ)){
			d.temp										<-ARRAY_GenCovACGTZ[,n.col,] 	;d.temp
			d.temp										<-d.temp+2						;d.temp# added in 07.2015
			d.temp										<-log(d.temp,10)					;d.temp				
			log_ARRAY_GenCovACGTZ[,n.col,]				<-d.temp}
	}else if(use.full.cov.spectrum.at.na[]=="no"){
		n.col<-1
		for(n.col in 1:ncol(ARRAY_GenCovACGTZ)){
			d.temp										<-ARRAY_GenCovACGTZ[,n.col,] 	;d.temp
			d.temp[which(ARRAY_GenCovACGTZ[,n.col,]<min.coverage.accepted)]<-0			;d.temp
			d.temp										<-d.temp+2						;d.temp# added in 07.2015
			d.temp										<-log(d.temp,10)				;d.temp				
			log_ARRAY_GenCovACGTZ[,n.col,]				<-d.temp}			
	}else{}					
		#	aa<-log_ARRAY_GenCovACGTZ; colnames(aa)<-1:ncol(aa); rownames(aa)<-1:nrow(aa); summary(aa)
		#	hist(aa[,,],breaks=200,col="yellow")
		#	hist(ARRAY_GenCovACGTZ_original[,,],breaks=200,col="yellow")


#   5. empty elements:
# 	-------------------------------------------------------------------------------------------------------
	ARRAY_pbinom.x.log10.x.Weight_ACGTZ					<-ARRAY_pbinom_ACGTZ#for.size 
	ARRAY_pbinom.x.log10.x.Weight_ACGTZ[,,]				<-0
	dimnames(ARRAY_pbinom.x.log10.x.Weight_ACGTZ)[[3]]  <-c("A","C","G","T","z")
	mode(ARRAY_pbinom.x.log10.x.Weight_ACGTZ)			<-"numeric"
	# ------ 			
	ARRAY_pbinom.x.log10.x.Weight.div.BetaCorr_ACGTZ		<-ARRAY_pbinom.x.log10.x.Weight_ACGTZ		
	ARRAY_log10.x.Weigth_ACGTZ							<-ARRAY_pbinom.x.log10.x.Weight_ACGTZ				
	ARRAY_log10.x.Weigth.div.BetaCorr_ACGTZ				<-ARRAY_pbinom.x.log10.x.Weight_ACGTZ
	MATRIX_Beta.Corr_ACGTZ								<-ARRAY_pbinom.x.log10.x.Weight_ACGTZ[,,1]


#   6. New Beta correction
# 	-------------------------------------------------------------------------------------------------------
	cat("[ ",DATA_LIST_NUMBER," ] - New Beta correction calculation","\n")
	MATRIX_Beta.Corr_ACGTZ[,]	<-1		
	is.nr<-1		
	for(is.nr in 1:length(reff.list)){		
		names(reff.list[[1]])	
		repl.nr		<-reff.list[[is.nr]][[1]]
		is.pos		<-which((reff.list[[is.nr]][[2]])[]==1)
		na.pos		<-reff.list[[is.nr]][[6]]
		if(repl.nr[]==1){
			MATRIX_Beta.Corr_ACGTZ[,is.pos]	<-repl.nr
		}else{
				na.pos.sum	<-rep(0,length(na.pos[[1]]))
							for(r.nr in 1:length(is.pos)){na.pos.sum<-c(na.pos.sum+na.pos[[r.nr]])}	
			# we dont chnage positions with na at all replicates (cpu time)	
				na.pos.sum[which(na.pos.sum[]==repl.nr)]					<-0
				noNA.pos.sum												<-as.numeric(repl.nr)-na.pos.sum
				noNA.pos.sum[which(noNA.pos.sum[]==as.numeric(repl.nr))]	<-1
			# introduce these values accodrdingly		
				r.nr<-1
				for(r.nr in 1:length(is.pos)){
					one.s.pos						<-is.pos[r.nr]
					na.pos.in.repl					<-na.pos[[r.nr]]
					###
					na.pos.in.repl					<-which(na.pos.in.repl[]==1)
					Beta.for.na.pos.in.repl			<-noNA.pos.sum[na.pos.in.repl]
					#names(Beta.for.na.pos.in.repl)	<-na.pos.in.repl	
					###
					MATRIX_Beta.Corr_ACGTZ[na.pos.in.repl,one.s.pos]	<-Beta.for.na.pos.in.repl
					#as.numeric(MATRIX_Beta.Corr_ACGTZ[,one.s.pos])
					}
					#aa<-MATRIX_Beta.Corr_ACGTZ[,is.pos]		
					#rownames(aa)<-1:nrow(aa)	
					#colnames(aa)<-1:ncol(aa)
					#aa	
			}#end else		
		}	#aa<-MATRIX_Beta.Corr_ACGTZ[,]; rownames(aa)<-1:nrow(aa); 	colnames(aa)<-1:ncol(aa); aa[1:100,1:20]		
			#colnames(MATRIX_Beta.Corr_ACGTZ)				
			
			
#   7a. Partial data for collapse array (No beta correction)
# 	-------------------------------------------------------------------------------------------------------
	# - 1 - ARRAY_log10.x.Weigth_ACGTZ
	for(n.dim in 1:4){
	ARRAY_log10.x.Weigth_ACGTZ[,,n.dim]<-(log_ARRAY_GenCovACGTZ[,,n.dim,drop=FALSE]*ARRAY_ValidationACGTZ[,,n.dim,drop=FALSE])}
	#	aa<-ARRAY_log10.x.Weigth_ACGTZ; colnames(aa)<-1:ncol(aa); rownames(aa)<-1:nrow(aa); aa[,,5]
	#	summary(ARRAY_log10.x.Weigth_ACGTZ); min(ARRAY_log10.x.Weigth_ACGTZ)	

	# - 2 - ARRAY_pbinom.x.log10.x.Weight_ACGTZ
	for(n.dim in 1:4){
	ARRAY_pbinom.x.log10.x.Weight_ACGTZ[,,n.dim]<-(ARRAY_log10.x.Weigth_ACGTZ[,,n.dim,drop=FALSE]*ARRAY_pbinom_ACGTZ[,,n.dim,drop=FALSE])}
	#	aa<-ARRAY_pbinom.x.log10.x.Weight_ACGTZ; colnames(aa)<-1:ncol(aa); rownames(aa)<-1:nrow(aa); aa[,,5]
	#	summary(ARRAY_pbinom.x.log10.x.Weight_ACGTZ[,,5]); min(ARRAY_pbinom.x.log10.x.Weight_ACGTZ[,,5])

	# - 3 - correction for "z"
	ARRAY_log10.x.Weigth_ACGTZ[,,5]			<-Zero.cov.MATRIX
	ARRAY_pbinom.x.log10.x.Weight_ACGTZ[,,5]<-Zero.cov.MATRIX
	#	#aa<-ARRAY_pbinom.x.log10.x.Weight_ACGTZ; colnames(aa)<-1:ncol(aa); rownames(aa)<-1:nrow(aa); aa
	#	#aa[1:100,1:5,5]; l<-53; round(aa[l,1:5,],digits=2)
	#	#min(aa)

#   7b. Partial data for collapse array (With beta correction)
# 	-------------------------------------------------------------------------------------------------------
	# - 1 - ARRAY_log10.x.Weigth_ACGTZ
	for(n.dim in 1:4){
	ARRAY_log10.x.Weigth.div.BetaCorr_ACGTZ[,,n.dim]<-ARRAY_log10.x.Weigth_ACGTZ[,,n.dim]/MATRIX_Beta.Corr_ACGTZ}
	#	aa<-ARRAY_log10.x.Weigth.div.BetaCorr_ACGTZ; colnames(aa)<-1:ncol(aa); rownames(aa)<-1:nrow(aa); round(aa,digits=3)
	#	min(MATRIX_Beta.Corr_ACGTZ)
	#	summary(ARRAY_log10.x.Weigth.div.BetaCorr_ACGTZ[,,1:4]); min(ARRAY_log10.x.Weigth.div.BetaCorr_ACGTZ[,,1:4])

	# - 2 - ARRAY_pbinom.x.log10.x.Weight.div.BetaCorr_ACGTZ
	for(n.dim in 1:4){
	ARRAY_pbinom.x.log10.x.Weight.div.BetaCorr_ACGTZ[,,n.dim]<-(ARRAY_log10.x.Weigth.div.BetaCorr_ACGTZ[,,n.dim,drop=FALSE]*ARRAY_pbinom_ACGTZ[,,n.dim,drop=FALSE])}
	
	# - 3 - correction for "z"
	ARRAY_log10.x.Weigth.div.BetaCorr_ACGTZ[,,5]			 <-Zero.cov.MATRIX
	ARRAY_pbinom.x.log10.x.Weight.div.BetaCorr_ACGTZ[,,5]<-Zero.cov.MATRIX
	#	#aa<-ARRAY_pbinom.x.log10.x.Weight_ACGTZ; colnames(aa)<-1:ncol(aa); rownames(aa)<-1:nrow(aa); aa
	#	#aa[1:100,1:5,5]; l<-53; round(aa[l,1:5,],digits=2)
	#	#min(aa)
	# 	dim(Zero.cov.MATRIX)

#   8. collapse.info
# 	-------------------------------------------------------------------------------------------------------
	setwd(Res.directory);
	###
	f.name			<-(strsplit(temp_Data_set_name,split="_ARRAY")[[1]])[1]; f.name
	binom.name		<-(strsplit(I.took.pbinom,split="_")[[1]])[3];binom.name
	collapse.info.acgtz			<-list(collapse.group.names,collapse.group.r.nr,collapse.group.pos)
	names(collapse.info.acgtz)	<-c("collapse.group.names","collapse.group.r.nr","collapse.group.pos")
	save.name					<-paste(f.name,additional.name,"_collapse.info.acgtz_RG.RData",sep=""); save.name								 
								save(collapse.info.acgtz,file=save.name)
								cat("[ ",DATA_LIST_NUMBER," ] ---",save.name," --- saved","\n")
								cat("\n")
								cat("!!!!    the next step can take a while !!!!","\n")
								cat("\n")
								
								
#   9. Collapse arrays - 2 long scripts for beta correctionand without
# 	-------------------------------------------------------------------------------------------------------

	# -------------------------------------------------------------------
	# A)	ARRAY collapse - No correction.
	# -------------------------------------------------------------------
	if(run.part.without.beta.correction[]=="yes"){
	setwd(Res.directory);
	###
	Ln.x.Wg.arr			<-ARRAY_log10.x.Weigth_ACGTZ
	pb.x.Ln.x.Wg.arr		<-ARRAY_pbinom.x.log10.x.Weight_ACGTZ
	z.cov.mat			<-Zero.cov.MATRIX
	group.names			<-collapse.group.names
	group.pos			<-collapse.group.pos	
	min.for.p			<-when.weight.is.zero.put	
	ARRAY_collapsed_pbinom_ACGTZ		<-GAMMA_collapse_acgtz_arrays_function20150725(
									Ln.x.Wg.arr,
									pb.x.Ln.x.Wg.arr,
									z.cov.mat,
									group.names,
									group.pos,
									min.for.p)
																		
	save.name	<-paste(f.name,additional.name,"_p",binom.name,"_NoCorr_ARRAY_collapsed_pbinom_ACGTZ",Name_add.method,"_RG.RData",sep=""); save.name
				save(ARRAY_collapsed_pbinom_ACGTZ,file=save.name)
				cat("[ ",DATA_LIST_NUMBER," ] - I A - no beta correction -",save.name,date(),"\n")	
	# -------------------------------------------------------------------
	if(prepare.info.arrays[]=="yes"){
	pbinom.arr		<-ARRAY_collapsed_pbinom_ACGTZ
	cov.arr			<-ARRAY_CovACGTZ
	gen.cov.arr		<-ARRAY_GenCovACGTZ_original
	z.cov.mat		<-Zero.cov.MATRIX
	group.names		<-collapse.group.names
	group.pos		<-collapse.group.pos	
	use.as.separator.in.info.arr			<-";"										
	ARRAY_collapsed_info_ACGTZ<-New_INFO_for_Collapse_acgtz_arrays__function_20160601(
											pbinom.arr,
											cov.arr	,
											gen.cov.arr,
											z.cov.mat,
											group.names,
											group.pos,
											use.as.separator.in.info.arr)
	###add. original sample names and repl nr						
	for(i in 1:length(reff.list)){
		gr.name			<-names(reff.list)[i];	gr.name	
		gr.pos			<-which(grepl(gr.name,colnames(ARRAY_collapsed_info_ACGTZ))[]==TRUE)						
		gr.repl.nr		<-reff.list[[i]][[1]][1]
		s.n				<-reff.list[[i]][[5]][[3]]
		gr.n				<-strsplit(s.n[1],split="__")[[1]][1]			
		for(ii in 1:length(s.n)){s.n[ii]<-strsplit(s.n[ii],split="__")[[1]][2]}
		org.s.names	<-paste(s.n,collapse=";")
		###
		dimnames(ARRAY_collapsed_info_ACGTZ)[[3]]
		ARRAY_collapsed_info_ACGTZ[,gr.pos,9]		<-gr.n
		ARRAY_collapsed_info_ACGTZ[,gr.pos,10]		<-as.numeric(gr.repl.nr)				
		ARRAY_collapsed_info_ACGTZ[,gr.pos,11]		<-org.s.names				
		}								
	save.name	<-paste(f.name,additional.name,"_p",binom.name,"_NoCorr_ARRAY_collapsed_info_ACGTZ",Name_add.method,"_RG.RData",sep=""); save.name
				save(ARRAY_collapsed_info_ACGTZ,file=save.name)
				cat("[ ",DATA_LIST_NUMBER," ] - I B - no beta correction -",save.name,date(),"\n")														
				#	no.corr<-ARRAY_collapsed_info_ACGTZ[,,]
	}else{
	cat("[ ",DATA_LIST_NUMBER," ] - I B - ARRAY INFO without beta corr. -"," HAS NOT BEEN DONE - TRY CHNAGES IN F. SETTINGS","\n")
	}
	
	}else{
	cat("[ ",DATA_LIST_NUMBER," ] - I A+B - ARRAY and ARRAY INFO without beta corr. -"," HAS NOT BEEN DONE - TRY CHNAGES IN F. SETTINGS","\n")	
	}

	# -------------------------------------------------------------------
	# B)		ARRAY collapse - Beta correction.
	# -------------------------------------------------------------------
	if(run.part.with.beta.correction[]=="yes"){
	setwd(Res.directory);
	###
	Ln.x.Wg.arr						<-ARRAY_log10.x.Weigth.div.BetaCorr_ACGTZ
	pb.x.Ln.x.Wg.arr					<-ARRAY_pbinom.x.log10.x.Weight.div.BetaCorr_ACGTZ
	z.cov.mat						<-Zero.cov.MATRIX
	group.names						<-collapse.group.names
	group.pos						<-collapse.group.pos	
	min.for.p						<-when.weight.is.zero.put								
	ARRAY_collapsed_pbinom_ACGTZ		<-GAMMA_collapse_acgtz_arrays_function20150725(
								Ln.x.Wg.arr,
								pb.x.Ln.x.Wg.arr,
								z.cov.mat,
								group.names,
								group.pos,
								min.for.p)			
	save.name	<-paste(f.name,additional.name,"_p",binom.name,"_BetaCorr_ARRAY_collapsed_pbinom_ACGTZ",Name_add.method,"_RG.RData",sep=""); save.name
				save(ARRAY_collapsed_pbinom_ACGTZ,file=save.name)
				cat("[ ",DATA_LIST_NUMBER," ] - II A - with beta corr -",save.name,date(),"\n")
	# -------------------------------------------------------------------
	if(prepare.info.arrays[]=="yes"){	
	pbinom.arr					<-ARRAY_collapsed_pbinom_ACGTZ
	cov.arr						<-ARRAY_CovACGTZ
	gen.cov.arr					<-ARRAY_GenCovACGTZ_original
	z.cov.mat					<-Zero.cov.MATRIX
	group.names					<-collapse.group.names
	group.pos					<-collapse.group.pos	
	use.as.separator.in.info.arr			<-";"										
	ARRAY_collapsed_info_ACGTZ<-New_INFO_for_Collapse_acgtz_arrays__function_20160601(
											pbinom.arr,
											cov.arr	,
											gen.cov.arr,
											z.cov.mat,
											group.names,
											group.pos,
											use.as.separator.in.info.arr)
	###add. original sample names and repl nr
	for(i in 1:length(reff.list)){
		gr.name			<-names(reff.list)[i];	gr.name	
		gr.pos			<-which(grepl(gr.name,colnames(ARRAY_collapsed_info_ACGTZ))[]==TRUE)						
		gr.repl.nr		<-reff.list[[i]][[1]][1]
		s.n				<-reff.list[[i]][[5]][[3]]
		gr.n				<-strsplit(s.n[1],split="__")[[1]][1]			
		for(ii in 1:length(s.n)){s.n[ii]<-strsplit(s.n[ii],split="__")[[1]][2]}
		org.s.names	<-paste(s.n,collapse=";")
		###
		dimnames(ARRAY_collapsed_info_ACGTZ)[[3]]
		ARRAY_collapsed_info_ACGTZ[,gr.pos,9]		<-gr.n
		ARRAY_collapsed_info_ACGTZ[,gr.pos,10]		<-as.numeric(gr.repl.nr)				
		ARRAY_collapsed_info_ACGTZ[,gr.pos,11]		<-org.s.names				
		}																	
	save.name	<-paste(f.name,additional.name,"_p",binom.name,"_BetaCorr_ARRAY_collapsed_info_ACGTZ",Name_add.method,"_RG.RData",sep=""); save.name
				save(ARRAY_collapsed_info_ACGTZ,file=save.name)
				cat("[ ",DATA_LIST_NUMBER," ] -  II B - with beta corr -",save.name,date(),"\n")
				#	beta.corr<-ARRAY_collapsed_info_ACGTZ[,,]																	
	}else{
	cat("[ ",DATA_LIST_NUMBER," ] - II B - ARRAY INFO with beta corr. -"," HAS NOT BEEN DONE - TRY CHNAGES IN F. SETTINGS","\n")
	}
	
	}else{
	cat("[ ",DATA_LIST_NUMBER," ] - II A+B - ARRAY and ARRAY INFO with beta corr. -"," HAS NOT BEEN DONE - TRY CHNAGES IN F. SETTINGS","\n")	
	}

}#for(DATA_LIST_NUMBER....
}



































###	---------------------------------------------------------------------------------------------------
###	7. Conservative_method__COLLAPSE_ACGTZ_ARRAYS__GAMMA__function__May2016	1500		
###	---------------------------------------------------------------------------------------------------
Conservative_method__COLLAPSE_ACGTZ_ARRAYS__GAMMA__function__May2016<-function(
	#	collapse.info
		Valid.dir,
		collapse.info,	
	#	arrays and results:
		acgtz.ar.dir,
		Res.directory,
	#	array lists:		
		pbinom.ACGTZ,
		valid.ACGTZ,
		cov.ACGTZ,
		gen.cov.ACGTZ,
	#	for threshold-based methods	
		filt.cov,
		binary.ACGTZ	,
		Name_add.method,	
	#	parameters:
		additional.name,
		pbinom.to.take,
		min.coverage.accepted,
		when.weight.is.zero.put,
		use.full.cov.spectrum.at.na,
	#	set.name
		run.part.without.beta.correction,
		run.part.with.beta.correction,
		prepare.info.arrays){		
	#	DATA_LIST_NUMBER<-2


# ---------------------------------------------------------------------------------------------------------
	for(DATA_LIST_NUMBER in 1:length(pbinom.ACGTZ)){			
# ---------------------------------------------------------------------------------------------------------
  
  
  
# 	1. load data:
# 	-------------------------------------------------------------------------------------------------------
	setwd(acgtz.ar.dir);dir()
 	temp_Data_set_name	<-pbinom.ACGTZ[DATA_LIST_NUMBER]; temp_Data_set_name
 	temp_Data_set_name	<-strsplit(temp_Data_set_name,split=".RData")
 	temp_Data_set_name	<-temp_Data_set_name[[1]]
 	temp_Data_set_name	<-temp_Data_set_name[1]; temp_Data_set_name		 	

	# - a - pbinom:		
		setwd(acgtz.ar.dir);load(pbinom.ACGTZ[DATA_LIST_NUMBER])
		ARRAY_pbinom_ACGTZ<-ARRAY_pbinom_ACGTZ_LIST[[pbinom.to.take]]
		I.took.pbinom<-names(ARRAY_pbinom_ACGTZ_LIST)[pbinom.to.take]
		mode(ARRAY_pbinom_ACGTZ)		<-"numeric"
	# - b - cov, validation, :	
		#setwd(acgtz.ar.dir);load(valid.ACGTZ[DATA_LIST_NUMBER]); 	dim(ARRAY_ValidationACGTZ)
		setwd(acgtz.ar.dir);load(cov.ACGTZ[DATA_LIST_NUMBER]);		dim(ARRAY_CovACGTZ)
		setwd(acgtz.ar.dir);load(gen.cov.ACGTZ[DATA_LIST_NUMBER]);	dim(ARRAY_GenCovACGTZ)
		setwd(acgtz.ar.dir);load(filt.cov[DATA_LIST_NUMBER]);		dim(ARRAY_FiltCovACGTZ)
		setwd(acgtz.ar.dir);load(binary.ACGTZ[DATA_LIST_NUMBER]);	dim(ARRAY_BinaryACGTZ)
		###
		#mode(ARRAY_ValidationACGTZ)		<-"numeric"
		mode(ARRAY_CovACGTZ) 			<-"numeric"
		mode(ARRAY_GenCovACGTZ)			<-"numeric"
		mode(ARRAY_FiltCovACGTZ)			<-"numeric"#needed for threshold based methods
		mode(ARRAY_BinaryACGTZ)			<-"numeric"
	# - c - collapse info:
		setwd(Valid.dir);load(collapse.info[DATA_LIST_NUMBER]);	
		names(Samples.Cross.Refference.essencials.for.substitution)
		reff.list							<-Samples.Cross.Refference.essencials.for.substitution; 	
		collapse.group.names				<-as.list(names(reff.list)); names(collapse.group.names)<-names(reff.list)
		collapse.group.pos					<-collapse.group.names
		collapse.group.r.nr					<-collapse.group.names
		collpase.group.const.pos				<-collapse.group.names
		collpase.group.n.zone.per.sample		<-collapse.group.names
		###
		for(i in 1:length(reff.list)){collapse.group.r.nr[[i]]				<-(reff.list[[i]])[[1]]; collapse.group.pos[[i]]<-(reff.list[[i]])[[2]]}
		for(i in 1:length(reff.list)){collpase.group.const.pos[[i]]			<-as.numeric((reff.list[[i]])[[8]])}
		for(i in 1:length(reff.list)){collpase.group.n.zone.per.sample[[i]]	<-(reff.list[[i]])[[6]]}
	 
	 # - d - info	
		cat("      -   ........................................................","\n")
		cat("[ ",DATA_LIST_NUMBER," ] - ",temp_Data_set_name,"\n")
		cat("[ ",DATA_LIST_NUMBER," ] -  with:",I.took.pbinom,"\n")
		cat("[ ",DATA_LIST_NUMBER," ] - ",date(),"\n")
		cat("[ ... ] -  elements are:","\n")	
		cat("[ ... ] -  ",cov.ACGTZ[DATA_LIST_NUMBER],"\n")
		cat("[ ... ] -  ",gen.cov.ACGTZ[DATA_LIST_NUMBER],"\n")
		cat("[ ... ] -  ",collapse.info[DATA_LIST_NUMBER],"\n")
		cat("[ ... ] -  ",filt.cov[DATA_LIST_NUMBER],"\n")	
		cat("[ ... ] -  ",binary.ACGTZ[DATA_LIST_NUMBER],"\n")	
		cat("[ ... ] \n")	
		cat("[ ... ] -  "," !!!!!!   Conservative_Threshold_method   !!!!!!!","\n")		
		cat("[...] -   ........................................................","\n")


#   3. Zero.cov.MATRIX (Binary 1 == 0 coverage) from general coverage matrix ! - it is important
# 	-------------------------------------------------------------------------------------------------------
	ARRAY_GenCovACGTZ_original<-ARRAY_GenCovACGTZ[,,1,drop=FALSE]					
	#ARRAY_ValidationACGTZ[which(ARRAY_ValidationACGTZ[,,1:5]==0)]<-as.numeric(when.weight.is.zero.put)	
					#summary(ARRAY_ValidationACGTZ)		
	Zero.cov.MATRIX	<-ARRAY_GenCovACGTZ[,,5,drop=FALSE]
					Zero.cov.MATRIX[which(ARRAY_GenCovACGTZ[,,5]==0)]<-1
					Zero.cov.MATRIX[which(ARRAY_GenCovACGTZ[,,5]>0)]<-0	
			
#   4. empty elements:
# 	-------------------------------------------------------------------------------------------------------
	ARRAY_pbinom.x.log10.x.Weight_ACGTZ					<-ARRAY_pbinom_ACGTZ#for.size 
	ARRAY_pbinom.x.log10.x.Weight_ACGTZ[,,]				<-0
	dimnames(ARRAY_pbinom.x.log10.x.Weight_ACGTZ)[[3]]  <-c("A","C","G","T","z")
	mode(ARRAY_pbinom.x.log10.x.Weight_ACGTZ)			<-"numeric"
	# ------ 			
	ARRAY_pbinom.x.log10.x.Weight.div.BetaCorr_ACGTZ		<-ARRAY_pbinom.x.log10.x.Weight_ACGTZ		
	ARRAY_log10.x.Weigth_ACGTZ							<-ARRAY_pbinom.x.log10.x.Weight_ACGTZ				
	ARRAY_log10.x.Weigth.div.BetaCorr_ACGTZ				<-ARRAY_pbinom.x.log10.x.Weight_ACGTZ
	MATRIX_Beta.Corr_ACGTZ								<-ARRAY_pbinom.x.log10.x.Weight_ACGTZ[,,1]

#   5. New Beta correction
# 	-------------------------------------------------------------------------------------------------------
	cat("[ ",DATA_LIST_NUMBER," ] - New Beta correction calculation","\n")
	MATRIX_Beta.Corr_ACGTZ[,]	<-1		
	is.nr<-1		
	for(is.nr in 1:length(reff.list)){		
		names(reff.list[[1]])	
		repl.nr		<-reff.list[[is.nr]][[1]]
		is.pos		<-which((reff.list[[is.nr]][[2]])[]==1)
		na.pos		<-reff.list[[is.nr]][[6]]
		if(repl.nr[]==1){
			MATRIX_Beta.Corr_ACGTZ[,is.pos]	<-repl.nr
		}else{
				na.pos.sum	<-rep(0,length(na.pos[[1]]))
							for(r.nr in 1:length(is.pos)){na.pos.sum<-c(na.pos.sum+na.pos[[r.nr]])}	
			# we dont chnage positions with na at all replicates (cpu time)	
				na.pos.sum[which(na.pos.sum[]==repl.nr)]					<-0
				noNA.pos.sum												<-as.numeric(repl.nr)-na.pos.sum
				noNA.pos.sum[which(noNA.pos.sum[]==as.numeric(repl.nr))]	<-1
			# introduce these values accodrdingly		
				r.nr<-1
				for(r.nr in 1:length(is.pos)){
					one.s.pos						<-is.pos[r.nr]
					na.pos.in.repl					<-na.pos[[r.nr]]
					###
					na.pos.in.repl					<-which(na.pos.in.repl[]==1)
					Beta.for.na.pos.in.repl			<-noNA.pos.sum[na.pos.in.repl]
					#names(Beta.for.na.pos.in.repl)	<-na.pos.in.repl	
					###
					MATRIX_Beta.Corr_ACGTZ[na.pos.in.repl,one.s.pos]	<-Beta.for.na.pos.in.repl
					#as.numeric(MATRIX_Beta.Corr_ACGTZ[,one.s.pos])
					}
					#aa<-MATRIX_Beta.Corr_ACGTZ[,is.pos]		
					#rownames(aa)<-1:nrow(aa)	
					#colnames(aa)<-1:ncol(aa)
					#aa	
			}#end else		
		}	#aa<-MATRIX_Beta.Corr_ACGTZ[,]; rownames(aa)<-1:nrow(aa); 	colnames(aa)<-1:ncol(aa); aa[1:100,1:20]		
			#colnames(MATRIX_Beta.Corr_ACGTZ)				

	
#   6. collapse.info
# 	-------------------------------------------------------------------------------------------------------
	setwd(Res.directory);
	###
	f.name			<-(strsplit(temp_Data_set_name,split="_ARRAY")[[1]])[1]; f.name
	binom.name		<-(strsplit(I.took.pbinom,split="_")[[1]])[3];binom.name
	collapse.info.acgtz			<-list(collapse.group.names,collapse.group.r.nr,collapse.group.pos)
	names(collapse.info.acgtz)	<-c("collapse.group.names","collapse.group.r.nr","collapse.group.pos")
					
#   7. Collapse arrays -Conservative method
# 	-------------------------------------------------------------------------------------------------------

	# A)	ARRAY collapse - Conservative method
	# -------------------------------------------------------------------
	setwd(Res.directory); dir()
	####
	bin.arr				<-ARRAY_BinaryACGTZ
	cov.arr				<-ARRAY_CovACGTZ
	beta.corr.arr		<-MATRIX_Beta.Corr_ACGTZ
	filt.cov.arr			<-ARRAY_FiltCovACGTZ
	###
	z.cov.mat			<-Zero.cov.MATRIX
	group.names			<-collapse.group.names
	group.pos			<-collapse.group.pos	
	min.for.p			<-when.weight.is.zero.put	
	###
	cont.pos			<-collpase.group.const.pos	
	nz.pos				<-collpase.group.n.zone.per.sample		
	###			
	ARRAY_collapsed_pbinom_ACGTZ		<-Collapse_array__Conservative_method(
		bin.arr,		
		cov.arr,				
		beta.corr.arr,	
		filt.cov.arr,		
		z.cov.mat,			
		group.names,		
		group.pos,			
		min.for.p,			
		cont.pos,			
		nz.pos)		
																		
	save.name	<-paste(f.name,additional.name,"_p",binom.name,"_NoCorr_ARRAY_collapsed_pbinom_ACGTZ","_",Name_add.method,"_RG.RData",sep=""); save.name
				save(ARRAY_collapsed_pbinom_ACGTZ,file=save.name)
				cat("[ ",DATA_LIST_NUMBER," ] - I A - collapsed array is done-",save.name,date(),"\n")
				cat("[ ",DATA_LIST_NUMBER," ] - I A - ",Name_add.method,"\n")

	# B)	ARRAY info - Conservative method
	# -------------------------------------------------------------------	
	if(prepare.info.arrays[]=="yes"){
	pbinom.arr		<-ARRAY_collapsed_pbinom_ACGTZ
	cov.arr			<-ARRAY_CovACGTZ
	gen.cov.arr		<-ARRAY_GenCovACGTZ_original
	z.cov.mat		<-Zero.cov.MATRIX
	group.names		<-collapse.group.names
	group.pos		<-collapse.group.pos	
	use.as.separator.in.info.arr			<-";"										
	ARRAY_collapsed_info_ACGTZ<-New_INFO_for_Collapse_acgtz_arrays__function_20160601(
											pbinom.arr,
											cov.arr	,
											gen.cov.arr,
											z.cov.mat,
											group.names,
											group.pos,
											use.as.separator.in.info.arr)
	###add. original sample names and repl nr						
	for(i in 1:length(reff.list)){
		gr.name			<-names(reff.list)[i];	gr.name	
		gr.pos			<-which(grepl(gr.name,colnames(ARRAY_collapsed_info_ACGTZ))[]==TRUE)						
		gr.repl.nr		<-reff.list[[i]][[1]][1]
		s.n				<-reff.list[[i]][[5]][[3]]
		gr.n				<-strsplit(s.n[1],split="__")[[1]][1]			
		for(ii in 1:length(s.n)){s.n[ii]<-strsplit(s.n[ii],split="__")[[1]][2]}
		org.s.names	<-paste(s.n,collapse=";")
		###
		dimnames(ARRAY_collapsed_info_ACGTZ)[[3]]
		ARRAY_collapsed_info_ACGTZ[,gr.pos,9]		<-gr.n
		ARRAY_collapsed_info_ACGTZ[,gr.pos,10]		<-as.numeric(gr.repl.nr)				
		ARRAY_collapsed_info_ACGTZ[,gr.pos,11]		<-org.s.names				
		}								
	save.name	<-paste(f.name,additional.name,"_p",binom.name,"_NoCorr_ARRAY_collapsed_info_ACGTZ","_",Name_add.method,"_RG.RData",sep=""); save.name
				save(ARRAY_collapsed_info_ACGTZ,file=save.name)
				cat("[ ",DATA_LIST_NUMBER," ] - I B - INFO ARRAY DONE:  -",save.name,date(),"\n")														
				#	no.corr<-ARRAY_collapsed_info_ACGTZ[,,]
	}else{}
				#ARRAY_collapsed_info_ACGTZ[53,,1:8]
				#ARRAY_collapsed_pbinom_ACGTZ[53,,]
}#for(DATA_LIST_NUMBER....
}






























###	---------------------------------------------------------------------------------------------------
###	8. Site_elimiantion_method__COLLAPSE_ACGTZ_ARRAYS__GAMMA__function__May2016		
###	---------------------------------------------------------------------------------------------------
Site_elimiantion_method__COLLAPSE_ACGTZ_ARRAYS__GAMMA__function__May2016<-function(
	#	collapse.info
		Valid.dir,
		collapse.info,	
	#	arrays and results:
		acgtz.ar.dir,
		Res.directory,
	#	array lists:		
		pbinom.ACGTZ,
		valid.ACGTZ,
		cov.ACGTZ,
		gen.cov.ACGTZ,
	#	for threshold-based methods	
		filt.cov,
		binary.ACGTZ	,
		Name_add.method,	
	#	parameters:
		additional.name,
		pbinom.to.take,
		min.coverage.accepted,
		when.weight.is.zero.put,
		use.full.cov.spectrum.at.na,
	#	set.name
		run.part.without.beta.correction,
		run.part.with.beta.correction,
		prepare.info.arrays,
		na.max.nr.per.locus,
		eliminate.loci.with.more.than.na.perc,
		ARRAY_CovACGTZ.bis){		
	#	DATA_LIST_NUMBER<-2



DATA_LIST_NUMBER<-1
# ---------------------------------------------------------------------------------------------------------
for(DATA_LIST_NUMBER in 1:length(pbinom.ACGTZ)){			
# ---------------------------------------------------------------------------------------------------------
    
# 	1. load data:
# 	-------------------------------------------------------------------------------------------------------
	setwd(acgtz.ar.dir);dir()
 	temp_Data_set_name	<-pbinom.ACGTZ[DATA_LIST_NUMBER]; temp_Data_set_name
 	temp_Data_set_name	<-strsplit(temp_Data_set_name,split=".RData")
 	temp_Data_set_name	<-temp_Data_set_name[[1]]
 	temp_Data_set_name	<-temp_Data_set_name[1]; temp_Data_set_name		 	

	# - a - pbinom:		
		setwd(acgtz.ar.dir);load(pbinom.ACGTZ[DATA_LIST_NUMBER])
		ARRAY_pbinom_ACGTZ<-ARRAY_pbinom_ACGTZ_LIST[[pbinom.to.take]]
		I.took.pbinom<-names(ARRAY_pbinom_ACGTZ_LIST)[pbinom.to.take]
		mode(ARRAY_pbinom_ACGTZ)		<-"numeric"
	# - b - cov, validation, :	
		#setwd(acgtz.ar.dir);load(valid.ACGTZ[DATA_LIST_NUMBER]); 	dim(ARRAY_ValidationACGTZ)
		setwd(acgtz.ar.dir);load(cov.ACGTZ[DATA_LIST_NUMBER]);		dim(ARRAY_CovACGTZ)
		setwd(acgtz.ar.dir);load(gen.cov.ACGTZ[DATA_LIST_NUMBER]);	dim(ARRAY_GenCovACGTZ)
		setwd(acgtz.ar.dir);load(filt.cov[DATA_LIST_NUMBER]);		dim(ARRAY_FiltCovACGTZ)
		setwd(acgtz.ar.dir);load(binary.ACGTZ[DATA_LIST_NUMBER]);	dim(ARRAY_BinaryACGTZ)
		###
				
		ARRAY_CovACGTZ					<-ARRAY_CovACGTZ.bis
				
		#mode(ARRAY_ValidationACGTZ)		<-"numeric"
		mode(ARRAY_CovACGTZ) 			<-"numeric"
		mode(ARRAY_GenCovACGTZ)			<-"numeric"
		mode(ARRAY_FiltCovACGTZ)			<-"numeric"#needed for threshold based methods
		mode(ARRAY_BinaryACGTZ)			<-"numeric"
	# - c - collapse info:
		setwd(Valid.dir);load(collapse.info[DATA_LIST_NUMBER]);	
		names(Samples.Cross.Refference.essencials.for.substitution)
		reff.list							<-Samples.Cross.Refference.essencials.for.substitution; 	
		collapse.group.names				<-as.list(names(reff.list)); names(collapse.group.names)<-names(reff.list)
		collapse.group.pos					<-collapse.group.names
		collapse.group.r.nr					<-collapse.group.names
		collpase.group.const.pos				<-collapse.group.names
		collpase.group.n.zone.per.sample		<-collapse.group.names
		###
		for(i in 1:length(reff.list)){collapse.group.r.nr[[i]]				<-(reff.list[[i]])[[1]]; collapse.group.pos[[i]]<-(reff.list[[i]])[[2]]}
		for(i in 1:length(reff.list)){collpase.group.const.pos[[i]]			<-as.numeric((reff.list[[i]])[[8]])}
		for(i in 1:length(reff.list)){collpase.group.n.zone.per.sample[[i]]	<-(reff.list[[i]])[[6]]}
	 
	 # - d - info	
		cat("      -   ........................................................","\n")
		cat("[ ",DATA_LIST_NUMBER," ] - ",temp_Data_set_name,"\n")
		cat("[ ",DATA_LIST_NUMBER," ] -  with:",I.took.pbinom,"\n")
		cat("[ ",DATA_LIST_NUMBER," ] - ",date(),"\n")
		cat("[ ... ] -  elements are:","\n")	
		cat("[ ... ] -  ",cov.ACGTZ[DATA_LIST_NUMBER],"\n")
		cat("[ ... ] -  ",gen.cov.ACGTZ[DATA_LIST_NUMBER],"\n")
		cat("[ ... ] -  ",collapse.info[DATA_LIST_NUMBER],"\n")
		cat("[ ... ] -  ",filt.cov[DATA_LIST_NUMBER],"\n")	
		cat("[ ... ] -  ",binary.ACGTZ[DATA_LIST_NUMBER],"\n")	
		cat("[ ... ] \n")	
		cat("[ ... ] -  "," !!!!!!   Site elimination method  !!!!!!!","\n")		
		cat("[...] -   ........................................................","\n")


#   3. Zero.cov.MATRIX (Binary 1 == 0 coverage) from general coverage matrix ! - it is important
# 	-------------------------------------------------------------------------------------------------------
	ARRAY_GenCovACGTZ_original<-ARRAY_GenCovACGTZ[,,1,drop=FALSE]					
	#ARRAY_ValidationACGTZ[which(ARRAY_ValidationACGTZ[,,1:5]==0)]<-as.numeric(when.weight.is.zero.put)	
					#summary(ARRAY_ValidationACGTZ)		
	Zero.cov.MATRIX	<-ARRAY_GenCovACGTZ[,,5,drop=FALSE]
					Zero.cov.MATRIX[which(ARRAY_GenCovACGTZ[,,5]==0)]<-1
					Zero.cov.MATRIX[which(ARRAY_GenCovACGTZ[,,5]>0)]<-0	
			
#   4. empty elements:
# 	-------------------------------------------------------------------------------------------------------
	ARRAY_pbinom.x.log10.x.Weight_ACGTZ					<-ARRAY_pbinom_ACGTZ#for.size 
	ARRAY_pbinom.x.log10.x.Weight_ACGTZ[,,]				<-0
	dimnames(ARRAY_pbinom.x.log10.x.Weight_ACGTZ)[[3]]  <-c("A","C","G","T","z")
	mode(ARRAY_pbinom.x.log10.x.Weight_ACGTZ)			<-"numeric"
	# ------ 			
	ARRAY_pbinom.x.log10.x.Weight.div.BetaCorr_ACGTZ		<-ARRAY_pbinom.x.log10.x.Weight_ACGTZ		
	ARRAY_log10.x.Weigth_ACGTZ							<-ARRAY_pbinom.x.log10.x.Weight_ACGTZ				
	ARRAY_log10.x.Weigth.div.BetaCorr_ACGTZ				<-ARRAY_pbinom.x.log10.x.Weight_ACGTZ
	MATRIX_Beta.Corr_ACGTZ								<-ARRAY_pbinom.x.log10.x.Weight_ACGTZ[,,1]
	
#   6. collapse.info
# 	-------------------------------------------------------------------------------------------------------
	setwd(Res.directory);
	###
	f.name			<-(strsplit(temp_Data_set_name,split="_ARRAY")[[1]])[1]; f.name
	binom.name		<-(strsplit(I.took.pbinom,split="_")[[1]])[3];binom.name
	collapse.info.acgtz			<-list(collapse.group.names,collapse.group.r.nr,collapse.group.pos)
	names(collapse.info.acgtz)	<-c("collapse.group.names","collapse.group.r.nr","collapse.group.pos")
					
#   7. Collapse arrays -Conservative method
# 	-------------------------------------------------------------------------------------------------------

	# A)	ARRAY collapse - site elimination method
	# -------------------------------------------------------------------
	setwd(Res.directory); 
	####
	bin.arr				<-ARRAY_BinaryACGTZ
	cov.arr				<-ARRAY_CovACGTZ
	beta.corr.arr		<-MATRIX_Beta.Corr_ACGTZ
	filt.cov.arr			<-ARRAY_FiltCovACGTZ
	###
	z.cov.mat			<-Zero.cov.MATRIX
	group.names			<-collapse.group.names
	group.pos			<-collapse.group.pos	
	min.for.p			<-when.weight.is.zero.put	
	###
	cont.pos				<-collpase.group.const.pos	
	nz.pos				<-collpase.group.n.zone.per.sample		
	###
	max.na.perc			<-eliminate.loci.with.more.than.na.perc
	max.na.pos.nr		<-na.max.nr.per.locus
	
	###				
	ARRAY_collapsed_pbinom_ACGTZ<-Collapse_array__Site_elimination_method(
		bin.arr,		
		cov.arr,				
		beta.corr.arr,	
		filt.cov.arr,		
		z.cov.mat,			
		group.names,		
		group.pos,			
		min.for.p,			
		cont.pos,			
		nz.pos,
		max.na.perc,
		max.na.pos.nr)
															
	save.name	<-paste(f.name,additional.name,"_p",binom.name,"_NoCorr_ARRAY_collapsed_pbinom_ACGTZ","_",Name_add.method,"_RG.RData",sep=""); save.name
				save(ARRAY_collapsed_pbinom_ACGTZ,file=save.name)
				cat("[ ",DATA_LIST_NUMBER," ] - I A - collapsed array is done-",save.name,date(),"\n")
				cat("[ ",DATA_LIST_NUMBER," ] - I A - ",Name_add.method,"\n")

	# B)	ARRAY info - Conservative method
	# -------------------------------------------------------------------	
	if(prepare.info.arrays[]=="yes"){
	pbinom.arr		<-ARRAY_collapsed_pbinom_ACGTZ
	cov.arr			<-ARRAY_CovACGTZ
	gen.cov.arr		<-ARRAY_GenCovACGTZ_original
	z.cov.mat		<-Zero.cov.MATRIX
	group.names		<-collapse.group.names
	group.pos		<-collapse.group.pos		
	use.as.separator.in.info.arr			<-";"									
	ARRAY_collapsed_info_ACGTZ<-New_INFO_for_Collapse_acgtz_arrays__function_20160601(
											pbinom.arr,
											cov.arr	,
											gen.cov.arr,
											z.cov.mat,
											group.names,
											group.pos,
											use.as.separator.in.info.arr)
	###add. original sample names and repl nr						
	for(i in 1:length(reff.list)){
		gr.name			<-names(reff.list)[i];	gr.name	
		gr.pos			<-which(grepl(gr.name,colnames(ARRAY_collapsed_info_ACGTZ))[]==TRUE)						
		gr.repl.nr		<-reff.list[[i]][[1]][1]
		s.n				<-reff.list[[i]][[5]][[3]]
		gr.n				<-strsplit(s.n[1],split="__")[[1]][1]			
		for(ii in 1:length(s.n)){s.n[ii]<-strsplit(s.n[ii],split="__")[[1]][2]}
		org.s.names	<-paste(s.n,collapse=";")
		###
		dimnames(ARRAY_collapsed_info_ACGTZ)[[3]]
		ARRAY_collapsed_info_ACGTZ[,gr.pos,9]		<-gr.n
		ARRAY_collapsed_info_ACGTZ[,gr.pos,10]		<-as.numeric(gr.repl.nr)				
		ARRAY_collapsed_info_ACGTZ[,gr.pos,11]		<-org.s.names				
		}								
	save.name	<-paste(f.name,additional.name,"_p",binom.name,"_NoCorr_ARRAY_collapsed_info_ACGTZ","_",Name_add.method,"_RG.RData",sep=""); save.name
				save(ARRAY_collapsed_info_ACGTZ,file=save.name)
				cat("[ ",DATA_LIST_NUMBER," ] - I B - INFO ARRAY DONE:  -",save.name,date(),"\n")														
				#	no.corr<-ARRAY_collapsed_info_ACGTZ[,,]
	}else{}
				#ARRAY_collapsed_info_ACGTZ[53,,1:8]
				#ARRAY_collapsed_pbinom_ACGTZ[53,,]
}#for(DATA_LIST_NUMBER....
}


























###	---------------------------------------------------------------------------------------------------
###	9. Fixed_Thresholds_method__COLLAPSE_ACGTZ_ARRAYS__GAMMA__function__May2016		2000
###	---------------------------------------------------------------------------------------------------
Fixed_Thresholds_method__COLLAPSE_ACGTZ_ARRAYS__GAMMA__function__May2016<-function(
	#	collapse.info
		Valid.dir,
		collapse.info,	
	#	arrays and results:
		acgtz.ar.dir,
		Res.directory,
	#	array lists:		
		pbinom.ACGTZ,
		valid.ACGTZ,
		cov.ACGTZ,
		gen.cov.ACGTZ,
	#	for threshold-based methods	
		filt.cov,
		binary.ACGTZ	,
		Name_add.method,	
	#	parameters:
		additional.name,
		pbinom.to.take,
		min.coverage.accepted,
		when.weight.is.zero.put,
		use.full.cov.spectrum.at.na,
	#	set.name
		prepare.info.arrays,
		eliminate.loci.with.more.than.na.perc,
		eliminate.allels.with.less.than.reads.perc,
		minimum.nr.of.replicates.with.same.allele,
		ARRAY_CovACGTZ.bis
		){		
			

			
# DATA_LIST_NUMBER<-1
# ---------------------------------------------------------------------------------------------------------
for(DATA_LIST_NUMBER in 1:length(pbinom.ACGTZ)){			
# ---------------------------------------------------------------------------------------------------------
    
# 	1. load data:
# 	-------------------------------------------------------------------------------------------------------
	setwd(acgtz.ar.dir);dir()
 	temp_Data_set_name	<-pbinom.ACGTZ[DATA_LIST_NUMBER]; temp_Data_set_name
 	temp_Data_set_name	<-strsplit(temp_Data_set_name,split=".RData")
 	temp_Data_set_name	<-temp_Data_set_name[[1]]
 	temp_Data_set_name	<-temp_Data_set_name[1]; temp_Data_set_name		 	

	# - a - pbinom:		
		setwd(acgtz.ar.dir);load(pbinom.ACGTZ[DATA_LIST_NUMBER])
		ARRAY_pbinom_ACGTZ<-ARRAY_pbinom_ACGTZ_LIST[[pbinom.to.take]]
		I.took.pbinom<-names(ARRAY_pbinom_ACGTZ_LIST)[pbinom.to.take]
		mode(ARRAY_pbinom_ACGTZ)		<-"numeric"
	# - b - cov, validation, :	
		#setwd(acgtz.ar.dir);load(valid.ACGTZ[DATA_LIST_NUMBER]); 	dim(ARRAY_ValidationACGTZ)
		setwd(acgtz.ar.dir);load(cov.ACGTZ[DATA_LIST_NUMBER]);		dim(ARRAY_CovACGTZ)
		setwd(acgtz.ar.dir);load(gen.cov.ACGTZ[DATA_LIST_NUMBER]);	dim(ARRAY_GenCovACGTZ)
		setwd(acgtz.ar.dir);load(filt.cov[DATA_LIST_NUMBER]);		dim(ARRAY_FiltCovACGTZ)
		setwd(acgtz.ar.dir);load(binary.ACGTZ[DATA_LIST_NUMBER]);	dim(ARRAY_BinaryACGTZ)
		###
		
		ARRAY_CovACGTZ					<-ARRAY_CovACGTZ.bis
		
		#mode(ARRAY_ValidationACGTZ)		<-"numeric"
		mode(ARRAY_CovACGTZ) 			<-"numeric"
		mode(ARRAY_GenCovACGTZ)			<-"numeric"
		mode(ARRAY_FiltCovACGTZ)			<-"numeric"#needed for threshold based methods
		mode(ARRAY_BinaryACGTZ)			<-"numeric"
	# - c - collapse info:
		setwd(Valid.dir);load(collapse.info[DATA_LIST_NUMBER]);	
		names(Samples.Cross.Refference.essencials.for.substitution)
		reff.list							<-Samples.Cross.Refference.essencials.for.substitution; 	
		collapse.group.names				<-as.list(names(reff.list)); names(collapse.group.names)<-names(reff.list)
		collapse.group.pos					<-collapse.group.names
		collapse.group.r.nr					<-collapse.group.names
		collpase.group.const.pos				<-collapse.group.names
		collpase.group.n.zone.per.sample		<-collapse.group.names
		###
		for(i in 1:length(reff.list)){collapse.group.r.nr[[i]]				<-(reff.list[[i]])[[1]]; collapse.group.pos[[i]]<-(reff.list[[i]])[[2]]}
		for(i in 1:length(reff.list)){collpase.group.const.pos[[i]]			<-as.numeric((reff.list[[i]])[[8]])}
		for(i in 1:length(reff.list)){collpase.group.n.zone.per.sample[[i]]	<-(reff.list[[i]])[[6]]}
	 
	 # - d - info	
		cat("      -   ........................................................","\n")
		cat("[ ",DATA_LIST_NUMBER," ] - ",temp_Data_set_name,"\n")
		cat("[ ",DATA_LIST_NUMBER," ] -  with:",I.took.pbinom,"\n")
		cat("[ ",DATA_LIST_NUMBER," ] - ",date(),"\n")
		cat("[ ... ] -  elements are:","\n")	
		cat("[ ... ] -  ",cov.ACGTZ[DATA_LIST_NUMBER],"\n")
		cat("[ ... ] -  ",gen.cov.ACGTZ[DATA_LIST_NUMBER],"\n")
		cat("[ ... ] -  ",collapse.info[DATA_LIST_NUMBER],"\n")
		cat("[ ... ] -  ",filt.cov[DATA_LIST_NUMBER],"\n")	
		cat("[ ... ] -  ",binary.ACGTZ[DATA_LIST_NUMBER],"\n")	
		cat("[ ... ] \n")	
		cat("[ ... ] -  "," !!!!!!   Fixed Thresholds method  !!!!!!!","\n")		
		cat("[...] -   ........................................................","\n")


#   3. Zero.cov.MATRIX (Binary 1 == 0 coverage) from general coverage matrix ! - it is important
# 	-------------------------------------------------------------------------------------------------------
	ARRAY_GenCovACGTZ_original<-ARRAY_GenCovACGTZ[,,1,drop=FALSE]					
	#ARRAY_ValidationACGTZ[which(ARRAY_ValidationACGTZ[,,1:5]==0)]<-as.numeric(when.weight.is.zero.put)	
					#summary(ARRAY_ValidationACGTZ)		
	Zero.cov.MATRIX	<-ARRAY_GenCovACGTZ[,,5,drop=FALSE]
					Zero.cov.MATRIX[which(ARRAY_GenCovACGTZ[,,5]==0)]<-1
					Zero.cov.MATRIX[which(ARRAY_GenCovACGTZ[,,5]>0)]<-0	
			
#   4. empty elements:
# 	-------------------------------------------------------------------------------------------------------
	ARRAY_pbinom.x.log10.x.Weight_ACGTZ					<-ARRAY_pbinom_ACGTZ#for.size 
	ARRAY_pbinom.x.log10.x.Weight_ACGTZ[,,]				<-0
	dimnames(ARRAY_pbinom.x.log10.x.Weight_ACGTZ)[[3]]  <-c("A","C","G","T","z")
	mode(ARRAY_pbinom.x.log10.x.Weight_ACGTZ)			<-"numeric"
	# ------ 			
	ARRAY_pbinom.x.log10.x.Weight.div.BetaCorr_ACGTZ		<-ARRAY_pbinom.x.log10.x.Weight_ACGTZ		
	ARRAY_log10.x.Weigth_ACGTZ							<-ARRAY_pbinom.x.log10.x.Weight_ACGTZ				
	ARRAY_log10.x.Weigth.div.BetaCorr_ACGTZ				<-ARRAY_pbinom.x.log10.x.Weight_ACGTZ
	MATRIX_Beta.Corr_ACGTZ								<-ARRAY_pbinom.x.log10.x.Weight_ACGTZ[,,1]
	
#   6. collapse.info
# 	-------------------------------------------------------------------------------------------------------
	setwd(Res.directory);
	###
	f.name			<-(strsplit(temp_Data_set_name,split="_ARRAY")[[1]])[1]; f.name
	binom.name		<-(strsplit(I.took.pbinom,split="_")[[1]])[3];binom.name
	collapse.info.acgtz			<-list(collapse.group.names,collapse.group.r.nr,collapse.group.pos)
	names(collapse.info.acgtz)	<-c("collapse.group.names","collapse.group.r.nr","collapse.group.pos")
					
#   7. Collapse arrays -Conservative method
# 	-------------------------------------------------------------------------------------------------------

	# A)	ARRAY collapse - site elimination method
	# -------------------------------------------------------------------
	setwd(Res.directory); 
	####
	bin.arr				<-ARRAY_BinaryACGTZ
	cov.arr				<-ARRAY_CovACGTZ
	beta.corr.arr		<-MATRIX_Beta.Corr_ACGTZ
	filt.cov.arr			<-ARRAY_FiltCovACGTZ
	###
	z.cov.mat			<-Zero.cov.MATRIX
	group.names			<-collapse.group.names
	group.pos			<-collapse.group.pos	
	min.for.p			<-when.weight.is.zero.put	
	###
	cont.pos				<-collpase.group.const.pos	
	nz.pos				<-collpase.group.n.zone.per.sample		
	###
	max.na.perc				<-eliminate.loci.with.more.than.na.perc
	min.allele.adm			<-eliminate.allels.with.less.than.reads.perc
	min.repl.nr.with.allele	<-minimum.nr.of.replicates.with.same.allele		
	###				
	ARRAY_collapsed_pbinom_ACGTZ<-Collapse_array__Fixed_Thresholds_method(
		bin.arr,		
		cov.arr,				
		beta.corr.arr,	
		filt.cov.arr,				
		###
		z.cov.mat,			
		group.names,		
		group.pos,			
		min.for.p,			
		cont.pos,			
		nz.pos,
		max.na.perc,
		min.allele.adm,
		min.repl.nr.with.allele)
	###														
	save.name	<-paste(f.name,additional.name,"_p",binom.name,"_NoCorr_ARRAY_collapsed_pbinom_ACGTZ","_",Name_add.method,"_RG.RData",sep=""); save.name
				save(ARRAY_collapsed_pbinom_ACGTZ,file=save.name)
				cat("[ ",DATA_LIST_NUMBER," ] - I A - collapsed array is done-",save.name,date(),"\n")
				cat("[ ",DATA_LIST_NUMBER," ] - I A - ",Name_add.method,"\n")

	# B)	ARRAY info - Conservative method
	# -------------------------------------------------------------------	
	if(prepare.info.arrays[]=="yes"){
	pbinom.arr		<-ARRAY_collapsed_pbinom_ACGTZ
	cov.arr			<-ARRAY_CovACGTZ
	gen.cov.arr		<-ARRAY_GenCovACGTZ_original
	z.cov.mat		<-Zero.cov.MATRIX
	group.names		<-collapse.group.names
	group.pos		<-collapse.group.pos		
	use.as.separator.in.info.arr			<-";"									
	ARRAY_collapsed_info_ACGTZ<-New_INFO_for_Collapse_acgtz_arrays__function_20160601(
											pbinom.arr,
											cov.arr	,
											gen.cov.arr,
											z.cov.mat,
											group.names,
											group.pos,
											use.as.separator.in.info.arr)
	###add. original sample names and repl nr						
	for(i in 1:length(reff.list)){
		gr.name			<-names(reff.list)[i];	gr.name	
		gr.pos			<-which(grepl(gr.name,colnames(ARRAY_collapsed_info_ACGTZ))[]==TRUE)						
		gr.repl.nr		<-reff.list[[i]][[1]][1]
		s.n				<-reff.list[[i]][[5]][[3]]
		gr.n				<-strsplit(s.n[1],split="__")[[1]][1]			
		for(ii in 1:length(s.n)){s.n[ii]<-strsplit(s.n[ii],split="__")[[1]][2]}
		org.s.names	<-paste(s.n,collapse=";")
		###
		dimnames(ARRAY_collapsed_info_ACGTZ)[[3]]
		ARRAY_collapsed_info_ACGTZ[,gr.pos,9]		<-gr.n
		ARRAY_collapsed_info_ACGTZ[,gr.pos,10]		<-as.numeric(gr.repl.nr)				
		ARRAY_collapsed_info_ACGTZ[,gr.pos,11]		<-org.s.names				
		}								
	save.name	<-paste(f.name,additional.name,"_p",binom.name,"_NoCorr_ARRAY_collapsed_info_ACGTZ","_",Name_add.method,"_RG.RData",sep=""); save.name
				save(ARRAY_collapsed_info_ACGTZ,file=save.name)
				cat("[ ",DATA_LIST_NUMBER," ] - I B - INFO ARRAY DONE:  -",save.name,date(),"\n")														
				#	no.corr<-ARRAY_collapsed_info_ACGTZ[,,]
	}else{}
				#ARRAY_collapsed_info_ACGTZ[53,,1:8]
				#ARRAY_collapsed_pbinom_ACGTZ[53,,]
}#for(DATA_LIST_NUMBER....
}


































###	---------------------------------------------------------------------------------------------------
###	10. Improved_GAMMA_Create_ARRAY_collapsed_short_info_ACGTZ_function20150725	
###	---------------------------------------------------------------------------------------------------
Improved_GAMMA_Create_ARRAY_collapsed_short_info_ACGTZ_function20150725<-function(
dir.list,
In.dir,
Res.dir,
split.by){
	####
for(i in 1:length(dir.list)){
	setwd(In.dir)
	load(dir.list[i])	#i<-66
	###
	f.name		<-(strsplit(dir.list[i],split=split.by)[[1]])[1]; f.name
	end.name		<-(strsplit(dir.list[i],split=split.by)[[1]])[2]; end.name
	save.name	<-paste(f.name,"_ARRAY_collapsed_short_info_ACGTZ",end.name,sep=""); save.name
	###
	cat("...............................................","\n")
	cat("from:",i,"-",dir.list[i],"\n")
	cat("to -> ",i,"-",save.name,"\n")
	cat("   ->",i,"-",date(),"\n"); 
	cat("in :",i,"-",Res.dir,"\n");cat("\n")	
	#:....
	dimnames(ARRAY_collapsed_info_ACGTZ)[[3]]
	#:....
	ARRAY_collapsed_short_info_ACGTZ				<-ARRAY_collapsed_info_ACGTZ[,,1:7]	#this must be modified!!!!!!!!!!!!!
	ARRAY_collapsed_short_info_ACGTZ[,,]			<-0
	dimnames(ARRAY_collapsed_short_info_ACGTZ)[[3]]	<-c("Probability","SNP.presence","Reads.per.snp","General.coverage", "group.name","repl.nr","sample.names")
	ARRAY_collapsed_short_info_ACGTZ[,,4]			<-ARRAY_collapsed_info_ACGTZ[,,8,drop=FALSE]
	ARRAY_collapsed_short_info_ACGTZ[,,3]			<-ARRAY_collapsed_info_ACGTZ[,,7,drop=FALSE]
	ARRAY_collapsed_short_info_ACGTZ[,,2]			<-ARRAY_collapsed_info_ACGTZ[,,6,drop=FALSE]
	####
	ARRAY_collapsed_short_info_ACGTZ[,,5]			<-ARRAY_collapsed_info_ACGTZ[,,9,drop=FALSE]
	ARRAY_collapsed_short_info_ACGTZ[,,6]			<-ARRAY_collapsed_info_ACGTZ[,,10,drop=FALSE]
	ARRAY_collapsed_short_info_ACGTZ[,,7]			<-ARRAY_collapsed_info_ACGTZ[,,11,drop=FALSE]	
	#:....	
	snp.names<-c("A","C","G","T","z")	#	locus.nr<-1; s.nr<-1
	for(locus.nr in 1:nrow(ARRAY_collapsed_info_ACGTZ)){					
		one.locus	<-ARRAY_collapsed_info_ACGTZ[locus.nr,,1:5,drop=FALSE]				
		for(s.nr in 1:ncol(one.locus)){			
			s.o.l	<-as.vector(one.locus[,s.nr,]); s.o.l
			s.o.l	<-paste((paste(snp.names,s.o.l,sep=":")),collapse=";"); s.o.l		
			ARRAY_collapsed_short_info_ACGTZ[locus.nr,s.nr,1]<-s.o.l}	
	}#end for...
	#:....
	setwd(Res.dir);getwd()
	save(ARRAY_collapsed_short_info_ACGTZ, file=save.name)
}#for(i in
}



























###	---------------------------------------------------------------------------------------------------
###	11. gamma_validation_v_two__function_20150723		
###	---------------------------------------------------------------------------------------------------
gamma_validation_v_two__function_20150723<-function(n.ar,unique.count.protocol,positions.to.validate){

###	part.1 - preparations:	###  for all samples it the same:)
		snp.names					<-c("A","C","G","T","z")
	# .	cosmetics:		
		dimnames(n.ar)[[3]]				<-c(snp.names,"n")		
		temp.pos							<-positions.to.validate			
	# .	empty trust.perc.table
		sample.names.original.names		<-colnames(n.ar)
		n.sample							<-ncol(n.ar)
		trust.perc.table				<-matrix(rep(0,n.sample*6),ncol=6)
		rownames(trust.perc.table)		<-sample.names.original.names
		colnames(trust.perc.table)		<-c(snp.names,"n")	
		trust.perc.table		
	# . other.results:
		error.perc.table				<-trust.perc.table;			
		trust.perc.sd.table				<-trust.perc.table;	
		error.perc.sd.table				<-trust.perc.table;	
	# . Unique and predicted:
		Predicted.pos.nr.table			<-trust.perc.table;		
		Unique.pos.nr.table				<-trust.perc.table;	
		Unique.pos.perc.table			<-trust.perc.table;	
	# . remaining.data.set (i.e. common for remaing set but not found in our one sample)
		Predicted.pos.nr.when.sample.is.removed						<-trust.perc.table;
		Common.pos.nr.when.sample.is.removed							<-trust.perc.table;	
		Common.pos.perc.when.sample.is.removed						<-trust.perc.table;
	# . Missing.data:
		Missing.pos.nr.which.are.common.for.rest.of.samples			<-trust.perc.table;
		Not.Missing.pos.nr.which.are.common.for.rest.of.samples		<-trust.perc.table;
		Missing.pos.perc.which.are.common.for.rest.of.samples		<-trust.perc.table;
		Not.Missing.pos.perc.which.are.common.for.rest.of.samples	<-trust.perc.table;	
	# . New validation data added at 2015.07.20 - done also earlie but separately
		Perc.of.NoNA_valid											<-trust.perc.table		
		Perc.of.NA_valid												<-trust.perc.table
		# individual is a stat calculating how ofthen i have na in my sample, when it is something in the other samples, check one by one
		# global is when we compare our sample to collective infomration from all other samples
		NA_in_s.set.when.noNA_in.sample								<-trust.perc.table	
		noNA_in_s.set.when.noNA_in.sample							<-trust.perc.table	
	
	###	part.2 - Validate.samples:	###
	sample.to.validate<-2
	for(sample.to.validate in 1:ncol(n.ar)){				
	### small trick to have nice data with 1,2, or mnore replicates:
	### i.e. it creates double sample when I work only with one remaining sample
	### it is only temporarly and it will be removed
	### and only for my convenience in script vriting		
		if(ncol(n.ar)[]>2){
			sample.out.ar		<-n.ar[,sample.to.validate,,drop=FALSE]	;dim(sample.out.ar)
			sample.rest.ar		<-n.ar[,-sample.to.validate,,drop=FALSE]	;dim(sample.rest.ar)	
		}else if(ncol(n.ar)[]==2){
			sample.out.ar		<-n.ar[,sample.to.validate,,drop=FALSE]	;dim(sample.out.ar)
			sample.rest.ar		<-n.ar[,-sample.to.validate,,drop=FALSE]	;dim(sample.rest.ar)		
			# -
			t.ar.dim								<-dim(sample.out.ar)
			t.ar.dim[2]							<-t.ar.dim[2]+1
			t.ar									<-array(dim=t.ar.dim)
			t.ar[,1,]							<-sample.rest.ar
			t.ar[,2,]							<-sample.rest.ar
			dimnames(n.ar)[[3]]					<-c(snp.names,"n")	
			sample.rest.ar						<-t.ar
		}else if(ncol(n.ar)[]==1){
			sample.out.ar		<-n.ar[,sample.to.validate,,drop=FALSE]	;dim(sample.out.ar)
			sample.rest.ar		<-n.ar[,sample.to.validate,,drop=FALSE]	;dim(sample.rest.ar)		
			# -
			t.ar.dim								<-dim(sample.out.ar)
			t.ar.dim[2]							<-t.ar.dim[2]+1
			t.ar									<-array(dim=t.ar.dim)
			t.ar[,1,]							<-sample.rest.ar
			t.ar[,2,]							<-sample.rest.ar
			dimnames(n.ar)[[3]]					<-c(snp.names,"n")	
			sample.rest.ar						<-t.ar
		}else{}
		#trust.and.error.mat - for.one.validated.sample
		if(ncol(n.ar)[]>2){
			trust.mat.perc			<-matrix(rep(0,dim(sample.out.ar)[3]*(ncol(n.ar)-1)),ncol=6)
			colnames(trust.mat.perc)<-c(snp.names,"n")	
			error.mat.perc			<-trust.mat.perc
			na.value.mat				<-trust.mat.perc
			No.na.value.mat			<-trust.mat.perc	
			NA_set.badness			<-trust.mat.perc	# i.e. how often snp from my site finds na in the rest of samples
			noNA_set.badness		<-trust.mat.perc	
		}else{
			trust.mat.perc			<-matrix(rep(0,dim(sample.out.ar)[3]*2),ncol=6)
			colnames(trust.mat.perc)<-c(snp.names,"n")	
			error.mat.perc			<-trust.mat.perc
			na.value.mat				<-trust.mat.perc
			No.na.value.mat			<-trust.mat.perc
			NA_set.badness			<-trust.mat.perc	# i.e. how often snp from my site finds na in the rest of samples
			noNA_set.badness			<-trust.mat.perc}
						
		# . Cross.validation.loop . the most important part:)
		sample.in.number			<-2
		for(sample.in.number in 1:ncol(sample.rest.ar)){	
			one.sample.rest.ar	<-sample.rest.ar[,sample.in.number,,drop=FALSE]; dim(one.sample.rest.ar)	; dim(sample.rest.ar)	
			# this llop takes each time, one sample from remaining sample set 
			# and compares it to sample which we are validating														
			#.	trust and error calculation:
				### first we have to fiind observed and predicted positions
				observed.pos.in.one.sample.rest.ar								<-one.sample.rest.ar
				observed.pos.in.one.sample.rest.ar[which(sample.out.ar[,,]==0)]	<-0; dim(observed.pos.in.one.sample.rest.ar)	
				predicted.pos.nr													<-rep(0,dim(sample.out.ar)[3])
				observed.pos.nr													<-rep(0,dim(sample.out.ar)[3])
				### for each snp
				for(trd.dim in 1:dim(sample.out.ar)[3]){
					predicted.pos.nr[trd.dim]									<-sum(sample.out.ar[,,trd.dim,drop=FALSE])
					observed.pos.nr[trd.dim]										<-sum(observed.pos.in.one.sample.rest.ar[,,trd.dim,drop=FALSE])
					}
				### trust:
				trust.perc														<-observed.pos.nr/predicted.pos.nr
				trust.mat.perc[sample.in.number,]								<-trust.perc
				###	error:
				error.perc														<-(predicted.pos.nr-observed.pos.nr)/predicted.pos.nr
				error.mat.perc[sample.in.number,]								<-error.perc
				error.mat.perc[which(as.numeric(predicted.pos.nr)[]==0)]			<-as.numeric(0)
				###	predicted.pos.nr
				Predicted.pos.nr.table[sample.to.validate,]						<-predicted.pos.nr									
			#.	Sample quality
				###	by checking how many positions with info in other samples is paired with na iin our validated sample
				### I dont want to do this for each snp
				### but the loop exits already:)		
				### fist done for each snp indiviually
					rest.sample.NoNA													<-one.sample.rest.ar[,,1:4,drop=TRUE]; dim(rest.sample.NoNA)
					rest.sample.NoNA.silenced										<-one.sample.rest.ar[,,1:4,drop=TRUE]; dim(rest.sample.NoNA)
					valid.sample.na.pos												<-as.numeric(sample.out.ar[,,6,drop=TRUE])
					rest.sample.NoNA.silenced[which(valid.sample.na.pos[]==1),]		<-0								
					opposite.predicted.pos.nr										<-apply(rest.sample.NoNA,2,function(x){sum(x)})
					opposite.predicted.pos.nr.with.info.at.sample.out				<-apply(rest.sample.NoNA.silenced,2,function(x){sum(x)})
					pred.op															<-opposite.predicted.pos.nr
					found.in														<-opposite.predicted.pos.nr.with.info.at.sample.out								
				### secondly done for na in a sample agains any snp in other sample 
					#sum(rest.sample.NoNA[which(rest.sample.NoNA[,]==1)])
					all.informative.pos.opp.s										<-as.numeric(apply(rest.sample.NoNA,1,function(x){sum(x)})[]>0)
					na.pos.s.out													<-as.numeric(sample.out.ar[,,6,drop=TRUE])
					pred.op.tot														<-sum(all.informative.pos.opp.s)
					all.informative.pos.opp.s[which(na.pos.s.out[]==1)]				<-0					
					found.in.tot													<-sum(all.informative.pos.opp.s)
				### combine:		
					found.in.all													<-c(found.in,1,found.in.tot)
					pred.op.all														<-c(pred.op,1,pred.op.tot)
					# I put 1 because i dont make any predictions of this kind for zeros :D
				###	perc.
					perc.of.op.pred.NoNA												<-found.in.all/pred.op.all									
					perc.of.op.pred.NoNA[which(pred.op.all[]==0)]					<-0
				### 
					na.value.mat	[sample.in.number,]								<-c(1-perc.of.op.pred.NoNA)	
					No.na.value.mat[sample.in.number,]								<-perc.of.op.pred.NoNA									
			#. 	Newest part
				#### I am comparing how often snp form my validated sample found any information within a dataset
					one.s.out.NoNA	<-as.numeric(as.numeric(sample.out.ar[,,6])[]==0)
					one.s.rest.NA	<-as.numeric(as.numeric(one.sample.rest.ar[,,6])[]==1)
					res.for.na		<-1-(((sum(one.s.out.NoNA)-sum((one.s.out.NoNA+one.s.rest.NA)[]==2))/sum(one.s.out.NoNA)))
					res.for.NoNa	<-1-res.for.na
									NA_set.badness[sample.in.number,]<-rep(res.for.na,6)
									noNA_set.badness[sample.in.number,]<-rep(res.for.NoNa,6)	
			}
			
		# . trust; error, which was calculated individually
			trust.perc.table[sample.to.validate,]		<-apply(trust.mat.perc,2,function(x){mean(x)})	
			error.perc.table[sample.to.validate,]		<-apply(error.mat.perc,2,function(x){mean(x)})		
			Perc.of.NoNA_valid[sample.to.validate,]		<-apply(No.na.value.mat,2,function(x){mean(x)})	
			Perc.of.NA_valid[sample.to.validate,]		<-apply(na.value.mat,2,function(x){mean(x)})	
			trust.perc.sd.table[sample.to.validate,]		<-apply(trust.mat.perc,2,function(x){sd(x)})	
			error.perc.sd.table[sample.to.validate,]		<-apply(error.mat.perc,2,function(x){sd(x)})		
			NA_in_s.set.when.noNA_in.sample[sample.to.validate,]			<-apply(NA_set.badness,2,function(x){mean(x)})
			noNA_in_s.set.when.noNA_in.sample[sample.to.validate,]		<-apply(noNA_set.badness,2,function(x){mean(x)})
		
		# . sample.sum:		
			remaining.s.nr		<-dim(sample.rest.ar)[2]
			sum.sample.rest.ar	<-sample.rest.ar[,1,,drop=FALSE]; dim(sum.sample.rest.ar)
			for(trd.d in 2:ncol(sample.rest.ar)){sum.sample.rest.ar<-sum.sample.rest.ar+sample.rest.ar[,trd.d,,drop=FALSE]}#sum.samples.not.snp's!
			sum.sample.rest.ar[which(sum.sample.rest.ar[,,]>0)]<-1	
		# . unique.positions:
			#unique.pos.nr											<-rep(0,dim(sample.out.ar)[3])
			#sample.out.ar.unique										<-sample.out.ar
			#sample.out.ar.unique[which(sum.sample.rest.ar[,,]>0)]	<-0
			#for(trd.dim in 1:dim(sample.out.ar.unique)[3]){
			#	unique.pos.nr[trd.dim]<-sum(sample.out.ar.unique[,,trd.dim,drop=FALSE])}
			#Unique.pos.nr.table[sample.to.validate,]					<-unique.pos.nr	
			#Unique.pos.perc.table[sample.to.validate,]				<-unique.pos.nr/predicted.pos.nr

		# . remaining.data.set (i.e. common for remaing set but not found in our one sample)
			predicted.rest.ar.pos.nr	<-rep(0,dim(sample.out.ar)[3])
			common.rest.ar.pos.nr		<-rep(0,dim(sample.out.ar)[3])
			# -
			predicted.rest.ar												<-sum.sample.rest.ar	
			predicted.rest.ar[which(sum.sample.rest.ar[,,]>0)]				<-1	
			common.rest.ar													<-sum.sample.rest.ar				
			common.rest.ar[which(sum.sample.rest.ar[,,]<remaining.s.nr)]		<-0
			common.rest.ar[which(sum.sample.rest.ar[,,]==remaining.s.nr)]	<-1		
			# -	
			for(trd.d in 1:dim(predicted.rest.ar)[3]){	predicted.rest.ar.pos.nr[trd.d]	<-sum(predicted.rest.ar[,,trd.d,drop=FALSE])}		
			for(trd.d in 1:dim(common.rest.ar)[3]){		common.rest.ar.pos.nr[trd.d]	<-sum(common.rest.ar[,,trd.d,drop=FALSE])}	
			# -			
			Predicted.pos.nr.when.sample.is.removed[sample.to.validate,]		<-predicted.rest.ar.pos.nr;
			Common.pos.nr.when.sample.is.removed	[sample.to.validate,]		<-common.rest.ar.pos.nr;	
			Common.pos.perc.when.sample.is.removed[sample.to.validate,]		<-common.rest.ar.pos.nr/predicted.rest.ar.pos.nr		
		
					
			#unique continuation:	
			unique.pos.nr											<-rep(0,dim(sample.out.ar)[3])
			sample.out.ar.unique										<-sample.out.ar
			sample.out.ar.unique[which(sum.sample.rest.ar[,,]>0)]	<-0
			for(trd.dim in 1:dim(sample.out.ar.unique)[3]){
				unique.pos.nr[trd.dim]<-sum(sample.out.ar.unique[,,trd.dim,drop=FALSE])}
			Unique.pos.nr.table[sample.to.validate,]					<-unique.pos.nr	
			
			if(unique.count.protocol[]==1){	
				total.al	<-c(sum(sum.sample.rest.ar[,,1:4])+Common.pos.nr.when.sample.is.removed[sample.to.validate,5]+sum(unique.pos.nr[1:4]))
				Unique.pos.perc.table[sample.to.validate,]				<-unique.pos.nr/	total.al
			}else if(unique.count.protocol[]==2){				
				total.al	<-c(apply(sum.sample.rest.ar[,,1:4],2,function(x){sum(x)})+unique.pos.nr[1:4]+Common.pos.nr.when.sample.is.removed[sample.to.validate,5])
				total.al<-c(total.al,1,1)
				Unique.pos.perc.table[sample.to.validate,]				<-unique.pos.nr/	total.al
			}else{}
			# and finally, no predictions for na and z
			Unique.pos.perc.table[sample.to.validate,c(5,6)]			<-0
			
		
		# . Missing.data:			
			not.missing.common.pos.in.s.out.with.common.pos.rest.nr					<-rep(0,dim(sample.out.ar)[3])		
			# -		
			not.missing.common.pos.in.s.out.with.common.pos.rest								<-sample.out.ar
			not.missing.common.pos.in.s.out.with.common.pos.rest[which(common.rest.ar[,,]==0)]	<-0
			# -
			for(trd.d in 1:dim(not.missing.common.pos.in.s.out.with.common.pos.rest)[3]){		
				not.missing.common.pos.in.s.out.with.common.pos.rest.nr[trd.d]<-sum(not.missing.common.pos.in.s.out.with.common.pos.rest[,,trd.d,drop=FALSE])}			
			# -
			missing.common.pos.in.s.out.with.common.pos.rest.nr		<-common.rest.ar.pos.nr-not.missing.common.pos.in.s.out.with.common.pos.rest.nr
			# -
			Not.Missing.pos.nr.which.are.common.for.rest.of.samples[sample.to.validate,]	<-not.missing.common.pos.in.s.out.with.common.pos.rest.nr
			Missing.pos.nr.which.are.common.for.rest.of.samples[sample.to.validate,]		<-missing.common.pos.in.s.out.with.common.pos.rest.nr
			# -
			Not.Missing.pos.perc.which.are.common.for.rest.of.samples[sample.to.validate,]	<-not.missing.common.pos.in.s.out.with.common.pos.rest.nr/common.rest.ar.pos.nr
			Missing.pos.perc.which.are.common.for.rest.of.samples[sample.to.validate,]		<-missing.common.pos.in.s.out.with.common.pos.rest.nr/common.rest.ar.pos.nr			
		}#for(sample.to.validate 	
		# .....................................
		
		# some empty elements:
		Conditional.trust	<-rep(0,6)								#this form is better to be done
		Conditional.NA		<-rep(0,6)	
		
		# ... list:		
		one.isolate.cross.validation.results<-list(
			trust.perc.table,
			error.perc.table,						
			trust.perc.sd.table,				
			error.perc.sd.table,				
			# . Unique and predicted:
			Predicted.pos.nr.table,				
			Unique.pos.nr.table,				
			Unique.pos.perc.table,
			# . remaining.data.set
			Predicted.pos.nr.when.sample.is.removed,	
			Common.pos.nr.when.sample.is.removed,
			Common.pos.perc.when.sample.is.removed,
			# . Missing.data:
			Missing.pos.nr.which.are.common.for.rest.of.samples,
			Not.Missing.pos.nr.which.are.common.for.rest.of.samples,
			Missing.pos.perc.which.are.common.for.rest.of.samples,
			Not.Missing.pos.perc.which.are.common.for.rest.of.samples,
			# . gamma protocol
			Perc.of.NoNA_valid,
			Perc.of.NA_valid,
			NA_in_s.set.when.noNA_in.sample,
			noNA_in_s.set.when.noNA_in.sample,
			Conditional.trust,
			Conditional.NA)			
		names(one.isolate.cross.validation.results)<-c(	
			"trust.perc.table",
			"error.perc.table",						
			"trust.perc.sd.table",				
			"error.perc.sd.table",				
			# . Unique and predicted:
			"Predicted.pos.nr.table",				
			"Unique.pos.nr.table",				
			"Unique.pos.perc.table",
			# . remaining.data.set
			"Predicted.pos.nr.when.sample.is.removed",	
			"Common.pos.nr.when.sample.is.removed",
			"Common.pos.perc.when.sample.is.removed",
			# . Missing.data:
			"Missing.pos.nr.which.are.common.for.rest.of.samples",
			"Not.Missing.pos.nr.which.are.common.for.rest.of.samples",
			"Missing.pos.perc.which.are.common.for.rest.of.samples",
			"Not.Missing.pos.perc.which.are.common.for.rest.of.samples",
			# . gamma protocol
			"Perc.of.NoNA_valid",
			"Perc.of.NA_valid",
			"NA_in_s.set.when.noNA_in.sample",
			"noNA_in_s.set.when.noNA_in.sample",
			"Conditional.trust",
			"Conditional.NA")
		
	# ... remove.problem.withdiv.by.zero(Inf)	and zer div by zero (NaN)		
	for(table.nr in 1:length(one.isolate.cross.validation.results)){
		# temp.table<-(Samples.Cross.validation.Six.Categories[[4]])[[14]]
		#table.nr<-1; a<-temp.table; a[,]<-0; temp.table<-temp.table/a
		temp.table<-one.isolate.cross.validation.results[[table.nr]]
		temp.table[which(temp.table[]==Inf)]<-0
		temp.table[is.na(temp.table)]<-0
		mode(temp.table)<-"numeric"
		one.isolate.cross.validation.results[[table.nr]]<-temp.table
		}
			#and finally we are caylculating the stuff used in cross validation:
				est.var								<-one.isolate.cross.validation.results
				trust								<-1
				unique								<-7
				noNA.perc							<-15
				NA.perc								<-16
				NA_in_s.set.when.noNA_in.sample		<-17
				###
				#Conditional.trust	<-est.var[[trust]]-(est.var[[trust]]*est.var[[NA.perc]])
				Conditional.trust	<-est.var[[trust]]*est.var[[noNA.perc]]							#this form is better to be done
				#and:
				Conditional.NA		<-est.var[[unique]]-(est.var[[unique]]*est.var[[NA_in_s.set.when.noNA_in.sample]])		
				###
				one.isolate.cross.validation.results[[19]]	<-Conditional.trust
				one.isolate.cross.validation.results[[20]]	<-Conditional.NA	
				
	# .....................................
	return(one.isolate.cross.validation.results)
	# .....................................
}
























###	---------------------------------------------------------------------------------------------------
###	12. GAMMA_collapse_acgtz_arrays_function20150725		
###	---------------------------------------------------------------------------------------------------
GAMMA_collapse_acgtz_arrays_function20150725<-function(
											Ln.x.Wg.arr,
											pb.x.Ln.x.Wg.arr,
											z.cov.mat,
											group.names,
											group.pos,
											min.for.p){												
			# deveoper
			# -------------------------------------------------------------------
			#Ln.x.Wg.arr		<-ARRAY_log10.x.Weigth.div.BetaCorr_ACGTZ
			#pb.x.Ln.x.Wg.arr<-ARRAY_pbinom.x.log10.x.Weight.div.BetaCorr_ACGTZ
			#group.names		<-collapse.group.names
			#group.pos			<-collapse.group.pos										
			# -------------------------------------------------------------------
			
# - 1 - empty.elements:
	ARRAY_collapsed_pbinom_ACGTZ	<-pb.x.Ln.x.Wg.arr[,1:length(group.names),,drop=FALSE]	
	ARRAY_collapsed_pbinom_ACGTZ[,,]<-0
	colnames(ARRAY_collapsed_pbinom_ACGTZ)<-names(group.names)
	#dimnames(ARRAY_collapsed_pbinom_ACGTZ)
	ARRAY_CollWg_ACGTZ				<-ARRAY_collapsed_pbinom_ACGTZ
	ARRAY_CollWg.x.pbinom_ACGTZ		<-ARRAY_collapsed_pbinom_ACGTZ
	Coll_ZERO_pbinom_matrix			<-ARRAY_collapsed_pbinom_ACGTZ[,,5,drop=FALSE]


# - 2 - ARRAY_CollWg_ACGTZ; s.nr<-1; d.dim<-1
	for(s.nr in 1:length(group.names)){
		one.gr.pos	<-group.pos[[s.nr]]
		temp.arr	<-Ln.x.Wg.arr[,which(one.gr.pos[]==1),,drop=FALSE]#dim(ARRAY_CollWg_ACGTZ); dim(temp.arr)
		for(d.dim in 1:4){
			temp.mat<-temp.arr[,,d.dim,drop=FALSE]; dim(temp.mat)
			ARRAY_CollWg_ACGTZ[,s.nr,d.dim]<-apply(temp.mat,1,function(x){sum(x)})}												
		}				

# - 3 - ARRAY_CollWg.x.pbinom_ACGTZ; s.nr<-1; d.dim<-1
	for(s.nr in 1:length(group.names)){
		one.gr.pos	<-group.pos[[s.nr]]
		temp.arr	<-pb.x.Ln.x.Wg.arr[,which(one.gr.pos[]==1),,drop=FALSE]#dim(ARRAY_CollWg_ACGTZ); dim(temp.arr)
		for(d.dim in 1:4){
			temp.mat<-temp.arr[,,d.dim,drop=FALSE]; dim(temp.mat)
			ARRAY_CollWg.x.pbinom_ACGTZ[,s.nr,d.dim]<-apply(temp.mat,1,function(x){sum(x)})}										
		}													
	#aa<-ARRAY_CollWg.x.pbinom_ACGTZ; colnames(aa)<-1:ncol(aa); rownames(aa)<-1:nrow(aa); round(aa,digits=3)
	#summary(aa)

# - 4 - Coll_ZERO_pbinom_matrix	
	for(s.nr in 1:length(group.names)){
		one.gr.pos	<-group.pos[[s.nr]]
		repl.nr		<-sum(as.numeric(group.pos[[s.nr]]))
		temp.arr	<-z.cov.mat[,which(one.gr.pos[]==1),,drop=FALSE]#dim(z.cov.mat); dim(temp.arr)
		sum.temp.arr<-as.numeric((apply(temp.arr,1,function(x){sum(x)}))[]==repl.nr)
		Coll_ZERO_pbinom_matrix[,s.nr,]<-sum.temp.arr
		}
	#aa<-Coll_ZERO_pbinom_matrix; colnames(aa)<-1:ncol(aa); rownames(aa)<-1:nrow(aa); round(aa,digits=3)
	
# - 5 - ARRAY_collapsed_pbinom_ACGTZ
	for(i.dim in 1:4){
		p.temp	<-ARRAY_CollWg.x.pbinom_ACGTZ[,,i.dim,drop=FALSE]; dim(p.temp)
		w.temp	<-ARRAY_CollWg_ACGTZ[,,i.dim,drop=FALSE]; dim(w.temp)
		for(c.dim in 1:ncol(p.temp)){
			c.p.temp<-p.temp[,c.dim,,drop=FALSE]; dim(c.p.temp)
			c.w.temp<-w.temp[,c.dim,,drop=FALSE]; dim(c.w.temp)
			ARRAY_collapsed_pbinom_ACGTZ[,c.dim,i.dim]<-as.numeric(c.p.temp/c.w.temp)
			ARRAY_collapsed_pbinom_ACGTZ[which(c.p.temp[]==0),c.dim,i.dim]<-0
			ARRAY_collapsed_pbinom_ACGTZ[which(c.w.temp[]==0),c.dim,i.dim]<-0}
	}		
	ARRAY_collapsed_pbinom_ACGTZ[,,5]<-Coll_ZERO_pbinom_matrix
	#aa<-ARRAY_collapsed_pbinom_ACGTZ; colnames(aa)<-1:ncol(aa); rownames(aa)<-1:nrow(aa); round(aa,digits=3)
	#summary(aa)
	#round(aa[1356,,],digits=2)
# - 5 - final correction
	ARRAY_collapsed_pbinom_ACGTZ[which(ARRAY_collapsed_pbinom_ACGTZ[,,]==0)]<-as.numeric(min.for.p)	
# - 6 - return	
	return(ARRAY_collapsed_pbinom_ACGTZ)
}



###	---------------------------------------------------------------------------------------------------
###	13. New_INFO_for_ollapse_acgtz_arrays__function_20150725		
###	---------------------------------------------------------------------------------------------------
# typically:  use.as.separator<-";" 

New_INFO_for_Collapse_acgtz_arrays__function_20160601<-function(
											pbinom.arr,
											cov.arr	,
											gen.cov.arr,
											z.cov.mat,
											group.names,
											group.pos,
											use.as.separator.in.info.arr
											){												
use.as.separator	<-use.as.separator.in.info.arr# because i was always forgetting where it is comming from
use.as.separator<-";"	
# - 1 - empty.elements:
		#:..array:					
			ar.colnames				<-names(group.names)
			ar.rownames				<-rownames(pbinom.arr)
			ar.depthnam				<-c(dimnames(pbinom.arr)[[3]],"snp.presence","snp.freq","gen.cov","group.name","repl.nr","sample.names")	
			ar.dimnames				<-list(ar.rownames,ar.colnames,ar.depthnam)	
			ar.dims					<-c(length(ar.rownames),length(ar.colnames),length(ar.depthnam))	
			Empty_matrix				<-matrix(nrow=ar.dims[1],ncol=ar.dims[2])
			Empty_matrix[,]			<-0
			mode(Empty_matrix)		<-"numeric" 
			array_matrices			<-c(rep(Empty_matrix,ar.dims[3]))
			ar						<-array(array_matrices,dim=ar.dims)
			dimnames(ar)				<-ar.dimnames	
			mode(ar)					<-"numeric"	;#dim(ar)			
		#:..new.empty.arrays		
			ARRAY_collapsed_info_ACGT	<-ar
			rm(ar)		
# - 2 - pbinom: old version
	#	ARRAY_collapsed_info_ACGT[,,1:5]<-round(pbinom.arr[,,,drop=FALSE],digits=2)
# - 2 - pbinom: 2016.06.01, no round ! 
	ARRAY_collapsed_info_ACGT[,,1:5]<-pbinom.arr[,,,drop=FALSE]
	
# - 3 - gen.cov:
	for(s.nr in 1:length(group.names)){
		one.gr.pos	<-group.pos[[s.nr]]
		temp.arr	<-gen.cov.arr[,which(one.gr.pos[]==1),,drop=FALSE]; dim(temp.arr)
		for(n.locus in 1:nrow(temp.arr)){
			cov.locus	<-as.vector(temp.arr[n.locus,,]); cov.locus
			cov.locus	<-paste(cov.locus,collapse=use.as.separator); cov.locus
			ARRAY_collapsed_info_ACGT[n.locus,s.nr,8]<-cov.locus
			}
		}	
# - 4 -"snp.presence":
	SNP.names	<-c("A","C","G","T","z")
	cov.arr[,,5]<-z.cov.mat
	for(s.nr in 1:length(group.names)){
		one.gr.pos	<-group.pos[[s.nr]]
		temp.arr	<-cov.arr[,which(one.gr.pos[]==1),,drop=FALSE]; dim(temp.arr)
		for(n.locus in 1:nrow(temp.arr)){
			one.locus.arr								<-temp.arr[n.locus,,,drop=FALSE];one.locus.arr
			one.locus.arr[which(one.locus.arr[,,]>0)]	<-1;one.locus.arr
			one.locus.sum								<-rep(0,ncol(one.locus.arr))
			for(i.dim in 1:ncol(one.locus.arr)){
				el	<-one.locus.arr[,i.dim,]
				if(sum(as.numeric((el[]>0)))[]>0){one.locus.sum[i.dim]<-paste(SNP.names[which(el[]>0)],collapse="/")
				}else{one.locus.sum[i.dim]<-"-"}
			}
			one.locus.sum<-paste(one.locus.sum,collapse=use.as.separator)	
			ARRAY_collapsed_info_ACGT[n.locus,s.nr,6]<-one.locus.sum	
		}		
		}				
# - 5 -"snp.freq":
	SNP.names	<-c("A","C","G","T","z")
	cov.arr[,,5]<-z.cov.mat
	for(s.nr in 1:length(group.names)){
		one.gr.pos	<-group.pos[[s.nr]]
		temp.arr	<-cov.arr[,which(one.gr.pos[]==1),,drop=FALSE]; dim(temp.arr)
		for(n.locus in 1:nrow(temp.arr)){
			one.locus.arr								<-temp.arr[n.locus,,,drop=FALSE];one.locus.arr
			one.locus.arr[which(one.locus.arr[,,5]>0)]	<-"z"
			one.locus.sum								<-rep(0,ncol(one.locus.arr))
			for(i.dim in 1:ncol(one.locus.arr)){
				el	<-one.locus.arr[,i.dim,]	
				if(sum(as.numeric((el[]>0)))[]>0){
					one.locus.sum[i.dim]<-paste(el[which(el[]!=0)],collapse="/")
				}else{one.locus.sum[i.dim]<-"-"}
			}
			one.locus.sum[which(one.locus.sum[]=="z/1")]<-"z"
			one.locus.sum<-paste(one.locus.sum,collapse=use.as.separator)	
			ARRAY_collapsed_info_ACGT[n.locus,s.nr,7]<-one.locus.sum	
		}		
		}				
	#ARRAY_collapsed_info_ACGT[1356,,]
	return(ARRAY_collapsed_info_ACGT)
}




























###	---------------------------------------------------------------------------------------------------
###	14. Collapse_array__Conservative_method		2880		
###	---------------------------------------------------------------------------------------------------
Collapse_array__Conservative_method<-function(
		bin.arr,		
		cov.arr,				
		beta.corr.arr,	
		filt.cov.arr,		
		z.cov.mat,			
		group.names,		
		group.pos,			
		min.for.p,			
		cont.pos,			
		nz.pos){			
																		
	# - 1 - empty.elements:
	ARRAY_collapsed_pbinom_ACGTZ				<-cov.arr[,1:length(group.names),,drop=FALSE]	
	ARRAY_collapsed_pbinom_ACGTZ[,,]		<-0
	colnames(ARRAY_collapsed_pbinom_ACGTZ)	<-names(group.names)
	dimnames(ARRAY_collapsed_pbinom_ACGTZ)
	Coll_ZERO_pbinom_matrix					<-ARRAY_collapsed_pbinom_ACGTZ[,,5,drop=FALSE]
	
	# - 2 speciffic for conservative method	
	s.nr	<-1
	for(s.nr in 1:length(group.names)){
			one.gr.pos			<-group.pos[[s.nr]]
			col.to.take			<-which(one.gr.pos[]==1)[1]# because we take only the same elemenets
			const.loci.in.gr		<-cont.pos[[s.nr]]
			####
			temp					<-bin.arr[,col.to.take,,drop=FALSE]; dim(temp)
			temp[which(const.loci.in.gr[]==0),,]<-0			#to remove all variable positions
			####
			ARRAY_collapsed_pbinom_ACGTZ	[,s.nr,]<-temp	
			}
	#	apply(ARRAY_collapsed_pbinom_ACGTZ,2,function(x){sum(x)})

	# - 3 - Coll_ZERO_pbinom_matrix	
	for(s.nr in 1:length(group.names)){
		one.gr.pos	<-group.pos[[s.nr]]
		repl.nr		<-sum(as.numeric(group.pos[[s.nr]]))
		temp.arr	<-z.cov.mat[,which(one.gr.pos[]==1),,drop=FALSE]#dim(z.cov.mat); dim(temp.arr)
		sum.temp.arr<-as.numeric((apply(temp.arr,1,function(x){sum(x)}))[]==repl.nr)
		Coll_ZERO_pbinom_matrix[,s.nr,]<-sum.temp.arr
		}

	# - 4 - final correction		
	ARRAY_collapsed_pbinom_ACGTZ[,,5]<-Coll_ZERO_pbinom_matrix
	ARRAY_collapsed_pbinom_ACGTZ[which(ARRAY_collapsed_pbinom_ACGTZ[,,]==0)]<-as.numeric(min.for.p)	
		#aa<-Coll_ZERO_pbinom_matrix; colnames(aa)<-1:ncol(aa); rownames(aa)<-1:nrow(aa); round(aa,digits=3)
		#aa<-ARRAY_collapsed_pbinom_ACGTZ; colnames(aa)<-1:ncol(aa); rownames(aa)<-1:nrow(aa); round(aa,digits=3)
		#summary(aa)
		#round(aa[1356,,],digits=2)
	
	# - 5 - return	
	return(ARRAY_collapsed_pbinom_ACGTZ)
	}























###	---------------------------------------------------------------------------------------------------
###	15. Collapse_array__Site_elimination_method		2960		
###	---------------------------------------------------------------------------------------------------
Collapse_array__Site_elimination_method<-function(
		bin.arr,		
		cov.arr,				
		beta.corr.arr,	
		filt.cov.arr,		
		z.cov.mat,			
		group.names,		
		group.pos,			
		min.for.p,			
		cont.pos,			
		nz.pos,
		max.na.perc,
		max.na.pos.nr){			
																		
	# - 1 - empty.elements:
	ARRAY_collapsed_pbinom_ACGTZ				<-cov.arr[,1:length(group.names),,drop=FALSE]	
	ARRAY_collapsed_pbinom_ACGTZ[,,]		<-0
	colnames(ARRAY_collapsed_pbinom_ACGTZ)	<-names(group.names)
	dimnames(ARRAY_collapsed_pbinom_ACGTZ)
	Coll_ZERO_pbinom_matrix					<-ARRAY_collapsed_pbinom_ACGTZ[,,5,drop=FALSE]
	
	# - 2 speciffic for site elimination method	(by na admixture level)
	ARRAY_CovACGTZ							<-cov.arr
	mode(ARRAY_CovACGTZ)						<-"numeric"
	cover.mat								<-ARRAY_CovACGTZ[,,1]+ARRAY_CovACGTZ[,,2]+ARRAY_CovACGTZ[,,3]+ARRAY_CovACGTZ[,,4]					
	zdl.mat									<-cover.mat
	zdl.mat[which(cover.mat[,]>0)]			<-0
	zdl.mat[which(cover.mat[,]==0)]			<-1
	
	
	s.nr	<-1
	for(s.nr in 1:length(group.names)){
		 	one.gr.pos			<-group.pos[[s.nr]]
			one.gr.pos
			if(sum(as.numeric(one.gr.pos))[]==1){
				#taken directly from conservative method
				col.to.take				<-which(one.gr.pos[]==1)[1]# because we take only the same elemenets
				const.loci.in.gr		<-cont.pos[[s.nr]]
				temp					<-bin.arr[,col.to.take,,drop=FALSE]; dim(temp)
				temp[which(const.loci.in.gr[]==0),,]<-0			#to remove all variable positions
				ARRAY_collapsed_pbinom_ACGTZ	[,s.nr,]<-temp	
			}else{
				### find how many z's you have per locus
					nz.mat									<-zdl.mat[,which(one.gr.pos[]==1)]
					mode(nz.mat)								<-"numeric"			
					nz.per.locus								<-apply(nz.mat,1,function(x){sum(x)})					
					nz.per.locus		
					###
					max.na.nr								<-as.numeric(max.na.perc*sum(as.numeric(one.gr.pos)))				
					max.na.nr
					if(max.na.pos.nr[]<=max.na.nr){						
					}else{max.na.nr<-max.na.pos.nr}
					###
					tt										<-nz.per.locus				
					take.loci								<-as.numeric(which(nz.per.locus[]<=max.na.nr));length(take.loci)
					dont.take.loci							<-as.numeric(which(nz.per.locus[]>max.na.nr));length(dont.take.loci)# be careful, it is sth differen then in fixed tr method
				###
				temp.arr								<-bin.arr[,which(one.gr.pos[]==1),,drop=FALSE]; dim(temp.arr); #sum(temp.arr)
				temp.arr[dont.take.loci,,]			<-0# eliminates all the info from unwanted loci
				
				### secondary filter; speciffic to this method
				### i.e. remve the positions that have more than na in a given nr of replicates
				###	it was removed  - for convenience
				
				###
				snp.bin.mat							<-matrix(ncol=5,nrow=nrow(temp.arr))
				snp.bin.mat[,]						<-0
				mode(snp.bin.mat)					<-"numeric"
				# last i.dim if for z
				###
				i.dim<-1
				for(i.dim in 1:4){
						i.temp.arr			<-temp.arr[,,i.dim,drop=TRUE]; dim(i.temp.arr)
					 	snp.pos				<-as.numeric( apply(i.temp.arr,1,function(x){sum(x)}) ); snp.pos	
					 	# binary presence
					 	snp.pos[which(snp.pos[]>0)]<-1
					 	#
					 	snp.bin.mat[,i.dim]	<-snp.pos
					 	}	
				### put it into collapsed array
				ARRAY_collapsed_pbinom_ACGTZ[,s.nr,]<-snp.bin.mat				
				}
			}
			#	sum(as.numeric(apply(snp.bin.mat,1,function(x){sum(x)})[]>1))/sum(as.numeric(apply(snp.bin.mat,1,function(x){sum(x)})[]>0))
			#	apply(ARRAY_collapsed_pbinom_ACGTZ,2,function(x){sum(x)})

	# - 3 - Coll_ZERO_pbinom_matrix	
	for(s.nr in 1:length(group.names)){
		one.gr.pos	<-group.pos[[s.nr]]
		repl.nr		<-sum(as.numeric(group.pos[[s.nr]]))
		temp.arr	<-z.cov.mat[,which(one.gr.pos[]==1),,drop=FALSE]#dim(z.cov.mat); dim(temp.arr)
		sum.temp.arr<-as.numeric((apply(temp.arr,1,function(x){sum(x)}))[]==repl.nr)
		Coll_ZERO_pbinom_matrix[,s.nr,]<-sum.temp.arr
		}

	# - 4 - final correction		
	ARRAY_collapsed_pbinom_ACGTZ[,,5]<-Coll_ZERO_pbinom_matrix
	ARRAY_collapsed_pbinom_ACGTZ[which(ARRAY_collapsed_pbinom_ACGTZ[,,]==0)]<-as.numeric(min.for.p)	
		#aa<-Coll_ZERO_pbinom_matrix; colnames(aa)<-1:ncol(aa); rownames(aa)<-1:nrow(aa); round(aa,digits=3)
		#aa<-ARRAY_collapsed_pbinom_ACGTZ; colnames(aa)<-1:ncol(aa); rownames(aa)<-1:nrow(aa); round(aa,digits=3)
		#summary(aa)
		#round(aa[1356,,],digits=2)
	
	# - 5 - return	
	return(ARRAY_collapsed_pbinom_ACGTZ)
	}


















###	---------------------------------------------------------------------------------------------------
###	16. Collapse_array__Fixed_Thresholds_method		3080	 / upgrade at 2016.06.02 3am	
###	---------------------------------------------------------------------------------------------------
Collapse_array__Fixed_Thresholds_method<-function(
		bin.arr,		
		cov.arr,				
		beta.corr.arr,	
		filt.cov.arr,				
		###
		z.cov.mat,			
		group.names,		
		group.pos,			
		min.for.p,			
		cont.pos,			
		nz.pos,
		max.na.perc,
		min.allele.adm,
		min.repl.nr.with.allele){	


	# - 1 - empty.elements:
	ARRAY_collapsed_pbinom_ACGTZ				<-cov.arr[,1:length(group.names),,drop=FALSE]	
	ARRAY_collapsed_pbinom_ACGTZ[,,]		<-0
	colnames(ARRAY_collapsed_pbinom_ACGTZ)	<-names(group.names)
	#dimnames(ARRAY_collapsed_pbinom_ACGTZ)
	Coll_ZERO_pbinom_matrix					<-ARRAY_collapsed_pbinom_ACGTZ[,,5,drop=FALSE]
	
	# - 2 speciffic for site elimination method	(by na admixture level)
	
	ARRAY_CovACGTZ							<-cov.arr
	mode(ARRAY_CovACGTZ)						<-"numeric"
	cover.mat								<-ARRAY_CovACGTZ[,,1]+ARRAY_CovACGTZ[,,2]+ARRAY_CovACGTZ[,,3]+ARRAY_CovACGTZ[,,4]					
	zdl.mat									<-cover.mat
	zdl.mat[which(cover.mat[,]>0)]			<-0
	zdl.mat[which(cover.mat[,]==0)]			<-1
	
	s.nr	<-1
	for(s.nr in 1:length(group.names)){
		 	one.gr.pos			<-group.pos[[s.nr]]
			one.gr.pos
						
			if(sum(as.numeric(one.gr.pos))[]==1){	
				#taken directly from conservative method
				col.to.take								<-which(one.gr.pos[]==1)[1]# because we take only the same elemenets
				const.loci.in.gr						<-cont.pos[[s.nr]]
				#temp									<-bin.arr[,col.to.take,,drop=FALSE]; dim(temp)							
				### first filter. see the min.allele admixture at each locus :)
				###
					temp.cov.arr								<-cov.arr[,which(one.gr.pos[]==1),,drop=FALSE];
					temp.filt.cov.arr						<-filt.cov.arr[,which(one.gr.pos[]==1),,drop=FALSE];				
					temp.arr								<-bin.arr[,which(one.gr.pos[]==1),,drop=FALSE]; dim(temp.arr); #sum(temp.arr)
					temp.arr[,,]								<-0
					###
					i.dim<-1
					for(i.dim in 1:4){
						one.tcr				<-temp.cov.arr[,,i.dim];		
						one.tfcr				<-temp.filt.cov.arr[,,i.dim]	;
						###
						rel.fr				<-one.tcr/one.tfcr
						rel.fr[which(one.tfcr[]==0)]				<-0; #hist(rel.fr,ylim=c(0,50),xlim=c(0.1,1),breaks=50,col="black")		
						rel.fr[which(rel.fr[]<min.allele.adm)]	<-0;	 #hist(rel.fr,ylim=c(0,50),xlim=c(0.1,1),breaks=50,col="yellow",add=TRUE)	
						rel.fr[which(rel.fr[]>0)]				<-1; #hist(rel.fr,ylim=c(0,50),xlim=c(0.1,1),breaks=50,col="darkgray",add=TRUE)
						temp.arr[,,i.dim]						<-rel.fr
						}
						#	sum(temp.arr)
						#	summary(temp.arr)
				ARRAY_collapsed_pbinom_ACGTZ	[,s.nr,]<-temp.arr				
			}else{
				### find how many z's you have per locus
				###
					nz.mat									<-zdl.mat[,which(one.gr.pos[]==1)]
					mode(nz.mat)								<-"numeric"			
					nz.per.locus								<-apply(nz.mat,1,function(x){sum(x)})					
					nz.per.locus		
					###
					max.na.nr								<-as.numeric(max.na.perc*sum(as.numeric(one.gr.pos)))				
					max.na.nr
					###
					tt										<-nz.per.locus				
					take.loci								<-as.numeric(which(nz.per.locus[]<max.na.nr));length(take.loci)
					dont.take.loci							<-as.numeric(which(nz.per.locus[]>=max.na.nr));length(dont.take.loci)										
				### first filter. see the min.allele admixture at each locus :)
				###
					temp.cov.arr								<-cov.arr[,which(one.gr.pos[]==1),,drop=FALSE];
					temp.filt.cov.arr						<-filt.cov.arr[,which(one.gr.pos[]==1),,drop=FALSE];				
					temp.arr								<-bin.arr[,which(one.gr.pos[]==1),,drop=FALSE]; dim(temp.arr); #sum(temp.arr)
					temp.arr[,,]								<-0
					###
					i.dim<-2
					for(i.dim in 1:4){
						one.tcr				<-temp.cov.arr[,,i.dim];			dim(one.tcr)
						one.tfcr			<-temp.filt.cov.arr[,,i.dim]	; 	dim(one.tfcr)
						###
						rel.fr				<-one.tcr/one.tfcr
						rel.fr[which(one.tfcr[,]==0)]			<-0; #hist(rel.fr,ylim=c(0,50),xlim=c(0.1,1),breaks=50,col="black")		
						rel.fr[which(rel.fr[]<min.allele.adm)]	<-0;	 #hist(rel.fr,ylim=c(0,50),xlim=c(0.1,1),breaks=50,col="yellow",add=TRUE)	
						rel.fr[which(rel.fr[]>0)]				<-1; #hist(rel.fr,ylim=c(0,50),xlim=c(0.1,1),breaks=50,col="darkgray",add=TRUE)
						temp.arr[,,i.dim]						<-rel.fr
						}
						#sum(temp.arr)				
				
				###	second filter: loci with more than speciefied perc of missing data
				###
					temp.arr[dont.take.loci,,]			<-0# eliminates all the info from unwanted loci
					#sum(temp.arr)						
				
				###	third filter: loci with minimum number of replicates among samples, but no less then total nr of replicates
				###													
					if(sum(one.gr.pos)<min.repl.nr.with.allele){
							min.a.nr		<-sum(one.gr.pos)
					}else{	min.a.nr		<-min.repl.nr.with.allele}
					### min.repl.nr.with.allele
					snp.bin.mat							<-matrix(ncol=5,nrow=nrow(temp.arr))
					snp.bin.mat[,]						<-0
					mode(snp.bin.mat)					<-"numeric"
					# last i.dim if for z
					###
					i.dim<-1
					for(i.dim in 1:4){
						i.temp.arr			<-temp.arr[,,i.dim,drop=TRUE]; dim(i.temp.arr)
					 	tls					<-as.numeric( apply(i.temp.arr,1,function(x){sum(x)}) ); 
					 	el.pos				<-rep(0,length(tls))
					 	# binary presence
					 	el.pos[which(tls[]>=as.numeric(min.a.nr))]<-1
					 	#
					 	snp.bin.mat[,i.dim]<-el.pos
					 	}	
						
					### put it into collapsed array
					ARRAY_collapsed_pbinom_ACGTZ[,s.nr,]<-snp.bin.mat				
				}#end if.else			
			
			
			}
			#	sum(as.numeric(apply(snp.bin.mat,1,function(x){sum(x)})[]>1))/sum(as.numeric(apply(snp.bin.mat,1,function(x){sum(x)})[]>0))
			#	apply(ARRAY_collapsed_pbinom_ACGTZ,2,function(x){sum(x)})

	# - 3 - Coll_ZERO_pbinom_matrix	
	for(s.nr in 1:length(group.names)){
		one.gr.pos	<-group.pos[[s.nr]]
		repl.nr		<-sum(as.numeric(group.pos[[s.nr]]))
		temp.arr	<-z.cov.mat[,which(one.gr.pos[]==1),,drop=FALSE]#dim(z.cov.mat); dim(temp.arr)
		sum.temp.arr<-as.numeric((apply(temp.arr,1,function(x){sum(x)}))[]==repl.nr)
		Coll_ZERO_pbinom_matrix[,s.nr,]<-sum.temp.arr
		}

	# - 4 - final correction		
	ARRAY_collapsed_pbinom_ACGTZ[,,5]<-Coll_ZERO_pbinom_matrix
	ARRAY_collapsed_pbinom_ACGTZ[which(ARRAY_collapsed_pbinom_ACGTZ[,,]==0)]<-as.numeric(min.for.p)	
		#aa<-Coll_ZERO_pbinom_matrix; colnames(aa)<-1:ncol(aa); rownames(aa)<-1:nrow(aa); round(aa,digits=3)
		#aa<-ARRAY_collapsed_pbinom_ACGTZ; colnames(aa)<-1:ncol(aa); rownames(aa)<-1:nrow(aa); round(aa,digits=3)
		#summary(aa)
		#round(aa[1356,,],digits=2)
	
	# - 5 - return	
	return(ARRAY_collapsed_pbinom_ACGTZ)
	}

























































###	---------------------------------------------------------------------------------------------------
###	17. COLLAPSE_InDel_MatrixNA__function__May2016		3300		
###	---------------------------------------------------------------------------------------------------
COLLAPSE_InDel_MatrixNA__function__May2016<-function(
	In.dir,
	Res.dir,
	collapse.info,
	matrix.na){


 	# 0. info	
	cat("[ .+. ] - ",date(),"\n")
	cat("[ .+. ] -  elements are:","\n")	
	cat("[ ... ] -  ",matrix.na,"\n")
	cat("[ ... ] -  ",collapse.info,"\n")		

	# 1. New object name
	setwd(In.dir);dir()
	temp_Data_set_name	<-matrix.na
	temp_Data_set_name	<-strsplit(temp_Data_set_name,split="__")[[1]][1]; temp_Data_set_name
	new.name				<-paste(temp_Data_set_name,"__","Collapsed_InDel_MatrixNA_RG.RData",sep=""); 
	new.name


	# 2. load
 	setwd(In.dir)
 	load(matrix.na)
 	load(collapse.info)
 	indel.matrixNA		<-matrixNA_LIST$indel_matrixNA
 	mode(indel.matrixNA)	<-"numeric"

 	
 	# 3. collapse info	
	names(Samples.Cross.Refference.essencials.for.substitution)
	reff.list							<-Samples.Cross.Refference.essencials.for.substitution; 	
	collapse.group.names					<-as.list(names(reff.list)); names(collapse.group.names)<-names(reff.list)
	collapse.group.pos					<-collapse.group.names
	collapse.group.r.nr					<-collapse.group.names
	###	
	for(i in 1:length(reff.list)){collapse.group.r.nr[[i]]<-(reff.list[[i]])[[1]]}
	for(i in 1:length(reff.list)){collapse.group.pos[[i]]<-(reff.list[[i]])[[2]]}

	
	# 4. new matrix	
	Collapsed_InDel_MatrixNA				<-matrix(ncol=length(collapse.group.names),nrow=nrow(indel.matrixNA)); dim(Collapsed_InDel_MatrixNA)
	colnames(Collapsed_InDel_MatrixNA)	<-names(collapse.group.names)
	rownames(Collapsed_InDel_MatrixNA)	<-rownames(indel.matrixNA)
	Collapsed_InDel_MatrixNA[,]			<-0	


	# 5. fill in
	gr.nr<-1;	
	for(gr.nr in 1:length(collapse.group.pos)){
		gr.pos							<-collapse.group.pos[[gr.nr]];gr.pos
		Collapsed_InDel_MatrixNA[,gr.nr]<-as.numeric(as.numeric(apply(indel.matrixNA[,which(gr.pos[]==1),drop=FALSE],1,function(x){sum(x)}))[]>0)
	}
	
		
	# 6. save
	setwd(Res.dir)
	save(Collapsed_InDel_MatrixNA,file=new.name)
	dir()


	# 7. info
	cat("[ ... ] -  ",sum(Collapsed_InDel_MatrixNA),"InDels has beed detected","\n")	
	cat("[ .s. ] -  ","Collapsed_InDel_MatrixNA SAVED","\n")
	cat("[ .s. ] -  ",date(),"\n")											
	cat("\n")						
	}







































###	* ---------------------------------------------------------------------------------------------------
###	 Run.GapJumper.five.methods   -update 2016.08.28	3420
###	* ---------------------------------------------------------------------------------------------------
Run.GapJumper.five.methods<-function(
		GapJumper_Data.dir,
		GapJumper_Intermediate_Files.dir,
		GapJumper_Results.dir,
		GapJumper_RgTable.dir,
		Project.Name	,
		RgTable.Name	,
		Position_Names,
		Positions_to.use,
		Positions_to.exclude,
		z.classifier	,
		Run.GapJumper.method,			
		Run.GapJumper.beta.method,
		Run.fixed.threshold.method,
		Run.site.elimination.method,
		Run.conservative.method,
		GapJumper.detection.limit,
		GapJumper.prior.p,
		FixedThresholds.detection.limit,				
		FixedThresholds.min.allele.perc,						
		FixedThresholds.min.nr.of.samples.with.same.allele,	
		FixedThresholds.max.na.positions.perc,					
		SiteElimination__na.max.nr.per.locus	,					
		SiteElimination__eliminate.pos.with.more.than.na.perc,	
		value.for.zero,				
		Run.process.process.subVCF,
		Run.array.creator_and.add.RG	,
		Run.process.sample.cross.validation,
		Run.process.ACGTZ.array.creator,	
		make.info.arrays,
		find.indels.in.each.samaple,
		make.multiSNV.files,
		put.replicate.names.into.multiSNV,
		Run.GapJumper.alfa.method
		){			
			
###	* ---------------------------------------------------------------------------------------------------
###	GapJumper
###	* ---------------------------------------------------------------------------------------------------	

	# 1. because I had to many lines of code:
		
		# /1. preset.variables
		vcf.data.dir									<-GapJumper_Data.dir
		gapjumper.working.dir						<-GapJumper_Intermediate_Files.dir
		gapjumper.results.dir						<-GapJumper_Results.dir
		rg.table.dir									<-GapJumper_RgTable.dir
		
		# /2. if else loops for omiting finished parts ; fill free to use them
		#	Run.process.process.subVCF					<-"yes"
		#	Run.array.creator_and.add.RG					<-"yes"	
		#	Run.process.sample.cross.validation			<-"yes"	
		#	Run.process.ACGTZ.array.creator				<-"yes"	
		#	make.info.arrays							<-"yes"
		
		# /3. other parts: # best dont touch this !
		run.main.set.only							<-"yes"		
		if(Run.GapJumper.beta.method[]=="y"){
				run.beta.correction					<-"yes"
		}else{	run.beta.correction					<-"no"}
		run.part.with.beta.correction				<-run.beta.correction	# it is other version of GapJumper	
																			# if you get here, best contact me via email.
		if(Run.GapJumper.alfa.method[]=="y"){
		run.par.without.beta.correction	<-"yes"														
		}else{
		run.par.without.beta.correction	<-"no"												
		}									
																											
		# /4. and more:		
		rg.table.name								<-RgTable.Name	
		Project.name									<-Project.Name	

		# /5. some data for final arrays
		#use.as.separator.in.info.arr		<-";"
		# it is inside one fo the functions


	# 2. loading the basic
		
		# /1. GapJumper - outside:
		#	setwd(GapJumper.Module.dir)
		#	source("GapJumper_v100_Module_May2016.R")	
	
		# /2. RG table
		setwd(GapJumper_RgTable.dir); dir()
		RG_sample_table			<-as.matrix(read.table(RgTable.Name	,sep=";",header=TRUE)) # RG_sample_table
								#  it has a header and two columns
								#	1st - .simplifiedVCF files
								#	2nd	- .multiVCF file name
			
		# /3. vcf.type
		#if(grepl("simplifiedVCF",RG_sample_table[1,1])[]==TRUE){vcf.type<-"simplifiedVCF"}else{vcf.type<-"subVCF"}; vcf.type
		setwd(vcf.data.dir); 
		if(sum(grepl("simplified_VCF",dir()))[]>0){vcf.type<-"simplified_VCF"}else{vcf.type<-"subVCF"}; vcf.type


	# 3. General info about the positions		fill free to modify

		# /1. load one file
		setwd(vcf.data.dir); 
		dir.list							<-dir(); 
		vcf.list							<-dir.list[grep(vcf.type,dir.list,fixed=T)]			
		any.vcf							<-read.table(vcf.list[1], sep="\t",as.is=T,strip.white=T,colClasses="character") ;vcf.list[1]
 		any.vcf							<-as.matrix(any.vcf)  				# it has loci, but not headers
		dim(any.vcf)

		# /2. positions.names
		if(Position_Names[1]=="n"){
		Position_names			<-1:nrow(any.vcf	)
		}else{Position_names	<-Position_Names}
		positions.names			<-Position_names
		#Position_names
		
		# /3. positions.to.use		
		if(Positions_to.use[1]=="n"){
		Positions_to.use				<-rep(1,length(positions.names))
		}else{}
		# binary, [0,1], zero will be excluded from the project
		
		# /4. positions.to.exclude
		if(Positions_to.exclude[1]=="n"){
		Positions_to.exclude			<-rep(0,length(positions.names))
		}else{}
		#exclude, alternative set, because sometimes  we have positions that are always exluded		
		# e.g. mitochondrial seq. from the genome 																		
		# in this case 1 will be excluded (binary, [0,1])

		# /5. z.classifier		
		if(z.classifier[1]=="n"){
		z.classifier						<-rep(0,length(positions.names))
		}else{}
		# binary, [0,1], positions with the same classifier with be cross-validated together
		
		RADloci.custom.names		<-as.character(Position_names)			
		RADloci.to.use				<-Positions_to.use	
		RADloci.to.remove			<-Positions_to.exclude	

		# remove what must be removed from z. clasiffier - it was causing many many problems
		remove.from.z									<-rep(0,length(z.classifier))
		remove.from.z[which(Positions_to.use[]>0)]		<-1
		remove.from.z[which(Positions_to.exclude[]>0)]	<-0
		z.classifier										<-z.classifier[which(remove.from.z[]==1)]
					
	# 4. Info:			
		cat("...","\n"); cat("******************************************************************************","\n")
		cat("***","GapJumper v.1.0.0","\n")
		cat("***","Pawel Rosikiewicz","\n")
		cat("***","University of Lausanne, Switzerland","\n")
		cat("***","May 2016","\n")
		cat("******************************************************************************","\n"); cat("...","\n")
		cat("...","prior p       			  :",GapJumper.prior.p,"\n")
		cat("...","dl     		 			  :",GapJumper.detection.limit,"\n")
		cat("...","Rg table		 			  :",rg.table.name,"\n")
		cat("...","Project Name  			  :",Project.name,"\n"); cat("...","\n")			
		cat("...","replicates.nr  			  :",length(table(RG_sample_table[,1])),"\n")
		cat("...","nr of consensus to build 	  :",length(table(RG_sample_table[,2])),"\n")
		cat("...","total nr of replicates	  :",nrow(RG_sample_table),"\n"); cat("...","\n")			
		cat("...","positions nr           	  :",length(Position_names),"\n"); cat("...","\n")		
		cat("...","Run GapJumper method		  :",Run.GapJumper.method,"\n")
		cat("...","Run fixed threshold method :",Run.fixed.threshold.method,"\n")
		cat("...","Run site elimination method:",Run.site.elimination.method	,"\n")
		cat("...","Run conservative method    :",Run.conservative.method	,"\n"); cat("...","\n")			
		cat("******************************************************************************","\n"); cat("...","\n"); cat("..","\n")
			
			
	#   .........................................................................
	###	-> if(1) ------>
	if(Run.process.process.subVCF[]=="yes"){
	#   .........................................................................

	# ---------------------------------------------------------------------------------------------------
	### 5. Run Process subVCF function		
	
	### info:
		cat("\n"); cat(".........................................................................","\n" )
		cat("***","\n" )
		cat("***","step 1. Process subVCF","   >>>  ",date(),"\n" )
		cat("***","        processing subVCF or simplifiedVCF into R format","\n" )
		cat("***","\n" )
		cat("\n"); cat(".........................................................................","\n" )	
	
	### parameters:
		In.directory					<-vcf.data.dir
		Res.directory				<-gapjumper.working.dir	
		subVCF.list					<-RG_sample_table[,1]; subVCF.list	
		#RADloci.custom.names		<-Position_names			
		#RADloci.to.use				<-Positions_to.use	
		#RADloci.to.remove			<-Positions_to.exclude		
		Set.name						<-Project.name
	
	### run:
		process_subVCF_May2016_function(
							In.directory	,
							Res.directory,	
							subVCF.list,
							RADloci.custom.names,
							RADloci.to.use,	
							RADloci.to.remove,
							Set.name	)	

	#   .........................................................................
	#if	(1) ------*
	}else{
		cat("\n"); cat(".........................................................................","\n" )
		cat("***","\n" )
		cat("***","THESE STEPS HAVE NOT BEED DONE:","   >>>  ",date(),"\n" )
		cat("***"," processing subVCF or simplifiedVCF into R format","\n" )
		cat("***","\n" )
		cat(".........................................................................","\n" ); cat("\n")		
	}
	
	
		
	###	-> if(2 and 3) ------>
	if(Run.array.creator_and.add.RG	[]=="yes"){
	#   .........................................................................	
		
	# ---------------------------------------------------------------------------------------------------
	### 6. Creating AGTZ coverage array	
	
	###	info:
		cat("\n"); cat(".........................................................................","\n" )
		cat("***","\n" )
		cat("***","step 2. Creating AGTZ coverage array","   >>>  ",date(),"\n" )
		cat("***","\n" )
		cat(".........................................................................","\n" ); cat("\n")	
	
	### parameters:
		setwd(gapjumper.working.dir)	;dir()
		In.directory				<-gapjumper.working.dir
		Res.directory			<-gapjumper.working.dir
		data_list				<-dir()[which(grepl(paste(Project.name,"__matrixNA_LIST.RData",sep=""),dir())[]==TRUE)]
		data_list
	
	### warnings:
		if(length(data_list)[]>1){
			cat("ERROR","\n" );
			cat("ERROR","      Multiple elements matrixNA_LIST fit to the desription","\n");
			cat("ERROR","      GapJumper can use only one !","\n");  	
			cat("ERROR","      Change the project name so the names of matrices are not confused with each other !","\n");  
			cat("ERROR","  !   Only the first Matrix set has been used	", "\n" );cat("ERROR","\n" )
		}else if(length(data_list)[]==0){
			cat("ERROR","\n" );
			cat("ERROR","      No elements matrixNA_LIST fit to the desription","\n");
			cat("ERROR","      GapJumper can use only one !","\n");  cat("ERROR","\n" )	
		}else{}
	
	### run:
		simple.coverageARRAY.creator_function20150719(
			data_list[1],
			In.directory,
			Res.directory)
		# first it calculates filtered coverage - takes a time 
		# filt cov takes also a cov for indels
		dir()

	# ---------------------------------------------------------------------------------------------------
	### 7. Add RG group to colnames in ACGT array (very fast, 1-10sec)
	###	info:
		cat("\n"); cat(".........................................................................","\n" )
		cat("***","\n" )
		cat("***","step 3.  Add Refference Group (RG) to sample in ACGT array","   >>>  ",date(),"\n" )
		cat("***","samples with the same RG will be crossvalidated and collapsed together in a following steps","\n" )
		cat("***","\n" )
		cat(".........................................................................","\n" ); cat("\n")	
	
	### parameters:
		setwd(gapjumper.working.dir)	;dir()
		In.directory				<-gapjumper.working.dir
		Res.directory			<-gapjumper.working.dir
		data_list				<-dir()[which(grepl(paste(Project.name,"__ARRAY_CovACGT.RData",sep=""),dir())[]==TRUE)]
		data_list
	
	### function:
		New_Add_RG_group_to_Sample_name__function<-
			function(data.list,RG_Group_table,In.directory,Res.directory,rg.separator){
			for(set.nr in 1:length(data_list)){		
			setwd(In.directory); load(data_list[set.nr]);data_list[set.nr]
			set.name		<-strsplit(data_list[set.nr],split=".RData")[[1]][1]	;set.name				
			save.name	<-paste(set.name,"_RG.RData",sep="")
			cat("\n"); cat("***","creating RG groups for: - ",set.nr," - ",save.name, "\n")
			s.names.rg<-paste(RG_Group_table[,2],rg.separator,RG_Group_table[,1],sep="")
			colnames(ARRAY_CovACGT)<-s.names.rg
			setwd(Res.directory); save(ARRAY_CovACGT,file=save.name)}}
		
	### run:
		rg.separator<-"__"
		New_Add_RG_group_to_Sample_name__function(data.list,RG_sample_table,In.directory,Res.directory,rg.separator)
		dir()
		#	test
		#	load("GenSt_Run1_TestSet__ARRAY_CovACGT_RG.RData")		
		#	colnames(ARRAY_CovACGT)

	#   .........................................................................
	#if	(2 and 3) ------*
	}else{
		cat("\n"); cat(".........................................................................","\n" )
		cat("***","\n" )
		cat("***","THESE STEPS HAVE NOT BEED DONE:","   >>>  ",date(),"\n" )
		cat("***","Creating AGTZ coverage array and adding RG group, ie sample name to each replicate name","\n" )
		cat("***","\n" )
		cat(".........................................................................","\n" ); cat("\n")		
	}	
	###	-> if(4 and 5) ------>
	if(Run.process.sample.cross.validation[]=="yes"){
	#   .........................................................................		

	# ---------------------------------------------------------------------------------------------------
	### 8. Cross-Validation - Additional elements (ARRAY_CovACGT_RG_additional_elements_for_CV):
	###	info:
		cat("\n"); cat(".........................................................................","\n" )
		cat("***","\n" )
		cat("***","step 4.  Cross-validation","   >>>  ",date(),"\n" )
		cat("***","a)       Additional elements","\n" )
		cat("***","         samples/replicates with the same RG will be cross-validated with each other","\n" )
		cat("***","         only variable positions will be cross-validated ","\n" )
		cat("***","         the stable positions have the weight equal to 1","\n" )
		cat("***","\n" ); 
		cat(".........................................................................","\n" ); cat("\n")	
	
	### parameters:
		setwd(gapjumper.working.dir)	;dir()
		In.directory				<-gapjumper.working.dir
		Res.directory			<-gapjumper.working.dir
		Res.two.dir				<-gapjumper.working.dir
		data_list				<-dir()[which(grepl(paste(Project.name,"__ARRAY_CovACGT_RG.RData",sep=""),dir())[]==TRUE)]
		data_list
		z.classifier				<-z.classifier	
	###	run:
		Additional_elements_for_Cross_Validation_function__June2016(In.directory,Res.directory,data_list[1],z.classifier)

	# ---------------------------------------------------------------------------------------------------
	### 9. Cross-validation itself
	###	info:
		cat("\n"); cat(".........................................................................","\n" )
		cat("***","\n" )
		cat("***","step 4.  Cross-validation","   >>>  ",date(),"\n" )
		cat("***","b)       The CV process","\n" )
		cat("***","\n" )
		cat(".........................................................................","\n" ); cat("\n")		
		
	###	data_list:	
		setwd(gapjumper.working.dir)	;dir()
		In.directory				<-gapjumper.working.dir
		Res.directory			<-gapjumper.working.dir
		Res.two.dir				<-gapjumper.working.dir
		###
		data.list				<-dir()[which(grepl(paste(Project.name,"__ARRAY_CovACGT_RG.RData",sep=""),dir())[]==TRUE)]
		data.list		
		###		
		validation.list			<-dir()[which(grepl(paste(Project.name,"__ARRAY_CovACGT_RG_additional_elements_for_CV.RData",sep=""),dir())[]==TRUE)]
		validation.list
								# this is a data.list done for one sample set				
	
	###	function.parameters:
		additional.name			<-""#best start with _			
		In.directory				<-Res.two.dir
		In.validation.directory	<-Res.two.dir
		Res.directory			<-Res.two.dir
		min.coverage.accepted	<-1;	#GapJumper.detection.limit, we use 1 to get better CV values, but we use filter later
		###
		fast.protocol			<-"yes"
		#fast.protocol			<-"no"
		# 		yes	-			with this we wont get cross validation to all pos and rep and not repeated
		# 						there are only 3 remaining categories i.e. variable pos with repeats 
		#						and not repeats and finally stable where I should get only 0 and 1	
		# 		no -			with this protocol I do all categories - it is alsmost 3x longer
	
	###
		#unique.count.protocol.nr	<-1	
		unique.count.protocol.nr		<-2 
		#		1 - 			smaler values, done with all known positions 
		# 		2 -			higher values, done only with each snp data
	
	###	run
		fast.set<-1:length(data.list)# to choose some in case I want
		#fast.set<-2
		gamma_cross_validation_procedure__function__May2016(
			In.directory	,
			In.validation.directory,
			Res.directory,
			data.list[1],
			validation.list[1],
			min.coverage.accepted,
			fast.protocol,
			unique.count.protocol.nr,
			additional.name)

	#   .........................................................................
	#if	(4 and 5) ------*
	}else{
		cat("\n"); cat(".........................................................................","\n" )
		cat("***","\n" )
		cat("***","THESE STEPS HAVE NOT BEED DONE:","   >>>  ",date(),"\n" )
		cat("***"," Cross-validation","\n" )
		cat("***","\n" )
		cat(".........................................................................","\n" ); cat("\n")		
	}	
	###	-> if(5) ------>
	if(Run.process.ACGTZ.array.creator[]=="yes"){
	#   .........................................................................

	# ---------------------------------------------------------------------------------------------------
	### 10. ACGTz array creator	
	###	   can be a long step 
	###	   5 - min for large sample set, 100 samples with 100.000 snp's
	###	info:
		cat("\n"); cat(".........................................................................","\n" )
		cat("***","\n" )
		cat("***","step 5. ACGTz array creator","   >>>  ",date(),"\n" )
		cat("***","\n" )
		cat(".........................................................................","\n" ); cat("\n")		
		
	###	data list:	
		setwd(gapjumper.working.dir)	;dir()
		In.directory				<-gapjumper.working.dir
		Res.directory			<-gapjumper.working.dir
		Res.two.dir				<-gapjumper.working.dir
		Res.dir					<-gapjumper.working.dir
		###
		data.list				<-dir()[which(grepl(paste(Project.name,"__ARRAY_CovACGT_RG.RData",sep=""),dir())[]==TRUE)]
		data.list		
		###		
		validation.list			<-dir()[which(grepl(paste(Project.name,"__ARRAY_CovACGT_RG_additional_elements_for_CV.RData",sep=""),dir())[]==TRUE)]
		validation.list
		###
		cv.essentials.list		<-dir()[grepl(paste(Project.name,"__GammaValidationProcedure_Samples.Cross.Refference.essencials.for.substitution.RData",sep=""),dir())]
		cv.essentials.list		

	### function parameters
		array.data.list								<-data.list[1]
		cross.validation.essentials.for.array.list	<-cv.essentials.list[1]
		
		In.validation.directory	<-gapjumper.working.dir
		In.Valid.essentials.dir	<-gapjumper.working.dir
		Res.directory			<-gapjumper.working.dir
		Res.dir					<-gapjumper.working.dir
		Res.two.dir				<-gapjumper.working.dir
		# ###
		additional.name			<-""#best start with _
		min.coverage.accepted	<-GapJumper.detection.limit						
		parameters.for.pbinom	<-GapJumper.prior.p
		Noise.value.for.zero.cov	<-value.for.zero
		coverage.for.pbinom		<-"filtered"
		#coverage.for.pbinom	<-"general"; gives lower Pr() to nucleotides at positions where was a lot of bad quality reads
		# ###
		# this must a very small value 
		# it is used when there is no data within a position
		# and it ensure that we have no later problems with the script (it shouldnt be but it is much better)
		# the value must be much smaller than noise level, you wish to use later
		# eg. noise = 0.01 and this site would be = 0.001 or even 0.0001 

	### in case you wish to make a subset
		array.data.list								<-data.list[1]
		cross.validation.essentials.for.array.list	<-cv.essentials.list
		fast.set										<-1:length(array.data.list)


		New.ACGTZ.array.creator__with.gamma.validation__function__May2016(
			array.data.list, 
			cross.validation.essentials.for.array.list, 
			parameters.for.pbinom,
			Noise.value.for.zero.cov,
			In.directory	,
			In.Valid.essentials.dir,
			Res.directory,
			additional.name,
			coverage.for.pbinom)

		setwd(Res.two.dir); dir()


	#   .........................................................................
	#if	(4 and 5) ------*
	}else{
		cat("\n"); cat(".........................................................................","\n" )
		cat("***","\n" )
		cat("***","THESE STEPS HAVE NOT BEED DONE:","   >>>  ",date(),"\n" )
		cat("***","step 5","Run.process.ACGTZ.array.creator, format better for GapJumper, 3D matrix","\n" )
		cat("***","\n" )
		cat(".........................................................................","\n" ); cat("\n")		
	}	

	# ---------------------------------------------------------------------------------------------------
	####	####################################################################			
	####	11.a.		COLLAPSE ARRAYS	- very long step 5 - min for large sample set, 100 samples with 100.000 snp's
	####			GapJumper Method
	####	####################################################################

	if(Run.GapJumper.method[]=="y"){
	### info
		cat("\n"); cat(".........................................................................","\n" )
		cat("***","\n" )
		cat("***","step 7. COLLAPSE ARRAYS - GapJumper Method","\n" )
		cat("***","\n" )
		cat(".........................................................................","\n" ); cat("\n")		
	
	###	collaspe info list
		setwd(gapjumper.working.dir)	;dir()
		In.directory				<-gapjumper.working.dir
		Res.directory			<-gapjumper.working.dir
		Res.two.dir				<-gapjumper.working.dir
		Res.dir					<-gapjumper.working.dir		
		Valid.dir				<-gapjumper.working.dir
		acgtz.ar.dir				<-gapjumper.working.dir
		Res.directory			<-gapjumper.working.dir
		###
	
	### 5x arrays (4 in use)													
		setwd(gapjumper.working.dir)	;dir()
		file_list				<-dir()
		cov.ACGTZ				<-file_list[grep(paste(Project.name,"__ARRAY_CovACGTZ_RG.RData",sep=""),file_list,fixed=T)]
		pbinom.ACGTZ				<-file_list[grep(paste(Project.name,"__ARRAY_pbinom_ACGTZ_LIST_RG.RData",sep=""),file_list,fixed=T)]
		valid.ACGTZ				<-file_list[grep(paste(Project.name,"__ARRAY_ValidationACGTZ_RG.RData",sep=""),file_list,fixed=T)]
		gen.cov.ACGTZ			<-file_list[grep(paste(Project.name,"__ARRAY_GenCovACGTZ_RG.RData",sep=""),file_list,fixed=T)]
		filt.cov					<-file_list[grep(paste(Project.name,"__ARRAY_FiltCovACGTZ_RG.RData",sep=""),file_list,fixed=T)]
		binary.ACGTZ				<-file_list[grep(paste(Project.name,"__ARRAY_BinaryACGTZ_RG.RData",sep=""),file_list,fixed=T)]	
		collapse.info			<-dir()[grepl(paste(Project.name,"__GammaValidationProcedure_Samples.Cross.Refference.essencials.for.substitution.RData",sep=""),dir())]
		collapse.info		
		
	### TEST:						
				# ---------------------
				collapse.info
				# ---------------------
				pbinom.ACGTZ
				valid.ACGTZ
				cov.ACGTZ
				gen.cov.ACGTZ
				# ---------------------
		setwd(gapjumper.working.dir); #dir()		
		load(cov.ACGTZ[1])	
		mode(ARRAY_CovACGTZ)<-"numeric"				
				
	### TEST:					
		additional.name					<-""		
		pbinom.to.take					<-1	#0.1# old version of softwear could use multiple p values, 
		min.coverage.accepted			<-GapJumper.detection.limit	
		when.weight.is.zero.put			<-value.for.zero
		use.full.cov.spectrum.at.na		<-"no"	#can be no i.e. all sites are with wieght equal to zero coverage
											#yes - slightly more cobnservative method"
	### elements.to.run
		run.part.without.beta.correction<-run.par.without.beta.correction
		run.part.with.beta.correction	<-run.beta.correction
		prepare.info.arrays				<-make.info.arrays
		Name_add.method					<-"_GapJumper_method"		

			
	### run
		fast.set<-1:length(collapse.info)
		#fast.set<-2
		COLLAPSE_ACGTZ_ARRAYS__GAMMA__function_May2016(
			#	collapse.info
				Valid.dir,
				collapse.info[fast.set],	
			#	arrays and results:
				acgtz.ar.dir,
				Res.directory,
			#	array lists:		
				pbinom.ACGTZ[fast.set],
				valid.ACGTZ[fast.set],
				cov.ACGTZ[fast.set],
				gen.cov.ACGTZ[fast.set],
			#	parameters:
				additional.name,
				pbinom.to.take,
				min.coverage.accepted,
				when.weight.is.zero.put,
				use.full.cov.spectrum.at.na,
			#	set.name
				run.part.without.beta.correction,
				run.part.with.beta.correction,
				prepare.info.arrays,
			#
				Name_add.method)				
	}else{
	### info
		cat("\n"); cat(".........................................................................","\n" )
		cat("***","\n" )
		cat("***","step 7. COLLAPSE ARRAYS - GapJumper Method - NOT RUN","\n" )
		cat("***","\n" )
		cat(".........................................................................","\n" ); cat("\n")	
	}


	# ---------------------------------------------------------------------------------------------------
	####	####################################################################			
	####	6.b. 	COLLAPSE ARRAYS		- very long step 5 - min for large sample set, 100 samples with 100.000 snp's
	####			CONSERVATIVE METHOD
	####	####################################################################
	if(Run.conservative.method[]=="y"){
	### info
		cat("\n" ); cat(".........................................................................","\n" )
		cat("***","\n" )
		cat("***","step 7.B CONSERVATIVE METHOD","\n" )
		cat("***","\n" )
		cat(".........................................................................","\n" ); cat("\n")		
	
	###	collaspe info list
		setwd(gapjumper.working.dir)	;dir()
		In.directory				<-gapjumper.working.dir
		Res.directory			<-gapjumper.working.dir
		Res.two.dir				<-gapjumper.working.dir
		Res.dir					<-gapjumper.working.dir		
		Valid.dir				<-gapjumper.working.dir
		acgtz.ar.dir				<-gapjumper.working.dir
		Res.directory			<-gapjumper.working.dir
		###
	
	### 5x arrays (4 in use)													
		setwd(gapjumper.working.dir)	;dir()
		file_list				<-dir()
		cov.ACGTZ				<-file_list[grep(paste(Project.name,"__ARRAY_CovACGTZ_RG.RData",sep=""),file_list,fixed=T)]
		pbinom.ACGTZ				<-file_list[grep(paste(Project.name,"__ARRAY_pbinom_ACGTZ_LIST_RG.RData",sep=""),file_list,fixed=T)]
		valid.ACGTZ				<-file_list[grep(paste(Project.name,"__ARRAY_ValidationACGTZ_RG.RData",sep=""),file_list,fixed=T)]
		gen.cov.ACGTZ			<-file_list[grep(paste(Project.name,"__ARRAY_GenCovACGTZ_RG.RData",sep=""),file_list,fixed=T)]
		filt.cov					<-file_list[grep(paste(Project.name,"__ARRAY_FiltCovACGTZ_RG.RData",sep=""),file_list,fixed=T)]
		binary.ACGTZ				<-file_list[grep(paste(Project.name,"__ARRAY_BinaryACGTZ_RG.RData",sep=""),file_list,fixed=T)]	
		collapse.info			<-dir()[grepl(paste(Project.name,"__GammaValidationProcedure_Samples.Cross.Refference.essencials.for.substitution.RData",sep=""),dir())]
		collapse.info		
		### TEST:						
				# ---------------------
				collapse.info
				# ---------------------
				pbinom.ACGTZ
				valid.ACGTZ
				cov.ACGTZ
				gen.cov.ACGTZ
				# ---------------------

	### TEST:					
		additional.name					<-""		
		pbinom.to.take					<-1	#0.1# old version of softwear could use multiple p values, 
		min.coverage.accepted			<-GapJumper.detection.limit	
		when.weight.is.zero.put			<-0
		use.full.cov.spectrum.at.na		<-"no"	#can be no i.e. all sites are with wieght equal to zero coverage
											#yes - slightly more cobnservative method"
	### elements.to.run
		prepare.info.arrays				<-make.info.arrays			
		Name_add.method					<-"Conservative_method"
	
	### run
		fast.set<-1:length(collapse.info)	
		#fast.set<-2
		Conservative_method__COLLAPSE_ACGTZ_ARRAYS__GAMMA__function__May2016(
		#	collapse.info
		Valid.dir,
		collapse.info[fast.set],	
		#	arrays and results:
		acgtz.ar.dir,
		Res.directory,
		#	array lists:		
		pbinom.ACGTZ[fast.set],
		valid.ACGTZ[fast.set],
		cov.ACGTZ[fast.set],
		gen.cov.ACGTZ[fast.set],
		#	for threshold-based methods	
		filt.cov[fast.set],
		binary.ACGTZ	[fast.set],
		Name_add.method,	
		#	parameters:
		additional.name,
		pbinom.to.take,
		min.coverage.accepted,
		when.weight.is.zero.put,
		use.full.cov.spectrum.at.na,
		#	set.name
		run.part.without.beta.correction,
		run.part.with.beta.correction,
		prepare.info.arrays)				
	}else{
	### info
		cat("\n" ); cat(".........................................................................","\n" )
		cat("***","\n" )
		cat("***","step 7.B CONSERVATIVE METHOD - NOT RUN","\n" )
		cat("***","\n" )
		cat(".........................................................................","\n" ); cat("\n")	
	}			


	# ---------------------------------------------------------------------------------------------------
	####	####################################################################			
	####	6.c. 	COLLAPSE ARRAYS		- very long step 5 - min for large sample set, 100 samples with 100.000 snp's
	####			SITE ELIMIATION
	####	####################################################################
	if(Run.site.elimination.method[]=="y"){
	### info
		cat("\n" ); cat(".........................................................................","\n" )
		cat("***","\n" )
		cat("***","step 7.C SITE ELIMIATION METHOD","\n" )
		cat("***","\n" )
		cat(".........................................................................","\n" ); cat("\n")		
	###	collaspe info list
		setwd(gapjumper.working.dir)	;dir()
		In.directory				<-gapjumper.working.dir
		Res.directory			<-gapjumper.working.dir
		Res.two.dir				<-gapjumper.working.dir
		Res.dir					<-gapjumper.working.dir		
		Valid.dir				<-gapjumper.working.dir
		acgtz.ar.dir				<-gapjumper.working.dir
		Res.directory			<-gapjumper.working.dir
		###
	
	### 5x arrays (4 in use)													
		setwd(gapjumper.working.dir)	;dir()
		file_list				<-dir()
		cov.ACGTZ				<-file_list[grep(paste(Project.name,"__ARRAY_CovACGTZ_RG.RData",sep=""),file_list,fixed=T)]
		pbinom.ACGTZ				<-file_list[grep(paste(Project.name,"__ARRAY_pbinom_ACGTZ_LIST_RG.RData",sep=""),file_list,fixed=T)]
		valid.ACGTZ				<-file_list[grep(paste(Project.name,"__ARRAY_ValidationACGTZ_RG.RData",sep=""),file_list,fixed=T)]
		gen.cov.ACGTZ			<-file_list[grep(paste(Project.name,"__ARRAY_GenCovACGTZ_RG.RData",sep=""),file_list,fixed=T)]
		filt.cov					<-file_list[grep(paste(Project.name,"__ARRAY_FiltCovACGTZ_RG.RData",sep=""),file_list,fixed=T)]
		binary.ACGTZ				<-file_list[grep(paste(Project.name,"__ARRAY_BinaryACGTZ_RG.RData",sep=""),file_list,fixed=T)]	
		collapse.info			<-dir()[grepl(paste(Project.name,"__GammaValidationProcedure_Samples.Cross.Refference.essencials.for.substitution.RData",sep=""),dir())]
		collapse.info		
		### TEST:						
				# ---------------------
				collapse.info
				# ---------------------
				pbinom.ACGTZ
				valid.ACGTZ
				cov.ACGTZ
				gen.cov.ACGTZ
				# ---------------------
		setwd(gapjumper.working.dir); #dir()		
		load(cov.ACGTZ[1])	
		mode(ARRAY_CovACGTZ)<-"numeric"				
		ARRAY_CovACGTZ.bis	<-ARRAY_CovACGTZ		

	### TEST:					
		additional.name					<-""		
		pbinom.to.take					<-1	#0.1# old version of softwear could use multiple p values, 
		min.coverage.accepted			<-GapJumper.detection.limit	
		when.weight.is.zero.put			<-0
		use.full.cov.spectrum.at.na		<-"no"	#can be no i.e. all sites are with wieght equal to zero coverage
											#yes - slightly more cobnservative method"
	### elements.to.run
		prepare.info.arrays							<-make.info.arrays	
		Name_add.method								<-"SiteElimination_method"
		na.max.nr.per.locus							<-SiteElimination__na.max.nr.per.locus
		eliminate.loci.with.more.than.na.perc		<-SiteElimination__eliminate.pos.with.more.than.na.perc

	### run
	fast.set<-1:length(collapse.info)
	#fast.set<-2
	Site_elimiantion_method__COLLAPSE_ACGTZ_ARRAYS__GAMMA__function__May2016(
		#	collapse.info
		Valid.dir,
		collapse.info[fast.set],	
		#	arrays and results:
		acgtz.ar.dir,
		Res.directory,
		#	array lists:		
		pbinom.ACGTZ[fast.set],
		valid.ACGTZ[fast.set],
		cov.ACGTZ[fast.set],
		gen.cov.ACGTZ[fast.set],
		#	for threshold-based methods	
		filt.cov[fast.set],
		binary.ACGTZ[fast.set]	,
		Name_add.method,	
		#	parameters:
		additional.name,
		pbinom.to.take,
		min.coverage.accepted,
		when.weight.is.zero.put,
		use.full.cov.spectrum.at.na,
		#	set.name
		run.part.without.beta.correction,
		run.part.with.beta.correction,
		prepare.info.arrays,
		na.max.nr.per.locus,
		eliminate.loci.with.more.than.na.perc,
		ARRAY_CovACGTZ.bis)
	}else{
	### info
		cat("\n" ); cat(".........................................................................","\n" )
		cat("***","\n" )
		cat("***","step 7.C SITE ELIMIATION METHOD - NOT RUN","\n" )
		cat("***","\n" )
		cat(".........................................................................","\n" ); cat("\n")	
	}




	# ---------------------------------------------------------------------------------------------------
	####	####################################################################			
	####	6.d. 	COLLAPSE ARRAYS		- very long step 5 - min for large sample set, 100 samples with 100.000 snp's
	####			FIXED THRESHOLDS METHOD
	####	####################################################################
	if(Run.fixed.threshold.method[]=="y"){
	### info
		cat("\n" ); cat(".........................................................................","\n" )
		cat("***","\n" )
		cat("***","step 7.D FIXED THRESHOLDS METHOD","\n" )
		cat("***","\n" ) 
		cat(".........................................................................","\n" ); cat("\n")		
	###	collaspe info list
		setwd(gapjumper.working.dir)	;dir()
		In.directory				<-gapjumper.working.dir
		Res.directory			<-gapjumper.working.dir
		Res.two.dir				<-gapjumper.working.dir
		Res.dir					<-gapjumper.working.dir		
		Valid.dir				<-gapjumper.working.dir
		acgtz.ar.dir				<-gapjumper.working.dir
		Res.directory			<-gapjumper.working.dir
		###
	
	### 5x arrays (4 in use)													
		setwd(gapjumper.working.dir)	;dir()
		file_list				<-dir()
		cov.ACGTZ				<-file_list[grep(paste(Project.name,"__ARRAY_CovACGTZ_RG.RData",sep=""),file_list,fixed=T)]
		pbinom.ACGTZ				<-file_list[grep(paste(Project.name,"__ARRAY_pbinom_ACGTZ_LIST_RG.RData",sep=""),file_list,fixed=T)]
		valid.ACGTZ				<-file_list[grep(paste(Project.name,"__ARRAY_ValidationACGTZ_RG.RData",sep=""),file_list,fixed=T)]
		gen.cov.ACGTZ			<-file_list[grep(paste(Project.name,"__ARRAY_GenCovACGTZ_RG.RData",sep=""),file_list,fixed=T)]
		filt.cov					<-file_list[grep(paste(Project.name,"__ARRAY_FiltCovACGTZ_RG.RData",sep=""),file_list,fixed=T)]
		binary.ACGTZ				<-file_list[grep(paste(Project.name,"__ARRAY_BinaryACGTZ_RG.RData",sep=""),file_list,fixed=T)]	
		collapse.info			<-dir()[grepl(paste(Project.name,"__GammaValidationProcedure_Samples.Cross.Refference.essencials.for.substitution.RData",sep=""),dir())]
		collapse.info		
		### TEST:						
				# ---------------------
				collapse.info
				# ---------------------
				pbinom.ACGTZ
				valid.ACGTZ
				cov.ACGTZ
				gen.cov.ACGTZ
				# ---------------------
		setwd(gapjumper.working.dir); #dir()		
		load(cov.ACGTZ[1])	
		mode(ARRAY_CovACGTZ)<-"numeric"		
		ARRAY_CovACGTZ.bis<-ARRAY_CovACGTZ

	### TEST:					
		additional.name					<-""		
		pbinom.to.take					<-1	#0.1# old version of softwear could use multiple p values, 
		min.coverage.accepted			<-GapJumper.detection.limit	
		when.weight.is.zero.put			<-0
		use.full.cov.spectrum.at.na		<-"no"	#can be no i.e. all sites are with wieght equal to zero coverage
											#yes - slightly more cobnservative method"
	### elements.to.run
		prepare.info.arrays								<-make.info.arrays	
		Name_add.method									<-"FixedThresholds_method"
		eliminate.loci.with.more.than.na.perc			<-FixedThresholds.max.na.positions.perc
		###
		eliminate.allels.with.less.than.reads.perc		<-FixedThresholds.min.allele.perc	
		minimum.nr.of.replicates.with.same.allele		<-FixedThresholds.min.nr.of.samples.with.same.allele

	### run	
	fast.set<-1:length(collapse.info)
	#fast.set<-2
	Fixed_Thresholds_method__COLLAPSE_ACGTZ_ARRAYS__GAMMA__function__May2016(
		#	collapse.info
			Valid.dir,
			collapse.info[fast.set],	
		#	arrays and results:
			acgtz.ar.dir,
			Res.directory,
		#	array lists:		
			pbinom.ACGTZ[fast.set],
			valid.ACGTZ[fast.set],
			cov.ACGTZ[fast.set],
			gen.cov.ACGTZ[fast.set],
		#	for threshold-based methods	
			filt.cov[fast.set],
			binary.ACGTZ	[fast.set],
			Name_add.method,	
		#	parameters:
			additional.name,
			pbinom.to.take,
			min.coverage.accepted,
			when.weight.is.zero.put,
			use.full.cov.spectrum.at.na,
		#	set.name
			prepare.info.arrays,
			eliminate.loci.with.more.than.na.perc,
			eliminate.allels.with.less.than.reads.perc,
			minimum.nr.of.replicates.with.same.allele,
			ARRAY_CovACGTZ.bis)
	}else{
	### info
		cat("\n" ); cat(".........................................................................","\n" )
		cat("***","\n" )
		cat("***","step 7.D FIXED THRESHOLDS METHOD - NOT DONE","\n" )
		cat("***","\n" )
		cat(".........................................................................","\n" ); cat("\n")
	}




	# ---------------------------------------------------------------------------------------------------
	####	####################################################################			
	####	7.	 	COLLAPSE matrix		
	####			with InDel presence/absence
	####	####################################################################
	
	
	if(find.indels.in.each.samaple[]=="yes"){	
	### info
		cat("\n" ); cat(".........................................................................","\n" )
		cat("***","\n" )
		cat("***","step 7 COLLAPSE matrix with InDel presence/absence","\n" )
		cat("***","\n" ) 
		cat(".........................................................................","\n" ); cat("\n")		
	###	collaspe info list
		setwd(gapjumper.working.dir)	;dir()
		In.dir				<-gapjumper.working.dir
		Res.dir				<-gapjumper.working.dir
	
	### 5x arrays (4 in use)													
		setwd(gapjumper.working.dir)	;#	dir()
		file_list				<-dir()
		matrix.na				<-file_list[grep(paste(Project.name,"__matrixNA_LIST.RData",sep=""),file_list,fixed=T)]
		collapse.info			<-dir()[grepl(paste(Project.name,"__GammaValidationProcedure_Samples.Cross.Refference.essencials.for.substitution.RData",sep=""),dir())]
		### TEST:						
				# ---------------------
				collapse.info
				matrix.na	
				# ---------------------

	### Run
		COLLAPSE_InDel_MatrixNA__function__May2016(
			In.dir,
			Res.dir,
			collapse.info[1],
			matrix.na[1])
	}else{
		cat("\n" ); cat(".........................................................................","\n" )
		cat("***","\n" )
		cat("***","step 7 COLLAPSE matrix with InDel presence/absence - NOT DONE","\n" )
		cat("***","\n" )
		cat(".........................................................................","\n" ); cat("\n")
	}
	

	####	####################################################################			
	####	8.	 	.multiSNV file format		
	####	####################################################################

	if(make.multiSNV.files[]=="yes"){	
		# /1. info
		cat("\n" ); cat(".........................................................................","\n" )
		cat("***","\n" )
		cat("***","step 8 preparing multiSNV","\n" )
		cat("***","Caution, this is the step where AllRecorded_Nucleotides.multiSNV files are created","\n" )
		cat("***","as a bonus to any other method","\n" )
		cat("***","\n" ) 
		cat(".........................................................................","\n" ); cat("\n")		
	
		# /2. DIR
		In.dir					<-GapJumper_Data.dir	
		In.two.dir				<-GapJumper_Intermediate_Files.dir
		Res.dir					<-GapJumper_Results.dir
			
		# /3. lists:
		setwd(In.dir)		;# 	dir()
		s.vcf.list				<-dir()[which(grepl(vcf.type,dir())[]==TRUE)]				
		
		setwd(In.two.dir)	;#	dir()
		file.list				<-dir()[which(grepl(Project.name,dir())[]==TRUE)]		
		info.arr.list			<-file.list[which(grepl("ARRAY_collapsed_info_ACGTZ",file.list)[]==TRUE)]; info.arr.list	
		indel.mat.list			<-file.list[which(grepl("Collapsed_InDel_MatrixNA_RG",file.list)[]==TRUE)]; indel.mat.list	
		
		# /4. Load.indel.mat 
		setwd(In.two.dir)
		load(indel.mat.list[1])
		if(length(indel.mat.list)[]!=1){
			cat("[ ERROR ] -  ","Project name is not unique","\n")
			cat("[ ERROR ] -  ","i.e. it can find more than one project","\n")
			cat("[ ERROR ] -  ","please change it or make it longer","\n")
			}else{}

		# /5. loop
		set.nr<-1 # for each method and project
		for(set.nr in 1:length(info.arr.list)){
		 	setwd(In.two.dir)
		 	load(info.arr.list[set.nr])
		 	
		 	# 0. info				
			cat("\n")
			cat("[ .",set.nr,". ] - ","    - ","\n")
			cat("[ .",set.nr,". ] - ",date(),"\n")
			cat("[ .",set.nr,". ] - ","elements are:","\n")	
			cat("[ .",set.nr,". ] - ",info.arr.list[set.nr],"\n")
			cat("[ .",set.nr,". ] - ",indel.mat.list[1],"\n")	
			cat("[ .",set.nr,". ] - ","it has",length(dimnames(ARRAY_collapsed_info_ACGTZ)[[2]]),".mutilSNV files to create with, each with one consensus","\n")	
			cat("\n")
					
			# 1. One file for each consensus
			s.nr<-1; # for each consensus that was build with one method and tstroed together
			for(s.nr in 1:dim(ARRAY_collapsed_info_ACGTZ)[2]){	 								
				#array
				one.s.arr				<-ARRAY_collapsed_info_ACGTZ[,s.nr,,drop=TRUE]; colnames(one.s.arr)
				pos.nrs					<-as.numeric(rownames(ARRAY_collapsed_info_ACGTZ)); pos.nrs[1:10]	
				first.repl.name			<-strsplit(one.s.arr[1,11],split=";")[[1]][1]
				sample.name				<-strsplit(colnames(ARRAY_collapsed_info_ACGTZ)[s.nr],split="__")[[1]][1]
				cat("[ .",set.nr,". ] - [ .",s.nr,". ] - ",sample.name,"\n")	
				colnames(one.s.arr)
				
				#indel position	
				one.s.indel				<-Collapsed_InDel_MatrixNA[,which(grepl(colnames(ARRAY_collapsed_info_ACGTZ)[s.nr],colnames(Collapsed_InDel_MatrixNA))[]==TRUE)]
				one.s.indel				<-as.numeric(one.s.indel)
				one.s.indel[which(one.s.indel[]==1)]<-"y"
				one.s.indel[which(one.s.indel[]==0)]<-"n"
				
					
				# s.vcf as source of postion information	
				setwd(In.dir)
				first.repl.mat			<-read.table(s.vcf.list[which(grepl(first.repl.name,s.vcf.list)[]==TRUE)], sep="\t",as.is=T,strip.white=T,colClasses="character")
				first.repl.mat			<-first.repl.mat[pos.nrs,]
				dim(first.repl.mat)
				head(first.repl.mat[,c(2:9)])
				
				
				Alleles					<-first.repl.mat[,1]
				Alleles[]				<-"na"
				
				
				# glue it:		value.for.zero
				multiSNV	<-cbind(	rep(sample.name,nrow(first.repl.mat)),
								first.repl.mat[,c(2:9)],
								one.s.arr[,c(1:5)],
								Alleles,
								one.s.arr[,c(6:8)],	
								one.s.indel,
								one.s.arr[,c(10)])	
				
				# new colnames:
				colnames(multiSNV)[1]		<-"Consensus.name"	
				colnames(multiSNV)[2:9]		<-c(colnames(first.repl.mat)[2:5],"x1","x2","x3","Ref")				
				colnames(multiSNV)[10:14]	<-c("M(A)","M(C)","M(G)","M(T)","na")				
				colnames(multiSNV)[15:20]	<-c("snp.presence","snp.presence.in.each.replicate",
				"snp.freq.in.each.replicate","position.coverage.in.each.replicate","InDel","Nr.of.pooled.replicates")
								
				# the method:
				if(grepl("Conservative_method",info.arr.list[set.nr])[]==TRUE){method<-"Conservative_method"; use.pcr5.method<-"no"
				}else if(grepl("SiteElimination_method",info.arr.list[set.nr])[]==TRUE){method<-"SiteElimination_method"; use.pcr5.method<-"no"	
				}else if(grepl("FixedThresholds_method",info.arr.list[set.nr])[]==TRUE){method<-"FixedThresholds_method"; use.pcr5.method<-"no"
				}else if(grepl("GapJumper_method_RG",info.arr.list[set.nr])[]==TRUE){
					use.pcr5.method<-"yes"				
					if(grepl("_NoCorr_",info.arr.list[set.nr])[]==TRUE){method<-"Probabilistic_method"	
					}else{method<-"Probabilistic_beta_method"}		
				};	method		


				# remove zero values
				multiSNV					<-as.matrix(multiSNV)
				ttt						<-as.matrix(multiSNV)[,c(10:13)]
				mode(ttt)				<-"numeric"
				ttt[which(ttt[,]==value.for.zero)]<-0
				multiSNV[,c(10:13)]		<-ttt
				
				# I want na column to be binary; it is specifically for probabilistic methods				
				na.col.NR			<-14
				na.col						<-as.numeric(multiSNV[,c(na.col.NR)])		
				na.col[which(na.col[]!=1)]	<-0	
				multiSNV[,c(na.col.NR)]		<-na.col	
								
				# now make allele presence in 		
				snps.n<-c("A","C","G","T")
				snp.pres<-apply(multiSNV[,c(10:13)],1,function(x){paste(snps.n[which(x[]>0)],collapse="/")})
				multiSNV[,15]<-snp.pres
				
				# postions with na for other methods must be updated
				pos.with.snps				<-as.matrix(multiSNV[,c(10:13)])
				mode(pos.with.snps)			<-"numeric"
				multiSNV[,c(14)]			<-as.numeric(apply(pos.with.snps,1,function(x){sum(x)[]>0})[]==FALSE)
						
				# put "na"
				multiSNV[which(multiSNV[,15]==""),15]<-"na"

				# column for genotype quality (Pr(genotype))									
				empty.column										<-rep(1,nrow(multiSNV))	
				empty.column[which(as.numeric(multiSNV[,14])==1)]	<-0
				multiSNV[,14]										<-empty.column							
				colnames(multiSNV)[14]								<-"M(genotype)"				
				multiSNV											<-as.matrix(multiSNV)				
								
				# caluclate genotype quality for positions (only for probabilistic methods)
				
				if(use.pcr5.method[]=="yes"){	
					pos.to.use											<-which(empty.column[]>0)
					temp.mat											<-multiSNV[pos.to.use,c(10:13)]
					mode(temp.mat)										<-"numeric"
					
					## only for snps's, future 
					##	v2 will be adapted for InDels
					##	I had to finish my PhD and I didnt have enought time ti make it nicely :)				
					
					pos.qualities<-as.numeric(apply(temp.mat,1,function(x){	
				if(sum((x)[]>0)[]==1){	
					Pr	<-	x[which(x[]>0)[1]]; #	m1	
					return(Pr)
				}else if(sum((x)[]>0)[]==2){									
					a1	<-	x[which(x[]>0)[1]]; #	m1	
					a2	<-	x[which(x[]>0)[2]]; #	m2
					m1	<-a1; m2<-a2
					Pr	<-m1*m2+((m1^2)*(1-m2))/(m1+(1-m2))+((m2^2)*(1-m1))/(m2+(1-m1))
					return(Pr)
				}else if(sum((x)[]>0)[]==3){		
					a1	<-	x[which(x[]>0)[1]]; #	m1	
					a2	<-	x[which(x[]>0)[2]]; #	m2
					a3	<-	x[which(x[]>0)[3]]; #	
					##
					m1	<-a1; m2<-a2
					a12	<-m1*m2+((m1^2)*(1-m2))/(m1+(1-m2))+((m2^2)*(1-m1))/(m2+(1-m1))
					##
					m1	<-a12; m2<-a3
					Pr	<-m1*m2+((m1^2)*(1-m2))/(m1+(1-m2))+((m2^2)*(1-m1))/(m2+(1-m1))
					return(Pr)
				}else if(sum((x)[]>0)[]==4){		
					a1	<-	x[which(x[]>0)[1]]; #	m1	
					a2	<-	x[which(x[]>0)[2]]; #	m2
					a3	<-	x[which(x[]>0)[3]]; #	
					a4	<-	x[which(x[]>0)[4]]; #
					##
					m1	<-a1; m2<-a2
					a12	<-m1*m2+((m1^2)*(1-m2))/(m1+(1-m2))+((m2^2)*(1-m1))/(m2+(1-m1))
					##
					m1	<-a12; m2<-a3
					a123<-m1*m2+((m1^2)*(1-m2))/(m1+(1-m2))+((m2^2)*(1-m1))/(m2+(1-m1))
					##
					m1	<-a123; m2<-a4
					Pr	<-m1*m2+((m1^2)*(1-m2))/(m1+(1-m2))+((m2^2)*(1-m1))/(m2+(1-m1))
					return(Pr)					
				}
				}))	
					multiSNV[pos.to.use,14]<-pos.qualities
				}else{}
									
				# save:	
				setwd(Res.dir)
				save.name		<-paste(sample.name,"__",method,".multiSNV",sep="")
				write.table(multiSNV	, file=save.name,sep="\t",quote=FALSE,row.names=FALSE)
				cat("[ .",set.nr,". ] - [ .",s.nr,". ] - ","saved as: ",save.name,"\n")	
				
				#additional method - for free						
				if(set.nr[]==1){
					info.mat				<-as.character(multiSNV[,c(16),drop=TRUE])
					arr.snp.bin				<-as.matrix(multiSNV[,c(10:13)]); dim(arr.snp.bin)
					arr.snp.bin[,]			<-0
					mode(arr.snp.bin)		<-"numeric"
									
					# transform
					snps					<-c("A","C","G","T")
					snp.nr<-1; is.nr<-1; 
					for(snp.nr in 1:length(snps)){
						arr.snp.bin[,snp.nr]<-as.numeric(grepl(snps[snp.nr],info.mat))}
					mode(arr.snp.bin)		<-"numeric"
					#cbind(arr.snp.bin,info.mat)	
					
					# put back
					multiSNV[,c(10:13)]		<-arr.snp.bin
									
					# postions with na
					pos.with.snps			<-as.matrix(multiSNV[,c(10:13)])
					mode(pos.with.snps)		<-"numeric"
					multiSNV[,c(14)]		<-as.numeric(apply(pos.with.snps,1,function(x){sum(x)[]>0})[]==TRUE)
						
					# snp.presence in a consensus
					snps.n<-c("A","C","G","T")
					snp.pres<-apply(multiSNV[,c(10:13)],1,function(x){paste(snps.n[which(x[]>0)],collapse="/")})
					multiSNV[,15]<-snp.pres				
					
					# put "na"
					multiSNV[which(multiSNV[,15]==""),15]<-"na"		
														
					# save.name
					setwd(Res.dir)
					save.name		<-paste(sample.name,"__AllRecorded_Nucleotides.multiSNV",sep="")
					write.table(multiSNV	, file=save.name,sep="\t",quote=FALSE,row.names=FALSE)
					cat("[ .",set.nr,". ] - plus - ","saved as: ",save.name,"\n")
															
					### in addition we produce small file with names of replicates that were used to build each consensus
					ReplicateNames	<-as.matrix(strsplit(one.s.arr[1,11],split=";")[[1]])
					###
					setwd(Res.dir)
					save.name		<-paste(sample.name,"__Info.ReplicateNames",sep="")
					write.table(ReplicateNames, file=save.name,sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
					cat("[ .",set.nr,". ] - additionally - ","saved as: ",save.name,"\n");cat("\n")						
		
				}else{cat("\n")}# done:				
		
			}#s.nr
		}#set.nr
	}else{	
		# /1. info
		cat("\n" ); cat(".........................................................................","\n" )
		cat("***","\n" )
		cat("***","step 8 preparing multiSNV - NOT DONE","\n" )
		cat("***","Caution, this is the step where AllRecorded_Nucleotides.multiSNV files are created","\n" )
		cat("***","as a bonus to any other method","\n" )
		cat("***",date(),"\n" )
		cat("***","\n" ) 
		cat(".........................................................................","\n" ); cat("\n")	
	}

	####	####################################################################
		cat("\n")
		cat("DONE :)","\n")
		cat("have a nice day","\n")
		cat("PS - Improve, change, do better and faster:)","\n")
	####	####################################################################
}
###	* ---------------------------------------------------------------------------------------------------
###	 END GapJumper
###	* ---------------------------------------------------------------------------------------------------



























###	---------------------------------------------------------------------------------------------------
###	19. GapJumper.plus.all.methods		4650 / upgrade at 2016.06.03 again 3am	
###	---------------------------------------------------------------------------------------------------
Run.All.methods<-function(
						prj.n,tab.n,
						vcf.dir,int.dir,res.dir,tab.dir,
						pos.y,pos.n,
						gp.dl,gp.p,
						se.na.nr,se.na.rt,
						ft.na.rt,ft.nr,ft.ap){
	
	# /1.Parameters defined by the user:			
			
		# main:	
			GapJumper_Data.dir					<-vcf.dir	# directory, with .simplified.vcf
			GapJumper_Intermediate_Files.dir	<-int.dir	# directory, for intermediate files
			GapJumper_Results.dir				<-res.dir	# directory for .multiSNV
			GapJumper_RgTable.dir				<-tab.dir	# directory with the table
			Project.Name							<-prj.n		# project name, make sure it is uniquie and has no "__"
			RgTable.Name							<-tab.n		# table name, make sure it is uniquie and has no "__"
			Positions_to.use					<-pos.y		# or vector, [0,1], length = numbe of positions in .simplifiedVCF
			Positions_to.exclude					<-pos.n		# or vector, [0,1], length = numbe of positions in .simplifiedVCF
			z.classifier							<-pos.z		# or vector, [0,1], length = numbe of positions in .simplifiedVCF
			
		# for each method	
			GapJumper.detection.limit								<-gp.dl
			GapJumper.prior.p										<-gp.p
			FixedThresholds.detection.limit							<-gp.dl# the same			
			FixedThresholds.min.allele.perc							<-ft.ap
			FixedThresholds.min.nr.of.samples.with.same.allele		<-ft.nr
			FixedThresholds.max.na.positions.perc					<-ft.na.rt
			SiteElimination__na.max.nr.per.locus						<-se.na.nr
			SiteElimination__eliminate.pos.with.more.than.na.perc	<-se.na.rt	
				
	# /2.Preset parameters:	
		#	methods (in total it is 6 methods); all recorder positions .multiSNV is always done.
			Run.GapJumper.method							<-"y"	# or "n" ; if n, alfa and beta method wont be run				
			Run.GapJumper.alfa.method					<-"y"	# or "n"
			Run.GapJumper.beta.method					<-"y"	# or "n"
			Run.fixed.threshold.method					<-"y"	# or "n"	
			Run.site.elimination.method					<-"y"	# or "n"
			Run.conservative.method						<-"y"	# or "n"		
		
		#	When R = 0, where R is number of reads for a given nucleotide
			value.for.zero								<-0.000000000001# i put this value because somewhere was a bug, i.e. dividing by zero
			# this value will be removed automatically from a consensus 
				
		#  	parts that should be run each time, but they dont have to if you plan large project
		#	or you have made some mistakes		
			Run.process.process.subVCF					<-"yes"
			Run.array.creator_and.add.RG					<-"yes"	
			Run.process.sample.cross.validation			<-"yes"	
			Run.process.ACGTZ.array.creator				<-"yes"	
			make.info.arrays							<-"yes"
			find.indels.in.each.samaple					<-"yes"
			make.multiSNV.files							<-"yes"
			put.replicate.names.into.multiSNV			<-"yes"	# if no, the last column will be missing in .multiSNV
																# saves a lot of memory
		# 	Positions
			Position_Names								<-"n" 
			# it is best to not change it; becuase the last step wont work properly 	 	
			# it was done for me, 
			# I was not using .multiSNV but the data stored in intermediate files
			# called collased ACGTZ info array
	
	# /3.Run genral function														
		Run.GapJumper.five.methods(
			GapJumper_Data.dir,
			GapJumper_Intermediate_Files.dir,
			GapJumper_Results.dir,
			GapJumper_RgTable.dir,
			Project.Name	,
			RgTable.Name	,
			Position_Names,
			Positions_to.use,
			Positions_to.exclude,
			z.classifier	,
			Run.GapJumper.method,			
			Run.GapJumper.beta.method,
			Run.fixed.threshold.method,
			Run.site.elimination.method,
			Run.conservative.method,
			GapJumper.detection.limit,
			GapJumper.prior.p,
			FixedThresholds.detection.limit,				
			FixedThresholds.min.allele.perc,						
			FixedThresholds.min.nr.of.samples.with.same.allele,	
			FixedThresholds.max.na.positions.perc,					
			SiteElimination__na.max.nr.per.locus	,					
			SiteElimination__eliminate.pos.with.more.than.na.perc,	
			value.for.zero,				
			Run.process.process.subVCF,
			Run.array.creator_and.add.RG	,
			Run.process.sample.cross.validation,
			Run.process.ACGTZ.array.creator,	
			make.info.arrays,
			find.indels.in.each.samaple,
			make.multiSNV.files,
			put.replicate.names.into.multiSNV,
			Run.GapJumper.alfa.method)
}# end six method






























###	---------------------------------------------------------------------------------------------------
###	20. Run.Probabilistic.method		4780 / upgrade at 2016.06.03 again 3am	
###	---------------------------------------------------------------------------------------------------

Run.Probabilistic.method<-function(
						prj.n,tab.n,
						vcf.dir,int.dir,res.dir,tab.dir,
						pos.y,pos.n,pos.z,
						gp.dl,gp.p){


	# /0. parampeters which are not used in this function. 
		# /3. GapJumper method
		ft.dl					<-"n"
		ft.ap					<-"n"
		ft.nr					<-"n"
		ft.na.rt					<-"n"
		se.na.nr					<-"n"
		se.na.rt					<-"n"

	# /1.Parameters defined by the user:			
			
		# main:	
			GapJumper_Data.dir					<-vcf.dir	# directory, with .simplified.vcf
			GapJumper_Intermediate_Files.dir	<-int.dir	# directory, for intermediate files
			GapJumper_Results.dir				<-res.dir	# directory for .multiSNV
			GapJumper_RgTable.dir				<-tab.dir	# directory with the table
			Project.Name							<-prj.n		# project name, make sure it is uniquie and has no "__"
			RgTable.Name							<-tab.n		# table name, make sure it is uniquie and has no "__"
			Positions_to.use					<-pos.y		# or vector, [0,1], length = numbe of positions in .simplifiedVCF
			Positions_to.exclude					<-pos.n		# or vector, [0,1], length = numbe of positions in .simplifiedVCF
			z.classifier							<-pos.z		# or vector, [0,1], length = numbe of positions in .simplifiedVCF
			
		# for each method	
			GapJumper.detection.limit								<-gp.dl
			GapJumper.prior.p										<-gp.p
			FixedThresholds.detection.limit							<-gp.dl# the same			
			FixedThresholds.min.allele.perc							<-ft.ap
			FixedThresholds.min.nr.of.samples.with.same.allele		<-ft.nr
			FixedThresholds.max.na.positions.perc					<-ft.na.rt
			SiteElimination__na.max.nr.per.locus						<-se.na.nr
			SiteElimination__eliminate.pos.with.more.than.na.perc	<-se.na.rt	
				
	# /2.Preset parameters:	
		#	methods (in total it is 6 methods); all recorder positions .multiSNV is always done.
			Run.GapJumper.method							<-"y"	# or "n" ; if n, alfa and beta method wont be run				
			Run.GapJumper.alfa.method					<-"y"	# or "n"
			Run.GapJumper.beta.method					<-"n"	# or "n"
			Run.fixed.threshold.method					<-"n"	# or "n"	
			Run.site.elimination.method					<-"n"	# or "n"
			Run.conservative.method						<-"n"	# or "n"		
		
		#	When R = 0, where R is number of reads for a given nucleotide
			value.for.zero								<-0.000000000001# i put this value because somewhere was a bug, i.e. dividing by zero
	
		#  	parts that should be run each time, but they dont have to if you plan large project
		#	or you have made some mistakes		
			Run.process.process.subVCF					<-"yes"
			Run.array.creator_and.add.RG					<-"yes"	
			Run.process.sample.cross.validation			<-"yes"	
			Run.process.ACGTZ.array.creator				<-"yes"	
			make.info.arrays							<-"yes"
			find.indels.in.each.samaple					<-"yes"
			make.multiSNV.files							<-"yes"
			put.replicate.names.into.multiSNV			<-"yes"	# if no, the last column will be missing in .multiSNV
																# saves a lot of memory
		# 	Positions
			Position_Names								<-"n" 
			# it is best to not change it; becuase the last step wont work properly 	 	
			# it was done for me, 
			# I was not using .multiSNV but the data stored in intermediate files
			# called collased ACGTZ info array
	
	# /3.Run genral function														
		Run.GapJumper.five.methods(
			GapJumper_Data.dir,
			GapJumper_Intermediate_Files.dir,
			GapJumper_Results.dir,
			GapJumper_RgTable.dir,
			Project.Name	,
			RgTable.Name	,
			Position_Names,
			Positions_to.use,
			Positions_to.exclude,
			z.classifier	,
			Run.GapJumper.method,			
			Run.GapJumper.beta.method,
			Run.fixed.threshold.method,
			Run.site.elimination.method,
			Run.conservative.method,
			GapJumper.detection.limit,
			GapJumper.prior.p,
			FixedThresholds.detection.limit,				
			FixedThresholds.min.allele.perc,						
			FixedThresholds.min.nr.of.samples.with.same.allele,	
			FixedThresholds.max.na.positions.perc,					
			SiteElimination__na.max.nr.per.locus	,					
			SiteElimination__eliminate.pos.with.more.than.na.perc,	
			value.for.zero,				
			Run.process.process.subVCF,
			Run.array.creator_and.add.RG	,
			Run.process.sample.cross.validation,
			Run.process.ACGTZ.array.creator,	
			make.info.arrays,
			find.indels.in.each.samaple,
			make.multiSNV.files,
			put.replicate.names.into.multiSNV,
			Run.GapJumper.alfa.method)
}# end six method









































###	---------------------------------------------------------------------------------------------------
###	21. Run.Probabilistic.beta.method		4930 / upgrade at 2016.06.03 again 3am	
###	---------------------------------------------------------------------------------------------------
Run.Probabilistic.beta.method<-function(
						prj.n,tab.n,
						vcf.dir,int.dir,res.dir,tab.dir,
						pos.y,pos.n,pos.z,
						gp.dl,gp.p){


	# /0. parampeters which are not used in this function. 
		# /3. GapJumper method
		ft.dl					<-"n"
		ft.ap					<-"n"
		ft.nr					<-"n"
		ft.na.rt					<-"n"
		se.na.nr					<-"n"
		se.na.rt					<-"n"


	# /1.Parameters defined by the user:			
			
		# main:	
			GapJumper_Data.dir					<-vcf.dir	# directory, with .simplified.vcf
			GapJumper_Intermediate_Files.dir	<-int.dir	# directory, for intermediate files
			GapJumper_Results.dir				<-res.dir	# directory for .multiSNV
			GapJumper_RgTable.dir				<-tab.dir	# directory with the table
			Project.Name							<-prj.n		# project name, make sure it is uniquie and has no "__"
			RgTable.Name							<-tab.n		# table name, make sure it is uniquie and has no "__"
			Positions_to.use					<-pos.y		# or vector, [0,1], length = numbe of positions in .simplifiedVCF
			Positions_to.exclude					<-pos.n		# or vector, [0,1], length = numbe of positions in .simplifiedVCF
			z.classifier							<-pos.z		# or vector, [0,1], length = numbe of positions in .simplifiedVCF
			
		# for each method	
			GapJumper.detection.limit								<-gp.dl
			GapJumper.prior.p										<-gp.p
			FixedThresholds.detection.limit							<-gp.dl# the same		
			FixedThresholds.min.allele.perc							<-ft.ap
			FixedThresholds.min.nr.of.samples.with.same.allele		<-ft.nr
			FixedThresholds.max.na.positions.perc					<-ft.na.rt
			SiteElimination__na.max.nr.per.locus						<-se.na.nr
			SiteElimination__eliminate.pos.with.more.than.na.perc	<-se.na.rt	
				
	# /2.Preset parameters:	
		#	methods (in total it is 6 methods); all recorder positions .multiSNV is always done.
			Run.GapJumper.method							<-"y"	# or "n" ; if n, alfa and beta method wont be run				
			Run.GapJumper.alfa.method					<-"n"	# or "n"
			Run.GapJumper.beta.method					<-"y"	# or "n"
			Run.fixed.threshold.method					<-"n"	# or "n"	
			Run.site.elimination.method					<-"n"	# or "n"
			Run.conservative.method						<-"n"	# or "n"		
		
		#	When R = 0, where R is number of reads for a given nucleotide
			value.for.zero								<-0.000000000001# i put this value because somewhere was a bug, i.e. dividing by zero
	
		#  	parts that should be run each time, but they dont have to if you plan large project
		#	or you have made some mistakes		
			Run.process.process.subVCF					<-"yes"
			Run.array.creator_and.add.RG					<-"yes"	
			Run.process.sample.cross.validation			<-"yes"	
			Run.process.ACGTZ.array.creator				<-"yes"	
			make.info.arrays							<-"yes"
			find.indels.in.each.samaple					<-"yes"
			make.multiSNV.files							<-"yes"
			put.replicate.names.into.multiSNV			<-"yes"	# if no, the last column will be missing in .multiSNV
																# saves a lot of memory
		# 	Positions
			Position_Names								<-"n" 
			# it is best to not change it; becuase the last step wont work properly 	 	
			# it was done for me, 
			# I was not using .multiSNV but the data stored in intermediate files
			# called collased ACGTZ info array
	
	# /3.Run genral function														
		Run.GapJumper.five.methods(
			GapJumper_Data.dir,
			GapJumper_Intermediate_Files.dir,
			GapJumper_Results.dir,
			GapJumper_RgTable.dir,
			Project.Name	,
			RgTable.Name	,
			Position_Names,
			Positions_to.use,
			Positions_to.exclude,
			z.classifier	,
			Run.GapJumper.method,			
			Run.GapJumper.beta.method,
			Run.fixed.threshold.method,
			Run.site.elimination.method,
			Run.conservative.method,
			GapJumper.detection.limit,
			GapJumper.prior.p,
			FixedThresholds.detection.limit,				
			FixedThresholds.min.allele.perc,						
			FixedThresholds.min.nr.of.samples.with.same.allele,	
			FixedThresholds.max.na.positions.perc,					
			SiteElimination__na.max.nr.per.locus	,					
			SiteElimination__eliminate.pos.with.more.than.na.perc,	
			value.for.zero,				
			Run.process.process.subVCF,
			Run.array.creator_and.add.RG	,
			Run.process.sample.cross.validation,
			Run.process.ACGTZ.array.creator,	
			make.info.arrays,
			find.indels.in.each.samaple,
			make.multiSNV.files,
			put.replicate.names.into.multiSNV,
			Run.GapJumper.alfa.method)
}# end six method





















###	---------------------------------------------------------------------------------------------------
###	22. GapJumper.plus.cons.method		5060 / upgrade at 2016.06.03 again 3am	
###	---------------------------------------------------------------------------------------------------

Run.Conservative.method	<-function(
	prj.n,tab.n,
	vcf.dir,int.dir,res.dir,tab.dir,
	pos.y,pos.n,
	gp.dl){
	
	# /0. parampeters which are not used in this function.
		gp.p					<-as.numeric(0.5)
		pos.z					<-"n"
	
	# /0. parampeters which are not used in this function. 
		# /3. GapJumper method
		ft.dl					<-"n"
		ft.ap					<-"n"
		ft.nr					<-"n"
		ft.na.rt					<-"n"
		se.na.nr					<-"n"
		se.na.rt					<-"n"

	# /1.Parameters defined by the user:			
			
		# main:	
			GapJumper_Data.dir					<-vcf.dir	# directory, with .simplified.vcf
			GapJumper_Intermediate_Files.dir	<-int.dir	# directory, for intermediate files
			GapJumper_Results.dir				<-res.dir	# directory for .multiSNV
			GapJumper_RgTable.dir				<-tab.dir	# directory with the table
			Project.Name							<-prj.n		# project name, make sure it is uniquie and has no "__"
			RgTable.Name							<-tab.n		# table name, make sure it is uniquie and has no "__"
			Positions_to.use					<-pos.y		# or vector, [0,1], length = numbe of positions in .simplifiedVCF
			Positions_to.exclude					<-pos.n		# or vector, [0,1], length = numbe of positions in .simplifiedVCF
			z.classifier							<-pos.z		# or vector, [0,1], length = numbe of positions in .simplifiedVCF
			
		# for each method	
			GapJumper.detection.limit								<-gp.dl
			GapJumper.prior.p										<-gp.p
			FixedThresholds.detection.limit							<-gp.dl# the same			
			FixedThresholds.min.allele.perc							<-ft.ap
			FixedThresholds.min.nr.of.samples.with.same.allele		<-ft.nr
			FixedThresholds.max.na.positions.perc					<-ft.na.rt
			SiteElimination__na.max.nr.per.locus						<-se.na.nr
			SiteElimination__eliminate.pos.with.more.than.na.perc	<-se.na.rt	
				
	# /2.Preset parameters:	
		#	methods (in total it is 6 methods); all recorder positions .multiSNV is always done.
			Run.GapJumper.method							<-"n"	# or "n" ; if n, alfa and beta method wont be run				
			Run.GapJumper.alfa.method					<-"n"	# or "n"
			Run.GapJumper.beta.method					<-"n"	# or "n"
			Run.fixed.threshold.method					<-"n"	# or "n"	
			Run.site.elimination.method					<-"n"	# or "n"
			Run.conservative.method						<-"y"	# or "n"		
		
		#	When R = 0, where R is number of reads for a given nucleotide
			value.for.zero								<-0.000000000001# i put this value because somewhere was a bug, i.e. dividing by zero
	
		#  	parts that should be run each time, but they dont have to if you plan large project
		#	or you have made some mistakes		
			Run.process.process.subVCF					<-"yes"
			Run.array.creator_and.add.RG					<-"yes"	
			Run.process.sample.cross.validation			<-"yes"	
			Run.process.ACGTZ.array.creator				<-"yes"	
			make.info.arrays							<-"yes"
			find.indels.in.each.samaple					<-"yes"
			make.multiSNV.files							<-"yes"
			put.replicate.names.into.multiSNV			<-"yes"	# if no, the last column will be missing in .multiSNV
																# saves a lot of memory
		# 	Positions
			Position_Names								<-"n" 
			# it is best to not change it; becuase the last step wont work properly 	 	
			# it was done for me, 
			# I was not using .multiSNV but the data stored in intermediate files
			# called collased ACGTZ info array
	
	# /3.Run genral function														
		Run.GapJumper.five.methods(
			GapJumper_Data.dir,
			GapJumper_Intermediate_Files.dir,
			GapJumper_Results.dir,
			GapJumper_RgTable.dir,
			Project.Name	,
			RgTable.Name	,
			Position_Names,
			Positions_to.use,
			Positions_to.exclude,
			z.classifier	,
			Run.GapJumper.method,			
			Run.GapJumper.beta.method,
			Run.fixed.threshold.method,
			Run.site.elimination.method,
			Run.conservative.method,
			GapJumper.detection.limit,
			GapJumper.prior.p,
			FixedThresholds.detection.limit,				
			FixedThresholds.min.allele.perc,						
			FixedThresholds.min.nr.of.samples.with.same.allele,	
			FixedThresholds.max.na.positions.perc,					
			SiteElimination__na.max.nr.per.locus	,					
			SiteElimination__eliminate.pos.with.more.than.na.perc,	
			value.for.zero,				
			Run.process.process.subVCF,
			Run.array.creator_and.add.RG	,
			Run.process.sample.cross.validation,
			Run.process.ACGTZ.array.creator,	
			make.info.arrays,
			find.indels.in.each.samaple,
			make.multiSNV.files,
			put.replicate.names.into.multiSNV,
			Run.GapJumper.alfa.method)
}# end six method
















































###	---------------------------------------------------------------------------------------------------
###	23. GapJumper.plus.SiteElim.method		5220 / upgrade at 2016.06.03 again 3am	
###	---------------------------------------------------------------------------------------------------
Run.Site.Elimination.method<-function(
						prj.n,tab.n,
						vcf.dir,int.dir,res.dir,tab.dir,
						pos.y,pos.n,
						gp.dl,
						se.na.nr,se.na.rt)	{

	# /0. parampeters which are not used in this function.
		gp.p					<-as.numeric(0.5)
		pos.z				<-"n"

	# /0. parampeters which are not used in this function. 
		# /3. GapJumper method
		ft.dl					<-"n"
		ft.ap					<-"n"
		ft.nr					<-"n"
		ft.na.rt					<-"n"

	# /1.Parameters defined by the user:			
			
		# main:	
			GapJumper_Data.dir					<-vcf.dir	# directory, with .simplified.vcf
			GapJumper_Intermediate_Files.dir	<-int.dir	# directory, for intermediate files
			GapJumper_Results.dir				<-res.dir	# directory for .multiSNV
			GapJumper_RgTable.dir				<-tab.dir	# directory with the table
			Project.Name							<-prj.n		# project name, make sure it is uniquie and has no "__"
			RgTable.Name							<-tab.n		# table name, make sure it is uniquie and has no "__"
			Positions_to.use					<-pos.y		# or vector, [0,1], length = numbe of positions in .simplifiedVCF
			Positions_to.exclude					<-pos.n		# or vector, [0,1], length = numbe of positions in .simplifiedVCF
			z.classifier							<-pos.z		# or vector, [0,1], length = numbe of positions in .simplifiedVCF
			
			
		# for each method	
			GapJumper.detection.limit								<-gp.dl
			GapJumper.prior.p										<-gp.p
			FixedThresholds.detection.limit							<-gp.dl# the same		
			FixedThresholds.min.allele.perc							<-ft.ap
			FixedThresholds.min.nr.of.samples.with.same.allele		<-ft.nr
			FixedThresholds.max.na.positions.perc					<-ft.na.rt
			SiteElimination__na.max.nr.per.locus						<-se.na.nr
			SiteElimination__eliminate.pos.with.more.than.na.perc	<-se.na.rt	
				
	# /2.Preset parameters:	
		#	methods (in total it is 6 methods); all recorder positions .multiSNV is always done.
			Run.GapJumper.method							<-"y"	# or "n" ; if n, alfa and beta method wont be run				
			Run.GapJumper.alfa.method					<-"y"	# or "n"
			Run.GapJumper.beta.method					<-"n"	# or "n"
			Run.fixed.threshold.method					<-"n"	# or "n"	
			Run.site.elimination.method					<-"n"	# or "n"
			Run.conservative.method						<-"n"	# or "n"		
		
		#	When R = 0, where R is number of reads for a given nucleotide
			value.for.zero								<-0.000000000001# i put this value because somewhere was a bug, i.e. dividing by zero
	
		#  	parts that should be run each time, but they dont have to if you plan large project
		#	or you have made some mistakes		
			Run.process.process.subVCF					<-"yes"
			Run.array.creator_and.add.RG					<-"yes"	
			Run.process.sample.cross.validation			<-"yes"	
			Run.process.ACGTZ.array.creator				<-"yes"	
			make.info.arrays							<-"yes"
			find.indels.in.each.samaple					<-"yes"
			make.multiSNV.files							<-"yes"
			put.replicate.names.into.multiSNV			<-"yes"	# if no, the last column will be missing in .multiSNV
																# saves a lot of memory
		# 	Positions
			Position_Names								<-"n" 
			# it is best to not change it; becuase the last step wont work properly 	 	
			# it was done for me, 
			# I was not using .multiSNV but the data stored in intermediate files
			# called collased ACGTZ info array
	
	# /3.Run genral function														
		Run.GapJumper.five.methods(
			GapJumper_Data.dir,
			GapJumper_Intermediate_Files.dir,
			GapJumper_Results.dir,
			GapJumper_RgTable.dir,
			Project.Name	,
			RgTable.Name	,
			Position_Names,
			Positions_to.use,
			Positions_to.exclude,
			z.classifier	,
			Run.GapJumper.method,			
			Run.GapJumper.beta.method,
			Run.fixed.threshold.method,
			Run.site.elimination.method,
			Run.conservative.method,
			GapJumper.detection.limit,
			GapJumper.prior.p,
			FixedThresholds.detection.limit,				
			FixedThresholds.min.allele.perc,						
			FixedThresholds.min.nr.of.samples.with.same.allele,	
			FixedThresholds.max.na.positions.perc,					
			SiteElimination__na.max.nr.per.locus	,					
			SiteElimination__eliminate.pos.with.more.than.na.perc,	
			value.for.zero,				
			Run.process.process.subVCF,
			Run.array.creator_and.add.RG	,
			Run.process.sample.cross.validation,
			Run.process.ACGTZ.array.creator,	
			make.info.arrays,
			find.indels.in.each.samaple,
			make.multiSNV.files,
			put.replicate.names.into.multiSNV,
			Run.GapJumper.alfa.method)
}# end six method



















###	---------------------------------------------------------------------------------------------------
###	24. GapJumper.plus.FixedTr.method		5350 / upgrade at 2016.06.03 again 3am	
###	---------------------------------------------------------------------------------------------------

Run.Fixed.Thresholds.method<-function(
						prj.n,tab.n,
						vcf.dir,int.dir,res.dir,tab.dir,
						pos.y,pos.n,
						gp.dl,
						ft,na.rt,ft.nr,ft.ap)	{


	# /0. parampeters which are not used in this function.
		gp.p					<-as.numeric(0.5)# just for making a name
		pos.z				<-"n"

	# /0. parampeters which are not used in this function. 
		# /3. GapJumper method
		se.na.nr					<-"n"
		se.na.rt					<-"n"

	# /1.Parameters defined by the user:			
			
		# main:	
			GapJumper_Data.dir					<-vcf.dir	# directory, with .simplified.vcf
			GapJumper_Intermediate_Files.dir	<-int.dir	# directory, for intermediate files
			GapJumper_Results.dir				<-res.dir	# directory for .multiSNV
			GapJumper_RgTable.dir				<-tab.dir	# directory with the table
			Project.Name							<-prj.n		# project name, make sure it is uniquie and has no "__"
			RgTable.Name							<-tab.n		# table name, make sure it is uniquie and has no "__"
			Positions_to.use					<-pos.y		# or vector, [0,1], length = numbe of positions in .simplifiedVCF
			Positions_to.exclude					<-pos.n		# or vector, [0,1], length = numbe of positions in .simplifiedVCF
			z.classifier							<-pos.z		# or vector, [0,1], length = numbe of positions in .simplifiedVCF
			
		# for each method	
			GapJumper.detection.limit								<-gp.dl
			GapJumper.prior.p										<-gp.p
			FixedThresholds.detection.limit							<-gp.dl# the same	
			FixedThresholds.min.allele.perc							<-ft.ap
			FixedThresholds.min.nr.of.samples.with.same.allele		<-ft.nr
			FixedThresholds.max.na.positions.perc					<-ft.na.rt
			SiteElimination__na.max.nr.per.locus						<-se.na.nr
			SiteElimination__eliminate.pos.with.more.than.na.perc	<-se.na.rt	
				
	# /2.Preset parameters:	
		#	methods (in total it is 6 methods); all recorder positions .multiSNV is always done.
			Run.GapJumper.method							<-"y"	# or "n" ; if n, alfa and beta method wont be run				
			Run.GapJumper.alfa.method					<-"y"	# or "n"
			Run.GapJumper.beta.method					<-"n"	# or "n"
			Run.fixed.threshold.method					<-"n"	# or "n"	
			Run.site.elimination.method					<-"n"	# or "n"
			Run.conservative.method						<-"n"	# or "n"		
		
		#	When R = 0, where R is number of reads for a given nucleotide
			value.for.zero								<-0.000000000001# i put this value because somewhere was a bug, i.e. dividing by zero
	
		#  	parts that should be run each time, but they dont have to if you plan large project
		#	or you have made some mistakes		
			Run.process.process.subVCF					<-"yes"
			Run.array.creator_and.add.RG					<-"yes"	
			Run.process.sample.cross.validation			<-"yes"	
			Run.process.ACGTZ.array.creator				<-"yes"	
			make.info.arrays							<-"yes"
			find.indels.in.each.samaple					<-"yes"
			make.multiSNV.files							<-"yes"
			put.replicate.names.into.multiSNV			<-"yes"	# if no, the last column will be missing in .multiSNV
																# saves a lot of memory
		# 	Positions
			Position_Names								<-"n" 
			# it is best to not change it; becuase the last step wont work properly 	 	
			# it was done for me, 
			# I was not using .multiSNV but the data stored in intermediate files
			# called collased ACGTZ info array
	
	# /3.Run genral function														
		Run.GapJumper.five.methods(
			GapJumper_Data.dir,
			GapJumper_Intermediate_Files.dir,
			GapJumper_Results.dir,
			GapJumper_RgTable.dir,
			Project.Name	,
			RgTable.Name	,
			Position_Names,
			Positions_to.use,
			Positions_to.exclude,
			z.classifier	,
			Run.GapJumper.method,			
			Run.GapJumper.beta.method,
			Run.fixed.threshold.method,
			Run.site.elimination.method,
			Run.conservative.method,
			GapJumper.detection.limit,
			GapJumper.prior.p,
			FixedThresholds.detection.limit,				
			FixedThresholds.min.allele.perc,						
			FixedThresholds.min.nr.of.samples.with.same.allele,	
			FixedThresholds.max.na.positions.perc,					
			SiteElimination__na.max.nr.per.locus	,					
			SiteElimination__eliminate.pos.with.more.than.na.perc,	
			value.for.zero,				
			Run.process.process.subVCF,
			Run.array.creator_and.add.RG	,
			Run.process.sample.cross.validation,
			Run.process.ACGTZ.array.creator,	
			make.info.arrays,
			find.indels.in.each.samaple,
			make.multiSNV.files,
			put.replicate.names.into.multiSNV,
			Run.GapJumper.alfa.method)
}# end six method








































###	---------------------------------------------------------------------------------------------------
###	25. Calculate.Noise.threshold		5500
###	---------------------------------------------------------------------------------------------------
	Calculate.Noise.threshold<-function(res.dir,prob.multiSNV.name,all.Rec.nucl.multiSNV.name,retain.nucleotides){
	
	# load:
	setwd(res.dir); 
	prob.multi			<-as.matrix(read.table(prob.multiSNV.name,header=TRUE,sep="\t",as.is=T,strip.white=T,colClasses="character"))
	all.nucl.multi		<-as.matrix(read.table(all.Rec.nucl.multiSNV.name,header=TRUE,sep="\t",as.is=T,strip.white=T,colClasses="character")	)
	dim.with.snps		<-c(11:14);colnames(prob.multi)[dim.with.snps]
	dim.with.indel.pres	<-19;colnames(prob.multi)[dim.with.indel.pres]
	
	# exclude pos. with indels
	ind.pres				<-as.numeric(grepl("y",prob.multi[,dim.with.indel.pres]))
	prob.multi			<-prob.multi[which(ind.pres[]==0),]	
	all.nucl.multi		<-all.nucl.multi[which(ind.pres[]==0),]	


	# take out the positions with nucleotides
	prob.snps			<-prob.multi[,dim.with.snps]
	all.nucl.snps		<-all.nucl.multi[,dim.with.snps]
	mode(all.nucl.snps)	<-"numeric"
	snps.probabilities	<-sort(as.numeric(prob.snps[which(all.nucl.snps[,]==1)]))
	
	# caluclate noise threshold
	if(retain.nucleotides[]<1){
		remove			<-round(length(snps.probabilities)*c(1-retain.nucleotides),digits=0)		
		if(remove[]==0){
		Ns.Tr			<-c(snps.probabilities[1]-c(0.1*snps.probabilities[1]))					
		}else{Ns.Tr<-snps.probabilities[remove]}
	}else{ 
		Ns.Tr			<-c(snps.probabilities[1]-c(0.1*snps.probabilities[1]))	
	}# end if else
	
	# check iy Ns.Tr is not lower than arbitrarly choosen value for empty pos.
	if(Ns.Tr[]<=0.000001){
		Ns.Tr			<-0.000002# becaus eit is an absolute zero in my system
	}else{}
	
	# number of nucleotides that are retaind or excluded with this threshold
	exc.nucl.nr	<-length(which(snps.probabilities[]<Ns.Tr))
	ret.nucl.nr	<-length(which(snps.probabilities[]>=Ns.Tr))
	
	# hist
	hist(snps.probabilities,col="gray",border="gray",ylab="Nulcleotide Frequency",xlab="Pr(snp)",
		main=paste("Noise Threshold:",Ns.Tr,sep=""),sub=paste("Retained: ",ret.nucl.nr,";  Excluded: ",exc.nucl.nr,sep=""))
	abline(v=Ns.Tr)
	
	# give back the value
	return(as.numeric(Ns.Tr))	
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
###	---------------------------------------------------------------------------------------------------
###	26. Test.Noise.threshold		5570
###	---------------------------------------------------------------------------------------------------	
	Test.Noise.threshold<-function(res.dir,prob.multiSNV.name,all.Rec.nucl.multiSNV.name,Noise.Threshold	){
	
	# load:
	setwd(res.dir); 
	prob.multi			<-as.matrix(read.table(prob.multiSNV.name,header=TRUE,sep="\t",as.is=T,strip.white=T,colClasses="character"))
	all.nucl.multi		<-as.matrix(read.table(all.Rec.nucl.multiSNV.name,header=TRUE,sep="\t",as.is=T,strip.white=T,colClasses="character")	)
	dim.with.snps		<-c(11:14);colnames(prob.multi)[dim.with.snps]
	dim.with.indel.pres	<-19;colnames(prob.multi)[dim.with.indel.pres]
	
	# exclude pos. with indels
	ind.pres				<-as.numeric(grepl("y",prob.multi[,dim.with.indel.pres]))
	prob.multi			<-prob.multi[which(ind.pres[]==0),]	
	all.nucl.multi		<-all.nucl.multi[which(ind.pres[]==0),]	


	# take out the positions with nucleotides
	prob.snps			<-prob.multi[,dim.with.snps]
	all.nucl.snps		<-all.nucl.multi[,dim.with.snps]
	mode(all.nucl.snps)	<-"numeric"
	snps.probabilities	<-sort(as.numeric(prob.snps[which(all.nucl.snps[,]==1)]))
	
	# caluclate noise threshold
	Ns.Tr				<-Noise.Threshold
	
	# number of nucleotides that are retaind or excluded with this threshold
	exc.nucl.nr	<-length(which(snps.probabilities[]<Ns.Tr))
	ret.nucl.nr	<-length(which(snps.probabilities[]>=Ns.Tr))
	
	# hist
	hist(snps.probabilities,col="gray",border="gray",ylab="Nulcleotide Frequency",xlab="Pr(snp)",
		main=paste("Noise Threshold:",Ns.Tr,sep=""),sub=paste("Retained: ",ret.nucl.nr,";  Excluded: ",exc.nucl.nr,sep=""))
	abline(v=Ns.Tr)
	
	}	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
###	---------------------------------------------------------------------------------------------------
###	27. Remove.Noise.Global		5630
###	---------------------------------------------------------------------------------------------------	

# the function used only if
#	value.for.zero > 0; i.e. when he re is no reads and each of these alleles obtains p>0 !


Remove.Noise.Global<-function(res.dir){			
	# load:
	setwd(res.dir); dir()
	search.for			<-"method.multiSNV"
	dir.list				<-dir()[which(grepl(search.for,dir())[]==TRUE)];dir.list	
	###
	search.for			<-"Probabilistic_"
	gp.met.list			<-dir.list[which(grepl(search.for,dir.list)[]==TRUE)]
	###
	search.for			<-"AllRecorded_Nucleotides"
	all.nucl.list		<-dir()[which(grepl(search.for,dir())[]==TRUE)]
	
	set.nr<-1
	for(set.nr in 1:length(gp.met.list)){
		
		# info
		cat("[ ",set.nr," ] Removing noise - Global - old function ","\n")
		cat("[ ",set.nr," ] left just in case ","\n")
		cat("[ ",set.nr," ] it removes all qualities for all postions that had no reads in any of the replicates","\n")
		cat("[ ",set.nr," ] the qualities are made with DST evidenctial theory, so this nucleotides already are == 0 at these positons ","\n")
		cat("[ ... ] loading:",gp.met.list[set.nr],"\n")	
		
		# load data
		setwd(res.dir)	
		prob.multi			<-as.matrix(read.table(gp.met.list[set.nr],header=TRUE,sep="\t",as.is=T,strip.white=T,colClasses="character"))
		set.name				<-strsplit(gp.met.list[set.nr],split="__")[[1]][1]
		next.set.to.load	<-all.nucl.list[which(grepl(set.name,all.nucl.list)[]==TRUE)]
		
		# error message:
		if(length(next.set.to.load)[]!=1){
		cat("ERROR  - you should have one (and only one !) .multiSNV file with AllRecorded_Nucleotides","\n")
		cat("ERROR  - that corresponds to your .multiSNV file prepared with the probabilitic method","\n")
		}else{}

		# load all nucleotides
		all.nucl.multi		<-as.matrix(read.table(next.set.to.load,header=TRUE,sep="\t",as.is=T,strip.white=T,colClasses="character")	)
		cat("[ ... ] loading:",next.set.to.load,"\n")	
			
		# exctract		
		dim.with.snps		<-c(10:13);colnames(prob.multi)[dim.with.snps]
		prob.snps			<-prob.multi[,dim.with.snps]
		mode(prob.snps)		<-"numeric"
		snps.presence		<-all.nucl.multi[,dim.with.snps]
		mode(snps.presence)	<-"numeric"
		
		# remove noise - globaly
		prob.snps[which(snps.presence[,]==0)]<-0
		prob.multi[,dim.with.snps]	<-prob.snps
							
		# snp.presence in a consensus
		snps.n<-c("A","C","G","T")
		snp.pres<-apply(prob.multi[,c(10:13)],1,function(x){paste(snps.n[which(x[]>0)],collapse="/")})
		prob.multi[,15]<-snp.pres				
					
		# put "na"
		prob.multi[which(prob.multi[,15]==""),14]<-1		
		prob.multi[which(prob.multi[,15]==""),15]<-"na"		
		
		# restore the name bcause i made so many mistakes ;)
		multiSNV<-prob.multi
				
		# save:
		setwd(res.dir)
		save.name		<-paste(strsplit(gp.met.list[set.nr],split=".multiSNV")[[1]][1],"_NsFiltered.multiSNV",sep="")
		write.table(multiSNV, file=save.name,sep="\t",quote=FALSE,row.names=FALSE)
		cat("[ ... ] saved as:   ",save.name,"\n");cat("\n")	
		
	}# end	
}# end function		


































###	---------------------------------------------------------------------------------------------------
###	28. Remove.Noise.using.Noise.Threshold		5750
###	---------------------------------------------------------------------------------------------------	
Remove.Noise.using.Noise.Threshold<-function(in.dir,res.dir,file.name,Ns,per.replicate,max.cutoff){			
		
	
	Ns.Tr<-Ns
	
	# search for:
	setwd(in.dir); dir()
	search.for			<-"multiSNV"
	dir.list				<-dir()[which(grepl(search.for,dir())[]==TRUE)];dir.list	
	###
	search.for			<-file.name
	gp.met.list			<-dir.list[which(grepl(search.for,dir.list)[]==TRUE)]
	gp.met.list	
	
	# just in case
	if(max.cutoff[]>=1){
		max.cutoff<-1
		}else{}
	
	# info:
	cat("***","-","----------------","-","***","\n") 
	cat("***","-","Noise filtration","-","***","\n") 
	cat("***","-",date(),"\n") 
	cat("***","-","----------------","-","\n") 
	cat("---","-","      Ns            =",Ns.Tr,"\n")
	cat("---","-","      per replicate =",per.replicate,"\n")
	cat("---","-","      max cutoff    =",max.cutoff,"\n")
	cat("   ","-","----------------","-","\n") 	
	cat("---","-","Ns.Tr is the percentage of the worst quality nucletides that will be removed","\n") 
	cat("---","-","from each consensus identified with file.name that you provided","\n")
	cat("---","-","if. per replicate == yes","\n")
	cat("---","-","ns.tr * number of replicates used to build a consensus","\n")
	cat("   ","-","----------------","-","\n") 	
	cat("---","-","file identified with:",file.name,"\n") 
	cat("---","-","the number of files =",length(gp.met.list),"\n") 
	cat("---","-","           found in =",in.dir,"\n") 
	cat("   ","-","----------------","-","\n") 	
	cat("---","-","CAUTION","\n") 
	cat("---","-","setting up too low max.cutoff threshold can reduce the % of removed nucleotides","\n") 
	cat("---","-","sometimes many nucleotides can haave the same qualiy score","\n") 
	cat("---","-","than the percentage of removed positions can be higher","\n") 
	cat("   ","-","----------------","-","\n") 
	cat("\n") ;cat("\n") ;cat("\n") 


	# big if.
	if(length(gp.met.list)[]==0){
		cat("ERROR","\n");cat("ERROR","\n");cat("ERROR","\n")
		cat("No files were found in in.dir","\n") 		
		cat("please control if the file.name is correct and or in.dir is correct and not empty","\n") 	
		cat("\n") 	
	}else{
	
	#second if. 
	if(Ns.Tr[]>0){
				
	set.nr<-1
	for(set.nr in 1:length(gp.met.list)){
					
		# PART 1				- Noise filtration
						
		# load data
		setwd(res.dir)	
		prob.multi			<-as.matrix(read.table(gp.met.list[set.nr],header=TRUE,sep="\t",as.is=T,strip.white=T,colClasses="character"))

		# Ns.tr
		repl.nr				<-as.numeric(prob.multi[1,ncol(prob.multi)])
		if(per.replicate[]=="yes"){
			ns.threshold		<-Ns.Tr*repl.nr
			if(ns.threshold[]>=1){ns.threshold[1]<-1}else{}
		}
				
		# exctract		
		dim.with.snps		<-c(10:13);colnames(prob.multi)[dim.with.snps]
		prob.snps			<-prob.multi[,dim.with.snps]
		mode(prob.snps)		<-"numeric"
		
		# calculate cutoff value
		all.nucl.probs		<-as.numeric(prob.snps)
		all.nucl.probs		<-all.nucl.probs[which(all.nucl.probs[]>0)]
		pos.with.cutoff		<-round(length(all.nucl.probs)*ns.threshold,digits=0)	
		cutoff				<-as.numeric(sort(all.nucl.probs,decreasing = FALSE)[pos.with.cutoff])
		
		# check for cutoff
		if(cutoff[]>=as.numeric(max.cutoff)){
			cutoff<-max.cutoff; cutoff
			}else{}
			
		# check how many positions will be removed
			nr.of.removed.nucl				<-length(which(as.numeric(all.nucl.probs)[]<=cutoff))
			perc.of.removed.nucl				<-round(nr.of.removed.nucl/length(all.nucl.probs),digits=4)
			
		# info
			cat("-",set.nr,"-","loading:",gp.met.list[set.nr],"\n")	
			cat("-",set.nr,"-",ns.threshold," - of the worst postions was set to be removed","\n")	
			cat("-",set.nr,"-","the cutoff value =",cutoff,"\n")					
			cat("-",set.nr,"-","which removes    =",nr.of.removed.nucl,"out off",length(all.nucl.probs),"nucleotides at that consensus","\n")
			cat("-",set.nr,"-","                  ","i.e = ",perc.of.removed.nucl,"of nucleotides in this consensus","\n")	
		
		# removing these positions and rebuild consensus	
			prob.snps[which(prob.snps[]<=cutoff)]	<-0
			prob.multi[,c(10:13)]					<-prob.snps
			
						
		# PART 2				- positions with na identification and calulation of qaluty scores for the positon	
			multiSNV			<-prob.multi
								
			# now make allele presence in 		
				snps.n									<-c("A","C","G","T")
				snp.pres									<-apply(multiSNV[,c(10:13)],1,function(x){paste(snps.n[which(x[]>0)],collapse="/")})
				multiSNV[,15]							<-snp.pres
				multiSNV[which(multiSNV[,15]==""),15]	<-"na"
								
			# column for genotype quality (Pr(genotype))									
				empty.column										<-rep(1,nrow(multiSNV))	
				empty.column[which(multiSNV[,15]=="na")]				<-0
				multiSNV[,14]										<-empty.column							
				colnames(multiSNV)[14]								<-"M(genotype)"				
				multiSNV												<-as.matrix(multiSNV)				
			
			# first check if recaltulation is needed
				prob.snps[,]<-1	
				m.values.in.a.consensus	<-as.numeric(names(table(as.numeric(prob.snps)[which(as.numeric(prob.snps)[]>0)])))
				if(sum(m.values.in.a.consensus[]!=1)[]>0){use.pcr5.method<-"no"
				}else{				use.pcr5.method	<-"yes"}
																
			#yes 
			#	caluclate genotype quality for positions 
			#	(only for probabilistic methods)							
				if(use.pcr5.method[]=="yes"){	
					pos.to.use											<-which(empty.column[]>0)
					temp.mat											<-multiSNV[pos.to.use,c(10:13)]
					mode(temp.mat)										<-"numeric"
					
					## only for snps's, future 
					##	v2 will be adapted for InDels
					##	I had to finish my PhD and I didnt have enought time ti make it nicely :)				
					
					pos.qualities<-as.numeric(apply(temp.mat,1,function(x){	
				if(sum((x)[]>0)[]==1){	
					Pr	<-	x[which(x[]>0)[1]]; #	m1	
					return(Pr)
				}else if(sum((x)[]>0)[]==2){									
					a1	<-	x[which(x[]>0)[1]]; #	m1	
					a2	<-	x[which(x[]>0)[2]]; #	m2
					m1	<-a1; m2<-a2
					Pr	<-m1*m2+((m1^2)*(1-m2))/(m1+(1-m2))+((m2^2)*(1-m1))/(m2+(1-m1))
					return(Pr)
				}else if(sum((x)[]>0)[]==3){		
					a1	<-	x[which(x[]>0)[1]]; #	m1	
					a2	<-	x[which(x[]>0)[2]]; #	m2
					a3	<-	x[which(x[]>0)[3]]; #	
					##
					m1	<-a1; m2<-a2
					a12	<-m1*m2+((m1^2)*(1-m2))/(m1+(1-m2))+((m2^2)*(1-m1))/(m2+(1-m1))
					##
					m1	<-a12; m2<-a3
					Pr	<-m1*m2+((m1^2)*(1-m2))/(m1+(1-m2))+((m2^2)*(1-m1))/(m2+(1-m1))
					return(Pr)
				}else if(sum((x)[]>0)[]==4){		
					a1	<-	x[which(x[]>0)[1]]; #	m1	
					a2	<-	x[which(x[]>0)[2]]; #	m2
					a3	<-	x[which(x[]>0)[3]]; #	
					a4	<-	x[which(x[]>0)[4]]; #
					##
					m1	<-a1; m2<-a2
					a12	<-m1*m2+((m1^2)*(1-m2))/(m1+(1-m2))+((m2^2)*(1-m1))/(m2+(1-m1))
					##
					m1	<-a12; m2<-a3
					a123<-m1*m2+((m1^2)*(1-m2))/(m1+(1-m2))+((m2^2)*(1-m1))/(m2+(1-m1))
					##
					m1	<-a123; m2<-a4
					Pr	<-m1*m2+((m1^2)*(1-m2))/(m1+(1-m2))+((m2^2)*(1-m1))/(m2+(1-m1))
					return(Pr)					
				}
				}))	
					multiSNV[pos.to.use,14]<-pos.qualities
				}else{}
				
		# PART 3			- SAVING
						
			# save:
			setwd(res.dir)
			save.name		<-paste(strsplit(gp.met.list[set.nr],split=".multiSNV")[[1]][1],"_NsFiltered.multiSNV",sep="")
			write.table(multiSNV, file=save.name,sep="\t",quote=FALSE,row.names=FALSE)
			cat("[ ... ] saved as:   ",save.name,"\n");cat("\n")			
		}# end	
		
		
		}else{
			cat("ERROR ??????  ","\n")
			cat("Noise threshold is equal to 0","\n") 		
			cat("no nucleotides will b removed","\n") 	
			cat("the program procedes, saving each file with new name","\n") 	
			cat("\n") 	

			set.nr<-1
			for(set.nr in 1:length(gp.met.list)){
				# load data
				setwd(res.dir)	
				multiSNV				<-as.matrix(read.table(gp.met.list[set.nr],header=TRUE,sep="\t",as.is=T,strip.white=T,colClasses="character"))
						
				# save:
				setwd(res.dir)
				save.name		<-paste(strsplit(gp.met.list[set.nr],split=".multiSNV")[[1]][1],"_NsFiltered.multiSNV",sep="")
				write.table(multiSNV, file=save.name,sep="\t",quote=FALSE,row.names=FALSE)
				cat("[ ... ] ","Ns.Tr == 0","\n")
				cat("[ ... ] saved as:   ",save.name,"\n");cat("\n")			
			}# end		
		}# e.d second if else
						
	}# end if else error message

}

































###	---------------------------------------------------------------------------------------------------
###	29. Remove.Noise.using.Acceptance.Threshold			6000
###	---------------------------------------------------------------------------------------------------	
Remove.Noise.using.Acceptance.Threshold	<-function(in.dir,res.dir,file.name,Acc){			
		
	cutoff.value		<-Acc
	max.cutoff		<-1# always
		
	# search for:
	setwd(in.dir); dir()
	search.for			<-"multiSNV"
	dir.list				<-dir()[which(grepl(search.for,dir())[]==TRUE)];dir.list	
	###
	search.for			<-file.name
	gp.met.list			<-dir.list[which(grepl(search.for,dir.list)[]==TRUE)]
	gp.met.list	
	
	
	# info:
	cat("***","-","----------------","-","***","\n") 
	cat("***","-","Noise filtration","-","***","\n") 
	cat("***","-",date(),"\n") 
	cat("***","-","----------------","-","\n") 
	cat("---","-","         Acc     =",cutoff.value,"\n")
	cat("   ","-","----------------","-","\n") 	
	cat("---","-","all nucleotides with this value or below will be removed","\n") 
	cat("---","-","from each consensus identified with file.name that you provided","\n")
	cat("   ","-","----------------","-","\n") 	
	cat("---","-","file identified with:",file.name,"\n") 
	cat("---","-","the number of files =",length(gp.met.list),"\n") 
	cat("---","-","           found in =",in.dir,"\n") 
	cat("   ","-","----------------","-","\n") 	
	cat("---","-","CAUTION","\n") 
	cat("---","-","setting up too low max.cutoff threshold can reduce the % of removed nucleotides","\n") 
	cat("---","-","sometimes many nucleotides can haave the same qualiy score","\n") 
	cat("---","-","than the percentage of removed positions can be higher","\n")
	cat("   ","-","----------------","-","\n") 
	cat("\n") ;cat("\n") ;cat("\n") 

	# big if.
	if(length(gp.met.list)[]==0){
		cat("ERROR","\n");cat("ERROR","\n");cat("ERROR","\n")
		cat("No files were found in in.dir","\n") 		
		cat("please control if the file.name is correct and or in.dir is correct and not empty","\n") 	
		cat("\n") 	
	}else{
	
	set.nr<-1
	for(set.nr in 1:length(gp.met.list)){
					
		# PART 1				- Noise filtration
						
		# load data
		setwd(res.dir)	
		prob.multi			<-as.matrix(read.table(gp.met.list[set.nr],header=TRUE,sep="\t",as.is=T,strip.white=T,colClasses="character"))
				
		# exctract		
		dim.with.snps		<-c(10:13);colnames(prob.multi)[dim.with.snps]
		prob.snps			<-prob.multi[,dim.with.snps]
		mode(prob.snps)		<-"numeric"
		
		# calculate cutoff value
		#ns.threshold			<-0
		all.nucl.probs		<-as.numeric(prob.snps)
		all.nucl.probs		<-all.nucl.probs[which(all.nucl.probs[]>0)]
		#pos.with.cutoff		<-round(length(all.nucl.probs)*ns.threshold,digits=0)	
		cutoff				<-as.numeric(cutoff.value)
		
			
		# check how many positions will be removed
			nr.of.removed.nucl				<-length(which(as.numeric(all.nucl.probs)[]<=cutoff))
			perc.of.removed.nucl				<-round(nr.of.removed.nucl/length(all.nucl.probs),digits=4)
			
		# info
			cat("-",set.nr,"-","loading:",gp.met.list[set.nr],"\n")	
			cat("-",set.nr,"-","the cutoff value =",cutoff,"\n")					
			cat("-",set.nr,"-","which removes    =",nr.of.removed.nucl,"out off",length(all.nucl.probs),"nucleotides at that consensus","\n")
			cat("-",set.nr,"-","                  ","i.e = ",perc.of.removed.nucl,"of nucleotides in this consensus","\n")	
		
		# removing these positions and rebuild consensus	
			prob.snps[which(prob.snps[]<=cutoff)]	<-0
			prob.multi[,c(10:13)]					<-prob.snps
			
						
		# PART 2				- positions with na identification and calulation of qaluty scores for the positon	
			multiSNV			<-prob.multi
								
			# now make allele presence in 		
				snps.n									<-c("A","C","G","T")
				snp.pres									<-apply(multiSNV[,c(10:13)],1,function(x){paste(snps.n[which(x[]>0)],collapse="/")})
				multiSNV[,15]							<-snp.pres
				multiSNV[which(multiSNV[,15]==""),15]	<-"na"
								
			# column for genotype quality (Pr(genotype))									
				empty.column										<-rep(1,nrow(multiSNV))	
				empty.column[which(multiSNV[,15]=="na")]				<-0
				multiSNV[,14]										<-empty.column							
				colnames(multiSNV)[14]								<-"M(genotype)"				
				multiSNV												<-as.matrix(multiSNV)				
			
			# first check if recaltulation is needed
				prob.snps[,]<-1	
				m.values.in.a.consensus	<-as.numeric(names(table(as.numeric(prob.snps)[which(as.numeric(prob.snps)[]>0)])))
				if(sum(m.values.in.a.consensus[]!=1)[]>0){use.pcr5.method<-"no"
				}else{				use.pcr5.method	<-"yes"}
																
			#yes 
			#	caluclate genotype quality for positions 
			#	(only for probabilistic methods)							
				if(use.pcr5.method[]=="yes"){	
					pos.to.use											<-which(empty.column[]>0)
					temp.mat											<-multiSNV[pos.to.use,c(10:13)]
					mode(temp.mat)										<-"numeric"
					
					## only for snps's, future 
					##	v2 will be adapted for InDels
					##	I had to finish my PhD and I didnt have enought time ti make it nicely :)				
					
					pos.qualities<-as.numeric(apply(temp.mat,1,function(x){	
				if(sum((x)[]>0)[]==1){	
					Pr	<-	x[which(x[]>0)[1]]; #	m1	
					return(Pr)
				}else if(sum((x)[]>0)[]==2){									
					a1	<-	x[which(x[]>0)[1]]; #	m1	
					a2	<-	x[which(x[]>0)[2]]; #	m2
					m1	<-a1; m2<-a2
					Pr	<-m1*m2+((m1^2)*(1-m2))/(m1+(1-m2))+((m2^2)*(1-m1))/(m2+(1-m1))
					return(Pr)
				}else if(sum((x)[]>0)[]==3){		
					a1	<-	x[which(x[]>0)[1]]; #	m1	
					a2	<-	x[which(x[]>0)[2]]; #	m2
					a3	<-	x[which(x[]>0)[3]]; #	
					##
					m1	<-a1; m2<-a2
					a12	<-m1*m2+((m1^2)*(1-m2))/(m1+(1-m2))+((m2^2)*(1-m1))/(m2+(1-m1))
					##
					m1	<-a12; m2<-a3
					Pr	<-m1*m2+((m1^2)*(1-m2))/(m1+(1-m2))+((m2^2)*(1-m1))/(m2+(1-m1))
					return(Pr)
				}else if(sum((x)[]>0)[]==4){		
					a1	<-	x[which(x[]>0)[1]]; #	m1	
					a2	<-	x[which(x[]>0)[2]]; #	m2
					a3	<-	x[which(x[]>0)[3]]; #	
					a4	<-	x[which(x[]>0)[4]]; #
					##
					m1	<-a1; m2<-a2
					a12	<-m1*m2+((m1^2)*(1-m2))/(m1+(1-m2))+((m2^2)*(1-m1))/(m2+(1-m1))
					##
					m1	<-a12; m2<-a3
					a123<-m1*m2+((m1^2)*(1-m2))/(m1+(1-m2))+((m2^2)*(1-m1))/(m2+(1-m1))
					##
					m1	<-a123; m2<-a4
					Pr	<-m1*m2+((m1^2)*(1-m2))/(m1+(1-m2))+((m2^2)*(1-m1))/(m2+(1-m1))
					return(Pr)					
				}
				}))	
					multiSNV[pos.to.use,14]<-pos.qualities
				}else{}
				
		# PART 3			- SAVING
						
			# save:
			setwd(res.dir)
			save.name		<-paste(strsplit(gp.met.list[set.nr],split=".multiSNV")[[1]][1],"_NsFiltered.multiSNV",sep="")
			write.table(multiSNV, file=save.name,sep="\t",quote=FALSE,row.names=FALSE)
			cat("[ ... ] saved as:   ",save.name,"\n");cat("\n")			
		}# end	
											
	}# end if else error message

}#end Remove.Noise.using.Cutoff.value


















###	---------------------------------------------------------------------------------------------------
###	30. fuse.consensus			6200
###	---------------------------------------------------------------------------------------------------	
fuse.consensus<-function(in.dir,res.dir,file.name,tr.m.name,pr.m.name){

	setwd(in.dir); 	#dir()
	dir.list			<-dir()[which(grepl(file.name,dir())[]==TRUE)]	
	tr.m.list		<-dir.list[which(grepl(tr.m.name,dir.list)[]==TRUE)]
	pr.m.list		<-dir.list[which(grepl(pr.m.name,dir.list)[]==TRUE)]	
	
	if(length(dir.list)[]==0){
		cat("ERROR","no files were found using - tr.m.name ","\n"); to.stop<-"yes"
		}else{	to.stop<-"no"}
	if(length(tr.m.list)[]==0){
		cat("ERROR","no files were found using - tr.m.name ","\n"); to.stop<-"yes"
		}else{	to.stop<-"no"}
	if(length(pr.m.list)[]==0){
		cat("ERROR","no files were found using - tr.m.name ","\n"); to.stop<-"yes"
		}else{	to.stop<-"no"}

	# big.if
	if(to.stop[]=="yes"){
	cat("the work couldnt be done - STOP")
	}else{		
	set.nr<-1	;tr.m.list[set.nr]
	for(set.nr in 1:length(tr.m.list)){
							
		# names
		setwd(in.dir)
		tr.m.to.load		<-tr.m.list[set.nr]
		search.name		<-strsplit(tr.m.list[set.nr],split="__")[[1]][1]		
		pr.m.to.load		<-pr.m.list[which(grepl(search.name,pr.m.list)[]==TRUE)]
		
		# small if:
		if(length(pr.m.to.load)[]>1){
			cat("[ .",set.nr,". ]","-","ERRROR !!!!!!","\n")
			cat("[ .",set.nr,". ]","-","ERRROR !!!!!!","this file:",tr.m.to.load	,"\n")
			cat("[ .",set.nr,". ]","-","ERRROR !!!!!!","has more than one target","\n")
			cat("[ .",set.nr,". ]","-","ERRROR !!!!!!","please check whether the proper file is in in.dir","\n")
			cat("[ .",set.nr,". ]","-","ERRROR !!!!!!","if you used __ in the name it can also cause problems","\n")
		}else if(length(pr.m.to.load)[]==0){
			cat("[ .",set.nr,". ]","-","ERRROR !!!!!!","\n")
			cat("[ .",set.nr,". ]","-","ERRROR !!!!!!","this file:",tr.m.to.load	,"\n")
			cat("[ .",set.nr,". ]","-","ERRROR !!!!!!","has more than NO target","\n")
			cat("[ .",set.nr,". ]","-","ERRROR !!!!!!","please check whether the other file is in in.dir","\n")
			cat("[ .",set.nr,". ]","-","ERRROR !!!!!!","if you used __ in the name it can also cause problems","\n")	
		}else{
			
			#PART 1 . masking
			
			# load data
			#.1.	
			setwd(in.dir)
			tr.multiSNV			<-as.matrix(read.table(tr.m.to.load,header=TRUE,sep="\t",as.is=T,strip.white=T,colClasses="character"))
			set.name				<-strsplit(tr.m.to.load	,split="__")[[1]][1]
			#.2.
			setwd(in.dir)
			pr.multiSNV			<-as.matrix(read.table(pr.m.to.load,header=TRUE,sep="\t",as.is=T,strip.white=T,colClasses="character"))
			set.name				<-strsplit(tr.m.to.load	,split="__")[[1]][1]			
			
			#info:
			
			cat("[ .",set.nr,". ]     fusing:   ",tr.m.to.load,"\n")		
			cat("[ .",set.nr,". ]  	  with:   ",pr.m.to.load,"\n")
			cat("[ .",set.nr,". ]  	      :   ",date(),"\n")	

			# new multi
			new.multiSNV			<-as.matrix(tr.multiSNV)
			
			# extract
			prob.snps.qual		<-as.matrix(pr.multiSNV[,c(10:13)])	
			mode(prob.snps.qual)<-"numeric"
			tr.snp.pres			<-as.matrix(tr.multiSNV[,c(10:13)])	
			mode(tr.snp.pres)	<-"numeric"
			
			# mask
			prob.snps.qual[which(tr.snp.pres[]==0)]<-0	
			new.multiSNV	[,c(10:13)]	<-prob.snps.qual	
						
			
			#PART 2			- positions with na identification and calulation of qaluty scores for the positon	
			multiSNV			<-new.multiSNV	# ido it because i dont want to cha ge my old script :)
			prob.snps		<-multiSNV			
			
			# now make allele presence in 		
				snps.n									<-c("A","C","G","T")
				snp.pres								<-apply(multiSNV[,c(10:13)],1,function(x){paste(snps.n[which(x[]>0)],collapse="/")})
				multiSNV[,15]							<-snp.pres
				multiSNV[which(multiSNV[,15]==""),15]	<-"na"
								
			# column for genotype quality (Pr(genotype))									
				empty.column										<-rep(1,nrow(multiSNV))	
				empty.column[which(multiSNV[,15]=="na")]				<-0
				multiSNV[,14]										<-empty.column							
				colnames(multiSNV)[14]								<-"M(genotype)"				
				multiSNV												<-as.matrix(multiSNV)				
			
			# first check if recaltulation is needed
				prob.snps[,]<-1	
				m.values.in.a.consensus	<-as.numeric(names(table(as.numeric(prob.snps)[which(as.numeric(prob.snps)[]>0)])))
				if(sum(m.values.in.a.consensus[]!=1)[]>0){use.pcr5.method<-"no"
				}else{				use.pcr5.method	<-"yes"}
																
			#yes 
			#	caluclate genotype quality for positions 
			#	(only for probabilistic methods)							
				if(use.pcr5.method[]=="yes"){	
					pos.to.use											<-which(empty.column[]>0)
					temp.mat											<-multiSNV[pos.to.use,c(10:13)]
					mode(temp.mat)										<-"numeric"
					
					## only for snps's, future 
					##	v2 will be adapted for InDels
					##	I had to finish my PhD and I didnt have enought time ti make it nicely :)				
					
					pos.qualities<-as.numeric(apply(temp.mat,1,function(x){	
				if(sum((x)[]>0)[]==1){	
					Pr	<-	x[which(x[]>0)[1]]; #	m1	
					return(Pr)
				}else if(sum((x)[]>0)[]==2){									
					a1	<-	x[which(x[]>0)[1]]; #	m1	
					a2	<-	x[which(x[]>0)[2]]; #	m2
					m1	<-a1; m2<-a2
					Pr	<-m1*m2+((m1^2)*(1-m2))/(m1+(1-m2))+((m2^2)*(1-m1))/(m2+(1-m1))
					return(Pr)
				}else if(sum((x)[]>0)[]==3){		
					a1	<-	x[which(x[]>0)[1]]; #	m1	
					a2	<-	x[which(x[]>0)[2]]; #	m2
					a3	<-	x[which(x[]>0)[3]]; #	
					##
					m1	<-a1; m2<-a2
					a12	<-m1*m2+((m1^2)*(1-m2))/(m1+(1-m2))+((m2^2)*(1-m1))/(m2+(1-m1))
					##
					m1	<-a12; m2<-a3
					Pr	<-m1*m2+((m1^2)*(1-m2))/(m1+(1-m2))+((m2^2)*(1-m1))/(m2+(1-m1))
					return(Pr)
				}else if(sum((x)[]>0)[]==4){		
					a1	<-	x[which(x[]>0)[1]]; #	m1	
					a2	<-	x[which(x[]>0)[2]]; #	m2
					a3	<-	x[which(x[]>0)[3]]; #	
					a4	<-	x[which(x[]>0)[4]]; #
					##
					m1	<-a1; m2<-a2
					a12	<-m1*m2+((m1^2)*(1-m2))/(m1+(1-m2))+((m2^2)*(1-m1))/(m2+(1-m1))
					##
					m1	<-a12; m2<-a3
					a123<-m1*m2+((m1^2)*(1-m2))/(m1+(1-m2))+((m2^2)*(1-m1))/(m2+(1-m1))
					##
					m1	<-a123; m2<-a4
					Pr	<-m1*m2+((m1^2)*(1-m2))/(m1+(1-m2))+((m2^2)*(1-m1))/(m2+(1-m1))
					return(Pr)					
				}
				}))	
					multiSNV[pos.to.use,14]<-pos.qualities
				}else{}
				
			#PART 3		- SAVING						
			
			#probability.method.name
			if(grepl("Probabilistic_beta_method",pr.m.to.load)){
						probability.method.name	<-"ProbBeta_QualityScores"
				}else{	probability.method.name	<-"Prob_QualityScores"}
						
			# save:
			setwd(res.dir)
			save.name		<-paste(strsplit(tr.m.to.load,split=".multiSNV")[[1]][1],"_with_",probability.method.name,".multiSNV",sep="")
			write.table(multiSNV, file=save.name,sep="\t",quote=FALSE,row.names=FALSE)
			cat("[ .",set.nr,". ]  saved as:   ",save.name,"\n");cat("\n")			
		
		
		}# end small if. else	
		} # end set.nr
	}#	# end big if else
}# end function
		
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
###	---------------------------------------------------------------------------------------------------
###	31. Compare.samples.using.multiSNV.files		6400
###	---------------------------------------------------------------------------------------------------			
Compare.samples.using.multiSNV.files		<-function(in.dir,res.dir,new.file.name,sample.list,remove.na){
	
	to.stop<-"no"
	if(length(sample.list)[]<2){
		cat("ERROR","you need at least two samples to perfomr a comparison","\n")
		to.stop<-"yes"
		}else{}		
		
	if(length(table(sample.list	))<length(sample.list)){
		cat("ERROR","You have duplicates in sample list, i.e. the same sample was put more than once","\n")
		to.stop<-"yes"
		}else{}
		
		#STEP 1
		# load each file and extract basics
			
			# list of elements:
			multiSNV.LIST			<-as.list(sample.list)
			names(multiSNV.LIST)		<-sample.list
			
			s.nr<-1
			for(s.nr in 1:length(sample.list)){				
				setwd(in.dir)
				multiSNV		<-as.matrix(read.table(sample.list[s.nr],header=TRUE,sep="\t",as.is=T,strip.white=T,colClasses="character"))
				###
				if(s.nr[]==1){
					first.columns		<-multiSNV[,c(1:9),drop=FALSE]
					first.columns[,1]	<-new.file.name
				}else{}
				multiSNV.LIST[[s.nr]]	<-multiSNV[,]
			}#end s.nr
				
		#STEP 2
		# find which positions to compare				
			take.col		<-14; colnames(multiSNV.LIST[[1]])[take.col]; colnames(multiSNV.LIST[[1]])
			temp								<-multiSNV.LIST[[1]][,take.col]
			for(i in 2:length(multiSNV.LIST)){temp<-cbind(temp,multiSNV.LIST[[i]][,take.col,drop=FALSE])}
			colnames(temp)					<-names(multiSNV.LIST)
			pos.qual.in.each.sample			<-temp					
			mode(pos.qual.in.each.sample	)	<-"numeric"
			BINARY_pos.qual.in.each.sample	<-pos.qual.in.each.sample
			BINARY_pos.qual.in.each.sample[which(pos.qual.in.each.sample[,]>0)]<-1
			noNA.pos							<-rowSums(BINARY_pos.qual.in.each.sample)
				
			
			pos.to.compare					<-rep(1,nrow(BINARY_pos.qual.in.each.sample))			
			if(remove.na[]=="no"){
			pos.to.compare[which(noNA.pos[]<=1)]<-0# if no only positions without data or with only one sample with data are removed
			}else{pos.to.compare[which(noNA.pos[]<length(multiSNV.LIST))]<-0}	
			###
			# one is also not acceptable because it is an informtaion in only one consensus
						
			
			cat("-----------------------------------------------------------------","\n")	
			cat("[ .-. ]     file name is",new.file.name,"\n");cat("\n")		
			cat("[ .-. ]      it contains",length(multiSNV.LIST),"amples that must be comapred: ","\n")		
			for(i in 1:length(multiSNV.LIST)){				
			cat("[ ... ]                 ",i,names(multiSNV.LIST)[i],"\n")}
			cat("[ ... ]         strated at:",date(),"\n")
			cat("-----------------------------------------------------------------","\n")					
			cat("[ ... ]            STATS","\n")
			cat("[ ... ]     samples.nr       =",length(multiSNV.LIST),"\n")
			cat("[ ... ]   total.pos.nr       =",length(pos.to.compare),"\n")
			cat("[ ... ] nr of pos.to compare =",sum(pos.to.compare),"\n")
			cat("[ ... ]","\n")
			cat("[ ... ] pos with NA were removed =",remove.na,"\n")
			cat("[ ... ]","\n")
			cat("[ ... ] Comment:","\n")
			cat("[ ... ] NA is a position that has no data in an entire consensus","\n")
			cat("[ ... ]","\n")
			
		#STEP 3		
		# 	SNP.combined.quality_mat	
			SNP.combined.quality_mat		<-multiSNV.LIST[[1]][,c(10:14)]
			SNP.combined.quality_mat[,]	<-0
		
		snp.nr<-2	
		for(snp.nr in 1:5){
			
			columns.with.snps		<-10:14
			take.col					<-columns.with.snps[snp.nr]; colnames(multiSNV.LIST[[1]])[columns.with.snps]
			temp						<-multiSNV.LIST[[1]][,take.col]
			for(i in 2:length(multiSNV.LIST)){temp<-cbind(temp,multiSNV.LIST[[i]][,take.col])}
			colnames(temp)			<-names(multiSNV.LIST)
			snp.quality				<-temp	
	
			#[first pair]
				temp.mat							<-snp.quality[,c(1:2)]			
				mode(temp.mat)						<-"numeric"				
				combined.quality.pos<-as.numeric(apply(temp.mat,1,function(x){	
				if(sum((x)[]>0)[]==0){
					Pr	<-	0
					return(Pr)					
				}else if(sum((x)[]>0)[]==1){	
					Pr	<-	x[which(x[]>0)[1]]; #	m1	
					return(Pr)
				}else {									
					a1	<-	x[which(x[]>0)[1]]; #	m1	
					a2	<-	x[which(x[]>0)[2]]; #	m2
					m1	<-a1; m2<-a2
					Pr	<-m1*m2+((m1^2)*(1-m2))/(m1+(1-m2))+((m2^2)*(1-m1))/(m2+(1-m1))
					return(Pr)			
				}
				}))							
			#[the rest]					
			if(length(multiSNV.LIST)[]>2){
			for(sample.nr in 1:length(multiSNV.LIST)){
				temp.mat							<-cbind(combined.quality.pos,snp.quality[,sample.nr])						
				mode(temp.mat)						<-"numeric"
				combined.quality.pos				<-as.numeric(apply(temp.mat,1,function(x){	
				if(sum((x)[]>0)[]==0){
					Pr	<-	0
					return(Pr)					
				}else if(sum((x)[]>0)[]==1){	
					Pr	<-	x[which(x[]>0)[1]]; #	m1	
					return(Pr)
				}else {									
					a1	<-	x[which(x[]>0)[1]]; #	m1	
					a2	<-	x[which(x[]>0)[2]]; #	m2
					m1	<-a1; m2<-a2
					Pr	<-m1*m2+((m1^2)*(1-m2))/(m1+(1-m2))+((m2^2)*(1-m1))/(m2+(1-m1))
					return(Pr)			
				}
				}))					
			}# end for(sample.nr	
			}else{}#end if else >2 samples	
												
		SNP.combined.quality_mat	[,snp.nr]<-combined.quality.pos
		}#end		
				
		# remove postions with na - if necessary
		SNP.combined.quality_mat[pos.to.compare[]==0,]<-0
		#	sum(as.numeric(SNP.combined.quality_mat[,5])[]>0)
		
				
		#STEP 4		
		shared.SNV						<-rep("-",nrow(SNP.combined.quality_mat))
		variable.SNV						<-rep("-",nrow(SNP.combined.quality_mat))
		SNPs								<-c("A","C","G","T")				
		
		# 	a) pos.with.no.change
			pos.with.no.change			<-as.numeric(apply(SNP.combined.quality_mat[,c(1:4)],1,function(x){sum(x[]>0)})[]<2)
			temp.mat						<-SNP.combined.quality_mat[which(pos.with.no.change[]==1),c(1:4)]
			shared.SNV[which(pos.with.no.change[]==1)]<-	apply(temp.mat,1,function(x){
				if(sum(x[]>0)[]==0){return("-")
				}else{	return(SNPs[which(x[]>0)])}	
				})

		if(length(which(pos.with.no.change[]==0))[]==0){	
		}else{	

		#	b) binary list
			binary.snp_LIST			<-as.list(names(multiSNV.LIST))
			names(binary.snp_LIST)	<-names(multiSNV.LIST)
			for(i in 1:length(multiSNV.LIST)){
				take.col					<-c(10:13); colnames(multiSNV.LIST[[1]])[take.col]
				temp						<-multiSNV.LIST[[i]][which(pos.with.no.change[]==0),take.col,drop=FALSE]
				mode(temp)				<-"numeric"
				temp[which(temp[,]>0)]	<-1
				binary.snp_LIST[[i]]	<-temp
			}# it can be still position with no chnage , i.e. more than one allele
				
		#	c) dissect each position
			short.no.change.snp			<-rep("-",length(which(pos.with.no.change[]==0)))
			short.variable.snp			<-rep("-",length(which(pos.with.no.change[]==0)))
		
			locus.nr<-1
			for(locus.nr in 1:nrow(binary.snp_LIST[[1]])){			
				
				# one locus mat
				one.locus			<-binary.snp_LIST[[1]][locus.nr,,drop=FALSE]
				for(s.nr in 2:length(binary.snp_LIST)){
					one.locus		<-rbind(one.locus,binary.snp_LIST[[s.nr]][locus.nr,,drop=FALSE])}; 
				mode(one.locus)		<-"numeric"
				rownames(one.locus)	<-names(binary.snp_LIST)
				one.locus
								
				# only missiing data
				samples.with.only.missing.data	<-as.numeric(rowSums(one.locus)[]==0)				
				one.locus						<-one.locus[which(samples.with.only.missing.data[]==0),,drop=FALSE]
				one.locus
				
				# check.whether.each.allele.is.in.each.remaining.sample
				snp.counts			<-colSums(one.locus)
				
				# shared.snps/common.snps
				shared.snps			<-SNPs[which(snp.counts[]==nrow(one.locus))];shared.snps	
				if(length(shared.snps)[]==0){
					}else{
					short.no.change.snp[locus.nr]<-paste(shared.snps,collapse="/")}		
				
				# variable snps			
				variable.snps			<-SNPs[which(snp.counts[]<nrow(one.locus)&snp.counts[]>0)];variable.snps
				if(length(variable.snps)[]==0){
					}else{
					short.variable.snp[locus.nr]<-paste(variable.snps,collapse="/")}						
				}
					
		#	d) reintroduce these position			
			shared.SNV[which(pos.with.no.change[]==0)]		<-short.no.change.snp	
			variable.SNV	[which(pos.with.no.change[]==0)]	<-short.variable.snp
					
		}# end if else for variable positions
				
		#STEP 5
		#		info.about.each.sample
			take.col					<-16; colnames(multiSNV.LIST[[1]])[take.col]
			temp						<-multiSNV.LIST[[1]][,take.col]
			for(i in 2:length(multiSNV.LIST)){temp<-cbind(temp,multiSNV.LIST[[i]][,take.col])}
			colnames(temp)			<-names(multiSNV.LIST)
			info.about.each.sample	<-temp	
						
						
		#STEP 6
		#		info.about.each.sample
			take.col					<-19; colnames(multiSNV.LIST[[1]])[take.col]; colnames(multiSNV.LIST[[1]])
			temp						<-multiSNV.LIST[[1]][,take.col]
			for(i in 2:length(multiSNV.LIST)){temp<-cbind(temp,multiSNV.LIST[[i]][,take.col])}
			colnames(temp)			<-names(multiSNV.LIST)
			indel.presence.in.each.sample	<-temp	
			
			indel.presence.in.each.sample[which(indel.presence.in.each.sample[,]=="n")]<-0			
			indel.presence.in.each.sample[which(indel.presence.in.each.sample[,]=="y")]<-1
			mode(indel.presence.in.each.sample)<-"numeric"
			InDel			<-as.numeric(rowSums(indel.presence.in.each.sample)[]>0)
			InDel[which(InDel[]==1)]<-"y"
			InDel[which(InDel[]==0)]<-"n"
			

		#STEP 7		
		# 	build a CombinedSNV file			
			SNP.combined.quality_mat[which(pos.to.compare[]==0),]	<-0
			shared.SNV[which(pos.to.compare[]==0)]					<-"-"
			variable.SNV	[which(pos.to.compare[]==0)]					<-"-"


			
		#STEP 8
		#	collect info about variable and stable postions
		
		#[1]
		snp.qual.ft			<-SNP.combined.quality_mat
		snp.qual.ft[,5]		<-shared.SNV
		shared.SNV.Quality	<-apply(snp.qual.ft,1,function(x){
		snps<-c("A","C","G","T");
		if(x[5]=="-"	){
			pos.snp.qual	<-0
		}else{
			snp.to.take	<-strsplit(x[5],split="/")[[1]];			snp.to.take
			for(i in 1:length(snp.to.take)){snp.to.take[i]<-which(snps[]==snp.to.take[i])}
			mode(snp.to.take)<-"numeric"	
			if(length(snp.to.take)[]==1){
				pos.snp.qual		<-as.numeric(x[snp.to.take])			
			}else if(length(snp.to.take)[]==2){							
					all.q			<-as.numeric(x[snp.to.take])
					a1	<-all.q[1]	
					a2	<-all.q[2]
					m1	<-a1; m2<-a2
					pos.snp.qual		<-m1*m2+((m1^2)*(1-m2))/(m1+(1-m2))+((m2^2)*(1-m1))/(m2+(1-m1))			
			}else if(length(snp.to.take)[]==3){	
					all.q			<-as.numeric(x[snp.to.take])
					a1	<-all.q[1]	
					a2	<-all.q[2]
					a3	<-all.q[3]
					##
					m1	<-a1; m2<-a2
					a12	<-m1*m2+((m1^2)*(1-m2))/(m1+(1-m2))+((m2^2)*(1-m1))/(m2+(1-m1))
					##
					m1	<-a12; m2<-a3
					pos.snp.qual		<-m1*m2+((m1^2)*(1-m2))/(m1+(1-m2))+((m2^2)*(1-m1))/(m2+(1-m1))
			}else if(length(snp.to.take)[]==4){	
					all.q			<-as.numeric(x[snp.to.take])
					a1	<-all.q[1]	
					a2	<-all.q[2]
					a3	<-all.q[3]
					a4	<-all.q[4]
					##
					m1	<-a1; m2<-a2
					a12	<-m1*m2+((m1^2)*(1-m2))/(m1+(1-m2))+((m2^2)*(1-m1))/(m2+(1-m1))
					##
					m1	<-a12; m2<-a3
					a123<-m1*m2+((m1^2)*(1-m2))/(m1+(1-m2))+((m2^2)*(1-m1))/(m2+(1-m1))
					##
					m1	<-a123; m2<-a4
					pos.snp.qual		<-m1*m2+((m1^2)*(1-m2))/(m1+(1-m2))+((m2^2)*(1-m1))/(m2+(1-m1))			
			}
		}	
		return(pos.snp.qual)
		})

		#[2]
		snp.qual.ft			<-SNP.combined.quality_mat
		snp.qual.ft[,5]		<-variable.SNV
		variable.SNV.Quality	<-apply(snp.qual.ft,1,function(x){
		snps<-c("A","C","G","T");
		if(x[5]=="-"	){
			pos.snp.qual	<-0
		}else{
			snp.to.take	<-strsplit(x[5],split="/")[[1]];			snp.to.take
			for(i in 1:length(snp.to.take)){snp.to.take[i]<-which(snps[]==snp.to.take[i])}
			mode(snp.to.take)<-"numeric"	
			if(length(snp.to.take)[]==1){
				pos.snp.qual		<-as.numeric(x[snp.to.take])			
			}else if(length(snp.to.take)[]==2){							
					all.q			<-as.numeric(x[snp.to.take])
					a1	<-all.q[1]	
					a2	<-all.q[2]
					m1	<-a1; m2<-a2
					pos.snp.qual		<-m1*m2+((m1^2)*(1-m2))/(m1+(1-m2))+((m2^2)*(1-m1))/(m2+(1-m1))			
			}else if(length(snp.to.take)[]==3){	
					all.q			<-as.numeric(x[snp.to.take])
					a1	<-all.q[1]	
					a2	<-all.q[2]
					a3	<-all.q[3]
					##
					m1	<-a1; m2<-a2
					a12	<-m1*m2+((m1^2)*(1-m2))/(m1+(1-m2))+((m2^2)*(1-m1))/(m2+(1-m1))
					##
					m1	<-a12; m2<-a3
					pos.snp.qual		<-m1*m2+((m1^2)*(1-m2))/(m1+(1-m2))+((m2^2)*(1-m1))/(m2+(1-m1))
			}else if(length(snp.to.take)[]==4){	
					all.q			<-as.numeric(x[snp.to.take])
					a1	<-all.q[1]	
					a2	<-all.q[2]
					a3	<-all.q[3]
					a4	<-all.q[4]
					##
					m1	<-a1; m2<-a2
					a12	<-m1*m2+((m1^2)*(1-m2))/(m1+(1-m2))+((m2^2)*(1-m1))/(m2+(1-m1))
					##
					m1	<-a12; m2<-a3
					a123<-m1*m2+((m1^2)*(1-m2))/(m1+(1-m2))+((m2^2)*(1-m1))/(m2+(1-m1))
					##
					m1	<-a123; m2<-a4
					pos.snp.qual		<-m1*m2+((m1^2)*(1-m2))/(m1+(1-m2))+((m2^2)*(1-m1))/(m2+(1-m1))			
			}
		}	
		return(pos.snp.qual)
		})		
		
						
		#STEP 9	
		# 	build a CombinedSNV file			
			CombinedSNV	<-cbind(first.columns,
								SNP.combined.quality_mat[,5,drop=FALSE],
								shared.SNV.Quality,
								variable.SNV.Quality,
								shared.SNV,
								variable.SNV,
								info.about.each.sample,
								InDel)		
		#STEP 10	
		# 	saving
			cat("-----------------------------------------------------------------","\n")	
			cat("-----------------------------------------------------------------","\n")	
			setwd(res.dir)
			save.name		<-paste(new.file.name,".CombinedSNV",sep=""); 	save.name
			write.table(CombinedSNV, file=save.name,sep="\t",quote=FALSE,row.names=FALSE)
			cat("[ .-. ] saved as:   ",save.name,"\n");cat("\n")	
			cat("[ ... ]       in:   ",res.dir,"\n")
			cat("[ ... ]       at:   ",date(),"\n")
			cat("-----------------------------------------------------------------","\n")		
		
		
		# test:
		#	tm			<-CombinedSNV[,c(10:18)]
		#	ql			<-as.numeric(tm[,5])
		#	names(ql)	<-1:length(ql)
		#	new.order	<-as.numeric(names(sort(ql,decreasing=TRUE)))		
		#	colnames(tm)<-1:ncol(tm)
		#	n.tm		<-tm[new.order,]	
		#	n.tm[which(n.tm[,7]!="-"),]		
		
	}# end function	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
###	---------------------------------------------------------------------------------------------------
###	32. Calculate.Acc.using.Ns	6800
###	---------------------------------------------------------------------------------------------------	
Calculate.Acc.using.Ns<-function(in.dir,file.name,Ns,per.replicate,max.cutoff){			
		
	
	Ns.Tr<-Ns
	
	# search for:
	setwd(in.dir); dir()
	search.for			<-"multiSNV"
	dir.list				<-dir()[which(grepl(search.for,dir())[]==TRUE)];dir.list	
	###
	search.for			<-file.name
	gp.met.list			<-dir.list[which(grepl(search.for,dir.list)[]==TRUE)]
	gp.met.list	
	
	# just in case
	if(max.cutoff[]>=1){
		max.cutoff<-1
		}else{}
	
	# info:
	cat("***","-","----------------","-","***","\n") 
	cat("***","-","Noise filtration","-","***","\n") 
	cat("***","-",date(),"\n") 
	cat("***","-","----------------","-","\n") 
	cat("---","-","      Ns            =",Ns.Tr,"\n")
	cat("---","-","      per replicate =",per.replicate,"\n")
	cat("---","-","      max cutoff    =",max.cutoff,"\n")
	cat("   ","-","----------------","-","\n") 	
	cat("---","-","file identified with:",file.name,"\n") 
	cat("---","-","the number of files =",length(gp.met.list),"\n") 
	cat("---","-","           found in =",in.dir,"\n") 
	cat("   ","-","----------------","-","\n") 
	cat("\n") ;cat("\n") ;cat("\n") 


	# big if.
	if(length(gp.met.list)[]==0){
		cat("ERROR","\n");cat("ERROR","\n");cat("ERROR","\n")
		cat("No files were found in in.dir","\n") 		
		cat("please control if the file.name is correct and or in.dir is correct and not empty","\n") 	
		cat("\n") 	
	}else{
	
	#second if. 
	if(Ns.Tr[]>0){
				
	set.nr<-1		
	Result.list					<-matrix(nrow=length(gp.met.list)+1,ncol=4)
	colnames(Result.list	)		<-c("tot.nr.of.nucl","removed.nucl","%.of.removed.nucl","Acc")
	rownames(Result.list	)		<-c(gp.met.list,"MEAN")
	Result.list[,]				<-0
	mode(Result.list	)			<-"numeric"
	
		
	for(set.nr in 1:length(gp.met.list)){
					
		# PART 1				- Noise filtration
						
		# load data
		setwd(res.dir)	
		prob.multi			<-as.matrix(read.table(gp.met.list[set.nr],header=TRUE,sep="\t",as.is=T,strip.white=T,colClasses="character"))

		# Ns.tr
		repl.nr				<-as.numeric(prob.multi[1,ncol(prob.multi)])
		if(per.replicate[]=="yes"){
			ns.threshold		<-Ns.Tr*repl.nr
			if(ns.threshold[]>=1){ns.threshold[1]<-1}else{}
		}
				
		# exctract		
		dim.with.snps		<-c(10:13);colnames(prob.multi)[dim.with.snps]
		prob.snps			<-prob.multi[,dim.with.snps]
		mode(prob.snps)		<-"numeric"
		
		# calculate cutoff value
		all.nucl.probs		<-as.numeric(prob.snps)
		all.nucl.probs		<-all.nucl.probs[which(all.nucl.probs[]>0)]
		pos.with.cutoff		<-round(length(all.nucl.probs)*ns.threshold,digits=0)	
		cutoff				<-as.numeric(sort(all.nucl.probs,decreasing = FALSE)[pos.with.cutoff])
		
		# check for cutoff
		if(cutoff[]>=as.numeric(max.cutoff)){
			cutoff<-max.cutoff; cutoff
			}else{}
			
		# check how many positions will be removed
			nr.of.removed.nucl				<-length(which(as.numeric(all.nucl.probs)[]<=cutoff))
			perc.of.removed.nucl				<-round(nr.of.removed.nucl/length(all.nucl.probs),digits=4)
			
		# fill in:
			Result.list[set.nr,1]		<-length(all.nucl.probs)
			Result.list[set.nr,2]		<-nr.of.removed.nucl	
			Result.list[set.nr,3]		<-perc.of.removed.nucl	
			Result.list[set.nr,4]		<-cutoff
		
		}# end	
		Result.list[nrow(Result.list),]<-apply(Result.list[c(1:c(nrow(Result.list)-1)),,drop=FALSE],2,function(x){mean(x)})		
		return(Result.list)
	}else{
	cat("ERROR ??????  ","\n")
	cat("Noise threshold is equal to 0","\n") 		
	cat("no nucleotides will b removed","\n") 	
	cat("the program procedes, saving each file with new name","\n") 	
	cat("\n") 	
	}# end		
	}# e.d second if else					

}
	
	
		
