#	OBJECTS

#	METHODS

setGeneric("dust",
			function(seqs, k=7L, e=3L, shift=1L) standardGeneric("dust"),
			signature=c("seqs"))
setMethod("dust", signature("DNAStringSet"),
	function(seqs, k=7L, e=3L, shift=1L) {
		gNames <- names(seqs)
		dna_strs <- c()
		ne <- k-e+1
		Ne <- 4**e
		TeE_ne.Ne <- as.integer(ne/Ne)+shift
		print(TeE_ne.Ne)
		pb <- txtProgressBar(0, length(seqs), style=3)
		for (i in 1:length(seqs)) {	#F1
	
			dna <- toString(seqs[[i]])
			dna <- strsplit(dna, split="")[[1]]
			tmp_dna <- c()
			
			for (j in seq(1, length(dna), k)) {	#F2
				kmer <- dna[j:(j+k-1)]
				kmer <- kmer[!is.na(kmer)]
				elemList <- list()
				
				if (length(kmer) > e){
					lenKmer <- length(kmer)-e
				}
				else{
					lenKmer <- length(kmer)
				}

				for (m in 1:lenKmer) {					#F3
					lmer <- paste(kmer[m:(m+e-1)], collapse="")
					if(lmer %in% names(elemList)) {
						elemList[[lmer]] <- elemList[[lmer]] + 1
					}
					else {
						elemList[[lmer]] <- 1
					}
				}												#EF3
				elemList <- unlist(elemList)
				if (sum(elemList > TeE_ne.Ne) > 0) {
					tmp_dna <- c(tmp_dna, rep("N", k))
				}
				else {
					tmp_dna <- c(tmp_dna, kmer)
				}
			}											#EF2
			dna_strs <- c(dna_strs, paste(tmp_dna, collapse=""))
			setTxtProgressBar(pb, i)
		}							#EF1
		names(dna_strs) <- gNames
		dna_result <- DNAStringSet(dna_strs, use.names=T)
		dna_result
})

setMethod("dust", signature("AAStringSet"),
	function(seqs) {

})




setGeneric("purge",
			function(seqs) standardGeneric("purge"),
			signature=c("seqs"))
setMethod("purge", signature("DNAStringSet"),
	function(seqs) {
		
})
setMethod("purge", signature("AAStringSet"),
	function(seqs) {

})
