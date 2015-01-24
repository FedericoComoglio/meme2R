#Federico Comoglio, ETH Zurich, 17 Feb 2014

library(XML)

extractPWM <- function( rootM, motifID ) {
	prob <- rootM[[ motifID ]][[ 'probabilities' ]]
	n <- length( prob[[ 1 ]] )
    	PWM <- matrix( Inf, nrow = 4, ncol = n )
	rownames( PWM ) <- c( 'A', 'C', 'G', 'T' )

	for ( i in seq_len( n ) ) {
        	for(j in seq_len( 4 ) ) {
            		PWM[ j, i ] <- as.numeric( xmlValue( prob[[ 1 ]][[ i ]][[ j ]]) )
        	}
    	}
	return( PWM  )
}

ConsensusFromPWM <- function( PWM ) {
	  nucleotides <- c('A', 'C', 'G', 'T')
	  pos <- apply( PWM, 2, which.max )
	  paste( nucleotides[ pos ], collapse = '' )
}

meme2R <- function( filename ) {   	
	file <- xmlTreeParse( file = filename, getDTD = FALSE )
    	root <- xmlRoot( file )
	rootM <- root[[ 'motifs' ]]
    	n <- length( rootM )
	
	attrs <- sapply( 1 : n, function( i ) xmlAttrs( rootM[[ i ]] ) )
	regExp <- sapply( 1 : n, function( i ) xmlValue( rootM[[ i ]][[ 'regular_expression' ]] ) )
	consensus <- sapply( 1 : n, function( i ) ConsensusFromPWM( extractPWM( rootM, motifID = i ) ) )

	motifs <- data.frame( t( attrs ), regExp, consensus )	
	colnames(motifs) <- c('id', 'name', 'width', 'sites', 
			      'ic', 're', 'llr',
			      'e_value', 'bayes_threshold',
		              'elapsed_time', 'regexp', 'consensus')
	
	return(motifs)
}
