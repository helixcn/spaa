sp.pair <- function(matr){
      if (any(is.na(matr))) {
         matr <- na.omit(matr)
         cat("Rows containing NA have been removed.\n")
      }
      
      pearson  <- cor(matr, method = "pearson")
      spearman <- cor(matr, method = "spearman")
      
      N <- nrow(matr)     # number of plots
      matr[matr > 1] <- 1 # presence-absence matrix
      
      empty_mat <- matrix(data = rep(NA, ncol(matr)*ncol(matr)), 
                          nrow = ncol(matr), 
                          ncol = ncol(matr), 
                          dimnames = list(colnames(matr), colnames(matr)))
      colnames(empty_mat) <- colnames(matr)
      rownames(empty_mat) <- colnames(matr)    
      
      ##Yate's correction for chisq
      chisq   <- empty_mat
      V       <- empty_mat
      Ochiai  <- empty_mat
      Dice    <- empty_mat
      Jaccard <- empty_mat
      PCC     <- empty_mat
      dd      <- empty_mat
      AC      <- empty_mat
      
      for (i in 1:ncol(matr)) {
         if(i < ncol(matr)){
            for (j in (i+1):ncol(matr)){
               temp_occurrence_sp1 <- matr[,i]
               temp_occurrence_sp2 <- matr[,j]
               temp_contigent_tab  <- table(temp_occurrence_sp1, temp_occurrence_sp2)
               
               a <- temp_contigent_tab[2,2] # number of plots in which both sp1 and sp2 are present
               b <- temp_contigent_tab[1,2] # number of plots in which sp1 is present, sp2 is absent
               c <- temp_contigent_tab[2,1] # number of plots in which sp1 is absent, sp2 is present
               d <- temp_contigent_tab[1,1] # number of plots in which both sp1 and sp2 are absent
               
               # Thanks for Ms Xueni Zhang for pointing out the error
               ### j is the row. i is the column
               ### always use the lower matrix
               ### j is dependent to i, j is always greater than i
               
               # V>0,positive, V<0,negative 
               V[j,i]       <-  ((a+d)-(b+c))/(a+b+c +d)
               Ochiai[j,i]  <-  a/sqrt((a+b)*(a+c))      #Ochiai index
               Dice[j,i]    <-  2*a/(2*a+b+c)            #Dice index
               Jaccard[j,i] <-  a/(a+b+c)                #Jaccard index
               PCC[j,i]     <-  (a*d-b*c)/((a+b)*(a+c)*(c+d)*(b+d)) ##Percentage cooccurance
               dd[j,i]      <-  a*d-b*c
               
               # Chi square
               chisq[j,i] <- ((((abs(a*d-b*c)-0.5*N) ^ 2)*N)/
                                  ((a+b)*(a+c)*(b+d)*(c+d)))
               
               # the Association index AC
               # -Wang Bosun, Peng Shaolin. Studies on the Measuring Techniques 
               # of Interspecific Association of Lower-Subtropical Evergreen-
               # Broadleaved Forests. I. The Exploration and the Revision 
               # on the Measuring Formulas of Interspecific Association[J]. 
               # Chin J Plan Ecolo, 1985, 9(4): 274-279.
               # -Hurlbert, S. H. (1969). A coefficient of interspecific association. 
               # Ecology, 50(1), 1-9.
               
               if (a*d >= b*c) {
                  AC[j,i] <- (a*d-b*c) /((a+b)*(b+d))
               } else { # a*d < b*c
                  if( a <= d ){
                     AC[j,i] <- (a*d-b*c) /((a+b)*(a+c))
                  } else { # a > d
                     AC[j,i] <- (a*d-b*c) /((b+d)*(c+d))
                  }
                  
               }
            }
         }
      }
      result <-
         list(
            chisq    = as.dist(chisq),
            chisqass = as.dist(dd),
            V        = as.dist(V),
            Ochiai   = as.dist(Ochiai),
            Dice     = as.dist(Dice),
            Jaccard  = as.dist(Jaccard),
            Pearson  = pearson,
            Spearman = spearman,
            PCC      = as.dist(PCC),
            AC       = as.dist(AC)
         )
      return(result)
}
