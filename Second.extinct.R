#################
#second.extinctM#
#################


#Description: This function is an update to the second.extinct function in the bipartite package to add a method for removing species according to its topological role in the network based on Modularity analysis. The method removes first network hubs, next module hubs, next connectors and finnally peripherial species. Species are removed randomly withinn a role.


#Load 2 Functions, extinctionM and second.extinctM

#Load libraries
library(bipartite)

###########
extinctionM <- function (web, participant = "both", method = "random", pol_ntw_hubs = NA, pol_mod_hubs = NA, pol_conn = NA, pol_periph = NA, plant_ntw_hubs = NA, plant_mod_hubs = NA, plant_conn = NA, plant_periph = NA) 
{
    partis <- c("lower", "higher", "both")
    partis.match <- pmatch(participant, partis)
    if (is.na(partis.match)) 
        stop("Choose participant: lower/higher/both.\n")
    meths <- c("random", "abundance", "degree", "modularity")
    meths.match <- pmatch(method, meths)
    if (is.na(meths.match)) 
        stop("Choose extinction method: random/abundance/degree/modularity.\n")
    nr <- NROW(web)
    nc <- NCOL(web)
    if (partis.match == 3) 
        partis.match <- sample(2, 1)
    if (meths.match == 1) {
        rexcl <- sample(nr, 1)
        cexcl <- sample(nc, 1)
        if (partis.match == 1) 
            web[rexcl, ] <- 0
        if (partis.match == 2) 
            web[, cexcl] <- 0
    }
    if (meths.match == 2) {
        rseq <- order(rowSums(web))
        cseq <- order(colSums(web))
        if (partis.match == 1) 
            web[rseq[1], ] <- 0
        if (partis.match == 2) 
            web[, cseq[1]] <- 0
    }
    if (meths.match == 3) {
        if (partis.match == 1) {
            sequ <- rowSums(web > 0)
            which.ex <- which(sequ == max(sequ))
            if (length(which.ex) > 1) 
                ex <- sample(which.ex, size = 1)
            else ex <- which.ex
            web[ex, ] <- 0
        }
        if (partis.match == 2) {
            sequ <- colSums(web > 0)
            which.ex <- which(sequ == max(sequ))
            if (length(which.ex) > 1) 
                ex <- sample(which.ex, size = 1)
            else ex <- which.ex
            web[, ex] <- 0
        }
    }
##new method (modularity), starts here 
	if (meths.match == 4) {        
        #exclude ntw_hubs first
        if(length(which(colnames(web) %in% pol_ntw_hubs) == TRUE) > 0) cexcl <- sample(colnames(web)[colnames(web) %in% pol_ntw_hubs], 1)
        if(length(which(rownames(web) %in% plant_ntw_hubs) == TRUE) > 0) rexcl <- sample(rownames(web)[rownames(web) %in% plant_ntw_hubs], 1)
        #exclude mod_hubs if there is no ntw_hubs remaining
        if(length(which(colnames(web) %in% pol_ntw_hubs) == TRUE) == 0){     
        	if(length(which(colnames(web) %in% pol_mod_hubs) == TRUE) > 0) cexcl <- sample(colnames(web)[colnames(web) %in% pol_mod_hubs], 1)
        	}
        if(length(which(rownames(web) %in% plant_ntw_hubs) == TRUE) == 0){     
       		if(length(which(rownames(web) %in% plant_mod_hubs) == TRUE) > 0) rexcl <- sample(rownames(web)[rownames(web) %in% plant_mod_hubs], 1)
       	}
       	#exclude conn ...
       	if(length(which(colnames(web) %in% c(pol_ntw_hubs, pol_mod_hubs)) == TRUE) == 0){         
       		if(length(which(colnames(web) %in% pol_conn) == TRUE) > 0) cexcl <- sample(colnames(web)[colnames(web) %in% pol_conn], 1)
       	}
       	if(length(which(rownames(web) %in% c(plant_ntw_hubs, plant_mod_hubs)) == TRUE) == 0){ 		        			if(length(which(rownames(web) %in% plant_conn) == TRUE) > 0) rexcl <- sample(rownames(web)[rownames(web) %in% plant_conn], 1)
		}
       	#exclude periph ...
		if(length(which(colnames(web) %in% c(pol_ntw_hubs, pol_mod_hubs,pol_conn)) == TRUE) == 0){  		
			if(length(which(colnames(web) %in% pol_periph) == TRUE) > 0) cexcl <- sample(colnames(web)[colnames(web) %in% pol_periph], 1)
		}
		if(length(which(rownames(web) %in% c(plant_ntw_hubs, plant_mod_hubs, plant_conn)) == TRUE) == 0){		        if(length(which(rownames(web) %in% plant_periph) == TRUE) > 0) rexcl <- sample(rownames(web)[rownames(web) %in% plant_periph], 1)
        }
		#remove excluded
        if (partis.match == 1) 
            web[rexcl, ] <- 0
        if (partis.match == 2) 
            web[, cexcl] <- 0
    }   
    
    return(web)
}





second.extinctM <- function (web, participant = "higher", method = "abun", nrep = 10, details = FALSE, pol_ntw_hubs = NA, pol_mod_hubs = NA, pol_conn = NA, pol_periph = NA, plant_ntw_hubs = NA, plant_mod_hubs = NA, plant_conn = NA, plant_periph = NA) 
{
    if (details == FALSE & pmatch(participant, c("both", "lower", 
        "higher")) == 1) {
        warning("\nFor random extinctions of both participants extinction sequences\n            will differ in length. Simply averaging sequences can hence not be used. Thus,\n            option 'details' will be set to FALSE internally.\n")
        details <- TRUE
    }
    one.second.extinct <- function(web = web, participant = participant, method = method) {
        dead <- matrix(nrow = 0, ncol = 3)
        colnames(dead) <- c("no", "ext.lower", "ext.higher")
        m2 <- web
        i = 1
        while (min(dim(m2)) > 1) {
        	#substitued extinction, for extinctionM
            n <- extinctionM(m2, participant = participant, method = method, pol_ntw_hubs = pol_ntw_hubs, pol_mod_hubs = pol_mod_hubs, pol_conn = pol_conn, pol_periph = pol_periph, plant_ntw_hubs = plant_ntw_hubs, plant_mod_hubs = plant_mod_hubs, plant_conn = plant_conn, plant_periph = plant_periph)
            dead <- rbind(dead, c(i, attributes(m2 <- empty(n, 
                count = TRUE))$empty))
            i <- i + 1
        }
        dead2 <- rbind(dead, c(NROW(dead) + 1, attributes(empty(m2, 
            count = TRUE))$empty))
        if (nrow(dead) != nrow(dead2)) {
            ci <- pmatch(participant, c("both", "lower", "higher"))
            dead2 <- matrix(0, nrow = ifelse(ci == 2, nrow(web), 
                ncol(web)), ncol = 3)
            colnames(dead2) <- colnames(dead)
            dead2[, 1] <- 1:nrow(dead2)
            dead2[1:nrow(dead), 2:3] <- dead[1:nrow(dead), 2:3]
        }
        else dead2 = dead
        dead2
    }
    if (is.vector(method)) 
        sequence = method
    if (pmatch(method, c("abundance", "random", "degree", "modularity")) == 
        1) {
        out <- one.second.extinct(web = web, participant = participant, method = method)
    }
    else {
        o <- replicate(nrep, one.second.extinct(web = web, participant = participant, method = method), simplify = FALSE)
        if (details) 
            out <- o
        else {
            z <- o[[1]]
            for (k in 2:length(o)) z <- z + o[[k]]
            out <- z/length(o)
        }
    }
    class(out) <- "bipartite"
    attr(out, "exterminated") <- c("both", "lower", "higher")[pmatch(participant, 
        c("both", "lower", "higher"))]
    out
}


####More examples:

#load data example
data(Safariland)
#need to define network roles (here I did that randomly)
pol_ntw_hubs <- colnames(web)[1:3]
pol_mod_hubs <- colnames(web)[4:10]
pol_conn <- colnames(web)[11:18]
pol_periph <- colnames(web)[19:27]
plant_ntw_hubs <- rownames(web)[1:3]
plant_mod_hubs <- rownames(web)[4:6]
plant_conn <- rownames(web)[7:8]
plant_periph <- rownames(web)[9]


#use of extinctionM with method random (same output as extinction)
extinctionM(web, participant = "both", method = "random", pol_ntw_hubs = NA, pol_mod_hubs = NA, pol_conn = NA, pol_periph = NA, plant_ntw_hubs = NA, plant_mod_hubs = NA, plant_conn = NA, plant_periph = NA)

#use of extinctionM with method modularity 
extinctionM(web, participant = "higher", method = "modularity", pol_ntw_hubs = pol_ntw_hubs, pol_mod_hubs = pol_mod_hubs, pol_conn = pol_conn, pol_periph = pol_periph, plant_ntw_hubs = plant_ntw_hubs, plant_mod_hubs = plant_mod_hubs, plant_conn = plant_conn, plant_periph = plant_periph)

#use of second.extinctM with method random (same output as second.extinct)
x <- second.extinctM(web, participant = "higher", method = "random", nrep = 10, details = FALSE, pol_ntw_hubs = NA, pol_mod_hubs = NA, pol_conn = NA, pol_periph = NA, plant_ntw_hubs = NA, plant_mod_hubs = NA, plant_conn = NA, plant_periph = NA)
robustness(x)

#use of second.extinctM with method modularity
x <- second.extinctM(web, participant = "lower", method = "modularity", nrep = 10, details = FALSE, pol_ntw_hubs = pol_ntw_hubs, pol_mod_hubs = pol_mod_hubs, pol_conn = pol_conn, pol_periph = pol_periph, plant_ntw_hubs = plant_ntw_hubs, plant_mod_hubs = plant_mod_hubs, plant_conn = plant_conn, plant_periph = plant_periph)
robustness(x)



