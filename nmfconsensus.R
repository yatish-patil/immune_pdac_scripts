nmfconsensus <- function(input.ds, k.init, k.final, num.clusterings, maxniter, error.function, rseed=123456789, stopconvergence = 40, stopfrequency = 10, non.interactive.run = F, doc.string = "", ...) {

#
#  GenePattern Methodology for:
#
#  Metagenes and Molecular Pattern Discovery using Matrix Factorization
#  Jean-Philippe Brunet, Pablo Tamayo, Todd. R. Golub, and Jill P. Mesirov
# 
#  Author:  Pablo Tamayo (tamayo@genome.wi.mit.edu)
#
#  Based on the original matlab version written by Jean-Philippe Brunet (brunet@broad.mit.edu) and
#  with additional contributions from: Ted Liefeld (liefeld@broad.mit.edu)   
#  Date:  November 27, 2003
#
#  Last change March 3, 2005: modifications to make the output more readable.
#
#  Execute from an R console window with this command:
#  source("<this file>", echo = TRUE)
#  E.g. someoutput <- mynmf2(input.ds="c:\\nmf\\all_aml.res",k.init=2,k.final=5,num.clusterings=20,maxniter=500, error.function="euclidean") 
#
#  For details on the method see:
#
#  Proc. Natl. Acad. Sci. USA 2004 101: 4164-4169
#  http://www.broad.mit.edu/cgi-bin/cancer/publications/pub_paper.cgi?mode=view&paper_id=89
#
#  Input parameters
#
#   input.ds
#                       input gene expression dataset in GCT or RES format
#   k.init
#                       initial value of k
#   k.final
#                       final value of k
#   num.clusterings
#                       number of NMF clusterings to build consensus matrix
#   maxniter
#                       maximum number of NMF iterations
#   error.function
#                       NMF error function: "divergence" of "euclidean"
#   rseed
#                       random number generator seed
#   stopconv
#                       how many no change checks are needed to stop NMF iterations (convergence)
#   stopfreq
#                       frequency (NMF iterations) of "no change" checks 
#   non.interactive.run 
#                       flag controlling if the plots are produced interatively (Rgui and saved) or only saved in files
#   doc.string
#                       prefix to be added to the output files
#
#  Output files are (prefix with "doc.string")
#
#   params.txt 
#                       run parameters and time of execution
#   membership.gct		
#			membership results for samples at all values of K
#   cophenetic.txt 
#			cophenetic values for each K
#   cophenetic.plot.jpeg
#			plot of cophenetic for each value of K		
#   consensus.k#.gct (for each value of K)
#			consensus matrix for k=#
#   consensus.plot.k#.jpeg (for each value of K)
#			plot of consensus matrix for k=#
#   graphs.k#.jpeg (for each value of K)

# save input parameters
directory = ""
time.string <- as.character(as.POSIXlt(Sys.time(),"GMT"))
if (non.interactive.run == "F"){ 
	filename <- paste(directory, doc.string, ".params.txt", sep="", collapse="")  
	write(paste("Run of NMF on ", time.string), file=filename)

	write(paste("input.ds =", input.ds, sep=" "), file=filename, append=T) 
	write(paste("k.init = ", k.init, sep=" "), file=filename, append=T) 
	write(paste("k.final =", k.final, sep=" "), file=filename, append=T) 
	write(paste("num.clusterings =", num.clusterings, sep=" "), file=filename, append=T) 
	write(paste("maxniter =", maxniter, sep=" "), file=filename, append=T) 
	write(paste("error.function =", error.function, sep=" "), file=filename, append=T) 
	write(paste("rseed =", rseed, sep=" "), file=filename, append=T) 
	write(paste("directory =", directory, sep=" "), file=filename, append=T) 
	write(paste("stopconv =", stopconvergence, sep=" "), file=filename, append=T) 
	write(paste("stopfreq =", stopfrequency, sep=" "), file=filename, append=T)
	write(paste("non.interctive.run =", non.interactive.run, sep=" "), file=filename, append=T) 
	write(paste("doc.string =", doc.string, sep=" "), file=filename, append=T) 
}


k.init<-as.integer(k.init)
k.final<-as.integer(k.final)
num.clusterings<-as.integer(num.clusterings)
n.iter<-as.integer(maxniter)
if (!is.na(rseed)){
     seed <- as.integer(rseed)
}
stopfreq <- as.integer(stopfrequency)
stopconv <- as.integer(stopconvergence)
# library(mva)
# library(MASS)
# library(GenePattern)

D <- read.dataset(input.ds)
A <- data.matrix(D)

# Threshold negative values to small quantity 

eps <- .Machine$double.eps
A[A < 0] <- eps



cols <- length(A[1,])
rows <- length(A[,1])

col.names <- names(D)

num.k <- k.final - k.init + 1

rho <- vector(mode = "numeric", length = num.k)
k.vector <- vector(mode = "numeric", length = num.k)

k.index <- 1

connect.matrix.ordered <- array(0, c(num.k, cols, cols))

for (k in k.init:k.final) {

   if (non.interactive.run == F) {
         if (.Platform$OS.type == "windows") {
             filename <- paste(directory, doc.string, "graphs.k", k, sep="", collapse="")
             windows(width = 9, height = 11)
         } else if (.Platform$OS.type == "unix") {
             filename <- paste(directory, doc.string,  "graphs.k", k, ".pdf", sep="", collapse="")
             pdf(file=filename, width = 9, height = 11)
         }
   } else {
         if (.Platform$OS.type == "unix") {
             filename <- paste(directory, doc.string,  "graphs.k", k, ".pdf", sep="", collapse="")
             pdf(file=filename, width = 9, height = 11)
         } else if (.Platform$OS.type == "windows") {
             filename <- paste(directory, doc.string, "graphs.k", k, ".pdf", sep="", collapse="")
             pdf(file=filename, width = 9, height = 11)
         }
   }

   nf <- layout(matrix(c(1,2,3,4,5,6,7,8), 4, 2, byrow=T), c(1, 1, 1, 1), c(1, 1), TRUE)
   assign <- matrix(0, nrow = num.clusterings, ncol = cols)

   for (i in 1:num.clusterings) {
	  
	if (non.interactive.run == "F"){
        print(paste("Computing clustering number=", i, " for k=", k, sep="", collapse=" "))
	}

        if (error.function == "divergence"){
	    NMF.out <- NMF.div(V = A, k = k, maxniter = n.iter, seed = seed + i, stopconv = stopconv, stopfreq = stopfreq)
	} else if (error.function == "euclidean"){
	    NMF.out <- NMF(V = A, k = k, maxniter = n.iter, seed = seed + i, stopconv = stopconv, stopfreq = stopfreq)
	} else {
            stop(paste("Un-supported error function=", error.function, sep=""))
        }
	if (non.interactive.run == "F"){
        print(paste(NMF.out$t, " NMF iterations performed", sep="", collapse=" "))
	}

        for (j in 1:cols) { # Find membership
            class <- order(NMF.out$H[,j], decreasing=T)
            assign[i, j] <- class[1]
        }

	if (i == 3) {  # '2''20' example for first clustering iteration
            H.saved <- NMF.out$H
            sub.string <- paste(doc.string, " k=", k, sep="")
            plot(1:NMF.out$t, NMF.out$error.v[1:NMF.out$t], pch = 20, cex = 1.5, col = 1, xlab="time", ylab="NMF error", sub=sub.string, main=paste("Example of NMF convergence plot k=", k, sep=""))


            if (rows < 1000) {
               W <- NMF.out$W
            } else {
               W <- NMF.out$W[sample(x = 1:rows, size = 1000),]
            }
            sub.string <- paste(doc.string, " k=", k, sep="")
            matrix.abs.plot(W, sub = sub.string, log = F, main = "Example W matrix (orig. ordering)", ylab = "genes", xlab ="metasamples")
            matrix.abs.plot(H.saved, sub = sub.string, log = F, main = "Example H matrix (orig. ordering)", ylab = "metagenes", xlab ="samples")
            metagene.plot(H = H.saved, main = "Example metagenes (orig. ordering)", sub = sub.string, xlab = "samples", ylab = "metagenes")

        }

        # Anguraj added to find the metagenes
         filename <- paste(directory, doc.string,  "metagenes.k.",k,".",i, ".txt", sep="", collapse="")
     write.table(NMF.out$W, filename,sep="\t")
     
     filename1 <- paste(directory, doc.string,  "metasamples.k.",k,".",i, ".txt", sep="", collapse="")
          write.table(NMF.out$H, filename1,sep="\t")

        rm(NMF.out)

     }  ## end  for (i in 1:num.clusterings)

   
     # compute consensus matrix
     connect.matrix <- matrix(0, nrow = cols, ncol = cols)

     for (i in 1:num.clusterings) {
       for (j in 1:cols) {
          for (p in 1:cols) {
             if (j != p) {
                  if (assign[i, j] == assign[i, p]) {
                    connect.matrix[j, p] <- connect.matrix[j, p] + 1
                  } 
              } else {
                    connect.matrix[j, p] <- connect.matrix[j, p] + 1
              }
           }
       }
     }

     connect.matrix <- connect.matrix / num.clusterings

     dist.matrix <- 1 - connect.matrix
     dist.matrix <- as.dist(dist.matrix)
     HC <- hclust(dist.matrix, method="average")

     dist.coph <- cophenetic(HC)
     k.vector[k.index] <- k
     rho[k.index] <- cor(dist.matrix, dist.coph)
     rho[k.index] <- signif(rho[k.index], digits = 4)
   
#     connect.matrix.ordered <- matrix(0, nrow=cols, ncol = cols)

     for (i in 1:cols) {
        for (j in 1:cols) {
           connect.matrix.ordered[k.index, i, j] <- connect.matrix[HC$order[i], HC$order[j]]
         }
     }

     # compute consensus clustering membership

     membership <- cutree(HC, k = k)

     max.k <- max(membership)
     items.names.ordered <- col.names[HC$order]
     membership.ordered <- membership[HC$order]
     results <- data.frame(cbind(membership.ordered, items.names.ordered))
     
     if (k > k.init){
          all.membership <- cbind(all.membership, membership);
     } else {
          all.membership <- cbind(membership);
     }

     sub.string <- paste(doc.string, " k=", k, sep="")
     matrix.abs.plot(connect.matrix.ordered[k.index,,], sub=sub.string, log = F, main = "Ordered Consensus Matrix", ylab = "samples", xlab ="samples")
     plot(HC, xlab="samples", cex = 0.75, labels = col.names, sub = sub.string, col = "blue", main = paste("Ordered Linkage Tree. Coph=", rho[k.index]))

     resultsGct <- data.frame(membership.ordered)
     row.names(resultsGct) <- items.names.ordered
     
     filename <- paste(directory, doc.string,  "consensus.k.",k, ".gct", sep="", collapse="")
     write.gct(resultsGct, filename)

     H.sorted <- H.saved[,HC$order]
     sub.string <- paste(doc.string, " k=", k, sep="")
   
     matrix.abs.plot(H.sorted, sub = sub.string, log = F, main = "Example H matrix (ordered)", ylab = "metagenes", xlab ="samples")
     metagene.plot(H = H.sorted, sub = sub.string, main = "Example metagenes (ordered)", xlab = "samples", ylab = "metagenes")

     if (non.interactive.run == F) {  
           if (.Platform$OS.type == "windows") {
               savePlot(filename = filename, type ="jpeg", device = dev.cur())
           } else if (.Platform$OS.type == "unix") {
               dev.off()
           }
      } else {
           dev.off()
      }

   if (non.interactive.run == F) {
         if (.Platform$OS.type == "windows") {
             filename <- paste(directory, doc.string,  "consensus.plot.k", k, sep="", collapse="")
             windows(width = 8.5, height = 11)
         } else if (.Platform$OS.type == "unix") {
             filename <- paste(directory, doc.string,  "consensus.plot.k", k, ".pdf", sep="", collapse="")
             pdf(file=filename, width = 8.5, height = 11)
         }
   } else {
         if (.Platform$OS.type == "unix") {
            filename <- paste(directory, doc.string,  "consensus.plot.k", k, ".pdf", sep="", collapse="")
             pdf(file=filename, width = 8.5, height = 11)
         } else if (.Platform$OS.type == "windows") {
            filename <- paste(directory, doc.string,  "consensus.plot.k", k, ".pdf", sep="", collapse="")
            pdf(file=filename, width = 8.5, height = 11)
         }
   }

     nf <- layout(matrix(c(1), 1, 1, byrow=T), c(1, 1), c(1, 1), TRUE)

     conlabel <- paste("Consensus k =", k, sep=" ", collapse="")

     sub.string <- paste("Consensus matrix k=", k, "; dataset= ", input.ds, sep="")
     ConsPlot(connect.matrix.ordered[k.index,,], col.labels = membership.ordered, col.names = items.names.ordered, main = " ", sub=sub.string, xlab=" ", ylab=" ")

      if (non.interactive.run == F) {  
           if (.Platform$OS.type == "windows") {
               savePlot(filename = filename, type ="jpeg", device = dev.cur())
           } else if (.Platform$OS.type == "unix") {
               dev.off()
           }
      } else {
           dev.off()
      }
  
     k.index <- k.index + 1

} # end of loop over k


# Save consensus matrices in one file

   if (non.interactive.run == F) {
         if (.Platform$OS.type == "windows") {
             filename <- paste(directory, doc.string, "consensus.all.k.plot", sep="")

             windows(width = 8.5, height = 11)
         } else if (.Platform$OS.type == "unix") {

             filename <- paste(directory, doc.string,  "consensus.all.k.plot.pdf", sep="")
             pdf(file=filename, width = 8.5, height = 11)
         }
   } else {
         if (.Platform$OS.type == "unix") {
            filename <- paste(directory, doc.string,  "consensus.all.k.plot.pdf", sep="")
            pdf(file=filename, width = 8.5, height = 11)
         } else if (.Platform$OS.type == "windows") {
            filename <- paste(directory, doc.string,  "consensus.all.k.plot.pdf", sep="")
            pdf(file=filename, width = 8.5, height = 11)
         }
   }

   nf <- layout(matrix(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16), 4, 4, byrow=T), c(1, 1, 1, 1), c(1, 1, 1, 1), TRUE)

  for (k in 1:num.k) { 
     matrix.abs.plot(connect.matrix.ordered[k,,], log = F, main = paste("k=", k.vector[k]), 
                          sub = paste("Cophenetic coef.=", rho[k]), ylab = "samples", xlab ="samples")
  }
   
  y.range <- c(1 - 2*(1 - min(rho)), 1)
  plot(k.vector, rho, main ="Cophenetic Coefficient", xlim=c(k.init, k.final), ylim=y.range, 
                          xlab = "k", ylab="Cophenetic correlation", type = "n")


  lines(k.vector, rho, type = "l", col = "black")
  points(k.vector, rho, pch=22, type = "p", cex = 1.25, bg = "black", col = "black")

  if (non.interactive.run == F) {  
        if (.Platform$OS.type == "windows") {
            savePlot(filename = filename, type ="jpeg", device = dev.cur())
        } else if (.Platform$OS.type == "unix") {
            dev.off()
        }
   } else {
        dev.off()
   }

   if (non.interactive.run == F) {
         if (.Platform$OS.type == "windows") {
             filename <- paste(directory, doc.string,  "cophenetic.plot", sep="")
             windows(width = 8.5, height = 11)
         } else if (.Platform$OS.type == "unix") {
             filename <- paste(directory, doc.string,  "cophenetic.plot.pdf", sep="")
             pdf(file=filename, width = 8.5, height = 11)
         }
   } else {
         if (.Platform$OS.type == "unix") {
            filename <- paste(directory, doc.string,  "cophenetic.plot.pdf", sep="")
             pdf(file=filename, width = 8.5, height = 11)
         } else if (.Platform$OS.type == "windows") {
            filename <- paste(directory, doc.string,  "cophenetic.plot.pdf", sep="")
             pdf(file=filename, width = 8.5, height = 11)
         }
   }


# Write the membership matrix

resultsmembership <- data.frame(all.membership)
row.names(resultsmembership) <- col.names
filename <- paste(directory, doc.string,  "membership", ".gct", sep="", collapse="")
write.gct(resultsmembership , filename)

y.range <- c(1 - 2*(1 - min(rho)), 1)
plot(k.vector, rho, main ="Cophenetic Coefficient", xlim=c(k.init, k.final), ylim=y.range, xlab = "k", ylab="Cophenetic correlation", type = "n")
lines(k.vector, rho, type = "l", col = "black")
points(k.vector, rho, pch=22, type = "p", cex = 1.25, bg = "black", col = "black")

      if (non.interactive.run == F) {  
           if (.Platform$OS.type == "windows") {
               savePlot(filename = filename, type ="jpeg", device = dev.cur())
           } else if (.Platform$OS.type == "unix") {
               dev.off()
           }
      } else {
           dev.off()
      }

xx <- cbind(k.vector, rho)
write(xx, file= paste(directory, doc.string, ".", "cophenetic.txt", sep=""))
}


read.dataset <- function(file) {
	result <- regexpr(paste(".txt","$",sep=""), tolower(file))

	if(result[[1]] != -1)
		return(read.gct(file))
	result <- regexpr(paste(".res","$",sep=""), tolower(file))
	if(result[[1]] != -1)
		return(read.res(file))
	stop("Input is not a res or gct file.")	
}

matrix.abs.plot <- function(V, axes = F, log = F, norm = T, transpose = T, matrix.order = T, max.v = 1, min.v = 0, main = " ", sub = " ", xlab = " ", ylab = "  ") {
      rows <- length(V[,1])
      cols <- length(V[1,])
      if (log == T) {
         V <- log(V)
      }
      B <- matrix(0, nrow=rows, ncol=cols)
	for (i in 1:rows) {
           for (j in 1:cols) {
                if (matrix.order == T) {
                   k <- rows - i + 1
                } else {
                   k <- i
                }
                if (norm == T) {
                  if ((max.v == 1) && (min.v == 0)) {
                     max.val <- max(V)
                     min.val <- min(V)
                  } else {
		     	   max.val = max.v
                     min.val = min.v
                  }
               }
	     B[k, j] <-  max.val - V[i, j] + min.val
           }
      }
	if (transpose == T) {
	  B <- t(B)
        }
	if (norm == T) {
#removed gamma = 1.5
          image(z = B, zlim = c(min.val, max.val), axes = axes, col = rainbow(100, s = 1.0, v = 0.75, start = 0.0, end = 0.75), main = main, sub = sub, xlab = xlab, ylab = ylab) 
      } else {
            image(z = B, axes = axes, col = rainbow(100, s = 1, v = 0.6, start = 0.1, end = 0.9, gamma = 1), main = main, sub = sub, xlab = xlab, ylab = ylab) 
      }
      return(list(B, max.val, min.val))
}

metagene.plot <- function(H, main = " ", sub = " ", xlab = "samples ", ylab = "amplitude") {
	k <- length(H[,1])
	S <- length(H[1,])
	index <- 1:S
	maxval <- max(H)
        minval <- min(H)
	plot(index, H[1,], xlim=c(1, S), ylim=c(minval, maxval), main = main, sub = sub, ylab = ylab, xlab = xlab, type="n")
	for (i in 1:k) {
	    lines(index, H[i,], type="l", col = i, lwd=2)
        }
}

NMF <- function(V, k, maxniter = 2000, seed = 123456, stopconv = 40, stopfreq = 10) {
        N <- length(V[,1])
        M <- length(V[1,])
        set.seed(seed)
        W <- matrix(runif(N*k), nrow = N, ncol = k)  # Initialize W and H with random numbers
        H <- matrix(runif(k*M), nrow = k, ncol = M)


        VP <- matrix(nrow = N, ncol = M)


        error.v <- vector(mode = "numeric", length = maxniter)
        new.membership <- vector(mode = "numeric", length = M)
        old.membership <- vector(mode = "numeric", length = M)
        eps <- .Machine$double.eps
        for (t in 1:maxniter) {

              VP = W %*% H

              H <- H * (crossprod(W, V)/crossprod(W, VP)) + eps
		  VP = W %*% H
              H.t <- t(H)
              W <- W * (V %*% H.t)/(VP %*% H.t) + eps
              error.v[t] <- sqrt(sum((V - VP)^2))/(N * M)
               if (t %% stopfreq == 0) {
                    for (j in 1:M) {
                        class <- order(H[,j], decreasing=T)
                        new.membership[j] <- class[1]
                     }
                     if (sum(new.membership == old.membership) == M) {
                        no.change.count <- no.change.count + 1
                     } else {
                        no.change.count <- 0
                     }
                     if (no.change.count == stopconv) break
                     old.membership <- new.membership
               }
        }
       
        return(list(W = W, H = H, t = t, error.v = error.v))
}

NMF.div <- function(V, k, maxniter = 2000, seed = 123456, stopconv = 40, stopfreq = 10) {

        N <- length(V[,1])
        M <- length(V[1,])
        set.seed(seed)
        W <- matrix(runif(N*k), nrow = N, ncol = k)  # Initialize W and H with random numbers
        H <- matrix(runif(k*M), nrow = k, ncol = M)
        VP <- matrix(nrow = N, ncol = M)
        error.v <- vector(mode = "numeric", length = maxniter)
        new.membership <- vector(mode = "numeric", length = M)
        old.membership <- vector(mode = "numeric", length = M)
        no.change.count <- 0
        eps <- .Machine$double.eps
        for (t in 1:maxniter) {
                VP = W %*% H
                W.t <- t(W)
                H <- H * (W.t %*% (V/VP)) + eps
                norm <- apply(W, MARGIN=2, FUN=sum)
                for (i in 1:k) {
                    H[i,] <- H[i,]/norm[i]
                }
                VP = W %*% H
                H.t <- t(H)
                W <- W * ((V/VP) %*% H.t) + eps
                norm <- apply(H, MARGIN=1, FUN=sum)
                for (i in 1:k) {
                    W[,i] <- W[,i]/norm[i]
                }
               error.v[t] <- sum(V * log((V + eps)/(VP + eps)) - V + VP)/(M * N)
               if (t %% stopfreq == 0) {

                    for (j in 1:M) {
                        class <- order(H[,j], decreasing=T)
                        new.membership[j] <- class[1]
                     }
                     if (sum(new.membership == old.membership) == M) {
                        no.change.count <- no.change.count + 1
                     } else {
                        no.change.count <- 0
                     }
                     if (no.change.count == stopconv) break
                     old.membership <- new.membership
               }
        }
        return(list(W = W, H = H, t = t, error.v = error.v))
}

ConsPlot <- function(V, col.labels, col.names, main = " ", sub = " ", xlab=" ", ylab=" ") {

# Plots a heatmap plot of a consensus matrix

     cols <- length(V[1,])
     B <- matrix(0, nrow=cols, ncol=cols)
     max.val <- max(V)
     min.val <- min(V)
     for (i in 1:cols) {
         for (j in 1:cols) {
             k <- cols - i + 1
	     B[k, j] <-  max.val - V[i, j] + min.val
          }
     }

     col.names2 <- rev(col.names)
     col.labels2 <- rev(col.labels)
     D <- matrix(0, nrow=(cols + 1), ncol=(cols + 1))

     col.tag <- vector(length=cols, mode="numeric")
     current.tag <- 0
     col.tag[1] <- current.tag
     for (i in 2:cols) {
        if (col.labels[i] != col.labels[i - 1]) {
             current.tag <- 1 - current.tag
        }
        col.tag[i] <- current.tag
     }
     col.tag2 <- rev(col.tag)
     D[(cols + 1), 2:(cols + 1)] <- ifelse(col.tag %% 2 == 0, 1.02, 1.01)
     D[1:cols, 1] <- ifelse(col.tag2 %% 2 == 0, 1.02, 1.01)
     D[(cols + 1), 1] <- 1.03
     D[1:cols, 2:(cols + 1)] <- B[1:cols, 1:cols]

     #gamma = 1.5 removed
     col.map <- c(rainbow(100, s = 1.0, v = 0.75, start = 0.0, end = 0.75), "#BBBBBB", "#333333", "#FFFFFF")
     image(1:(cols + 1), 1:(cols + 1), t(D), col = col.map, axes=FALSE, main=main, sub=sub, xlab= xlab, ylab=ylab)
     for (i in 1:cols) {
         col.names[i]  <- paste("      ", substr(col.names[i], 1, 12), sep="")
         col.names2[i] <- paste(substr(col.names2[i], 1, 12), "     ", sep="")
     }

     axis(2, at=1:cols, labels=col.names2, adj= 0.5, tick=FALSE, las = 1, cex.axis=0.50, font.axis=1, line=-1)
     axis(2, at=1:cols, labels=col.labels2, adj= 0.5, tick=FALSE, las = 1, cex.axis=0.65, font.axis=1, line=-1)

     axis(3, at=2:(cols + 1), labels=col.names, adj= 1, tick=FALSE, las = 3, cex.axis=0.50, font.axis=1, line=-1)
     axis(3, at=2:(cols + 1), labels=as.character(col.labels), adj = 1, tick=FALSE, las = 1, cex.axis=0.65, font.axis=1, line=-1)

     return()
}

read.res <- function(filename = "NULL") { 
#
# Reads a gene expression dataset in RES format and converts it into an R data frame
#
   header.cont <- readLines(filename, n = 1)
   temp <- unlist(strsplit(header.cont, "\t"))
   colst <- length(temp)
   header.labels <- temp[seq(3, colst, 2)]
   ds <- read.delim(filename, header=F, row.names = 2, sep="\t", quote="", skip=3, blank.lines.skip=T, comment.char="", as.is=T)
   colst <- length(ds[1,])
   cols <- (colst - 1)/2
   rows <- length(ds[,1])
   A <- matrix(nrow=rows - 1, ncol=cols)
   A <- ds[1:rows, seq(2, colst, 2)]
   table1 <- data.frame(A)
   names(table1) <- header.labels
   return(table1)
}

read.gct <- function(filename = "NULL") { 
#
# Reads a gene expression dataset in GCT format and converts it into an R data frame
#
ds <- read.delim(filename, header=T,  quote="",   row.names=1, blank.lines.skip=T, comment.char="", as.is=T,na.strings="NA")
#   ds <- read.delim(filename, header=T,  quote="",  skip=2, row.names=1, blank.lines.skip=T, comment.char="", as.is=T)
  # ds <- ds[-1]
   return(ds)
}

write.gct <- function (gct, filename) 
{
    f <- file(filename, "w")

    gct_1<-cbind(rownames(gct),gct)
    colnames(gct_1)[1]<-"Name"

    write.table(gct_1,file = f, row.names=FALSE, sep="\t")
#    cat("#1.2", "\n", file = f, append = TRUE, sep = "")
 #   cat(dim(gct)[1], "\t", dim(gct)[2], "\n", file = f, append = TRUE, sep = "")

#    cat("Name", file = f, append = TRUE, sep = "")
   # cat("Description", file = f, append = TRUE, sep = "")

 #   names <- names(gct)
  #  for (j in 1:length(names)) {
  #      cat("\t", names[j], file = f, append = TRUE, sep = "")
  #  }
   # cat("\n", file = f, append = TRUE, sep = "")
   # oldWarn <- options(warn = -1)

   # m <- matrix(nrow = dim(gct)[1], ncol = dim(gct)[2] +  1)
   # m[, 1] <- row.names(gct)
   # m[, 2] <- row.names(gct)
   # index <- 2
   # for (i in 1:dim(gct)[2]) {
    #    m[, index] <- gct[, i]
    #    index <- index + 1
   # }
   # write.table(m, file = f, append = TRUE, quote = FALSE, sep = "\t", eol = "\n", col.names = FALSE, row.names = FALSE)
   # close(f)
   # options(warn = 0)
   return(gct)
}

