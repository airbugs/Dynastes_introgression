data.loci <- scan('GAntillesALLtaxaALLvar.loci', what = 'character', sep = '\n');#change name of the .loci file

break.lines <- grep('//', data.loci);

var.vec <- NULL;

for(i in 1:length(break.lines)){
	
	test <- data.loci[break.lines[i]]
	temp <- unlist(strsplit(test,''))
	n.var.sites <- sum(temp == '-') + sum(temp == '*')
	var.vec <- c(var.vec,n.var.sites)
	
}


#variable site position

var.site.vec <- NULL;

for(i in 1:length(break.lines)){
	
	test <- data.loci[break.lines[i]];
	temp <- unlist(strsplit(test,''));
	s1 <- which(temp == '-');
	sM <- which(temp == '*');
	s1ed <- s1 - 18;#change 17 to the number of space&character before the first DNA character
	sMed <- sM - 18;#change 17 to the number of space&character before the first DNA character
	var.site.vec <- c(var.site.vec, s1ed, sMed);
	
}

par(mfrow = c(2,1));
hist(var.vec, breaks = c(seq(-1, 100, by = 1)), xlim = c(-1,50), main = 'Number of segregating sites per locus', xlab = 'Number of segregating sites');
abline(v = quantile(var.vec, 0.025), lwd = 1.5, lty = 2, col = 'red');
abline(v = quantile(var.vec, 0.5), lwd = 2, col = 'red');
abline(v = quantile(var.vec, 0.975), lwd = 1.5, lty = 2, col = 'red');
hist(var.site.vec, xlim = c(-1,150), breaks = c(-1:150), xlab = 'Position along the loci', main = 'The position of segregating sites');
abline(h = 320, col = 'red', lwd = 2);
#abline(v = 110, col = 'red', lwd = 2);
#abline(v = 100, col = 'red', lty = 2, lwd = 1.5);
abline(v = 130, col = 'red', lty = 2, lwd = 1.5);

#the following is for plotting the after edited data. including the after ed statistics
stat.file <- read.csv('afterchopped_EDloci_stat.csv');

par(mfrow = c(2,2));
hist(var.vec, breaks = c(seq(-1, 100, by = 1)), xlim = c(-1,20), main = 'Number of segregating sites per locus', xlab = 'Number of segregating sites');
abline(v = quantile(var.vec, 0.025), lwd = 1.5, lty = 2, col = 'red');
abline(v = quantile(var.vec, 0.5), lwd = 2, col = 'red');
abline(v = quantile(var.vec, 0.975), lwd = 1.5, lty = 2, col = 'red');
hist(var.site.vec, xlim = c(-1,150), breaks = c(0:150), xlab = 'Position along the loci', main = 'The position of segregating sites');

#then the summaries dis of theta & pd
#par(mfrow = c(2,1));
hist(stat.file[,1], main = 'Number of taxa per locus', xlab = 'Number of taxa', xlim = c(3,9), breaks = c(-1:20));
hist(stat.file[,2]*100, main = 'Maximum p-distance', xlab = '% pairwise divergence (interspecific)', xlim = c(-0.5,5.5), breaks = c(seq(-0.5, 15.5, by = 1)));