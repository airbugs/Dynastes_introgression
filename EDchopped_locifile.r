#editing .loci file based on per site theta within taxon
#& also % pairwise divergence between taxa
#Additionally, remove loci that have less than 4 taxa & those that are not variable
#need a library names pegas (the function theta.s)
#library(pegas);
chopped.loci <- scan('filename.loci', what = 'character', sep = '\n');#change file name
break.lines <- grep('//', chopped.loci);
index.lines <- c(0, break.lines);

loci.seq.length <- 110;#how long are the aligned sequences (in the Dynastes example it is set to 110 basepairs
#max.theta.intra <- 0.02;#dont need this here
max.div.inter <- 0.15;#change number here, maximum p distance between species
null.data <- rep(0, length(break.lines)*2);#for creating a data matrix
stat.mat <- matrix(null.data, nrow = length(break.lines), ncol = 2);#data matrix for sum stats
colnames(stat.mat) <- c('N_taxa', 'max_p_dis');

for(i in 1:length(break.lines)){
	start.line <- index.lines[i] + 1;
	end.line <- break.lines[i] - 1;
	temp.alignment <- chopped.loci[start.line:end.line];
	
	s.pattern <- '>(.+)_.+[[:space:]]+(.+)';#
	sp.vec <- gsub(s.pattern, '\\1', temp.alignment);
	seq.vec <- gsub(s.pattern, '\\2', temp.alignment);
	uni.sp.vec <- unique(sp.vec);
	stat.mat[i,1] <- length(uni.sp.vec);	

# if you want to filter loci by the theta value (e.g. population genetic studies) you may want to use the following:
#	theta.vec <- NULL;#calculate theta first
#	for(iter in 1:length(uni.sp.vec)){
#		temp.sp <- uni.sp.vec[iter];
#		temp.sample.index <- grep(temp.sp, sp.vec);
#		if(length(temp.sample.index) > 1){#at least 2 individual in a taxon
#			var.sites <- NULL;
#			for(ii in 1:(length(temp.sample.index)-1)){
#				tt.seq1 <- unlist(strsplit(seq.vec[temp.sample.index[ii]], ''));
#				tt.seq2 <- unlist(strsplit(seq.vec[temp.sample.index[ii + 1]], ''));
#				temp.var.sites <- which(tt.seq1 != tt.seq2 & tt.seq1 != 'N' & tt.seq1 != '-' & tt.seq2 != 'N' & tt.seq2 != '-');
#				var.sites <- c(var.sites, temp.var.sites);
#			}
#			n.segregating.sites <- length(unique(var.sites));
#			temp.theta <- theta.s(n.segregating.sites, length(temp.sample.index))/loci.seq.length;
#			theta.vec <- c(theta.vec, temp.theta);
#		}else{#only one ind for a given taxon, just fill in 0
#			theta.vec <- c(theta.vec, 0);
#		}
#	}
#	max.theta <- max(theta.vec);
#	stat.mat[i, 2] <- max.theta;
		
	p.dis.vec <- NULL;#p.dis across samples from the same or diff taxa
	for(iii in 1:(length(seq.vec)-1)){
		ttt.seq1 <- unlist(strsplit(seq.vec[iii], ''));
		for(itt in (iii+1):length(seq.vec)){
			ttt.seq2 <- unlist(strsplit(seq.vec[itt], ''));
			dis.sites <- which(ttt.seq1 != ttt.seq2 & ttt.seq1 != 'N' & ttt.seq2 != 'N' & ttt.seq1 != '-' & ttt.seq2 != '-');
			p.dis <- length(dis.sites)/loci.seq.length;
			p.dis.vec <- c(p.dis.vec, p.dis);
			#cat(iii, '\t', itt, '\n');
		}
	}
	max.p.dis <- max(p.dis.vec);
	stat.mat[i, 2] <- max.p.dis;					
	cat(i, '\n');	
}

write.table(stat.mat, file = 'stat_chopped_loci.csv', sep = ',', quote = FALSE, row.names = FALSE);#save it for future use

less.4.taxa <- which(stat.mat[,1] < 6);#those loci that have less than 4 taxa;
#bad.theta <- which(stat.mat[,2] > max.theta.intra & stat.mat[,2] != 'NA');
bad.pd <- which(stat.mat[,2] > max.div.inter & stat.mat[,2] != 'NA');
non.var <- which(stat.mat[,2] == 0);#those that have no variable sites between taxa (so also within taxon);
unwanted.loci <- unique(c(less.4.taxa, bad.pd, non.var));

new.stat.mat <- stat.mat[-unwanted.loci,];
write.table(new.stat.mat, file = 'afterchopped_EDloci_stat.csv', sep = ',', quote = FALSE, row.names = FALSE);#save it for future use

unwanted.lines <- NULL
for(ite in 1:length(unwanted.loci)){
	loci.id <- unwanted.loci[ite];
	start.line <- index.lines[loci.id] + 1;
	end.line <- break.lines[loci.id];
	range.lines <- start.line:end.line;
	unwanted.lines <- c(unwanted.lines, range.lines);
	cat(ite, '\n')
}

new.loci <- chopped.loci[-unwanted.lines];
write.table(new.loci, file = 'DyEDpd15AllSP.loci', quote = FALSE, sep = '\n', row.names = FALSE, col.names = FALSE);#change file name to your preference
