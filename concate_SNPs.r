#from unlinked SNPs to bpp inputs
data.loci <- scan('DyEDpd15AllSP.loci', what = 'character', sep = '\n');
#this .loci file should be already edited with a max theta within taxon & a max pd between taxa
n.ind <- 61;#change this number to the sample size of yours

break.lines <- grep('//', data.loci);
index.lines <- c(0, break.lines);

#create a matrix to store the SNP data
missing <- rep('?', n.ind*length(break.lines));#the first number is the n of samples the 2nd is # of loci
data.matrix <- matrix(missing, nrow = n.ind, ncol = length(break.lines));
stats <- read.csv('samplelist.csv');#this is for my data only
sample.list <- stats$taxa;#a vector of all your sample names. you need to somehow create this vector by yourself
rownames(data.matrix) <- sample.list;

#a for loop to fill in the matrix with the first variable site per locus
for(iter in 1:length(break.lines)){	
	n.samples <- (index.lines[iter + 1] - index.lines[iter]) - 1;
	id.line <- break.lines[iter];
	temp.snp.n <- unlist(strsplit(data.loci[id.line], ''));
	m.snp.id <- which(temp.snp.n == '*');
	s.snp.id <- which(temp.snp.n == '-');
	if(length(m.snp.id) > 0){
		for(i in 1:n.samples){
			temp.index <- m.snp.id[1] - 17;
			sample.index <- index.lines[iter] + i;
			temp.sample <- data.loci[sample.index];
			s.pattern <- '>.+_([[:alnum:]]+)[[:space:]]+(.+)';
			sample.ind <- sub(s.pattern, '\\1', temp.sample);
			seq.con <- sub(s.pattern, '\\2', temp.sample);
			seq.separate <- unlist(strsplit(seq.con, ''));
			temp.SNP <- seq.separate[temp.index];
			data.matrix[sample.ind, iter] <- temp.SNP;
		}
	}else if(length(s.snp.id) > 0){
		for(i in 1:n.samples){
			temp.index <- s.snp.id[1] - 17;
			sample.index <- index.lines[iter] + i;
			temp.sample <- data.loci[sample.index];
			s.pattern <- '>.+_([[:alnum:]]+)[[:space:]]+(.+)';
			sample.ind <- sub(s.pattern, '\\1', temp.sample);
			seq.con <- sub(s.pattern, '\\2', temp.sample);
			seq.separate <- unlist(strsplit(seq.con, ''));
			temp.SNP <- seq.separate[temp.index];
			data.matrix[sample.ind, iter] <- temp.SNP;
		}			
	}	
}

write.table(data.matrix, file='DyEDunlinksnppd15AllSP.txt', quote = FALSE, sep = '', col.names = FALSE);#you may want to specify the file name