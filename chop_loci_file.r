#chop off the sites after position ### in all loci
data.loci <- scan('filename.loci', what = 'character', sep = '\n');
break.lines <- grep('//', data.loci);

for(i in 1:length(data.loci)){#the length of the .loci file
	temp <- data.loci[i];
	s.pattern <- '(.{118}).+';#this xx number need to be changed according to your need
	#it is the sequence length plus the number of space before the sequence
	r.pattern <- '\\1';
	t.txt <- sub(s.pattern, r.pattern, temp);
	data.loci[i] <- t.txt;
}
for(iter in 1:length(break.lines)){
	line.id <- break.lines[iter];
	data.loci[line.id] <- paste(data.loci[line.id], '|', sep = '');
}

write.table(data.loci, file = 'outputfilename.loci', quote = FALSE, sep = '\n', row.names = FALSE, col.names = FALSE);#change file name to your preference

