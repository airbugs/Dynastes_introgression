#chop off the 1st site and the sites after ### in all loci
data.loci <- scan('GLAntilles.loci', what = 'character', sep = '\n');
break.lines <- grep('//', data.loci);

for(i in 1:length(data.loci)){#the length of the .loci file
	temp <- data.loci[i];
	s.pattern <- '(.{148}).+';#this xx number need to be changed according to your need
	r.pattern <- '\\1';
	t.txt <- sub(s.pattern, r.pattern, temp);
	data.loci[i] <- t.txt;
}
for(iter in 1:length(break.lines)){
	line.id <- break.lines[iter];
	data.loci[line.id] <- paste(data.loci[line.id], '|', sep = '');
}

write.table(data.loci, file = 'GAntilles130.loci', quote = FALSE, sep = '\n', row.names = FALSE, col.names = FALSE);#change file name to your preference




#below is a old r scrip. slower than the upper one.
#number <- length(data.loci);
#newdata.loci <- NULL;
#for(i in 1:number){#the length of the .loci file
# 	temp <- data.loci[i];
# 	temp.split <- unlist(strsplit(temp, ''));
#	ed.txt <- temp.split[1:129];#xx is the length of seq plus the name of the seq and spaces
#	new.txt <- paste(ed.txt, collapse = '');
#	newdata.loci <- c(newdata.loci, new.txt);
#}

#write.table(newdata.loci, file = 'chop110.loci', quote = FALSE, sep = '\n', row.names = FALSE, col.names = FALSE);