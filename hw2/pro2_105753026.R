######################################
# the reference code of program2 
######################################

######################################
# initial
######################################
library("Biostrings",verbose=F,quietly=T)

# read parameters
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("USAGE: Rscript pro2_105753026.R --input test.fasta --score pam250.txt --aln global --gap_open -10 --gap_extend -2 --output result.fasta", call.=FALSE)
}

# parse parameters
i<-1 
while(i < length(args))
{
  if(args[i] == "--input"){
    i_f<-args[i+1]
    i<-i+1
  }else if(args[i] == "--score"){
    s_f<-args[i+1]
    i<-i+1
  }else if(args[i] == "--aln"){
    aln_mode <- args[i+1]
    i<-i+1
  }else if(args[i] == "--gap_open"){
    g_o<- as.integer(args[i+1]) 
    i<-i+1
  }else if(args[i] == "--gap_extend"){
    g_e<-as.integer(args[i+1])
    i<-i+1    
  }else if(args[i] == "--output"){
    o_f<-args[i+1]
    i<-i+1
  }else{
    stop(paste("Unknown flag", args[i]), call.=FALSE)
  }
  i<-i+1
}

print("PARAMETERS")
print(paste("input file         :", i_f))
print(paste("output file        :", o_f))
print(paste("score file         :", s_f))
print(paste("gap open penalty   :", g_o))
print(paste("gap extend penalty :", g_e))

######################################
# main
######################################
# read fasta file
ff <- readAAStringSet(i_f)
seq_name = names(ff)
sequence = paste(ff)

# aln length
aln_length <- nchar(sequence[1])

# read score file
s_m<-read.table(s_f)
s_m<-as.matrix(s_m)

aln_score<-0
for(i in 1:aln_length)
{
  a<-substring(sequence[1], i, i)
  b<-substring(sequence[2], i, i)
  
  if((a != "-")&&(b != "-"))
  {
    print(paste(a, "-", b, "=", s_m[a,b]))
    aln_score = aln_score + s_m[a,b]
  }
  else{
    aln_score = aln_score + g_o
  }
}

print(aln_score)

# global alignment

# create seqence matrix of seq1 and seq2
seq1 <- unlist(strsplit(paste("*", sequence[1], sep=""), split=""))
seq2 <- unlist(strsplit(paste("*", sequence[2], sep=""), split=""))
seq_matrix <- matrix(data=NA, nrow=length(seq2), ncol=length(seq1))
dimnames(seq_matrix) <- list(seq2, seq1)

# set default score of sequence matrix
seq_matrix[1,1] = 0

for(i in 2:length(seq1))
  seq_matrix[1,i] = as.numeric(g_o) + (i-2)*as.numeric(g_e)
for(j in 2:length(seq2))
  seq_matrix[j,1] = as.numeric(g_o) + (j-2)*as.numeric(g_e)

# calculate sequence matrix value
gap_open <- TRUE

for(i in 2:length(seq1)) {
  for(j in 2:length(seq2)) {
    match <- s_m[seq1[i],seq2[j]] + seq_matrix[j-1,i-1]
    if(gap_open == TRUE) {
      gap1 <- seq_matrix[j-1,i] + as.numeric(g_o)
      gap2 <- seq_matrix[j,i-1] + as.numeric(g_o)
      seq_matrix[j,i] <- max(match, gap1, gap2)
      if(seq_matrix[j,i] != match)
        gap_open = FALSE
    }
    else {
      gap1 <- seq_matrix[j-1,i] + as.numeric(g_e)
      gap2 <- seq_matrix[j,i-1] + as.numeric(g_e)
      seq_matrix[j,i] <- max(match, gap1, gap2)
      if(seq_matrix[j,i] == match)
        gap_open = TRUE
    }
  }
}

# trace back path

count1 = count2 = 2
aln_seq1 = aln_seq2 = ""

while(count1 <= length(seq1) || count2 <= length(seq2)) {
  if(count1 == length(seq1)) {
    aln_seq1 <- paste(aln_seq1, "-", sep="")
    aln_seq2 <- paste(aln_seq2, seq2[count2], sep="")
    count2 = count2+1
    if(count2 > length(seq2))
      break
  }
  else if(count2 == length(seq2)) {
    aln_seq1 <- paste(aln_seq1, seq1[count1], sep="")
    aln_seq2 <- paste(aln_seq2, "-", sep="")
    count1 = count1+1
    if(count1 > length(seq1))
      break
  }
  else {
    aln_max = max(seq_matrix[count2,count1+1], seq_matrix[count2+1,count1], seq_matrix[count2+1, count1+1])
    if(aln_max == seq_matrix[count2+1,count1+1]) {
      aln_seq1 <- paste(aln_seq1, seq1[count1], sep="")
      aln_seq2 <- paste(aln_seq2, seq2[count2], sep="")
      count1 = count1+1
      count2 = count2+1
    }
    else if(aln_max == seq_matrix[count2,count1+1]) {
      aln_seq1 <- paste(aln_seq1, seq1[count1], sep="")
      aln_seq2 <- paste(aln_seq2, "-", sep="")
      count1 = count1+1
    }
    else if(aln_max == seq_matrix[count2+1,count1]) {
      aln_seq1 <- paste(aln_seq1, "-", sep="")
      aln_seq2 <- paste(aln_seq2, seq2[count2], sep="")
      count2 = count2+1
    }
  }
}

# output
writeXStringSet(ff, o_f)
write(aln_seq1, file=o_f, append=FALSE)
write(aln_seq2, file=o_f, append=TRUE)
print(aln_seq1)
print(aln_seq2)
