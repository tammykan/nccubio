# read PAM1 from data
pam1 <- read.table("pam1.txt", header=TRUE)

# check PAM1 data
dim(pam1)
str(pam1)

# construct PAM250 from PAM1
#library('expm')
pam1 <- as.matrix(pam1)/10000
#pam250 <- as.matrix(pam1) %^% 250

pam250 <- as.matrix(pam1)

for(x in 1:250) {
  pam250 <- pam250 %*% pam1
}

pam250 <- as.matrix(pam250)*100
pam250 <- round(pam250, 0)

# output PAM250 as a file
data <- iris
write.table(pam250, file="pam250.txt", sep=",")
