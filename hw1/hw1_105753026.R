# read PAM1 from data
pam1 <- read.table(file.choose(), header=TRUE)

# check PAM1 data
dim(pam1)
str(pam1)

# construct PAM250 from PAM1
library('expm')
pam1 <- as.matrix(pam1)/10000
pam250 <- as.matrix(pam1) %^% 250
pam250 <- as.matrix(pam250)*100


# output PAM250 as a file
data <- iris
write.table(data, file="pam250.txt")
