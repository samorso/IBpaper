# Show how we build the data from "student-mat.csv"
student_math <- read.table("data-raw/student-mat.csv", sep=",", header=TRUE)

# Construct the response variable (0 if very low alcohol on workday and low on weekend, 1 otherwise)
y <- ifelse(student_math$Walc <= 2 & student_math$Dalc == 1, 0, 1)

# Construct the design matrix
# Subset original dataset
var_names <- colnames(student_math)
var_names <- var_names[-c(7,8,27,28,29)] # remove Medu, Fedu, Dalc, Walc and health
tmp <- subset(student_math, select=var_names)

# Re-level some factors
tmp$traveltime <- factor(student_math$traveltime, levels = c("1","2","3","4"), labels = c("1","2","3","3")) # too few "4"
tmp$studytime <- as.factor(student_math$studytime)
tmp$freetime <- as.factor(student_math$freetime)
tmp$failures <- factor(student_math$failures, levels = c("0","1","2","3"), labels = c("0","1","2","2"))
tmp$famrel <- factor(student_math$famrel, levels = c("1","2","3","4","5"), labels = c("2","2","3","4","5"))
tmp$goout <- as.factor(student_math$goout)

# Scale numerical variables
# tmp$age <- scale(student_math$age)
# tmp$absences <- scale(student_math$absences)
# tmp$G1 <- scale(student_math$G1)
# tmp$G2 <- scale(student_math$G2)
# tmp$G3 <- scale(student_math$G3)

# Save the data frame
student <- as.data.frame(cbind(y,tmp))
colnames(student)[1] <- "alc"
usethis::use_data(student, overwrite = TRUE)
