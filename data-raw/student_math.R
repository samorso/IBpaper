# Show how we build the data from "student-mat.csv"
student_math <- read.table("data-raw/student-mat.csv", sep=",", header=TRUE)

# Construct the response variable (0 if very low alcohol on workday and low on weekend, 1 otherwise)
y <- ifelse(student_math$Walc <= 2 & student_math$Dalc == 1, 0, 1)

# Construct the design matrix
# Subset original dataset
var_names <- colnames(student_math)
var_names <- var_names[-c(7,8,27,28,29)] # remove Medu, Fedu, Dalc, Walc and health
tmp <- subset(student_math, select=var_names)

# Define and re-level some factors
tmp$traveltime <- factor(student_math$traveltime, levels = c("1","2","3","4"), labels = c("1","2","3","3")) # too few "4"
tmp$studytime <- as.factor(student_math$studytime)
tmp$failures <- factor(student_math$failures, levels = c("0","1","2","3"), labels = c("0","1","2","2"))
tmp$freetime <- as.factor(student_math$freetime)
tmp$Mjob[tmp$Mjob == "other"] <- "aaa" # rename to make "other" the first contrast
tmp$Fjob[tmp$Fjob == "other"] <- "aaa" # rename to make "other" the first contrast
tmp$reason[tmp$reason == "other"] <- "aaa" # rename to make "other" the first contrast
tmp$guardian[tmp$guardian == "other"] <- "aaa" # rename to make "other" the first contrast

# Scale numerical variables
tmp$age <- scale(student_math$age)
tmp$famrel <- scale(student_math$famrel)
tmp$absences <- scale(student_math$absences)
tmp$goout <- scale(student_math$goout)
tmp$G1 <- scale(student_math$G1)
tmp$G2 <- scale(student_math$G2)
tmp$G3 <- scale(student_math$G3)

# Define the design matrix
x <- model.matrix(~., tmp)
x <- x[,!grepl("Intercept",colnames(x))] # remove intercept

# Save the data frame
student <- as.data.frame(cbind(y,x))
colnames(student)[1] <- "alc"
usethis::use_data(student, overwrite = TRUE)
