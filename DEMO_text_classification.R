library(data.table)
library(Matrix)
library(xgboost)
library(caret)
library(stringr)
library(tm)
library(syuzhet) 
library(bigmemory)
library(tibble)
library(tidyr)
library(dplyr)
library(readr)
memory.limit(40000)
options(java.parameters = "-Xmx15000m")
library(rJava)


# [0] initial setting & functions -------------------------------------------------------------
# mac
file_path <- "/Users/bee/Google Drive/R/Projects/mykaggle/"
# windows
file_path <- "D:/Google Drive/R/Projects/mykaggle/"


# LabelCount Encoding function
labelCountEncoding <- function(column){
  return(match(column,levels(column)[order(summary(column,maxsum=nlevels(column)))]))
}

# a matrix which is zero except for the column corresponding to the class.
#The function is currently defined as
class.ind <- function(cl, cap)
{
  n <- length(cl)
  cl <- as.factor(cl)
  x <- matrix(0, n, length(levels(cl)) )
  x[(1:n) + n*(unclass(cl)-1)] <- 1
  dimnames(x) <- list(names(cl), paste0(cap , levels(cl)) )
  x
}

# [1] Data pre-processing -------------------------------------------------------------

# Load CSV files 
cat("Read data")

# read training_variants
train_var <- read.table(file = paste0(file_path,"data/training_variants"), header = T, sep = ",")

# read training_text
train_txt_dump <- tibble(text = read_lines(paste0(file_path,"data/training_text"),skip = 1))
train_txt <- train_txt_dump %>%
  separate(text, into = c("ID", "txt"), sep = "\\|\\|")
train_txt <- train_txt %>%
  mutate(ID = as.integer(ID))

# read test_variants
test_var <- read.table(file = paste0(file_path,"data/test_variants"), header = T, sep = ",")

# read testing_text
test_txt_dump <- tibble(text = read_lines(paste0(file_path,"data/test_text"),skip = 1))
test_txt <- test_txt_dump %>%
  separate(text, into = c("ID", "txt"), sep = "\\|\\|")
test_txt <- test_txt %>%
  mutate(ID = as.integer(ID))

# remove useless data
train <- merge(train_var, train_txt, by="ID")
test <- merge(test_var, test_txt, by="ID")
rm(train_txt, train_txt_dump, train_var, test_txt, test_txt_dump, test_var);gc()

test$Class <- -1
data <- rbind(train,test) %>% as.data.table()
rm(train,test);gc()
names(data) <- c("ID", "Gene", "Variation", "Class", "Text")

# Basic text features
cat("Basic text features")
data$nchar <- as.numeric(nchar(data$Text))
data$nwords <- as.numeric(str_count(data$Text, "\\S+"))

# LabelCount Encoding for Gene and Variation 
# We can do more advanced feature engineering later, e.g. char-level n-grams
#data$Gene <- labelCountEncoding(data$Gene)
#data$Variation <- labelCountEncoding(data$Variation)
data <- data %>% 
  extract(Variation, into=c("First_Letter", "Var_2"), 
          regex = "^(?=.{1,7}$)([a-zA-Z]+)([0-9].*)$", 
          remove = FALSE)

# Split number and last letter for typical variations
data <- data %>% 
  extract(Var_2, into=c("Gene_Location", "Last_Letter"),
          regex = "^([0-9]+)([a-zA-Z]|.*)$",
          remove = TRUE)


data$Gene_Location[is.na(data$Gene_Location)] <- "0"
data$Gene_Location <- as.numeric(data$Gene_Location)
data[is.na(data)] <- "IamNA"
data$First_Letter <- as.factor(data$First_Letter)
data$Last_Letter <- as.factor(data$Last_Letter)

# factor to dummy coding (First & last letter)
# data <- as.data.frame(data)
# data <-  cbind.data.frame(data[,-c(4,6)], 
#                           class.ind(data$First_Letter, "First_"),
#                           class.ind(data$Last_Letter , "Last_") 
#                           )

# Compute Text length - nchar & nwords
data$nchar <- as.numeric(nchar(data$Text))
data$nwords <- as.numeric(str_count(data$Text, "\\S+"))

#labelCountEncoding(data$Gene)
#match(data$Gene,levels(data$Gene)[order(summary(data$Gene,maxsum=nlevels(data$Gene)))])


# Identify and encode deletions
data$is_del <- ifelse(grepl("del", data$Variation, ignore.case = T), 1, 0) 

# Identify and encode insertions
data$is_ins <- ifelse(grepl("ins", data$Variation, ignore.case = T), 1, 0)

# Identify and encode fusions
data$is_fus <- ifelse(grepl("fus", data$Variation, ignore.case = T), 1, 0) 

# Identify and encode truncation
data$is_trunc <- ifelse(grepl("trunc", data$Variation, ignore.case = T), 1, 0) 

# Identify and encode methylations
data$is_methyl <- ifelse(grepl("methyl", data$Variation, ignore.case = T), 1, 0) 

# Identify and encode amplifications
data$is_amp <- ifelse(grepl("amp", data$Variation, ignore.case = T), 1, 0) 

# Identify and encode silencing
data$is_sil <- ifelse(grepl("sil", data$Variation, ignore.case = T), 1, 0) 

# Identify and encode overexpression
data$is_expr <- ifelse(grepl("expr", data$Variation, ignore.case = T), 1, 0) 

# Identify and encode splicing
data$is_splice <- ifelse(grepl("splice", data$Variation, ignore.case = T), 1, 0) 

# Identify and encode exon variations
data$is_exon <- ifelse(grepl("exon", data$Variation, ignore.case = T), 1, 0) 

# Identify and encode dup variations
data$is_dup <- ifelse(grepl("dup", data$Variation, ignore.case = T), 1, 0) 

# Identify and encode Fusion variations
data$is_fusion <- ifelse(grepl("Fusion", data$Variation, ignore.case = T), 1, 0)

# Identify and encode BRAF variations
data$is_BRAF <- ifelse(grepl("BRAF", data$Variation, ignore.case = T), 1, 0) 

# Identify and encode BCOR variations
data$is_BCOR <- ifelse(grepl("BCOR", data$Variation, ignore.case = T), 1, 0) 

# Identify and encode BCR variations
data$is_BCR <- ifelse(grepl("BCR", data$Variation, ignore.case = T), 1, 0) 

# Identify and encode Mutation variations
data$is_Mutation <- ifelse(grepl("Mutation", data$Variation, ignore.case = T), 1, 0) 

# Identify and encode * variations
data$is_star <- ifelse(grepl("*", data$Variation, ignore.case = T), 1, 0) 

# Identify and encode Deletion variations
data$is_Deletion <- ifelse(grepl("Deletion", data$Variation, ignore.case = T), 1, 0) 

# Identify and encode DNA variations
data$is_DNA <- ifelse(grepl("DNA", data$Variation, ignore.case = T), 1, 0) 

# Identify and encode ESR variations
data$is_ESR <- ifelse(grepl("ESR", data$Variation, ignore.case = T), 1, 0) 

# Identify and encode EWSR variations
data$is_EWSR <- ifelse(grepl("EWSR", data$Variation, ignore.case = T), 1, 0) 

# Identify and encode Exon variations
data$is_Exon <- ifelse(grepl("Exon", data$Variation, ignore.case = T), 1, 0) 

# [2] feature extraction & selection -------------------------------------------------------------

# [2-1] TF-IDF---------------------------------------------------- 
cat("TF-IDF")
txt <- Corpus(VectorSource(data$Text))
txt <- tm_map(txt, stripWhitespace)
txt <- tm_map(txt, content_transformer(tolower))
txt <- tm_map(txt, removePunctuation)
txt <- tm_map(txt, removeWords, stopwords("english"))
txt <- tm_map(txt, stemDocument, language="english")
txt <- tm_map(txt, removeNumbers)

#dtm <- DocumentTermMatrix(txt, control = list(weighting = weightTfIdf))
#dtm <- DocumentTermMatrix(txt)
dtm <- DocumentTermMatrix(txt, control = list(weighting = function(x) weightSMART(x, spec = "ntc")))

removeSparseTerms(dtm, 0.99) %>% dim()
dtm <- removeSparseTerms(dtm, 0.99)
data <- cbind(data, as.matrix(dtm))
data_backup <- data
data <- data_backup

# [2-2] Information Gain, IG--------------------------------------------
library(FSelector)
data$Class <- data$Class %>% as.character() %>% as.factor()
# error message: duplicated name 'naiv' in data frame using '.'
# grep("naiv",names(data) ) -> 3711
tmp_data <- data[c(1:3321), ] %>% as.data.frame()
#tmp_data <- tmp_data[ , -c(1:6, 8:18,  grep("naiv",names(data) ))]
tmp_data <- tmp_data[ , -c(1:4, 6, 
                           grep("accord",names(data)), 
                           grep("naiv",names(data)),
                           grep("activ",names(data))
)]

#IG <- information.gain(Class~. , data=tmp_data, unit = "log2")
IG <- information.gain(Class~. , data=tmp_data)

# pick IG top 500 words
IG_top500 <- head(IG %>% mutate(index=1:nrow(IG))
                  %>% arrange(desc(attr_importance)), 1000)
# make some plot
plot(IG_top500)

IG %>% 
  arrange(desc(attr_importance)) %>% 
  mutate(index=1:nrow(IG)) %>% 
  plot

# combine to new data
data <- as.data.frame(data)[ , -c(grep("accord",names(data)), 
                                  grep("naiv",names(data)),
                                  grep("activ",names(data)))]
#data <- data[ ,c(1:18,(IG_top500$index+18))] %>% as.data.table 
data <- data[ ,c(1:6,(IG_top500$index+6))] %>% as.data.table 

# [3] XGB model construction -------------------------------------------------------------

# XG boost 
# Sentiment analysis
#cat("Sentiment analysis")
#sentiment <- get_nrc_sentiment(data$Text) 
#data <- cbind(data,sentiment) 

# Set seed
set.seed(123456) 

# cross validation setting
data$Class <- data$Class %>% as.character() %>% as.integer()
data <- as.data.table(data)
cvFoldsList <- createFolds(data$Class[data$Class > -1], k=5, list=TRUE, returnTrain=FALSE)

# To sparse matrix
cat("Create sparse matrix")
varnames <- setdiff(colnames(data), c("ID", "Class","Variation", "Text"))
train_sparse <- Matrix(as.matrix(sapply(data[Class > -1, varnames, with=FALSE],as.numeric)), sparse=TRUE)
test_sparse <- Matrix(as.matrix(sapply(data[Class == -1, varnames, with=FALSE],as.numeric)), sparse=TRUE)
y_train <- data[Class > -1,Class]-1
test_ids <- data[Class == -1,ID]
dtrain <- xgb.DMatrix(data=train_sparse, label=y_train)
dtest <- xgb.DMatrix(data=test_sparse)

# Params for xgboost
param <- list(booster = "gbtree",
              objective = "multi:softprob", 
              eval_metric = "mlogloss",
              num_class = 9,
              eta = .2,
              gamma = 1,
              max_depth = 6, 
              min_child_weight = 1,
              subsample = 0.7, 
              colsample_bytree = 0.7 
)

# Cross validation - determine CV scores & optimal amount of rounds
cat("XGB cross validation")
xgb_cv <- xgb.cv(data = dtrain,
                 params = param,
                 nrounds = 3000,
                 maximize = FALSE,
                 prediction = TRUE,
                 folds = cvFoldsList,
                 print.every.n = 10,
                 early.stop.round = 100
)
rounds <- xgb_cv$best_iteration
#rounds <- which.min(xgb_cv$dt[, test.mlogloss.mean])

# Train model
cat("XGB training")
xgb_model <- xgb.train(data = dtrain,
                       params = param,
                       watchlist = list(train = dtrain),
                       nrounds = rounds,
                       verbose = 1,
                       print.every.n = 20
)

# [4] revise model, and repeat step 3 -------------------------------------------------------------

# Feature importance
cat("Plotting feature importance")
names <- dimnames(train_sparse)[[2]]
importance_matrix <- xgb.importance(names, model=xgb_model)
xgb.plot.importance(importance_matrix = importance_matrix, top_n = 1000)

# top 1000 importance
imp <- importance_matrix$Feature[1:1000]
imp_index <- sapply((1:length(imp)), function(i) which( names(data) == imp[i]) ) 

# imp_index
data <- data %>% as.data.frame() 
data <- data[ ,c(1, grep("Class",names(data)), imp_index)]  %>% as.data.table 


# [5] Output -------------------------------------------------------------

# Predict and output csv
cat("Predictions")
preds <- as.data.table(t(matrix(predict(xgb_model, dtest), nrow=9, ncol=nrow(dtest))))
colnames(preds) <- c("class1","class2","class3","class4","class5","class6","class7","class8","class9")

write.table(data.table(ID=test_ids, preds), paste0(file_path,"output.csv"), sep=",", dec=".", quote=FALSE, row.names=FALSE)
