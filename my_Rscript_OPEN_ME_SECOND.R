# A hashtag creates a 'comment' - anything you write here won't be 'read' or processed as code 
# Comments are useful for explaining what your code does, for example:

#set the folder (directory) you want to work in 
setwd("../Documents/")

#list the files found in current directory 
list.files()

# Rstudio can easily split your code into sections that are easy to navigate
# To do this, add 4 hashes to the end of a comment, like so: 
## Code section 1 ####

## Code section 2 ####

## Creating a variable ####
#create a variable - an 'object' - using <- or =

#variables can contain different datatypes

#character: "a", "swc"
my_char <- "hello, I'm a character"

#numeric: 2, 15.5
my_num <- 2

#logical: TRUE, FALSE
my_log <- TRUE

## Data structures ####
#datatypes can be organised into different data 'structures':
# atomic vector
# list
# matrix
# data frame ( commonly we will use dataframes - think of it as a table or excel sheet containing data 
#   that can be words or numbers)
# factors

#a vector contains the same datatype, eg.:
my_vector <- c(1,2,3)
#c() concatenates things together 

#a list can contain a variety of datatypes 
my_list <- list(2,"two",FALSE)

#a matrix is like a table containing only one type of data - eg. numbers only
my_matrix <- matrix(nrow = 2, ncol = 2)

#a dataframe is like a table which can hold a variety of different datatypes (eg. numbers, text)
my_df <- data.frame(id = letters[1:10], x = 1:10, y = 11:20)

#to access a column, use $:
my_df$id
#elements within a dataframe can be accessed using their row and column index in square brackets:
# [row,column]
my_df[1,3]

#column names:
colnames(my_df)

#rownames:
rownames(my_df)

## reading data in from a file ####
metadata <- read.table(file="C:/Users/Public/Mycobiome_tutorial/metadata.csv")

## Printing a plot to file ####
# first you tell R what file you want to make
#   here we are going to make a png image file
png(filename="my_plot.png", width = 500, height = 500, bg = "white")
#once you have told R to open such a file, 
# we then need to tell it what to put in the file
# so, now we make our plot:
plot(my_df$x, my_df$y)
#then we tell R that we are finished, 
# and it will create the image file for us
dev.off()
#If you check the folder the Documents folder on your computer,
# there should now be a "my_plot.png" file with our graph inside!
