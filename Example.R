library("mixor")
library("reshape2")
library("data.table")
library("mstate")
library("mvQuad")
library("readxl")
library("tidyverse")

Example_Data <- read_excel("Example_Data.xlsx")

Initial(Example_Data)
jmol(Example_Data)
