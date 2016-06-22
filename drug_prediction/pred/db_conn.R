library(RMySQL)
mysql_drvr <-dbDriver("MySQL")
db_con <- dbConnect(mysql_drvr,user="app_proj_reposit",password="4repositioning",host="buttelab-db1",dbname="proj_repositioning")

#nosology_con <- dbConnect(mysql_drvr,user="app_proj_nosolog",password="4nosology",host="buttelab-db1",dbname="proj_nosology")

