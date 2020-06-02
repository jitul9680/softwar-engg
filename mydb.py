import mysql.connector


mydb=mysql.connector.connect(
	host="localhost",
	user="root",
	passwd="metasulphanol",
	database="geneinfo",
	)

my_cursor=mydb.cursor()

#my_cursor.execute("CREATE DATABASE geneinfo")
#my_cursor.execute("SHOW DATABASES")

#for db in my_cursor:
	#print(db[0])

#my_cursor.execute("CREATE TABLE gdata(Gene_Name VARCHAR(50) PRIMARY KEY,Gene_Sequence VARCHAR(8000),lenofgeneseq INTEGER(10),Amino_Acid VARCHAR(4000),lenofProtein INTEGER(10),Remarks VARCHAR(10) )")
'''
#my_cursor.execute("SHOW TABLES")
#for table in my_cursor:
#	print(table[0])
'''
sql_records="INSERT INTO gdata(Gene_Name,Gene_Sequence,lenofgeneseq,Amino_Acid,lenofProtein,Remarks) VALUES(%s,%s,%s,%s,%s,%s)"

