from mydb import *
from collections import Counter

def acid_analysis():
	
	my_cursor.execute("SELECT Amino_Acid FROM gdata where Remarks='Good'")
	my_result=my_cursor.fetchall()
	val=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
	for row in my_result:
		#print(str(row[0])[:-4])
		d=Counter(list(str(row[0])[:-4]))
		for key,value in d.items():
			#print(key,value)
			if key=="A":
				val[0]=val[0]+d[key]
			elif key=="C":
				val[1]=val[1]+d[key]
			elif key=="D":
				val[2]=val[2]+d[key]
			elif key=="E":
				val[3]=val[3]+d[key]
			elif key=="F":
				val[4]=val[4]+d[key]
			elif key=="G":
				val[5]=val[5]+d[key]
			elif key=="H":
				val[6]=val[6]+d[key]
			elif key=="I":
				val[7]=val[7]+d[key]
			elif key=="K":
				val[8]=val[8]+d[key]
			elif key=="L":
				val[9]=val[9]+d[key]
			elif key=="M":
				val[10]=val[10]+d[key]
			elif key=="N":
				val[11]=val[11]+d[key]
			elif key=="P":
				val[12]=val[12]+d[key]
			elif key=="Q":
				val[13]=val[13]+d[key]
			elif key=="R":
				val[14]=val[14]+d[key]
			elif key=="S":
				val[15]=val[15]+d[key]
			elif key=="T":
				val[16]=val[16]+d[key]
			elif key=="V":
				val[17]=val[17]+d[key]
			elif key=="W":
				val[18]=val[18]+d[key]
			elif key=="Y":
				val[19]=val[19]+d[key]


			
				

		key_max=max(d.keys(),key=(lambda k:d[k]))
		key_min=min(d.keys(),key=(lambda k:d[k]))
		#print("Most frequent Amino Acid ={} {}".format(key_max,d[key_max]))
		#print("Least frequent Amino Acid = {} {}".format(key_min,d[key_min]))
		
	entry_into_db="INSERT INTO amino_acids(Alanine,Cysteine,Aspartic_acid,Glutamic_acid,Phenylalanine,Glycine,Histidine,Isoleucine,Lysine,Leucine,Methionine,Asparagine,Proline,Glutamine,Arginine,Serine,Threonine,Valine,Tryptophan,Tyrosine) VALUES(%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)"
	my_cursor.execute(entry_into_db,(val[0],val[1],val[2],val[3],val[4],val[5],val[6],val[7],val[8],val[9],val[10],val[11],val[12],val[13],val[14],val[15],val[16],val[17],val[18],val[19]))
	mydb.commit()