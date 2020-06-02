import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
def start_program():
	with open('Ecol_K12_MG1655_.wgs') as fp:
		ecoli_gene_sequence=''
		for i in fp:
			if (i[0]!='>'):
				ecoli_gene_sequence=ecoli_gene_sequence+str(i[:-1])

		parser_ecoli_gene_sequence(ecoli_gene_sequence)
		#print(gc_list)




	#print(len(ecoli_gene_sequence))
def parser_ecoli_gene_sequence(a):
	count=0
	gc_list=[]
	ch=""
	values=[]
	j=1
	for i in range(len(a)):
		ch=ch+a[i]
		count=count+1
		if (count==99):
			values=counting_the_no_of_atgc(ch)
			gc_list.append((values[2]-values[3])/(values[3]+values[2]))
			#print('fragment_no{} ---------->count of a,t,g,c respectively are {} {} {} {}'.format(j,values[0],values[1],values[2],values[3]))
			j=j+1
			count=0
			ch=''
		else:
			continue
	gc_cummulative=np.cumsum(gc_list)
	#print(gc_list[0],gc_list[1],gc_list[2],gc_list[3],gc_list[4],gc_list[5],gc_list[6])
	#print(gc_cummulative)
	print(len(gc_cummulative))
	c=max(gc_cummulative)
	d=min(gc_cummulative)
	for x,y in enumerate(gc_cummulative):
		if y==c:
			print("fragment no:{} and peak-value={}".format(x,y))
			co=x
		if y==d:
			print("fragment no:{} and peak-value={}".format(x,y))
			c01=x
	print('negative percentage : {}'.format((c01-co)*100/(len(gc_cummulative))))
	print('positive percentage : {}'.format(100-(c01-co)*100/(len(gc_cummulative))))
	plt.plot(range(1,46866),gc_cummulative,color='green')
	plt.xlabel('fragments')
	plt.ylabel('cummulative sum')
	plt.title('gc-skew diagram')
	plt.show()

	



            
    
def counting_the_no_of_atgc(parsed_sequence):

	l=[0,0,0,0]

	for x in range(len(parsed_sequence)):
		if parsed_sequence[x]=='a':
			l[0]=l[0]+1;
		if parsed_sequence[x]=='t':
			l[1]=l[1]+1;
		if parsed_sequence[x]=='g':
			l[2]=l[2]+1;
		if parsed_sequence[x]=='c':
			l[3]=l[3]+1;

	return (l)

start_program()