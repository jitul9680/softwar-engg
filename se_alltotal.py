from mydb import *
import re
from collections import Counter
import matplotlib.pyplot as plt
from scipy import signal
import numpy as np
from tkinter import *
from tkinter import filedialog
root=Tk()


def parser(a):
    count=0
    ch=""
    protein=""
    protein_table={'ttt':'F','tct':'S','tat':'Y','tgt':'C',
                   'ttc':"F",'tcc':'S','tac':'Y','tgc':'C',
                   'tta':'L',"tca":'S',"taa":'STOP','tga':"STOP",
                   'ttg':'L','tcg':'S','tag':'STOP','tgg':"W",
                   'ctt':'L',"cct":'P','cat':'H','cgt':'R',
                   "ctc":'L','ccc':'P','cac':'H','cgc':'R',
                   'cta':'L','cca':'P','caa':'Q','cga':'R',
                   'ctg':'L','ccg':'P','cag':'Q','cgg':'R',
                    'att':'I','act':'T','aat':'N','agt':'S',
                    'atc':'I','acc':'T','aac':'N','agc':'S',
                    'ata':'I','aca':'T','aaa':'K','aga':'R',
                    'atg':'M','acg':'T','aag':'K','agg':'R',
                    'gtt':'V','gct':'A','gat':'D','ggt':'G',
                    'gtc':'V','gcc':'A','gac':'D','ggc':'G',
                    'gta':'V','gca':'A','gaa':'E','gga':'G',
                    'gtg':'V','gcg':'A','gag':'E','ggg':'G'
                  }
    for i in range(len(a)):
        ch=ch+a[i]
        count=count+1
        if count==3:
            
            protein=protein+protein_table[ch]
            count=0
            ch=''
            
        else:
            
            continue
            
    return protein

def start(x):

  

  with open(x,'r') as my_file:
    genename=''
    gene_name=''
    geneseq=''
    protein=''
    for i in my_file:
        if(i[0]=='>'):
          genename=str(i[:-1])
          #print(genen)
          if(geneseq!=''):

            if len(geneseq)%3==0:
              protein=parser(geneseq)
              pattern="STOP"
              match=re.findall(pattern,protein)
              if (match.count("STOP")!=2):
                #print(genename,'\n',geneseq,'\n',len(geneseq),'\n',protein,'\n',len(protein[:-3]),"Good")
                my_cursor.execute(sql_records,(gene_name,geneseq,len(geneseq),protein,len(protein[:-3]),"Good"))
              else:
                #print(genename,'\n',geneseq,'\n',len(geneseq),'\n',protein,'\n',len(protein[:-3]),"Good")
                my_cursor.execute(sql_records,(gene_name,geneseq,len(geneseq),protein,len(protein[:-3]),"BAD"))



              #print(genename,'\n',geneseq,'\n',len(geneseq),'\n',protein,'\n',len(protein[:-3]),"Good")
              
            else:
              protein=parser(geneseq)
              #print(genename,'\n',geneseq,'\n',len(geneseq),'\n',protein,'\n',len(protein[:-3]),"BAD")
              my_cursor.execute(sql_records,(gene_name,geneseq,len(geneseq),protein,len(protein[:-3]),"BAD"))


          geneseq=''
          protein=''
        else:
          geneseq=geneseq+str(i[:-1])
          gene_name=genename;


    #print(genename,'\n',geneseq,'\n',len(geneseq),'\n',parser(geneseq),len(parser(geneseq)[:-3]),"Good")
    #my_cursor.execute(sql_records,(genename,geneseq,len(geneseq),parser(geneseq),len(parser(geneseq)[:-3]),"Good"))
    mydb.commit()


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
  names_of_aminoacids=["Alanine","Cysteine","Aspartic_acid  ","Glutamic_acid  ","Phenylalanine  ","Glycine","Histidine","Isoleucine","Lysine","Leucine","Methionine  ","Asparagine","Proline","Glutamine","Arginine","Serine","Threonine","Valine","Tryptophan","Tyrosine"]
  numeric_values=[val[0]*100/sum(val),val[1]*100/sum(val),val[2]*100/sum(val),val[3]*100/sum(val),val[4]*100/sum(val),val[5]*100/sum(val),val[6]*100/sum(val),val[7]*100/sum(val),val[8]*100/sum(val),val[9]*100/sum(val),val[10]*100/sum(val),val[11]*100/sum(val),val[12]*100/sum(val),val[13]*100/sum(val),val[14]*100/sum(val),val[15]*100/sum(val),val[16]*100/sum(val),val[17]*100/sum(val),val[18]*100/sum(val),val[19]*100/sum(val)]
  plt.bar(names_of_aminoacids,numeric_values)
  plt.xlabel("****************Amino Acids********************")
  plt.ylabel("****************Quantities***********************")
  #plt.setp(plt.xticklabels(), rotation=30, horizontalalignment='right')
  plt.xticks(rotation=45)
  plt.show()


def show_databases():

  my_cursor.execute("select * from gdata")
  result=my_cursor.fetchall();

  for i in result:
    print(i)

def start_program(x):
  with open(x) as fp:
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




    
print("************************************************menu***********************************")


while(True):
  print("welcome, \n 1. Upload the file \n 2. show database \n 3. show histogram \n 4. show skew-graph \n 5. Help")

  choice=input()


  if choice=='1':
    

    root.filename=filedialog.askopenfilename(filetypes=(("png files",".png"),("all files","*.*")))
    start(root.filename)

  if choice=='2':
    show_databases()

  if choice=='3':
    acid_analysis()

  if choice=='4':
    root.filename=filedialog.askopenfilename(filetypes=(("png files",".png"),("all files","*.*")))
    start_program(root.filename)

  if choice=='5':
    print("NOTE**:if the desired file is not in the same directory then you have to enter the location of the file")
    print("submitted by \n CSB17072 \n CSB17060 \n CSB17056 \n CSB17018")
