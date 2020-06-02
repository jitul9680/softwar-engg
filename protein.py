from mydb import *
import re

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

def start():

  

  with open("Ecol_K12_MG1655_.txt",'r') as my_file:
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
    
start()