from Bio import SeqIO
from collections import defaultdict
import collections
import numpy as np
import math

dicodonList = []
aaDict = {
    'ttt': 'Phe', 'tct': 'Ser', 'tat': 'Tyr', 'tgt': 'Cys',
    'ttc': 'Phe', 'tcc': 'Ser', 'tac': 'Tyr', 'tgc': 'Cys',
    'tta': 'Leu', 'tca': 'Ser', 'taa': '*', 'tga': '*',
    'ttg': 'Leu', 'tcg': 'Ser', 'tag': '*', 'tgg': 'Trp',
    'ctt': 'Leu', 'cct': 'Pro', 'cat': 'His', 'cgt': 'Arg',
    'ctc': 'Leu', 'ccc': 'Pro', 'cac': 'His', 'cgc': 'Arg',
    'cta': 'Leu', 'cca': 'Pro', 'caa': 'Gln', 'cga': 'Arg',
    'ctg': 'Leu', 'ccg': 'Pro', 'cag': 'Gln', 'cgg': 'Arg',  #dictionary kuri sutapatina kodonus su jais atitinkančiomis amino rūgštimis
    'att': 'Ile', 'act': 'Thr', 'aat': 'Asn', 'agt': 'Ser',
    'atc': 'Ile', 'acc': 'Thr', 'aac': 'Asn', 'agc': 'Ser',
    'ata': 'Ile', 'aca': 'Thr', 'aaa': 'Lys', 'aga': 'Arg',
    'atg': 'Met', 'acg': 'Thr', 'aag': 'Lys', 'agg': 'Arg', 
    'gtt': 'Val', 'gct': 'Ala', 'gat': 'Asp', 'ggt': 'Gly',
    'gtc': 'Val', 'gcc': 'Ala', 'gac': 'Asp', 'ggc': 'Gly',
    'gta': 'Val', 'gca': 'Ala', 'gaa': 'Glu', 'gga': 'Gly',
    'gtg': 'Val', 'gcg': 'Ala', 'gag': 'Glu', 'ggg': 'Gly'
   }
orderedNames = ['Phe','Leu','Ile','Met','Val','Ser','Pro','Thr','Ala','Tyr','*','His','Gln','Asn','Lys','Asp','Glu','Cys','Trp','Arg','Gly']

#visų amino rugščių sąrašas

s={'ttt','ttc','tta','ttg','ctt','ctc','cta','ctg','att','atc','ata','atg','gtt',
   'gtc','gta','gtg','tct','tcc','tca','tcg','cct','ccc','cca','ccg','act','acc', #visų kodonų sąrašas
   'aca','acg','gct','gcc','gca','gcg','tat','tac','taa','tag','cat','cac','caa',
   'cag','aat','aac','aaa','aag','gat','gac','gaa','gag','tgt','tgc','tga',
   'tgg','cgt','cgc','cga','cgg','agt','agc','aga','agg','ggt','ggc','gga','ggg'}

dicodonArray=[];
i=0
for a in s:
 for b in s:
     c=(a+b)
     dicodonArray.append(c)








def calculate(filename):
 dicodonDict=dict.fromkeys(dicodonArray,0)
 fname=filename
 virus_name=""
 for seq_record in SeqIO.parse(fname, "fasta"):
  print(seq_record.id)
  print(repr(seq_record.seq))
  print(len(seq_record))
  virus_name=repr(str(seq_record.name))
 n = 3
 out = [(seq_record.seq[i:i+n]) for i in range(0, len(seq_record.seq), n)] 
 length = len(out) 
 i = 0
 fragmentcount=0
 f= open(fname+"codonresults.txt","w")
 finalstring=""
# Suskirstom seką į start-stop kodonų poras
 while i < length:
     if out[i]=="ATG":
        currentsequence=""
        j=i
        while j<length:
            currentsequence+=out[j]
            if out[j]=="TAG" or out[j]=="TAA" or out[j]=="TGA":
                if len(str(currentsequence))>=100: #atfiltruojam visas sekas trumpesnes už 100
                 finalstring += currentsequence
                 f.write(str(currentsequence))
                 f.write("\n")
                i=j
                break
            j+=1
     i += 1
 f.close()
 finalstring=finalstring.lower()
 codons = [finalstring[i:i + 3] for i in range(0, len(finalstring), 3)]
 dicodons=[str(finalstring[i:i + 6]) for i in range(0, len(finalstring)-3, 3)]
 codon_count = collections.Counter(aaDict[c] for c in codons) #skaičiuojam visų  tipų kodonų skaičių


 dicodon_count= collections.Counter(dicodons) #skaičiuojam visų rušių dikodonų skaičių
 updatedDict=dicodonDict
 updatedDict.update(dicodon_count) #sujungiam musų tuščią dikodonų dictionary, kur visi values yra 0 su dikodonų dict kurioje yra kiekvieno dikodono skaičius.
 
 

 codon_number=len(codons)
 dicodon_number=len(dicodons)
 for key, value in updatedDict.items():
    updatedDict[key] = value/dicodon_number #dikodonų skaičių pakeičiame su jo dažniu toje sekoje.

 dicodonList.append(updatedDict)
 f=open("codon frequency.txt","a+")
 dataMatrix = np.array([round(codon_count[i]/codon_number,4) for i in orderedNames]) #surašome kodonų pasikartojimo dažnius į txt failą
 for i in range(0, len(dataMatrix)):
  f.write(str(dataMatrix[i]))
  f.write(" ")
 f.write("\n")
 f.close()

 
 f=open(fname+" "+virus_name+" codon frequency.txt","w")          
 for ele in codon_count:
  f.write(ele+" ")                                 #suskaičiuojame ir surašome kiekvieno kodono dažnį tam tikro viruso DNR.
  f.write(str(round(codon_count[ele]/codon_number,4)))
  f.write("\n")
 f.close() 
 f=open(fname+" "+virus_name+" codon count.txt","w") #surašom suskaičioutą info i txt failą
 f.write(str(codon_count))
 f.close()
 f=open(fname+" "+virus_name+" dicodon count.txt","w")
 f.write(str(dicodon_count))
 f.close()
 f=open(fname+" "+virus_name+" dicodon frequency.txt","w")
 for ele in dicodon_count:
  f.write(ele+" ")
  f.write(str((dicodon_count[ele]/dicodon_number)))
  f.write("\n")
 f.close()


def calculatereverse(filename): #funkcija, kuri skaičiuoja info susijusią su atvirkštiniu komplementu
 fname=filename
 virus_name=""
 for seq_record in SeqIO.parse(fname, "fasta"):
  print(seq_record.id)
  print(repr(seq_record.seq))
  print(len(seq_record))
  virus_name=repr(str(seq_record.name))
 n = 3
 rc_record=seq_record.reverse_complement().seq
 out = [(rc_record[i:i+n]) for i in range(0, len(rc_record), n)]
 length = len(out)
 i = 0
 fragmentcount=0
 f= open(fname+"reverse codonresults.txt","w")
 finalstring=""
# Iterating using while loop 
 while i < length:
     if out[i]=="ATG":
        currentsequence=""
        j=i
        while j<length:
            currentsequence+=out[j]
            if out[j]=="TAG" or out[j]=="TAA" or out[j]=="TGA":
                if len(str(currentsequence))>=100:
                 finalstring += currentsequence
                 f.write(str(currentsequence))
                 f.write("\n")
                i=j
                break
            j+=1
     i += 1
 f.close()
 finalstring=finalstring.lower()
 codons = [finalstring[i:i + 3] for i in range(0, len(finalstring), 3)] #suskaidom visų kodonų sąrašą į array po 3 raides
 dicodons=[str(finalstring[i:i + 6]) for i in range(0, len(finalstring)-3, 3)] #suskaidom visą kodonų sąrašą į array po 6 raides
 codon_count = collections.Counter(aaDict[c] for c in codons) #skaičiuojam kiek kokių yra kodonų
 dicodon_count= collections.Counter(dicodons) #skaičiuojam kiek kokių yra dikodonų
 codon_number=len(codons)
 dicodon_number=len(dicodons)


 f=open(fname+" "+virus_name+"reverse codon frequency.txt","w")

 
 for ele in codon_count:
  f.write(ele+" ")
  f.write(str((codon_count[i]/codon_number,4)))           
  f.write("\n")
 f.close() 
 f=open(fname+" "+virus_name+"reverse codon count.txt","w")  #rašom į faila kodonų skaičius atvirkštiniam komplementui
 f.write(str(codon_count))
 f.close()
 f=open(fname+" "+virus_name+"reverse dicodon count.txt","w") #rašom į failą dikodonu skaičius atvirkštiniam komplementui
 f.write(str(dicodon_count))
 f.close()
 f=open(fname+" "+virus_name+"reverse dicodon frequency.txt","w") #rašom į faila
 for ele in dicodon_count:
  f.write(ele+" ")
  f.write(str((dicodon_count[ele]/dicodon_number)))
  f.write("\n")
 f.close()

 
f=open("distance matrix codons.txt","w+")
for i in range(0, len(orderedNames)):
 f.write(orderedNames[i])
 f.write(" ")
f.write("\n")
f.close()
calculate("bacterial1.fasta")
calculate("bacterial2.fasta")
calculate("bacterial3.fasta")
calculate("bacterial4.fasta")
calculate("mamalian1.fasta")
calculate("mamalian2.fasta")
calculate("mamalian3.fasta")
calculate("mamalian4.fasta")

finalresult=0;
ivalue=0;
jvalue=0;
ikey=""

f=open("Dicodon Distance Matrix.txt","w+")  #Čia sudaroma dikodonu astumų matrica
for i in dicodonList:
    for j in dicodonList:
      for key, value in i.items():
       ivalue=value
       ikey=key
       for keyj, valuej in j.items():
           if(ikey==keyj):
            jvalue=valuej
            finalresult+=(ivalue-jvalue)**2
            break
      f.write(str(round(finalresult,4))) #rašom dikodonu matricą į failą
      f.write(" ")
      finalresult=0
    f.write("\n")
f.close()

calculatereverse("bacterial1.fasta")
calculatereverse("bacterial2.fasta")
calculatereverse("bacterial3.fasta")
calculatereverse("bacterial4.fasta")
calculatereverse("mamalian1.fasta")
calculatereverse("mamalian2.fasta")
calculatereverse("mamalian3.fasta")
calculatereverse("mamalian4.fasta")
