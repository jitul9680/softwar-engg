import matplotlib.pyplot as plt
from acid_analysis import * 




names_of_aminoacids=["Alanine","Cysteine","Aspartic_acid  ","Glutamic_acid  ","Phenylalanine  ","Glycine","Histidine","Isoleucine","Lysine","Leucine","Methionine  ","Asparagine","Proline","Glutamine","Arginine","Serine","Threonine","Valine","Tryptophan","Tyrosine"]
numeric_values=[val[0]*100/sum(val),val[1]*100/sum(val),val[2]*100/sum(val),val[3]*100/sum(val),val[4]*100/sum(val),val[5]*100/sum(val),val[6]*100/sum(val),val[7]*100/sum(val),val[8]*100/sum(val),val[9]*100/sum(val),val[10]*100/sum(val),val[11]*100/sum(val),val[12]*100/sum(val),val[13]*100/sum(val),val[14]*100/sum(val),val[15]*100/sum(val),val[16]*100/sum(val),val[17]*100/sum(val),val[18]*100/sum(val),val[19]*100/sum(val)]
plt.bar(names_of_aminoacids,numeric_values)
plt.xlabel("****************Amino Acids********************")
plt.ylabel("****************Quantities***********************")
#plt.setp(plt.xticklabels(), rotation=30, horizontalalignment='right')
plt.xticks(rotation=45)
plt.show()