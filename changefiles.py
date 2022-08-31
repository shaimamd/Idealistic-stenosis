import os
import natsort 
editFiles = []
for file in os.listdir("/scratch/smdsouza/opensource/DNS/matlabpostprocess/"):
    if file.startswith("s11pipe0.f000"):
        editFiles.append(file)
print(editFiles)
        #os.rename(old_file, new_file)
      

sortfiles=natsort.natsorted(editFiles,reverse=False)
print(sortfiles)
k=1
for i in sortfiles:
    converted_num = str(k)
    os.rename(i,'s11pipe0.f'+converted_num.rjust(5, '0'))
    k=k+1
 
