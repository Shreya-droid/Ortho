#Pipeline to collect orthologs
#Date:21-12-2020
#Version:1.0
#Author:Shreya Sharma

import glob
import os
import re
import ast

def Convert(string):
    l = list(string.split(" "))
    return l

def Replace(line):
    line_strpd = line.replace("[","").replace("]","").replace(",","").replace("'","").rstrip('\n')
    return line_strpd

def stripline(line):
    xd = line.split(",")
    result = []
    for item in xd:
        result.append(re.sub(r"[\'\[\]]","",item).strip())
    return result

def checksort(myfile):
    mydict = dict()
    with open(myfile) as thefile:
        x = thefile.readline()
        while x != "":
            x = stripline(x)
        
            if bool(mydict.get(x[0])):
                mydict[x[0]] = mydict[x[0]].union(set(x))
            else:
                mydict[x[0]] = set(x)
            x = thefile.readline()
    return mydict

directory = os.chdir(input("enter the dir path: "))
identity = float(input("enter your desired cutoff identity: "))
coverage = float(input("enter your desired qlen(coverage): "))
eValue = float(input("enter your desired evalue: "))

###collecting only qualified ones
result=""
for file in os.listdir(directory):
    if file.endswith(".blast"):
        with open(file,'r') as f:
            lines = f.readlines()
            for i in range(len(lines)):
                q = lines[i].split('\t')
                if ( float(q[2]) >= identity) and (float(q[3]) >= coverage) and (float(q[10]) <= eValue): 
                    result+=lines[i]
                i+=1       
            my_outfile = "".join(file.split('.')[:-1]) + '_qualified.txt'  
            output_file = open(my_outfile, "w")
            output_file.write(result)
            output_file.close()

              
with open("rbh_output.txt", "w") as fa: 
    with open(my_outfile, "r") as filehandle:
        lines = filehandle.readlines()
        for i in range(len(lines)):
            l = lines[i].split('\t')
            splitted_first_ele_l = l[0].split("_")
            splitted_second_ele_l =l[1].split("_")
            if ((l[0]!=l[1]) and ((splitted_first_ele_l[1])!=(splitted_second_ele_l[1]))):
                collect = l[:2]
                fa.write("%s\n" %collect)
                i+=1
            else:
                continue
fa.close()      


f = open("rbh_output.txt", "r")
Lines = f.readlines()
i = 0
j = 0
outputlines = ''
skipNumber = ''
# extract the line not to be printed in file in skipnumbers
for line in Lines:
    if i==0:
        outputlines = (outputlines + line).replace("[","").replace("]","").replace(",","").replace("'","")
        i = i+1
    else :
        line1 = line.replace("[","").replace("]","").replace(",","").replace("'","").rstrip('\n')
        for l in outputlines.splitlines():
            if l == line1 or l == ' '.join(reversed(line1.split(' '))) :
                skipNumber = skipNumber + str(i) +'\n'
        outputlines = (outputlines + line).replace("[","").replace("]","").replace(",","").replace("'","")
        i = i+1
i = 0
#write the lines into file except skip number lines from file
file = open("rbh.txt", "w")
for line in Lines:
    if j==0:
        outputlines = line
        j = j+1
    else :
        for k in skipNumber.splitlines():
            if k==str(j):
                i=1
        if i==0:
            outputlines = outputlines + line
        i = 0
        j = j+1
file.write(outputlines)
file.close()

f = open('rbh.txt','r')
lines = f.readlines()
i = 0
fileLines = ''
os.mkdir('temp')
for line in lines:
    line_li = Convert(Replace(line))
    i+=1
    for element in line_li:
        if len(element) != 0 :
            for line in lines:
                if element in line:
                    if line not in fileLines :
                        fileLines = fileLines + line
    fi = open(os.path.join('temp/', line_li[0] +'-'+ line_li[1]+".txt"), "w")
    fi.write(fileLines)
    fi.close()
    fileLines = ''
print("Processing....")

final_file = open("final.txt","a+")
arr = os.listdir('temp')
for fileName in arr:
    fileN = open(os.path.join('temp/', fileName),"r")
    lines=fileN.readlines()
    element_A = []
    first_line = fileName.split('-')
    first_element = first_line[0]
    second_element = (first_line[1]).rstrip(".txt")
    element_A.append(first_element)
    element_A.append(second_element)
    for line in lines:
        line = Convert(Replace(line))
        if line!=first_line:
            if(first_element in line):
                if(first_element == line[0]):
                    ele_to_be_append = line[1]
                    element_A.append(ele_to_be_append)
                else:
                    ele_to_be_append = line[0]
                    element_A.append(ele_to_be_append)
            if(second_element in line):
                if(second_element == line[0]):
                    ele_to_be_append = line[1]
                    element_A.append(ele_to_be_append)
                else:
                    ele_to_be_append = line[0]
                    element_A.append(ele_to_be_append)
    final_file.write("%s\n" %element_A)
    element_A.clear()
final_file.close()


lines_seen=set()
out_f = open("yay_final.txt","w")
in_f = open("final.txt","r")
for line in in_f:
    if line not in lines_seen:
        out_f.write(line)
        lines_seen.add(line)
out_f.close()


file_w = open("yay2.txt","w")
file_l = open("Final.txt","r")
file_lines=file_l.readlines()
for line_x in file_lines:
    line_x = Convert(Replace(line_x))
    line_x = list(dict.fromkeys(line_x))
    file_w.write("%s\n" %line_x)
file_w.close()     

res = checksort("yay2.txt")
filee = open("res.txt", "w")
data = str(res)
filee.write(data)
filee.close()

ast.literal_eval(data)
all_keys = []
for key in res:
    all_keys.append(key)
final_file = open("orthologous_gene.txt","w")
temp_list = []
rem_ele = []
rem_genome = []
for element_as_key in all_keys:
    values_of_key = (res[element_as_key])
    for value_in in values_of_key:
        if value_in not in all_keys:
            rem_ele.append(value_in)
     
        else:
            if (set(res[element_as_key]) == set(res[value_in])):
                genome = value_in.split("_")[1]
                temp_list.append(genome)
                for s in temp_list:
                    coun = temp_list.count(s)
                    if(coun == 1):
                        continue
                    else:
                        if value_in not in rem_genome:
                            rem_genome.append(value_in)
        remove_from = set(rem_ele).union(set(rem_genome))
        dataset = set(values_of_key)-set(remove_from)
    final_file.write("\n")    
    final_file.write(str(dataset))
    temp_list.clear()
    rem_ele.clear()
    rem_genome.clear()
final_file.close()

file_end = open("orthologs_from_rbh.txt","w")
aa = open("orthologous_gene.txt","r")
aaa = aa.readlines()
file_y = []

for line in aaa:
    if line not in file_y:
        file_y.append(line)
        line_to_app = str(line)
    file_end.write(line_to_app)
file_end.close()


###cleanUp

#for CleanUp in glob.glob(directory):
    #print(CleanUp)
    #if not CleanUp.endswith("orthologs_from_rbh.txt"):
        #os.remove(CleanUp)



