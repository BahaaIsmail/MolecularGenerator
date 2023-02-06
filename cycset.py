
### approved for cyclic atoms
from os import walk ,getcwd , system , remove
from random import uniform

d = {1 : {'a':{1:'C' , 2:'H'} , 'b':{1:[1,2,1]} , 'c':[1,2]} , 
     2 : {'a':{1:'C' , 2:'H'} , 'b':{1:[1,2,1]} , 'c':[2,1]} , 
     3 : {'a':{1:'C' , 2:'O'} , 'b':{1:[1,2,2]} , 'c':[1,1]} ,
     4 : {'a':{1:'C' , 2:'C' , 3:'H' , 4:'H'} , 'b':{1:[1,2,2] , 2:[2,3,1] , 3:[2,4,1]} , 'c':[1,1]} ,
     5 : {'a':{1:'N+', 2:'H'} , 'b':{1:[1,2,1]} , 'c':[1,2]} , 
     6 : {'a':{1:'N+', 2:'H'} , 'b':{1:[1,2,1]} , 'c':[2,1]} , 
     7 : {'a':{1:'N+', 2:'O'} , 'b':{1:[1,2,2]} , 'c':[1,1]} ,
     8 : {'a':{1:'N'}, 'b':{} , 'c':[1,2]} , 
     9 : {'a':{1:'N'}, 'b':{} , 'c':[2,1]} , 
     10: {'a':{1:'C' , 2:'H' , 3:'H'} , 'b':{1:[1,2,1] , 2:[1,3,1]} , 'c':[1,1]} ,
     12: {'a':{1:'N' , 2:'H'} , 'b':{1:[1,2,1]} , 'c':[1,1]} , 
     13: {'a':{1:'O'}, 'b':{} , 'c':[1,1]} , 
     14: {'a':{1:'S'}, 'b':{} , 'c':[1,1]}}

dl = {1:{1 : {'a':{1:'C'} , 'b':{} , 'c':[1,2]} , 
         2 : {'a':{1:'C'} , 'b':{} , 'c':[2,1]} , 
         5 : {'a':{1:'N+'} , 'b':{} , 'c':[1,2]} , 
         6 : {'a':{1:'N+'} , 'b':{} , 'c':[2,1]} , 
         12: {'a':{1:'N'} , 'b':{} , 'c':[1,1]} } , 

      2:{3 : {'a':{1:'C'} , 'b':{} , 'c':[1,1]} ,
         7 : {'a':{1:'N+'} , 'b':{} , 'c':[1,1]}}}
         

def new(L,pool) :
    s = {1:1, 2:1, 3:3, 4:4, 5:5, 6:5, 7:7, 8:8, 9:8, 10:10, 12:12, 13:13, 14:14}
    M = [s[i] for i in L]
    for i in range(1,len(M)+1) :
        if M[i:]+ M[:i] in pool :  
            return 0   
    M = [s[L[-i]] for i in range(1,len(L)+1)]
    for i in range(1,len(M)+1) :
        if M[i:]+ M[:i] in pool :  
            return 0          
    pool += [M]
    return pool

    
def link(i,j):
    return d[i]['c'][1] == d[j]['c'][0]
        

def aromatic(L): 
    aroma = 0
    for i in L : 
        if i == 10 : 
            return 0
        elif i in [1,2,5,6,8,9]: 
            aroma += 1
        elif i in [12,13,14] : 
            aroma += 2
    if aroma == 6 : 
        return 1
            

def valid(g) : 
    cord = {'C':4 , 'N+':4 , 'N':3 , 'O':2 , 'S':2 , 'O-':1 , 'S-':1 , 'H':1 , 'F':1 , 'BR':1 ,'CL':1 , 'I':1}
    for i in g['a'] : 
        nbs = 0
        for j in g['b'] : 
            if i in g['b'][j][0:2] : 
                nbs += g['b'][j][2]
        if nbs != cord[g['a'][i]]:
            return 0
    return 1
            
          

def export(L,x,c): 
    print
    print c , '  _______________________________________________________________' 
    if 1 :
        g = {'a':{} , 'b':{}}
        n = len(L)
        a = 1
        b = 1
        al = []
        for h in range(n) :  
            f = L[h]
            ai = a
            if x and h in [x[0],x[1]]: 
                al += [ai]       
            for i in f['a'] :             
                g['a'][a] = f['a'][i]
                a += 1
            for i in f['b'] :
                g['b'][b] = []
                for j in f['b'][i][0:2] :
                    g['b'][b] += [j+ai-1]            
                g['b'][b] += [f['b'][i][2]]             
                b += 1
            if h < n-1 : 
                g['b'][b] = [ai,a,f['c'][1]]
                b += 1
        g['b'][b] = [1,ai,L[0]['c'][0]] 
        if x : 
            g['b'][b+1] = [al[0] ,al[1],x[2] ]
            b += 1
            
        #if valid(g) : 
        mol2 = 'cyc_'+str(c)+'.mol2'
        tf = open(mol2,'w')
        tf.write('@<TRIPOS>MOLECULE\n*****\n')
        tf.write(str(a-1) + '\t' + str(b) + '  0  0  0\n')
        tf.write('SMALL\n')
        tf.write('molecule\n\n')
        tf.write('@<TRIPOS>ATOM\n')
        x = [0,0,0]
        for i in g['a'] : 
            j = str(i)
            at = g['a'][i].strip('+-') 
            x = [x[0]+uniform(-2,2),x[1]+uniform(-2,2),x[2]+uniform(-2,2)]
            tf.write(j + '\t' + at+j + '\t' + str(x[0])+  '\t' + str(x[1])+ '\t' + str(x[2]) + '\t' + at + '\n') 
        tf.write('@<TRIPOS>BOND\n')
        for i in g['b'] : 
            b = g['b'][i]
            tf.write( str(i) + '\t' +str(b[0]) + '\t' + str(b[1]) + '\t' + str(b[2]) + '\n' )     
        tf.close() 
    c = c + 1
    return c

c = 1   

### 3-mem cycles
pool = []
for i in [3,10] : 
    for j in [3,10] :  
        for k in [3,10] : 
            if new ([i,j,k],pool):    
                c = export([d[i],d[j],d[k]],0,c)


### 4-mem cycles
pool = []
for i in [3,10,12,13] : 
    for j in [3,10,12,13] :  
        for k in [3,10,12,13] : 
            for l in [3,10,12,13]: 
                if new ([i,j,k,l],pool):    
                    c = export([d[i],d[j],d[k],d[l]],0,c)

### one 5-mem cycle 
pool = []    
for i in d :
    if link(3,i) : 
        for j in d : 
            if link(i,j) : 
                for k in d : 
                    if link(j,k): 
                        for l in d : 
                            if link(k,l): 
                                if link(l,3) :
                                    if new([i,j,k,l,3],pool) : 
                                        c = export([d[i],d[j],d[k],d[l],d[3]],0,c) 

### one 6-mem aromatic ring
pool = []    
for i in d :
    if link(2,i) : 
        for j in d : 
            if link(i,j) : 
                for k in d : 
                    if link(j,k): 
                        for l in d : 
                            if link(k,l): 
                                if link(l,1) :
                                    if aromatic([i,j,k,l,1,2]):
                                        if new ([i,j,k,l,1,2],pool):                                              
                                            c = export([d[i],d[j],d[k],d[l],d[1],d[2]],0,c)

### one 7-mem aromatic ring
pool = []    
for i in d :
    if link(3,i) : 
        for j in d : 
            if link(i,j) : 
                for k in d : 
                    if link(j,k): 
                        for l in d : 
                            if link(k,l): 
                                if link(l,1) :
                                    if aromatic([i,j,k,l,1,2]):
                                        if new ([i,j,k,l,1,2],pool):                                              
                                            c = export([d[i],d[j],d[k],d[l],d[1],d[2]],0,c)
                                        if new ([i,j,k,l,1,2,3],pool):    
                                                c = export([d[i],d[j],d[k],d[l],d[1],d[2],d[3]],0,c)


### fused cycles
def aromatic(L): 
    aroma = 0
    for j in range(len(L)): 
        i = L[j] 
        if i == 10: 
            return 0
        elif j == 2 and i in [2,6,9]:
            aroma += 1            
        elif j == 3 and i in [1,5,8]:            
            aroma += 1            
        elif j not in [2,3] and i in [1,2,5,6,8,9]: 
            aroma += 1
        elif i in [12,13,14] : 
            aroma += 2
    return aroma




### fused two cycles 5/6 , 5/7 and 6/7
pool = []    
for i in d :
    for j in d : 
        if link(i,j) : 
            for k in dl[1] : 
                if link(j,k): 
                    for l in d : 
                        if link(k,l): 
                            for m in d :
                                if link(l,m): 
                                    for n in d : 
                                        if link(m,n) : 
                                            for o in dl[1] : 
                                                if link(n,o): 
                                                    for p in d : 
                                                        if link(o,p): 
                                                            for q in d : 
                                                                if link(p,q) and link(q,i):
                                                                    if aromatic([i,j,k,o,p,q]) == 6: 
                                                                        #if new([i,j,k,l,m,n,o,p,q],pool) :                                                                        
                                                                        c = export([d[i],d[j],dl[1][k],d[l],d[m],d[n],dl[1][o],d[p],d[q]],[2,6,1],c)
                                                                        #if link(q,3) and link(3,i): 
                                                                         #   if new([i,j,k,l,m,n,o,p,q,3],pool) :
                                                                          #      c = export([d[i],d[j],dl[1][k],d[l],d[m],d[n],dl[1][o],d[p],d[q],d[3]],[2,6,1],c)

pool = []    
for i in d :
    for j in d : 
        if link(i,j) : 
            for k in dl[2] : 
                if link(j,k): 
                    for l in d : 
                        if link(k,l): 
                            for m in d :
                                if link(l,m): 
                                    for n in d : 
                                        if link(m,n) : 
                                            for o in dl[2] : 
                                                if link(n,o): 
                                                    for p in d : 
                                                        if link(o,p): 
                                                            for q in d : 
                                                                if link(p,q) and link(q,i): 
                                                                    if aromatic([i,j,k,o,p,q]) == 4: 
                                                                        #if new([i,j,k,l,m,n,o,p,q],pool) :                                                                        
                                                                        c = export([d[i],d[j],dl[2][k],d[l],d[m],d[n],dl[2][o],d[p],d[q]],[2,6,2],c)
                                                                        #if link(q,3) and link(3,i): 
                                                                         #   if new([i,j,k,l,m,n,o,p,q,3],pool) :
                                                                          #       c = export([d[i],d[j],dl[2][k],d[l],d[m],d[n],dl[2][o],d[p],d[q],d[3]],[2,6,2],c)

###############################################################################
###############################################################################
###############################################################################
### execute the cgenff for many files 

from os import walk ,getcwd , system

count = 1
### getting the mol2 file
mol2 = 0
path = getcwd()
for root , dirs, files in walk(path):
    for f in files :
        if '.mol2' in f :
            print count , '===================  ' , f , '  ===================\n'
            srt = f.replace('mol2','srt')
            system('./cgenff -a ' + f + ' > ' + srt)
            print '\n\n\n'
            count += 1



print '***************************************'
print 'ENDED SUCCESSFULLY'
                        
                
           




         
'''







                                                                             


                                                                          
# fused 6/7
pool = []    
for i in d :
    if link(3,i):
        for j in d : 
            if link(i,j) : 
                for k in d : 
                    if link(j,k): 
                        for l in dl[1] : 
                            if link(k,l): 
                                for m in d :
                                    if link(l,m): 
                                        for n in d : 
                                            if link(m,n) : 
                                                for o in d : 
                                                    if link(n,o): 
                                                        for p in d : 
                                                            if link(o,p): 
                                                                for q in dl[1] : 
                                                                    if link(p,q):
                                                                        for r in d : 
                                                                            if link(q,r) and link(r,3):
                                                                                if new([i,j,k,l,m,n,o,p,q,r,3],pool) :         
                                                                                    if aromatic([j,k,l,q,r,i]) == 6 and aromatic ([o,p,q,l,m,n]) == 6 :                                                                                                                                        
                                                                                        c = export([d[i],d[j],d[k],dl[1][l],d[m],d[n],d[o],d[p],dl[1][q],d[r],d[3]],[3,8,1],c)
 
pool = []    
for i in d :
    if link(3,i):
        for j in d : 
            if link(i,j) : 
                for k in d : 
                    if link(j,k): 
                        for l in dl[2] : 
                            if link(k,l): 
                                for m in d :
                                    if link(l,m): 
                                        for n in d : 
                                            if link(m,n) : 
                                                for o in d : 
                                                    if link(n,o): 
                                                        for p in d : 
                                                            if link(o,p): 
                                                                for q in dl[2] : 
                                                                    if link(p,q):
                                                                        for r in d : 
                                                                            if link(q,r) and link(r,3):
                                                                                if new([i,j,k,l,m,n,o,p,q,r,3],pool) :         
                                                                                    if aromatic([j,k,l,q,r,i]) == 6 and aromatic ([o,p,q,l,m,n]) == 6 :                                                                                                                                        
                                                                                        c = export([d[i],d[j],d[k],dl[2][l],d[m],d[n],d[o],d[p],dl[2][q],d[r],d[3]],[3,8,2],c)
  

'''