#!/usr/bin/env python
import commands
x=[]
y=[]

def write_file(xx,yy):
      fileName="data.big.%d.%d.json"%(int(xx),int(yy));
      x=open(fileName,"w")
      x.write('{\n')
      x.write('"d1" : 0.004,\n')
      x.write('"d2" : 25,\n')
      x.write('"d3" : 25,\n')
      x.write('"d4" : 25,\n')
      x.write('"d5" : 25,\n')
      x.write('"filename" : "/scr1/bob/data.big.json.%d.%d.dat",\n'%(int(xx),int(yy)))
      x.write('"label1" : "Undefined",\n')
      x.write('"label2" : "Undefined",\n')
      x.write('"label3" : "Undefined",\n')
      x.write('"label4" : "Undefined",\n')
      x.write('"label5" : "Undefined",\n')
      x.write('"n1" : 1500,\n')
      x.write('"n2" : 240,\n')
      x.write('"n3" : 240,\n')
      x.write('"n4" : 1,\n')
      x.write('"n5" : 1,\n')
      x.write('"o1" : 0,\n')
      x.write('"o2" : 0,\n')
      x.write('"o3" : 0,\n')
      x.write('"o4" : %d,\n'%xx)
      x.write('"o5" : %d\n}\n'%yy)
      x.close()
      par="model.%d.%d.P"%(xx,yy);
      y=open(par,"w");
      y.write("{\n")
      y.write('"data": "data.big.%d.%d.json",\n'%(xx,yy))
      y.write('"image": "refl.big.json\",\n')
      y.write('"velocity": "vel.big.json\",\n')
      y.write('"wavelet": "wavelet.json"\n}\n')
      y.close()
  
def copy_file(x,y):
  print "hould copy"
  commands.getoutput("cp /scr1/bob/data.big.json.dat /scr1/bob/data.big.%d.%d.json.dat"%(x,y));
  
  
def run_job(x,y):
  xx="Model3D json=model.%d.%d.P"%(x,y)
  print xx;
  commands.getoutput("Model3D json=model.%d.%d.P"%(x,y))
	
  





for i2 in range(38):
  for i1 in range(38):
     x.append(i2*200.);
     y.append(i1*200.);
     



for i in range(len(x)):
   print "writing file",x[i],y[i]
   write_file(x[i],y[i]);
   print "copying file",x[i],y[i]
   copy_file(int(x[i]),int(y[i]))
   run_job(int(x[i]),int(y[i]));
   