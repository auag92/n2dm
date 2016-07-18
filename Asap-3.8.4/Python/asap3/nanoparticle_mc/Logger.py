import datetime

class Log:
	def __init__(self,name = "simulation_log"):
		d = datetime.date.today()
		t = datetime.datetime.now().time()
		
		self.date = str(d.day)+":"+str(d.month)+":"+str(d.year)
		self.time = str(t.hour)+":"+str(t.minute)
		
		self.filename = name+"_"+str(self.date)+"_"+str(self.time)+".txt"
		f = open(self.filename,"a")
		f.write("Log file containing parameters\n \n")
		f.close()
		
		
	def dumpstring(self,dumpstr):
		f = open(self.filename,"a")
		f.write(dumpstr+"\n")
		f.close()
		
	def dumpparamdict(self,paramdict):
		for p in paramdict:
			self.dumpstring(p+":"+str(paramdict[p]))
	
	
	def load2append(self,filetoload):
		f = open(filetoload,"r")
		str_in = f.read()
		self.dumpstring(str_in)
		
