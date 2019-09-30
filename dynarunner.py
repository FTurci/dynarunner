import subprocess
import os
import numpy as np
import os
import xml.etree.ElementTree as ET
import string
 
def discrete_std(a,b):
  return np.sqrt(((b-a+1)**2-1.)/12.)
#A helpful function to load compressed or uncompressed XML files
def loadXMLFile(filename):
    #Check if the file is compressed or not, and 
    if (os.path.splitext(filename)[1][1:].strip() == "bz2"):
        import bz2
        f = bz2.BZ2File(filename)
        doc = ET.parse(f)
        f.close()
        return doc
    else:
        return ET.parse(filename)

def read_pressure(filename):
  doc = loadXMLFile(filename)
  p=doc.find(".//Misc")
  return float(p.find("Pressure").get("Avg"))

def read_temperature(filename):
  doc = loadXMLFile(filename)
  p=doc.find(".//Misc")
  return float(p.find("Temperature").get("Mean") )

def conf_to_xyz(filename):
  XMLDoc=loadXMLFile(filename)

  #We can create a list of all particle tags using an xpath expression
  #(xpath expressions always return lists)
  PtTags = XMLDoc.findall(".//Pt")
  #print the x, y, and z positions of each particle
  with open(filename+".xyz", 'w') as fw:
    fw.write(f"{len(PtTags)}\nAtoms\n")
    for PtElement in PtTags:
        PosTag = PtElement.find("P")
        x,y,z=PosTag.get("x"),PosTag.get("y"),PosTag.get("z")
        fw.write(f"A {x} {y} {z}\n")

def conf_to_atom(filename,timestep):
  XMLDoc=loadXMLFile(filename)
  types=list(string.ascii_uppercase)

  SizeTag = XMLDoc.find(".//SimulationSize")
  sizes = float(SizeTag.get("x")),float(SizeTag.get("y")),float(SizeTag.get("z"))

  box =np.array([[-l*0.5, l*0.5] for l in sizes])
  xlo,xhi = box[0,0],box[0,1]
  ylo,yhi = box[1,0],box[1,1]
  zlo,zhi = box[2,0],box[2,1]

  PtTags = XMLDoc.findall(".//Pt")
  N=len(PtTags)
  #print the x, y, and z positions of each particle
  with open(filename+".atom", 'w') as fw:
    fw.write(f"""ITEM: TIMESTEP
{timestep}
ITEM: NUMBER OF ATOMS
{N}
ITEM: BOX BOUNDS pp pp pp
{xlo} {xhi}
{ylo} {yhi}
{zlo} {zhi}\n""")

    # check if there are diameters
    p0=PtTags[0]

    if p0.get("D") is not None:
      print(":: Found a diameter. Saving radius information.")
      radius = True
      fw.write("ITEM: ATOMS type x y z radius\n")
    else:
      print(":: No diameters. Assuming single species.",flush=True)
      radius = False
      fw.write("ITEM: ATOMS type x y z\n")
  
    for PtElement in PtTags:
      PosTag = PtElement.find("P")
      x,y,z=PosTag.get("x"),PosTag.get("y"),PosTag.get("z")
      
      if radius:
        diam = float(PtElement.get("D"))
        fw.write(f"A {x} {y} {z} {diam/2}\n")
      else:
        fw.write(f"A {x} {y} {z}\n")


def create_polydisperse_from_conf(confin, confout,sizes, probabilities):
  #Load the XML file
  
  XMLFile = ET.parse(f"{confin}")

  #Add the diameter and mass tags
  for ParticleTag in XMLFile.findall("./ParticleData/Pt"):
      #Calculate a random diameter
      diameter = np.random.choice(sizes,p=probabilities)
      #Add attributes for the diameter
      ParticleTag.attrib["D"] = str(diameter)
      #Add a mass attribute which scales with the particle volume (for D=1 M=1)
      ParticleTag.attrib["M"] = str(diameter * diameter * diameter)

  #Tell DynamO about these properties
  massProperty = ET.SubElement(XMLFile.findall("./Properties")[0], "Property")
  massProperty.attrib["Type"] = "PerParticle"
  massProperty.attrib["Units"] = "Mass"
  massProperty.attrib["Name"] = "M"
  diamProperty = ET.SubElement(XMLFile.findall("./Properties")[0], "Property")
  diamProperty.attrib["Type"] = "PerParticle"
  diamProperty.attrib["Units"] = "Length"
  diamProperty.attrib["Name"] = "D"

  #Change the Interaction and Species to use the properties
  XMLFile.findall(".//Interaction")[0].attrib["Diameter"] = "D"
  XMLFile.findall(".//Species")[0].attrib["Mass"] = "M"

  #Write out XML file
  XMLFile.write(f"{confout}",xml_declaration=True)

class DynamoMD:
  def __init__(self,root,output="./output"):
    self.root = root
    self.output= output
    os.system(f"mkdir {self.output} ")
    self.run_counter = 0
    
  def configure(self, density,C,polydispersity=None, N=None):
    self.conf_density= density
    self.C = C
    self.start = self.output+"/config.start.xml"
    os.system(f"{self.root}/dynamod -m 0 -d {density} -C {C} --i1 1 -o {self.start} -r 1 -T 1.")
    if N is not None:
      self.take_n(self.start, N)

    if polydispersity is not None:
      self.polydispersity = polydispersity
      # assume 5 component system with assigned polydispersity
      # uni
      n = 5
      p = 1.0/n
      mean = 1.0
      radii = mean+np.arange(-(n-1)/2, (n+1)/2)/discrete_std(0,n-1)*polydispersity/mean
      probs = np.ones(n)*p
      create_polydisperse_from_conf(self.start,self.start,radii,probs)

  def take_n(self,conf,n):
    doc = loadXMLFile(conf)
    docroot = doc.getroot()
    print("Docroot")
    particledata=docroot.find('.//ParticleData')
    for Pt in particledata.findall('.//Pt'):      
      Id = int(Pt.get('ID'))
      if Id >=n:
        particledata.remove(Pt)
    doc.write(conf)

  def compress(self, target_packing, conf=None):
    print(":: Compression", flush=True)
    if conf==None:
      conf= self.start

    self.compressed = f"{self.output}/compression.target.{target_packing}.conf.xml.bz2"
    os.system(f"{self.root}/dynarun --engine=3 --target-pack-frac {target_packing} --config-file {conf}  -o {self.compressed} --out-data-file {self.output}/compression.target.{target_packing}.log.xml")

  def equilibrate(self,finaltime,conf=None):
    print(":: Equilibrate", flush=True)
    if conf == None:
      conf = self.compressed

    self.equilibrated = f"{self.output}/equilibration.conf.xml.bz2"

    os.system(f"{self.root}/dynarun -E -f {finaltime} --engine=1 --config-file {conf} -o {self.equilibrated} --out-data-file {self.output}/equilibration.log.xml")
    conf_to_atom(self.equilibrated, finaltime)

  def run(self,conf,finaltime,snapshot=10):
    self.run_counter+=1
    os.system(f"{self.root}/dynarun -f {finaltime} --engine=1 --config-file {conf} --snapshot {snapshot} -o {self.output}/config.end.xml.bz2 --out-data-file {self.output}/output.end.xml.bz2")
    os.system(f"mkdir {self.output}/run{self.run_counter}")
    os.system(f"mv Snapshot* {self.output}/run{self.run_counter}")
