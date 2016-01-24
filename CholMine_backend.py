#! /usr/bin/env python
# CLR & CHD prediction server
#v 1.2, Python 2.7
#
# Nan Liu
#10/01/13
#The server has two options: CLR and CHD. The two options cannot be used at the same time. The current version only support one binding site input, given the protein file, the ligand residue name, the ligand residue number and the ligand chain ID.
# support multiple models and the first model will be chosen.
# the bug about counting ATOM and HETATM lines has been fixed. 
# fix the CLR and CHD connect record problem
# let usr define cleft residues
# fill in water in the cleft defined and generate fake ligand
# add judgement of cleft box volume, which ranges [0, 10000]
# add judgement of cleft residue numbers cutoff below 25
# add judgement of fake ligand, if file empty, return
# add judgement of miniPDB as subset of original regular PDB
# add sitemap points count judgement
# fix bug to valid PDB if pdb contains no HETATM
# add tail explanation

from __future__ import division
import argparse
import os
import sys
import numpy as np
from ASCbasePy.utils import pdb
import csv
from datetime import datetime
 
def usage_args():   
   parser = argparse.ArgumentParser(prog="CholMine",usage='%(prog)s [options]',description='version1.00 server',formatter_class=argparse.RawDescriptionHelpFormatter,epilog='''
   Typically,
   CholMine -CLR -opse -projout /usr/SimSite3D/CholMine -pfile xxxx  -lnum xxx -lchain x
   CholMine -CHD -opse -projout /usr/SimSite3D/CholMine -pfile xxxx  -lnum xxx -lchain x 
   CholMine -CLR -projout /usr/SimSite3D/CholMine -pfile xxxx -keyresidue xxx''')
   parser.add_argument("-CLR", "--CLR", help="CLR predictor turned on",action='store_true')
   parser.add_argument("-CHD", "--CHD", help ="CHD predictor turned on", action='store_true')
   parser.add_argument("-opse", "--opse", help="output PyMOL pse file",action='store_true')
   parser.add_argument("-v", help = "debug, output all the resulting files", action='store_true')
   parser.add_argument("-projout", "--projout", help = "usr defined pathway for the output directory, such as /usr/SimSite3D/CholMine")
   parser.add_argument("-pfile", "--pfile",  help='protein file name such as 2RH1')
#   parser.add_argument("-lname", "--lname", help='ligand residue name such as CLR')
   parser.add_argument("-lnum", "--lnum", help='ligand residue number such as 412')
   parser.add_argument("-lchain", "--lchain",  help='ligand chainID such as A. If the chain ID is not supplied, the default chain ID will be _.')
   parser.add_argument("-keyresidue", "--keyresidue",  help='keyresidue file name such as keyresidue')
   args = parser.parse_args()
   return args
def usage_help():
   parser = argparse.ArgumentParser(prog="CholMine",usage='%(prog)s [options]',description='version1.00 server',formatter_class=argparse.RawDescriptionHelpFormatter,epilog='''
   Typically
   CholMine -CLR -projout /usr/SimSite3D/CholMine -pfile xxxx -lnum xxx -lchain x
   CholMine -CHD -projout /usr/SimSite3D/CholMine -pfile xxxx -lnum xxx -lchain x 
   CholMine -CLR -projout /usr/SimSite3D/CholMine -pfile xxxx -keyresidue xxx''')
   parser.add_argument("-CLR", "--CLR", help="CLR predictor turned on",action='store_true')
   parser.add_argument("-CHD", "--CHD", help ="CHD predictor turned on", action='store_true')
   parser.add_argument("-opse", "--opse", help="output PyMOL pse file",action='store_true')
   parser.add_argument("-v", "--v", help = "debug, output all the resulting files", action='store_true')
   parser.add_argument("-projout", "--projout", help = "usr defined pathway for the output directory, such as /usr/SimSite3D/CholMine")
   parser.add_argument("-pfile", "--pfile",  help='protein file name such as 2RH1')
 #  parser.add_argument("-lname", "--lname", help='ligand residue name such as CLR')
   parser.add_argument("-lnum", "--lnum", help='ligand residue number such as 412')
   parser.add_argument("-lchain", "--lchain",  help='ligand chainID such as A. If the chain ID is not supplied, the default chain ID will be _.')
   parser.add_argument("-keyresidue", "--keyresidue",  help='keyresidue file name such as keyresidue')
   help = parser.print_help()
   return help

def directory(projpath, pdbfile):
   os.system("rm -rf " + projpath + "/" + pdbfile)
   os.system("mkdir "  +  projpath + "/" + pdbfile)
   os.system("mkdir " +  projpath + "/" + pdbfile + "/initialdata")
   dest = projpath + "/" + pdbfile + "/"
   
   os.system("mkdir " + dest + "dataset")
   os.system("mkdir " + dest + "dataset/ligands")
   os.system("mkdir " + dest + "dataset/dbase")
   return dest

def prepare(dest,wholename,pdbfilefull,linput,logfile,v ):
   #prepare initial protein files
   input = open(dest + pdbfilefull, "r")
   output = open(dest + "initialdata/" + wholename + "_p.pdb", "w")
   for line in input:
      line = line.strip()
      if line.startswith("ATOM"):
         print >>output, line
      if line.startswith("CONECT"):
         print >> output, line
      if line.startswith("MODEL"):
         line = line.split()
         if float(line[1])==2:
            break
   input.close()
   output.close()
   #prepare initial ligand file
   input = open(dest + pdbfilefull, "r")
   output = open(dest + "initialdata/" + wholename + "_l.pdb", "w")
   Lig=[]
   for line in input:
      line = line.strip()
      if line.startswith("HETATM"):
        # lname = line[17:20]
        # lname = lname.strip()
         lindex = line[22:26]
         lindex = lindex.strip()
         cid = line[21]
         cid = cid.strip()
         if not cid:
            cid = ""
        # linfo = lname + lindex + cid
         linfo = cid + "_" + lindex
         if cmp(linput[4:].lower(), linfo.lower())==0:
            print >>output, line
            Lig.append(line)
      if line.startswith("CONECT"):
         print >> output, line
      if line.startswith("MODEL"):
         line = line.split()
         if float(line[1])==2:
            break
   input.close()
   output.close()
   if not Lig:
      print >>logfile, "ERROR: This ligand site cannot be found in the PDB file. Please recheck the chain ID and residue number or the format of the protein file and ensure that the ligand is labeled as HETATM rather than ATOM."
      print "ERROR: This ligand site cannot be found in the PDB file. Please recheck the chain ID and residue number or the format of the protein file and ensure that the ligand is labeled as HETATM rather than ATOM."
      if  not v:
         os.system("rm -rf " + dest + "initialdata/")
         os.system("rm -rf " + dest + "dataset")
         os.system("rm -rf " + dest + "*.pdb")
         os.system("rm -rf " + dest + pdbfilefull)
         os.system("rm -rf " + dest + "center.out")
      sys.exit(13)

#check whether the mini PDB is subset of the original regular PDB file
def checksubset(dest, keyresidue, pdbfilefull, v):
   f = open(dest + keyresidue, "r")
   count = 0
   list = []
   for l in f:
      if l.startswith("ATOM"):
         count = count +1
         l = l.strip()
         info = l[17:65]
         ref = open(dest + pdbfilefull, "r")
         for line in ref:
            line = line.strip()
            if (info in line):
               list.append("yes")
               break
         ref.close()
   f.close()
   if len(list) < count:
      print "ERROR: The mini PDB need to be a subset of the regular PDB file."
      if not v:
         os.system("rm -rf " + dest + "initialdata/")
         os.system("rm -rf " + dest + "dataset")
         os.system("rm -rf " + dest + "*.pdb")
         os.system("rm -rf " + dest + pdbfilefull)
         os.system("rm -rf " + dest + "center.out")
      sys.exit(26) 
 
def countResidue(dest, keyresidue,pdbfilefull, v):
   f = open(dest + keyresidue, "r")
   residues = []
   for l in f:
      if l.startswith("ATOM"):
         l = l.strip()
         chainID = l[21]
         chainID = chainID.strip()
         lindex = l[22:26]
         lindex = lindex.strip()
         info = chainID + lindex
         inlist = 1
         if len(residues) == 0:
            residues.append(info)
         for index in residues:
            if cmp(index, info) == 0:
               inlist=0
               break
         if inlist ==1:
            residues.append(info)
   f.close()
   count = len(residues)
   if (count > 25):
      if  not v:
         os.system("rm -rf " + dest + "initialdata/")
         os.system("rm -rf " + dest + "dataset")
         os.system("rm -rf " + dest + "*.pdb")
         os.system("rm -rf " + dest + pdbfilefull)

      sys.exit(23)
      print "The recommended cleft number should be under 25. Please reduce the number of residues."
      

def preparekeyresidue(dest,wholename,pdbfilefull,keyresidue, linput,logfile,v ):
   #prepare initial protein files
   input = open(dest + pdbfilefull, "r")
   output = open(dest + "initialdata/" + wholename + "_p.pdb", "w")
   for line in input:
      line = line.strip()
      if line.startswith("ATOM"):
         print >>output, line
      if line.startswith("CONECT"):
         print >> output, line
      if line.startswith("MODEL"):
         line = line.split()
         if float(line[1])==2:
            break
   input.close()
   output.close()
   #prepare initial ligand file
    
   output = open(dest + "center.out", "w")
   file = open(dest + keyresidue, "r") 
   xlist = []
   ylist = []
   zlist = []
   for line in file:
      line = line.strip()
      if line.startswith("ATOM"):
         x = float(line[30:38])
         y = float(line[38:46])
         z = float(line[46:54])
         xlist.append(x)
         ylist.append(y)
         zlist.append(z)
   file.close()
   box = []
   minxyz = []
   maxxyz = []
   minxyz.append(min(xlist))
   minxyz.append(min(ylist))
   minxyz.append(min(zlist))
   maxxyz.append(max(xlist))
   maxxyz.append(max(ylist))
   maxxyz.append(max(zlist))
   box.append(minxyz)
   box.append(maxxyz)
   dis = [(b-a) for (a,b) in zip(*box)]
   alocate = [a for (a,b) in zip(*box)]
   blocate = [b for (a,b) in zip(*box)]
   nx = int(dis[0]/1)
   ny = int(dis[1]/1)
   nz = int(dis[2]/1)
   # judgement of volume
   if ((nx*ny*nz)>20000):
      if  not v:
         os.system("rm -rf " + dest + "initialdata/")
         os.system("rm -rf " + dest + "dataset")
         os.system("rm -rf " + dest + "*.pdb")
         os.system("rm -rf " + dest + pdbfilefull)
         os.system("rm -rf " + dest + "center.out")
      sys.exit(21)
      print "The defined cleft is too big for typical ligand binding site. Please reduce the number of residues in the cleft file." 
     # nwater = nx*ny*nz
   if ((nx*ny*nz)==0):
      if  not v:
         os.system("rm -rf " + dest + "initialdata/")
         os.system("rm -rf " + dest + "dataset")
         os.system("rm -rf " + dest + "*.pdb")
         os.system("rm -rf " + dest + pdbfilefull)
         os.system("rm -rf " + dest + "center.out")
      sys.exit(22)
      print "The defined cleft has one dimention 0. Please check the coordinates of residues in the cleft file."
   f = open(dest + "initialdata/" + wholename + "_linitial.pdb", "w")
     # f2 = open(initial_path + index + "twopoints_l.pdb", "w") 
   formatstr = 'HETATM%6d O   HOH%6d    %8.3f%8.3f%8.3f1.00100.00           O  '
     # print >> f2, formatstr % (1, 1, alocate[0], alocate[1], alocate[2])
     # print >> f2, formatstr % (2, 2, blocate[0], blocate[1], blocate[2])
     # f2.close()

   for i in range(nx+1):
      coorx = (blocate[0])-1*(i)
      for j in range(ny+1):
         coory = (blocate[1])-1*(j)
         for k in range(nz+1):
            coorz = (blocate[2])-1*(k)
            print >> f, formatstr %((i+1)*(j+1)*(k+1),(i+1)*(j+1)*(k+1), coorx, coory, coorz)
   print >> f, "END"
   f.close()
   print >> output, wholename + " " + str(nx) + " " + str(ny) + " " + str(nz)
   output.close()

   f = open(dest + "initialdata/" + wholename + "_linitial.pdb", "r")
   fnew = open(dest + "initialdata/" + wholename + "_l.pdb", "w")
   for line in f:
      line = line.strip()
      if line.startswith("HETATM"):
         x = float(line[30:38])
         y = float(line[38:46])
         z = float(line[46:54])
         in_list =1
         ref = open(dest + "initialdata/" + wholename + "_p.pdb", "r")
         for l in ref:
            l = l.strip()
            if l.startswith("ATOM"):
               xref = float(l[30:38])
               yref = float(l[38:46])
               zref = float(l[46:54])
               dis = ((x-xref)**2 + (y-yref)**2 + (z-zref)**2)**0.5
               if (dis > 3.4):
                  continue
               else:
          #     print dis
               #   print xref, yref, zref
                  in_list=0
                  break
         if (in_list ==1):
            print >> fnew, line
         else:
            continue
         ref.close()
   f.close()
   fnew.close()
   filepath = dest + 'initialdata/' + wholename + '_l.pdb'
   if os.stat(filepath).st_size ==0:
      if not v:
         os.system("rm -rf " + dest + "initialdata")
         os.system("rm -rf " + dest + "dataset")
         os.system("rm -rf " + dest + "*.pdb")
         os.system("rm -rf " + dest + pdbfilefull)
         os.system("rm -rf " + dest + "center.out")
      #   os.system("rm -rf " + dest + "initialdata")
      sys.exit(24)
      print "This binding site is too small for a possible CLR binding site. Please enlarge the cleft defined."

def gen_points(dest, wholename, mol2_dbase_path, sitemap_dbase_path, logfile, v, pdbfilefull):
   os.system("/soft/linux64/openeye/bin/molcharge -in " + dest + "initialdata/" + wholename + "_l.pdb -out " + dest + "initialdata/" + wholename + "_l.mol2")
   os.system("mv " + dest + "initialdata/" + wholename + "_l.mol2 " + mol2_dbase_path)
   if os.system("/soft/linux64/SimSite3D_v4.5/bin/gen_points -p " + dest + "initialdata/" + wholename + "_p.pdb -l " + mol2_dbase_path + wholename + "_l.mol2 --allow_small_site_maps --no_normalization --include_metals --msms_surf " + dest +  wholename + "_s.csv")== 0:
      os.system("cp " + dest + wholename + "_s.csv " + sitemap_dbase_path)
      os.system("cp " + dest + wholename + "_s.pdb " + sitemap_dbase_path)
      os.system("cp " + dest + wholename + "_a.pdb " + sitemap_dbase_path)
      os.system("cp " + dest + wholename + "_rad.pdb " + sitemap_dbase_path)
      os.system("cp " + dest + wholename + "_surf.face " + sitemap_dbase_path)
      os.system("cp " + dest + wholename + "_surf.vert " + sitemap_dbase_path)
      print >> logfile, "##         The binding site shape and chemical properties have been extracted."
      print >> logfile, "## ..........................................................................."
   else:
      print >> logfile, "ERROR: Sitemap representation not generated for crstalID: " + wholename + " because of CholMine internal problems."
      print >> logfile, "..........................................................................."
      if  not v:
         os.system("rm -rf " + dest + "initialdata/")
         os.system("rm -rf " + dest + "dataset")
         os.system("rm -rf " + dest + "3KDP_CLR3001D_results/")
         os.system("rm -rf " + dest + "2DYR_CHD525C_results/")
         os.system("rm -rf " + dest + "*.pdb")
         os.system("rm -rf " + dest + "*csv")
         os.system("rm -rf " + dest + "*face")
         os.system("rm -rf " + dest + "*vert")
         os.system("rm -rf " + dest + "rotateddataset")
         os.system("rm -rf " + dest + pdbfilefull)
         os.system("rm -rf " + dest + "center.out")
      sys.exit(14)

def countsitemap(dest, wholename, v, pdbfilefull):
   sfile = open(dest + wholename + "_s.pdb", "r")
   counthphilic = 0
   counts = 0
   for line in sfile:
      line = line.strip()
      if line.startswith("HETATM"):
         counts = counts + 1   
         property = float(line[60:65])
         if (property!=100.00):
            counthphilic = counthphilic + 1
   sfile.close()
   if (counts<20 or (counthphilic < 6 and counts < 30)):
      if  not v:
         os.system("rm -rf " + dest + "initialdata")
         os.system("rm -rf " + dest + "dataset")
         os.system("rm -rf " + dest + "3KDP_CLR3001D_results/")
         os.system("rm -rf " + dest + "2DYR_CHD525C_results/")
         os.system("rm -rf " + dest + "*.pdb")
         os.system("rm -rf " + dest + "*csv")
         os.system("rm -rf " + dest + "*face")
         os.system("rm -rf " + dest + "*vert")
         os.system("rm -rf " + dest + "rotateddataset")
         os.system("rm -rf " + dest + pdbfilefull)
         os.system("rm -rf " + dest + "center.out")
      print "the sitemap points are not enough, please include more residues."
      sys.exit(24)
          
def CLR_searchsites(dest, sitemap_dbase_path, mol2_dbase_path, logfile, wholename, v, pdbfilefull):
   if os.system("search_sitemaps --score_threshold 0.0  --lig_frag_size 0 --prot_lig_scoring --fine_tune_tier2 --dbase_sites " + sitemap_dbase_path + " --dbase_ligs " + mol2_dbase_path + " --proj_output " + dest + "3KDP_CLR3001D_results ./hardcodedfiles/3KDP_CLR3001D_s.csv")==0:
      print >> logfile, "##         The binding site has been compared with the CholMine cholesterol site."
      print >> logfile, "## ..........................................................................."
   else:
      print >> logfile, "ERROR: Site comparison failed for crstalID: " + wholename + " because of CholMine internal problems."
      print >> logfile, "..........................................................................."
      if  not v:
         os.system("rm -rf " + dest + "initialdata/")
         os.system("rm -rf " + dest + "dataset")
         os.system("rm -rf " + dest + "3KDP_CLR3001D_results/")
         os.system("rm -rf " + dest + "2DYR_CHD525C_results/")
         os.system("rm -rf " + dest + "*.pdb")
         os.system("rm -rf " + dest + "*csv")
         os.system("rm -rf " + dest + "*face")
         os.system("rm -rf " + dest + "*vert")
         os.system("rm -rf " + dest + "rotateddataset")
         os.system("rm -rf " + dest + pdbfilefull)
         os.system("rm -rf " + dest + "center.out")
      sys.exit(15)

      
def CHD_searchsites(dest, sitemap_dbase_path, mol2_dbase_path, logfile, wholename,v, pdbfilefull):
   if os.system("search_sitemaps --score_threshold 0.0  --lig_frag_size 0 --prot_lig_scoring --fine_tune_tier2 --dbase_sites " + sitemap_dbase_path + " --dbase_ligs " + mol2_dbase_path + " --proj_output " + dest + "2DYR_CHD525C_results ./hardcodedfiles/2DYR_CHD525C_s.csv")==0:
      print >> logfile, "##         The binding site has been compared with CholMine cholate site."
      print >> logfile, "## ..........................................................................."
   else:
      print >> logfile, "Error: Site comparison failed for crstalID: " + wholename + " because of CholMine internal problems."
      print >> logfile, "..........................................................................."
      if  not v:
         os.system("rm -rf " + dest + "initialdata/")
         os.system("rm -rf " + dest + "dataset")
         os.system("rm -rf " + dest + "3KDP_CLR3001D_results/")
         os.system("rm -rf " + dest + "2DYR_CHD525C_results/") 
         os.system("rm -rf " + dest + "*.pdb")
         os.system("rm -rf " + dest + "*csv")
         os.system("rm -rf " + dest + "*face")
         os.system("rm -rf " + dest + "*vert")
         os.system("rm -rf " + dest + "rotateddataset")
         os.system("rm -rf " + dest + pdbfilefull)
         os.system("rm -rf " + dest + "center.out")
      sys.exit(15)
def CLRmatchprint(dest):
   matchprint_file = open(dest + '3KDP_CLR3001D_results/3KDP_CLR3001D_matchprint.out', 'w')
   res_file = open(dest + '3KDP_CLR3001D_results/3KDP_CLR3001D_blast.out', 'r')
   for line in res_file:
      if(line.startswith("#")): continue
      else:
         toks = line.split("|")
         matchprint = toks[4]
         print >> matchprint_file, matchprint, toks[0], toks[1]
   matchprint_file.close()
   res_file.close()

def CHDmatchprint(dest):
   matchprint_file = open(dest + '2DYR_CHD525C_results/2DYR_CHD525C_matchprint.out', 'w')
   res_file = open(dest + '2DYR_CHD525C_results/2DYR_CHD525C_blast.out', 'r')
   for line in res_file:
      if(line.startswith("#")): continue
      else:
         toks = line.split("|")
         matchprint = toks[4]
         print >> matchprint_file, matchprint, toks[0], toks[1]
   matchprint_file.close()
   res_file.close()

def matchprint_analysis(file_path, file_path3, logfile, pdbfile, chainid, ligandindex, linput, dest, v, pdbfilefull):
   matrix = []
   input = open(file_path, 'r')
   for line in input:
      line = line.strip()
      toks = line.split()
      matrix.append(toks[0])
   input.close()
   #figure out the columns which satifies 85% conservation and put (column number-1) to columnlist
   columnlist = []
   columnlist2 = []
   for i in range(len(matrix[0])):
      sum = 0
      for j in range(len(matrix)):
         sum = sum + int(matrix[j][i])
      if (sum < len(matrix)*0.70): # conservation 85%
         continue
      else:
         columnlist.append(i)
         number = i + 1
         columnlist2.append(number)#need to plus 1 for column number starting from 1 instead of 0

   #test set
   input3 = open(file_path3, 'r')# open test matchprint file
   matrix_test = []
   proteinidtest = []
   idtest = []
   for line in input3:
      if line == "":
         break

      line = line.strip()
      toks = line.split()
      matrix_test.append(toks[0])
      proteinidtest.append(toks[1])
   row_testlist = []
   for j in range(len(matrix_test)):
      sum = 0
      for index in columnlist:
         sum = sum + int(matrix_test[j][index])
      if (sum < len(columnlist)*0.70):
         continue
      else:
         row_testlist.append(j)
   input3.close()
  # print >> logfile, "How many sites are left after screening: ", len(row_testlist)
   marker = chainid + " " + ligandindex
   marker = marker.strip()
   if len(row_testlist) ==0:
         print >> logfile, "## Step 4: Prediction result:"
         if cmp(linput[:3], "CLR")==0:
            print >> logfile, "##         CholMine reports that the site at the " + marker  + " in " + pdbfile + " was not predicted as a cholesterol binding site after checking conserved interactions with the prototypic site."
         if cmp(linput[:3], "CHD")==0:
            print >> logfile, "##         CholMine reports that the site at the " + marker  + " in " + pdbfile + " was not predicted as a cholate binding site after checking conserved interactions with the prototypic site."
         if not v:
            os.system("rm -rf " + dest + "initialdata")
            os.system("rm -rf " + dest + "dataset")
            os.system("rm -rf " + dest + "3KDP_CLR3001D_results/")
            os.system("rm -rf " + dest + "2DYR_CHD525C_results/")
            os.system("rm -rf " + dest + "*.pdb")
            os.system("rm -rf " + dest + "*csv")
            os.system("rm -rf " + dest + "*face")
            os.system("rm -rf " + dest + "*vert")
            os.system("rm -rf " + dest + "rotateddataset")
            os.system("rm -rf " + dest + pdbfilefull)
            os.system("rm -rf " + dest + "center.out")
         if cmp(linput[:3], "CLR")==0:
            sys.exit(16)
         if cmp(linput[:3], "CHD")==0:
            sys.exit(17)
   print >> logfile, "## Step 4: Prediction result:"
   if cmp(linput[:3], "CLR")==0:
      print >> logfile, "##         Cholesterol site predicted!"
      print >> logfile, "##         Based on similarities in shape and conserved interactions, CholMine reports that the site at the " + marker  + " in " + pdbfile+" matches known cholesterol site features."
   if cmp(linput[:3], "CHD")==0:
      print >> logfile, "##         Cholate site predicted!"
      print >> logfile, "##         Based on similarities in shape and conserved interactions, CholMine reports that the site at the " + marker  + " in " + pdbfile+" matches known cholate site features."


#   print >> logfile, "##         ID of unknown site picked by CholMine: "
   for i in range(len(row_testlist)):
 #     print >> logfile, "##         " + proteinidtest[row_testlist[i]][:-13]
      idtest.append(proteinidtest[row_testlist[i]])
   return idtest 
  
def rotatedresidues(dest, input,idtest,logfile,wholename,label,linput):
   from ASCbasePy.utils import pdb
   os.system("mkdir " + dest + "rotateddataset")
   output_path =  dest + "rotateddataset/"
   coormatrix = []
   for line in input:
      line = line.strip()
      if line.startswith("HETATM"):
         keycoor = []
         x = float(line[30:38])
         y = float(line[38:46])
         z = float(line[46:54])
         keycoor.append(x)
         keycoor.append(y)
         keycoor.append(z)
         coormatrix.append(keycoor)
   input.close()

   for index in idtest:
      input = open(dest + label + "_results/" + label + "_blast.out", "r")

      for line in input:
         if(line.startswith("#")): continue
         else:
            toks = line.split("|")
            if cmp(index, toks[0])==0:
               M = toks[2] + " " + toks[3]
               M = [float(s) for s in M.split()]
               R = np.array(M[:9]).reshape(3,3)
               T = np.array(M[9:])
               pdb.transform(dest + wholename + "_s.pdb", output_path + index[:-13] + "_rotateds.pdb", R.T, T)
               pdb.transform(dest + "initialdata/" + wholename + "_p.pdb", output_path + index[:-13] + "_rotatedp.pdb", R.T, T)
  
   keypointindex = []
   for index in idtest:
      keyeach = []
      for i in range(len(coormatrix)):
         dis = []
         dex = []
         chem = []
         file = open(output_path + index[:-13] + "_rotateds.pdb", "r")
         for line in file:
            line = line.strip()
            if line.startswith("HETATM"):
               x = float(line[30:38])
               y = float(line[38:46])
               z = float(line[46:54])
      #         type = float(line[60:66])
             #  if type == 100.00:
               distance =((x-float(coormatrix[i][0]))**2 + (y -float(coormatrix[i][1]))**2 + (z- float(coormatrix[i][2]))**2)**0.5
                 #print distance
               if (distance > 1.5):
                  continue
               else:
          #        print distance               
                  dis.append(distance)
                  pointindexall = int(line[20:26])
                  type = float(line[60:66])
                  dex.append(pointindexall)
                  chem.append(type)
         file.close()
         if len(dis) == 0:
            continue
         else:
            locate = dis.index(min(dis))
            pointindex = dex[locate]
            chemtype = chem[locate]
            if chemtype == 100.00:
               innum =1
               if (len(keyeach) == 0):
                  keyeach.append(pointindex)
                  continue
               for ins in keyeach:
                  if (ins == pointindex):
                     innum = 0
                     break
               if (innum == 1):
                  keyeach.append(pointindex)
      keypointindex.append(keyeach)
   resimatrix = []
   for i in range(len(idtest)):
      resieach = []
      f = open (dest + wholename + "_s.csv", "rb")
      reader = csv.reader(f, delimiter='|')
      for line in reader:
         if line == []:
            continue
         else:
            if len(line) == 3:
#            print line
               num = line[0]
               for j in range(len(keypointindex[i])):
                  numrefer = keypointindex[i][j]
                  if cmp(num, str(numrefer)) == 0:
                     resi = line[1]
                     resieach.append(resi)
      resimatrix.append(resieach)
      f.close()

   #newmatrix contains residues and atoms making conserved interactions
   newmatrix = []
   for i in range(len(idtest)):
      residues = []
      for j in range(len(resimatrix[i])):
         test = resimatrix[i][j].split("_")
         #print test
         for k in range(int(len(test)/3)):
            resi = test[(k+1)*3-3] + "_"+ test[(k+1)*3-2] + "_" + test[(k+1)*3-1]
            in_list =1
            if len(residues) == 0:
               residues.append(resi)
            for index in residues:
               if cmp(index, resi) == 0:
                  in_list=0
                  break
            if in_list ==1:
               residues.append(resi)
      newmatrix.append(residues)
   #print newmatrix

   #newmatrix2 only contains residues making conserved interactions
   newmatrix2 = []
   for i in range(len(idtest)):
      residues2 = []
      for j in range(len(resimatrix[i])):
         test = resimatrix[i][j].split("_")
         for k in range(int(len(test)/3)):
            resi = test[(k+1)*3-3]+ "_"+ test[(k+1)*3-2]
            if test[(k+1)*3-2] == "":
               resi = test[(k+1)*3-3]+ "_"+ " "
            in_list =1
            if len(residues2) == 0:
               residues2.append(resi)
            for index in residues2:
               if cmp(index, resi) == 0:
                  in_list=0
                  break
            if in_list ==1:
               residues2.append(resi)
      newmatrix2.append(residues2)
   #print newmatrix2

   print >> logfile, "## ..........................................................................."
   if cmp(linput[:3], "CLR")==0:
      print >> logfile, "## The residues making conserved interactions with cholesterol:"
   if cmp(linput[:3], "CHD")==0:
      print >> logfile, "## The residues making conserved interactions with cholate:"
   print >> logfile, "## SiteID: Residues"

#print >> finalfile, "#ID|important residues|"
   for i in range(len(idtest)):
      for j in range(len(newmatrix[i])):
         test = newmatrix[i][j]
       #  all =  str(test) + " " + all
         print >> logfile, idtest[i][:-13] + ": " + test
   print >> logfile, "## Note that in the l.pdb file (predicted ligand position) and PyMOL .pse file (if requested), the position of the steroid tail is arbitrary, as it is outside the conserved recognition motif."
   logfile.close()
#read rotated pdb file

   for i in range(len(idtest)):
      pdb = open(output_path + idtest[i][:-13] + "_rotatedp.pdb", "r")
      pdbout = open(output_path + idtest[i][:-13] + "_rotatedkeyresidue.pdb", "w")
      for line in pdb:
         line=line.strip()
         if line.startswith("ATOM"):
            resiname = line[17:20]
            chainID = line[21]
            resinum = int(line[22:26])
            resiwhole = resiname + str(resinum) + "_" + chainID
       # print resiwhole
            for j in range(len(newmatrix2[i])):
               test = newmatrix2[i][j]
               if cmp(resiwhole, test)==0:
                  print >> pdbout, line

         if line.startswith("CONECT"):
            print  >> pdbout, line
      pdb.close()
      pdbout.close()


def rotatedback(in_fname, out_fname, R, T):
   """
  Assumptions 
    Assumes R is for premultiplication by coordinates and 
    R and T are 3x3 and 3x1 numpy arrays 
    apply transformation to all ATOM and HETATM records
    transformation is x' =R.T*x + T
  According to the equation, if the coordinates are rotated back, then x = (x'-T)*(1/(R.T))
  """
   infile = open(in_fname, "r")
   out = open(out_fname, "w+")
   fmat = '%8.3f%8.3f%8.3f'
   for line in infile:
      if(line.startswith("ATOM  ") or line.startswith("HETATM")):
         position= [float(line[30:38]), float(line[38:46]), float(line[46:54])]
         diff = position-T
         position= np.dot(diff, R.I)
       #  position =np.dot((position-T), R.I)
         x = np.array(position)[0][0]
         y =  np.array(position)[0][1]
         z =   np.array(position)[0][2]
         out.write(line[0:30]+fmat%(x, y, z)+line[54:])
      else:
         out.write(line)

def CLR_rotatedback_proteins(dest, input,idtest,logfile,wholename,label):
   for index in idtest:
      input = open(dest + label + "_results/" + label + "_blast.out", "r")
      output_path = dest + "rotateddataset/"
      for line in input:
         if(line.startswith("#")): continue
         else:
            toks = line.split("|")
            if cmp(index, toks[0])==0:
               M = toks[2] + " " + toks[3]
               M = [float(s) for s in M.split()]
               R = np.matrix([[M[0],M[1],M[2]],[M[3],M[4],M[5]],[M[6],M[7],M[8]]])
               T = np.matrix([M[9],M[10],M[11]])
              # os.system("cat " + output_path + index[:-13] + "_rotatedp.pdb ./hardcodedfiles/3KDP_CLR3001D_l.pdb > " + output_path + index[:-13] + "_addligandp.pdb" )
               rotatedback(output_path + index[:-13] + "_rotatedp.pdb",output_path + index[:-13] + "_p.pdb", R.T, T)
              # rotatedback(output_path + index[:-13] + "_addligandp.pdb",output_path + index[:-13] + "_p.pdb", R.T, T)
               rotatedback("./hardcodedfiles/3KDP_CLR3001D_l.pdb",output_path + index[:-13] + "_l.pdb", R.T, T)
               rotatedback(output_path + index[:-13] + "_rotatedkeyresidue.pdb",output_path + index[:-13] + "_keyres.pdb", R.T, T)
           #    rotatedback("./hardcodedfiles/conservedsitemappoints3kdp.pdb", output_path + index[:-13] + "_keypoints.pdb", R.T, T)
def CHD_rotatedback_proteins(dest, input,idtest,logfile,wholename,label):
   for index in idtest:
      input = open(dest + label + "_results/" + label + "_blast.out", "r")
      output_path = dest + "rotateddataset/"
      for line in input:
         if(line.startswith("#")): continue
         else:
            toks = line.split("|")
            if cmp(index, toks[0])==0:
               M = toks[2] + " " + toks[3]
               M = [float(s) for s in M.split()]
               R = np.matrix([[M[0],M[1],M[2]],[M[3],M[4],M[5]],[M[6],M[7],M[8]]])
               T = np.array(M[9:])
               os.system("cat " + output_path + index[:-13] + "_rotatedp.pdb ./hardcodedfiles/2DYR_CHD525C_l.pdb > " + output_path + index[:-13] + "_addligandp.pdb" )
               rotatedback(output_path + index[:-13] + "_rotatedp.pdb",output_path + index[:-13] + "_p.pdb", R.T, T)
              # rotatedback(output_path + index[:-13] + "_addligandp.pdb",output_path + index[:-13] + "_p.pdb", R.T, T)
               rotatedback("./hardcodedfiles/2DYR_CHD525C_l.pdb",output_path + index[:-13] + "_l.pdb", R.T, T)
               rotatedback(output_path + index[:-13] + "_rotatedkeyresidue.pdb",output_path + index[:-13] + "_keyres.pdb", R.T, T)
            #   rotatedback("./hardcodedfiles/conservedsitemappoints2DYR.pdb", output_path + index[:-13] + "_keypoints.pdb", R.T, T)

def pseout(dest, logread,wholename,resn):
   import __main__
   __main__.pymol_argv = ['pymol','-qc'] # Pymol: quiet and no GUI
   from time import sleep
   import pymol
   pymol.finish_launching()
   pymol.cmd.do("load " + dest + wholename + ".pdb")
   pymol.cmd.do("load " + dest + wholename + "_l.pdb")
   pymol.cmd.do("hide ")
   pymol.cmd.do("show cartoon, " + wholename)
   pymol.cmd.do("cartoon loop")
   pymol.cmd.do("show sticks, " + wholename + "_l")
   pymol.cmd.do("zoom resn " + resn) 
  # pymol.cmd.do("show surface, " + wholename)
  # pymol.cmd.do("set surface_cavity_mode, 1")
  # pymol.cmd.do("set cartoon_oval_quality, 2")
  # pymol.cmd.do("set cartoon_tube_radius, 0.075")
  # pymol.cmd.do("set cartoon_loop_radius, 0.075")
   pymol.cmd.do("bg_color white")
  # pymol.cmd.do("load " + dest + wholename + "_keypoints.pdb")
  # pymol.cmd.do("color magenta, " + wholename + "_keypoints")
  # pymol.cmd.do("show spheres, " + wholename + "_keypoints")
   pymol.cmd.do("load " + dest + wholename + "_keyres.pdb")
   pymol.cmd.do("show sticks, " + wholename + "_keyres")
   pymol.cmd.do("show surface, " + wholename + "_keyres")
   pymol.cmd.do("set line_width, 2")
   pymol.cmd.do("util.cbaw " + wholename + "_keyres")
 #  pymol.cmd.do("zoom " + wholename + "_keyres") 
   pymol.cmd.do("set transparency, 0.5")
   for line in logread:
      line = line.strip()
      if line.startswith("##"):
         continue
      toks = line.split(":")
      id = toks[0]
      toks[1] = toks[1].strip()
      info = toks[1].split("_")
      chain = info[1]
      rname = info[0][:3]
      rnum = info[0][3:]
      atype = info[2]
      pymol.cmd.do("select a, /" + id +  "_keyres//" + chain + "/" + rname + "`" + rnum+ "/"+atype)
      pymol.cmd.do("show spheres, /" + id +  "_keyres//" + chain + "/" + rname + "`" + rnum+ "/"+atype)
      pymol.cmd.do("color green, /" + id +  "_keyres//" + chain + "/" + rname + "`" +rnum+ "/"+atype)
    #  ref = open(dest + wholename + "_keypoints.pdb", "r")
    #  i =0
    #  for line2 in ref:
     #    i = i+1
      #   if line2.startswith("HETATM"):
       #     Onum = line2[24:26].strip()
        #    print Onum
         #   pymol.cmd.distance("distance"+str(i), "a", "/" + wholename + "_keypoints///HOH`" + Onum + "/O", 4.0)
     # ref.close()
   logread.close()
  # pymol.cmd.do("hide labels, distance*")
  # pymol.cmd.do("set dash_gap, 0.5")
  # pymol.cmd.do("set dash_radius, 0.1")
   pymol.cmd.do("set sphere_scale, 0.3")
   pymol.cmd.do("set surface_quality, 1")
   pymol.cmd.do("set solvent_radius, 2")
   pymol.cmd.do("delete a")
   pymol.cmd.save(dest + wholename + ".pse")     
   
def validPDB(pdbfilefull, ligandindex):
   #Judge whether the pdb file is valid pdb file
   input = open(pdbfilefull, "r")
   atom=0
   hetatm = 0
   for line in input:
      line = line.strip()
      toks = line.split()
      if len(toks)==0:
         continue
      firstword=toks[0]
      if firstword.startswith("ATOM"):
         atom=atom+1
      if firstword.startswith("HETATM"):
         hetatm=hetatm+1
   input.close()
   if atom==0:
      print "ERROR: Not a valid PDB file; missing protein information (ATOM lines)."
      sys.exit(6)
   if hetatm==0:
      if cmp(ligandindex, "cleft")!=0:
         print "ERROR: Not a valid PDB file; missing ligand information (HETATM lines)."
         sys.exit(7)
   input = open(pdbfilefull, "r")
   for line in input:
      line=line.strip()
      if line.startswith("ATOM"):
         if len(line) < 66:
            print "ERROR: Not a valid PDB file; missing protein information (ATOM lines)."
            sys.exit(6)
      if line.startswith("HETATM"):
         if len(line) < 66:
            print "ERROR: Not a valid PDB file; missing ligand information (HETATM lines)."
            sys.exit(7)
    
def validCleft(keyresidue):
   input = open(keyresidue, "r")
   atom=0
   for line in input:
      line = line.strip()
      toks = line.split()
      if len(toks)==0:
         continue
      firstword=toks[0]
      if firstword.startswith("ATOM"):
         atom=atom+1
   input.close()
   if atom==0:
      print "ERROR: Not a valid Cleft PDB file; missing cleft residue information (ATOM lines)."
      sys.exit(25)
   input = open(keyresidue, "r")
   for line in input:
      line=line.strip()
      if line.startswith("ATOM"):
         if len(line) < 66:
            print "ERROR: Not a valid Cleft PDB file; missing cleft residue information (ATOM lines)."
            sys.exit(25)

def main():
   if len(sys.argv) < 6:
      print "ERROR: The command line has no enough arguments. Please check the usage."
      print usage_help()
      sys.exit(19)
   args=usage_args()
  # if not args.pname or not args.lname or not args.lnum:
   if not args.pfile:
      print "ERROR: The command line is missing an argument parameter for protein file. Please upload protein file."
      print usage_help()
      sys.exit(1)
   if not args.lnum and not args.keyresidue:
      print "ERROR: The command line is missing an argument parameter for ligand residue number or cleft pdb file. Please enter ligand residue number."
      sys.exit(2)
   if not args.CHD and not args.CLR:
      print "ERROR: Please specify cholesterol or cholate predictor."
      print usage_help()
      sys.exit(3)
   pdbfilefull= args.pfile
  # ligandcode= args.lname
   ligandindex=args.lnum
   projpath = args.projout
   keyresidue = args.keyresidue
   if args.lnum and args.keyresidue:
      print "lnum and keysidue can not be turned on at the same time."
      sys.exit(20)
   if args.CLR:
      if args.lnum:
         if not args.lchain or args.lchain ==" ":
            chainid = ""
            linput = "CLR_" + "" + "_" + ligandindex.strip()
         else:
            chainid=args.lchain
            if not chainid.isalpha():
               print "ERROR: Not a valid chain ID."
               sys.exit(4)
            linput = "CLR_" + chainid.strip() + "_" + ligandindex.strip()
      if args.keyresidue:
         chainid = "" 
         ligandindex="cleft"
         linput = "CLR_" + chainid.strip() + "_" + ligandindex.strip()
   if args.CHD:
      if args.lnum:
         if not args.lchain or args.lchain ==" ":
            chainid = ""
    #  linput = ligandcode.strip() + ligandindex.strip()+ "A"
            linput = "CHD_" + "" + "_" + ligandindex.strip()
         else:
            chainid=args.lchain
            if not chainid.isalpha():
               print "ERROR: Not a valid chain ID."
               sys.exit(4)
     # linput = ligandcode.strip() + ligandindex.strip() + chainid.strip()
            linput = "CHD_" + chainid.strip() + "_" + ligandindex.strip()
      if args.keyresidue:
         chainid = ""
         ligandindex="cleft"
         linput = "CHD_" + "" + "_" + ligandindex.strip()
   num = pdbfilefull.count(".")
   if num ==0:
      pdbfile = pdbfilefull
      wholename = pdbfilefull + "_" + linput
   elif num>=1:
      for i in range(len(pdbfilefull)):
         if cmp(pdbfilefull[i], ".")==0:
            indexnum =i  
      pdbfile = pdbfilefull[0:indexnum]
      wholename = pdbfilefull[0:indexnum] + "_" + linput
   
   # Judge whether the file exist or not"
   v = args.v
   if os.path.isfile(pdbfilefull) == False:
      print "ERROR: The PDB file doesn't exist."
      sys.exit(5)

   #Judge whether the pdb file is valid pdb file
   print ligandindex
   validPDB(pdbfilefull, ligandindex)
   if args.keyresidue:
      validCleft(keyresidue)
   marker = chainid + " " + ligandindex 
   marker = marker.strip() 
   #build directory 
   dest = directory(projpath, pdbfile )
   mol2_dbase_path = dest + "dataset/ligands/"
   sitemap_dbase_path =  dest + "dataset/dbase/"
     
   logfile = open(dest + wholename + "_results.txt", "w")
  # print >> logfile, str(len(sys.argv))
   os.system("cp " +  pdbfilefull + " " + dest)
   if args.keyresidue:
      os.system("cp " +  keyresidue + " " + dest)
   if args.CLR:
      if (args.CLR and args.CHD):
         print "ERROR: Cholesterol and cholate predictor cannot be turned on at the same time; do successive runs."
         sys.exit(8)
      print "CLR predictor is turned on."
      print >> logfile, "## " + str(datetime.now())
      print >> logfile, "## CholMine server 1.00"
      print >> logfile, "## CLR predictor is turned on."
      print >> logfile, "## The binding site at the " + marker + " in " + pdbfile + " is being analyzed."
      if args.lnum:
         print >> logfile, "## Step 1: Read protein and ligand information and generate initial files for the binding site."
         print >> logfile, "## ..........................................................................."
         prepare(dest,wholename,pdbfilefull,linput,logfile,v)
         print >> logfile, "##         Protein and ligand information have been read."
      if args.keyresidue:
         print >> logfile, "## Step 1: Read protein and cleft information and generate initial files for the binding site."
         print >> logfile, "## ..........................................................................."
         checksubset(dest, keyresidue, pdbfilefull, v)
         countResidue(dest, keyresidue,pdbfilefull, v)
         preparekeyresidue(dest,wholename,pdbfilefull,keyresidue, linput,logfile,v)
         print >> logfile, "##         Protein and cleft information have been read."

      print >> logfile, "## ..........................................................................."
      print >> logfile, "## Step 2: Extract the shape and chemical properties in the binding site of " + wholename + "."
      print >> logfile, "## ..........................................................................."
      gen_points(dest, wholename, mol2_dbase_path, sitemap_dbase_path, logfile, v, pdbfilefull)
      countsitemap(dest, wholename, v, pdbfilefull)
      print >> logfile, "## Step 3: Compare the binding site shape and chemical properties with CholMine prototypic cholesterol site."
      print >> logfile, "## ..........................................................................."
      CLR_searchsites(dest, sitemap_dbase_path, mol2_dbase_path, logfile, wholename,v, pdbfilefull)
      CLRmatchprint(dest)
      file_path = './hardcodedfiles/3KDP_CLR3001D_matchprint.out'
      file_path3 = dest + '3KDP_CLR3001D_results/3KDP_CLR3001D_matchprint.out'
      if os.path.isfile(file_path3) == False:
         print >> logfile,  "ERROR: The CholMine cholesterol matchprint file doesn't exist."
         if  not args.v:
            os.system("rm -rf " + dest + "initialdata")
            os.system("rm -rf " + dest + "dataset")
            os.system("rm -rf " + dest + "3KDP_CLR3001D_results/")
            os.system("rm -rf " + dest + "2DYR_CHD525C_results/")
            os.system("rm -rf " + dest + "*.pdb")
            os.system("rm -rf " + dest + "*csv")
            os.system("rm -rf " + dest + "*face")
            os.system("rm -rf " + dest + "*vert")
            os.system("rm -rf " + dest + "rotateddataset")
            os.system("rm -rf " + dest + pdbfilefull)
            os.system("rm -rf " + dest + "center.out")
         sys.exit(9)
      if os.stat(file_path3).st_size ==0:
         print >> logfile, "## Step 4: Prediction result:"
         print >> logfile, "##         CholMine reports that the ligand site at the " + marker + " in " + pdbfile + " was not predicted as a cholesterol binding site after site alignment to the CholMine prototypic site."
         if  not args.v:
            os.system("rm -rf " + dest + "initialdata")
            os.system("rm -rf " + dest + "dataset")
            os.system("rm -rf " + dest + "3KDP_CLR3001D_results/")
            os.system("rm -rf " + dest + "2DYR_CHD525C_results/")
            os.system("rm -rf " + dest + "*.pdb")
            os.system("rm -rf " + dest + "*csv")
            os.system("rm -rf " + dest + "*face")
            os.system("rm -rf " + dest + "*vert")
            os.system("rm -rf " + dest + "rotateddataset")
            os.system("rm -rf " + dest + pdbfilefull)
            os.system("rm -rf " + dest + "center.out")
         sys.exit(10)
      idtest = matchprint_analysis(file_path, file_path3, logfile, pdbfile, chainid, ligandindex, linput,dest,v, pdbfilefull)
      input = open("./hardcodedfiles/conservedsitemappoints3kdp.pdb", "r")
      label="3KDP_CLR3001D"
      rotatedresidues(dest, input,idtest,logfile,wholename, label, linput)
      CLR_rotatedback_proteins(dest, input,idtest,logfile,wholename,label)
   elif args.CHD:
      if (args.CLR and args.CHD):
         print "Cholesterol and cholate predictor cannot be turned on at the same time."
         sys.exit(8)
      print "CHD predictor turned on."
      print >> logfile, "## " + str(datetime.now())
      print >> logfile, "## CholMine server 1.00"
      print >> logfile, "## CHD predictor is turned on."
      print >> logfile, "## The binding site at the " + marker + " in " + pdbfile + " is being analyzed."
      if args.lnum:
         print >> logfile, "## Step 1: Read protein and ligand information and generate initial files for the binding site."
         print >> logfile, "## ..........................................................................."
         prepare(dest,wholename,pdbfilefull,linput,logfile,v)
         print >> logfile, "##         Protein and ligand information have been read."
      if args.keyresidue:
         print >> logfile, "## Step 1: Read protein and cleft information and generate initial files for the binding site."
         print >> logfile, "## ..........................................................................."
         checksubset(dest, keyresidue, pdbfilefull, v)
         countResidue(dest, keyresidue,pdbfilefull, v)
         preparekeyresidue(dest,wholename,pdbfilefull,keyresidue, linput,logfile,v)
         print >> logfile, "##         Protein and cleft information have been read."
      print >> logfile, "## ..........................................................................."
      print >> logfile, "## Step 2: Extract the shape and chemical properties in the binding site of " + wholename + "."
      print >> logfile, "## ..........................................................................."
      gen_points(dest, wholename, mol2_dbase_path, sitemap_dbase_path, logfile,v, pdbfilefull)
      countsitemap(dest, wholename, v, pdbfilefull)
      print >> logfile, "## Step 3: Compare the binding site shape and chemical properties with CholMine prototypic cholate site."
      print >> logfile, "## ..........................................................................."
      CHD_searchsites(dest, sitemap_dbase_path, mol2_dbase_path, logfile, wholename, v, pdbfilefull)
      CHDmatchprint(dest)
      file_path = './hardcodedfiles/2DYR_CHD525C_matchprint.out'
      file_path3 = dest + '2DYR_CHD525C_results/2DYR_CHD525C_matchprint.out'
      if os.path.isfile(file_path3) == False:
         print >> logfile,  "ERROR: the CholMine cholate matchprint file doesn't exist."
         if  not args.v:
            os.system("rm -rf " + dest + "initialdata")
            os.system("rm -rf " + dest + "dataset")
            os.system("rm -rf " + dest + "3KDP_CLR3001D_results/")
            os.system("rm -rf " + dest + "2DYR_CHD525C_results/")
            os.system("rm -rf " + dest + "*.pdb")
            os.system("rm -rf " + dest + "*csv")
            os.system("rm -rf " + dest + "*face")
            os.system("rm -rf " + dest + "*vert")
            os.system("rm -rf " + dest + "rotateddataset")
            os.system("rm -rf " + dest + pdbfilefull)
            os.system("rm -rf " + dest + "center.out")
         sys.exit(11)
      if os.stat(file_path3).st_size ==0:
         print >> logfile, "## Step 4: Prediction result:"
         print >> logfile, "##         CholMine reports that the ligand site at the " + marker + " in " + pdbfile + " was not predicted as a cholate binding site after site alignment to the CholMine prototypic site."
         if  not args.v:
            os.system("rm -rf " + dest + "initialdata")
            os.system("rm -rf " + dest + "dataset")
            os.system("rm -rf " + dest + "3KDP_CLR3001D_results/")
            os.system("rm -rf " + dest + "2DYR_CHD525C_results/")
            os.system("rm -rf " + dest + "*.pdb")
            os.system("rm -rf " + dest + "*csv")
            os.system("rm -rf " + dest + "*face")
            os.system("rm -rf " + dest + "*vert")
            os.system("rm -rf " + dest + "rotateddataset")
            os.system("rm -rf " + dest + pdbfilefull)
            os.system("rm -rf " + dest + "center.out")
         sys.exit(12)
      idtest = matchprint_analysis(file_path, file_path3, logfile, pdbfile, chainid, ligandindex, linput,dest,v,pdbfilefull)
      input = open("./hardcodedfiles/conservedsitemappoints2DYR.pdb", "r")
      label="2DYR_CHD525C"
      rotatedresidues(dest, input,idtest,logfile,wholename,label,linput)
      CHD_rotatedback_proteins(dest, input,idtest,logfile,wholename,label)
   logfile.close()
   if  not args.v:
      os.system("rm -rf " + dest + "initialdata")
      os.system("rm -rf " + dest + "dataset")
      os.system("rm -rf " + dest + "3KDP_CLR3001D_results/")
      os.system("rm -rf " + dest + "2DYR_CHD525C_results/")
      os.system("rm -rf " + dest + "*.pdb")
      os.system("rm -rf " + dest + "*csv")
      os.system("rm -rf " + dest + "*face")
      os.system("rm -rf " + dest + "*vert")
      os.system("rm -rf " + dest + pdbfilefull)
      os.system("rm -rf " + dest + "center.out")
      os.system("cp " + dest + "rotateddataset/"+"*_p.pdb " + dest + wholename + ".pdb")
      os.system("cp " + dest + "rotateddataset/"+"*_l.pdb " + dest + wholename + "_l.pdb")
      os.system("cp " + dest + "rotateddataset/"+"*_keyres.pdb " + dest)
  #    os.system("cp " + dest + "rotateddataset/"+"*_keypoints.pdb " + dest)
      os.system("rm -rf " + dest + "rotateddataset")
   if args.opse:
      logread = open(dest + wholename + "_results.txt", "r")
      if args.CLR:
         resn = "CLR"
         pseout(dest, logread,wholename,resn)
         sys.exit(0)
      if args.CHD:
         resn = "CHD"
         pseout(dest, logread,wholename,resn)
         sys.exit(18)
   if not args.opse:
      if args.CLR:
         sys.exit(0)
      if args.CHD:
         sys.exit(18)

main()







