
import argparse
import sys
import json
import re
import datetime
import os
import math, csv
from   os import  path
from   Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import operator
from   collections import defaultdict


arguments = argparse.ArgumentParser(description='Summarize selection analysis results.')

arguments.add_argument('-s', '--site',   help = 'Site to extract', required = True, type = int, nargs = '*')
arguments.add_argument('-d', '--dir',   help = 'Directory of result files', required = True, type = str)
import_settings = arguments.parse_args()

dir_path = import_settings.dir

def make_path (file):
    return os.path.join (dir_path, file)

def load_slac (file, sites, results, tag):
    def extract_site_info (site, results):
        residues = {}
        subs = {}
        for branch, tag in slac_data['tested']["0"].items():
            if tag == 'test':
                aa = slac_data["branch attributes"]["0"][branch]["amino-acid"][0][site-1]
                if aa in residues:
                    residues[aa] += 1
                else:
                    residues[aa] = 1
        results['info'] = slac_data["MLE"]["content"]["0"]["by-site"]["RESOLVED"][site-1]
        results['residues'] = residues

    with open (file, 'r') as fh:
        slac_data = json.load (fh)
        for s in sites:
            if not s in results:
                results[s] = {}
            if not tag in results[s]:
                results[s][tag] = {}
            results[s][tag]['slac'] = {}
            extract_site_info (s, results[s][tag]['slac'])
        
    return results
    
def load_fel (file, sites, results, tag):

    with open (file, 'r') as fh:
        fel_data = json.load (fh)
        for s in sites:
            if not s in results:
                results[s] = {}
            if not tag in results[s]:
                results[s][tag] = {}
            results[s][tag]['fel'] = fel_data["MLE"]["content"]["0"][s-1]
        
    return results

def load_meme (file, sites, results, tag):

    with open (file, 'r') as fh:
        meme_data = json.load (fh)
        for s in sites:
            if not s in results:
                results[s] = {}
            if not tag in results[s]:
                results[s][tag] = {}
            results[s][tag]['meme'] = meme_data["MLE"]["content"]["0"][s-1]
        
    return results
 
def load_cfel (file, sites, results):
    with open (file, 'r') as fh:
        cfel_data = json.load (fh)
        for s in sites:
            if not s in results:
                results[s] = {}
            results[s]['cfel'] = {
                'overall' : cfel_data["MLE"]["content"]["0"][s-1][10],
                'human' : cfel_data["MLE"]["content"]["0"][s-1][2],
                'avian' : cfel_data["MLE"]["content"]["0"][s-1][4],
                'mammalian' : cfel_data["MLE"]["content"]["0"][s-1][5],
                'between' : cfel_data["MLE"]["content"]["0"][s-1][3],
                'HZ' : cfel_data["MLE"]["content"]["0"][s-1][12],
                'HA' : cfel_data["MLE"]["content"]["0"][s-1][13],
                'HM' : cfel_data["MLE"]["content"]["0"][s-1][14],
                'ZA' : cfel_data["MLE"]["content"]["0"][s-1][15],
                'ZM' : cfel_data["MLE"]["content"]["0"][s-1][16],
                'AM' : cfel_data["MLE"]["content"]["0"][s-1][17]
            }
            
    return results
            
           
    
def load_fade (file, sites, results, tag):
    
    with open (file, 'r') as fh:
        fade_data = json.load (fh)
        for s in sites:
            if not s in results:
                results[s] = {}
            if not tag in results[s]:
                results[s][tag] = {}
            results[s][tag]['fade'] = {}
            
            for res, rows in fade_data[ "MLE"]["content"].items():
                if rows["0"][s-1][3] >= 10:
                    results[s][tag]['fade'][res] = rows["0"][s-1][3]
        
    return results

def print_ratio (n,d):
    if d == 0.:
        if n == 0.: return "N/A"
        return "inf"
    n = n/d
    if n > 1000.:
        return ">1000"
    return "%.3g" % n

def print_sig (p):
    if p > 0.05: return ''
    if p > 0.01: return '*'
    if p > 0.001: return '**'
    return '***'
    
def print_fade (f):
    if len (f):
        return " ".join (["%s (%.2g)" % (k,v) for k, v in f.items()])
    return ""  
    
def print_cfel (f):
    res = []
    if f['HA'] < 0.05: res.append ("H%sA %s" % ('<' if f["avian"] > f["human"] else '>', print_sig (f['HA'])))
    if f['HM'] < 0.05: res.append ("H%sM %s" % ('<' if f["mammalian"] > f["human"] else '>', print_sig (f['HM'])))
    if f['AM'] < 0.05: res.append ("A%sM %s" % ('<' if f["mammalian"] > f["avian"] else '>', print_sig (f['AM'])))
    if f['overall'] < 0.05: res.append ("overall %s" % print_sig(f['overall']))
    return ", ".join (res)
    
results = {}
site_list = [int (k) for k in import_settings.site]

load_slac (make_path('SLAC-human.json'), site_list, results, "Human")
load_slac (make_path('SLAC-avian.json'), site_list, results, "Avian")
load_slac (make_path('SLAC-mammals.json'), site_list, results, "Mammals")
load_slac (make_path('SLAC-zoonotic.json'), site_list, results, "Zoonotic")

load_fel (make_path('FEL-human.json'), site_list, results, "Human")
load_fel (make_path('FEL-avian.json'), site_list, results, "Avian")
load_fel (make_path('FEL-mammalian.json'), site_list, results, "Mammals")
load_fel (make_path('FEL-zoonotic.json'), site_list, results, "Zoonotic")

load_meme (make_path('MEME-human.json'), site_list, results, "Human")
load_meme (make_path('MEME-avian.json'), site_list, results, "Avian")
load_meme (make_path('MEME-mammalian.json'), site_list, results, "Mammals")
load_meme (make_path('MEME-zoonotic.json'), site_list, results, "Zoonotic")

load_fade (make_path('HA.prot.fas.FADE.json'), site_list, results, "Human")

load_cfel (make_path('CFEL.json'), site_list, results)


for s,d in results.items():
    #print (s, d['cfel'])
    #continue
    print ("\t".join ([str (k) for k in [
        s-16 - (1 if s < 60 else 0),
        "%d:%d" % (d['Avian']['slac']['info'][2] ,d['Avian']['slac']['info'][3]),
        "%d:%d" % (d['Mammals']['slac']['info'][2] ,d['Mammals']['slac']['info'][3]),
        "%d:%d" % (d['Human']['slac']['info'][2] ,d['Human']['slac']['info'][3]),
        "%s %s" % (print_ratio (d['Avian']['fel'][1], d['Avian']['fel'][0]), print_sig(d['Avian']['fel'][4])),
        "%s %s" % (print_ratio (d['Mammals']['fel'][1], d['Mammals']['fel'][0]), print_sig(d['Mammals']['fel'][4])),
        "%s %s" % (print_ratio (d['Human']['fel'][1], d['Human']['fel'][0]), print_sig(d['Human']['fel'][4])),
        "%d %s" % (d['Avian']['meme'][7], print_sig(d['Avian']['meme'][6])),
        "%d %s" % (d['Mammals']['meme'][7], print_sig(d['Mammals']['meme'][6])),
        "%d %s" % (d['Human']['meme'][7], print_sig(d['Human']['meme'][6])),
        "%s" % print_fade (d['Human']['fade']),
        "%s" % print_cfel (d['cfel'])
   ]]))

