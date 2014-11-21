# -*- coding: utf-8 -*-
"""
Created on Thu Nov 20 16:14:02 2014

@author: dima
"""

import subprocess

def runcmd(x):
   assert (type(x) == str), 'not a string!'
   p = subprocess.Popen( x , shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
   out, err = p.communicate()
   return (out.decode('utf-8'), err.decode('utf-8'))