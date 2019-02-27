#!/usr/bin/env python

import random
from collections import defaultdict as dd
from collections import Counter as cc
import sys
import os
from math import log 
from random import random



class dot_log:
	#def __init__(self,VERBOSE,msg=' New_Module ',PARENT=False):
	def __init__(self,args,prefix=None,command_line=None):

		if args.verbose:	self.active = True 
		else:			self.active = False
		self.prefix = prefix

		if self.prefix:  self.log_file = 	 self.prefix+'.logfile'
		else: 		 self.log_file = 'progress_log.logfile'

		self.fout = open(self.log_file,'w') 
		self.out = sys.stderr
#		self.active = True

		self.intro = 'Module: '
		self.log = [[]] 

		self.out_write("\nRAGE COMMAND: "+"".join(" ".join(command_line).split("\'"))+'\n')


		if args.test: self.out_write('Rage Begins (test run): '+'\n\n')
		else: 	      self.out_write('Rage Begins: '+'\n\n')

#		if self.active: 
#			if args.test: self.out.write('Rage Begins (test run): '+'\n\n')
#			else: 	      self.out.write('Rage Begins: '+'\n\n')

		self.blank = ' '
		self.topic, self.subtopic = None, None 
		self.topics = [] 
		self.events = {} 
		self.notes = dd(list) 
		self.status = None 


	def reset(self): 
		if self.prefix:  self.log_file = 	 self.prefix+'.logfile'
		else: 		 self.log_file = 'progress_log.logfile'

		self.fout = open(self.log_file,'w') 
		self.out = sys.stderr


	def out_write(self,out_str):
	
		if self.active: 
			sys.stderr.write(out_str) 

		self.fout.write(out_str) 


	def start_major(self,topic,msg=None,subtopics=[]): 


		if self.status == 'minor': self.end() 

		self.status = 'major'


		self.out_write(self.intro+topic+'\n') 
		self.blank = ''.join([' ' for x in range(len(self.intro))])

		if msg: 
			self.out_write(self.blank+msg+' '+', '.join([str(s) for s in subtopics])+'\n') 
			self.sub_blank = ''.join([' ' for x in range(len(self.intro+msg))])
		else:
			self.sub_blank = self.blank	








		self.log.append([topic])
		self.topic, self.subtopic = topic, None 
		self.topics.append(topic) 
		self.events[self.topic] = dd(list) 
		self.topic_counter = 0 
		


	def start_minor(self,topic,block_len=10,ENUMERATE=False):

		if self.status == 'minor': self.end() 

		self.status = 'minor'
		self.subtopic = topic
		self.log[-1].append(topic)
		self.counter = 0 
		self.ENUMERATE=ENUMERATE
		self.block_len,self.mp,self.dots,self.numbers = block_len,block_len,0,0 
		if block_len > 100: self.mp = 100 

		
		#if self.active: self.out.write(self.sub_blank+topic+'...')
		#if self.active: 
		self.out_write(self.sub_blank+topic+'...')



	def mark(self,dotrate=0):
		self.counter +=1 

		if dotrate >= 1:  self.out_write('.') 
		elif dotrate > 0: 
			if random() < dotrate: self.out_write('.') 
		else:

			if self.counter % self.mp == 0: 
				self.out_write('.') 
				self.dots +=1
				if self.dots == 100:  self.mp = self.mp * self.mp 
				elif self.dots == 50: self.mp *= 4
				elif self.dots == 25: self.mp *= 2  
				elif self.dots == 10: self.mp += self.mp 

				if self.dots > 10 and self.dots % 10 == 0: self.mp += self.dots 

			

			if self.counter % self.block_len == 0 and self.ENUMERATE: 

				self.out_write('.'+str(self.counter)+'.')
				self.numbers += 1 
				if self.numbers == 4: self.block_len *= 2 
				elif self.numbers == 8: self.block_len *= 2 
				elif self.numbers % 10 == 0: self.block_len *= 2 
	










				
	def notate(self,note):
		self.out_write(self.sub_blank+note)
		self.notes[self.topic].append(note) 	
					


	def end(self): 

		self.out_write('Complete\n') 
			
		self.events[self.topic][self.subtopic] = [self.counter]
		self.subtopic = None 

