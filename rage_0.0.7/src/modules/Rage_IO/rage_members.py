import sys
from collections import defaultdict as dd
from collections import Counter as cc
from collections import MutableSequence
import numpy as np
from math import log
import mmap 
from sklearn.preprocessing import MinMaxScaler
from random import shuffle
import matplotlib.pyplot as plt
from matplotlib import cm as cm
from matplotlib.lines import Line2D as Line
from matplotlib.patches import Circle, Wedge, Polygon
import copy 





COLORS_1 = [ 'indigo', 'gold', 'hotpink', 'firebrick', 'indianred', 'sage', 'yellow', 'mistyrose', 'darkolivegreen', 'olive', 'darkseagreen', 'pink', 'tomato', 'lightcoral', 'orangered', 'navajowhite', 'lime', 'palegreen', 'darkslategrey', 'greenyellow', 'burlywood', 'seashell', 'mediumspringgreen', 'fuchsia', 'papayawhip', 'blanchedalmond', 'chartreuse', 'dimgray', 'black', 'peachpuff', 'springgreen', 'aquamarine', 'white', 'orange', 'lightsalmon', 'darkslategray', 'brown', 'ivory', 'dodgerblue', 'peru', 'darkgrey', 'lawngreen', 'chocolate', 'crimson', 'forestgreen', 'slateblue', 'lightseagreen', 'cyan', 'mintcream', 'silver', 'antiquewhite']

COLORS_2 = [ 'mediumorchid', 'skyblue', 'gray', 'darkturquoise', 'goldenrod', 'darkgreen', 'floralwhite', 'darkviolet', 'darkgray', 'moccasin', 'saddlebrown', 'grey', 'darkslateblue', 'lightskyblue', 'lightpink', 'mediumvioletred', 'slategrey', 'red', 'deeppink', 'limegreen', 'darkmagenta', 'palegoldenrod', 'plum', 'turquoise', 'lightgrey', 'lightgoldenrodyellow', 'darkgoldenrod', 'lavender', 'maroon', 'yellowgreen', 'sandybrown', 'thistle', 'violet', 'navy', 'magenta', 'dimgrey', 'tan', 'rosybrown', 'olivedrab', 'blue', 'lightblue', 'ghostwhite', 'honeydew', 'cornflowerblue', 'linen', 'darkblue', 'powderblue', 'seagreen', 'darkkhaki']

COLORS_3 = [ 'snow', 'sienna', 'mediumblue', 'royalblue', 'lightcyan', 'green', 'mediumpurple', 'midnightblue', 'cornsilk', 'paleturquoise', 'bisque', 'slategray', 'darkcyan', 'khaki', 'wheat', 'teal', 'darkorchid', 'deepskyblue', 'salmon', 'darkred', 'steelblue', 'palevioletred', 'lightslategray', 'aliceblue', 'lightslategrey', 'lightgreen', 'orchid', 'gainsboro', 'mediumseagreen', 'lightgray', 'mediumturquoise', 'darksage', 'lemonchiffon', 'cadetblue', 'lightyellow', 'lavenderblush', 'coral', 'purple', 'aqua', 'lightsage', 'whitesmoke', 'mediumslateblue', 'darkorange', 'mediumaquamarine', 'darksalmon', 'beige', 'blueviolet', 'azure', 'lightsteelblue', 'oldlace']

COMMON = ['red','blue','green','yellow','orange','purple','lime','cyan','k']
COMMON  = ['blue','green', 'red', 'orange','purple']+['lime','cyan',  'magenta','yellow']+['gold','peru',  'salmon',  'olive',    'khaki','indigo','crimson']
ALL_COLORS = COMMON+list(set([x for x in [c for c in COLORS_1+COLORS_2+COLORS_3] if x not in COMMON]))








def get_color_list(list_size,MISSING,rank=None,OFFSET=0):


	

	basics  = ['blue','green', 'red', 'orange','purple']
	brights = ['lime','cyan',  'magenta','yellow']
	pastels = ['gold','peru',  'salmon',  'olive',    'khaki','indigo','crimson']





	if MISSING: return [c for c in ALL_COLORS[OFFSET:OFFSET+list_size] if c != 'black'],'black'
	else: 	    return ALL_COLORS[OFFSET:OFFSET+list_size],None





def get_marker_list(list_size,MISSING,rank=None): 


	marks = ["o","v","^","<",">","8","s","p","H","D","d"]+["o","v","^","<",">","8","s","p","H","D","d"]+["o","v","^","<",">","8","s","p","H","D","d"]+["o","v","^","<",">","8","s","p","H","D","d"]
	
	if MISSING: return [m for m in marks if m != "s"][0:list_size],"s"
	else:       return marks[0:list_size], None 















def member_warning(msg):
        sys.stderr.write('RageMembersWarning: '+msg+'\n')



def member_error(msg):
	sys.stderr.write('\n')
        sys.stderr.write('RageMembersError: '+msg+'\n')
        sys.exit()


def command_line_error(msg):
	sys.stderr.write('\n')
        sys.stderr.write('RageCommandLineError: '+msg+'\n')
        sys.exit()


class cnt_key(dd):
    	def __missing__(self, key):
        	return 0.0

	def update(self,mylist):
		for (a,b) in mylist: self[a] = b 
		return self


class attribute_key(dd):
    	def __missing__(self, key):
        #value = self.default_factory(key)
        	self[key] = 'NA'
        	return 'NA'


MIN_GRP_SIZE=20





		
		
			
class Member:
	def __init__(self,name,idx):
		self.name = name 
		self.idx  = idx 
		self.cnts = cnt_key(None) 
		self.attributes = attribute_key(None)
		self.preds =  dd(float)  
		self.regVariables = dd(list) 
		self.label = None 
		self.notes = dd(bool)
		self.null = False  
		self.hide = False
		self.norms = dd(lambda: cnt_key(None)) 

 
	def update(self):
		self.cnt_total = sum(self.cnts.values()) 	
		self.len = len(self.cnts) 
		return self.name 


	def add_line_cnts(self,cnts):
		for idx in [i for i in range(len(cnts)) if cnts[i] > 0]: self.cnts[idx] = cnts[idx]  
		self.cnt_total = sum(self.cnts.values()) 	
		self.len = len(self.cnts) 
		return self






class Members(MutableSequence):

	"""A container for manipulating lists of hosts"""
	def __init__(self, data=None):
		"""Initialize the class"""
		super(Members, self).__init__()
		if (data is not None):
			self._list = [Member(data[i],i) for i in range(len(data))]
			self.len = len(self._list) 
			self.removed, self.attributes,self.attribute_type,self.attribute_class = [], [], [] , {} 
		else:
	    		self._list = list()
			self.len = 0 
			self.removed, self.attributes,self.attribute_type,self.attribute_class = [], [], [] , {} 

		self.keys = {} 
		self.notes = {} 
		self.color_offset = 0 
		self.min_group_size = 5 
		self.max_group_members = 10 


	def collate(self,label): 
		self.label = label
		self.names = [x.update() for x in self._list]
		self.lookup = {member.name : member.idx for member in self._list}
		self.len = len(self._list)
#		print len(self._list)
#		print len(self.__list) 
#		self.__list = self.__list[0:20] 

	def add_attribute(self,opt,optType):
		self.attributes.append(opt) 
		self.attribute_type.append(optType) 
		self.attribute_class[opt] = optType

#	def get_attributes(self,attributes):

#		self.key  = {} 
#		s_key = dd(lambda: {}) 
#		for a in attributes: 
#			s_key[a] = [s.attributes[a] for i,s in enumerate(self._list)]
#		return s_key


	def __repr__(self):
		return "<{0} {1}>".format(self.__class__.__name__, self._list)

	def __len__(self):
		"""List length"""
		return len(self._list)

	def __getitem__(self, ii):
		"""Get a list item"""
		return self._list[ii]

	def __delitem__(self, ii):
		"""Delete an item"""
		del self._list[ii]

	def __setitem__(self, ii, val):
		# optional: self._acl_check(val)
		return self._list[ii]

	def __str__(self):
		return str(self._list)

	def insert(self, ii, val):
		# optional: self._acl_check(val)
		self._list.insert(ii, val)

	def append(self, member):
		self.insert(len(self._list), member)

			
	def segregate(self,attribute,preds=[]):

		

		if attribute in self.attribute_class and self.attribute_class[attribute] != 'binary': 
			slist= sorted([(s.attributes[attribute],i) for i,s in enumerate(self._list)])
			slist= [g for g in slist if g[0] != 'NA']
			gIdx = len(slist)/2 
			g1 = slist[0:gIdx]


			while slist[gIdx][0] == g1[-1][0]: 
				g1.append(slist[gIdx])
				gIdx+=1 
				if gIdx == len(slist):
					seg =  {'HI': [g[1] for g in slist if g[0] !='NA'], 'LO': [g[1] for g in slist if g[0] == 'NA']} 
					return seg,{sa: len(seg[sa]) for sa in seg.keys()}

			g2 = slist[gIdx::]  
			seg =  {'HI': [g[1] for g in g2], 'LO': [g[1] for g in g1]} 

		else:
			seg = dd(list) 
			self.attribute_class[attribute] = 'binary'
			for i,s in enumerate(self._list):
				seg[s.attributes[attribute]].append(i) 
		return seg,{sa: len(seg[sa]) for sa in seg.keys()}


	def prune_cnts(self,min_val):

		self._list = [f for f in self._list if len(f.cnts) > min_val]
		self.len = len(self._list) 
		return self	
 


	def transpose(self): 
		self._list = [[self._list[i][j] for i in range(len(self._list))] for j in range(len(self._list[0]))]
		self.n_len,self.m_len,self.n_rng,self.m_rng = self.m_len,self.n_len,self.m_rng,self.n_rng
		return self

	def transform(self,transform):
		hide = False
		if   transform == 'log':	self._list = [[log(self._list[i][j]+1,2) for j in range(self.n_len)] for i in range(self.m_len)]
		elif transform == 'hide':  hide = True
		else:
			print transform,'wtf'
			sys.exit()
		return hide,self

	
	def center(self):
		self._list = [[self._list[i][j] - sum(self._list[i])/float(self.n_len) for j in range(self.n_len)] for i in range(self.m_len)]
		return self


	def filter(self,attribute_list,attribute_prune_list=[]):



		unknown_attributes = [p.split('=')[0] for p in attribute_prune_list if p.split('=')[0].split('!')[0] not in self.attributes]+[a.split('=')[0] for a in attribute_list if a.split('=')[0] not in self.attributes]
		if len(unknown_attributes)>0: 	member_error('Unknown Predictors: '+",".join(unknown_attributes))



		self.removed = [] 
		for a in attribute_list:


			if self.attribute_class[a.split('=')[0]] != 'binary': 	
				self.removed.extend([s.idx for s in self._list if s.attributes[a] == 'NA'])	
			else:


				if len(a.split('=')) > 1: 
					self.removed.extend([s.idx for s in self._list if s.attributes[a.split('=')[0]].split('~')[-1] not in a.split('=')[1].split(',')])
				else:
					self.removed.extend([s.idx for s in self._list if s.attributes[a].split('~')[-1] == 'NA' or 'NA' in s.attributes[a].split('~')[-1].split('-')])	
					
			continue 


			if len(self._list[0].name.split('~'))>2 and self._list[0].name.split("~")[-2] == a: 
				for s in self._list:
					s.attributes[a] = s.name.split('~')[-1] 
			else:
				
				if self.attribute_class[a.split('=')[0]] != 'binary': 	
					self.removed.extend([s.idx for s in self._list if s.attributes[a] == 'NA'])	
				else:


					if len(a.split('=')) > 1: 
						self.removed.extend([s.idx for s in self._list if s.attributes[a.split('=')[0]].split('~')[-1] not in a.split('=')[1].split(',')])
					else:
						self.removed.extend([s.idx for s in self._list if s.attributes[a].split('~')[-1] == 'NA' or 'NA' in s.attributes[a].split('~')[-1].split('-')])	
						

















		for p in attribute_prune_list:
			pA,pR,pX = p.split('=')[0], p.split('=')[0].split('!')[0],p.split('=')[-1].split(',') 
			if pR != pA: 	self.removed.extend([s.idx for s in self._list if s.attributes[pR].split('~')[-1] in pX])
			else:		self.removed.extend([s.idx for s in self._list if s.attributes[pA].split('~')[-1] not in pX])

		self.removed = sorted(list(set(self.removed)))
		
		if len(self._list) - len(self.removed) < 2: 
			member_error('Incongrous Covariates/Prune Requests: All samples removed') 



		self.swap = {}  
		self._list = [s for s in self._list if s.idx not in self.removed]
		for i in range(len(self._list)):	
			self.swap[self._list[i].idx] = i 
			self._list[i].idx = i 
			
		self.len = len(self._list) 
		return self
	



	def colorize_names(self,key):


		for i,s in enumerate(self._list):


			
			self._list[i].label 


			if s.name.split('=')[-1] in key: 
			
				self._list[i].label = key[s.name.split('=')[-1]]

				#s.notes['color'] = key[s.name.split('=')[-1]]
			elif s.name.split('~')[-1] in key:
				
				#s.notes['color'] = key[s.name.split('~')[-1]]
				self._list[i].label = key[s.name.split('~')[-1]]

			elif "~".join(s.name.split('~')[1::]) in key:
				self._list[i].label = key["~".join(s.name.split('~')[1::])]

			else: 
				s.notes['color'] = 'silver'









	
	def set_regression_labels(self,args):
		
		colors, missed, marker, size  = [], [] , None, None 




 
		p_bin   = [p for p  in [a.split('=')[0] for b in [P.split('*') for P in args.predictors] for a in b] if self.attribute_class[p] == 'binary']
		p_cont  = [p for p  in [a.split('=')[0] for b in [P.split('*') for P in args.predictors] for a in b] if self.attribute_class[p] != 'binary']
		c_bin   = [p for p  in [a.split('=')[0] for b in [P.split('*') for P in args.covariates] for a in b] if self.attribute_class[p] == 'binary']
		c_cont  = [p for p  in [a.split('=')[0] for b in [P.split('*') for P in args.covariates] for a in b] if self.attribute_class[p] != 'binary']
		variables = p_bin + p_cont + c_bin + c_cont 


		

		for r in variables: 
			if (r in colors) or (r == marker) or (r == size): continue 
			if len(colors) == 0: colors.append(r) 

			elif self.attribute_class[r] == 'binary': 



				if marker == None:  marker = r 
				elif len(colors) < 3: colors.append(r) 
				else: 		     missed.append(r) 
			else: 
				if size == None:     	size = r 
				elif len(colors) < 3:  	colors.append(r)       

		return colors, marker, size





	def create_plot_labels(self,args,SIZE=0):



		if args.mingroupsize > 0: self.min_group_size = args.mingroupsize 
		else: 			  self.min_group_size = 10 

		if args.maxgroups > 0:    self.max_group_members = args.maxgroups
		else: 			  self.max_group_members = 10 


		self.notes['colors'], label_notes, self.colorize_key, binary_colors, continuous_colors, binary_markers = dd(list), {}, {}, {},{}, {} 
		colors, missed, marker, size, left_key, right_key, edge_key, size_labels, marker_labels = [], [], None, None, None,None,None, None, None 


		if 'color' in args and len(args.color) > 0:
			colors, marker, size = args.color, args.marker, args.size
		elif args.option == 'regression':		
			colors, marker, size = self.set_regression_labels(args) 
		else:
			colors = [self.attributes[0]]


		for i,a in enumerate(colors):
			if self.attribute_class[a] == 'binary': 	
				binary_colors[a] = self.binary_labels(a,'color',rank=(i+1,len(colors))) 
			else:						
				continuous_colors[a] = self.continuous_labels(a,'color')

		if len(binary_colors) > 0: 
			left_key, right_key  = binary_colors.keys()[0],  binary_colors.keys()[-1] 
			if len(continuous_colors)> 0: 	edge_key  = continuous_colors.keys()[0]
		elif len(continuous_colors) > 0:
			left_key, right_key  = continuous_colors.keys()[0],  continuous_colors.keys()[-1] 

		if marker != None: 	binary_markers = self.binary_labels(marker,'marker') 

		if SIZE != 0: size_labels = [SIZE for s in self._list]
		elif size != None: 
			if self.attribute_class[size] == 'binary': 
				print 'binary size'
				sys.exit() 
			else:
				s_mean = np.mean([s.attributes[size] for s in self._list if s.attributes[size] != 'NA']) 
				s_data = np.array([s.attributes[size] if s.attributes[size] != 'NA' else s_mean for s in self._list]).reshape(-1,1)
				scaler = MinMaxScaler().fit_transform(s_data)
				size_labels = [(scaler[j][0]*7)+3 for j in range(len(scaler))]
		if left_key: 
			if left_key == right_key: self.notes['labels'] = {'color': left_key,'shape': marker, 'size': size} 
			else: 			  self.notes['labels'] = {'left_color': left_key, 'right_color': right_key, 'shape': marker, 'size': size}



		for i,s in enumerate(self._list):
			p = PlotLabel()

			if left_key:				
				try: p.set_left_color(binary_colors[left_key][1][i], binary_colors[left_key][0][i])
				except KeyError: p.set_left_color(continuous_colors[left_key][1][i], continuous_colors[left_key][0][i])
			if right_key and right_key != left_key:	
				try: p.set_right_color(binary_colors[right_key][1][i], binary_colors[right_key][0][i])
				except KeyError: p.set_right_color(continuous_colors[right_key][1][i], continuous_colors[right_key][0][i])
			if edge_key:				p.set_edge_color(continuous_colors[edge_key][1][i], continuous_colors[edge_key][0][i]) #  edge_labels[i], edge_colors[i])
			if size_labels: 			p.set_size_label(size_labels[i]) 
			if marker: 				p.set_marker_label(binary_markers[1][i], binary_markers[0][i]) #  marker_labels[i],markers[i]) 
			self._list[i].label = p.check() 
			self.colorize_key[p.id] =  self._list[i].label  
			





	def continuous_labels(self,attribute,label='color'):

		sIds = [s.attributes[attribute] for s in self._list if s.attributes[attribute] != 'NA']
		sAvg = np.mean(sIds)
		sImpute = [s.attributes[attribute] if s.attributes[attribute] != 'NA' else sAvg for s in self._list]
		if label == 'color':
			norm = plt.Normalize(min(sIds),max(sIds))
			return [plt.cm.hot(norm(x)) for x in sImpute], [attribute for x in sImpute] 








	def binary_labels(self,attribute,label='color',rank=None):

		
		sAll,sValid = [s.attributes[attribute] for s in self._list], [s.attributes[attribute] for s in self._list if s.attributes[attribute] != 'NA']


			


		if len(set(sValid)) > 0:	
			try: minSize = min(self.min_group_size,sorted(cc(sValid).values())[-2])
			except IndexError: minSize = 10
			sGrps, sIdx, sMissing = list(set(sValid)) , [], [] 
			if len(sGrps) > 2:  sGrps = sorted([a for (a,b) in cc(sValid).items() if b > minSize],reverse=True)
			if len(sGrps) > self.max_group_members:  sGrps = sGrps[0:self.max_group_members]



			for s in self._list:
				if s.attributes[attribute] in sGrps:	
					sIdx.append(sGrps.index(s.attributes[attribute]))
				else:
					sIdx.append(len(sGrps))
					sMissing.append(s.attributes[attribute])
			if label == 'color': 
				label_list, missing_label = get_color_list(len(sGrps),len(sMissing)!=0,rank,OFFSET=self.color_offset)
				self.color_offset+=len(label_list) 

			elif label == 'marker':

 
				label_list, missing_label = get_marker_list(len(sGrps),len(sMissing)!=0,rank)



			if len(sMissing) > 0:
				sGrps.append(attribute+'='+",".join([sM.split('~')[-1] for sM in list(set(sMissing))]))
				label_list.append(missing_label) 




		
			label_vals, label_labels = [label_list[i] for i in sIdx], [sGrps[i] for i in sIdx]
			return label_vals, label_labels 
				











class PlotLabel:
	def __init__(self):


		self.left_color  = None
		self.right_color = None
		self.edge_color  = None
		self.marker      = None 
		self.size        = 10
		self.label_id    = [] 
		self.labels = [] 

	def set_size_label(self,size_label): 
		self.size = size_label

	def set_marker_label(self,label,marker): 

		self.marker = marker 
		self.marker_label = label
		self.label_id.append([label,marker]) 

	def set_left_color(self,label,color):

		self.left_color = color
		self.left_color_label = label 
		self.label_id.append([label,color]) 

	def set_right_color(self,label,color):

		self.right_color = color
		self.right_color_label = label 
		self.label_id.append([label,color]) 
		
	def set_edge_color(self,label,color):
		self.edge_color = color
		self.edge_color_label = label
		self.label_id.append([label,color]) 

	def update_size(self,size):
		self.vals['markersize'] = size


	def check(self,msize=10):


		self.id = '|'.join([L[0] for L in self.label_id]) 


		if self.left_color and self.right_color:	self.vals = {'color': self.left_color, 'fillstyle': 'left', 'marker': 'o', 'markerfacecoloralt': self.right_color, 'mew': 0.0}
		elif self.left_color:				self.vals = {'color': self.left_color, 'marker': 'o', 'mew': 0.0} #, 'markersize': self.size}
		else: 						self.vals = {'mew': 0.0} 
	

		if self.marker: self.vals['marker'] = self.marker 		
		self.proxy =  Line([0],[0],linestyle='none',markersize=msize,**self.vals)

		self.vals['markersize'] = self.size 


		return self























































	



























































	
