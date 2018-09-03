# Provides colormap classes.
from abc import ABCMeta, abstractmethod
from tempfile import NamedTemporaryFile
from os import remove
import numpy as np

try:
	from matplotlib.pyplot import get_cmap as mpl_get_cmap
	from matplotlib.colors import Normalize, Colormap
	from matplotlib.cm import ScalarMappable
	has_mpl = True
except:
	has_mpl = False

class Colormap(metaclass=ABCMeta):
	"""
	Base class for color maps.
	
	This class can be used to hold color map information and generate
	temporary color palette tables on the fly when needed.
	
	Although rather quick and dirty, this yields comfort somewhat
	similar to matplotlib color maps.
	"""
	def __init__(self, name, vmin=0, vmax=1):
		self._name = name
		self._cpt = None
		self._cpt_path = None
		self._vmin = vmin
		self._vmax = vmax
	
	def __del__(self):
		if self._cpt_path is not None:
			remove(self._cpt_path)
	
	def _create_cpt(self, vmin, vmax):
		# Write CPT to temporary file:
		self._cpt = NamedTemporaryFile(mode='w',suffix='.cpt',delete=False)
		self._cpt_path = self._cpt.name
		
		self._cpt.write("# Color table for color map '" + self._name + "'\n"
		                "# COLOR_MODEL = RGB\n"
		                "#----------------------------------------------------------\n")
		
		dv = vmax-vmin
		for i in range(self.cpt_table.shape[0]-1):
			x0 = vmin + dv*self.cpt_table[i,0]
			x1 = vmin + dv*self.cpt_table[i+1,0]
			self._cpt.write("%.5f\t%i/%i/%i\t%.5f\t%i/%i/%i\n" % 
			                (x0,255*self.cpt_table[i,1],
			                 255*self.cpt_table[i,2],255*self.cpt_table[i,3],
			                 x1,255*self.cpt_table[i+1,1],
			                 255*self.cpt_table[i+1,2],255*self.cpt_table[i+1,3]))
		
		self._cpt.write("#----------------------------------------------------------")
		self._cpt.close()
	
	def __str__(self):
		# Make sure we have created the temporary cpt file:
		if self._cpt is None:
			self._create_cpt(self._vmin,self._vmax)
		
		return self._cpt_path


class MPLWrapperColormap(Colormap):
	"""
	This class wraps matplotlib colormaps and creates color palette tables
	from them.
	"""
	def __init__(self, cmap, vmin=0, vmax=1):
		if not has_mpl:
			raise ImportError("Could not import all matplotlib dependencies!")
		
		if isinstance(cmap,str):
			Colormap.__init__(self,cmap,vmin,vmax)
			cmap = mpl_get_cmap(cmap)
		else:
			try:
				name = str(cmap.name)
			except:
				raise ValueError("The argument 'cmap' is of wrong type!")
			
			Colormap.__init__(self,name,vmin,vmax)
		
		# Obtain table to convert to .cpt:
		if 'colors' in cmap.__dict__.keys():
			if isinstance(cmap.colors,list):
				N = len(cmap.colors)
			elif isinstance(cmap.colors,np.ndarray):
				N = cmap.colors.shape[0]
			else:
				raise ValueError("Data type not understood: " + str(type(cmap.colors)))
			self.cpt_table = np.zeros((N,4))
			self.cpt_table[:,1:4] = cmap.colors
			self.cpt_table[:,0] = np.linspace(0,1,N)
		elif '_segmentdata' in cmap.__dict__.keys():
			# Obtain the segments for all three colors:
			segs = [np.array(cmap._segmentdata['red']),
			        np.array(cmap._segmentdata['green']),
			        np.array(cmap._segmentdata['blue'])]
			
			changes =   set([seg[0] for seg in segs[0]]) \
			          | set([seg[0] for seg in segs[1]]) \
			          | set([seg[0] for seg in segs[2]])
			
			changes = np.sort(np.array(list(changes)))[1:]
			
			colors = []
			colors += [[0, segs[0][0,1], segs[1][0,1], segs[2][0,1]]]
			for x in changes:
				c0 = [0,0,0]
				c1 = [0,0,0]
				for k in range(3):
					id_ = np.argwhere(segs[k][:,0] < x)[-1]
					if id_ < len(segs[k]):
						if segs[k][id_+1,0] == x:
							# Red value changes at x:
							c0[k] = segs[k][id_+1,1][0]
							c1[k] = segs[k][id_+1,2][0]
						else:
							# Red value does not change at x. Interpolate linearly:
							x0 = segs[k][id_,0][0]
							v0 = segs[k][id_,2][0]
							x1 = segs[k][id_+1,0][0]
							v1 = segs[k][id_+1,1][0]
							c0[k] = (x-x0)/(x1-x0)*v0 + (x1-x)/(x1-x0)*v1
							c1[k] = c0[k]
				colors += [[x,c0[0],c0[1],c0[2]]]
				if c0[0] != c1[0] or c0[0] != c1[0] or c0[0] != c1[0]:
					colors += [[x,c1[0],c1[1],c1[2]]]
			
			self.cpt_table = np.array(colors)
