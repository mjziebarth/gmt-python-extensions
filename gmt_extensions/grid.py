# A grid class.
#from scipy.interpolate import SmoothSphereBivariateSpline
from scipy.interpolate import LinearNDInterpolator,NearestNDInterpolator
from tempfile import NamedTemporaryFile
import numpy as np
import gmt
from os import remove

def _is_gridded(lon, lat):
	"""
	This method checks whether the coordinates given by lon
	and lat are already gridded.
	"""
	# Quick and dirty:
	nontrivial_dim_lon = np.count_nonzero([x > 1 for x in lon.shape])
	nontrivial_dim_lat = np.count_nonzero([x > 1 for x in lat.shape])
	if nontrivial_dim_lon == 1 and nontrivial_dim_lat == 1:
		return False
	elif nontrivial_dim_lon == 2 and nontrivial_dim_lat == 2:
		# Some more checks would be suited.
		if len(lon.shape) != len(lat.shape):
			raise TypeError("Error: lon and lat are not compatible in shape!")
		try:
			sameshape = all(lon.shape == lat.shape)
		except:
			sameshape = lon.shape == lat.shape
		if sameshape:
			return True
	print("success!")
	raise TypeError("Error: lon and lat are not compatible in shape!")
	

class Grid:
	
	@staticmethod
	def grid_data(lon, lat, val, gridlon, gridlat, rect, filename=None,
	              interpolator='nearest', maskfun=None, **kwargs):
		
		# Check input data:
		if not isinstance(gridlon,np.ndarray):
			gridlon = np.array(gridlon)
		if not isinstance(gridlat,np.ndarray):
			gridlat = np.array(gridlat)
		
		# Create grid, if gridlon and gridlat are not yet in
		# gridded shape:
		if not _is_gridded(gridlon, gridlat):
			gridlon, gridlat = np.meshgrid(gridlon, gridlat,
			                               indexing='ij')
		
		# Grid the data.
		# Use linear interpolation in 3D. Should be close enough to
		# real spherical interpolation at short distances.
		D2R = np.pi/180.0
		vec_dat =  np.array([
		              np.cos(D2R*lon)*np.cos(D2R*lat),
		              np.sin(D2R*lon)*np.cos(D2R*lat),
		              np.sin(D2R*lat)
		           ]).T
		vec_grid = np.array([
		              np.cos(D2R*gridlon)*np.cos(D2R*gridlat),
		              np.sin(D2R*gridlon)*np.cos(D2R*gridlat),
		              np.sin(D2R*gridlat)
		           ]).T
		if interpolator == 'nearest':
			interpolator = NearestNDInterpolator(vec_dat, val)
		elif interpolator == 'linear':
			# Not linear in spherical coordinates but close at
			# small distances.
			interpolator = LinearNDInterpolator(vec_dat, val)
		else:
			raise ValueError("Interpolator '" + str(interpolator) + "' unknown!")
		grd = interpolator(vec_grid)
		
		# Mask values:
		if maskfun is not None:
			try:
				# Try fast indexing first:
				grd[maskfun(grd)] = np.NaN
			except:
				# Iterate over elements:
				for g in grd:
					if maskfun(g):
						g = np.nan
		
		# Return grid:
		return Grid(gridlon, gridlat, grd, rect, filename=filename, **kwargs)
	
	def __init__(self, lon, lat, val, rect, filename=None):
		"""
		Parameters:
		
		- lon
		- lat
		- val
		
		Optional:
		
		- filename (default: None)
		  (Persistent) Cache file name
		"""
		self._lon = lon
		self._lat = lat
		self._val = val
		self._filename = filename
		self._remove_file = False
		self._has_grid = False
		self._rect = rect
		if not _is_gridded(lon, lat) and \
		   not (np.issorted(lon) and np.issorted(lat)):
		   # Q&D
		   raise ValueError("For non-gridded data, lon and lat "
		                    "have to be sorted!")
	
	def __del__(self):
		if hasattr(self,'_remove_file') and self._remove_file:
			remove(self._filename)
		
	
	def _create_grd_file(self):
		
		# Step 1: Create temporary csv:
		if not _is_gridded(self._lon, self._lat):
			lon, lat = np.meshgrid(self._lon, self._lat,
			                       indexing='ij')
			dlon = lon[1]-lon[0]
			dlat = lat[1]-lat[0]
		else:
			lon, lat = self._lon, self._lat
			lon_unq = np.sort(np.unique(lon))
			lat_unq = np.sort(np.unique(lat))
			dlon = lon_unq[1]-lon_unq[0]
			dlat = lat_unq[1]-lat_unq[0]
		
		csv_path = None
		try:
			# Create temporary CSV:
			nanvalue = None
			with NamedTemporaryFile(mode='w',suffix='.csv', delete=False) as fcsv:
				if np.any(np.isnan(self._val)):
					# Handle NaN:
					nanvalue = int(np.floor(self._val[~np.isnan(self._val)].min()-1))
					val = self._val.flatten()
					val[np.isnan(val)] = nanvalue
				else:
					val = self._val.flatten()
				csv_path = fcsv.name
				np.savetxt(fcsv, val, delimiter=",")
			
			
			# Step 2: Call xyz2grd:
			with gmt.clib.Session() as lib:
				# Make sure we have a filename:
				if self._filename is None:
					# Temporary file:
					with NamedTemporaryFile(suffix='.grd',delete=False) as f:
						self._filename=f.name
						self._remove_file = True
					# Now we hope that between the previous and the
					# next command nobody takes hold of the tempfile name.
				
				
				# Create the argument string:
				args = "%s -G%s -ZBL -I%.5f/%.5f -V -R" \
					   % (fcsv.name, self._filename, dlon, dlat) \
					   +str(self._rect)
				if nanvalue is not None:
					# Write NaN values to a value outside range:
					args += " -di%i" % (nanvalue)
				
				
				lib.call_module('xyz2grd',args)
			
				# Should probably do a check beforehand:
				self._has_grid = True
		finally:
			# Make sure temporary CSV is removed:
			if csv_path is not None:
				try:
					os.remove(csv_path)
				except:
					pass
	
	def __str__(self):
		# Make sure we have a .grd file:
		if not self._has_grid:
			self._create_grd_file()
		
		# Return path to .grd file:
		return self._filename
