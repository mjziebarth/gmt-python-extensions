# A grid class.
#from scipy.interpolate import SmoothSphereBivariateSpline
from scipy.interpolate import LinearNDInterpolator,NearestNDInterpolator
from tempfile import NamedTemporaryFile
import numpy as np
import gmt
from os import remove
from shutil import copyfile

def same_shape(array1, array2):
	return np.array_equal(array1.shape, array2.shape)

def _mask_array(array, maskfun):
	"""
	Masks values of array where maskfun(array)==True.
	This is done by setting array[maskfun(array)] = NaN
	"""
	# Mask values:
	if maskfun is not None:
		try:
			# Try fast indexing first:
			array[maskfun(array)] = np.NaN
		except:
			# Iterate over elements:
			for g in array:
				if maskfun(g):
					g = np.nan


def _compatible_shape(lon, lat):
	# Quick and dirty:
	nontrivial_dim_lon = np.count_nonzero([x > 1 for x in lon.shape])
	nontrivial_dim_lat = np.count_nonzero([x > 1 for x in lat.shape])
	if nontrivial_dim_lon == 1 and nontrivial_dim_lat == 1:
		return True
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
	raise TypeError("Error: lon and lat are not compatible in shape!")
	

def _is_gridded(lon, lat, tolerance=1e-5):
	"""
	This method checks whether the coordinates given by lon
	and lat are already gridded.
	"""
	
	# Make sure that shape of both arrays is compatible:
	if not _compatible_shape(lon, lat):
		return False
	
	# Inference is hard if grid is small:
	if lon.size < 9:
		raise TypeError("Error: lon and lat arrays too small! Interpolate to "
		                "larger grid!")
	
	lon_ = lon.copy().flatten()
	lat_ = lat.copy().flatten()
	lon_[lon_<lon_[0]] += 360.0
	
	# Test the grid itself.
	# Find out which index varies faster:
	grid_property = dict()
	grid_property["fast_index"] = -1
	if lon_[0] != lon_[1]:
		# Longitude varies faster:
		if lat_[0] != lat_[1]:
			# May still be a grid but tilted:
			# TODO!
			dlat_ = lat_[1] - lat_[0]
			DLAT = np.round((lat_[1:] - lat_[:-1]) / dlat_).astype(int)
			
			grid_property["is_gridded"] = False
			return grid_property
		grid_property["fast_index"] = "LON"
		dl1 = lon_[1]-lon_[0]
		dl2 = lon_[2]-lon_[1]
		# We may just be at a wraparound.
		if abs(dl1) < abs(dl2):
			dlon = dl1
		else:
			dlon = dl2
		
	elif lat_[0] != lat_[1]:
		# Latitude varies faster:
		grid_property["fast_index"] = "LAT"
		
		# We assume that there is no wraparound in latitude. #TODO?
		dlat = lat_[1] - lat_[0]
	else:
		raise ValueError("Error: Longitude and latitude do not form a grid!")
	
	
	# Now make sure we're on a well-defined grid:
	if grid_property["fast_index"] == "LON":
		# Determine dimension of longitudes:
		M=0
		while M < lon_.size-1:
			# Check if next longitude coordinate is where expected:
			if abs((lon_[M+1] - lon_[M]) / dlon - 1) > 0.01:
				break
			M += 1
		M += 1
		# Determine dimension of latitudes:
		N=int(lon_.size/M)
		dlat = lat_[M] - lat_[0]
	elif grid_property["fast_index"] == "LAT":
		N=0
		while N < lat_.size-1:
			# Check if next latitude coordinate is expected:
			if abs((lat_[N+1] - lat_[N]) / dlat - 1) > 0.01:
				break
			N += 1
		N += 1
		# Determine dimension of longitudes:
		M=int(lat_.size/N)
		dl1 = lon_[N]-lon_[0]
		dl2 = lon_[2*N]-lon_[N]
		# We may just be at a wraparound.
		if abs(dl1) < abs(dl2):
			dlon = dl1
		else:
			dlon = dl2
	
	grid_property["Nlon"] = M
	grid_property["Nlat"] = N
	
	# Check grid:
	if grid_property["fast_index"] == "LON":
		lon_t, lat_t = np.meshgrid(dlon*np.arange(M)+lon_[0],
		                           dlat*np.arange(N)+lat_[0],
		                           indexing='xy')
	else:
		lon_t, lat_t = np.meshgrid(dlon*np.arange(M)+lon_[0],
		                           dlat*np.arange(N)+lat_[0],
		                           indexing='ij')
	
	lon_t = lon_t.flatten()
	lat_t = lat_t.flatten()
	
	diff_lon = np.max(np.abs(lon_t-lon_))
	diff_lat = np.max(np.abs(lat_t-lat_))
	diff = max(diff_lon,diff_lat)
	
	grid_property["is_gridded"] = diff < tolerance
	return grid_property
				

class Grid:
	"""
	A class for the automatic conversion of numpy arrays to gmt grids.
	
	
	::Class methods::
	
	Grid(lon, lat, val, rect, filename=None)
	
	   Create a grid instance.
	
	
	Static methods:
	
	grid_if_needed(lon, lat, val, gridlon, gridlat, filename=None,
	               interpolator='nearest', maskfun=None, verbosity=0)
	
	grid_data(lon, lat, val, gridlon, gridlat, rect, filename=None,
	          interpolator='nearest', maskfun=None)
	
	__str__():
	   Returns the path of a temporary .grd file which contains the grid.
	   The path can be passed to gmt. Can be heavy on file operations the
	   first time the .grd file is created.
	
	"""
	
	
	@staticmethod
	def grid_if_needed(lon, lat, val, gridlon, gridlat, rect,
	                   filename=None, interpolator='nearest',
	                   maskfun=None, verbosity=0,
	                   **kwargs):
		# TODO : Remove this method? Better to have option regrid=True
		#        or interpolate=True in constructor?
		
		# Check input data:
		if not isinstance(gridlon,np.ndarray):
			gridlon = np.array(gridlon)
		if not isinstance(gridlat,np.ndarray):
			gridlat = np.array(gridlat)
	
		# See whether we even have to grid lon and lat:
		grid_property = _is_gridded(lon,lat)
		if grid_property["is_gridded"]:
			gridlon = lon.copy()
			gridlat = lat.copy()
			grd = val.copy()
			_mask_array(grd, maskfun)
#			if not np.any(np.isnan(grd)):
#				print("MASK DID NOT CHANGE ANYTHING!")
#				print("should change:",np.count_nonzero(~maskfun(grd)))
			if verbosity > 1:
				print("No regrids!")
			return Grid(lon, lat, grd, rect, filename=filename, _no_grid_check=True,
			            fast_index=grid_property["fast_index"],
			            Nlon=grid_property["Nlon"], Nlat=grid_property["Nlat"], **kwargs)
		else:
			if verbosity > 0:
				print("Regridding...")
			return Grid.grid_data(lon, lat, val, gridlon, gridlat, rect,
				                  filename, interpolator, maskfun, **kwargs)
	
	
	@staticmethod
	def grid_data(lon, lat, val, gridlon, gridlat, rect, filename=None,
	              interpolator='nearest', maskfun=None,
	              **kwargs):
		
		# Check input data:
		if not isinstance(gridlon,np.ndarray):
			gridlon = np.array(gridlon)
		if not isinstance(gridlat,np.ndarray):
			gridlat = np.array(gridlat)
		
		# Create grid, if gridlon and gridlat are not yet in
		# gridded shape:
		grid_property = _is_gridded(gridlon, gridlat)
		if not grid_property["is_gridded"]:
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
		_mask_array(grd, maskfun)
		
		
		# Return grid:
		return Grid(gridlon, gridlat, grd, rect, filename=filename,
		            registration=registration, _no_grid_check=True,
		            fast_index=grid_property["fast_index"],
		            Nlon=grid_property["Nlon"], Nlat=grid_property["Nlat"], **kwargs)
	
	def __init__(self, lon, lat, val, rect, filename=None, registration='gridline',
	             fast_index=None, Nlon=None, Nlat=None, _no_grid_check=False):
		"""
		Parameters:
		
		- lon
		- lat
		- val
		
		Optional:
		
		- filename (default: None)
		  (Persistent) Cache file name
		"""
		# TODO: It would be best to only store lon, lat, and val in np.ndarrays.
		#       Do conversion in the creation routines. This reduces code base
		#       size and makes error-free code easier to write.
		if not isinstance(lon,np.ndarray) or not isinstance(lat,np.ndarray) \
		    or not isinstance(val,np.ndarray):
			raise TypeError("Error: lon, lat, and val have to be given as "
			                "numpy arrays!")
		if not same_shape(lon, lat) or not same_shape(lon,val):
			raise TypeError("Error: lon, lat, and val have to be given in "
			                "same shape!")
		self._lon = lon
		self._lat = lat
		self._val = val
		self._filename = filename
		self._remove_file = False
		self._has_grid = False
		self._rect = rect
		self._registration = registration
		if registration not in ['gridline','pixel']:
			raise ValueError("Unknown registration '" + str(registration) + "'. "
			                 "Must be either 'gridline' or 'pixel'.")
		if not _no_grid_check or fast_index is None or Nlon is None or \
		    Nlat is None:
		   grid_property = _is_gridded(lon, lat)
		   fast_index = grid_property["fast_index"]
		   Nlon = grid_property["Nlon"]
		   Nlat = grid_property["Nlat"]
		   if not grid_property["is_gridded"] and \
		       not (np.issorted(lon) and np.issorted(lat)):
		       # Q&D
		       raise ValueError("For non-gridded data, lon and lat "
		                        "have to be sorted!")
		
		
		self._fast_index = fast_index
		self._Nlon = Nlon
		self._Nlat = Nlat
	
	def __del__(self):
		if hasattr(self,'_remove_file') and self._remove_file:
			remove(self._filename)
		
	
	def _create_grd_file(self):
		
		# Step 1: Create temporary csv:
		if not _is_gridded(self._lon, self._lat)["is_gridded"]:
			lon, lat = np.meshgrid(self._lon, self._lat,
			                       indexing='ij')
			dlon = self._lon[1]-self._lon[0]
			dlat = self._lat[1]-self._lat[0]
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
				if self._fast_index == 'LAT':
					if np.count_nonzero([x > 1 for x in self._val.shape]) == 1:
						val = self._val.reshape((self._Nlon,self._Nlat)).T.flatten()
					else:
						val = self._val.flatten()
				else:
					val = self._val.flatten()
				if np.any(np.isnan(self._val)):
					# Handle NaN:
					nanvalue = int(np.floor(self._val[~np.isnan(self._val)].min()-1))
					val[np.isnan(val)] = nanvalue
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
				# -F : Pixel registration.
				args = "%s -G%s -ZBL -I%.5f/%.5f -V -R" \
					   % (fcsv.name, self._filename, dlon, dlat) \
					   +str(self._rect)
				if nanvalue is not None:
					# Write NaN values to a value outside range:
					args += " -di%i" % (nanvalue)
				if self._registration is 'pixel':
					# Pixel registration:
					args += " -r "
				
				print("args: '",args,"'")
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
