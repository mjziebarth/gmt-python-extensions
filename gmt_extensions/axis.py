# Supply a simple axis.

class Axis:
	"""
	Make plotting a little bit like matplotlib.
	"""
	def __init__(self, fig, region, projection, frame=True):
		self._region = region
		self._projecion = projection
		self._fig = fig
		self._frame = frame
	
	def coast(self, land=None, water=None, shownationalborders=False,
	          showstateborders=False):
		# Preprocess:
		kwargs = dict()
		
		# Borders (-N):
		if not (shownationalborders is False and showstateborders is False):
			kwargs['borders'] = ''
			if shownationalborders is True:
				kwargs['borders'] += '1'
			if showstateborders is True:
				kwargs['borders'] += '2'
		self._fig.coast(region=self._region, projection=self._projection,
		                land=land, water=water, **kwargs)
	
	def grdimage(self, grid, cmap):
		if 
