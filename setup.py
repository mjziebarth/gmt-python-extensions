#!/bin/python
from setuptools import setup


extensions=[]

setup(
	name='gmt-extensions',
	packages=['gmt_extensions'],
#	requires=['numpy (>=1.8)', 'scipy (>=0.14)',
#		      'pythonigraph (>=0.7)'],
	provides=['gmt_extensions'],
	scripts=[],
	license='BSD',
	)
