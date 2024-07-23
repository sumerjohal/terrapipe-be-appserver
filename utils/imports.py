import azure
import sys
import os, uuid
import glob
import sys
import shapely.wkt
import shutil
import shapely.ops as ops
import folium
import shapely.wkt
import rasterio
import rasterio.features
import rasterio.warp
import contextily as ctx
import ftplib
import re
import ast
import pylab as P
import os
import pickle
import csv
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
#import shapefile
import urllib
import earthpy as et
import os, fnmatch
import math
import sentinelhub
import geojson
import requests
import json
import shapely.wkt
import csv
import xml.sax, xml.sax.handler
import base64
import pyproj
import zipfile
import multiprocessing
import math
import matplotlib.colors as colors
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from email import encoders
from email.mime.base import MIMEBase
import smtplib, ssl
from eto import ETo, datasets
import wget
import matplotlib.colors as cl
import azure
from azure.storage.blob import BlobServiceClient, BlobClient, ContainerClient
import tarfile
from PIL import Image
from sentinelhub import AwsTile
from sentinelsat import SentinelAPI
from sentinelsat.sentinel import SentinelAPI, read_geojson, geojson_to_wkt
from shapely.geometry import shape, Point, Polygon
from zipfile import ZipFile
from rasterio.enums import Resampling
from s2cloudless import S2PixelCloudDetector, CloudMaskRequest
from sentinelhub import WmsRequest, BBox, CRS, MimeType, CustomUrlParam, get_area_dates

from shapely import wkt
from azure.storage.blob import BlobServiceClient, BlobClient, ContainerClient
from shapely.geometry import mapping, Point, Polygon
from ftplib import FTP
from pathlib import Path
from functools import partial
from PIL import Image
from shapely.wkt import dumps, loads
from rasterio import plot
from rasterio import warp
from rasterio import mask
from matplotlib import cm
from rasterio import Affine as A
from rasterio.warp import reproject, Resampling
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from datetime import date
from datetime import datetime 
from datetime import timedelta  
from dateutil import parser
from rasterio import plot as rioplot
from area import area
from pandas import json_normalize
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from matplotlib import cm
from shapely.geometry import Point, mapping, shape
from fiona import collection
from fiona.crs import from_epsg
from pandas import ExcelWriter
from pandas import ExcelFile
from geojson import Polygon
from sentinelhub import AwsTile
from sentinelsat.sentinel import SentinelAPI, read_geojson, geojson_to_wkt
from datetime import date
from collections import OrderedDict
from sentinelsat import SentinelAPI
from rasterio import plot
from pyproj import Proj, transform
from PIL import Image
from zipfile import ZipFile
from os import listdir
from os.path import isfile, join
from time import gmtime, strftime
from matplotlib.pyplot import *
from rasterio.warp import calculate_default_transform, reproject, Resampling
from html.parser import HTMLParser

import plotly
import plotly.express as px
import plotly.graph_objects as go
import chart_studio.plotly as py
from plotly.subplots import make_subplots
import statsmodels
from scipy import stats
from shapely.ops import nearest_points
import xmltodict

import sqlalchemy
import decimal
from geoalchemy2.shape import to_shape
from sqlalchemy import create_engine
from geoalchemy2 import Geometry, WKTElement
from sqlalchemy import *

from bs4 import BeautifulSoup

#######
from raster2points import raster2df
import parquet_tools
import pyarrow
import pyarrow.parquet as pq
import gc
import pyarrow.dataset as ds
import dask.dataframe as dd 
import tarfile
import wget
import s2sphere as s2
import json
import fastparquet
import matplotlib.colors as cl
from geopandas.tools import sjoin
import warnings
warnings.filterwarnings('ignore')

######### Imports for H3

#import h3
#import h3.api.basic_int as h3int
#import h3pandas

# np.warnings.filterwarnings('ignore')  
gpd.io.file.fiona.drvsupport.supported_drivers['KML'] = 'rw'

_type_py2sql_dict = {
 'object': sqlalchemy.sql.sqltypes.Unicode,
 int: sqlalchemy.sql.sqltypes.BigInteger,
 str: sqlalchemy.sql.sqltypes.Unicode,
 float: sqlalchemy.sql.sqltypes.Float,
 'float64': sqlalchemy.sql.sqltypes.Float,
 decimal.Decimal: sqlalchemy.sql.sqltypes.Numeric,
 datetime: sqlalchemy.sql.sqltypes.DateTime,
 bytes: sqlalchemy.sql.sqltypes.LargeBinary,
 bool: sqlalchemy.sql.sqltypes.Boolean,
 date: sqlalchemy.sql.sqltypes.Date,
 time: sqlalchemy.sql.sqltypes.Time,
 timedelta: sqlalchemy.sql.sqltypes.Interval,
 list: sqlalchemy.sql.sqltypes.ARRAY,
 dict: sqlalchemy.sql.sqltypes.JSON
}
#datetime.datetime: sqlalchemy.sql.sqltypes.DateTime,
#datetime.timedelta: sqlalchemy.sql.sqltypes.Interval,
#datetime.date: sqlalchemy.sql.sqltypes.Date,
#datetime.time: sqlalchemy.sql.sqltypes.Time,
