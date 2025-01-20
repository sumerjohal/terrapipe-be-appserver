import sys
import os
#sys.path.append('/home/user/terrapipe')
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from utils.imports import *
from utils.settings import *
from utils.agls_utilities import *


from flask import Flask, render_template, jsonify, request, url_for
import hashlib
import h3
import plotly
import plotly.express as px
import plotly.graph_objects as go
import numpy as np
import urllib
import urllib.request
import urllib.parse as urlparse
from urllib.parse import urlencode
from flask import make_response
import folium
import geojson
import s2sphere as s2
import uuid
import pandas as pd
import geopandas as gpd
from sqlalchemy import create_engine
import shapely
import shapely.wkt as wkt
from geoalchemy2 import Geometry, WKTElement
import warnings
import time
import json
from datetime import datetime
from datetime import date, timedelta, datetime ,timezone
from dateutil.parser import parse
from shapely.geometry import mapping, Point, Polygon
import pyarrow.dataset as ds
import matplotlib.pyplot as plt
import matplotlib
from os import path, walk
import pandas as pd
import json
import s2sphere as s2
import s2sphere
from concurrent.futures import ThreadPoolExecutor
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler
import pyarrow as pa
from s2sphere import CellId
import threading
import logging
from s2sphere import CellId, Cell, LatLng
import shapely.geometry
from flask_caching import Cache
from functools import wraps
from flask_cors import CORS

import concurrent.futures
import pyarrow.parquet as pq
import signal
import pytz
from datetime import timedelta
from collections import Counter
from shapely.geometry import Polygon as ShapelyPolygon


warnings.filterwarnings('ignore')

################ GLOBALS
#global res_gdf
#global pa_filePath
from fiona.drvsupport import supported_drivers
supported_drivers['LIBKML'] = 'rw'

pa_filePath = '/network/SENTINEL/PARQUET_NDVI_L20_3/'
tileStr = '11SKA'
fieldDataConfigRoot = '/home/rnaura/terrapipe/config/'
fieldDataOutputRoot = '/home/rnaura/terrapipe/OUTPUT/'

geoidLookupFile = fieldDataConfigRoot+'geoid_lookup.json'

fieldDataConfigPath = fieldDataConfigRoot + tileStr+ '/'
resLevel = 20
s2_max_resolution = resLevel
maxResColName = 's2_index__L'+str(resLevel)

logging.basicConfig(
    filename='weather_data_processing.log',  # Log file name
    level=logging.INFO,  # Log level (INFO, DEBUG, ERROR, etc.)
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',  # Log message format
    datefmt='%Y-%m-%d %H:%M:%S'  # Date format
)# cache = Cache(app)


# config = {
#     "DEBUG": True,          # some Flask specific configs
#     "CACHE_TYPE": "SimpleCache",  # Flask-Caching related configs
#     "CACHE_DEFAULT_TIMEOUT": 300
# }
app = Flask(__name__)
CORS(app)
# tell Flask to use the above defined config
# app.config.from_mapping(config)
# cache = Cache(app)

app.config["TEMPLATES_AUTO_RELOAD"] = False


""" logs config """

log_directory = 'logs'
os.makedirs(log_directory, exist_ok=True)

# Configure logging
logging.basicConfig(level=logging.INFO)
general_log_file = os.path.join(log_directory, 'general.log')
error_log_file = os.path.join(log_directory, 'error.log')

# General log handler
general_handler = logging.FileHandler(general_log_file)
general_handler.setLevel(logging.INFO)
general_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))

# Error log handler
error_handler = logging.FileHandler(error_log_file)
error_handler.setLevel(logging.ERROR)
error_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))

# Adding handlers to the logger
logger = logging.getLogger(__name__)
logger.addHandler(general_handler)
logger.addHandler(error_handler)




################ DICTIONARIES
property_dict = {
    'bdod': 'Bulk density of the fine earth fraction',
    'cec': 'Cation Exchange Capacity of the soil',
    'cfvo': 'Volumetric fraction of coarse fragments (> 2 mm)',
    'clay': 'Proportion of clay particles (< 0.002 mm) in the fine earth fraction',
    'sand': 'Proportion of sand particles (> 0.05 mm) in the fine earth fraction',
    'silt': 'Proportion of silt particles (≥ 0.002 mm and ≤ 0.05 mm) in the fine earth fraction',
    'nitrogen': 'Total nitrogen (N)', 
    'ocd': 'Organic carbon density',
    'ocs': 'Organic carbon stocks',
    'phh2o': 'Soil pH',
    'soc': 'Soil organic carbon content in the fine earth fraction'
}
soilsList = [
'Acrisols',
'Alisols',
'Andosols',
'Anthrosols',
'Arenosols',
'Calcisols',
'Cambisols',
'Chernozem',
'Cryosols',
'Durisols',
'Ferralsols',
'Fluvisols',
'Gleysols',
'Gypsisols',
'Histosol',
'Kastanozems',
'Leptosols',
'Lixisols',
'Luvisols',
'Nitisols',
'Phaeozems',
'Planosols',
'Plinthosols',
'Podzols',
'Regosols',
'Retisols',
'Solonchaks',
'Solonetz',
'Stagnosols',
'Technosols',
'Umbrisols',
'Vertisols',
]
soils_dict = dict()
for s in soilsList:
    if s[-1]=='s':
        url = 'https://en.wikipedia.org/wiki/' + s[0:-1]
    else:
        url = 'https://en.wikipedia.org/wiki/' + s
    soils_dict.update({s: url})
    
# def getSoilPropertiesDF(propertyCodeStr, data):
    
#     allDepth_df = pd.DataFrame()

#     #first get the subset by name
#     for i in range(len(data['properties']['layers'])):
#         propName = data['properties']['layers'][i]['name']
#         if (propName==propertyCodeStr):
            
#             #get the unit measure
#             um = data['properties']['layers'][i]['unit_measure']
#             fac = um['d_factor']
#             unitStr_target = um['target_units']
#             unitStr_mapped = um['mapped_units']
            
            
#             data_code = data['properties']['layers'][i]
#             num_depths = len(data_code['depths'])
#             for d in range(num_depths):
#                 data_at_depth = data_code['depths'][d]
#                 row_label = data_at_depth['label']
#                 vals = data_at_depth['values']
#                 rng = data_at_depth['range']
                
#                 top_depth = rng['top_depth']
#                 bottom_depth = rng['bottom_depth']
#                 unit_depth = rng['unit_depth']
                
#                 df = pd.DataFrame(list(vals.values())).T
#                 df = df / fac
                
#                 df.columns = list(vals.keys())
#                 cols = ['Depth'] + ['Top_Depth', 'Botton_Depth', 'Units_Depth'] + df.columns.tolist()
                
#                 df['Depth'] = row_label
#                 df['Top_Depth'] = top_depth
#                 df['Botton_Depth'] = bottom_depth
#                 df['Units_Depth'] = unit_depth
                
#                 df = df[cols]
#                 allDepth_df = allDepth_df.append(df, ignore_index=True)
#         else:
#             continue
    
#     return [unitStr_mapped, fac, unitStr_target, propertyCodeStr, property_dict[propertyCodeStr]], allDepth_df


# def getSoilPropertiesDF2(propertyCodeStr, data):
    
#     allDepth_df = pd.DataFrame()

#     #first get the subset by name
#     for i in range(len(data['properties']['layers'])):
#         propName = data['properties']['layers'][i]['name']
#         if (propName==propertyCodeStr):
            
#             #get the unit measure
#             um = data['properties']['layers'][i]['unit_measure']
#             fac = um['d_factor']
#             unitStr_target = um['target_units']
#             unitStr_mapped = um['mapped_units']
            
            
#             data_code = data['properties']['layers'][i]
#             num_depths = len(data_code['depths'])
#             for d in range(num_depths):
#                 data_at_depth = data_code['depths'][d]
#                 row_label = data_at_depth['label']
#                 vals = data_at_depth['values']
#                 rng = data_at_depth['range']
                
#                 top_depth = rng['top_depth']
#                 bottom_depth = rng['bottom_depth']
#                 unit_depth = rng['unit_depth']
                
#                 df = pd.DataFrame(list(vals.values())).T
#                 df = df / fac
                
#                 df.columns = list(vals.keys())
#                 cols = ['Depth'] + ['Top_Depth', 'Botton_Depth', 'Units_Depth'] + df.columns.tolist()
                
#                 df['Depth'] = row_label
#                 df['Top_Depth'] = top_depth
#                 df['Botton_Depth'] = bottom_depth
#                 df['Units_Depth'] = unit_depth
                
#                 df = df[cols]
#                 allDepth_df = allDepth_df.append(df, ignore_index=True)
#         else:
#             continue
    
#     #fix the allDepth_df
#     df = allDepth_df[['Depth','mean']]
#     df = df.rename(columns={'mean':propertyCodeStr})
#     df = df.drop(columns=['Depth'], axis=1)
    
#     return [unitStr_mapped, fac, unitStr_target, propertyCodeStr, property_dict[propertyCodeStr]], df



# def getAWCMean(silt, clay, BD):
#     #https://journals.lww.com/soilsci/Abstract/1985/07000/Estimating_Available_Water_Holding_Capacity_of.7.aspx#:~:text=Available%20water%2Dstorage%20capacity%20(AWSC,(r2%20%3D%200.92)%3A%20AWSCcore%20%3D
#     if not silt.index.name=='Depth':
#         silt.set_index('Depth', drop=True, inplace=True)
#     silt = silt[['mean']]
#     silt = silt.astype('float')
    
#     if not clay.index.name=='Depth':
#         clay.set_index('Depth', drop=True, inplace=True)
#     clay = clay[['mean']]
#     clay = clay.astype('float')
    
#     if not BD.index.name=='Depth':
#         BD.set_index('Depth', drop=True, inplace=True)
#     BD = BD[['mean']]
#     BD = BD.astype('float')
    
#     AWSC = 14.01 + 0.03*(silt * clay) - 8.78*BD
#     AWSC.rename(columns={'mean':'AWSC'}, inplace=True)
#     AWSC = AWSC.round(2)
#     #units are %Volume
#     unitsStr = '%vol [volume-fraction]'
    
#     return unitsStr, AWSC, 'Available water holding capacity'

# def getAWC(silt, clay, BD):
#     #https://journals.lww.com/soilsci/Abstract/1985/07000/Estimating_Available_Water_Holding_Capacity_of.7.aspx#:~:text=Available%20water%2Dstorage%20capacity%20(AWSC,(r2%20%3D%200.92)%3A%20AWSCcore%20%3D
#     if not silt.index.name=='Depth':
#         silt.set_index('Depth', drop=True, inplace=True)
#     s = silt[['Q0.05','Q0.5','Q0.95','mean']]
#     s = s.astype('float')
    
#     if not clay.index.name=='Depth':
#         clay.set_index('Depth', drop=True, inplace=True)
#     c = clay[['Q0.05','Q0.5','Q0.95','mean']]
#     c = c.astype('float')
    
#     if not BD.index.name=='Depth':
#         BD.set_index('Depth', drop=True, inplace=True)
#     bd = BD[['Q0.05','Q0.5','Q0.95','mean']]
#     bd = bd.astype('float')
    
#     AWSC = 14.01 + 0.03*(s * c) - 8.78*bd
    
#     #AWSC.rename(columns={'mean':'AWSC'}, inplace=True)
#     AWSC = AWSC.round(2)
    
#     #add the columns for ['Top_Depth', 'Botton_Depth', 'Units_Depth']
#     cols = ['Top_Depth', 'Botton_Depth', 'Units_Depth']
#     AWSC[cols]=silt[cols]
    
#     return AWSC

# def getKsatMean(silt, sand, clay, BD):
#     #, OC, PD, MC):
#     #PD = Particle density
#     #MC = Moisture Content
#     #https://www.researchgate.net/figure/Saturated-hydraulic-conductivity-of-soil-as-influenced-by-per-cent-silt-th-clay-content-in_fig2_248885151

#     #https://www.tandfonline.com/doi/full/10.1080/24749508.2018.1481633
#     #EIR = -30,578.81–305.56(sand%)-306.16(silt%)-0.306.33(clay%)-5.18(BD%)+.34(MC%)+4.18(PD)+16.85(OC%)

#     #https://essd.copernicus.org/preprints/essd-2020-149/essd-2020-149.pdf
#     #log(Ksat) = b0 + b1 · BD + b2 · BD2 + b3 · CL + b4 · BD · CL + b5 · CL2 + b6 · SA + b7 · BD · SA + b8 · CL · SA + b9 · SA2
#     if not silt.index.name=='Depth':
#         silt.set_index('Depth', drop=True, inplace=True)
    
#     silt = silt[['mean']]
#     silt = silt.astype('float')

#     if not clay.index.name=='Depth':
#         clay.set_index('Depth', drop=True, inplace=True)
#     clay = clay[['mean']]
#     clay = clay.astype('float')

#     if not sand.index.name=='Depth':
#         sand.set_index('Depth', drop=True, inplace=True)
#     sand = sand[['mean']]
#     sand = sand.astype('float')

#     if not BD.index.name=='Depth':
#         BD.set_index('Depth', drop=True, inplace=True)
#     BD = BD[['mean']]
#     BD = BD.astype('float')

#     b0=2.17
#     b1=0.9387
#     b2=-0.8026
#     b3=0.0037
#     b4=-0.017
#     b5=0
#     b6=0.0025
#     b7=0
#     b8=0
#     b9=0
#     #Ksat is in cm/day, clay (CL) and sand (SA) are expressed in % and bulk density (BD) is in g/cm3 or kg/dm3
#     log_Ksat = b0 + b1*BD + b2*BD.pow(2) + b3*clay + b4*BD*clay + b5*clay.pow(2) + b6*sand + b7*BD*sand + b8*clay*sand + b9*sand.pow(2)
#     log_Ksat = log_Ksat[['mean']]
#     Ksat = log_Ksat.apply(lambda x: np.exp(x))
#     Ksat.rename(columns = {'mean':'Ksat'}, inplace=True)
#     units = 'cm/day'
    
#     #Convert to inches / hr
#     Ksat_inchesPerHr = Ksat * 0.0164042
#     units_inchesPerHr = 'in/hr'

#     return units_inchesPerHr, Ksat_inchesPerHr, 'Saturated hydraulic conductivity'

# def getKsat(silt, sand, clay, BD):
#     #, OC, PD, MC):
#     #PD = Particle density
#     #MC = Moisture Content
#     #https://www.researchgate.net/figure/Saturated-hydraulic-conductivity-of-soil-as-influenced-by-per-cent-silt-th-clay-content-in_fig2_248885151

#     #https://www.tandfonline.com/doi/full/10.1080/24749508.2018.1481633
#     #EIR = -30,578.81–305.56(sand%)-306.16(silt%)-0.306.33(clay%)-5.18(BD%)+.34(MC%)+4.18(PD)+16.85(OC%)

#     #https://essd.copernicus.org/preprints/essd-2020-149/essd-2020-149.pdf
#     #log(Ksat) = b0 + b1 · BD + b2 · BD2 + b3 · CL + b4 · BD · CL + b5 · CL2 + b6 · SA + b7 · BD · SA + b8 · CL · SA + b9 · SA2

#     if not silt.index.name=='Depth':
#         silt.set_index('Depth', drop=True, inplace=True)
#     st = silt[['Q0.05','Q0.5','Q0.95','mean']]
#     st = st.astype('float')

#     if not clay.index.name=='Depth':
#         clay.set_index('Depth', drop=True, inplace=True)
#     c = clay[['Q0.05','Q0.5','Q0.95','mean']]
#     c = c.astype('float')

    
#     if not sand.index.name=='Depth':
#         sand.set_index('Depth', drop=True, inplace=True)
#     sa = sand[['Q0.05','Q0.5','Q0.95','mean']]
#     sa = sa.astype('float')

#     if not BD.index.name=='Depth':
#         BD.set_index('Depth', drop=True, inplace=True)
#     bd = BD[['Q0.05','Q0.5','Q0.95','mean']]
#     bd = bd.astype('float')

#     b0=2.17
#     b1=0.9387
#     b2=-0.8026
#     b3=0.0037
#     b4=-0.017
#     b5=0
#     b6=0.0025
#     b7=0
#     b8=0
#     b9=0
#     #Ksat is in cm/day, clay (CL) and sand (SA) are expressed in % and bulk density (BD) is in g/cm3 or kg/dm3
#     log_Ksat = b0 + b1*bd + b2*bd.pow(2) + b3*c + b4*bd*c + b5*c.pow(2) + b6*sa + b7*bd*sa + b8*c*sa + b9*sa.pow(2)
#     log_Ksat = log_Ksat[['Q0.05','Q0.5','Q0.95','mean']]
#     Ksat = log_Ksat.apply(lambda x: np.exp(x))
#     #Ksat.rename(columns = {'mean':'Ksat'}, inplace=True)
#     units = 'cm/day'
    
#     #Convert to inches / hr
#     Ksat_inchesPerHr = Ksat * 0.0164042
#     units_inchesPerHr = 'in/hr'
    
#     cols = ['Top_Depth', 'Botton_Depth', 'Units_Depth']
#     Ksat_inchesPerHr[cols]=silt[cols]

#     return units_inchesPerHr, Ksat_inchesPerHr

# def getSoils(lat,lon):
#     number_classes=10
#     #url_loc = 'https://rest.soilgrids.org/soilgrids/v2.0/classification/query?'+'lon='+str(lon)+'&lat='+str(lat)+'&number_classes='+str(number_classes)
#     url_loc = 'https://rest.isric.org/soilgrids/v2.0/classification/query?'+'lon='+str(lon)+'&lat='+str(lat)+'&number_classes='+str(number_classes)
#     with urllib.request.urlopen(url_loc) as response:
#         loc_data = json.load(response)
    
#     soilsList = loc_data['wrb_class_probability']
#     colsList=['SOIL','PROB','INFO_URL']
#     soils_df = pd.DataFrame()
#     i=0
#     for s in soilsList:
#         soilName = s[0]
#         prob = s[1]
        
#         try:
#             urlStr = soils_dict[soilName]
#             soilLink = '<a href='+urlStr+'>'+sn+'</a>' 
#             #print(soilLink)
#         except:
#             #didn't find the link
#             if soilName[-1]=='s':
#                 sn = soilName[0:-1]
#             else:
#                 sn = soilName
#             urlStr = 'https://en.wikipedia.org/wiki/'+sn
#             soilLink = '<a href='+urlStr+'>'+sn+'</a>' 
#             #print(soilLink)
        
#         soils_df.loc[i,'SOIL']=soilName
#         soils_df.loc[i,'PROB']=prob
#         soils_df.loc[i,'INFO_URL']=soilLink
#         i=i+1
    
    
#     soils_df = soils_df[colsList]
    
#     newCols = ['Soil Type', 'Percentage', 'More Information']
#     soils_df.columns = newCols
    
#     return soils_df

# def getSoilsForGeoid(agstack_geoid):
#     configFile =  fieldDataConfigPath + agstack_geoid+'.json'
#     #NDays = (datetime.today() - datetime.strptime(dtStr,'%Y-%m-%d')).days
    
#     #Validate the date 
#     #1. Dt should be in the last 90 days from today
#     #2. Dt should be a valida date
    
#     s1_time_start = time.time()
#     with open(configFile, "r") as jsonfile:
#         fieldJSON = json.load(jsonfile)
#         field_geoid = fieldJSON['geoid']
#         field_wkt = fieldJSON['wkt']
#     fieldPoly = shapely.wkt.loads(field_wkt)
#     c=fieldPoly.centroid
#     lat=c.y
#     lon=c.x
    
#     number_classes=10
#     #url_loc = 'https://rest.soilgrids.org/soilgrids/v2.0/classification/query?'+'lon='+str(lon)+'&lat='+str(lat)+'&number_classes='+str(number_classes)
#     url_loc = 'https://rest.isric.org/soilgrids/v2.0/classification/query?'+'lon='+str(lon)+'&lat='+str(lat)+'&number_classes='+str(number_classes)
#     with urllib.request.urlopen(url_loc) as response:
#         loc_data = json.load(response)
    
#     soilsList = loc_data['wrb_class_probability']
#     colsList=['SOIL','PROB','INFO_URL']
#     soils_df = pd.DataFrame()
#     i=0
#     for s in soilsList:
#         soilName = s[0]
#         prob = s[1]
        
#         try:
#             urlStr = soils_dict[soilName]
#             soilLink = '<a href='+urlStr+'>'+sn+'</a>' 
#             #print(soilLink)
#         except:
#             #didn't find the link
#             if soilName[-1]=='s':
#                 sn = soilName[0:-1]
#             else:
#                 sn = soilName
#             urlStr = 'https://en.wikipedia.org/wiki/'+sn
#             soilLink = '<a href='+urlStr+'>'+sn+'</a>' 
#             #print(soilLink)
        
#         soils_df.loc[i,'SOIL']=soilName
#         soils_df.loc[i,'PROB']=prob
#         soils_df.loc[i,'INFO_URL']=soilLink
#         i=i+1
    
    
#     soils_df = soils_df[colsList]
    
#     newCols = ['Soil Type', 'Percentage', 'More Information']
#     soils_df.columns = newCols
    
#     return soils_df


# def getPropertiesDFForGeoid(agstack_geoid):
    
#     configFile =  fieldDataConfigPath + agstack_geoid+'.json'
#     #NDays = (datetime.today() - datetime.strptime(dtStr,'%Y-%m-%d')).days
    
#     #Validate the date 
#     #1. Dt should be in the last 90 days from today
#     #2. Dt should be a valida date
    
#     s1_time_start = time.time()
#     with open(configFile, "r") as jsonfile:
#         fieldJSON = json.load(jsonfile)
#         field_geoid = fieldJSON['geoid']
#         field_wkt = fieldJSON['wkt']
#     fieldPoly = shapely.wkt.loads(field_wkt)
#     c=fieldPoly.centroid
#     lat=c.y
#     lon=c.x
    
#     #API call for the details about a point
#     """
#     Query a single pixel point on the soilgrids stack, returning a GeoJSON
#     layer: soilgrids layer name to be queried
#     depth: specific depth to be queried
#     values: statistical values Optional[List] = LayerQuery
#     """
    
#     #API #2 is meta data
#     url_layers = 'https://rest.isric.org/soilgrids/v2.0/properties/layers'
#     with urllib.request.urlopen(url_layers) as response:
#         layer_data = json.load(response)
#     #get the list of properties
#     propertiesList = []
#     length = len(layer_data['layers'])
#     for idx in range(length):
#         p = layer_data['layers'][idx]['property']
#         propertiesList.append(p)
    
    
#     #to get the meta data for any property
#     depths=[]
#     values = []
#     #modify the property list
#     propertiesList = ['bdod',
#      'cec',
#      'cfvo',
#      'soc',
#      'nitrogen',
#      'ocd',
#      'phh2o',
#      'clay',
#      'sand',
#      'silt']
    
#     prop_to_find = 'nitrogen'
#     for idx in range(length):
#         p = layer_data['layers'][idx]['property']
#         if (p==prop_to_find):
#             info = layer_data['layers'][idx]['layer_structure']
#             for i in range(len(info)):
#                 depths.append(info[i]['range'])
#                 values.append(info[i]['values'])
#     valuesList = values[1]
#     prop_url = ''
#     for p in propertiesList:
#         prop_url = prop_url + '&property='+str(p)

#     value_url = ''
#     valuesList=['mean']
#     for v in valuesList:
#         value_url = value_url + '&value='+str(v)

#     depth_url = ''
#     depths=['0-5cm']
#     for d in depths:
#         depth_url = depth_url + '&depth='+str(d)

#     main_url = 'https://rest.isric.org/soilgrids/v2.0/properties/query?' + 'lon='+str(lon)+'&lat='+str(lat)
#     url_details = main_url + prop_url + depth_url + value_url
#     with urllib.request.urlopen(url_details) as response:
#         data = json.load(response) 

    

#     propertyResult = pd.DataFrame()
#     i=0
#     for p in propertiesList:
#         prop_mean_value = getSoilPropertiesDF2(p, data)[1].iloc[0][0]
#         prop_units = getSoilPropertiesDF2(p, data)[0][2]
#         prop_desc = getSoilPropertiesDF2(p, data)[0][4]
#         propertyResult.loc[i,'Name']=p
#         propertyResult.loc[i,'Mean Value in Top-Soil']=prop_mean_value
#         propertyResult.loc[i,'Units']=prop_units
#         propertyResult.loc[i,'Description']=prop_desc
#         i=i+1

#     #Lets; do KSat
#     ksat_tuple = getKsatMean(getSoilPropertiesDF('silt', data)[1], getSoilPropertiesDF('sand', data)[1], getSoilPropertiesDF('clay', data)[1], getSoilPropertiesDF('bdod', data)[1])
    
#     nameStr = ksat_tuple[1].columns.tolist()[0]
#     valStr = ksat_tuple[1].iloc[0][0]
#     unitsStr = ksat_tuple[0]
#     descStr = ksat_tuple[2]
#     propertyResult.loc[len(propertyResult.index)] = [nameStr, valStr, unitsStr, descStr]

#     #Now let's do AWC
#     awc_tuple = getAWCMean(getSoilPropertiesDF('silt', data)[1], getSoilPropertiesDF('clay', data)[1], getSoilPropertiesDF('bdod', data)[1])
#     nameStr = awc_tuple[1].columns.tolist()[0]
#     valStr = awc_tuple[1].iloc[0][0]
#     unitsStr = awc_tuple[0]
#     descStr = awc_tuple[2]
#     propertyResult.loc[len(propertyResult.index)] = [nameStr, valStr, unitsStr, descStr]
    
#     renamedCols = ['Acronym', 'Mean Value in Top-Soil', 'Units', 'Property']
#     propertyResult.columns=renamedCols
#     subsetCols = ['Property', 'Mean Value in Top-Soil', 'Units']
#     propertyResult = propertyResult[subsetCols]
    
#     #round to 2 decimal places for Mean Value
#     propertyResult['Mean Value in Top-Soil'] = propertyResult['Mean Value in Top-Soil'].round(2)
    
#     return propertyResult


# def getPropertiesDF(lat,lon):
    
#     #API call for the details about a point
#     """
#     Query a single pixel point on the soilgrids stack, returning a GeoJSON
#     layer: soilgrids layer name to be queried
#     depth: specific depth to be queried
#     values: statistical values Optional[List] = LayerQuery
#     """
    
#     #API #2 is meta data
#     url_layers = 'https://rest.isric.org/soilgrids/v2.0/properties/layers'
#     with urllib.request.urlopen(url_layers) as response:
#         layer_data = json.load(response)
#     #get the list of properties
#     propertiesList = []
#     length = len(layer_data['layers'])
#     for idx in range(length):
#         p = layer_data['layers'][idx]['property']
#         propertiesList.append(p)
    
    
#     #to get the meta data for any property
#     depths=[]
#     values = []
#     #modify the property list
#     propertiesList = ['bdod',
#      'cec',
#      'cfvo',
#      'soc',
#      'nitrogen',
#      'ocd',
#      'phh2o',
#      'clay',
#      'sand',
#      'silt']
    
#     prop_to_find = 'nitrogen'
#     for idx in range(length):
#         p = layer_data['layers'][idx]['property']
#         if (p==prop_to_find):
#             info = layer_data['layers'][idx]['layer_structure']
#             for i in range(len(info)):
#                 depths.append(info[i]['range'])
#                 values.append(info[i]['values'])
#     valuesList = values[1]
#     prop_url = ''
#     for p in propertiesList:
#         prop_url = prop_url + '&property='+str(p)

#     value_url = ''
#     valuesList=['mean']
#     for v in valuesList:
#         value_url = value_url + '&value='+str(v)

#     depth_url = ''
#     depths=['0-5cm']
#     for d in depths:
#         depth_url = depth_url + '&depth='+str(d)

#     main_url = 'https://rest.isric.org/soilgrids/v2.0/properties/query?' + 'lon='+str(lon)+'&lat='+str(lat)
#     url_details = main_url + prop_url + depth_url + value_url
#     with urllib.request.urlopen(url_details) as response:
#         data = json.load(response) 

    

#     propertyResult = pd.DataFrame()
#     i=0
#     for p in propertiesList:
#         prop_mean_value = getSoilPropertiesDF2(p, data)[1].iloc[0][0]
#         prop_units = getSoilPropertiesDF2(p, data)[0][2]
#         prop_desc = getSoilPropertiesDF2(p, data)[0][4]
#         propertyResult.loc[i,'Name']=p
#         propertyResult.loc[i,'Mean Value in Top-Soil']=prop_mean_value
#         propertyResult.loc[i,'Units']=prop_units
#         propertyResult.loc[i,'Description']=prop_desc
#         i=i+1

#     #Lets; do KSat
#     ksat_tuple = getKsatMean(getSoilPropertiesDF('silt', data)[1], getSoilPropertiesDF('sand', data)[1], getSoilPropertiesDF('clay', data)[1], getSoilPropertiesDF('bdod', data)[1])
    
#     nameStr = ksat_tuple[1].columns.tolist()[0]
#     valStr = ksat_tuple[1].iloc[0][0]
#     unitsStr = ksat_tuple[0]
#     descStr = ksat_tuple[2]
#     propertyResult.loc[len(propertyResult.index)] = [nameStr, valStr, unitsStr, descStr]

#     #Now let's do AWC
#     awc_tuple = getAWCMean(getSoilPropertiesDF('silt', data)[1], getSoilPropertiesDF('clay', data)[1], getSoilPropertiesDF('bdod', data)[1])
#     nameStr = awc_tuple[1].columns.tolist()[0]
#     valStr = awc_tuple[1].iloc[0][0]
#     unitsStr = awc_tuple[0]
#     descStr = awc_tuple[2]
#     propertyResult.loc[len(propertyResult.index)] = [nameStr, valStr, unitsStr, descStr]
    
#     renamedCols = ['Acronym', 'Mean Value in Top-Soil', 'Units', 'Property']
#     propertyResult.columns=renamedCols
#     subsetCols = ['Property', 'Mean Value in Top-Soil', 'Units']
#     propertyResult = propertyResult[subsetCols]
    
#     #round to 2 decimal places for Mean Value
#     propertyResult['Mean Value in Top-Soil'] = propertyResult['Mean Value in Top-Soil'].round(2)
    
#     return propertyResult

# ################  FUNCTIONS
def getLocalNDVIDatasets(s2_index__L8_list, s2_index__L10_list, pa_filePath):

    list_of_paths = []
    list_of_paths_L8 = [pa_filePath+'s2_index__L8='+x for x in s2_index__L8_list]
    for p in list_of_paths_L8:
        list_of_paths = list_of_paths + [p+'/s2_index__L10='+x for x in s2_index__L10_list]
        print(list_of_paths)
    
    datasets = []
    for p in list_of_paths:
        print(p)
        if os.path.exists(p):
            datasets.append( ds.dataset(p,format="parquet", partitioning="hive") )
        else:
            continue
    print(f'datasets--{datasets}')
    return datasets

def getNDVIFromParquetForDataSetsAndDate(datasets, dtStr, s2_index__L10_list, s2_index__L15_list, s2_index__L20_list):
    start_time = time.time()
    subset_dirs = []
    endDt = datetime.now()
    
    
    endDt =  datetime.strptime(dtStr, '%Y-%m-%d')
    startDt = endDt - timedelta(days=10)
    NDays = (endDt - startDt).days
    
    dtStrs_list = [datetime.strftime(x, '%Y_%m_%d') for x in pd.date_range(startDt, periods=NDays).tolist()]
    
    print(f'dtStrs_list--{dtStrs_list}')
    
    p_df = pd.DataFrame()
    for dtset in datasets:
        print(f"Processing dataset: {dtset}")

        p_subset_df = dtset.to_table(
            columns=['s2_index__L20', 'NDVI', 'UTC_DATETIME','YY_MM_DD'],
            filter=(
                ds.field('s2_index__L20').isin(s2_index__L20_list) &
                ds.field('YY_MM_DD').isin(dtStrs_list)
            )).to_pandas()
        print(f"Columns in p_subset_df: {p_subset_df.columns}")
        p_df = pd.concat([p_df, p_subset_df], ignore_index=True)

    #get the closest dt and screen for it
    
    print(p_df.columns)
    list_of_dates = pd.unique(p_df.YY_MM_DD).tolist()
    list_of_dates.sort()
    print(f'list_of_dates--{list_of_dates}')
    closest_dtStr = list_of_dates[-1]
    
    p_df = p_df[p_df.YYYY_MM_DD==closest_dtStr].reset_index(drop=True)
    
    if len(p_df)==0:
        print('Empty Dataset')
    else:
        p_df.reset_index(drop=True, inplace=True)
        
    end_time = time.time()
    time_elapsed = (end_time - start_time)
    
   

    return p_df, time_elapsed


# def getNDVIFromParquetForDataSets(datasets, NDays, s2_index__L10_list, s2_index__L15_list, s2_index__L20_list):
#     start_time = time.time()
#     subset_dirs = []
#     endDt = datetime.now()

#     startDt = endDt - timedelta(days=NDays) 
#     dtStrs_list = [datetime.strftime(x, '%Y-%m-%d') for x in pd.date_range(startDt, periods=NDays).tolist()]
    
#     p_df = pd.DataFrame()
#     for dtset in datasets:
#         try:
#             p_subset_df = dtset.to_table(
#                 columns=['s2_index__L20', 'NDVI', 'UTC_DATETIME','YY_MM_DD'],
#                 filter=(
#                     ds.field('s2_index__L20').isin(s2_index__L20_list) &
#                     ds.field('YY_MM_DD').isin(dtStrs_list)
#                 )).to_pandas()
#             p_df = p_df.append(p_subset_df, ignore_index=True)
#         except:
#             continue
    
#     if len(p_df)==0:
#         print('Empty Dataset')
#     else:
#         p_df.reset_index(drop=True, inplace=True)
        
#     end_time = time.time()
#     time_elapsed = (end_time - start_time)

#     return p_df, time_elapsed
# def getStats(dt_df,colName):
#     colName='NDVI'
#     dt_df[colName] = dt_df[colName].astype('float')
    
#     arr = np.array(dt_df[colName])
#     arr[arr>1]=np.nan
#     arr[arr<0]=np.nan

#     #stats
#     p5 = np.nanpercentile(arr,5)
#     p10 = np.nanpercentile(arr,10)
#     p25 = np.nanpercentile(arr,25)
#     p33 = np.nanpercentile(arr,33)
#     p50 = np.nanpercentile(arr,50)
#     p66 = np.nanpercentile(arr,66)
#     p75 = np.nanpercentile(arr,75)
#     p90 = np.nanpercentile(arr,90)
#     p95 = np.nanpercentile(arr,95)
#     m = np.nanmean(arr)
#     stdDev = np.nanstd(arr)
#     mn = np.nanmin(arr)
#     mx = np.nanmax(arr)
#     unif1 = p25/m
#     unif2 = stdDev/m
    
#     stats_dict = {
#         'p5': p5,
#         'p10': p10,
#         'p25': p25, 
#         'p33': p33, 
#         'p50': p50,
#         'p66': p66,
#         'p75': p75,
#         'p90': p90,
#         'p95': p95,
#         'Mean': m,
#         'Std': stdDev, 
#         'Min': mn,
#         'Max': mx,
#         'Uniformity-1': unif1,
#         'Uniformity-2': unif2
#     }
    
#     return stats_dict


# def getBoundaryCoverageNew(cids, poly, maxResColName):
#     s2_index__L19_list = []
#     p_gdf = gpd.GeoDataFrame()
#     idx=0
#     for cellid in cids:
#         c = s2.Cell(cellid)
#         level = cellid.level()
#         vertices = []
#         for i in range(0,4):
#             vertex = c.get_vertex(i)
#             latlng = s2.LatLng.from_point(vertex)
#             vertices.append((latlng.lng().degrees, latlng.lat().degrees))
#         geo = shapely.geometry.Polygon(vertices)
#         if poly.intersects(geo):
#             s2_index__L19_list.append(cellid.to_token())
#             p_gdf.loc[idx,maxResColName] = cellid.to_token()
#             p_gdf.loc[idx,'geometry']=geo
#         idx+=1
    
#     p_gdf.reset_index(drop=True, inplace=True)
#     return p_gdf


# def isValidPolygon(user_fieldWKT):
#     try:
#         poly = shapely.wkt.loads(user_fieldWKT)
#         if poly.geom_type == 'Polygon':
#             return True
#         else:
#             return False
#     except Exception as e:
#         print(e)
#         return False

# def getNDVIStatsFromGeoid(agstack_geoid, startDtStr, endDtStr):
    
#     #checks LATER
#     #1. to make sure startDtStr is before endDtStr
#     #2. to make sure they are valida dates and within the date scope
#     #3. Dt should be in the last 180 days from today
    
#     time_start = time.time()
    
#     NDays = (datetime.today() - datetime.strptime(startDtStr,'%Y-%m-%d')).days
    
    
#     configFile =  fieldDataConfigPath + agstack_geoid+'.json'
    
#     s1_time_start = time.time()
#     with open(configFile, "r") as jsonfile:
#         fieldJSON = json.load(jsonfile)
#         field_geoid = fieldJSON['geoid']
#         field_wkt = fieldJSON['wkt']
#         p20_gdf = gpd.read_file(fieldJSON['s2_L20_index_gdfJSON'], driver='GeoJSON')
#         p19_gdf = gpd.read_file(fieldJSON['s2_L19_index_gdfJSON'], driver='GeoJSON')
#         bound_gdf =  gpd.read_file(fieldJSON['boundary_gdfJSON'], driver='GeoJSON') 
#         field_tileStr = fieldJSON['sentinel-tile-list']
#         field_s2indices = fieldJSON['indices']
#     p_gdf = p20_gdf.copy()
    
#     #get the s2 indices
#     s2_index__L8_list = json.loads(field_s2indices['s2_index__L8_list'])
#     s2_index__L10_list = json.loads(field_s2indices['s2_index__L10_list'])
#     s2_index__L15_list = json.loads(field_s2indices['s2_index__L15_list'])
#     s2_index__L19_list = json.loads(field_s2indices['s2_index__L19_list'])
#     s2_index__L20_list = json.loads(field_s2indices['s2_index__L20_list'])
    
#     #get the NDVI dataset
#     datasets_ndvi = getLocalNDVIDatasets(s2_index__L8_list, s2_index__L10_list, pa_filePath)
#     time_end = time.time()
#     time_elapsed = time_end  - time_start
#     #print('load dataset directory: '+str(round(time_elapsed,2)))
    
#     #read the data
#     time_start = time.time()
#     s2_index__maxRes_list = eval('s2_index__L'+str(s2_max_resolution)+'_list')
#     maxResColName = 's2_index__L'+str(s2_max_resolution)
    
#     ndvi_df, t_elapsed = getNDVIFromParquetForDataSets(
#         datasets_ndvi,
#         NDays, 
#         s2_index__L10_list, 
#         s2_index__L15_list,
#         s2_index__L20_list
#     )
#     time_end = time.time()
#     time_elapsed = time_end  - time_start
#     #print('get NDVI data: '+str(round(time_elapsed,2)))
    
#     #Screen for the dates that are within the date range specified
#     n_df = ndvi_df[(ndvi_df.YY_MM_DD>=startDtStr) & (ndvi_df.YY_MM_DD<=endDtStr)].reset_index(drop=True)
    
#     #for each dt, get the stats
#     list_of_dates = pd.unique(n_df.YY_MM_DD).tolist()
#     list_of_dates.sort()
    
#     all_dates_df = pd.DataFrame()
#     for dtStr in list_of_dates:
#         dt_df = n_df[n_df.YY_MM_DD==dtStr].reset_index(drop=True)
#         stats_dict = getStats(dt_df,'NDVI')
#         stats_df = pd.DataFrame([stats_dict])
#         stats_df['Date'] = [dtStr]
#         all_dates_df = pd.concat([all_dates_df, stats_df], axis=0)
    
#     #all_dates_df.Date = pd.to_datetime(all_dates_df.Date)
#     #all_dates_df.Date = all_dates_df.Date.dt.strftime('%Y-%m-%d')
#     all_dates_df.set_index('Date', inplace=True, drop=True)

#     return all_dates_df


# #get available dates for NDVI
# def getAvailableListOfDatesForNDVI(agstack_geoid):
#     configFile =  fieldDataConfigPath + agstack_geoid+'.json'
#     #NDays = (datetime.today() - datetime.strptime(dtStr,'%Y-%m-%d')).days
    
#     #Validate the date 
#     #1. Dt should be in the last 90 days from today
#     #2. Dt should be a valida date
    
#     s1_time_start = time.time()
#     with open(configFile, "r") as jsonfile:
#         fieldJSON = json.load(jsonfile)
#         field_geoid = fieldJSON['geoid']
#         field_wkt = fieldJSON['wkt']
#         field_s2indices = fieldJSON['indices']
        
#     #get the s2 indices
#     s2_index__L8_list = json.loads(field_s2indices['s2_index__L8_list'])
#     s2_index__L10_list = json.loads(field_s2indices['s2_index__L10_list'])
#     s2_index__L20_list = json.loads(field_s2indices['s2_index__L20_list'])

#     #get the datasets
#     datasets_ndvi = getLocalNDVIDatasets(s2_index__L8_list, s2_index__L10_list, pa_filePath)
    
#     #get the NDVI and dts
#     p_df = pd.DataFrame()
#     for dtset in datasets_ndvi:
#         try:
#             p_subset_df = dtset.to_table(
#                 columns=['s2_index__L20', 'NDVI', 'UTC_DATETIME','YY_MM_DD'],
#                 filter=(
#                     ds.field('s2_index__L20').isin(s2_index__L20_list)
#                 )).to_pandas()
#             p_df = p_df.append(p_subset_df, ignore_index=True)
#         except:
#             continue
    
#     #get the closest dt and screen for it
#     list_of_dates = pd.unique(p_df.YY_MM_DD).tolist()
#     list_of_dates.sort()
    
#     return list_of_dates

# def getInterpolatedValueFromEdgeNeightbors_ndvi(tkn, res_gdf):
#     try:
#         neighboring_cells = s2.CellId.from_token(tkn).get_edge_neighbors()
#         arr_vals = [res_gdf[res_gdf['s2_index__L20']==x.to_token()].NDVI.values.tolist() for x in neighboring_cells]
#         est_val = np.nanmean([x[0] for x in arr_vals if len(x)>0])
#     except:
#         est_val = np.nan
       
#     return est_val


def fetchGeometryForS2Token(s2_tokens):
    """
    Convert S2 tokens to geometries (polygons).
    
    :param s2_tokens: List of S2 cell tokens as strings (e.g., ['4c4c', '7f8f'])
    :return: List of Polygon geometries for the S2 cells
    """
    polygons = []
    for s2_token in s2_tokens:
        try:
            # Convert the S2 token to a CellId
            cell_id = CellId.from_token(s2_token)
            
            # Initialize an S2Cell
            s2_cell = Cell(cell_id)

            # Get the four vertices of the S2 cell
            vertices = [s2_cell.get_vertex(i) for i in range(4)]
            
            # Convert vertices to LatLng
            lat_lngs = [LatLng.from_point(vertex) for vertex in vertices]
            
            # Extract latitude and longitude as floats
            lat_lng_coords = [(lat_lng.lat().degrees, lat_lng.lng().degrees) for lat_lng in lat_lngs]
                        
            # Convert LatLngs to Shapely Polygon
            polygon = ShapelyPolygon(lat_lng_coords)
            
            polygons.append(polygon)
        except Exception as e:
            print(f"Error processing S2 token {s2_token}: {e}")
    
    return polygons


def getBoundaryGeoDataFrame(s2_tokens_list):
    """
    Generate a GeoDataFrame containing boundary polygons for each S2 token.
    
    :param s2_tokens_list: List of S2 cell tokens
    :return: GeoDataFrame with geometry column containing S2 token boundary polygons
    """
    geometries = []
    for token in s2_tokens_list:
        geometry = fetchGeometryForS2Token(token)
        geometries.append(geometry)
    
    # Create a GeoDataFrame with the geometry
    gdf = gpd.GeoDataFrame(geometry=geometries)
    return gdf

def parse_date_to_utc(date_str):
    """
    Parse a date string in 'YYYY-MM-DD' format to a UTC datetime object.
    
    :param date_str: Date string in the format 'YYYY-MM-DD'
    :return: UTC localized datetime object
    """
    local_dt = datetime.strptime(date_str, '%Y-%m-%d')
    return pytz.utc.localize(local_dt)

from shapely.geometry import Polygon

def reorder_coordinates(polygon):
    coords = list(polygon.exterior.coords)
    reordered_coords = [(coord[1], coord[0]) for coord in coords]  # Swap lat, lon to lon, lat
    return Polygon(reordered_coords)
  
    
    
def getGeoDataFrameWithNDVI(s_df, polygons):
    if len(s_df) != len(polygons):
        raise ValueError("The number of NDVI values does not match the number of polygons")

    reordered_polygons = [reorder_coordinates(poly) for poly in polygons]

    # Check if the polygons are valid
    for i, poly in enumerate(reordered_polygons):
        if not poly.is_valid:
            print(f"Invalid polygon at index {i}: {poly}")

    geo_df = pd.DataFrame({
        'NDVI': s_df['average_NDVI'],
        'geometry': reordered_polygons
    })
    
    geo_gdf = gpd.GeoDataFrame(geo_df, geometry='geometry')
    
    # Check for invalid geometries
    invalid_geometries = geo_gdf[~geo_gdf.is_valid]
    if not invalid_geometries.empty:
        print("Invalid geometries found:", invalid_geometries)

    return geo_gdf


# NDVI IMAGE

def getNDVIgdfFromGeoid(agstack_geoid, dt_str):
    """
    Fetch NDVI data and geometries for a given geoid and date, then create a GeoDataFrame.
    """
    # Convert date string to UTC datetime
    utc_start_dt = parse_date_to_utc(dt_str)
    start_date_str = utc_start_dt.strftime('%Y-%m-%d')

    # Fetch S2 tokens and geometries
    s2_tokens_gdf, boundary_gdf = polygon_to_s2_tokens(agstack_geoid, 19)
    print(s2_tokens_gdf)
    # Perform spatial join to filter intersecting polygons
    intersecting_s19 = gpd.sjoin(s2_tokens_gdf, boundary_gdf, how="inner", predicate="intersects")

    # Extract unique S2 indices and geometries
    L8_token_list = set(intersecting_s19['s2_index__L8'])
    L10_token_list = set(intersecting_s19['s2_index__L10'])
    L19_token_list = set(intersecting_s19['s2_index__L19'])
    geometry = intersecting_s19[['s2_index__L19', 'geometry']].drop_duplicates()


    ndvi_img_gdf = pd.DataFrame()

    # Generate valid L8 paths
    list_of_L8_paths = {
        os.path.join(pa_filePath, f's2_index__L8={x}')
        for x in L8_token_list
        if os.path.exists(os.path.join(pa_filePath, f's2_index__L8={x}'))
    }
    
    # Generate valid L10 paths
    list_of_L10_paths = {
        os.path.join(L8_path, f's2_index__L10={x}')
        for L8_path in list_of_L8_paths
        for x in L10_token_list
        if os.path.exists(os.path.join(L8_path, f's2_index__L10={x}'))
    }

    if not list_of_L10_paths:
        # Return empty GeoDataFrame and boundary_gdf if no paths are found
        return gpd.GeoDataFrame(columns=['s2_index__L19', 'NDVI', 'geometry']), boundary_gdf

    # Load weather datasets
    weather_datasets = [
        ds.dataset(path, format="parquet", partitioning="hive")
        for path in list_of_L10_paths
    ]

    for weather_dataset in weather_datasets:
        # Filter data by date
        filtered_table = weather_dataset.to_table(
            columns=['s2_index__L19', 'NDVI'],
            filter=ds.field('YY_MM_DD') == start_date_str
        ).to_pandas()

        # Filter rows with valid S2 tokens in L10
        intersected_L19_filtered_table = filtered_table[filtered_table['s2_index__L19'].isin(L19_token_list)]

        # Ensure no duplicate S2 tokens and calculate average NDVI
        aggregated_df = (
            intersected_L19_filtered_table
            .drop_duplicates(subset='s2_index__L19', keep='first')
            .groupby('s2_index__L19', as_index=False)['NDVI']
            .mean()
        )

        # Merge NDVI data with geometry based on S2 token
        merged_gdf = pd.merge(aggregated_df, geometry, on='s2_index__L19', how='inner')
        # Convert to GeoDataFrame
        ndvi_img_gdf = gpd.GeoDataFrame(merged_gdf, geometry='geometry', crs="EPSG:4326")
        ndvi_img_gdf = ndvi_img_gdf.drop(columns=['s2_index__L19'])
        
        
        # # Ensure both GeoDataFrames are in the same CRS
        ndvi_img_gdf = ndvi_img_gdf.to_crs(epsg=4326)
        boundary_gdf = boundary_gdf.to_crs(epsg=4326)

        # Clip the ndvi_gdf with the boundary_gdf
        clipped_ndvi_gdf = gpd.clip(ndvi_img_gdf, boundary_gdf)
        
        
        # print("NDVI GeoDataFrame CRS:", ndvi_img_gdf.crs)
        
        # print("Boundary GeoDataFrame CRS:", boundary_gdf.crs)
        
        # if ndvi_img_gdf.empty:
        #     print("NDVI GeoDataFrame is empty.")
        # else:
        #     print(f"NDVI GeoDataFrame has {len(ndvi_img_gdf)} records.")

        # if boundary_gdf.empty:
        #     print("Boundary GeoDataFrame is empty.")
        # else:
        #     print(f"Boundary GeoDataFrame has {len(boundary_gdf)} records.")
            
        # print("Are all NDVI geometries valid?", ndvi_img_gdf.is_valid.all())
        # print("Are all boundary geometries valid?", boundary_gdf.is_valid.all())
        
        
        # print("Geometry types in NDVI GeoDataFrame:", ndvi_img_gdf.geometry.geom_type.unique())
        # print("Geometry types in Boundary GeoDataFrame:", boundary_gdf.geometry.geom_type.unique())
        
        # intersects = gpd.sjoin(ndvi_img_gdf, boundary_gdf, how='inner', predicate='intersects')
        # print(f"Found {len(intersects)} intersecting geometries.")
        
        # boundary_gdf['geometry'] = boundary_gdf.geometry.apply(lambda x: x.simplify(0.001))
        # ndvi_img_gdf['geometry'] = ndvi_img_gdf.geometry.apply(lambda x: x.simplify(0.001))

    return clipped_ndvi_gdf, boundary_gdf


def getTMFgdfFromGeoid(agstack_geoid):
    pa_filePath_jrc = '/network/TMF_UAD/'
    # Convert date string to UTC datetime
    # utc_start_dt = parse_date_to_utc(dt_str)
    # start_date_str = utc_start_dt.strftime('%Y-%m-%d')

    # Fetch S2 tokens and geometries
    s2_tokens_gdf, boundary_gdf = polygon_to_s2_tokens_tmf(agstack_geoid, 18)
    print(s2_tokens_gdf)
    # Perform spatial join to filter intersecting polygons
    intersecting_s18 = gpd.sjoin(s2_tokens_gdf, boundary_gdf, how="inner", predicate="intersects")

    # Extract unique S2 indices and geometries
    L3_token_list = set(intersecting_s18['s2_index__L3'])
    L8_token_list = set(intersecting_s18['s2_index__L8'])
    L18_token_list = set(intersecting_s18['s2_index__L18'])
    geometry = intersecting_s18[['s2_index__L18', 'geometry']].drop_duplicates()


    ndvi_img_gdf = pd.DataFrame()

    # Generate valid L3 paths
    list_of_L3_paths = {
        os.path.join(pa_filePath_jrc, f's2_index__L3={x}')
        for x in L3_token_list
        if os.path.exists(os.path.join(pa_filePath_jrc, f's2_index__L3={x}'))
    }
    print("=====list_of_L3_paths=====",list_of_L3_paths)
    # Generate valid L8 paths
    list_of_L8_paths = {
        os.path.join(L3_path, f's2_index__L8={x}')
        for L3_path in list_of_L3_paths
        for x in L8_token_list
        if os.path.exists(os.path.join(L3_path, f's2_index__L8={x}'))
    }
    print("=========list_of_L8_paths=======",list_of_L8_paths)

    if not list_of_L8_paths:
        # Return empty GeoDataFrame and boundary_gdf if no paths are found
        return gpd.GeoDataFrame(columns=['s2_index__L18', 'TMF_UAD', 'geometry']), boundary_gdf

    # Load weather datasets
    weather_datasets = [
        ds.dataset(path, format="parquet", partitioning="hive")
        for path in list_of_L8_paths
    ]
    for weather_dataset in weather_datasets:
        # Filter data by s2_index__L18
        filtered_table = weather_dataset.to_table(
            columns=['s2_index__L18', 'TMF_UAD']
        ).to_pandas()

        # Filter rows with valid S2 tokens in L18
        intersected_L18_filtered_table = filtered_table[filtered_table['s2_index__L18'].isin(L18_token_list)]

        # Ensure no duplicate S2 tokens and calculate average TMF_UAD
        aggregated_df = (
            intersected_L18_filtered_table
            .drop_duplicates(subset='s2_index__L18', keep='first')
            .groupby('s2_index__L18', as_index=False)['TMF_UAD']
            .mean()
        )

        # Merge TMF_UAD data with geometry based on S2 token
        merged_gdf = pd.merge(aggregated_df, geometry, on='s2_index__L18', how='inner')
        # Convert to GeoDataFrame
        ndvi_img_gdf = gpd.GeoDataFrame(merged_gdf, geometry='geometry', crs="EPSG:4326")
        ndvi_img_gdf = ndvi_img_gdf.drop(columns=['s2_index__L18'])

        # Ensure both GeoDataFrames are in the same CRS
        ndvi_img_gdf = ndvi_img_gdf.to_crs(epsg=4326)
        boundary_gdf = boundary_gdf.to_crs(epsg=4326)

        # Clip the ndvi_gdf with the boundary_gdf
        clipped_ndvi_gdf = gpd.clip(ndvi_img_gdf, boundary_gdf)

    return clipped_ndvi_gdf, boundary_gdf

def getAvailableDatesForGeoID(agstack_geoid):
    base_path="/network/SENTINEL/PARQUET_NDVI_L20_3/"
    try:
        # Fetch WKT polygon and extract latitude and longitde
        wkt_polygon = fetchWKT(agstack_geoid)
        lat, lon = extractLatLonFromWKT(wkt_polygon)
    except Exception as e:
        return {"error": f"Failed to fetch WKT or extract lat/lon: {str(e)}"}

    available_dates = set()

    try:
        # Fetch S2 cell IDs for levels 8 and 10
        s2_index_L8_list, _ = get_s2_cellids_and_token_list(8, [lat], [lon])
        s2_index_L10_list, _ = get_s2_cellids_and_token_list(10, [lat], [lon])

        for s2_index_L8 in s2_index_L8_list:
            L8_path = os.path.join(base_path, f"s2_index__L8={s2_index_L8}")
            if not os.path.exists(L8_path):
                continue

            for s2_index_L10 in s2_index_L10_list:
                L10_path = os.path.join(L8_path, f"s2_index__L10={s2_index_L10}")
                if not os.path.exists(L10_path):
                    continue

                for item in os.listdir(L10_path):
                    item_path = os.path.join(L10_path, item)
                    if os.path.isdir(item_path) and re.match(r"YY_MM_DD=\d{4}-\d{2}-\d{2}", item):
                        date = item.split("=")[1]
                        available_dates.add(date)

        return {"dates": sorted(list(available_dates))}

    except Exception as e:
        return {"error": f"Failed to process geoid {agstack_geoid}: {str(e)}"}




# def getNDVIgdfFromGeoid(agstack_geoid, start_date,end_date):
#     # Fetch WKT polygon and extract lat, lon
#     wkt_polygon = fetchWKT(agstack_geoid)
#     lat, lon = extractLatLonFromWKT(wkt_polygon)

#     # Fetch S2 tokens dynamically based on geoid and coordinates
#     s2_index__L5_list, L5_cids = get_s2_cellids_and_token_list(5, [lat], [lon])
#     s2_index__L8_list, L8_cids = get_s2_cellids_and_token_list(8, [lat], [lon])
#     s2_index__L10_list, L10_cids = get_s2_cellids_and_token_list(10, [lat], [lon])
#     s2_index__L15_list, L15_cids = get_s2_cellids_and_token_list(15, [lat], [lon])
#     s2_index__L19_list, L19_cids = get_s2_cellids_and_token_list(19, [lat], [lon])
#     s2_index__L20_list, L20_cids = get_s2_cellids_and_token_list(20, [lat], [lon])

#     # Start the timer for S2 token retrieval
#     s1_time_start = time.time()

#     s1_time_end = time.time()
#     s1_time = s1_time_end - s1_time_start

#     # Fetch the NDVI datasets using the dynamically fetched S2 indices
#     s2_time_start = time.time()
#     datasets_ndvi = getLocalNDVIDatasets(s2_index__L8_list, s2_index__L10_list, pa_filePath)
    
#     print(f'datasets_ndvi--{datasets_ndvi}')
#     s2_time_end = time.time()
#     s2_time = s2_time_end - s2_time_start

#     # Get NDVI data for the requested date
#     s3_time_start = time.time()
#     ndvi_df, t_elapsed = getNDVIFromParquetForDataSetsAndDate(
#         datasets_ndvi, 
#         start_date,
#         s2_index__L10_list, 
#         s2_index__L15_list,
#         s2_index__L20_list
#     )
#     ndvi_df = ndvi_df.drop_duplicates(subset=['s2_index__L20', 'YY_MM_DD']).reset_index(drop=True)
#     s3_time_end = time.time()
#     s3_time = s3_time_end - s3_time_start
#     print('get NDVI data: ' + str(round(s3_time, 2)))

#     # Create GeoDataFrame for the S2 L20 tiles dynamically using the fetched s2_index__L20_list
#     s5_time_start = time.time()
#     # Here we would need a function that generates or fetches the correct geometry for each token
    # p20_gdf = generateGeoDataFrameFromS2Tokens(s2_index__L20_list)  # Replace with your actual method

#     # Merge NDVI data with the field GeoDataFrame
#     bound_gdf = getBoundaryGeoDataFrame(agstack_geoid)  # Adjust this function to return your boundary gdf
#     bound_gdf.set_index('s2_index__L20', drop=True, inplace=True)
#     res_gdf = bound_gdf.join(ndvi_df)
#     s5_time_end = time.time()
#     s5_time = s5_time_end - s5_time_start
#     print('merge: ' + str(round(s5_time, 2)))

#     # Interpolation to handle missing NDVI values
#     s5_time_start = time.time()
#     res_gdf.reset_index(inplace=True)
#     res2_gdf = res_gdf.copy()
#     res2_gdf['NDVI'] = res_gdf[['s2_index__L20', 'NDVI']].apply(
#         lambda x: getInterpolatedValueFromEdgeNeightbors_ndvi(x['s2_index__L20'], res2_gdf) if pd.isnull(x.NDVI) else x.NDVI,
#         axis=1
#     )
#     s5_time_end = time.time()
#     s5_time = s5_time_end - s5_time_start
#     print('interpolate: ' + str(round(s5_time, 2)))

#     # Clip the data to the boundary
#     s6_time_start = time.time()
#     resClipped_gdf = gpd.clip(res2_gdf, bound_gdf)  # Clip the results based on the boundary
#     resClipped_gdf.reset_index(drop=True, inplace=True)
#     s6_time_end = time.time()
#     s6_time = s6_time_end - s6_time_start
#     print('Clip boundaries: ' + str(round(s6_time, 2)))

#     # Final result with NDVI and boundary
#     final_res_gdf = resClipped_gdf[['s2_index__L20', 'geometry', 'NDVI']]
#     return final_res_gdf, bound_gdf



# def getNDVIgdfFromGeoid(agstack_geoid, dtStr):
    
#     configFile =  fieldDataConfigPath + agstack_geoid+'.json'
#     #NDays = (datetime.today() - datetime.strptime(dtStr,'%Y-%m-%d')).days
    
#     #Validate the date 
#     #1. Dt should be in the last 90 days from today
#     #2. Dt should be a valida date
    
#     s1_time_start = time.time()
#     with open(configFile, "r") as jsonfile:
#         fieldJSON = json.load(jsonfile)
#         field_geoid = fieldJSON['geoid']
#         field_wkt = fieldJSON['wkt']
#         p20_gdf = gpd.read_file(fieldJSON['s2_L20_index_gdfJSON'], driver='GeoJSON')
#         p19_gdf = gpd.read_file(fieldJSON['s2_L19_index_gdfJSON'], driver='GeoJSON')
#         bound_gdf =  gpd.read_file(fieldJSON['boundary_gdfJSON'], driver='GeoJSON') 
#         field_tileStr = fieldJSON['sentinel-tile-list']
#         field_s2indices = fieldJSON['indices']

#     p_gdf = p20_gdf.copy()
    
#     #get the s2 indices
#     s2_index__L8_list = json.loads(field_s2indices['s2_index__L8_list'])
#     s2_index__L10_list = json.loads(field_s2indices['s2_index__L10_list'])
#     # s2_index__L13_list = json.loads(field_s2indices['s2_index__L13_list'])
#     # s2_index__L15_list = json.loads(field_s2indices['s2_index__L15_list'])
#     s2_index__L19_list = json.loads(field_s2indices['s2_index__L19_list'])
#     s2_index__L20_list = json.loads(field_s2indices['s2_index__L20_list'])
#     s1_time_end = time.time()
#     s1_time = s1_time_end  - s1_time_start
#     #print('load field data: '+str(round(s1_time,2)))
    
#     #get the NDVI dataset
#     s2_time_start = time.time()
#     datasets_ndvi = getLocalNDVIDatasets(s2_index__L8_list, s2_index__L10_list, pa_filePath)
#     s2_time_end = time.time()
#     s2_time = s2_time_end  - s2_time_start
#     #print('load dataset directory: '+str(round(s2_time,2)))
    
#     #read the data
#     s3_time_start = time.time()
#     ndvi_df, t_elapsed = getNDVIFromParquetForDataSetsAndDate(
#         datasets_ndvi, 
#         dtStr,
#         s2_index__L10_list, 
#         s2_index__L15_list,
#         s2_index__L20_list
#     )
#     ndvi_df = ndvi_df.drop_duplicates(subset=['s2_index__L20','YY_MM_DD']).reset_index(drop=True)
#     s3_time_end = time.time()
#     s3_time = s3_time_end  - s3_time_start
#     print('get NDVI data: '+str(round(s3_time,2)))
    
#     #merge with the acct_gdf to get the gdf
#     s5_time_start = time.time()
#     ndvi_df.set_index('s2_index__L20', drop=True, inplace=True)
#     returnedDtStr = ndvi_df.YY_MM_DD.iloc[0]
#     p_gdf.set_index('s2_index__L20', drop=True, inplace=True)
#     #set the value of res_gdf
#     res_gdf = p_gdf.join(ndvi_df)
#     s5_time_end = time.time()
#     s5_time = s5_time_end  - s5_time_start
#     print('merge: '+str(round(s5_time,2)))

#     #interpolate
#     s5_time_start = time.time()
#     res_gdf.reset_index(inplace=True)
#     res2_gdf = res_gdf.copy()
#     res2_gdf['NDVI'] = res_gdf[['s2_index__L20','NDVI']].apply(lambda x: getInterpolatedValueFromEdgeNeightbors_ndvi(x['s2_index__L20'], res2_gdf) if(pd.isnull(x.NDVI)) else x.NDVI, axis = 1)
#     s5_time_end = time.time()
#     s5_time = s5_time_end  - s5_time_start
#     print('interpolate: '+str(round(s5_time,2)))
    
#     #clip
#     s6_time_start = time.time()
#     resClipped_gdf = gpd.clip(res2_gdf, bound_gdf)
#     resClipped_gdf.reset_index(drop=True, inplace=True)
#     s6_time_end = time.time()
#     s6_time = s6_time_end  - s6_time_start
#     print('Clip boundaries: '+str(round(s6_time,2)))
    
#     final_res_gdf = resClipped_gdf[['s2_index__L20','geometry','NDVI']]
#     return final_res_gdf, bound_gdf, returnedDtStr


# def get_s2_cellids_and_token_list_for_wkt(resLevel, user_fieldWKT):
#     poly = shapely.wkt.loads(user_fieldWKT)
    
#     """
#     (ln, lt) = poly.boundary.xy
#     lons=ln.tolist()
#     lats=lt.tolist()
#     """
#     lns, lts = poly.exterior.coords.xy
#     lons = lns.tolist()
#     lats = lts.tolist()

#     min_level=resLevel
#     max_level=resLevel
#     r = s2.RegionCoverer()
#     r.min_level = min_level
#     r.max_level = max_level
#     lb_lat = min(lats)
#     ub_lat = max(lats)
#     lb_lon = min(lons)
#     ub_lon = max(lons)
#     lb = s2.LatLng.from_degrees(lb_lat, lb_lon)
#     ub = s2.LatLng.from_degrees(ub_lat, ub_lon)
#     cids = r.get_covering(s2.LatLngRect.from_point_pair(lb,ub))
#     s2_token_list = []
#     for cellid in cids:
#         s2_token_list.append(cellid.to_token())

#     return s2_token_list, cids

def get_level_11_tokens(level_13_tokens):
    level_11_tokens = []
    for token in level_13_tokens:
        cell_id = s2sphere.CellId.from_token(token)
        parent_cell_id = cell_id.parent(11)  # Get parent cell at level 11
        level_11_token = parent_cell_id.to_token()
        level_11_tokens.append(level_11_token)
        print(f'level_11_tokens {level_11_tokens}')
    return level_11_tokens

def get_level_9_tokens(level_9_tokens):
    # Using list comprehension for efficiency
    level_5_tokens = [
        s2sphere.CellId.from_token(token).parent(9).to_token()
        for token in level_9_tokens
    ]
    
    # Print outside the loop to avoid multiple print calls
    print(f'level_5_tokens {level_5_tokens}')
    
    return level_9_tokens


def get_level_7_tokens(level_9_tokens):
    # Using list comprehension for efficiency
    level_7_tokens = [
        s2sphere.CellId.from_token(token).parent(7).to_token()
        for token in level_9_tokens
    ]
    
    # Print outside the loop to avoid multiple print calls
    print(f'level_7_tokens {level_7_tokens}')
    
    return level_7_tokens

def get_level_5_tokens(level_9_tokens):
    # Using list comprehension for efficiency
    level_5_tokens = [
        s2sphere.CellId.from_token(token).parent(5).to_token()
        for token in level_9_tokens
    ]
    return level_5_tokens 

# world = gpd.read_file('/home/rnaura/terrapipe/ne_110m_admin_0_countries.shp')


# last_region = None
# def determine_region(lat, lon):
#     point = Point(lon, lat)
#     # Finding the row in the dataframe where the point lies within the country's boundary
#     country = world[world.contains(point)]
#     if not country.empty:
#         # Extracting the first row as the result
#         return country.iloc[0]['SOVEREIGNT']  # or another column name if different
#     else:
#         return 'Unknown'


def empty_response():
    weather_df_filtered = pd.DataFrame({
    "latitude": [],
    "longitude": [],
    "time": [],
    "PRMSL_meansealevel": [],
    "GUST_surface": [],
    "TSOIL_0M0D1mbelowground": [],
    "SOILW_0M0D1mbelowground": [],
    "TSOIL_0D1M0D4mbelowground": [],
    "SOILW_0D1M0D4mbelowground": [],
    "TSOIL_0D4M1mbelowground": [],
    "SOILW_0D4M1mbelowground": [],
    "TSOIL_1M2mbelowground": [],
    "SOILW_1M2mbelowground": [],
    "TMP_2maboveground": [],
    "SPFH_2maboveground": [],
    "DPT_2maboveground": [],
    "RH_2maboveground": [],
    "PRATE_surface": [],
    "CRAIN_surface": [],
    "SOTYP_surface": [],
    "TCDC_entireatmosphere": [],
    "DSWRF_surface": [],
    "DLWRF_surface": [],
    "USWRF_surface": [],
    "ULWRF_surface": [],
    "PRES_maxwind": [],
    "s2_token_L5": [],
    "s2_token_L7": [],
    "s2_token_L9": [],
    "Year": [],
    "Month": [],
    "Day": [],
    "timestamp": []
})
    
    return weather_df_filtered



def meta_data_response():
    metadata = {
        "CRAIN_surface": "Categorical Rain [-]",
        "DLWRF_surface": "W/m^2",
        "DPT_2maboveground": "K",
        "DSWRF_surface": "W/m^2",
        "Day": "integer",
        "GUST_surface": "m/s",
        "Month": "integer",
        "PRATE_surface": "kg/m^2/s",
        "PRES_maxwind": "Pa",
        "PRMSL_meansealevel": "Pa",
        "RH_2maboveground": "%",
        "SOILW_0D1M0D4mbelowground": "m^3/m^3",
        "SOILW_0D4M1mbelowground": "m^3/m^3",
        "SOILW_0M0D1mbelowground": "m^3/m^3",
        "SOILW_1M2mbelowground": "m^3/m^3",
        "SOTYP_surface": "unitless",
        "SPFH_2maboveground": "kg/kg",
        "TCDC_entireatmosphere": "%",
        "TMP_2maboveground": "K",
        "TSOIL_0D1M0D4mbelowground": "K",
        "TSOIL_0D4M1mbelowground": "K",
        "TSOIL_0M0D1mbelowground": "K",
        "TSOIL_1M2mbelowground": "K",
        "ULWRF_surface": "W/m^2",
        "USWRF_surface": "W/m^2",
        "Year": "integer",
        "index": "integer",
        "latitude": "degrees",
        "longitude": "degrees",
        "time": "datetime (UTC)",
        "T_min": "K",
        "T_max": "K",
        "T_mean": "K",
        "ETo_FAO_MM": "mm",
        "ETo_HAR_MM": "mm",
        "ETo_AVG_MM": "mm",
        "ETo_FAO_IN": "inches",
        "ETo_HAR_IN": "inches",
        "ETo_AVG_IN": "inches"
    }
    
    return metadata



def meta_data_with_description():
    metadata = {
        "CRAIN_surface": "Categorical Rain",
        "DLWRF_surface": "Downward Long-Wave Rad. Flux",
        "DPT_2maboveground": "Dew Point Temperature",
        "DSWRF_surface": "Downward Short-Wave Radiation Flux",
        "Day": "Day of the month",
        "GUST_surface": "Wind Speed (Gust)",
        "Month": "Month of the year",
        "PRATE_surface": "Precipitation Rate",
        "PRES_maxwind": "Pressure at the maximum wind level",
        "PRMSL_meansealevel": "Pressure Reduced to Mean Sea Level",
        "RH_2maboveground": "Relative Humidity",
        "SOILW_0D1M0D4mbelowground": "Volumetric Soil Moisture Content at 0.1-0.4m depth",
        "SOILW_0D4M1mbelowground": "Volumetric Soil Moisture Content at 0.4-1m depth",
        "SOILW_0M0D1mbelowground": "Volumetric Soil Moisture Content at 0-0.1m depth",
        "SOILW_1M2mbelowground": "Volumetric Soil Moisture Content at 1-2m depth",
        "SOTYP_surface": "Soil Type",
        "SPFH_2maboveground": "Specific Humidity",
        "TCDC_entireatmosphere": "Total Cloud Cover",
        "TMP_2maboveground": "Temperature at 2 meters above ground",
        "TSOIL_0D1M0D4mbelowground": "Soil Temperature at 0.1-0.4m depth",
        "TSOIL_0D4M1mbelowground": "Soil Temperature at 0.4-1m depth",
        "TSOIL_0M0D1mbelowground": "Soil Temperature at 0-0.1m depth",
        "TSOIL_1M2mbelowground": "Soil Temperature at 1-2m depth",
        "ULWRF_surface": "Upward Long-Wave Rad. Flux",
        "USWRF_surface": "Upward Short-Wave Radiation Flux",
        "Year": "Year of the record",
        "index": "Index of the record",
        "latitude": "Latitude",
        "longitude": "Longitude",
        "time": "Time of the record in UTC format",
        "T_min": "Minimum Temperature (2m above ground)",
        "T_max": "Maximum Temperature (2m above ground)",
        "T_mean": "Mean Temperature (2m above ground)",
        "ETo_FAO_MM": "FAO-based Evapotranspiration in mm",
        "ETo_HAR_MM": "Hargreaves-based Evapotranspiration in mm",
        "ETo_AVG_MM": "Average Evapotranspiration in mm",
        "ETo_FAO_IN": "FAO-based Evapotranspiration in inches",
        "ETo_HAR_IN": "Hargreaves-based Evapotranspiration in inches",
        "ETo_AVG_IN": "Average Evapotranspiration in inches"
    }
    
    return metadata




def get_noaa_metadata():
    metadata = {
        "ETo_FAO_IN": "inches",
        "ETo_HAR_IN": "inches",
        "ETo_AVG_IN": "inches",
        "ETo_FAO_MM": "millimetre",
        "ETo_HAR_MM": "millimetre",
        "ETo_AVG_MM": "millimetre",
        "TS": "date (UTC)",
        "T_max": "Celsius",
        "T_mean": "Celsius",
        "T_min": "Celsius",
        "lat": "degrees",
        "lon": "degrees",
    }

    return metadata


def get_noaa_forecast_metadata():
    metadata = {
        "Day": "day of the month (integer)",
        "ETo_FAO_inches": "inches",
        "ETo_HAR_inches": "inches",
        "ETo_average_inches": "inches",
        "Month": "month of the year (integer)",
        "RH_max": "percent (%)",
        "RH_mean": "percent (%)",
        "RH_min": "percent (%)",
        "R_n__MJpm2": "megajoules per square meter (MJ/m^2)",
        "R_s__MJpm2": "megajoules per square meter (MJ/m^2)",
        "TS": "timestamp (UTC)",
        "T_dew": "Celsius",
        "T_max": "Celsius",
        "T_mean": "Celsius",
        "T_min": "Celsius",
        "U_z": "meters per second (m/s)",
        "Year": "year (integer)",
        "index": "row index (integer)",
        "lat": "degrees",
        "lon": "degrees",
        "z_msl": "meters above mean sea level (m)",
    }

    return metadata


def get_noaa_description():
    
    metadata_noaa_description = {
        "ETo_FAO_IN": "Evapotranspiration using FAO method (inches)",
        "ETo_HAR_IN": "Evapotranspiration using Hargreaves method (inches)",
        "ETo_AVG_IN": "Average evapotranspiration (inches)",
        "ETo_FAO_MM": "Evapotranspiration using FAO method (millimetres)",
        "ETo_HAR_MM": "Evapotranspiration using Hargreaves method (millimetres)",
        "ETo_AVG_MM": "Average evapotranspiration (millimetres)",
        "TS": "Timestamp in UTC format",
        "T_max": "Maximum daily air temperature",
        "T_mean": "Mean daily air temperature",
        "T_min": "Minimum daily air temperature",
        "lat": "Latitude",
        "lon": "Longitude"
        }


    return metadata_noaa_description


def get_noaa_forecast_description():
    metadata_noaa_description = {
        "Day": "Day of the month (1-31)",
        "ETo_FAO_inches": "Evapotranspiration using FAO method (inches)",
        "ETo_HAR_inches": "Evapotranspiration using Hargreaves method (inches)",
        "ETo_average_inches": "Average daily evapotranspiration (inches)",
        "Month": "Month of the year (1-12)",
        "RH_max": "Maximum relative humidity (%)",
        "RH_mean": "Mean relative humidity (%)",
        "RH_min": "Minimum relative humidity (%)",
        "R_n__MJpm2": "Net radiation at the crop surface (MJ/m²)",
        "R_s__MJpm2": "Solar radiation (MJ/m²)",
        "TS": "Timestamp in UTC format (ISO 8601)",
        "T_dew": "Dew point temperature (°C)",
        "T_max": "Maximum daily air temperature (°C)",
        "T_mean": "Mean daily air temperature (°C)",
        "T_min": "Minimum daily air temperature (°C)",
        "U_z": "Wind speed at 2 meters (m/s)",
        "Year": "Year (e.g., 2024)",
        "index": "Row index (used for data processing)",
        "lat": "Latitude in decimal degrees",
        "lon": "Longitude in decimal degrees",
        "z_msl": "Elevation above mean sea level (m)"
    }

    return metadata_noaa_description




def get_GHCND_metdata():
    metadata = {
        "AWDR": "Degrees (°)",
        "AWND": "meters per second",
        "Day": "integer",
        "ETo_AVG_IN": "inches",
        "ETo_AVG_MM": "millimeters",
        "ETo_FAO_IN": "inches",
        "ETo_FAO_MM": "millimeters",
        "ETo_HAR_IN": "inches",
        "ETo_HAR_MM": "millimeters",
        "Month": "integer",
        "PRCP": "inches",
        "SNOW": "Millimeters (mm)",
        "SNWD": "Millimeters (mm)",
        "TAVG": "Celsius",
        "TMAX": "Celsius",
        "TMIN": "Celsius",
        "TOBS": "Celsius",
        "WDF2": "degrees",
        "WDF5": "degrees",
        "WSF2": "km/h",
        "WSF5": "km/h",
        "WESD": "Millimeters (mm)",
        "WSFI": "Miles per hour (mph)",  # Added
        "Year": "integer",
        "aifstime_utc": "datetime (UTC)",
        "elev": "meters",
        "index": "integer",
        "lat": "degrees",
        "lon": "degrees"
    }
    return metadata


def get_GHCND_description():
    metadata_global_description = {
        "AWDR": "Average wind direction over a specified time period",
        "AWND": "Average daily wind speed",
        "Day": "Day of the month",
        "ETo_AVG_IN": "Average evapotranspiration",
        "ETo_AVG_MM": "Average evapotranspiration",
        "ETo_FAO_IN": "Evapotranspiration using FAO method",
        "ETo_FAO_MM": "Evapotranspiration using FAO method",
        "ETo_HAR_IN": "Evapotranspiration using Hargreaves method",
        "ETo_HAR_MM": "Evapotranspiration using Hargreaves method",
        "Month": "Month of observation",
        "PRCP": "Precipitation",
        "SNOW": "The total snowfall accumulated over a specific time period, usually measured in a 24-hour interval.",
        "SNWD": "The snow depth at observation time, measured in millimeters (mm).",
        "TAVG": "Average temperature",
        "TMAX": "Maximum temperature",
        "TMIN": "Minimum temperature",
        "TOBS": "Temperature observed at a specific time, typically in Celsius.",
        "WDF2": "Direction of fastest 2-minute wind",
        "WDF5": "Direction of fastest 5-minute wind",
        "WSF2": "Fastest 2-minute wind speed",
        "WSF5": "Fastest 5-minute wind speed",
        "WESD": "Water equivalent of snow depth, measured in millimeters (mm).",
        "WSFI": "Fastest instantaneous wind speed, typically recorded in miles per hour (mph).",  # Added
        "Year": "Year of observation",
        "aifstime_utc": "UTC timestamp",
        "elev": "Elevation",
        "index": "Index or row number",
        "lat": "Latitude",
        "lon": "Longitude"
    }

    return metadata_global_description




def get_aus_metadata():
    metadata = {
        "Day": "integer",
        "ETo_AVG_IN": "inches",
        "ETo_FAO_IN": "inches",
        "ETo_HAR_IN": "inches",
        "Month": "integer",
        "RH_mean": "%",
        "T_dew": "Celsius",
        "T_max": "Celsius",
        "T_mean": "Celsius",
        "T_min": "Celsius",
        "U_z": "m/s",
        "Year": "integer",
        "air_temp": "Celsius",
        "apparent_t": "Celsius",
        "cloud": "string",
        "cloud_base_m": "meters",
        "cloud_oktas": "oktas",
        "cloud_type": "string",
        "cloud_type_id": "integer",
        "delta_t": "Celsius",
        "dewpt": "Celsius",
        "gust_kmh": "km/h",
        "gust_kt": "knots",
        "aifstime_utc":"datetime (UTC)",
        "history_product": "string",
        "index": "integer",
        "lat": "degrees",
        "local_date_time": "string",
        "local_date_time_full": "string",
        "lon": "degrees",
        "name": "string",
        "press": "hPa",
        "press_msl": "hPa",
        "press_qnh": "hPa",
        "press_tend": "string",
        "rel_hum": "%",
        "sea_state": "string",
        "sort_order": "float",
        "swell_dir_worded": "string",
        "swell_height": "float",
        "swell_period": "string",
        "vis_km": "km",
        "weather": "string",
        "wind_dir": "string",
        "wind_spd_kmh": "km/h",
        "wind_spd_kt": "knots",
        "wmo": "float"
    }
    return metadata

def get_aus_description():
    metadata_aus_description = {
        "Day": "Day of the month",
        "ETo_AVG_IN": "Average evapotranspiration",
        "ETo_FAO_IN": "Evapotranspiration using FAO method",
        "ETo_HAR_IN": "Evapotranspiration using Hargreaves method",
        "Month": "Month of the year",
        "RH_mean": "Mean relative humidity",
        "T_dew": "Dew point temperature",
        "T_max": "Maximum temperature of the day",
        "T_mean": "Mean temperature of the day",
        "T_min": "Minimum temperature of the day",
        "U_z": "Wind speed at reference height",
        "Year": "Year of the observation",
        "air_temp": "Air temperature",
        "apparent_t": "Apparent temperature",
        "cloud": "Cloud cover description",
        "cloud_base_m": "Cloud base height",
        "cloud_oktas": "Cloud cover in oktas",
        "cloud_type": "Type of clouds",
        "cloud_type_id": "Cloud type identifier",
        "delta_t": "Temperature difference between dry and wet bulbs",
        "dewpt": "Dew point temperature",
        "gust_kmh": "Wind gust speed",
        "gust_kt": "Wind gust speed",
        "history_product": "Historical weather product",
        "index": "Unique identifier for the record",
        "lat": "Latitude of the location",
        "local_date_time": "Local date and time",
        "local_date_time_full": "Full local date and time in format YYYYMMDDHHMMSS",
        "lon": "Longitude of the location",
        "aifstime_utc":"Timestamp in UTC format",
        "name": "Station or location name",
        "press": "Air pressure",
        "press_msl": "Mean sea level air pressure",
        "press_qnh": "Air pressure at the station",
        "press_tend": "Pressure tendency",
        "rel_hum": "Relative humidity",
        "sea_state": "State of the sea",
        "sort_order": "Sort order index",
        "swell_dir_worded": "Swell direction",
        "swell_height": "Swell height",
        "swell_period": "Swell period",
        "vis_km": "Visibility distance",
        "weather": "Weather description",
        "wind_dir": "Wind direction",
        "wind_spd_kmh": "Wind speed",
        "wind_spd_kt": "Wind speed",
        "wmo": "World Meteorological Organization station identifier"
    }
    return metadata_aus_description


def get_aus_forcasted_metadata():
    metadata = {
        "Day": "integer",
        "ETo_FAO_inches": "inches",
        "ETo_HAR_inches": "inches",
        "ETo_average_inches": "inches",
        "Month": "integer",
        "R_n__MJpm2": "MJ/m^2",
        "R_s__MJpm2": "MJ/m^2",
        "TS": "datetime (UTC)",
        "T_max": "Celsius",
        "T_mean": "Celsius",
        "T_min": "Celsius",
        "Year": "integer",
        "index": "integer",
        "lat": "degrees",
        "lon": "degrees",
        "z_msl": "meters"
    }
    return metadata

def get_aus_forcasted_description():
    metadata_aus_forecasted_description = {
        "Day": "Day of the month",
        "ETo_FAO_inches": "Evapotranspiration using FAO method",
        "ETo_HAR_inches": "Evapotranspiration using Hargreaves method",
        "ETo_average_inches": "Average evapotranspiration",
        "Month": "Month of the year",
        "R_n__MJpm2": "Net radiation at crop surface",
        "R_s__MJpm2": "Solar radiation at crop surface",
        "TS": "Timestamp in UTC format",
        "T_max": "Maximum temperature of the day",
        "T_mean": "Mean temperature of the day",
        "T_min": "Minimum temperature of the day",
        "Year": "Year of the observation",
        "index": "Unique identifier for the record",
        "lat": "Latitude of the station",
        "lon": "Longitude of the station",
        "z_msl": "Elevation above mean sea level"
    }
    return metadata_aus_forecasted_description


def satelite_metadata():
    metadata = {
    "NDVI": "dimensionless",
    "UTC_DATETIME": "datetime (UTC)",
    "s2_index__L19": "string (identifier)"
    }

    return metadata

def get_sentinel_data_description():
    metadata_sentinel_description = {
        "NDVI": "Normalized Difference Vegetation Index, indicating vegetation health",
        "UTC_DATETIME": "Date and time of data capture in UTC format",
        "s2_index__L19": "Unique identifier for the Sentinel-2 data product at level 19"
    }
    return metadata_sentinel_description


def eto_metadata():
    metadata = {
        "average_eto": "float (inches/day)",
        "date": "datetime (UTC)"
    }
    return metadata

def get_eto_data_description():
    metadata_eto_description = {
        "average_eto": "Average reference evapotranspiration (ETo) in inches per day, representing the rate of water loss through evaporation and transpiration from a reference surface.",
        "date":"Date and time of data capture in UTC format"
    }
    return metadata_eto_description


def get_satelite_stats_metadata():
<<<<<<< HEAD

    metadata = {
            "10th_prctile": "dimensionless",
            "25th_prctile": "dimensionless",
            "5th_prctile": "dimensionless",
            "75th_prctile": "dimensionless",
            "90th_prctile": "dimensionless",
            "95th_prctile": "dimensionless",
            "date": "datetime (GMT)",
            "avg": "dimensionless",
            "max": "dimensionless",
            "median": "dimensionless",
            "min": "dimensionless",
            "stddev": "dimensionless"
        }
=======
>>>>>>> Rajat

    metadata = {
            "10th_prctile": "dimensionless",
            "25th_prctile": "dimensionless",
            "5th_prctile": "dimensionless",
            "75th_prctile": "dimensionless",
            "90th_prctile": "dimensionless",
            "95th_prctile": "dimensionless",
            "date": "datetime (GMT)",
            "avg": "dimensionless",
            "max": "dimensionless",
            "median": "dimensionless",
            "min": "dimensionless",
            "stddev": "dimensionless"
        }

<<<<<<< HEAD
    return metadata


def get_satelite_desc():
    
    satelite_metadata_descriptions = {
        "10th_prctile": "10th percentile value of NDVI data",
        "25th_prctile": "25th percentile value of NDVI data",
        "5th_prctile": "5th percentile value of NDVI data",
        "75th_prctile": "75th percentile value of NDVI data",
        "90th_prctile": "90th percentile value of NDVI data",
        "95th_prctile": "95th percentile value of NDVI data",
        "date": "Date of data observation in GMT",
        "avg": "Average value of NDVI data",
        "max": "Maximum value of NDVI data",
        "median": "Median value of NDVI data",
        "min": "Minimum value of NDVI data",
        "stddev": "Standard deviation of NDVI data"
    }


    return satelite_metadata_descriptions

#cimis 

def get_cimis_metadata():
    metadata = {
        "Date": "datetime (UTC)",
        "Day": "integer",
        "HlyAirTmp": "degrees Fahrenheit",
        "HlyAsceEto": "inches",
        "HlyAsceEtr": "inches",
        "HlyDewPnt": "degrees Fahrenheit",
        "HlyEto": "inches",
        "HlyNetRad": "watts per square meter",
        "HlyPrecip": "inches",
        "HlyRelHum": "percentage",
        "HlyResWind": "miles per hour",
        "HlySoilTmp": "degrees Fahrenheit",
        "HlySolRad": "watts per square meter",
        "HlyVapPres": "hectopascals",
        "HlyWindDir": "degrees",
        "HlyWindSpd": "miles per hour",
        "Month": "integer",
        "Year": "integer",
        "elevation": "feet",
        "lat": "degrees",
        "lon": "degrees"
    }

    return metadata

def get_cimis_description():
    metadata_cimis_description = {
        "Date": "Timestamp in UTC format",
        "Day": "Day of the month",
        "HlyAirTmp": "Hourly air temperature",
        "HlyAsceEto": "Hourly evapotranspiration using ASCE method",
        "HlyAsceEtr": "Hourly reference evapotranspiration using ASCE method",
        "HlyDewPnt": "Hourly dew point temperature",
        "HlyEto": "Hourly evapotranspiration",
        "HlyNetRad": "Hourly net radiation",
        "HlyPrecip": "Hourly precipitation",
        "HlyRelHum": "Hourly relative humidity",
        "HlyResWind": "Hourly wind speed at reference height",
        "HlySoilTmp": "Hourly soil temperature",
        "HlySolRad": "Hourly solar radiation",
        "HlyVapPres": "Hourly vapor pressure",
        "HlyWindDir": "Hourly wind direction",
        "HlyWindSpd": "Hourly wind speed",
        "Month": "Month of the year",
        "Year": "Year of observation",
        "elevation": "Elevation of the location",
        "lat": "Latitude of the location",
        "lon": "Longitude of the location"
    }
    return metadata_cimis_description

def read_parquet_file(file):
    return pq.read_table(file)


def api_ticket_response():
    return {
        "NCEP": {
            "id": 1,
            "type": "Weather",
            "description": "https://www.ncei.noaa.gov/",
            "coverage": "Global",
            "frequency": "Hourly",
            "authors": "US National Centers for Environmental Prediction (NCEP)",
            "company": "US NCEP-NCEI",
            "create_date": "2024-10-19",
            "picture": "https://images.app.goo.gl/ytEe5jAFyQGj6rfC6",
            "dev_stage": "alpha",
            "price_scope_month": 1,
            "price_units": "USD",
            "root": "getNCEP",
            "meta_data": "{\"parameter_name\": \"star_date\", \"parameter_type\": \"Datetime\", \"parameter_example\": \"start_date=2022-10-10\", \"parameter_required\": True},{\"parameter_name\": \"end_date\", \"parameter_type\": \"Datetime\", \"parameter_example\": \"end_date=2022-10-20\", \"parameter_required\": False}, {\"parameter_name\": \"geoid\", \"parameter_type\": \"String\", \"parameter_example\": \"geoid=f2efb5783281811923424bd80fddbd04beb7ba36c22433031530c34911a3832f\", \"parameter_required\": True}"
        },
        "NLDAS": {
            "id": 3,
            "type": "Weather",
            "description": "https://www.weather.gov/",
            "coverage": "United States",
            "frequency": "Daily",
            "authors": "US National Weather Service",
            "company": "NASA NOAA",
            "create_date": "2024-10-19",
            "picture": "https://images.app.goo.gl/kN6xhQqhfqhHEENh9",
            "dev_stage": "alpha",
            "price_scope_month": 1,
            "price_units": "USD",
            "root": "getNLDAS",
            "meta_data": "{\"parameter_name\": \"star_date\", \"parameter_type\": \"Datetime\", \"parameter_example\": \"start_date=2022-10-10\", \"parameter_required\": True},{\"parameter_name\": \"end_date\", \"parameter_type\": \"Datetime\", \"parameter_example\": \"end_date=2022-10-20\", \"parameter_required\": False}, {\"parameter_name\": \"geoid\", \"parameter_type\": \"String\", \"parameter_example\": \"geoid=f2efb5783281811923424bd80fddbd04beb7ba36c22433031530c34911a3832f\", \"parameter_required\": True}"
        },
        "NOAA": {
            "id": 3,
            "type": "Weather",
            "description": "https://www.weather.gov/",
            "coverage": "United States",
            "frequency": "Daily",
            "authors": "US National Weather Service",
            "company": "NASA NOAA",
            "create_date": "2024-10-19",
            "picture": "https://images.app.goo.gl/kN6xhQqhfqhHEENh9",
            "dev_stage": "alpha",
            "price_scope_month": 1,
            "price_units": "USD",
            "root": "getNOAA",
            "meta_data": "{\"parameter_name\": \"star_date\", \"parameter_type\": \"Datetime\", \"parameter_example\": \"start_date=2022-10-10\", \"parameter_required\": True},{\"parameter_name\": \"end_date\", \"parameter_type\": \"Datetime\", \"parameter_example\": \"end_date=2022-10-20\", \"parameter_required\": False}, {\"parameter_name\": \"geoid\", \"parameter_type\": \"String\", \"parameter_example\": \"geoid=f2efb5783281811923424bd80fddbd04beb7ba36c22433031530c34911a3832f\", \"parameter_required\": True}"
        },
        "NOAAForecast" :{
            "id": 5,
            "type": "Weather",
            "description": "https://www.ngdc.noaa.gov/",
            "coverage": "Global",
            "frequency": "Daily",
            "authors": "NOAA",
            "company": "",
            "create_date": "2024-10-19",
            "picture": "https://images.app.goo.gl/Tu7GTZqGQq7UmSwa8",
            "dev_stage": "alpha",
            "price_scope_month": 1,
            "price_units": "USD",
            "root": "getNOAAForecast",
            "meta_data": "{\"parameter_name\": \"star_date\", \"parameter_type\": \"Datetime\", \"parameter_example\": \"start_date=2022-10-10\", \"parameter_required\": True},{\"parameter_name\": \"end_date\", \"parameter_type\": \"Datetime\", \"parameter_example\": \"end_date=2022-10-20\", \"parameter_required\": False}, {\"parameter_name\": \"geoid\", \"parameter_type\": \"String\", \"parameter_example\": \"geoid=f2efb5783281811923424bd80fddbd04beb7ba36c22433031530c34911a3832f\", \"parameter_required\": True}"
        },
        
        "GHCND": {
            "id": 5,
            "type": "Weather",
            "description": "https://www.ngdc.noaa.gov/",
            "coverage": "Global",
            "frequency": "Daily",
            "authors": "NOAA",
            "company": "",
            "create_date": "2024-10-19",
            "picture": "https://images.app.goo.gl/Tu7GTZqGQq7UmSwa8",
            "dev_stage": "alpha",
            "price_scope_month": 1,
            "price_units": "USD",
            "root": "getGHCND",
            "meta_data": "{\"parameter_name\": \"star_date\", \"parameter_type\": \"Datetime\", \"parameter_example\": \"start_date=2022-10-10\", \"parameter_required\": True},{\"parameter_name\": \"end_date\", \"parameter_type\": \"Datetime\", \"parameter_example\": \"end_date=2022-10-20\", \"parameter_required\": False}, {\"parameter_name\": \"geoid\", \"parameter_type\": \"String\", \"parameter_example\": \"geoid=f2efb5783281811923424bd80fddbd04beb7ba36c22433031530c34911a3832f\", \"parameter_required\": True}"
        
        }
    }


def getWeatherFromNLDAS(agstack_geoid, start_date,end_date):
    filePath = '/network/NLDAS/'

    tok = start_date.split('-')
    YYYY_str = tok[0]
    MM_str = tok[1]
    DD_str = tok[2]

    # Create a datetime object and localize it to UTC
    local_dt = datetime(int(YYYY_str), int(MM_str), int(DD_str))
    utc_dt = pytz.utc.localize(local_dt)  # Convert to UTC

    if end_date:
        # Parse end date and localize it to UTC if provided
        end_tok = end_date.split('-')
        end_YYYY_str = int(end_tok[0])
        end_MM_str = int(end_tok[1])
        end_DD_str = int(end_tok[2])
        end_local_dt = datetime(end_YYYY_str, end_MM_str, end_DD_str)
        utc_end_dt = pytz.utc.localize(end_local_dt)

=======

    return metadata


def get_satelite_desc():
    
    satelite_metadata_descriptions = {
        "10th_prctile": "10th percentile value of NDVI data",
        "25th_prctile": "25th percentile value of NDVI data",
        "5th_prctile": "5th percentile value of NDVI data",
        "75th_prctile": "75th percentile value of NDVI data",
        "90th_prctile": "90th percentile value of NDVI data",
        "95th_prctile": "95th percentile value of NDVI data",
        "date": "Date of data observation in GMT",
        "avg": "Average value of NDVI data",
        "max": "Maximum value of NDVI data",
        "median": "Median value of NDVI data",
        "min": "Minimum value of NDVI data",
        "stddev": "Standard deviation of NDVI data"
    }


    return satelite_metadata_descriptions

#cimis 

def get_cimis_metadata():
    metadata = {
        "Date": "datetime (UTC)",
        "Day": "integer",
        "HlyAirTmp": "degrees Fahrenheit",
        "HlyAsceEto": "inches",
        "HlyAsceEtr": "inches",
        "HlyDewPnt": "degrees Fahrenheit",
        "HlyEto": "inches",
        "HlyNetRad": "watts per square meter",
        "HlyPrecip": "inches",
        "HlyRelHum": "percentage",
        "HlyResWind": "miles per hour",
        "HlySoilTmp": "degrees Fahrenheit",
        "HlySolRad": "watts per square meter",
        "HlyVapPres": "hectopascals",
        "HlyWindDir": "degrees",
        "HlyWindSpd": "miles per hour",
        "Month": "integer",
        "Year": "integer",
        "elevation": "feet",
        "lat": "degrees",
        "lon": "degrees"
    }

    return metadata

def get_cimis_description():
    metadata_cimis_description = {
        "Date": "Timestamp in UTC format",
        "Day": "Day of the month",
        "HlyAirTmp": "Hourly air temperature",
        "HlyAsceEto": "Hourly evapotranspiration using ASCE method",
        "HlyAsceEtr": "Hourly reference evapotranspiration using ASCE method",
        "HlyDewPnt": "Hourly dew point temperature",
        "HlyEto": "Hourly evapotranspiration",
        "HlyNetRad": "Hourly net radiation",
        "HlyPrecip": "Hourly precipitation",
        "HlyRelHum": "Hourly relative humidity",
        "HlyResWind": "Hourly wind speed at reference height",
        "HlySoilTmp": "Hourly soil temperature",
        "HlySolRad": "Hourly solar radiation",
        "HlyVapPres": "Hourly vapor pressure",
        "HlyWindDir": "Hourly wind direction",
        "HlyWindSpd": "Hourly wind speed",
        "Month": "Month of the year",
        "Year": "Year of observation",
        "elevation": "Elevation of the location",
        "lat": "Latitude of the location",
        "lon": "Longitude of the location"
    }
    return metadata_cimis_description

def read_parquet_file(file):
    return pq.read_table(file)


def api_ticket_response():
    return {
        "NCEP": {
            "id": 1,
            "type": "Weather",
            "description": "https://www.ncei.noaa.gov/",
            "coverage": "Global",
            "frequency": "Hourly",
            "authors": "US National Centers for Environmental Prediction (NCEP)",
            "company": "US NCEP-NCEI",
            "create_date": "2024-10-19",
            "picture": "https://images.app.goo.gl/ytEe5jAFyQGj6rfC6",
            "dev_stage": "alpha",
            "price_scope_month": 1,
            "price_units": "USD",
            "root": "getNCEP",
            "meta_data": "{\"parameter_name\": \"star_date\", \"parameter_type\": \"Datetime\", \"parameter_example\": \"start_date=2022-10-10\", \"parameter_required\": True},{\"parameter_name\": \"end_date\", \"parameter_type\": \"Datetime\", \"parameter_example\": \"end_date=2022-10-20\", \"parameter_required\": False}, {\"parameter_name\": \"geoid\", \"parameter_type\": \"String\", \"parameter_example\": \"geoid=f2efb5783281811923424bd80fddbd04beb7ba36c22433031530c34911a3832f\", \"parameter_required\": True}"
        },
        "NLDAS": {
            "id": 3,
            "type": "Weather",
            "description": "https://www.weather.gov/",
            "coverage": "United States",
            "frequency": "Daily",
            "authors": "US National Weather Service",
            "company": "NASA NOAA",
            "create_date": "2024-10-19",
            "picture": "https://images.app.goo.gl/kN6xhQqhfqhHEENh9",
            "dev_stage": "alpha",
            "price_scope_month": 1,
            "price_units": "USD",
            "root": "getNLDAS",
            "meta_data": "{\"parameter_name\": \"star_date\", \"parameter_type\": \"Datetime\", \"parameter_example\": \"start_date=2022-10-10\", \"parameter_required\": True},{\"parameter_name\": \"end_date\", \"parameter_type\": \"Datetime\", \"parameter_example\": \"end_date=2022-10-20\", \"parameter_required\": False}, {\"parameter_name\": \"geoid\", \"parameter_type\": \"String\", \"parameter_example\": \"geoid=f2efb5783281811923424bd80fddbd04beb7ba36c22433031530c34911a3832f\", \"parameter_required\": True}"
        },
        "NOAA": {
            "id": 3,
            "type": "Weather",
            "description": "https://www.weather.gov/",
            "coverage": "United States",
            "frequency": "Daily",
            "authors": "US National Weather Service",
            "company": "NASA NOAA",
            "create_date": "2024-10-19",
            "picture": "https://images.app.goo.gl/kN6xhQqhfqhHEENh9",
            "dev_stage": "alpha",
            "price_scope_month": 1,
            "price_units": "USD",
            "root": "getNOAA",
            "meta_data": "{\"parameter_name\": \"star_date\", \"parameter_type\": \"Datetime\", \"parameter_example\": \"start_date=2022-10-10\", \"parameter_required\": True},{\"parameter_name\": \"end_date\", \"parameter_type\": \"Datetime\", \"parameter_example\": \"end_date=2022-10-20\", \"parameter_required\": False}, {\"parameter_name\": \"geoid\", \"parameter_type\": \"String\", \"parameter_example\": \"geoid=f2efb5783281811923424bd80fddbd04beb7ba36c22433031530c34911a3832f\", \"parameter_required\": True}"
        },
        "NOAAForecast" :{
            "id": 5,
            "type": "Weather",
            "description": "https://www.ngdc.noaa.gov/",
            "coverage": "Global",
            "frequency": "Daily",
            "authors": "NOAA",
            "company": "",
            "create_date": "2024-10-19",
            "picture": "https://images.app.goo.gl/Tu7GTZqGQq7UmSwa8",
            "dev_stage": "alpha",
            "price_scope_month": 1,
            "price_units": "USD",
            "root": "getNOAAForecast",
            "meta_data": "{\"parameter_name\": \"star_date\", \"parameter_type\": \"Datetime\", \"parameter_example\": \"start_date=2022-10-10\", \"parameter_required\": True},{\"parameter_name\": \"end_date\", \"parameter_type\": \"Datetime\", \"parameter_example\": \"end_date=2022-10-20\", \"parameter_required\": False}, {\"parameter_name\": \"geoid\", \"parameter_type\": \"String\", \"parameter_example\": \"geoid=f2efb5783281811923424bd80fddbd04beb7ba36c22433031530c34911a3832f\", \"parameter_required\": True}"
        },
        
        "GHCND": {
            "id": 5,
            "type": "Weather",
            "description": "https://www.ngdc.noaa.gov/",
            "coverage": "Global",
            "frequency": "Daily",
            "authors": "NOAA",
            "company": "",
            "create_date": "2024-10-19",
            "picture": "https://images.app.goo.gl/Tu7GTZqGQq7UmSwa8",
            "dev_stage": "alpha",
            "price_scope_month": 1,
            "price_units": "USD",
            "root": "getGHCND",
            "meta_data": "{\"parameter_name\": \"star_date\", \"parameter_type\": \"Datetime\", \"parameter_example\": \"start_date=2022-10-10\", \"parameter_required\": True},{\"parameter_name\": \"end_date\", \"parameter_type\": \"Datetime\", \"parameter_example\": \"end_date=2022-10-20\", \"parameter_required\": False}, {\"parameter_name\": \"geoid\", \"parameter_type\": \"String\", \"parameter_example\": \"geoid=f2efb5783281811923424bd80fddbd04beb7ba36c22433031530c34911a3832f\", \"parameter_required\": True}"
        
        }
    }


def getWeatherFromNLDAS(agstack_geoid, start_date,end_date):
    filePath = '/network/NLDAS/'

    tok = start_date.split('-')
    YYYY_str = tok[0]
    MM_str = tok[1]
    DD_str = tok[2]

    # Create a datetime object and localize it to UTC
    local_dt = datetime(int(YYYY_str), int(MM_str), int(DD_str))
    utc_dt = pytz.utc.localize(local_dt)  # Convert to UTC

    if end_date:
        # Parse end date and localize it to UTC if provided
        end_tok = end_date.split('-')
        end_YYYY_str = int(end_tok[0])
        end_MM_str = int(end_tok[1])
        end_DD_str = int(end_tok[2])
        end_local_dt = datetime(end_YYYY_str, end_MM_str, end_DD_str)
        utc_end_dt = pytz.utc.localize(end_local_dt)

>>>>>>> Rajat

    try:
            
        # Fetch WKT polygon and extract latitude and longitude
        wkt_polygon = fetchWKT(agstack_geoid)
        lat, lon = extractLatLonFromWKT(wkt_polygon)
        start_time = time.time()
        # Get the list of S2 indices and CIDs for the data point
        s2_index__L5_list, L5_cids = get_s2_cellids_and_token_list(5, [lat], [lon])
        print(f'token---{s2_index__L5_list}')
        
        list_of_5_paths = [filePath + 's2_token_L5=' + x for x in s2_index__L5_list
                        if os.path.exists(filePath + 's2_token_L5=' + x)]
        if not list_of_5_paths:
            return empty_response()
        
        weather_datasets = []
        for x in list_of_5_paths:
            weather_datasets.append(ds.dataset(x, format="parquet", partitioning="hive"))

        w_all = pd.DataFrame()
        w_ret = pd.DataFrame()
        for weatherDataset in weather_datasets:
            if end_date is not None:
<<<<<<< HEAD
                # Create integer representations for the start and end dates
                start_date_int = int(YYYY_str) * 10000 + int(MM_str) * 100 + int(DD_str)
                end_date_int = end_YYYY_str * 10000 + end_MM_str * 100 + end_DD_str

                # Filter between start and end date
                weather_df = weatherDataset.to_table(
                    filter=(
                        (ds.field('Year') * 10000 + ds.field('Month') * 100 + ds.field('Day')) >= start_date_int
                    ) & (
                        (ds.field('Year') * 10000 + ds.field('Month') * 100 + ds.field('Day')) <= end_date_int
                    )
                ).to_pandas()
            else:
                # Create integer representation for the start date
                start_date_int = int(YYYY_str) * 10000 + int(MM_str) * 100 + int(DD_str)

                # Filter for start date only
                weather_df = weatherDataset.to_table(
                    filter=(
                        (ds.field('Year') * 10000 + ds.field('Month') * 100 + ds.field('Day')) == start_date_int
                    )
                ).to_pandas()
=======
                # Filter between start and end date
                weather_df = weatherDataset.to_table(
                    filter=(
                        (ds.field('Year') >= int(YYYY_str)) & (ds.field('Year') <= end_YYYY_str) &
                        (ds.field('Month') >= int(MM_str)) & (ds.field('Month') <= end_MM_str) &
                        (ds.field('Day') >= int(DD_str)) & (ds.field('Day') <= end_DD_str)
                    )
                ).to_pandas()
            elif end_date is None:
                # Filter for start date only
                weather_df = weatherDataset.to_table(
                    filter=(
                        (ds.field('Year') == int(YYYY_str)) &
                        (ds.field('Month') == int(MM_str)) &
                        (ds.field('Day') == int(DD_str))
                    )
                ).to_pandas()

>>>>>>> Rajat
            # if not weather_df.empty:
            #     weather_df.dropna(inplace=True)
            w_all = pd.concat([w_all, weather_df], ignore_index=True)
        if len(w_all) > 0:
            w_all = w_all.fillna(0)
            # Convert 'time' column to datetime format if it isn't already
            w_all['time'] = pd.to_datetime(w_all['time'])

            # Ensure all numeric columns are converted to numeric type
            numeric_cols = w_all.select_dtypes(include=[np.number]).columns
            w_all[numeric_cols] = w_all[numeric_cols].apply(pd.to_numeric, errors='coerce')

            # Group by 'time' and calculate the mean for each group
            w_ret = w_all.groupby('time').mean(numeric_only=True).reset_index()
            
            # Ensure 'Year', 'Month', and 'Day' are included in the response
            w_ret['Year'] = w_ret['Year'].astype(int)
            w_ret['Month'] = w_ret['Month'].astype(int)
            w_ret['Day'] = w_ret['Day'].astype(int)
        
            w_ret['T_min'] = w_all['T_min']
            w_ret['T_max'] = w_all['T_max']
            w_ret['T_mean'] = w_all['T_mean']
            w_ret['ETo_FAO_MM'] = w_all['ETo_FAO_MM']
            w_ret['ETo_HAR_MM'] = w_all['ETo_HAR_MM']
            w_ret['ETo_AVG_MM'] = w_all['ETo_AVG_MM']
            w_ret['ETo_FAO_IN'] = w_all['ETo_HAR_IN']
            w_ret['ETo_AVG_IN'] = w_all['ETo_AVG_IN']
            # Ensure TS column is created   
    except Exception as e:
        print(e)
        w_ret = empty_response()
        
    end_time = time.time()
    time_elapsed = (end_time - start_time)

    return w_ret



# # Custom cache key function
# def make_cache_key():
#     args = request.args
#     dtStr = args.get('date')
#     agstack_geoid = args.get('geoid')
#     key_str = f"{dtStr}_{agstack_geoid}"
#     print(f'Cache key generated: {hashlib.md5(key_str.encode()).hexdigest()}')
#     return hashlib.md5(key_str.encode()).hexdigest()

# def make_key(*args, **kwargs):
#    """A function which is called to derive the key for a computed value.
#       The key in this case is the concat value of all the json request
#       parameters. Other strategy could to use any hashing function.
#    :returns: unique string for which the value should be cached.
#    """
#    user_data = request.args
#    return ",".join([f"{key}={value}" for key, value in user_data.items()])




# @cache.cached(timeout=5000000, make_cache_key=make_key)
# @memoize_ncep
# def getWeatherFromNCEP(agstack_geoid, dtStr, end_date=None):
#     filePath = '/mnt/md1/NCEP/PARQUET_S2/'
#     # Parse the date
#     tok = dtStr.split('-')
#     YYYY_str = int(tok[0])
#     MM_str = int(tok[1])
#     DD_str = int(tok[2])

#     # Create a datetime object for start date and localize it to UTC
#     local_dt = datetime(YYYY_str, MM_str, DD_str)
#     utc_dt = pytz.utc.localize(local_dt)  # Convert to UTC

#     if end_date:
#         # Parse end date and localize it to UTC if provided
#         end_tok = end_date.split('-')
#         end_YYYY_str = int(end_tok[0])
#         end_MM_str = int(end_tok[1])
#         end_DD_str = int(end_tok[2])
#         end_local_dt = datetime(end_YYYY_str, end_MM_str, end_DD_str)
#         utc_end_dt = pytz.utc.localize(end_local_dt)

    
#     w_ret = pd.DataFrame()

#     try:
        
#         # Fetch WKT polygon and extract latitude and longitude
#         wkt_polygon = fetchWKT(agstack_geoid)
#         lat, lon = extractLatLonFromWKT(wkt_polygon)
#         start_time = time.time()
        
#         list_of_L5_paths = []
#         list_of_L7_paths = []
#         list_of_L9_paths = []
        
#         # Get the list of S2 indices and CIDs for the data point
#         s2_index__L5_list, L5_cids = get_s2_cellids_and_token_list(5, [lat], [lon])

#         s2_index__L7_list, L7_cids = get_s2_cellids_and_token_list(7, [lat], [lon])
#         s2_index__L9_list, L9_cids = get_s2_cellids_and_token_list(9, [lat], [lon])
        
#         list_of_L5_paths += [filePath + 's2_token_L5=' + x for x in s2_index__L5_list
#                             if os.path.exists(filePath + 's2_token_L5=' + x)]
#         for p in list_of_L5_paths:
#             list_of_L7_paths+= [p+'/s2_token_L7='+x for x in s2_index__L7_list 
#                             if os.path.exists(p+'/s2_token_L7='+x)]
            
#         # for p in list_of_L7_paths:
#         #     list_of_L9_paths+= [p+'/s2_token_L9='+x for x in s2_index__L9_list 
#         #                     if os.path.exists(p+'/s2_token_L9='+x)]
        
#         if list_of_L7_paths == []:
#             return empty_response()

#         weather_datasets = []

#         for x in list_of_L5_paths:
#             weather_datasets.append(ds.dataset(x, format="parquet", partitioning="hive"))
            
#         w_all = pd.DataFrame()
#         for weatherDataset in weather_datasets:
#             if end_date:
#                 # Filter between start and end date
#                 weather_df = weatherDataset.to_table(
#                     filter=(
#                         (ds.field('Year') >= YYYY_str) & (ds.field('Year') <= end_YYYY_str) &
#                         (ds.field('Month') >= MM_str) & (ds.field('Month') <= end_MM_str) &
#                         (ds.field('Day') >= DD_str) & (ds.field('Day') <= end_DD_str)
#                     )
#                 ).to_pandas()
#             elif end_date is None:
#                 # Filter for start date only
#                 weather_df = weatherDataset.to_table(
#                     filter=(
#                         (ds.field('Year') == YYYY_str) &
#                         (ds.field('Month') == MM_str) &
#                         (ds.field('Day') == DD_str)
#                     )
#                 ).to_pandas()

#             if not weather_df.empty:
#                 weather_df.dropna(inplace=True)

#             # Concatenate to the main DataFrame
#             w_all = pd.concat([w_all, weather_df], ignore_index=True)
#             print(w_all)

#         if len(w_all) > 0:
#             # Convert 'time' column to datetime format if it isn't already
#             w_all['time'] = pd.to_datetime(w_all['time'])

#             # Ensure all numeric columns are converted to numeric type
#             numeric_cols = w_all.select_dtypes(include=[np.number]).columns
#             w_all[numeric_cols] = w_all[numeric_cols].apply(pd.to_numeric, errors='coerce')

#             # Group by 'time' and calculate the mean for each group
#             w_ret = w_all.groupby('time').mean(numeric_only=True).reset_index()

#             # Ensure 'Year', 'Month', and 'Day' are included in the response
#             w_ret['Year'] = w_ret['Year'].astype(int)
#             w_ret['Month'] = w_ret['Month'].astype(int)
#             w_ret['Day'] = w_ret['Day'].astype(int)
#             w_ret['T_min'] = w_all['T_min']
#             w_ret['T_max'] = w_all['T_max']
#             w_ret['T_mean'] = w_all['T_mean']
#             w_ret['ETo_FAO_MM'] = w_all['ETo_FAO_MM']
#             w_ret['ETo_HAR_MM'] = w_all['ETo_HAR_MM']
#             w_ret['ETo_AVG_MM'] = w_all['ETo_AVG_MM']
#             w_ret['ETo_FAO_IN'] = w_all['ETo_HAR_IN']
#             w_ret['ETo_AVG_IN'] = w_all['ETo_AVG_IN']
#     except:
#         w_ret = empty_response()

#     end_time = time.time()
#     time_elapsed = (end_time - start_time)

#     return w_ret

# def getWeatherFromNCEP(agstack_geoid, dtStr, end_date=None):
#     filePath = '/network/NCEP/'
    
#     # Parse the date
#     tok = dtStr.split('-')
#     YYYY_str = int(tok[0])
#     MM_str = int(tok[1])
#     DD_str = int(tok[2])

#     # Create a datetime object for start date
#     local_dt = datetime(YYYY_str, MM_str, DD_str)
#     utc_dt = pytz.utc.localize(local_dt)  # Convert to UTC

#     if end_date:
#         # Parse end date and localize it to UTC if provided
#         end_tok = end_date.split('-')
#         end_YYYY_str = int(end_tok[0])
#         end_MM_str = int(end_tok[1])
#         end_DD_str = int(end_tok[2])
#         end_local_dt = datetime(end_YYYY_str, end_MM_str, end_DD_str)
#         utc_end_dt = pytz.utc.localize(end_local_dt)


#     try:
            
#         # Fetch WKT polygon and extract latitude and longitude
#         wkt_polygon = fetchWKT(agstack_geoid)
#         lat, lon = extractLatLonFromWKT(wkt_polygon)
#         start_time = time.time()
#         # Get the list of S2 indices and CIDs for the data point
#         s2_index__L5_list, L5_cids = get_s2_cellids_and_token_list(5, [lat], [lon])
#         # s2_index__L5_list = ['4c4c']
#         print(f'token---{s2_index__L5_list}')
        
#         list_of_5_paths = [filePath + 's2_token_L5=' + x for x in s2_index__L5_list
#                         if os.path.exists(filePath + 's2_token_L5=' + x)]
#         if not list_of_5_paths:
#             return empty_response()
        
#         weather_datasets = []
#         for x in list_of_5_paths:
#             weather_datasets.append(ds.dataset(x, format="parquet", partitioning="hive"))

#         w_all = pd.DataFrame()
#         w_ret = pd.DataFrame()
#         for weatherDataset in weather_datasets:
#             if end_date is not None:
#                 # Filter between start and end date
#                 weather_df = weatherDataset.to_table(
#                     filter=(
#                         (ds.field('Year') >= int(YYYY_str)) & (ds.field('Year') <= end_YYYY_str) &
#                         (ds.field('Month') >= int(MM_str)) & (ds.field('Month') <= end_MM_str) &
#                         (ds.field('Day') >= int(DD_str)) & (ds.field('Day') <= end_DD_str)
#                     )
#                 ).to_pandas()
#             elif end_date is None:
#                 # Filter for start date only
#                 weather_df = weatherDataset.to_table(
#                     filter=(
#                         (ds.field('Year') == int(YYYY_str)) &
#                         (ds.field('Month') == int(MM_str)) &
#                         (ds.field('Day') == int(DD_str))
#                     )
#                 ).to_pandas()

#             # if not weather_df.empty:
#             #     weather_df.dropna(inplace=True)
#             w_all = pd.concat([w_all, weather_df], ignore_index=True)
#         if len(w_all) > 0:
#             w_all = w_all.fillna(0)
#             # Convert 'time' column to datetime format if it isn't already
#             w_all['time'] = pd.to_datetime(w_all['time'])

#             # Ensure all numeric columns are converted to numeric type
#             numeric_cols = w_all.select_dtypes(include=[np.number]).columns
#             w_all[numeric_cols] = w_all[numeric_cols].apply(pd.to_numeric, errors='coerce')

#             # Group by 'time' and calculate the mean for each group
#             w_ret = w_all.groupby('time').mean(numeric_only=True).reset_index()
            
#             # Ensure 'Year', 'Month', and 'Day' are included in the response
#             w_ret['Year'] = w_ret['Year'].astype(int)
#             w_ret['Month'] = w_ret['Month'].astype(int)
#             w_ret['Day'] = w_ret['Day'].astype(int)
        
#             w_ret['T_min'] = w_all['T_min']
#             w_ret['T_max'] = w_all['T_max']
#             w_ret['T_mean'] = w_all['T_mean']
#             w_ret['ETo_FAO_MM'] = w_all['ETo_FAO_MM']
#             w_ret['ETo_HAR_MM'] = w_all['ETo_HAR_MM']
#             w_ret['ETo_AVG_MM'] = w_all['ETo_AVG_MM']
#             w_ret['ETo_FAO_IN'] = w_all['ETo_HAR_IN']
#             w_ret['ETo_AVG_IN'] = w_all['ETo_AVG_IN']
#             # Ensure TS column is created   
#     except Exception as e:
#         print(e)
#         w_ret = empty_response()
        
#     end_time = time.time()
#     time_elapsed = (end_time - start_time)

#     return w_ret


def getWeatherFromNCEP(agstack_geoid, dtStr, end_date=None):
    filePath = '/network/NCEP/'
    
    # Parse the start date
    tok = dtStr.split('-')
    YYYY_str = int(tok[0])
    MM_str = int(tok[1])
    DD_str = int(tok[2])

    # Create a datetime object for the start date
    local_dt = datetime(YYYY_str, MM_str, DD_str)
    utc_dt = pytz.utc.localize(local_dt)  # Convert to UTC

    if end_date:
        # Parse end date and localize it to UTC if provided
        end_tok = end_date.split('-')
        end_YYYY_str = int(end_tok[0])
        end_MM_str = int(end_tok[1])
        end_DD_str = int(end_tok[2])
        end_local_dt = datetime(end_YYYY_str, end_MM_str, end_DD_str)
        utc_end_dt = pytz.utc.localize(end_local_dt)

    try:
        # Fetch WKT polygon and extract latitude and longitude
        wkt_polygon = fetchWKT(agstack_geoid)
        lat, lon = extractLatLonFromWKT(wkt_polygon)
        start_time = time.time()
        
        # Get the list of S2 indices and CIDs for the data point
        s2_index__L5_list, L5_cids = get_s2_cellids_and_token_list(5, [lat], [lon])
        print(f'token---{s2_index__L5_list}')
        
        list_of_5_paths = [filePath + 's2_token_L5=' + x for x in s2_index__L5_list
                           if os.path.exists(filePath + 's2_token_L5=' + x)]
        if not list_of_5_paths:
            return empty_response()
        
        weather_datasets = [ds.dataset(x, format="parquet", partitioning="hive") for x in list_of_5_paths]

        w_all = pd.DataFrame()
        for weatherDataset in weather_datasets:
            if end_date:
                weather_df = weatherDataset.to_table(
                    filter=(
                        (ds.field('Year') >= YYYY_str) & (ds.field('Year') <= end_YYYY_str) &
                        (ds.field('Month') >= MM_str) & (ds.field('Month') <= end_MM_str) &
                        (ds.field('Day') >= DD_str) & (ds.field('Day') <= end_DD_str)
                    )
                ).to_pandas()
            else:
                weather_df = weatherDataset.to_table(
                    filter=(
                        (ds.field('Year') == YYYY_str) &
                        (ds.field('Month') == MM_str) &
                        (ds.field('Day') == DD_str)
                    )
                ).to_pandas()

            if not weather_df.empty:
                print(f"Filtered rows: {len(weather_df)} for {YYYY_str}-{MM_str}-{DD_str}")
                w_all = pd.concat([w_all, weather_df], ignore_index=True)

        # Remove duplicates after concatenation
        w_all.drop_duplicates(inplace=True)

        # Process hourly or date-level aggregation if necessary
        if not w_all.empty:
            # Ensure 'time' is a datetime object
            w_all['time'] = pd.to_datetime(w_all['time'])
            w_all['Hour'] = w_all['time'].dt.hour
            w_all['date'] = pd.to_datetime(w_all[['Year', 'Month', 'Day']]).dt.date
            columns = ['s2_token_L7','s2_token_L9']
            w_all = w_all.drop(columns, axis=1)
            print(w_all)
            # Group by date and hour, converting to the expected response format
            grouped = w_all.groupby(['date', 'Hour'])
            result = {}

            for (date, hour), group in grouped:
                date_str = str(date)
                if date_str not in result:
                    result[date_str] = []
                # Add hourly data as a dictionary
                result[date_str].append(group.drop(columns=['date', 'time']).to_dict(orient='records')[0])

            end_time = time.time()
            time_elapsed = (end_time - start_time)
            print(f"Processing time: {time_elapsed} seconds")
            
        return result
    except Exception as e:
        print(e)
        return empty_response()

    

# @memoize_ncep_forecast
def fetchNCEPWeatherForecast(start_date, end_date, agstack_geoid):
    # Convert date strings to datetime objects
    start_date = pd.to_datetime(start_date)
    end_date = pd.to_datetime(end_date)
    
    # Ensure start_date is before end_date
    if start_date > end_date:
        raise ValueError("start_date must be earlier than or equal to end_date")

    global cached_weather_files_ncep, last_reload_time_ncep, cached_weather_df_ncep

    start_time_total = time.time()

    # Initialize an empty list to hold all the daily DataFrames
    all_dfs = []

    # Loop through each date in the range from start_date to end_date
    current_date = start_date
    while current_date <= end_date:
        YYYY_str = current_date.strftime('%Y')
        MM_str = current_date.strftime('%m')
        DD_str = current_date.strftime('%d')

        day_key = f'Day={DD_str}'

        # Reload data for the current day
        try:
            day_data = reload_data_ncep(YYYY_str, MM_str, DD_str)
            with concurrent.futures.ThreadPoolExecutor() as executor:
                tables = list(executor.map(read_parquet_file, day_data[day_key]))
            
            dfs = [table.to_pandas() for table in tables]
            daily_df = pd.concat(dfs)
            daily_df.dropna(inplace=True)
            daily_df.set_index(['s2_token_L9', 's2_token_L7', 's2_token_L5'], inplace=True)
            all_dfs.append(daily_df)
        except Exception as e:
            logging.error(f"Error reading Parquet files for {current_date}: {e}")

        # Move to the next day
        current_date += pd.Timedelta(days=1)

    # Combine all the daily DataFrames
    if all_dfs:
        cached_weather_df_ncep = pd.concat(all_dfs)
    else:
        cached_weather_df_ncep = pd.DataFrame()

    # Fetch WKT polygon and extract lat/lon
    wkt_polygon = fetchWKT(agstack_geoid)
    lat, lon = extractLatLonFromWKT(wkt_polygon)

    # Filter data by level 9 tokens
    s2_index__L9_list, _ = get_s2_cellids_and_token_list(9, [lat], [lon])
    s2_index__L9_list = [str(token) for token in s2_index__L9_list]

    # Initialize an empty DataFrame
    weather_df = pd.DataFrame()

    try:
        weather_df = cached_weather_df_ncep.loc[(s2_index__L9_list, slice(None), slice(None))]
    except KeyError:
        pass

    # Filter data by level 7 tokens if level 9 data is empty
    if weather_df.empty:
        level_7_tokens = get_level_7_tokens(s2_index__L9_list)
        try:
            weather_df = cached_weather_df_ncep.loc[(slice(None), level_7_tokens, slice(None))]
        except KeyError:
            pass

    # Filter data by level 5 tokens if level 7 data is empty
    if weather_df.empty:
        level_5_tokens = get_level_5_tokens(s2_index__L9_list)
        try:
            weather_df = cached_weather_df_ncep.loc[(slice(None), slice(None), level_5_tokens)]
        except KeyError:
            pass

    # Filter data by the date range
    weather_df = weather_df[
        (pd.to_datetime(weather_df['time']).dt.date >= start_date.date()) &
        (pd.to_datetime(weather_df['time']).dt.date <= end_date.date())
    ]

    print(weather_df.tail())

    # If multiple records exist, group by 'time' and calculate the mean
    if not weather_df.empty:
        weather_df['time'] = pd.to_datetime(weather_df['time'], errors='coerce')
        numeric_cols = weather_df.columns.drop('time')
        weather_df[numeric_cols] = weather_df[numeric_cols].apply(pd.to_numeric, errors='coerce')
        weather_df = weather_df.groupby('time').mean().reset_index()
        weather_df['time'] = weather_df['time'].dt.strftime('%Y-%m-%dT%H:%M:%S')

        # Convert Year, Month, Day columns to integers
        weather_df['Year'] = weather_df['Year'].astype(int)
        weather_df['Month'] = weather_df['Month'].astype(int)
        weather_df['Day'] = weather_df['Day'].astype(int)
        
    # Return empty response if no data is found
    if weather_df.empty:
        weather_df = pd.DataFrame(empty_response())

    print(f"NCEP Data for geoid {agstack_geoid} from {start_date} to {end_date}: {weather_df}")
    end_time_total = time.time()
    total_duration = end_time_total - start_time_total
    print(f'Total function execution took {total_duration:.2f} seconds')
    
    


    return weather_df



# @memoize_noaa
def getWeatherFromNOAA(agstack_geoid, start_date, end_date):
    # filePath = '/home/rajat/Downloads/rnaura_work/mnt/md0/NOAA/PARQUET_S2/'
    filePath = '/network/NOAA/PARQUET_S2/'
    # Parse the date
    tok = start_date.split('-')
    
    YYYY_str = tok[0]
    
    MM_str = tok[1]
    DD_str = tok[2]

    # Create a datetime object and localize it to UTC
    local_dt = datetime(int(YYYY_str), int(MM_str), int(DD_str))
    utc_dt = pytz.utc.localize(local_dt)  # Convert to UTC

    if end_date:
        # Parse end date and localize it to UTC if provided
        end_tok = end_date.split('-')
        end_YYYY_str = int(end_tok[0])
        end_MM_str = int(end_tok[1])
        end_DD_str = int(end_tok[2])
        end_local_dt = datetime(end_YYYY_str, end_MM_str, end_DD_str)
        utc_end_dt = pytz.utc.localize(end_local_dt)
    w_ret = pd.DataFrame()
    
<<<<<<< HEAD


    s1_time_start = time.time()
    # Fetch WKT polygon and extract latitude and longitude
    wkt_polygon = fetchWKT(agstack_geoid)
    lat, lon = extractLatLonFromWKT(wkt_polygon)
    start_time = time.time()
    # Get the list of S2 indices and CIDs for the data point
    s2_index__L3_list, L3_cids = get_s2_cellids_and_token_list(3, [lat], [lon])

    list_of_3_paths = [filePath + 's2_index__L3=' + x for x in s2_index__L3_list
            if os.path.exists(filePath + 's2_index__L3=' + x)]
    print(list_of_3_paths)
    if not list_of_3_paths:
        return empty_response()
    
    weather_datasets = []
    for x in list_of_3_paths:
        weather_datasets.append(ds.dataset(x, format="parquet", partitioning="hive"))

    w_all = pd.DataFrame()
    
    try:
        
        for weatherDataset in weather_datasets:
            if end_date is not None:
                # Combine YYYY, MM, and DD for filtering
                start_date_int = int(YYYY_str) * 10000 + int(MM_str) * 100 + int(DD_str)
                end_date_int = end_YYYY_str * 10000 + end_MM_str * 100 + end_DD_str

                weather_df = weatherDataset.to_table(
                    filter=(
                        (ds.field('YYYY') * 10000 + ds.field('MM') * 100 + ds.field('DD')) >= start_date_int
                    ) & (
                        (ds.field('YYYY') * 10000 + ds.field('MM') * 100 + ds.field('DD')) <= end_date_int
                    )
                ).to_pandas()
            else:
                start_date_int = int(YYYY_str) * 10000 + int(MM_str) * 100 + int(DD_str)

                weather_df = weatherDataset.to_table(
                    filter=(
                        (ds.field('YYYY') * 10000 + ds.field('MM') * 100 + ds.field('DD')) == start_date_int
                    )
                ).to_pandas()

            w_all = pd.concat([w_all, weather_df], ignore_index=True)

            print(w_all)

            # Ensure TS column is created
            if len(w_all) > 0:
                
                w_all = w_all.fillna(0)
                # Convert 'time' column to datetime format if it isn't already
                w_all['TS'] = pd.to_datetime(w_all['TS'])

                # Ensure all numeric columns are converted to numeric type
                numeric_cols = w_all.select_dtypes(include=[np.number]).columns
                w_all[numeric_cols] = w_all[numeric_cols].apply(pd.to_numeric, errors='coerce')

                # Group by 'time' and calculate the mean for each group
                w_ret = w_all.groupby('TS').mean(numeric_only=True).reset_index()

                # Ensure 'Year', 'Month', and 'Day' are included in the response
                # w_ret['Year'] = w_ret['YYYY'].astype(int)
                # w_ret['Month'] = w_ret['MM'].astype(int)
                # w_ret['Day'] = w_ret['DD'].astype(int)
                w_ret = w_ret.drop(columns=['YYYY','MM','DD'], axis=1)
                # weather_df['TS'] = pd.to_datetime(weather_df[['Year', 'Month', 'Day']].astype(str).agg('-'.join, axis=1)).dt.tz_localize('UTC')

                # numeric_cols = weather_df.select_dtypes(include=[np.number]).columns
                # w_df = pd.DataFrame(weather_df[numeric_cols].mean(axis=0)).T

                # # Retain the first TS value or handle TS appropriately
                # if 'TS' in weather_df.columns and not weather_df['TS'].empty:
                #     w_df['TS'] = weather_df['TS'].iloc[0]
                # else:
                #     w_df['TS'] = pd.NaT  # or some other default value

                # w_all = pd.concat([w_all, w_df], ignore_index=True)

        # if not w_all.empty:
        #     # Fill NaN values
        #     w_all = w_all.fillna(0)
        #     # Take the average of the columns
        #     w_ret = pd.DataFrame(w_all.groupby(['TS']).mean())
        #     w_ret.reset_index(inplace=True)
        # else:
        #     w_ret = pd.DataFrame()
        # except Exception as e:

        #     print(f'Error--{e}')
            # w_ret = empty_response()

        # Drop the index column
        w_ret = w_ret.reset_index(drop=True)
    except Exception as e:
        w_ret = pd.DataFrame()
    end_time = time.time()
    time_elapsed = (end_time - start_time)

    return w_ret


def getWeatherFromNOAAFORECASTED(agstack_geoid, start_date,end_date):
    filePath =  '/network/NOAA/FORECASTED/PARQUET_S2/'

    # Parse the date
    tok = start_date.split('-')
    
    YYYY_str = tok[0]
    
    MM_str = tok[1]
    DD_str = tok[2]

    # Create a datetime object and localize it to UTC
    local_dt = datetime(int(YYYY_str), int(MM_str), int(DD_str))
    utc_dt = pytz.utc.localize(local_dt)  # Convert to UTC

    if end_date:
        # Parse end date and localize it to UTC if provided
        end_tok = end_date.split('-')
        end_YYYY_str = int(end_tok[0])
        end_MM_str = int(end_tok[1])
        end_DD_str = int(end_tok[2])
        end_local_dt = datetime(end_YYYY_str, end_MM_str, end_DD_str)
        utc_end_dt = pytz.utc.localize(end_local_dt)
    w_ret = pd.DataFrame()
    try:


        s1_time_start = time.time()
        # Fetch WKT polygon and extract latitude and longitude
        wkt_polygon = fetchWKT(agstack_geoid)
        lat, lon = extractLatLonFromWKT(wkt_polygon)
        start_time = time.time()
        # Get the list of S2 indices and CIDs for the data point
        s2_index__L5_list, L5_cids = get_s2_cellids_and_token_list(5, [lat], [lon])
        print(s2_index__L5_list)
    
        list_of_5_paths = [filePath + 's2_tokens_l5=' + x for x in s2_index__L5_list
                if os.path.exists(filePath + 's2_tokens_l5=' + x)]
        
        if not list_of_5_paths:
            return empty_response()
        
        weather_datasets = []
        for x in list_of_5_paths:
            weather_datasets.append(ds.dataset(x, format="parquet", partitioning="hive"))

        w_all = pd.DataFrame()
        for weatherDataset in weather_datasets:
            if end_date is not None:
                # Filter between start and end date
            #     weather_df = weatherDataset.to_table(
            #         filter=(
            #             (ds.field('Year') >= int(YYYY_str)) & (ds.field('Year') <= end_YYYY_str) &
            #             (ds.field('Month') >= int(MM_str)) & (ds.field('Month') <= end_MM_str) &
            #             (ds.field('Day') >= int(DD_str)) & (ds.field('Day') <= end_DD_str)
            #         )
            #     ).to_pandas()
            # elif end_date is None:
            #     # Filter for start date only
            #     weather_df = weatherDataset.to_table(
            #         filter=(
            #             (ds.field('Year') == int(YYYY_str)) &
            #             (ds.field('Month') == int(MM_str)) &
            #             (ds.field('Day') == int(DD_str))
            #         )
            #     ).to_pandas()
            
            
                            # Combine Year, Month, and Day into a single datetime field and filter between start and end dates
                start_date_dt = datetime(int(YYYY_str), int(MM_str), int(DD_str))
                end_date_dt = datetime(end_YYYY_str, end_MM_str, end_DD_str)

                weather_df = weatherDataset.to_table(
                    filter=(
                        (ds.field('Year') * 10000 + ds.field('Month') * 100 + ds.field('Day')) >=
                        (start_date_dt.year * 10000 + start_date_dt.month * 100 + start_date_dt.day)
                    ) & (
                        (ds.field('Year') * 10000 + ds.field('Month') * 100 + ds.field('Day')) <=
                        (end_date_dt.year * 10000 + end_date_dt.month * 100 + end_date_dt.day)
                    )
                ).to_pandas()
            else:
                # Filter for start date only
                start_date_dt = datetime(int(YYYY_str), int(MM_str), int(DD_str))
                weather_df = weatherDataset.to_table(
                    filter=(
                        (ds.field('Year') * 10000 + ds.field('Month') * 100 + ds.field('Day')) ==
                        (start_date_dt.year * 10000 + start_date_dt.month * 100 + start_date_dt.day)
                    )
                ).to_pandas()
            w_all = pd.concat([w_all, weather_df], ignore_index=True)

            # Ensure TS column is created
            if len(w_all) > 0:
                
                w_all = w_all.fillna(0)
                # Convert 'time' column to datetime format if it isn't already
                w_all['TS'] = pd.to_datetime(w_all['TS'])

                # Ensure all numeric columns are converted to numeric type
                numeric_cols = w_all.select_dtypes(include=[np.number]).columns
                w_all[numeric_cols] = w_all[numeric_cols].apply(pd.to_numeric, errors='coerce')

                # Group by 'time' and calculate the mean for each group
                w_ret = w_all.groupby('TS').mean(numeric_only=True).reset_index()

                # Ensure 'Year', 'Month', and 'Day' are included in the response
                w_ret['Year'] = w_ret['Year'].astype(int)
                w_ret['Month'] = w_ret['Month'].astype(int)
                w_ret['Day'] = w_ret['Day'].astype(int)
                # weather_df['TS'] = pd.to_datetime(weather_df[['Year', 'Month', 'Day']].astype(str).agg('-'.join, axis=1)).dt.tz_localize('UTC')

                # numeric_cols = weather_df.select_dtypes(include=[np.number]).columns
                # w_df = pd.DataFrame(weather_df[numeric_cols].mean(axis=0)).T

                # # Retain the first TS value or handle TS appropriately
                # if 'TS' in weather_df.columns and not weather_df['TS'].empty:
                #     w_df['TS'] = weather_df['TS'].iloc[0]
                # else:
                #     w_df['TS'] = pd.NaT  # or some other default value

                # w_all = pd.concat([w_all, w_df], ignore_index=True)

        # if not w_all.empty:
        #     # Fill NaN values
        #     w_all = w_all.fillna(0)
        #     # Take the average of the columns
        #     w_ret = pd.DataFrame(w_all.groupby(['TS']).mean())
        #     w_ret.reset_index(inplace=True)
        # else:
        #     w_ret = pd.DataFrame()
    except Exception as e:

        print(f'Error--{e}')
        w_ret = empty_response()


    end_time = time.time()
    time_elapsed = (end_time - start_time)

    return w_ret

def getWeatherFromGHCND(agstack_geoid, start_date,end_date):

    filePath = '/network/GHCND/PARQUET_S2/'
    tok = start_date.split('-')
    YYYY_str = tok[0]
    MM_str = tok[1]
    DD_str = tok[2]

        # Create a datetime object and localize it to UTC
    local_dt = datetime(int(YYYY_str), int(MM_str), int(DD_str))
    utc_dt = pytz.utc.localize(datetime(int(YYYY_str), int(MM_str), int(DD_str)))

    if end_date:
        # Parse end date and localize it to UTC if provided
        end_tok = end_date.split('-')
        end_YYYY_str = int(end_tok[0])
        end_MM_str = int(end_tok[1])
        end_DD_str = int(end_tok[2])
        end_local_dt = datetime(end_YYYY_str, end_MM_str, end_DD_str)
        utc_end_dt = pytz.utc.localize(datetime(int(end_YYYY_str), int(end_MM_str), int(end_DD_str)))
    
    w_ret = pd.DataFrame()
    try:
        
        # Fetch WKT polygon and extract latitude and longitude
        wkt_polygon = fetchWKT(agstack_geoid)
        lat, lon = extractLatLonFromWKT(wkt_polygon)
        start_time = time.time()
        # Get the list of S2 indices and CIDs for the data point
        s2_index__L5_list, L5_cids = get_s2_cellids_and_token_list(5, [lat], [lon])
        print(f'token---{s2_index__L5_list}')
        
        list_of_5_paths = [filePath + 's2_token_L5=' + x for x in s2_index__L5_list
                        if os.path.exists(filePath + 's2_token_L5=' + x)]
        if not list_of_5_paths:
            return empty_response()
        
        weather_datasets = []
        for x in list_of_5_paths:
            weather_datasets.append(ds.dataset(x, format="parquet", partitioning="hive"))

        w_all = pd.DataFrame()
        
        for weatherDataset in weather_datasets:
            start_date_int = int(YYYY_str) * 10000 + int(MM_str) * 100 + int(DD_str)
            if end_date:
                end_date_int = end_YYYY_str * 10000 + end_MM_str * 100 + end_DD_str
                weather_df = weatherDataset.to_table(
                    filter=(
                        (ds.field('Year') * 10000 + ds.field('Month') * 100 + ds.field('Day')) >= start_date_int
                    ) & (
                        (ds.field('Year') * 10000 + ds.field('Month') * 100 + ds.field('Day')) <= end_date_int
                    )
                ).to_pandas()
            else:
                weather_df = weatherDataset.to_table(
                    filter=(
                        (ds.field('Year') * 10000 + ds.field('Month') * 100 + ds.field('Day')) == start_date_int
                    )
                ).to_pandas()
        

            w_all = pd.concat([w_all, weather_df], ignore_index=True)
        if len(w_all)> 0:
            w_all = w_all.fillna(0)
            # Convert 'time' column to datetime format if it isn't already
            w_all['aifstime_utc'] = pd.to_datetime(w_all['aifstime_utc'])

            # Ensure all numeric columns are converted to numeric type
            numeric_cols = w_all.select_dtypes(include=[np.number]).columns
            w_all[numeric_cols] = w_all[numeric_cols].apply(pd.to_numeric, errors='coerce')

            # Group by 'time' and calculate the mean for each group
            w_ret = w_all.groupby('aifstime_utc').mean(numeric_only=True).reset_index()
            

            # Ensure 'Year', 'Month', and 'Day' are included in the response
            w_ret['Year'] = w_ret['Year'].astype(int)
            w_ret['Month'] = w_ret['Month'].astype(int)
            w_ret['Day'] = w_ret['Day'].astype(int)
        #     if not weather_df.empty:
        #         weather_df['aifstime_utc'] = pd.to_datetime(
        #             weather_df[['Year', 'Month', 'Day']].astype(str).agg('-'.join, axis=1)
        #         ).dt.tz_localize('UTC')

        #         # Exclude 's2_token_L9' from numerical columns
        #         exclude_cols = ['s2_token_L9']
        #         numeric_cols = [
        #             col for col in weather_df.select_dtypes(include=[np.number]).columns if col not in exclude_cols
        #         ]
        #         w_df = pd.DataFrame(weather_df[numeric_cols].mean(axis=0)).T

        #         if 'aifstime_utc' in weather_df.columns and not weather_df['aifstime_utc'].empty:
        #             w_df['aifstime_utc'] = weather_df['aifstime_utc'].iloc[0]
        #         else:
        #             w_df['aifstime_utc'] = pd.NaT

        #         w_all = pd.concat([w_all, w_df], ignore_index=True)
            # print("========w_all======",w_all)

        # if not w_all.empty:
        #     w_all = w_all.fillna(0)
        #     # w_all = w_all.dropna(inplace=True)
        #     w_ret = pd.DataFrame(w_all.groupby(['aifstime_utc']).mean())
=======


    s1_time_start = time.time()
    # Fetch WKT polygon and extract latitude and longitude
    wkt_polygon = fetchWKT(agstack_geoid)
    lat, lon = extractLatLonFromWKT(wkt_polygon)
    start_time = time.time()
    # Get the list of S2 indices and CIDs for the data point
    s2_index__L3_list, L3_cids = get_s2_cellids_and_token_list(3, [lat], [lon])

    list_of_3_paths = [filePath + 's2_index__L3=' + x for x in s2_index__L3_list
            if os.path.exists(filePath + 's2_index__L3=' + x)]
    print(list_of_3_paths)
    if not list_of_3_paths:
        return empty_response()
    
    weather_datasets = []
    for x in list_of_3_paths:
        weather_datasets.append(ds.dataset(x, format="parquet", partitioning="hive"))

    w_all = pd.DataFrame()
    for weatherDataset in weather_datasets:
        if end_date is not None:
            # Filter between start and end date
            weather_df = weatherDataset.to_table(
                filter=(
                    (ds.field('YYYY') >= int(YYYY_str)) & (ds.field('YYYY') <= end_YYYY_str) &
                    (ds.field('MM') >= int(MM_str)) & (ds.field('MM') <= end_MM_str) &
                    (ds.field('DD') >= int(DD_str)) & (ds.field('DD') <= end_DD_str)
                )
            ).to_pandas()
        elif end_date is None:
            # Filter for start date only
            weather_df = weatherDataset.to_table(
                filter=(
                    (ds.field('YYYY') == int(YYYY_str)) &
                    (ds.field('MM') == int(MM_str)) &
                    (ds.field('DD') == int(DD_str))
                )
            ).to_pandas()
        w_all = pd.concat([w_all, weather_df], ignore_index=True)

        print(w_all)

        # Ensure TS column is created
        if len(w_all) > 0:
            
            w_all = w_all.fillna(0)
            # Convert 'time' column to datetime format if it isn't already
            w_all['TS'] = pd.to_datetime(w_all['TS'])

            # Ensure all numeric columns are converted to numeric type
            numeric_cols = w_all.select_dtypes(include=[np.number]).columns
            w_all[numeric_cols] = w_all[numeric_cols].apply(pd.to_numeric, errors='coerce')

            # Group by 'time' and calculate the mean for each group
            w_ret = w_all.groupby('TS').mean(numeric_only=True).reset_index()

            # Ensure 'Year', 'Month', and 'Day' are included in the response
            # w_ret['Year'] = w_ret['YYYY'].astype(int)
            # w_ret['Month'] = w_ret['MM'].astype(int)
            # w_ret['Day'] = w_ret['DD'].astype(int)
            w_ret = w_ret.drop(columns=['YYYY','MM','DD'], axis=1)
            # weather_df['TS'] = pd.to_datetime(weather_df[['Year', 'Month', 'Day']].astype(str).agg('-'.join, axis=1)).dt.tz_localize('UTC')

            # numeric_cols = weather_df.select_dtypes(include=[np.number]).columns
            # w_df = pd.DataFrame(weather_df[numeric_cols].mean(axis=0)).T

            # # Retain the first TS value or handle TS appropriately
            # if 'TS' in weather_df.columns and not weather_df['TS'].empty:
            #     w_df['TS'] = weather_df['TS'].iloc[0]
            # else:
            #     w_df['TS'] = pd.NaT  # or some other default value

            # w_all = pd.concat([w_all, w_df], ignore_index=True)

    # if not w_all.empty:
    #     # Fill NaN values
    #     w_all = w_all.fillna(0)
    #     # Take the average of the columns
    #     w_ret = pd.DataFrame(w_all.groupby(['TS']).mean())
    #     w_ret.reset_index(inplace=True)
    # else:
    #     w_ret = pd.DataFrame()
    # except Exception as e:

    #     print(f'Error--{e}')
        # w_ret = empty_response()

    # Drop the index column
    w_ret = w_ret.reset_index(drop=True)

    end_time = time.time()
    time_elapsed = (end_time - start_time)

    return w_ret


def getWeatherFromNOAAFORECASTED(agstack_geoid, start_date,end_date):
    filePath =  '/network/NOAA/FORECASTED/PARQUET_S2/'

    # Parse the date
    tok = start_date.split('-')
    
    YYYY_str = tok[0]
    
    MM_str = tok[1]
    DD_str = tok[2]

    # Create a datetime object and localize it to UTC
    local_dt = datetime(int(YYYY_str), int(MM_str), int(DD_str))
    utc_dt = pytz.utc.localize(local_dt)  # Convert to UTC

    if end_date:
        # Parse end date and localize it to UTC if provided
        end_tok = end_date.split('-')
        end_YYYY_str = int(end_tok[0])
        end_MM_str = int(end_tok[1])
        end_DD_str = int(end_tok[2])
        end_local_dt = datetime(end_YYYY_str, end_MM_str, end_DD_str)
        utc_end_dt = pytz.utc.localize(end_local_dt)
    w_ret = pd.DataFrame()
    try:


        s1_time_start = time.time()
        # Fetch WKT polygon and extract latitude and longitude
        wkt_polygon = fetchWKT(agstack_geoid)
        lat, lon = extractLatLonFromWKT(wkt_polygon)
        start_time = time.time()
        # Get the list of S2 indices and CIDs for the data point
        s2_index__L5_list, L5_cids = get_s2_cellids_and_token_list(5, [lat], [lon])
        print(s2_index__L5_list)
    
        list_of_5_paths = [filePath + 's2_tokens_l5=' + x for x in s2_index__L5_list
                if os.path.exists(filePath + 's2_tokens_l5=' + x)]
        
        if not list_of_5_paths:
            return empty_response()
        
        weather_datasets = []
        for x in list_of_5_paths:
            weather_datasets.append(ds.dataset(x, format="parquet", partitioning="hive"))

        w_all = pd.DataFrame()
        for weatherDataset in weather_datasets:
            if end_date is not None:
                # Filter between start and end date
            #     weather_df = weatherDataset.to_table(
            #         filter=(
            #             (ds.field('Year') >= int(YYYY_str)) & (ds.field('Year') <= end_YYYY_str) &
            #             (ds.field('Month') >= int(MM_str)) & (ds.field('Month') <= end_MM_str) &
            #             (ds.field('Day') >= int(DD_str)) & (ds.field('Day') <= end_DD_str)
            #         )
            #     ).to_pandas()
            # elif end_date is None:
            #     # Filter for start date only
            #     weather_df = weatherDataset.to_table(
            #         filter=(
            #             (ds.field('Year') == int(YYYY_str)) &
            #             (ds.field('Month') == int(MM_str)) &
            #             (ds.field('Day') == int(DD_str))
            #         )
            #     ).to_pandas()
            
            
                            # Combine Year, Month, and Day into a single datetime field and filter between start and end dates
                start_date_dt = datetime(int(YYYY_str), int(MM_str), int(DD_str))
                end_date_dt = datetime(end_YYYY_str, end_MM_str, end_DD_str)

                weather_df = weatherDataset.to_table(
                    filter=(
                        (ds.field('Year') * 10000 + ds.field('Month') * 100 + ds.field('Day')) >=
                        (start_date_dt.year * 10000 + start_date_dt.month * 100 + start_date_dt.day)
                    ) & (
                        (ds.field('Year') * 10000 + ds.field('Month') * 100 + ds.field('Day')) <=
                        (end_date_dt.year * 10000 + end_date_dt.month * 100 + end_date_dt.day)
                    )
                ).to_pandas()
            else:
                # Filter for start date only
                start_date_dt = datetime(int(YYYY_str), int(MM_str), int(DD_str))
                weather_df = weatherDataset.to_table(
                    filter=(
                        (ds.field('Year') * 10000 + ds.field('Month') * 100 + ds.field('Day')) ==
                        (start_date_dt.year * 10000 + start_date_dt.month * 100 + start_date_dt.day)
                    )
                ).to_pandas()
            w_all = pd.concat([w_all, weather_df], ignore_index=True)

            # Ensure TS column is created
            if len(w_all) > 0:
                
                w_all = w_all.fillna(0)
                # Convert 'time' column to datetime format if it isn't already
                w_all['TS'] = pd.to_datetime(w_all['TS'])

                # Ensure all numeric columns are converted to numeric type
                numeric_cols = w_all.select_dtypes(include=[np.number]).columns
                w_all[numeric_cols] = w_all[numeric_cols].apply(pd.to_numeric, errors='coerce')

                # Group by 'time' and calculate the mean for each group
                w_ret = w_all.groupby('TS').mean(numeric_only=True).reset_index()

                # Ensure 'Year', 'Month', and 'Day' are included in the response
                w_ret['Year'] = w_ret['Year'].astype(int)
                w_ret['Month'] = w_ret['Month'].astype(int)
                w_ret['Day'] = w_ret['Day'].astype(int)
                # weather_df['TS'] = pd.to_datetime(weather_df[['Year', 'Month', 'Day']].astype(str).agg('-'.join, axis=1)).dt.tz_localize('UTC')

                # numeric_cols = weather_df.select_dtypes(include=[np.number]).columns
                # w_df = pd.DataFrame(weather_df[numeric_cols].mean(axis=0)).T

                # # Retain the first TS value or handle TS appropriately
                # if 'TS' in weather_df.columns and not weather_df['TS'].empty:
                #     w_df['TS'] = weather_df['TS'].iloc[0]
                # else:
                #     w_df['TS'] = pd.NaT  # or some other default value

                # w_all = pd.concat([w_all, w_df], ignore_index=True)

        # if not w_all.empty:
        #     # Fill NaN values
        #     w_all = w_all.fillna(0)
        #     # Take the average of the columns
        #     w_ret = pd.DataFrame(w_all.groupby(['TS']).mean())
>>>>>>> Rajat
        #     w_ret.reset_index(inplace=True)
        # else:
        #     w_ret = pd.DataFrame()
    except Exception as e:
<<<<<<< HEAD
        print(e)
        w_ret = empty_response()

    return w_ret



# AUS
def getWeatherFromAUS(agstack_geoid, dtStr,end_date=None):
    # filePath = '/mnt/md1/NLDAS/PARQUETE_S2/'
    # filePath = '/home/rajat/Downloads/Rnaura_Work/mnt/md1/AUS/DAILY/PROCESSED/PARQUET_S2_TEST/'
    filePath = "/network/AUS/DAILY/PARQUET_S2/"

    tok = dtStr.split('-')
=======

        print(f'Error--{e}')
        w_ret = empty_response()


    end_time = time.time()
    time_elapsed = (end_time - start_time)

    return w_ret

def getWeatherFromGHCND(agstack_geoid, start_date,end_date):

    filePath = '/network/GHCND/PARQUET_S2/'
    tok = start_date.split('-')
>>>>>>> Rajat
    YYYY_str = tok[0]
    MM_str = tok[1]
    DD_str = tok[2]

<<<<<<< HEAD
    # Create a datetime object and localize it to UTC
    local_dt = datetime(int(YYYY_str), int(MM_str), int(DD_str))
    utc_dt = pytz.utc.localize(local_dt)  # Convert to UTC

    if end_date:
        # Parse end date and localize it to UTC if provided
        end_tok = end_date.split('-')
        end_YYYY_str = int(end_tok[0])
        end_MM_str = int(end_tok[1])
        end_DD_str = int(end_tok[2])
        end_local_dt = datetime(end_YYYY_str, end_MM_str, end_DD_str)
        utc_end_dt = pytz.utc.localize(end_local_dt)

    try:
        

        s1_time_start = time.time()
        # Fetch WKT polygon and extract latitude and longitude
        wkt_polygon = fetchWKT(agstack_geoid)
        lat, lon = extractLatLonFromWKT(wkt_polygon)
        start_time = time.time()
        # Get the list of S2 indices and CIDs for the data point
        s2_index__L5_list, L5_cids = get_s2_cellids_and_token_list(5, [lat], [lon])
        print(f'token---{s2_index__L5_list}')
        
        list_of_5_paths = [filePath + 's2_token_L5=' + x for x in s2_index__L5_list
                        if os.path.exists(filePath + 's2_token_L5=' + x)]
        if not list_of_5_paths:
            return empty_response()
        
        weather_datasets = []
        for x in list_of_5_paths:
            weather_datasets.append(ds.dataset(x, format="parquet", partitioning="hive"))

        w_all = pd.DataFrame()
        w_ret = pd.DataFrame()
        for weatherDataset in weather_datasets:
            if end_date is not None:
                # Create integer representations for the start and end dates
                start_date_int = int(YYYY_str) * 10000 + int(MM_str) * 100 + int(DD_str)
                end_date_int = end_YYYY_str * 10000 + end_MM_str * 100 + end_DD_str

                # Filter between start and end date
                weather_df = weatherDataset.to_table(
                    filter=(
                        (ds.field('Year') * 10000 + ds.field('Month') * 100 + ds.field('Day')) >= start_date_int
                    ) & (
                        (ds.field('Year') * 10000 + ds.field('Month') * 100 + ds.field('Day')) <= end_date_int
                    )
                ).to_pandas()
            else:
                # Create integer representation for the start date
                start_date_int = int(YYYY_str) * 10000 + int(MM_str) * 100 + int(DD_str)

                # Filter for start date only
                weather_df = weatherDataset.to_table(
                    filter=(
                        (ds.field('Year') * 10000 + ds.field('Month') * 100 + ds.field('Day')) == start_date_int
                    )
                ).to_pandas()


            # if not weather_df.empty:
            #     weather_df.dropna(inplace=True)
            print(weather_df.columns)
            w_all = pd.concat([w_all, weather_df], ignore_index=True)
        if len(w_all) > 0:
            w_all = w_all.fillna(0)
            # Convert 'time' column to datetime format if it isn't already
            # After constructing w_all DataFrame and filling Year, Month, and Day
            if 'Year' in w_all.columns and 'Month' in w_all.columns and 'Day' in w_all.columns:
                # Create 'aifstime_utc' based on Year, Month, and Day
                w_all['aifstime_utc'] = pd.to_datetime(
                    w_all[['Year', 'Month', 'Day']].astype(str).agg('-'.join, axis=1)
                ).dt.tz_localize('UTC')

            # Now you can convert the 'aifstime_utc' column
            w_all['aifstime_utc'] = pd.to_datetime(w_all['aifstime_utc'])
            # print("======w_all=====",w_all)
            # Ensure all numeric columns are converted to numeric type
            numeric_cols = w_all.select_dtypes(include=[np.number]).columns
            w_all[numeric_cols] = w_all[numeric_cols].apply(pd.to_numeric, errors='coerce')
            print("Columns in w_all before aggregation:", w_all.columns)

            # Group by 'time' and calculate the mean for each group
            w_ret = w_all.groupby('aifstime_utc').mean(numeric_only=True).reset_index()
            
            print("Non-numeric columns:", w_all.select_dtypes(exclude=[np.number]).columns)

            # print("=======w_ret========",w_ret)
=======
        # Create a datetime object and localize it to UTC
    local_dt = datetime(int(YYYY_str), int(MM_str), int(DD_str))
    utc_dt = pytz.utc.localize(datetime(int(YYYY_str), int(MM_str), int(DD_str)))

    if end_date:
        # Parse end date and localize it to UTC if provided
        end_tok = end_date.split('-')
        end_YYYY_str = int(end_tok[0])
        end_MM_str = int(end_tok[1])
        end_DD_str = int(end_tok[2])
        end_local_dt = datetime(end_YYYY_str, end_MM_str, end_DD_str)
        utc_end_dt = pytz.utc.localize(datetime(int(end_YYYY_str), int(end_MM_str), int(end_DD_str)))
    
    w_ret = pd.DataFrame()
    try:
        
        # Fetch WKT polygon and extract latitude and longitude
        wkt_polygon = fetchWKT(agstack_geoid)
        lat, lon = extractLatLonFromWKT(wkt_polygon)
        start_time = time.time()
        # Get the list of S2 indices and CIDs for the data point
        s2_index__L5_list, L5_cids = get_s2_cellids_and_token_list(5, [lat], [lon])
        print(f'token---{s2_index__L5_list}')
        
        list_of_5_paths = [filePath + 's2_token_L5=' + x for x in s2_index__L5_list
                        if os.path.exists(filePath + 's2_token_L5=' + x)]
        if not list_of_5_paths:
            return empty_response()
        
        weather_datasets = []
        for x in list_of_5_paths:
            weather_datasets.append(ds.dataset(x, format="parquet", partitioning="hive"))

        w_all = pd.DataFrame()
        for weatherDataset in weather_datasets:
            if end_date is not None:
                # Filter between start and end date
                weather_df = weatherDataset.to_table(
                    filter=(
                        (ds.field('Year') >= int(YYYY_str)) & (ds.field('Year') <= end_YYYY_str) &
                        (ds.field('Month') >= int(MM_str)) & (ds.field('Month') <= end_MM_str) &
                        (ds.field('Day') >= int(DD_str)) & (ds.field('Day') <= end_DD_str)
                    )
                ).to_pandas()
            elif end_date is None:
                # Filter for start date only
                weather_df = weatherDataset.to_table(
                    filter=(
                        (ds.field('Year') == int(YYYY_str)) &
                        (ds.field('Month') == int(MM_str)) &
                        (ds.field('Day') == int(DD_str))
                    )
                ).to_pandas()
        

            w_all = pd.concat([w_all, weather_df], ignore_index=True)
        if len(w_all)> 0:
            w_all = w_all.fillna(0)
            
            # Convert 'time' column to datetime format if it isn't already
            # w_all['aifstime_utc'] = pd.to_datetime(w_all['aifstime_utc'],utc = True)
            
            
            w_all['aifstime_utc'] = pd.to_datetime(w_all['aifstime_utc']).dt.tz_localize('UTC')

            # Ensure all numeric columns are converted to numeric type
            numeric_cols = w_all.select_dtypes(include=[np.number]).columns
            w_all[numeric_cols] = w_all[numeric_cols].apply(pd.to_numeric, errors='coerce')

            # Group by 'time' and calculate the mean for each group
            w_ret = w_all.groupby('aifstime_utc').mean(numeric_only=True).reset_index()

>>>>>>> Rajat
            # Ensure 'Year', 'Month', and 'Day' are included in the response
            w_ret['Year'] = w_ret['Year'].astype(int)
            w_ret['Month'] = w_ret['Month'].astype(int)
            w_ret['Day'] = w_ret['Day'].astype(int)
<<<<<<< HEAD
            w_ret['cloud'] = w_all['cloud']
            w_ret['cloud_type'] = w_all['cloud_type']
            w_ret['name'] = w_all['name']
            w_ret['press_tend'] = w_all['press_tend']
            # w_ret['rain_trace'] = w_all['rain_trace']
            w_ret['local_date_time'] = w_all['local_date_time']
            w_ret['local_date_time_full'] = w_all['local_date_time_full']
            w_ret['history_product'] = w_all['history_product']
            w_ret['local_date_time_full'] = w_all['local_date_time_full']
            w_ret['vis_km'] = w_all['vis_km']
            w_ret['weather'] = w_all['weather']
            w_ret['wind_dir'] = w_all['wind_dir']
            
            
            
            print(f'w_ret--{w_ret.columns}')
            
            # Ensure TS column is created
            # if not weather_df.empty:
            #     weather_df['aifstime_utc'] = pd.to_datetime(weather_df[['Year', 'Month', 'Day']].astype(str).agg('-'.join, axis=1)).dt.tz_localize('UTC')

            #     numeric_cols = weather_df.select_dtypes(include=[np.number]).columns
            #     w_df = pd.DataFrame(weather_df[numeric_cols].mean(axis=0)).T

            #     # Retain the first TS value or handle TS appropriately
            #     if 'aifstime_utc' in weather_df.columns and not weather_df['aifstime_utc'].empty:
            #         w_df['aifstime_utc'] = weather_df['aifstime_utc'].iloc[0]
            #     else:
            #         w_df['aifstime_utc'] = pd.NaT  # or some other default value

            #     w_all = pd.concat([w_all, w_df], ignore_index=True)

        # if not w_all.empty:
        #     w_all = w_all.fillna(0)
        #     # Take the average of the columns
=======
        #     if not weather_df.empty:
        #         weather_df['aifstime_utc'] = pd.to_datetime(
        #             weather_df[['Year', 'Month', 'Day']].astype(str).agg('-'.join, axis=1)
        #         ).dt.tz_localize('UTC')

        #         # Exclude 's2_token_L9' from numerical columns
        #         exclude_cols = ['s2_token_L9']
        #         numeric_cols = [
        #             col for col in weather_df.select_dtypes(include=[np.number]).columns if col not in exclude_cols
        #         ]
        #         w_df = pd.DataFrame(weather_df[numeric_cols].mean(axis=0)).T

        #         if 'aifstime_utc' in weather_df.columns and not weather_df['aifstime_utc'].empty:
        #             w_df['aifstime_utc'] = weather_df['aifstime_utc'].iloc[0]
        #         else:
        #             w_df['aifstime_utc'] = pd.NaT

        #         w_all = pd.concat([w_all, w_df], ignore_index=True)
            # print("========w_all======",w_all)

        # if not w_all.empty:
        #     w_all = w_all.fillna(0)
        #     # w_all = w_all.dropna(inplace=True)
>>>>>>> Rajat
        #     w_ret = pd.DataFrame(w_all.groupby(['aifstime_utc']).mean())
        #     w_ret.reset_index(inplace=True)
        # else:
        #     w_ret = pd.DataFrame()
    except Exception as e:
<<<<<<< HEAD
        print(e)      
        
    end_time = time.time()
    time_elapsed = (end_time - start_time)

    return w_ret


# def getEtoFromWeatherData(agstack_geoid, start_date, end_date):
#     # Define file paths
#     paths = {
#         "ncep": '/mnt/md1/NCEP/PARQUET_S2/', 
#         "nldas": '/mnt/md1/NLDAS/PARQUET_S2/',
#         "noaa": '/mnt/md1/NOAA/DAILY/PROCESSED/PARQUET_S2/',
#         "noaa_forecasts": '/mnt/md1/NOAA/FORECASTED/PARQUET_S2/',
#         "global": '/home/rnaura/mnt/md1/GHCND/DAILY/PROCESSED/PARQUET_S2/',
#         "aus": '/mnt/md1/AUS/DAILY/PROCESSED/PARQUET_S2/',
#         "aus_forecast": '/mnt/md1/AUS/DAILY/FORECASTED/PARQUET_S2/'
#     }
    
#     # Fetch WKT and extract latitude and longitude
#     wkt_polygon = fetchWKT(agstack_geoid)
#     lat, lon = extractLatLonFromWKT(wkt_polygon)

#     # Fetch ETo data from various sources
#     eto_val_aus = getEtoFromAUS(lat, lon, paths['aus'], start_date, end_date)
#     eto_val_aus_forecast = getEtoFromAUSForecast(lat, lon, paths['aus_forecast'], start_date, end_date)
#     eto_val_ncep = getEtoFromNCEP(lat, lon, paths['ncep'], start_date, end_date)
#     eto_val_nldas = getEtoFromNLDAS(lat, lon, paths['nldas'], start_date, end_date)
#     eto_val_noaa = getEtoFromNOAA(lat, lon, paths['noaa'], start_date, end_date)
#     eto_val_noaa_forecast = getEtoFromNOAAForecast(lat, lon, paths['noaa_forecasts'], start_date, end_date)
#     eto_val_global = getEtoFromGHCND(lat, lon, paths['global'], start_date, end_date)
    
#     # Helper function to extract the 'ETo_AVG_IN' value safely
#     def extract_eto_value(df):
#         return df['ETo_AVG_IN'].values[0] if not df.empty else None

#     # Extract the 'ETo_AVG_IN' values from each dataset
#     eto_values = [
#         extract_eto_value(eto_val_noaa),
#         extract_eto_value(eto_val_noaa_forecast),
#         extract_eto_value(eto_val_global),
#         extract_eto_value(eto_val_ncep),
#         extract_eto_value(eto_val_nldas),
#         extract_eto_value(eto_val_aus),
#         extract_eto_value(eto_val_aus_forecast)
#     ]
=======
        print(e)
        w_ret = empty_response()

    return w_ret

>>>>>>> Rajat

#     # Log the values (Optional)
#     for source, value in zip(["NOAA", "NOAA Forecast", "Global", "NCEP", "NLDAS", "AUS", "AUS Forecast"], eto_values):
#         print(f"{source} ETo: {value}")

<<<<<<< HEAD
#     # Filter out None, invalid data (extreme negative values), and non-numeric values
#     valid_eto_values = [val for val in eto_values if isinstance(val, (int, float)) and val > 0]
    
#     print(f"Valid ETo values: {valid_eto_values}")
    
#     # Calculate the average if there are valid values
#     average_eto = np.mean(valid_eto_values) if valid_eto_values else np.nan  # Using np.nan as a default value
    
    
    
#     return average_eto_df

def getEtoFromWeatherData(agstack_geoid, start_date, end_date):
    ncep_filePath = '/network/NCEP/' 
    nldas_filePath = '/network/NLDAS/'
    noaaPath = '/network/NOAA/PARQUET_S2/'
    noaa_forecasts_path = '/network/NOAA/FORECASTED/PARQUET_S2/'
    global_path = '/network/GHCND/PARQUET_S2/'
    aus_filePath = '/network/AUS/DAILY/PARQUET_S2/'
    aus_forecast_filePath = '/network/AUS/FORECASTED/PARQUET_S2/'
    cimis_filePath = '/network/CIMIS/HOURLY/PARQUET_S2/'

    wkt_polygon = fetchWKT(agstack_geoid)
    lat, lon = extractLatLonFromWKT(wkt_polygon)
    
    # Set date range
    date_range = [pd.to_datetime(start_date)] if end_date is None else pd.date_range(start=start_date, end=end_date, freq='D')

    # Initialize list to store daily ETo values
    daily_eto_values = []

    # Define data sources and their corresponding functions
    sources = [
        ('AUS', getEtoFromAUS, aus_filePath),
        ('AUS_FORECAST', getEtoFromAUSForecast, aus_forecast_filePath),
        ('NCEP', getEtoFromNCEP, ncep_filePath),
        ('NLDAS', getEtoFromNLDAS, nldas_filePath),
        ('NOAA', getEtoFromNOAA, noaaPath),
        ('NOAA_FORECAST', getEtoFromNOAAForecast, noaa_forecasts_path),
        ('GLOBAL', getEtoFromGHCND, global_path),
        ('CIMIS', getEToFromCIMIS, cimis_filePath),
    ]

    for date in date_range:
        date_str = date.strftime('%Y-%m-%d')

        eto_values = []
        for source_name, fetch_func, file_path in sources:
            # Fetch data and check if it's not empty
            eto_data = fetch_func(lat, lon, file_path, date_str, date_str)
            if not eto_data.empty:
                eto_value = eto_data['ETo_AVG_IN'].iloc[0]
                if 0 <= eto_value <= 1:  # Condition to filter valid ETo values
                    print(f"{source_name} gives valid ETo value {eto_value} for {date_str}")
                    eto_values.append(eto_value)
                else:
                    pass
                    # print(f"{source_name} gives invalid ETo value {eto_value} for {date_str}, skipping.")
            else:
                pass
                # print(f"{source_name} returned empty dataset for {date_str}")
        
        # Calculate the average ETo ignoring empty datasets
        valid_eto_values = [val for val in eto_values if val is not None]
        average_eto = round(np.mean(valid_eto_values), 2) if valid_eto_values else None
        
        # Store the result for the current day
        daily_eto_values.append({'date': date_str, 'average_eto': average_eto})

    if daily_eto_values:
        # Convert results to DataFrame
        daily_eto_df = pd.DataFrame(daily_eto_values)
        daily_eto_df['average_eto'] = daily_eto_df['average_eto'].fillna(0.0)
    else:
        daily_eto_df = pd.DataFrame()
        
    return daily_eto_df

def getWeatherFromSatelite(agstack_geoid, start_date, end_date=None):
    # Convert start and end dates using pd.to_datetime (much faster and cleaner)
    utc_start_dt = pd.to_datetime(start_date).tz_localize('UTC')
    
    # If end_date is None, set utc_end_dt to the same value as utc_start_dt
    if end_date is not None:    
        utc_end_dt = pd.to_datetime(end_date).tz_localize('UTC')
    else:
        utc_end_dt = utc_start_dt  # When end_date is None, use start_date as end_date

    w_ret = pd.DataFrame()

    s1_time_start = time.time()
    try:
        # Fetch WKT polygon and extract latitude and longitude
        wkt_polygon = fetchWKT(agstack_geoid)
        lat, lon = extractLatLonFromWKT(wkt_polygon)
    except Exception as e:
        print(e)
        return empty_response()
   
    try:
            # Fetch list of S2 paths and filter existing paths in one step
        s2_index_L8_list, _ = get_s2_cellids_and_token_list(8, [lat], [lon])
        s2_index_L10_list, _ = get_s2_cellids_and_token_list(10, [lat], [lon])

        list_of_L8_paths = [
            os.path.join(pa_filePath, f's2_index__L8={x}')
            for x in s2_index_L8_list if os.path.exists(os.path.join(pa_filePath, f's2_index__L8={x}'))
        ]
        list_of_L10_paths = [
            os.path.join(L8_path, f's2_index__L10={x}')
            for L8_path in list_of_L8_paths
            for x in s2_index_L10_list if os.path.exists(os.path.join(L8_path, f's2_index__L10={x}'))
        ]
        
        if not list_of_L10_paths:
            return empty_response()

        # Generate paths for L8 and L10
        weather_datasets = [ds.dataset(x, format="parquet", partitioning="hive") for x in list_of_L10_paths]
        
        # Generate date range for filtering
        date_range = pd.date_range(utc_start_dt, utc_end_dt, freq='D')
        
        # Create an empty list to accumulate DataFrames
        all_weather_data = []

        for weatherDataset in weather_datasets:
            for current_date in date_range:
                current_date_str = current_date.strftime('%Y-%m-%d')
                weather_df = weatherDataset.to_table(
                    filter=(
                        ds.field('YY_MM_DD') == current_date_str
                    )
                ).to_pandas()

                if not weather_df.empty:
                    all_weather_data.append(weather_df)

        if all_weather_data:
            w_all = pd.concat(all_weather_data, ignore_index=True)

            # Drop unnecessary columns and fill missing values
            w_all = w_all.drop(['s2_index__L20', 'latitude', 'longitude'], axis=1, errors='ignore')  # Avoid error if columns don't exist
            w_all = w_all.fillna(0)
            
            
            numeric_cols = w_all.select_dtypes(include=[np.number]).columns
            w_all[numeric_cols] = w_all[numeric_cols].apply(pd.to_numeric, errors='coerce')

            # Group by 'time' and calculate the mean for each group
            w_ret = w_all.groupby('UTC_DATETIME').mean(numeric_only=True).reset_index()
            w_ret['s2_index__L19'] = w_all['s2_index__L19']
    except Exception as e:
        print(e)
        w_ret = pd.DataFrame()
    print(w_ret)
    return w_ret


def getNDVIDatesForGeoid(agstack_geoid, start_date, end_date=None):
    
    filePath = "/network/SENTINEL/PARQUET_NDVI_L20/"
    tok = start_date.split('-')
    YYYY_str = tok[0]
    # print(YYYY_str)
    MM_str = tok[1]
    DD_str = tok[2]
    
    # Create a datetime object and localize it to UTC
    local_dt = datetime(int(YYYY_str), int(MM_str), int(DD_str))
    utc_start_dt = pytz.utc.localize(local_dt)  # Convert to UTC
    
    if end_date:
        # Parse end date and localize it to UTC if provided
        end_tok = end_date.split('-')
        end_YYYY_str = int(end_tok[0])
        end_MM_str = int(end_tok[1])
        end_DD_str = int(end_tok[2])
        end_local_dt = datetime(end_YYYY_str, end_MM_str, end_DD_str)
        utc_end_dt = pytz.utc.localize(end_local_dt)

    # Fetch WKT polygon and extract latitude and longitude
=======
# AUS
def getWeatherFromAUS(agstack_geoid, dtStr,end_date=None):
    # filePath = '/mnt/md1/NLDAS/PARQUETE_S2/'
    # filePath = '/home/rajat/Downloads/Rnaura_Work/mnt/md1/AUS/DAILY/PROCESSED/PARQUET_S2_TEST/'
    filePath = "/network/AUS/DAILY/PARQUET_S2/"

    tok = dtStr.split('-')
    YYYY_str = tok[0]
    MM_str = tok[1]
    DD_str = tok[2]

    # Create a datetime object and localize it to UTC
    local_dt = datetime(int(YYYY_str), int(MM_str), int(DD_str))
    utc_dt = pytz.utc.localize(local_dt)  # Convert to UTC

    if end_date:
        # Parse end date and localize it to UTC if provided
        end_tok = end_date.split('-')
        end_YYYY_str = int(end_tok[0])
        end_MM_str = int(end_tok[1])
        end_DD_str = int(end_tok[2])
        end_local_dt = datetime(end_YYYY_str, end_MM_str, end_DD_str)
        utc_end_dt = pytz.utc.localize(end_local_dt)

    try:
        

        s1_time_start = time.time()
        # Fetch WKT polygon and extract latitude and longitude
        wkt_polygon = fetchWKT(agstack_geoid)
        lat, lon = extractLatLonFromWKT(wkt_polygon)
        start_time = time.time()
        # Get the list of S2 indices and CIDs for the data point
        s2_index__L5_list, L5_cids = get_s2_cellids_and_token_list(5, [lat], [lon])
        print(f'token---{s2_index__L5_list}')
        
        list_of_5_paths = [filePath + 's2_token_L5=' + x for x in s2_index__L5_list
                        if os.path.exists(filePath + 's2_token_L5=' + x)]
        if not list_of_5_paths:
            return empty_response()
        
        weather_datasets = []
        for x in list_of_5_paths:
            weather_datasets.append(ds.dataset(x, format="parquet", partitioning="hive"))

        w_all = pd.DataFrame()
        w_ret = pd.DataFrame()
        for weatherDataset in weather_datasets:
            if end_date is not None:
                # Filter between start and end date
                weather_df = weatherDataset.to_table(
                    filter=(
                        (ds.field('Year') >= int(YYYY_str)) & (ds.field('Year') <= end_YYYY_str) &
                        (ds.field('Month') >= int(MM_str)) & (ds.field('Month') <= end_MM_str) &
                        (ds.field('Day') >= int(DD_str)) & (ds.field('Day') <= end_DD_str)
                    )
                ).to_pandas()
            elif end_date is None:
                # Filter for start date only
                weather_df = weatherDataset.to_table(
                    filter=(
                        (ds.field('Year') == int(YYYY_str)) &
                        (ds.field('Month') == int(MM_str)) &
                        (ds.field('Day') == int(DD_str))
                    )
                ).to_pandas()

            # if not weather_df.empty:
            #     weather_df.dropna(inplace=True)
            print(weather_df.columns)
            w_all = pd.concat([w_all, weather_df], ignore_index=True)
        if len(w_all) > 0:
            w_all = w_all.fillna(0)
            # Convert 'time' column to datetime format if it isn't already
            # After constructing w_all DataFrame and filling Year, Month, and Day
            if 'Year' in w_all.columns and 'Month' in w_all.columns and 'Day' in w_all.columns:
                # Create 'aifstime_utc' based on Year, Month, and Day
                w_all['aifstime_utc'] = pd.to_datetime(
                    w_all[['Year', 'Month', 'Day']].astype(str).agg('-'.join, axis=1)
                ).dt.tz_localize('UTC')

            # Now you can convert the 'aifstime_utc' column
            w_all['aifstime_utc'] = pd.to_datetime(w_all['aifstime_utc'])
            # print("======w_all=====",w_all)
            # Ensure all numeric columns are converted to numeric type
            numeric_cols = w_all.select_dtypes(include=[np.number]).columns
            w_all[numeric_cols] = w_all[numeric_cols].apply(pd.to_numeric, errors='coerce')
            print("Columns in w_all before aggregation:", w_all.columns)

            # Group by 'time' and calculate the mean for each group
            w_ret = w_all.groupby('aifstime_utc').mean(numeric_only=True).reset_index()
            
            print("Non-numeric columns:", w_all.select_dtypes(exclude=[np.number]).columns)

            # print("=======w_ret========",w_ret)
            # Ensure 'Year', 'Month', and 'Day' are included in the response
            w_ret['Year'] = w_ret['Year'].astype(int)
            w_ret['Month'] = w_ret['Month'].astype(int)
            w_ret['Day'] = w_ret['Day'].astype(int)
            w_ret['cloud'] = w_all['cloud']
            w_ret['cloud_type'] = w_all['cloud_type']
            w_ret['name'] = w_all['name']
            w_ret['press_tend'] = w_all['press_tend']
            # w_ret['rain_trace'] = w_all['rain_trace']
            w_ret['local_date_time'] = w_all['local_date_time']
            w_ret['local_date_time_full'] = w_all['local_date_time_full']
            w_ret['history_product'] = w_all['history_product']
            w_ret['local_date_time_full'] = w_all['local_date_time_full']
            w_ret['vis_km'] = w_all['vis_km']
            w_ret['weather'] = w_all['weather']
            w_ret['wind_dir'] = w_all['wind_dir']
            
            
            
            print(f'w_ret--{w_ret.columns}')
            
            # Ensure TS column is created
            # if not weather_df.empty:
            #     weather_df['aifstime_utc'] = pd.to_datetime(weather_df[['Year', 'Month', 'Day']].astype(str).agg('-'.join, axis=1)).dt.tz_localize('UTC')

            #     numeric_cols = weather_df.select_dtypes(include=[np.number]).columns
            #     w_df = pd.DataFrame(weather_df[numeric_cols].mean(axis=0)).T

            #     # Retain the first TS value or handle TS appropriately
            #     if 'aifstime_utc' in weather_df.columns and not weather_df['aifstime_utc'].empty:
            #         w_df['aifstime_utc'] = weather_df['aifstime_utc'].iloc[0]
            #     else:
            #         w_df['aifstime_utc'] = pd.NaT  # or some other default value

            #     w_all = pd.concat([w_all, w_df], ignore_index=True)

        # if not w_all.empty:
        #     w_all = w_all.fillna(0)
        #     # Take the average of the columns
        #     w_ret = pd.DataFrame(w_all.groupby(['aifstime_utc']).mean())
        #     w_ret.reset_index(inplace=True)
        # else:
        #     w_ret = pd.DataFrame()
    except Exception as e:
        print(e)
        
        
    end_time = time.time()
    time_elapsed = (end_time - start_time)

    return w_ret


# def getEtoFromWeatherData(agstack_geoid, start_date, end_date):
#     # Define file paths
#     paths = {
#         "ncep": '/mnt/md1/NCEP/PARQUET_S2/', 
#         "nldas": '/mnt/md1/NLDAS/PARQUET_S2/',
#         "noaa": '/mnt/md1/NOAA/DAILY/PROCESSED/PARQUET_S2/',
#         "noaa_forecasts": '/mnt/md1/NOAA/FORECASTED/PARQUET_S2/',
#         "global": '/home/rnaura/mnt/md1/GHCND/DAILY/PROCESSED/PARQUET_S2/',
#         "aus": '/mnt/md1/AUS/DAILY/PROCESSED/PARQUET_S2/',
#         "aus_forecast": '/mnt/md1/AUS/DAILY/FORECASTED/PARQUET_S2/'
#     }
    
#     # Fetch WKT and extract latitude and longitude
#     wkt_polygon = fetchWKT(agstack_geoid)
#     lat, lon = extractLatLonFromWKT(wkt_polygon)

#     # Fetch ETo data from various sources
#     eto_val_aus = getEtoFromAUS(lat, lon, paths['aus'], start_date, end_date)
#     eto_val_aus_forecast = getEtoFromAUSForecast(lat, lon, paths['aus_forecast'], start_date, end_date)
#     eto_val_ncep = getEtoFromNCEP(lat, lon, paths['ncep'], start_date, end_date)
#     eto_val_nldas = getEtoFromNLDAS(lat, lon, paths['nldas'], start_date, end_date)
#     eto_val_noaa = getEtoFromNOAA(lat, lon, paths['noaa'], start_date, end_date)
#     eto_val_noaa_forecast = getEtoFromNOAAForecast(lat, lon, paths['noaa_forecasts'], start_date, end_date)
#     eto_val_global = getEtoFromGHCND(lat, lon, paths['global'], start_date, end_date)
    
#     # Helper function to extract the 'ETo_AVG_IN' value safely
#     def extract_eto_value(df):
#         return df['ETo_AVG_IN'].values[0] if not df.empty else None

#     # Extract the 'ETo_AVG_IN' values from each dataset
#     eto_values = [
#         extract_eto_value(eto_val_noaa),
#         extract_eto_value(eto_val_noaa_forecast),
#         extract_eto_value(eto_val_global),
#         extract_eto_value(eto_val_ncep),
#         extract_eto_value(eto_val_nldas),
#         extract_eto_value(eto_val_aus),
#         extract_eto_value(eto_val_aus_forecast)
#     ]

#     # Log the values (Optional)
#     for source, value in zip(["NOAA", "NOAA Forecast", "Global", "NCEP", "NLDAS", "AUS", "AUS Forecast"], eto_values):
#         print(f"{source} ETo: {value}")

#     # Filter out None, invalid data (extreme negative values), and non-numeric values
#     valid_eto_values = [val for val in eto_values if isinstance(val, (int, float)) and val > 0]
    
#     print(f"Valid ETo values: {valid_eto_values}")
    
#     # Calculate the average if there are valid values
#     average_eto = np.mean(valid_eto_values) if valid_eto_values else np.nan  # Using np.nan as a default value
    
    
    
#     return average_eto_df

def getEtoFromWeatherData(agstack_geoid, start_date, end_date):
    ncep_filePath = '/network/NCEP/' 
    nldas_filePath = '/network/NLDAS/'
    noaaPath = '/network/NOAA/PARQUET_S2/'
    noaa_forecasts_path = '/network/NOAA/FORECASTED/PARQUET_S2/'
    global_path = '/network/GHCND/PARQUET_S2/'
    aus_filePath = '/network/AUS/DAILY/PARQUET_S2/'
    aus_forecast_filePath = '/network/AUS/FORECASTED/PARQUET_S2/'
    cimis_filePath = '/network/CIMIS/HOURLY/PARQUET_S2/'

>>>>>>> Rajat
    wkt_polygon = fetchWKT(agstack_geoid)
    lat, lon = extractLatLonFromWKT(wkt_polygon)
    
    # Set date range
    date_range = [pd.to_datetime(start_date)] if end_date is None else pd.date_range(start=start_date, end=end_date, freq='D')

    # Initialize list to store daily ETo values
    daily_eto_values = []

    # Define data sources and their corresponding functions
    sources = [
        ('AUS', getEtoFromAUS, aus_filePath),
        ('AUS_FORECAST', getEtoFromAUSForecast, aus_forecast_filePath),
        ('NCEP', getEtoFromNCEP, ncep_filePath),
        ('NLDAS', getEtoFromNLDAS, nldas_filePath),
        ('NOAA', getEtoFromNOAA, noaaPath),
        ('NOAA_FORECAST', getEtoFromNOAAForecast, noaa_forecasts_path),
        ('GLOBAL', getEtoFromGHCND, global_path),
        ('CIMIS', getEToFromCIMIS, cimis_filePath),
    ]

<<<<<<< HEAD
    # Get the list of S2 indices and CIDs for the data point
        
    s2_index__L8_list, L8_cids = get_s2_cellids_and_token_list(8, [lat], [lon])
    print(s2_index__L8_list)

    list_of_8_paths = [filePath + 's2_index__L8=' + x for x in s2_index__L8_list
                    if os.path.exists(filePath + 's2_index__L8=' + x)]
    if not list_of_8_paths:
        return empty_response()
    
    start_date_str = utc_start_dt.strftime('%Y-%m-%d')

    satellite_datasets = []
    for x in list_of_8_paths:
        satellite_datasets.append(ds.dataset(x, format="parquet", partitioning="hive"))

    s_all = pd.DataFrame()
    for SatelliteDataset in satellite_datasets:
        if end_date is not None:
            end_date_str = utc_end_dt.strftime('%Y-%m-%d')
            # Filter between start and end date
            s_all = SatelliteDataset.to_table(
                columns=['YY_MM_DD'],
                filter=(
                    (ds.field('YY_MM_DD') >= start_date_str) & (ds.field('YY_MM_DD') <= end_date_str)
                )
            ).to_pandas()
            
            print(s_all)

        elif end_date is None:
            # Filter for start date only
            
            s_all = SatelliteDataset.to_table(
                columns=['YY_MM_DD'],
                filter=(
                    (ds.field('YY_MM_DD') == start_date_str)
                )
            ).to_pandas()
        
    return s_all


def get_all_available_dates(filePath):
    """Extract all available dates from the directory structure."""
    available_dates = []
    
    # Traverse through the directory and collect dates from the YY_MM_DD folder names
    for root, dirs, files in os.walk(filePath):
        for dir_name in dirs:
            if dir_name.startswith("YY_MM_DD="):
                # Extract the date from the directory name
                date_str = dir_name.split("YY_MM_DD=")[-1]
                available_dates.append(date_str)
    print(f'available_dates--{available_dates}')
    
    return available_dates



def get_closest_historical_date(available_dates, target_date):
    """Find the closest historical date prior to the target date."""
    # Convert string dates to datetime objects for comparison
    available_dates = [datetime.strptime(date, '%Y-%m-%d').replace(tzinfo=timezone.utc) for date in available_dates]
    
    # Make sure the target_date is also aware
    if target_date.tzinfo is None:
        target_date = target_date.replace(tzinfo=timezone.utc)

    # Sort available dates
    available_dates.sort()

    # Find the most recent date that is less than or equal to the query date
    recent_date = None
    for date in available_dates:
        if date <= target_date:
            recent_date = date
        else:
            break

    return recent_date.strftime('%Y-%m-%d') if recent_date else None

def calculate_stats(weather_df):
    """Calculate required statistics for a given weather_df."""
    if 'NDVI' not in weather_df.columns:
        print("NDVI column not found in the DataFrame.")
        return None

    ndvi_data = weather_df['NDVI']
    
    
    stats = {
        "min": ndvi_data.min(),
        "max": ndvi_data.max(),
        "avg": ndvi_data.mean(),
        "median": ndvi_data.median(),
        "stddev": ndvi_data.std(),
        "5th_prctile": ndvi_data.quantile(0.05),
        "10th_prctile": ndvi_data.quantile(0.10),
        "25th_prctile": ndvi_data.quantile(0.25),
        "75th_prctile": ndvi_data.quantile(0.75),
        "90th_prctile": ndvi_data.quantile(0.90),
        "95th_prctile": ndvi_data.quantile(0.95),
    }
    
    return pd.Series(stats)


# Define an empty result structure
def empty_result():
    columns = ["Date", "NDVI", "min", "max", "median","stddev","5th_prctile","10th_prctile","25th_prctile","75th_prctile","90th_prctile","95th_prctile"]  # Adjust with the actual column names
    return pd.DataFrame(columns=columns).fillna("")

def getSateliteStatsFn(agstack_geoid, start_date, end_date=None):
    filePath = "/network/SENTINEL/PARQUET_NDVI_L20/"
    
    # Parse start date
    try:
        utc_start_dt = pytz.utc.localize(datetime.strptime(start_date, '%Y-%m-%d'))
    except Exception as e:
        print(e)
        return empty_result()
        
    # Parse end date
    if end_date:
        try:
            utc_end_dt = pytz.utc.localize(datetime.strptime(end_date, '%Y-%m-%d'))
        except (ValueError, IndexError) as e:
            print(f"Error parsing end_date: {e}")
            return empty_result()
    else:
        utc_end_dt = utc_start_dt

    # Fetch WKT and extract latitude and longitude
    try:
        wkt_polygon = fetchWKT(agstack_geoid)
        lat, lon = extractLatLonFromWKT(wkt_polygon)
    except Exception as e:
        print(f"Error fetching WKT or extracting lat/lon: {e}")
        return empty_result()
    
    list_of_L8_paths = []
    list_of_L10_paths = []
    # Get the list of S2 indices
    try:
        s2_index_L8_list, _ = get_s2_cellids_and_token_list(8, [lat], [lon])
        s2_index_L10_list, _ = get_s2_cellids_and_token_list(10, [lat], [lon])

        list_of_L8_paths = [
            os.path.join(pa_filePath, f's2_index__L8={x}')
            for x in s2_index_L8_list if os.path.exists(os.path.join(pa_filePath, f's2_index__L8={x}'))
        ]
        list_of_L10_paths = [
            os.path.join(L8_path, f's2_index__L10={x}')
            for L8_path in list_of_L8_paths
            for x in s2_index_L10_list if os.path.exists(os.path.join(L8_path, f's2_index__L10={x}'))
        ]
        
        if not list_of_L10_paths:
            print("No paths found.")
            return empty_result()
    except Exception as e:
        print(f"Error getting S2 indices or building paths: {e}")
        return empty_result()

    start_date_str = utc_start_dt.strftime('%Y-%m-%d')

    # Check date availability and get closest date
    try:
        available_dates = get_all_available_dates(list_of_L10_paths[0])
        date_counts = Counter(available_dates)
        print(f'Date counts: {date_counts}')
        
        if start_date_str not in available_dates:
            print(f"No data available for {start_date_str}. Checking for the closest historical date...")
            closest_date_str = get_closest_historical_date(available_dates, utc_start_dt)
            print(f'closest_date_str--{closest_date_str}')
            if not closest_date_str:
                print("No historical data available.")
                return empty_result()
        else:
            closest_date_str = start_date_str
    except Exception as e:
        print(f"Error checking date availability or finding closest date: {e}")
        return empty_result()

    # Process data for each dataset and date range
    stats_list = []
    try:
        date_range = pd.date_range(start=closest_date_str, end=utc_end_dt.date()).strftime('%Y-%m-%d').tolist() if end_date else [closest_date_str]
    except Exception as e:
        date_range = pd.date_range(start=start_date_str, end=utc_end_dt.date()).strftime('%Y-%m-%d').tolist() if end_date else [closest_date_str]
        print(f'An error occurred - {e}')

    try:
        weather_datasets = [ds.dataset(path, format="parquet", partitioning="hive") for path in list_of_L10_paths]
        
        for weatherDataset in weather_datasets:
            for date in date_range:
                try:
                    date_timestamp = pyarrow.scalar(datetime.strptime(date, '%Y-%m-%d'), type=pyarrow.timestamp('us'))
                    satelite_df = weatherDataset.to_table(
                        filter=(ds.field('UTC_DATETIME') == date_timestamp)
                    ).to_pandas()
                    
                    # Check if the DataFrame is empty
                    if satelite_df.empty:
                        print(f"No data available for date {date}.")
                        continue  # Skip to the next date

                    # Safely access UTC_DATETIME
                    utc_dt = satelite_df['UTC_DATETIME'].iloc[0]

                    # Calculate statistics
                    stats = calculate_stats(satelite_df)
                    if stats is not None:
                        stats['date'] = utc_dt
                        stats_list.append(stats)

                except Exception as e:
                    print(f"Error processing data for date {date}: {e}")
                    return empty_result()

    except Exception as e:
        print(f"Error setting up weather datasets: {e}")
        return empty_result()

    # Combine all stats into a single DataFrame
    if stats_list:
        try:
            s_all = pd.DataFrame(stats_list)
            s_all = s_all.dropna()
        except Exception as e:
            print(f"Error combining statistics: {e}")
            return empty_result()
    else:
        print("No statistics calculated.")
        return empty_result()

    # Process based on the number of unique dates
    if len(date_range) > 1:
        # More than 1 date, index by Date and resample
        s_all.set_index('date', inplace=True)

        # Ensure index is timezone-aware (convert if necessary)
        if s_all.index.tz is None:
            s_all.index = s_all.index.tz_localize('UTC')

        try:
            duplicates = s_all.index[s_all.index.duplicated()]
            if not duplicates.empty:
                print("Duplicate index labels found:", duplicates)
            s_all = s_all[~s_all.index.duplicated(keep='first')]  # Keep the first occurrence
            s_all.index = pd.to_datetime(s_all.index)

        except Exception as e:
            print(f'error -{e}')

=======
    for date in date_range:
        date_str = date.strftime('%Y-%m-%d')

        eto_values = []
        for source_name, fetch_func, file_path in sources:
            # Fetch data and check if it's not empty
            eto_data = fetch_func(lat, lon, file_path, date_str, date_str)
            if not eto_data.empty:
                eto_value = eto_data['ETo_AVG_IN'].iloc[0]
                if 0 <= eto_value <= 1:  # Condition to filter valid ETo values
                    print(f"{source_name} gives valid ETo value {eto_value} for {date_str}")
                    eto_values.append(eto_value)
                else:
                    pass
                    # print(f"{source_name} gives invalid ETo value {eto_value} for {date_str}, skipping.")
            else:
                pass
                # print(f"{source_name} returned empty dataset for {date_str}")
        
        # Calculate the average ETo ignoring empty datasets
        valid_eto_values = [val for val in eto_values if val is not None]
        average_eto = round(np.mean(valid_eto_values), 2) if valid_eto_values else None
        
        # Store the result for the current day
        daily_eto_values.append({'date': date_str, 'average_eto': average_eto})

    if daily_eto_values:
        # Convert results to DataFrame
        daily_eto_df = pd.DataFrame(daily_eto_values)
        daily_eto_df['average_eto'] = daily_eto_df['average_eto'].fillna(0.0)
    else:
        daily_eto_df = pd.DataFrame()
        
    return daily_eto_df

def getWeatherFromSatelite(agstack_geoid, start_date, end_date=None):
    # Convert start and end dates using pd.to_datetime (much faster and cleaner)
    utc_start_dt = pd.to_datetime(start_date).tz_localize('UTC')
    
    # If end_date is None, set utc_end_dt to the same value as utc_start_dt
    if end_date is not None:    
        utc_end_dt = pd.to_datetime(end_date).tz_localize('UTC')
    else:
        utc_end_dt = utc_start_dt  # When end_date is None, use start_date as end_date

    w_ret = pd.DataFrame()

    s1_time_start = time.time()
    try:
        # Fetch WKT polygon and extract latitude and longitude
        wkt_polygon = fetchWKT(agstack_geoid)
        lat, lon = extractLatLonFromWKT(wkt_polygon)
    except Exception as e:
        print(e)
        return empty_response()
   
    try:
            # Fetch list of S2 paths and filter existing paths in one step
        s2_index_L8_list, _ = get_s2_cellids_and_token_list(8, [lat], [lon])
        s2_index_L10_list, _ = get_s2_cellids_and_token_list(10, [lat], [lon])

        list_of_L8_paths = [
            os.path.join(pa_filePath, f's2_index__L8={x}')
            for x in s2_index_L8_list if os.path.exists(os.path.join(pa_filePath, f's2_index__L8={x}'))
        ]
        list_of_L10_paths = [
            os.path.join(L8_path, f's2_index__L10={x}')
            for L8_path in list_of_L8_paths
            for x in s2_index_L10_list if os.path.exists(os.path.join(L8_path, f's2_index__L10={x}'))
        ]
        
        if not list_of_L10_paths:
            return empty_response()

        # Generate paths for L8 and L10
        weather_datasets = [ds.dataset(x, format="parquet", partitioning="hive") for x in list_of_L10_paths]
        
        # Generate date range for filtering
        date_range = pd.date_range(utc_start_dt, utc_end_dt, freq='D')
        
        # Create an empty list to accumulate DataFrames
        all_weather_data = []

        for weatherDataset in weather_datasets:
            for current_date in date_range:
                current_date_str = current_date.strftime('%Y-%m-%d')
                weather_df = weatherDataset.to_table(
                    filter=(
                        ds.field('YY_MM_DD') == current_date_str
                    )
                ).to_pandas()

                if not weather_df.empty:
                    all_weather_data.append(weather_df)

        if all_weather_data:
            w_all = pd.concat(all_weather_data, ignore_index=True)

            # Drop unnecessary columns and fill missing values
            w_all = w_all.drop(['s2_index__L20', 'latitude', 'longitude'], axis=1, errors='ignore')  # Avoid error if columns don't exist
            w_all = w_all.fillna(0)
            
            
            numeric_cols = w_all.select_dtypes(include=[np.number]).columns
            w_all[numeric_cols] = w_all[numeric_cols].apply(pd.to_numeric, errors='coerce')

            # Group by 'time' and calculate the mean for each group
            w_ret = w_all.groupby('UTC_DATETIME').mean(numeric_only=True).reset_index()
            w_ret['s2_index__L19'] = w_all['s2_index__L19']
    except Exception as e:
        print(e)
        w_ret = pd.DataFrame()
    print(w_ret)
    return w_ret


def getNDVIDatesForGeoid(agstack_geoid, start_date, end_date=None):
    
    filePath = "/network/SENTINEL/PARQUET_NDVI_L20/"
    tok = start_date.split('-')
    YYYY_str = tok[0]
    # print(YYYY_str)
    MM_str = tok[1]
    DD_str = tok[2]
    
    # Create a datetime object and localize it to UTC
    local_dt = datetime(int(YYYY_str), int(MM_str), int(DD_str))
    utc_start_dt = pytz.utc.localize(local_dt)  # Convert to UTC
    
    if end_date:
        # Parse end date and localize it to UTC if provided
        end_tok = end_date.split('-')
        end_YYYY_str = int(end_tok[0])
        end_MM_str = int(end_tok[1])
        end_DD_str = int(end_tok[2])
        end_local_dt = datetime(end_YYYY_str, end_MM_str, end_DD_str)
        utc_end_dt = pytz.utc.localize(end_local_dt)

    # Fetch WKT polygon and extract latitude and longitude
    wkt_polygon = fetchWKT(agstack_geoid)
    lat, lon = extractLatLonFromWKT(wkt_polygon)

    # Get the list of S2 indices and CIDs for the data point
        
    s2_index__L8_list, L8_cids = get_s2_cellids_and_token_list(8, [lat], [lon])
    print(s2_index__L8_list)

    list_of_8_paths = [filePath + 's2_index__L8=' + x for x in s2_index__L8_list
                    if os.path.exists(filePath + 's2_index__L8=' + x)]
    if not list_of_8_paths:
        return empty_response()
    
    start_date_str = utc_start_dt.strftime('%Y-%m-%d')

    satellite_datasets = []
    for x in list_of_8_paths:
        satellite_datasets.append(ds.dataset(x, format="parquet", partitioning="hive"))

    s_all = pd.DataFrame()
    for SatelliteDataset in satellite_datasets:
        if end_date is not None:
            end_date_str = utc_end_dt.strftime('%Y-%m-%d')
            # Filter between start and end date
            s_all = SatelliteDataset.to_table(
                columns=['YY_MM_DD'],
                filter=(
                    (ds.field('YY_MM_DD') >= start_date_str) & (ds.field('YY_MM_DD') <= end_date_str)
                )
            ).to_pandas()
            
            print(s_all)

        elif end_date is None:
            # Filter for start date only
            
            s_all = SatelliteDataset.to_table(
                columns=['YY_MM_DD'],
                filter=(
                    (ds.field('YY_MM_DD') == start_date_str)
                )
            ).to_pandas()
        
    return s_all


def get_all_available_dates(filePath):
    """Extract all available dates from the directory structure."""
    available_dates = []
    
    # Traverse through the directory and collect dates from the YY_MM_DD folder names
    for root, dirs, files in os.walk(filePath):
        for dir_name in dirs:
            if dir_name.startswith("YY_MM_DD="):
                # Extract the date from the directory name
                date_str = dir_name.split("YY_MM_DD=")[-1]
                available_dates.append(date_str)
    print(f'available_dates--{available_dates}')
    
    return available_dates



def get_closest_historical_date(available_dates, target_date):
    """Find the closest historical date prior to the target date."""
    # Convert string dates to datetime objects for comparison
    available_dates = [datetime.strptime(date, '%Y-%m-%d').replace(tzinfo=timezone.utc) for date in available_dates]
    
    # Make sure the target_date is also aware
    if target_date.tzinfo is None:
        target_date = target_date.replace(tzinfo=timezone.utc)

    # Sort available dates
    available_dates.sort()

    # Find the most recent date that is less than or equal to the query date
    recent_date = None
    for date in available_dates:
        if date <= target_date:
            recent_date = date
        else:
            break

    return recent_date.strftime('%Y-%m-%d') if recent_date else None

def calculate_stats(weather_df):
    """Calculate required statistics for a given weather_df."""
    if 'NDVI' not in weather_df.columns:
        print("NDVI column not found in the DataFrame.")
        return None

    ndvi_data = weather_df['NDVI']
    
    
    stats = {
        "min": ndvi_data.min(),
        "max": ndvi_data.max(),
        "avg": ndvi_data.mean(),
        "median": ndvi_data.median(),
        "stddev": ndvi_data.std(),
        "5th_prctile": ndvi_data.quantile(0.05),
        "10th_prctile": ndvi_data.quantile(0.10),
        "25th_prctile": ndvi_data.quantile(0.25),
        "75th_prctile": ndvi_data.quantile(0.75),
        "90th_prctile": ndvi_data.quantile(0.90),
        "95th_prctile": ndvi_data.quantile(0.95),
    }
    
    return pd.Series(stats)


# Define an empty result structure
def empty_result():
    columns = ["Date", "NDVI", "min", "max", "median","stddev","5th_prctile","10th_prctile","25th_prctile","75th_prctile","90th_prctile","95th_prctile"]  # Adjust with the actual column names
    return pd.DataFrame(columns=columns).fillna("")

def getSateliteStatsFn(agstack_geoid, start_date, end_date=None):
    filePath = "/network/SENTINEL/PARQUET_NDVI_L20/"
    
    # Parse start date
    try:
        utc_start_dt = pytz.utc.localize(datetime.strptime(start_date, '%Y-%m-%d'))
    except Exception as e:
        print(e)
        return empty_result()
        
    # Parse end date
    if end_date:
        try:
            utc_end_dt = pytz.utc.localize(datetime.strptime(end_date, '%Y-%m-%d'))
        except (ValueError, IndexError) as e:
            print(f"Error parsing end_date: {e}")
            return empty_result()
    else:
        utc_end_dt = utc_start_dt

    # Fetch WKT and extract latitude and longitude
    try:
        wkt_polygon = fetchWKT(agstack_geoid)
        lat, lon = extractLatLonFromWKT(wkt_polygon)
    except Exception as e:
        print(f"Error fetching WKT or extracting lat/lon: {e}")
        return empty_result()
    
    list_of_L8_paths = []
    list_of_L10_paths = []
    # Get the list of S2 indices
    try:
        s2_index_L8_list, _ = get_s2_cellids_and_token_list(8, [lat], [lon])
        s2_index_L10_list, _ = get_s2_cellids_and_token_list(10, [lat], [lon])

        list_of_L8_paths = [
            os.path.join(pa_filePath, f's2_index__L8={x}')
            for x in s2_index_L8_list if os.path.exists(os.path.join(pa_filePath, f's2_index__L8={x}'))
        ]
        list_of_L10_paths = [
            os.path.join(L8_path, f's2_index__L10={x}')
            for L8_path in list_of_L8_paths
            for x in s2_index_L10_list if os.path.exists(os.path.join(L8_path, f's2_index__L10={x}'))
        ]
        
        if not list_of_L10_paths:
            print("No paths found.")
            return empty_result()
    except Exception as e:
        print(f"Error getting S2 indices or building paths: {e}")
        return empty_result()

    start_date_str = utc_start_dt.strftime('%Y-%m-%d')

    # Check date availability and get closest date
    try:
        available_dates = get_all_available_dates(list_of_L10_paths[0])
        date_counts = Counter(available_dates)
        print(f'Date counts: {date_counts}')
        
        if start_date_str not in available_dates:
            print(f"No data available for {start_date_str}. Checking for the closest historical date...")
            closest_date_str = get_closest_historical_date(available_dates, utc_start_dt)
            print(f'closest_date_str--{closest_date_str}')
            if not closest_date_str:
                print("No historical data available.")
                return empty_result()
        else:
            closest_date_str = start_date_str
    except Exception as e:
        print(f"Error checking date availability or finding closest date: {e}")
        return empty_result()

    # Process data for each dataset and date range
    stats_list = []
    try:
        date_range = pd.date_range(start=closest_date_str, end=utc_end_dt.date()).strftime('%Y-%m-%d').tolist() if end_date else [closest_date_str]
    except Exception as e:
        date_range = pd.date_range(start=start_date_str, end=utc_end_dt.date()).strftime('%Y-%m-%d').tolist() if end_date else [closest_date_str]
        print(f'An error occurred - {e}')

    try:
        weather_datasets = [ds.dataset(path, format="parquet", partitioning="hive") for path in list_of_L10_paths]
        
        for weatherDataset in weather_datasets:
            for date in date_range:
                try:
                    date_timestamp = pyarrow.scalar(datetime.strptime(date, '%Y-%m-%d'), type=pyarrow.timestamp('us'))
                    satelite_df = weatherDataset.to_table(
                        filter=(ds.field('UTC_DATETIME') == date_timestamp)
                    ).to_pandas()
                    
                    # Check if the DataFrame is empty
                    if satelite_df.empty:
                        print(f"No data available for date {date}.")
                        continue  # Skip to the next date

                    # Safely access UTC_DATETIME
                    utc_dt = satelite_df['UTC_DATETIME'].iloc[0]

                    # Calculate statistics
                    stats = calculate_stats(satelite_df)
                    if stats is not None:
                        stats['date'] = utc_dt
                        stats_list.append(stats)

                except Exception as e:
                    print(f"Error processing data for date {date}: {e}")
                    return empty_result()

    except Exception as e:
        print(f"Error setting up weather datasets: {e}")
        return empty_result()

    # Combine all stats into a single DataFrame
    if stats_list:
        try:
            s_all = pd.DataFrame(stats_list)
            s_all = s_all.dropna()
        except Exception as e:
            print(f"Error combining statistics: {e}")
            return empty_result()
    else:
        print("No statistics calculated.")
        return empty_result()

    # Process based on the number of unique dates
    if len(date_range) > 1:
        # More than 1 date, index by Date and resample
        s_all.set_index('date', inplace=True)

        # Ensure index is timezone-aware (convert if necessary)
        if s_all.index.tz is None:
            s_all.index = s_all.index.tz_localize('UTC')

        try:
            duplicates = s_all.index[s_all.index.duplicated()]
            if not duplicates.empty:
                print("Duplicate index labels found:", duplicates)
            s_all = s_all[~s_all.index.duplicated(keep='first')]  # Keep the first occurrence
            s_all.index = pd.to_datetime(s_all.index)

        except Exception as e:
            print(f'error -{e}')

>>>>>>> Rajat

        # Interpolate missing data
        s_resampled = s_all.resample('D').interpolate(method='linear')
        s_resampled.reset_index(inplace=True)

        # Extrapolation if the last date is less than the requested end_date
        last_date = s_resampled['date'].iloc[-1]
        if last_date < utc_end_dt:
            
            new_dates = pd.date_range(start=last_date + pd.Timedelta(days=1), end=utc_end_dt, freq='D')
            new_rows = pd.DataFrame({'date': new_dates})
            # new_rows.set_index('date', inplace=True)

            # Concatenate new rows with the existing DataFrame
            s_resampled = pd.concat([s_resampled, new_rows], axis=0)
            s_resampled.set_index('date', inplace=True)
            s_resampled = s_resampled.interpolate(method='linear')  # Extrapolate linearly
            s_resampled.reset_index(inplace=True)

            
        # Filter the DataFrame for the date range
        s_resampled = s_resampled[(s_resampled['date'].dt.date >= utc_start_dt.date()) & (s_resampled['date'].dt.date <= utc_end_dt.date())]

        return s_resampled  # Return the DataFrame directly
    else:
        # Single date, return the DataFrame
        return s_all




# AUS_FORECASTED
def getWeatherFromAUSFORECASTED(agstack_geoid, start_date,end_date=None):
    filePath = "/network/AUS/FORECASTED/PARQUET_S2/"
    tok = start_date.split('-')
    print(tok)
    YYYY_str = tok[0]
    print(YYYY_str)
    MM_str = tok[1]
    DD_str = tok[2]

    # Create a datetime object and localize it to UTC
    local_dt = datetime(int(YYYY_str), int(MM_str), int(DD_str))
    utc_dt = pytz.utc.localize(local_dt)  # Convert to UTC

    if end_date:
        # Parse end date and localize it to UTC if provided
        end_tok = end_date.split('-')
        end_YYYY_str = int(end_tok[0])
        end_MM_str = int(end_tok[1])
        end_DD_str = int(end_tok[2])
        end_local_dt = datetime(end_YYYY_str, end_MM_str, end_DD_str)
        utc_end_dt = pytz.utc.localize(end_local_dt)
    w_ret = pd.DataFrame()
    try:
<<<<<<< HEAD
        # Fetch WKT polygon and extract latitude and longitude
        wkt_polygon = fetchWKT(agstack_geoid)
        lat, lon = extractLatLonFromWKT(wkt_polygon)

        # Get S2 tokens for the given location
        s2_index__L5_list, L5_cids = get_s2_cellids_and_token_list(5, [lat], [lon])

        # Construct paths for relevant S2 partitions
        list_of_5_paths = [
            os.path.join(filePath, f's2_token_L5={s2_token}')
            for s2_token in s2_index__L5_list if os.path.exists(os.path.join(filePath, f's2_token_L5={s2_token}'))
        ]

        if not list_of_5_paths:
            print("No relevant partitions found for S2 tokens.")
            return empty_response()

        # Load datasets
        weather_datasets = [ds.dataset(path, format="parquet", partitioning="hive") for path in list_of_5_paths]

        w_all = pd.DataFrame()

        for weatherDataset in weather_datasets:
            # Create filter conditions for cross-year date ranges
            start_date_int = start_date_dt.year * 10000 + start_date_dt.month * 100 + start_date_dt.day
            end_date_int = end_date_dt.year * 10000 + end_date_dt.month * 100 + end_date_dt.day

            if end_date:  # Filter for a range of dates
                weather_df = weatherDataset.to_table(
                    filter=(
                        (ds.field('Year') * 10000 + ds.field('Month') * 100 + ds.field('Day')) >= start_date_int
                    ) & (
                        (ds.field('Year') * 10000 + ds.field('Month') * 100 + ds.field('Day')) <= end_date_int
                    )
                ).to_pandas()
            else:  # Filter for a single date
                weather_df = weatherDataset.to_table(
                    filter=(
                        (ds.field('Year') * 10000 + ds.field('Month') * 100 + ds.field('Day')) == start_date_int
                    )
                ).to_pandas()

            # Concatenate data
            w_all = pd.concat([w_all, weather_df], ignore_index=True)
        print("w_all=======",w_all)
        if not w_all.empty:
            # Fill missing values
            w_all = w_all.fillna(0)

            # Add 'aifstime_utc' if 'Year', 'Month', and 'Day' exist
            if {'Year', 'Month', 'Day'}.issubset(w_all.columns):
                w_all['aifstime_utc'] = pd.to_datetime(
                    w_all[['Year', 'Month', 'Day']].astype(str).agg('-'.join, axis=1)
                ).dt.tz_localize('UTC')

            # Convert 'TS' to datetime
            if 'TS' in w_all.columns:
                w_all['TS'] = pd.to_datetime(w_all['TS'])

            # Ensure numeric columns are properly typed
            numeric_cols = w_all.select_dtypes(include=[np.number]).columns
            w_all[numeric_cols] = w_all[numeric_cols].apply(pd.to_numeric, errors='coerce')

            # Group by 'TS' and calculate mean values
            w_ret = w_all.groupby('TS').mean(numeric_only=True).reset_index()

            # Retain 'Year', 'Month', and 'Day' as integers
            if {'Year', 'Month', 'Day'}.issubset(w_ret.columns):
=======


        s1_time_start = time.time()
        # Fetch WKT polygon and extract latitude and longitude
        wkt_polygon = fetchWKT(agstack_geoid)
        lat, lon = extractLatLonFromWKT(wkt_polygon)
        start_time = time.time()
        # Get the list of S2 indices and CIDs for the data point
        s2_index__L5_list, L5_cids = get_s2_cellids_and_token_list(5, [lat], [lon])
        
        list_of_5_paths = [filePath + 's2_token_L5=' + x for x in s2_index__L5_list
                        if os.path.exists(filePath + 's2_token_L5=' + x)]
        if not list_of_5_paths:
            return empty_response()
        
        weather_datasets = []
        for x in list_of_5_paths:
            weather_datasets.append(ds.dataset(x, format="parquet", partitioning="hive"))

        w_all = pd.DataFrame()
        for weatherDataset in weather_datasets:
            if end_date is not None:
                # Filter between start and end date
                weather_df = weatherDataset.to_table(
                    filter=(
                        (ds.field('Year') >= int(YYYY_str)) & (ds.field('Year') <= end_YYYY_str) &
                        (ds.field('Month') >= int(MM_str)) & (ds.field('Month') <= end_MM_str) &
                        (ds.field('Day') >= int(DD_str)) & (ds.field('Day') <= end_DD_str)
                    )
                ).to_pandas()
            elif end_date is None:
                # Filter for start date only
                weather_df = weatherDataset.to_table(
                    filter=(
                        (ds.field('Year') == int(YYYY_str)) &
                        (ds.field('Month') == int(MM_str)) &
                        (ds.field('Day') == int(DD_str))
                    )
                ).to_pandas()
            w_all = pd.concat([w_all, weather_df], ignore_index=True)
            print("========",w_all)
            # Ensure TS column is created
            if len(w_all) > 0:
                
                w_all = w_all.fillna(0)
                if 'Year' in w_all.columns and 'Month' in w_all.columns and 'Day' in w_all.columns:
                    # Create 'aifstime_utc' based on Year, Month, and Day
                    w_all['aifstime_utc'] = pd.to_datetime(
                        w_all[['Year', 'Month', 'Day']].astype(str).agg('-'.join, axis=1)
                    ).dt.tz_localize('UTC')
                
                # Convert 'time' column to datetime format if it isn't already
                w_all['TS'] = pd.to_datetime(w_all['TS'])
                # w_all['TS'] = pd.to_datetime(w_all[['Year', 'Month', 'Day']].astype(str).agg('-'.join, axis=1)).dt.tz_localize('UTC')


                # Ensure all numeric columns are converted to numeric type
                numeric_cols = w_all.select_dtypes(include=[np.number]).columns
                w_all[numeric_cols] = w_all[numeric_cols].apply(pd.to_numeric, errors='coerce')

                # Group by 'time' and calculate the mean for each group
                w_ret = w_all.groupby('TS').mean(numeric_only=True).reset_index()

                # Ensure 'Year', 'Month', and 'Day' are included in the response
>>>>>>> Rajat
                w_ret['Year'] = w_ret['Year'].astype(int)
                w_ret['Month'] = w_ret['Month'].astype(int)
                w_ret['Day'] = w_ret['Day'].astype(int)

    except Exception as e:
<<<<<<< HEAD
        print(f"Error: {e}")
        w_ret = empty_response()

=======

        print(f'Error--{e}')
        w_ret = empty_response()


    end_time = time.time()
    time_elapsed = (end_time - start_time)

>>>>>>> Rajat
    return w_ret

#cimis
def getWeatherFromCIMIS(agstack_geoid, start_date, end_date):
    filePath = '/network/CIMIS/HOURLY/PARQUET_S2/'

    
    # Parse start date
    tok = start_date.split('-')
    YYYY_str, MM_str, DD_str = tok[0], tok[1], tok[2]
    local_dt = datetime(int(YYYY_str), int(MM_str), int(DD_str))
    utc_dt = pytz.utc.localize(local_dt)

    # Parse end date if provided
    if end_date:
        end_tok = end_date.split('-')
        end_YYYY_str, end_MM_str, end_DD_str = int(end_tok[0]), int(end_tok[1]), int(end_tok[2])
        end_local_dt = datetime(end_YYYY_str, end_MM_str, end_DD_str)
        utc_end_dt = pytz.utc.localize(end_local_dt)

    w_ret = pd.DataFrame()
    
    try:
        wkt_polygon = fetchWKT(agstack_geoid)
        lat, lon = extractLatLonFromWKT(wkt_polygon)
        start_time = time.time()

        s2_index__L5_list, L5_cids = get_s2_cellids_and_token_list(5, [lat], [lon])
        print(f'token---{s2_index__L5_list}')

        list_of_5_paths = [filePath + 's2_token_L5=' + x for x in s2_index__L5_list
                           if os.path.exists(filePath + 's2_token_L5=' + x)]
        if not list_of_5_paths:
            return empty_response()

        weather_datasets = [ds.dataset(x, format="parquet", partitioning="hive") for x in list_of_5_paths]

        w_all = pd.DataFrame()
        for weatherDataset in weather_datasets:
            if end_date is not None:
<<<<<<< HEAD
                # Create integer representations for the start and end dates
                start_date_int = int(YYYY_str) * 10000 + int(MM_str) * 100 + int(DD_str)
                end_date_int = end_YYYY_str * 10000 + end_MM_str * 100 + end_DD_str

                # Filter between start and end date
                weather_df = weatherDataset.to_table(
                    filter=(
                        (ds.field('Year') * 10000 + ds.field('Month') * 100 + ds.field('Day')) >= start_date_int
                    ) & (
                        (ds.field('Year') * 10000 + ds.field('Month') * 100 + ds.field('Day')) <= end_date_int
                    )
                ).to_pandas()
            else:
                # Create integer representation for the start date
                start_date_int = int(YYYY_str) * 10000 + int(MM_str) * 100 + int(DD_str)

                # Filter for start date only
                weather_df = weatherDataset.to_table(
                    filter=(
                        (ds.field('Year') * 10000 + ds.field('Month') * 100 + ds.field('Day')) == start_date_int
                    )
                ).to_pandas()

=======
                weather_df = weatherDataset.to_table(
                    filter=(
                        (ds.field('Year') >= int(YYYY_str)) & 
                        (ds.field('Year') <= end_YYYY_str) &
                        (ds.field('Month') >= int(MM_str)) & 
                        (ds.field('Month') <= end_MM_str) &
                        (ds.field('Day') >= int(DD_str)) & 
                        (ds.field('Day') <= end_DD_str)
                    )
                ).to_pandas()
            else:
                weather_df = weatherDataset.to_table(
                    filter=(
                        (ds.field('Year') == int(YYYY_str)) &
                        (ds.field('Month') == int(MM_str)) &
                        (ds.field('Day') == int(DD_str))
                    )
                ).to_pandas()
>>>>>>> Rajat
            w_all = pd.concat([w_all, weather_df], ignore_index=True)

        if len(w_all) > 0:
            w_all = w_all.fillna(0)
<<<<<<< HEAD
            w_all['Date'] = pd.to_datetime(w_all['Date'])
=======
            # w_all['Date'] = pd.to_datetime(w_all['Date'])
            w_all['Date'] = pd.to_datetime(w_all['Date'], utc=True)
>>>>>>> Rajat

            numeric_cols = ['HlyAirTmp', 'HlyDewPnt', 'HlyEto', 'HlyNetRad', 'HlyAsceEto',
                            'HlyAsceEtr', 'HlyPrecip', 'HlyRelHum', 'HlyResWind',
                            'HlySoilTmp', 'HlySolRad', 'HlyVapPres', 'HlyWindDir',
                            'HlyWindSpd']
            
            for col in numeric_cols:
                if col in w_all.columns:
                    w_all[col] = pd.to_numeric(w_all[col], errors='coerce').fillna(0)

<<<<<<< HEAD
            # Group data by Date to get hourly averages
=======
>>>>>>> Rajat
            w_ret = w_all.groupby('Date', as_index=False).mean(numeric_only=True)
            w_ret['Year'] = w_ret['Date'].dt.year
            w_ret['Month'] = w_ret['Date'].dt.month
            w_ret['Day'] = w_ret['Date'].dt.day
<<<<<<< HEAD

    except Exception as e:
        print(e)
        w_ret = empty_response()
=======
            w_ret['Date'] = w_ret['Date'].dt.strftime("%a, %d %b %Y %H:%M:%S UTC")
    except Exception as e:
        print(e)
        w_ret = empty_response()

>>>>>>> Rajat
    return w_ret

def get_jrc_metadata():
    metadata = {
        "latitude": "float",
        "longitude": "float",
        "TMF_UAD": "integer",
        "s2_index__L13": "string",
        "s2_index__L18": "string"
    }
    return metadata

def get_jrc_description():
    metadata_jrc_description = {
        "latitude": "Latitude of the location",
        "longitude": "Longitude of the location",
        "TMF_UAD": "Tropical Moist Forest value",
        "s2_index__L13": "S2 token representing a geographic area at level 13 for spatial indexing.",
        "s2_index__L18": "S2 token representing a geographic area at level 18 for spatial indexing."
    }
    return metadata_jrc_description



# Function to get metadata for NDVI
def get_ndvi_metadata():
    """
    Returns metadata for the satellite data related to NDVI and boundary information.

    :return: Dictionary with metadata information
    """
    metadata = {
        "NDVI": "float (NDVI value for the geographical area, a value between -1 and 1)",
        "boundary_geoDataFrameDict": "GeoJSON FeatureCollection containing boundary geometries.",
        "coordinates": "list of lists (coordinates of the boundary geometry)",
        "coordinates_ndvi": "list of lists (coordinates of the NDVI geometry)",
        "geometry_type": "string (type of geometry, e.g., Polygon)",
        "geometry_type_ndvi": "string (type of geometry, e.g., Polygon)",
        "ndvi_img": "GeoJSON FeatureCollection containing NDVI values for specific geographical areas.",
        "properties": "dictionary (properties associated with the boundary, currently empty)"
    }

    return metadata


# Function to get metadata description for NDVI
def get_ndvi_metadata_description():
    """
    Returns a detailed description of the metadata fields for the satellite data.

    :return: Dictionary with metadata descriptions
    """
    metadata_sentinel_description = {
        "NDVI": "Normalized Difference Vegetation Index (NDVI) value for the defined geographical area. NDVI is a measure of vegetation health, with values ranging from -1 to 1, where positive values indicate healthy vegetation.",
        "boundary_geoDataFrameDict": "A GeoJSON FeatureCollection representing the boundary geometry of a specific geographical area.",
        "coordinates": "A list of lists that defines the coordinates of the boundary polygon. Each list represents a vertex of the polygon, and the first and last coordinates must be the same to close the polygon.",
        "coordinates_ndvi": "A list of lists that defines the coordinates of the NDVI polygon. Each list represents a vertex of the polygon, and the first and last coordinates must be the same to close the polygon.",
        "geometry_type": "The type of geometry used, which in this case is typically 'Polygon'.",
        "geometry_type_ndvi": "The type of geometry used, which is typically 'Polygon' for NDVI data.",
        "ndvi_img": "A GeoJSON FeatureCollection representing NDVI values within a specific geographical area.",
        "properties": "A dictionary containing properties associated with the boundary geometry. In this case, no properties are defined, but it could contain attributes like area, region name, etc."
    }
    return metadata_sentinel_description

def get_tmf_metadata():
    
    metadata = {
        "TMF_UAD": "float (TMF value for the geographical area, representing the Tropical Moist Forest)",
        "boundary_geoDataFrameDict": "GeoJSON FeatureCollection containing boundary geometries.",
        "coordinates": "list of lists (coordinates of the boundary geometry)",
        "coordinates_tmf": "list of lists (coordinates of the TMF geometry)",
        "geometry_type": "string (type of geometry, e.g., Polygon)",
        "geometry_type_tmf": "string (type of geometry, e.g., Polygon)",
        "tmf_img": "GeoJSON FeatureCollection containing TMF values for specific geographical areas.",
        "properties": "dictionary (properties associated with the boundary, currently empty)"
    }


    return metadata

def get_tmf_metadata_description():

    metadata_tmf_description = {
        "TMF_UAD": "Tropical Moist Forest (TMF) classification value for the defined geographical area. The value indicates forest cover status and can range from 0 to 3",
        "boundary_geoDataFrameDict": "A GeoJSON FeatureCollection representing the boundary geometry of a specific geographical area.",
        "coordinates": "A list of lists that defines the coordinates of the boundary polygon. Each list represents a vertex of the polygon, and the first and last coordinates must be the same to close the polygon.",
        "coordinates_tmf": "A list of lists that defines the coordinates of the TMF polygon. Each list represents a vertex of the polygon, and the first and last coordinates must be the same to close the polygon.",
        "geometry_type": "The type of geometry used, which in this case is typically 'Polygon'.",
        "geometry_type_tmf": "The type of geometry used, which is typically 'Polygon' for TMF data.",
        "tmf_img": "A GeoJSON FeatureCollection representing TMF values within a specific geographical area.",
        "properties": "A dictionary containing properties associated with the boundary geometry. In this case, no properties are defined, but it could contain attributes like area, region name, etc."
    }
    return metadata_tmf_description


# JRC
def getWeatherFromJRC(agstack_geoid):
    filePath = "/network/TMF_UAD/"
    
    # Fetch WKT polygon and extract lat/lon
    print("======agstack_geoid======",agstack_geoid)
    wkt_polygon = fetchWKT(agstack_geoid)
    lat, lon = extractLatLonFromWKT(wkt_polygon)
    print("=====lat lon======",lat,lon)

    # Get S2 token paths
    list_of_L3_paths = []
    list_of_L8_paths = []

    s2_index_L3_list, _ = get_s2_cellids_and_token_list(3, [lat], [lon])
    s2_index_L8_list, _ = get_s2_cellids_and_token_list(8, [lat], [lon])
    s2_index_L13_list,_ = get_s2_cellids_and_token_list(13, [lat], [lon])
    s2_index_L18_list,_ = get_s2_cellids_and_token_list(18, [lat], [lon])
    print(s2_index_L3_list)
    print(s2_index_L8_list)
    print(s2_index_L18_list)

    #Make L3 paths
    list_of_L3_paths = [filePath+'s2_index__L3='+x for x in s2_index_L3_list
                        if os.path.exists(filePath+'s2_index__L3='+x)]

    #Make L8 paths
    for p in list_of_L3_paths:
        list_of_L8_paths+= [p+'/s2_index__L8='+x for x in s2_index_L8_list 
                        if os.path.exists(p+'/s2_index__L8='+x)]

    if not list_of_L8_paths:
        return empty_response()
    
    print("=======list_of_L3_paths======",list_of_L3_paths)
    print("=======list_of_L8_paths======",list_of_L8_paths)
    datasets = []
    for x in list_of_L8_paths:
        datasets.append( ds.dataset(x,format="parquet", partitioning="hive") )

    #get the Data
    res_df = pd.DataFrame()
    for dtaset in datasets:
        df = dtaset.to_table().to_pandas()
        res_df = pd.concat([res_df, df], ignore_index=True)


    #screen for the L18
    res_df = res_df[res_df.s2_index__L18.isin(s2_index_L18_list)].reset_index(drop=True)

    # Format TMF_UAD column
    if "TMF_UAD" in res_df.columns:
        def format_tmf_uad(value):
            return int(value) if value >= 1 else value

        res_df["TMF_UAD"] = res_df["TMF_UAD"].apply(format_tmf_uad)

    return res_df
   


# def reload_gunicorn_workers():
#     print("Reloading Gunicorn workers...")
#     # Send the HUP signal to Gunicorn to gracefully reload the workers
#     os.kill(os.getppid(), signal.SIGHUP)  # Get Gunicorn's parent process ID and send SIGHUP

# def reload_data():
#     print("Reloading data...")
    
#     # Reload Gunicorn workers
#     reload_gunicorn_workers()



# def getWeatherFromNLDAS(agstack_geoid, dtStr):
#     filePath = '/mnt/md1/NLDAS/PARQUET_S2/'
#     tok=dtStr.split('-')
#     YYYY_str=tok[0]
#     MM_str=tok[1]
#     DD_str=tok[2]
    
#     configFile =  fieldDataConfigPath + agstack_geoid+'.json'
#     #NDays = (datetime.today() - datetime.strptime(dtStr,'%Y-%m-%d')).days
    
#     #Validate the date 
#     #1. Dt should be in the last 90 days from today
#     #2. Dt should be a valida date
    
#     s1_time_start = time.time()
#     with open(configFile, "r") as jsonfile:
#         fieldJSON = json.load(jsonfile)
#         field_geoid = fieldJSON['geoid']
#         field_wkt = fieldJSON['wkt']
#     fieldPoly = shapely.wkt.loads(field_wkt)
#     c=fieldPoly.centroid
#     lat=c.y
#     lon=c.x
        
    
#     lats=[lat]
#     lons=[lon]
#     start_time = time.time()
#     #get the list of S2 indeces and CIDs for the data point
#     s2_index__L8_list, L8_cids = get_s2_cellids_and_token_list(8, lats, lons)
#     list_of_L8_paths = [filePath+'s2_index__L8='+x for x in s2_index__L8_list 
#                         if os.path.exists(filePath+'s2_index__L8='+x)]
    
#     weather_datasets = []
#     for x in list_of_L8_paths:
#         weather_datasets.append( ds.dataset(x,format="parquet", partitioning="hive") )
    
#     #get the Data
#     w_all = pd.DataFrame()
#     for weatherDataset in weather_datasets:

#         yrStr=  YYYY_str
#         moStr = MM_str
#         dtStr = DD_str

#         YYYY_list=[yrStr]
#         MM_list=[moStr]
#         DD_list=[dtStr]
#         dt = datetime.strptime(yrStr+'-'+moStr+'-'+dtStr, '%Y-%m-%d')

#         weather_df = weatherDataset.to_table(
#             columns=['YYYY','MM','DD','ETo_AVG_IN','T_min', 'T_max', 'T_mean'],
#             filter=(
#                 ds.field('YYYY').isin(YYYY_list) &
#                 ds.field('MM').isin(MM_list) & 
#                 ds.field('DD').isin(DD_list)
#             )
#         ).to_pandas()

#         #get the DataFrame for the average of that day
#         weather_df['YYYY']=weather_df['YYYY'].astype(str)
#         weather_df['MM']=weather_df['MM'].astype(str).str.zfill(2)
#         weather_df['DD']=weather_df['DD'].astype(str).str.zfill(2)

#         weather_df['Date']=weather_df.YYYY+'-'+weather_df.MM+'-'+weather_df.DD
#         weather_df.Date = pd.to_datetime(weather_df.Date)

#         weather_df = weather_df.drop(columns=['YYYY','MM','DD'], axis=1)
#         w_df = pd.DataFrame(weather_df.mean(axis=0)).T
#         w_df['Date']=dt

#         cols = ['Date','ETo_AVG_IN','T_min', 'T_max', 'T_mean']
#         w_df = w_df[cols]
#         w_all = w_all.append(w_df, ignore_index=True)
    
#     if len(w_all)>0:
#         #take the average of the columns
#         w_ret = pd.DataFrame(w_all.groupby(['Date'])['ETo_AVG_IN','T_min', 'T_max', 'T_mean'].mean())
#         w_ret.reset_index(inplace=True)
#     else:
#         w_ret=pd.DataFrame()
    
#     end_time = time.time()
#     time_elapsed = (end_time - start_time)
    
#     return weather_df, w_ret, time_elapsed


#############

# def registerField(fieldWKT): 
#     """
#     Rules:
#     1. Registering only fields < 650 acres in area each
#     2. Find the list of L13s for this polygon
#     3. Find all other fields which have the one or more L13s from 
#     4. Overlay of L19 to provide overlay stats - if 95% overlap, then return the same id if one exists already
#     5. If doesn't exist, then create one and update both tables
#     """

#     fieldPoly = shapely.wkt.loads(fieldWKT)

#     #get the tile
#     #tileStr = S_gdf[S_gdf.intersects(fieldPoly)].TILE.iloc[0]
#     fieldDataConfigPath = fieldDataConfigRoot + 'geoids/'
#     #+ tileStr+ '/'
#     #if not os.path.exists(fieldDataConfigPath):
#     #    os.makedirs(fieldDataConfigPath)

#     #get the L13 indeces
#     s2_index__L13_list, L13_cids = get_s2_cellids_and_token_list_for_wkt(13, fieldWKT)
#     s2_index__L19_list, L19_cids = get_s2_cellids_and_token_list_for_wkt(19, fieldWKT)
#     s2_index__L20_list, L20_cids = get_s2_cellids_and_token_list_for_wkt(20, fieldWKT)

#     #generate the geoid
#     geoid = makeGeoID(s2_index__L20_list)

#     #read the lookup to get all the geoids for the L13s
#     geoidLookupFile = fieldDataConfigRoot+'geoid_lookup.json'
#     df = read_lookup(geoidLookupFile,13)

#     list_of_common_indices = list(set(s2_index__L13_list).intersection(set(df.s2_index__L13.tolist())))
#     list_of_existing_local_geoids = pd.unique(df[df.s2_index__L13.isin(list_of_common_indices)].geoid).tolist()

#     if len(list_of_existing_local_geoids)==0: #no match found; create it
#         print('\nFIELD NOT FOUND in registry ... creating and updating')
#         geoid = storeFieldInfo(fieldWKT, fieldDataConfigPath, fieldDataConfigRoot)

#         #update the registry for lookup
#         write_lookup(geoidLookupFile, s2_index__L15_list, geoid)
#         print('DONE!')

#         #construct the match_df dataframe
#         match_df = pd.DataFrame(columns=['new_geoid','existing_geoid', 'percent_overlap'])
#         match_df['new_geoid']=geoid
#         match_df['existing_geoid']=geoid
#         match_df['percent_overlap']=100.0
#     else:
#         match_df = getMatchingGeoids(geoid,s2_index__L19_list,list_of_existing_local_geoids, fieldDataConfigPath)
        
#     return match_df

# def getMatchingGeoids(geoid,s2_index__L19_list,list_of_existing_local_geoids, fieldDataConfigPath): #for each of the geoids find the match
#     df = pd.DataFrame(columns=['new_geoid','existing_geoid', 'percent_overlap'])
#     i=0
#     for gid in list_of_existing_local_geoids:
#         list_of_indices = getFieldIndices(fieldDataConfigPath, gid, 19)
#         percent_overlap = getPercentOverlap(list_of_indices,s2_index__L19_list)
#         if percent_overlap > 0:
#             df.loc[i,'new_geoid']=geoid
#             df.loc[i,'existing_geoid']=gid
#             df.loc[i,'percent_overlap']=percent_overlap
#             i=i+1
#     return df
                      
# def getFieldIndices(fieldDataConfigPath, geoid,resLevel):
#     configFile = fieldDataConfigPath + geoid+'.json'
#     with open(configFile, "r") as jsonfile:
#         fieldJSON = json.load(jsonfile)
    
#     #read the info
#     field_geoid = fieldJSON['geoid']
#     field_wkt = fieldJSON['wkt']
#     field_s2indices = fieldJSON['indices']
    
#     #get the s2 indices
#     s2_index__L8_list = json.loads(field_s2indices['s2_index__L8_list'])
#     s2_index__L10_list = json.loads(field_s2indices['s2_index__L10_list'])
#     s2_index__L13_list = json.loads(field_s2indices['s2_index__L13_list'])
#     s2_index__L15_list = json.loads(field_s2indices['s2_index__L15_list'])
#     s2_index__L19_list = json.loads(field_s2indices['s2_index__L19_list'])
#     s2_index__L20_list = json.loads(field_s2indices['s2_index__L20_list'])
    
#     list_of_indeces = eval('s2_index__L'+str(resLevel)+'_list')
    
#     return list_of_indeces

# def getPercentOverlap(list_of_indices,s2_index__L19_list):
#     list1 = list_of_indices
#     list2 = list(set(list_of_indices).intersection(set(s2_index__L19_list))) #index is in both lists
    
#     overlap_percent = round((len(list2)/len(list1))*100, 2)
    
#     return overlap_percent
    
# def getBoundaryCoverageNew(cids, poly, maxResColName):
#     s2_index__L19_list = []
#     p_gdf = gpd.GeoDataFrame()
#     idx=0
#     for cellid in cids:
#         c = s2.Cell(cellid)
#         level = cellid.level()
#         vertices = []
#         for i in range(0,4):
#             vertex = c.get_vertex(i)
#             latlng = s2.LatLng.from_point(vertex)
#             vertices.append((latlng.lng().degrees, latlng.lat().degrees))
#         geo = shapely.geometry.Polygon(vertices)
#         if poly.intersects(geo):
#             s2_index__L19_list.append(cellid.to_token())
#             p_gdf.loc[idx,maxResColName] = cellid.to_token()
#             p_gdf.loc[idx,'geometry']=geo
#         idx+=1
    
#     p_gdf.reset_index(drop=True, inplace=True)
#     return p_gdf


# def isValidPolygon(user_fieldWKT):
#     try:
#         poly = shapely.wkt.loads(user_fieldWKT)
#         if poly.geom_type == 'Polygon':
#             return True
#         else:
#             return False
#     except Exception as e:
#         print(e)
#         return False
    

# def makeGeoID(s2_index__L20_list):
#     s2_hash_method = 'S2 API + sha-256'
#     s2Tuple = tuple(s2_index__L20_list)
#     m = hashlib.sha256()
#     for s in s2Tuple:
#         m.update(s.encode())
#     geoid = m.hexdigest() #<-- geoid
#     return geoid
    
# #store config info on a field boundary
# def storeFieldInfo(user_fieldWKT, fieldDataConfigPath, fieldDataConfigRoot):
    
#     if isValidPolygon(user_fieldWKT):
#         #get the poly
#         poly = shapely.wkt.loads(user_fieldWKT)
#         polyShape = shapely.geometry.mapping(poly)
        
#         #Make the GeoDataFrame
#         acct_gdf = gpd.GeoDataFrame()
#         acct_gdf.loc[0,'name']='field'
#         acct_gdf.loc[0,'geometry']=shapely.wkt.loads(user_fieldWKT)
#         acct_gdf.crs = {'init':'epsg:4326'}
        
#         #get the indeces
#         print('Getting indeces for field boundary...')
#         s2_index__L20_list, L20_cids = get_s2_cellids_and_token_list_for_wkt(20, user_fieldWKT)
#         s2_index__L19_list, L19_cids = get_s2_cellids_and_token_list_for_wkt(19, user_fieldWKT)
#         s2_index__L15_list, L15_cids = get_s2_cellids_and_token_list_for_wkt(15, user_fieldWKT)
#         s2_index__L13_list, L13_cids = get_s2_cellids_and_token_list_for_wkt(13, user_fieldWKT)
#         s2_index__L10_list, L10_cids = get_s2_cellids_and_token_list_for_wkt(10, user_fieldWKT)
#         s2_index__L8_list, L8_cids = get_s2_cellids_and_token_list_for_wkt(8, user_fieldWKT)
#         s2_index__L5_list, L5_cids = get_s2_cellids_and_token_list_for_wkt(5, user_fieldWKT)
#         s2_index__L3_list, L3_cids = get_s2_cellids_and_token_list_for_wkt(3, user_fieldWKT)
#         print('\tdone!')
        
#         #s2_max_resolution
#         print('Getting p_gdf L20 coverage...')
#         p20_gdf = getBoundaryCoverageNew(L20_cids, poly, 's2_index__L20')
#         print('\tdone!')
        
#         #s2_L19_resolution
#         print('Getting p_gdf L19 coverage...')
#         p19_gdf = getBoundaryCoverageNew(L19_cids, poly, 's2_index__L19')
#         print('\tdone!')
        
#         #s2_L15_resolution
#         print('Getting p_gdf L15 coverage...')
#         p15_gdf = getBoundaryCoverageNew(L15_cids, poly, 's2_index__L15')
#         print('\tdone!')
        
#         #s2_L13_resolution
#         print('Getting p_gdf L13 coverage...')
#         p13_gdf = getBoundaryCoverageNew(L13_cids, poly, 's2_index__L13')
#         print('\tdone!')
        
#         #s2_L10_resolution
#         print('Getting p_gdf L10 coverage...')
#         p10_gdf = getBoundaryCoverageNew(L10_cids, poly, 's2_index__L10')
#         print('\tdone!')
        
#         #s2_L8_resolution
#         print('Getting p_gdf L8 coverage...')
#         p8_gdf = getBoundaryCoverageNew(L8_cids, poly, 's2_index__L8')
#         print('\tdone!')
        
#         #s2_L5_resolution
#         print('Getting p_gdf L5 coverage...')
#         p5_gdf = getBoundaryCoverageNew(L5_cids, poly, 's2_index__L5')
#         print('\tdone!')
        
#         #s2_L3_resolution
#         print('Getting p_gdf L3 coverage...')
#         p3_gdf = getBoundaryCoverageNew(L3_cids, poly, 's2_index__L3')
#         print('\tdone!')
        
        
#         #get the sentinel tile
#         print('getting tile...')
#         tile_gdf = S_gdf[S_gdf.intersects(poly)]
#         tile_gdf.reset_index(drop=True, inplace=True)
#         print('\tdone!')
        
#         #get the hashkey
#         geoid = makeGeoID(s2_index__L20_list)
        
        
        
#         fieldJSON = {
#             'geoid': geoid,
#             'wkt': user_fieldWKT,
#             's2_L20_index_gdfJSON': p20_gdf.to_json(drop_id=True),
#             's2_L19_index_gdfJSON': p19_gdf.to_json(drop_id=True),
#             's2_L15_index_gdfJSON': p15_gdf.to_json(drop_id=True),
#             's2_L13_index_gdfJSON': p13_gdf.to_json(drop_id=True),
#             's2_L10_index_gdfJSON': p10_gdf.to_json(drop_id=True),
#             's2_L8_index_gdfJSON': p8_gdf.to_json(drop_id=True),
#             's2_L5_index_gdfJSON': p5_gdf.to_json(drop_id=True),
#             's2_L3_index_gdfJSON': p3_gdf.to_json(drop_id=True),
#             'boundary_gdfJSON': acct_gdf.to_json(drop_id=True),
#             'sentinel-tile-list': tile_gdf.TILE.tolist(),
#             'indices': {
#                 's2_index__L3_list': json.dumps(s2_index__L3_list),
#                 's2_index__L5_list': json.dumps(s2_index__L5_list),
#                 's2_index__L8_list': json.dumps(s2_index__L8_list),
#                 's2_index__L10_list': json.dumps(s2_index__L10_list),
#                 's2_index__L13_list': json.dumps(s2_index__L13_list),
#                 's2_index__L15_list': json.dumps(s2_index__L15_list),
#                 's2_index__L19_list': json.dumps(s2_index__L19_list),
#                 's2_index__L20_list': json.dumps(s2_index__L20_list),
#             },
#         }
        
#         print('Writing JSON ...')
#         configFile = fieldDataConfigPath + geoid+'.json'
#         #Write the config file to disk
#         with open(configFile, "w") as jsonfile:
#             json.dump(fieldJSON, jsonfile)
#         print('\tdone!')
        
        
#         #append to the dictionary of geoids L13a
#         geoidLookupFile = fieldDataConfigRoot+'geoid_lookup.json'
#         write_lookup(geoidLookupFile, s2_index__L13_list, geoid)
        
#         return geoid
    
# def read_lookup(geoidLookupFile, resLevel):
#     colName = 's2_index__L'+str(resLevel)
#     df = pd.DataFrame(columns=[colName,'geoid'])
#     #read the file
#     try:
#         with open(geoidLookupFile, "r") as geoidfile:
#             df = pd.DataFrame.from_dict(json.loads(json.load(geoidfile)))
#     except:
#         return df
#     return df

# def write_lookup(geoidLookupFile, s2_index_list, geoid, resLevel):
#     df = read_lookup(geoidLookupFile, resLevel)
#     i=len(df)
#     colName = 's2_index__L'+str(resLevel)
#     for idx in s2_index_list:
#         df.loc[i,colName] = idx
#         df.loc[i,'geoid']=geoid
#         i=i+1
    
#     #drop all the duplicate rows
#     df = df.drop_duplicates(subset=[colName,'geoid'])
    
#     #Write the file
#     data = df.to_json()
#     with open(geoidLookupFile, 'w') as f:
#         json.dump(data, f)
#     return

################  ROUTES
# @app.route('/')
#@app.route('/index')
#@app.route('/home')
#def index():
#    return render_template('index.html')
#
cache_nldas = {}
@app.route('/')
#def hello_world():
 #   return 'Hello World!'
@app.route('/getNLDAS')
def getNLDASData():
    geoid = request.args['geoid']
    # s2_token = ["8094c20c","8094e9e4"]
    # date = request.args.get('date', '')

    start_date = request.args.get('start_date', '')
    end_date = request.args.get('end_date',None)
    
    if end_date is not None :

        weather_df = getWeatherFromNLDAS(geoid, start_date, end_date)
    else:
        weather_df = getWeatherFromNLDAS(geoid, start_date, end_date=None)
    
    try:
        weather_df.reset_index(inplace=True)
    except:

        pass

        

    # Convert DataFrame to a dictionary
    json_data = weather_df.to_dict(orient='records')

    response = {
        "data": json_data,
        "metadata": meta_data_response(),
        "metadata-description":meta_data_with_description()
    }

    return jsonify(response)


cache_ncep = {}

@app.route('/getNCEP')
def getNCEPWeatherData():
    geoid = request.args['geoid']
    start_date = request.args.get('start_date', '')
    end_date = request.args.get('end_date', None)

   
    # Fetch data from NCEP
    if end_date is not None:
        weather_df = getWeatherFromNCEP(geoid, start_date, end_date)
    else:
        weather_df = getWeatherFromNCEP(geoid, start_date)

    try:
        weather_df.reset_index(inplace=True)
    except:
        pass

    # Convert DataFrame to JSON
    json_data = weather_df.to_dict(orient='records')


    response = {
        "data": json_data,
        "metadata": meta_data_response(),
        "metadata-description": meta_data_with_description()
    }

    return jsonify(response)


@app.route('/getNCEPByDate')
def getNCEP():
    geoid = request.args['geoid']
    start_date = request.args.get('start_date', '')
    end_date = request.args.get('end_date', None)

    if end_date:
        weather_data = getWeatherFromNCEP(geoid, start_date, end_date)
    else:
        weather_data = getWeatherFromNCEP(geoid, start_date)

    # Ensure that the result from getWeatherFromNCEP is already a dictionary or JSON-serializable format
    response = {
        "data": weather_data if isinstance(weather_data, dict) else weather_data.to_dict(),
        "metadata": meta_data_response(),
        "metadata-description": meta_data_with_description()
    }

    return jsonify(response)



@app.route('/getWeatherFromNCEPForecast', methods=['GET'])
def getNCEPWeatherForecastData():
    geoid = request.args.get('geoid')
    start_date = request.args.get('start_date')
    end_date = request.args.get('end_date')

    if not geoid or not start_date or not end_date:
        return jsonify({"error": "Please provide geoid, start_date, and end_date"}), 400

    # Start timing for cache access
    start_time_cache_access = time.time()
    cache_key_ncep = ('NCEP', geoid, start_date, end_date)

    if cache_key_ncep in cache_ncep:
        # Log cache access time
        cache_access_duration = time.time() - start_time_cache_access
        logging.info(f"Returning cached data for geoid: {geoid} for the range {start_date} to {end_date} took {cache_access_duration:.2f} seconds")
        return jsonify({"data": cache_ncep[cache_key_ncep], "metadata": meta_data_response()})

    # Data not in cache, fetch from NCEP
    weather_df = fetchNCEPWeatherForecast(start_date, end_date, geoid)

    # Convert DataFrame to a dictionary and reset index if needed
    try:
        weather_df.reset_index(inplace=True)
    except Exception as e:
        logging.warning(f"Resetting index failed for dataframe with geoid: {geoid}, error: {e}")

    json_data = weather_df.to_dict(orient='records')

    # Cache the response
    cache_ncep[cache_key_ncep] = json_data

    # Log the completion of the request
    logging.info(f"Data for geoid: {geoid} from {start_date} to {end_date} fetched from NCEP and cached")

    response = {
        "data": json_data,
        "metadata": meta_data_response()
    }

    return jsonify(response)


#NOAA
cache_noaa = {}
@app.route('/getNOAA')
def getNOAAWeatherData():
    geoid = request.args['geoid']
    start_date = request.args.get('start_date', '')
    end_date = request.args.get('end_date',None)   

    if end_date is not None :

        weather_df = getWeatherFromNOAA(geoid, start_date, end_date)
    else:
        weather_df = getWeatherFromNOAA(geoid, start_date, end_date=None)
    
    try:
        weather_df.reset_index(inplace=True)
    except:
        pass

    # Convert DataFrame to a dictionary
    json_data = weather_df.to_dict(orient='records')


    response = {
        "data": json_data,
        "metadata": get_noaa_metadata(),
        "metadata-description":get_noaa_description()
    }


    return jsonify(response)

#NOAA_FORECASTED
@app.route('/getNOAAForecast')
def getNOAAFORECASTEDWeatherData():
    geoid = request.args['geoid']
    start_date = request.args.get('start_date', '')
    end_date = request.args.get('end_date',None)

    if end_date is not None :

        weather_df = getWeatherFromNOAAFORECASTED(geoid, start_date, end_date)
    else:
        weather_df = getWeatherFromNOAAFORECASTED(geoid, start_date, end_date=None)
    
    try:
        weather_df.reset_index(inplace=True)
    except:

        pass

    
    # Convert DataFrame to a dictionary
    json_data = weather_df.to_dict(orient='records')


    response = {
        "data": json_data,
        "metadata": get_noaa_forecast_metadata(),
        "metadata-description":get_noaa_forecast_description()
    }
    return jsonify(response)


#GHCND
cache_ghcnd = {}
@app.route('/getGHCND')
def getGHCNDWeatherData():
    geoid = request.args['geoid']
    start_date = request.args.get('start_date', '')
    end_date = request.args.get('end_date', None)


    if end_date is not None:
        weather_df = getWeatherFromGHCND(geoid, start_date,end_date)
    else:
        weather_df = getWeatherFromGHCND(geoid, start_date,end_date)
    
    try:
        weather_df.reset_index(inplace=True)
    except:
        pass

    # Convert DataFrame to a dictionary
    json_data = weather_df.to_dict(orient='records')

    response = {
        "data": json_data,
        "metadata": get_GHCND_metdata(),
        "metadata-description":get_GHCND_description()
    }

    return jsonify(response)


#AUS_FORECASTED
@app.route('/getAUSFORECASTED')
def getAUSFORECASTEDWeatherData():
    geoid = request.args['geoid']
    # date = request.args.get('date', '')

    # geoid = request.args['geoid']
    start_date = request.args.get('start_date', '')
    end_date = request.args.get('end_date',None)

    # Data not in cache, fetch from AUS FORECASTED
    # weather_df = getWeatherFromAUSFORECASTED(geoid, date)
    
    if end_date is not None:

        weather_df = getWeatherFromAUSFORECASTED(geoid, start_date,end_date)
    else:
        weather_df = getWeatherFromAUSFORECASTED(geoid, start_date,end_date=None)
    
    
    try:
        weather_df.reset_index(inplace=True)
    except:
        # logging.warning(f"Resetting index failed for dataframe with geoid: {geoid}")
        pass

    # Convert DataFrame to a dictionary
    json_data = weather_df.to_dict(orient='records')


    response = {
        "data": json_data,
        "metadata": get_aus_forcasted_metadata(),
        "metadata-description":get_aus_forcasted_description()
    }

    return jsonify(response)


#AUS
cache_aus = {}
@app.route('/getAUS')
def getAUSWeatherData():
    geoid = request.args['geoid']
    start_date = request.args.get('start_date', '')
    end_date = request.args.get('end_date',None)

    # Start timing for cache access
    # start_time_cache_access = time.time()
    # cache_key_aus = ('AUS',geoid, date)
    # if cache_key_aus in cache_aus:
    #     # Log cache access time
    #     cache_access_duration = time.time() - start_time_cache_access
    #     # logging.info(f"Returning cached data for geoid: {geoid} took {cache_access_duration:.2f} seconds")
    #     return jsonify(cache_aus[cache_key_aus])

    
    if end_date is not None:

        weather_df = getWeatherFromAUS(geoid, start_date,end_date)
    else:
        weather_df = getWeatherFromAUS(geoid, start_date,end_date=None)
    
    
    try:
        weather_df.reset_index(inplace=True)
    except:
        # logging.warning(f"Resetting index failed for dataframe with geoid: {geoid}")
        pass

    # Convert DataFrame to a dictionary
    json_data = weather_df.to_dict(orient='records')

    # Cache the response
    # cache_aus[cache_key_aus] = json_data

    # Log the completion of the request
    logging.info(f"Data for geoid: {geoid} fetched from AUS and cached")

    response = {
        "data": json_data,
        "metadata": get_aus_metadata(),
        "metadata-description":get_aus_description()
    }

    return jsonify(response)


@app.route('/getSatelite')
def getSateliteData():
    
    geoid = request.args['geoid']
    start_date = request.args.get('start_date', '')
    end_date = request.args.get('end_date',None)

    
        
    satelite_df = getWeatherFromSatelite(geoid,start_date,end_date)
    
    json_data = satelite_df.to_dict(orient='records')
    
    
    response = {
        "data": json_data,
        "metadata": satelite_metadata(),
        "metadata-description":get_sentinel_data_description()
    }
    
    return jsonify(response)

@app.route('/getSateliteStats')
def getSateliteStatsRoute():
    geoid = request.args['geoid']
    start_date = request.args.get('start_date', '')
    end_date = request.args.get('end_date', None)  # This retrieves the end_date correctly

    # Pass end_date directly without overriding it to None
    satelite_stats = getSateliteStatsFn(geoid, start_date, end_date)

    json_data = satelite_stats.to_dict(orient='records')

    response = {
        "NDVI": json_data,
        "metadata":get_satelite_stats_metadata(),
        "metadata-description":get_satelite_desc()
    }

    return response


#CIMIS
cache_cimis = {}
@app.route('/getCIMIS')
def getCIMISWeatherData():
    geoid = request.args['geoid']
    start_date = request.args.get('start_date', '')
    end_date = request.args.get('end_date', None)

    weather_df = getWeatherFromCIMIS(geoid, start_date,end_date)
    
    try:
        weather_df.reset_index(inplace=True)
    except:
        pass

    # Convert DataFrame to a dictionary
    json_data = weather_df.to_dict(orient='records')

    response = {
        "data": json_data,
        "metadata": get_cimis_metadata(),
        "metadata-description":get_cimis_description()
    }

    return jsonify(response)

@app.route('/getEto')
def getEtoFromWeather():
    agstack_geoid = request.args['geoid']
    start_date = request.args['start_date']
    end_date = request.args.get('end_date', None)

    # Call the function to get the average ETo value
    avg_eto__df = getEtoFromWeatherData(agstack_geoid, start_date, end_date)
    
    json_data = avg_eto__df.to_dict(orient='records')
    # Return the average ETo value in the response
    response = {
                'data': json_data,
                'metadata':eto_metadata(),
                'metadata-description':get_eto_data_description()
                
                }
    
    return jsonify(response)

#JRC
cache_jrc = {}
@app.route('/getJRC')
def getJRCWeatherData():
    geoid = request.args['geoid']

    weather_df = getWeatherFromJRC(geoid)

    # Convert DataFrame to a dictionary
    json_data = weather_df.to_dict(orient='records')

    response = {
        "data": json_data,
        "metadata": get_jrc_metadata(),
        "metadata-description":get_jrc_description()
    }

    return jsonify(response)


@app.route('/ticket-response')
def getTicketresponse():
    
    data = api_ticket_response()

    return data


@app.route('/')
def home():
    return 'hello world'

# @app.route('/getNDVIStats')
# def getNDVIStats():
#     startDtStr =  request.args.get('startDt', '')
#     endDtStr =  request.args.get('endDt', '')
#     agstack_geoid =  request.args.get('geoid', '')
#     #call the function
#     dates_df = getNDVIStatsFromGeoid(agstack_geoid, startDtStr, endDtStr)

#     return dates_df.to_json()

@app.route('/getNDVIImg')
def getNDVIImg():
    """
    Endpoint to fetch the NDVI image and boundary GeoJSON data.
    """
    date = request.args.get('date')
    agstack_geoid =  request.args.get('geoid', '')
    
    # Assuming `getNDVIgdfFromGeoid` is a function that returns GeoDataFrames for NDVI and boundary data
    final_res_gdf, bound_gdf = getNDVIgdfFromGeoid(agstack_geoid, date)
    
    imgDict = {
        'date': date,
        'ndvi_img': json.loads(final_res_gdf.to_json()),  # GeoDataFrame to GeoJSON object
        'boundary_geoDataFrameDict': json.loads(bound_gdf.to_json()),  # GeoDataFrame to GeoJSON object
        'metadata': get_ndvi_metadata(),  # Call metadata function
        'metadata-description': get_ndvi_metadata_description()  # Call metadata description function
    }
    
    return jsonify(imgDict)  # Return the response as JSON

@app.route('/getTMFImg')
def getTMFImg():
    """
    Endpoint to fetch the NDVI image and boundary GeoJSON data.
    """
    # date = request.args.get('date')
    agstack_geoid =  request.args.get('geoid', '')
    
    # Assuming `getNDVIgdfFromGeoid` is a function that returns GeoDataFrames for NDVI and boundary data
    final_res_gdf, bound_gdf = getTMFgdfFromGeoid(agstack_geoid)
    
    imgDict = {
        # 'date': date,
        'tmf_img': json.loads(final_res_gdf.to_json()),  # GeoDataFrame to GeoJSON object
        'boundary_geoDataFrameDict': json.loads(bound_gdf.to_json()),  # GeoDataFrame to GeoJSON object
        'metadata': get_tmf_metadata(),  # Call metadata function
        'metadata-description': get_tmf_metadata_description()  # Call metadata description function
    }
    
    return jsonify(imgDict)  # Return the response as JSON



@app.route('/getNDVIDatesForGeoid')
def getNDVIDatesForGeoId():
    geoid = request.args.get('geoid')
    if not geoid:
        return jsonify({"error": "Missing 'geoid' parameter"}), 400

    result = getAvailableDatesForGeoID(geoid)
    if "error" in result:
        return jsonify(result), 500

    return jsonify(result), 200

# @app.route('/regField')
# def regField():
#     fieldWKT = "POLYGON((0 0, 1 0, 1 1, 0 1, 0 0))"
#     fieldWKT =  request.args.get('fieldWKT', '')
#     #print(fieldWKT)
#     match_df = registerField(fieldWKT)
#     matchDict = {}
#     matchDict['type']='Feature'
#     matchDict['geoid_matches']=match_df.to_json()

#     return matchDict

# Directory watching and hot-reloading setup
base_path_nldas_reload = '/mnt/md1/NLDAS/PARQUET_S2/s2_tokens_L5/s2_tokens_L7/s2_tokens_L9/Year=2024/'
base_path_ncep_reload = '/mnt/md1/NCEP/PARQUET_S2/s2_tokens_L5/s2_tokens_L7/s2_tokens_L9/Year=2024/'
extra_dirs = [base_path_nldas_reload, base_path_ncep_reload]
extra_files = extra_dirs[:]

for extra_dir in extra_dirs:
    for dirname, dirs, files in walk(extra_dir):
        for filename in files:
            file_path = path.join(dirname, filename)
            if path.isfile(file_path):
                extra_files.append(file_path)

if __name__ == '__main__':
    # extra_files = [updated_data_available_file,]

    app.run(debug=True,extra_files=extra_files,port=5000)



"""
@app.route('/getNDVI')
def getNDVI():
    user_fieldWKT =  request.args.get('wkt', '')



    s2_max_resolution = resLevel
    #get the field boundary gdf
    acct_gdf = gpd.GeoDataFrame()
    acct_gdf.loc[0,'name']='field'
    acct_gdf.loc[0,'geometry']=shapely.wkt.loads(user_fieldWKT)

    #print(user_fieldWKT)

    if isValidPolygon(user_fieldWKT):
        #try:

        poly = shapely.wkt.loads(user_fieldWKT)
        polyShape = shapely.geometry.mapping(poly)
        

        s2_index__maxRes_list, maxRes_cids = get_s2_cellids_and_token_list_for_wkt(s2_max_resolution, user_fieldWKT)
        s2_index__L15_list, L15_cids = get_s2_cellids_and_token_list_for_wkt(15, user_fieldWKT)
        s2_index__L10_list, L10_cids = get_s2_cellids_and_token_list_for_wkt(10, user_fieldWKT)
        s2_index__L8_list, L8_cids = get_s2_cellids_and_token_list_for_wkt(8, user_fieldWKT)

        list_of_paths = [pa_filePath+'s2_index__L8='+x for x in s2_index__L8_list]
        datasets = []
        for p in list_of_paths:
            datasets.append( ds.dataset(p,format="parquet", partitioning="hive") )

        #s2_max_resolution=19
        maxResColName = 's2_index__L'+str(s2_max_resolution)

        p_gdf = getBoundaryCoverageNew(maxRes_cids, poly, maxResColName)

        p_df, t_elapsed = getnewNDVIFromParquetForDataSets_maxRes(s2_max_resolution, datasets,300, s2_index__L10_list, s2_index__L15_list, s2_index__maxRes_list)
        #print('Got p_df in '+str(t_elapsed))

        p_gdf = getFieldOutline(user_fieldWKT, s2_max_resolution, maxResColName)

        list_of_dts = pd.unique(p_df['YY_MM_DD']).tolist()
        if len(list_of_dts)>0:
            list_of_dts.sort()

            min_dtStr, max_dtStr, no_dtStrs = list_of_dts[0], list_of_dts[-1], len(list_of_dts)

            dtStr = max_dtStr

            r_gdf, r_stats = getNDVIAndStatsForDtStr_maxRes(s2_max_resolution, p_df, p_gdf, dtStr)

            polyDict = {}
            polyDict['type']='Feature'
            polyDict['properties']={
                'index': 'NDVI',
                'mean' :  r_stats[0],
                'min' :  r_stats[1],
                'p10' :  r_stats[2],
                'p25' :  r_stats[3],
                'p50' :  r_stats[4],
                'p75' :  r_stats[5],
                'p90' :  r_stats[6],
                'max' :  r_stats[7],
            }
            polyDict['date_UTC']=dtStr
            polyDict['extent']=json.dumps(getExtent(poly))
            polyDict['crs']={'init':'epsg:4326'}
            polyDict['geometry']=polyShape
            polyDict['extent']=json.dumps(getExtent(poly))
            
            #return [dtStr, r_gdf, p_gdf, polyDict]
            return json.loads(json.dumps(polyDict))
        else:
            return'empty dataset' 
    else:
        #print('Invalid Polygon')
        return 'Invalid Polygon'
"""
