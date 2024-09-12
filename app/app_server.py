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
from datetime import date, timedelta, datetime
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
import cProfile
import pstats
from io import StringIO
from functools import wraps
from line_profiler import LineProfiler
from flask_cors import CORS
import concurrent.futures
import pyarrow.parquet as pq
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler
from watchdog.observers import Observer

from functools import wraps
import multiprocessing
import signal
import pytz 
warnings.filterwarnings('ignore')


################ GLOBALS
#global res_gdf
#global pa_filePath
profiler = LineProfiler()

pa_filePath = '/mnt/md1/SENTINEL/PARQUET_S2_L20/'
tileStr = '11SKA'
fieldDataConfigRoot = '/home/rnaura/terrapipe/config/'
fieldDataOutputRoot = '/home/rnaura/terrapipe/OUTPUT/'

geoidLookupFile = fieldDataConfigRoot+'geoid_lookup.json'

fieldDataConfigPath = fieldDataConfigRoot + tileStr+ '/'
resLevel = 20
s2_max_resolution = resLevel
maxResColName = 's2_index__L'+str(resLevel)

# logging.basicConfig(
#     filename='weather_data_processing.log',  # Log file name
#     level=logging.INFO,  # Log level (INFO, DEBUG, ERROR, etc.)
#     format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',  # Log message format
#     datefmt='%Y-%m-%d %H:%M:%S'  # Date format
# )# cache = Cache(app)




app = Flask(__name__)
CORS(app)

app.config["TEMPLATES_AUTO_RELOAD"] = True

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
# def getLocalNDVIDatasets(s2_index__L8_list, s2_index__L10_list, pa_filePath):

#     list_of_paths = []
#     list_of_paths_L8 = [pa_filePath+'s2_index__L8='+x for x in s2_index__L8_list]
#     for p in list_of_paths_L8:
#         list_of_paths = list_of_paths + [p+'/s2_index__L10='+x for x in s2_index__L10_list]
    
#     datasets = []
#     for p in list_of_paths:
#         #print(p)
#         if os.path.exists(p):
#             datasets.append( ds.dataset(p,format="parquet", partitioning="hive") )
#         else:
#             continue
        
#     return datasets

# def getNDVIFromParquetForDataSetsAndDate(datasets, dtStr, s2_index__L10_list, s2_index__L15_list, s2_index__L20_list):
#     start_time = time.time()
#     subset_dirs = []
#     endDt = datetime.now()
    
    
#     endDt =  datetime.strptime(dtStr, '%Y-%m-%d')
#     startDt = endDt - timedelta(days=10)
#     NDays = (endDt - startDt).days
    
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
    
#     #get the closest dt and screen for it

#     list_of_dates = pd.unique(p_df.YY_MM_DD).tolist()
#     list_of_dates.sort()
#     closest_dtStr = list_of_dates[-1]
    
#     p_df = p_df[p_df.YY_MM_DD==closest_dtStr].reset_index(drop=True)
    
#     if len(p_df)==0:
#         print('Empty Dataset')
#     else:
#         p_df.reset_index(drop=True, inplace=True)
        
#     end_time = time.time()
#     time_elapsed = (end_time - start_time)
    
   

#     return p_df, time_elapsed


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
#     s2_index__L13_list = json.loads(field_s2indices['s2_index__L13_list'])
#     s2_index__L15_list = json.loads(field_s2indices['s2_index__L15_list'])
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


            

""" Find  s2 tokens code  for different levels"""

def get_level_11_tokens(level_13_tokens):
    level_11_tokens = []
    for token in level_13_tokens:
        cell_id = s2sphere.CellId.from_token(token)
        parent_cell_id = cell_id.parent(11)  # Get parent cell at level 11
        level_11_token = parent_cell_id.to_token()
        level_11_tokens.append(level_11_token)
        print(f'level_11_tokens {level_11_tokens}')
    return level_11_tokens

def get_level_20_tokens(level_9_tokens):
    # Using list comprehension for efficiency
    level_20_tokens = [
        s2sphere.CellId.from_token(token).parent(20).to_token()
        for token in level_9_tokens
    ]
    
    # Print outside the loop to avoid multiple print calls
    print(f'level_20_tokens {level_20_tokens}')
    
    return level_20_tokens

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
    print(f'level_5_tokens {level_5_tokens}')
    return level_5_tokens 

base_path_ncep = '/mnt/md1/NCEP/PARQUET_S2/'
base_path_nldas = '/mnt/md1/NLDAS/PARQUETE_S2/'
base_path_noaa = '/home/rnaura/mnt/md1/NOAA/DAILY/PROCESSED/S2_GEOJSON/PARQUET_S2/'
# base_path_noaa = '/home/rajat/Downloads/Rnaura/mnt/md1/NOAA/DAILY/PROCESSED/S2_GEOJSON/PARQUET_S2/'
# base_path_noaa_forecasted = '/home/rajat/Downloads/Rnaura/mnt/md1/NOAA/DAILY/FORECASTED/PARQUET_S2/'
base_path_noaa_forecasted = '/home/rnaura/mnt/md1/NOAA/DAILY/FORECASTED/PARQUET_S2/'
# base_path_ghcnd = '/home/rajat/Downloads/Rnaura/mnt/md1/GHCND/PARQUET_S2/'
base_path_ghcnd = '/home/rnaura/mnt/md1/GHCND/PARQUET_S2/'
# base_path_aus_forecasted = '/home/rajat/Downloads/Rnaura_Work/mnt/md1/AUS/DAILY/FORECASTED/PARQUET_S2/'
base_path_aus_forecasted = '/home/rnaura/mnt/md1/AUS/DAILY/FORECASTED/PARQUET_S2/'

# base_path_aus = '/home/rajat/Downloads/Rnaura_Work/mnt/md1/AUS/DAILY/PROCESSED/PARQUET_S2/'
base_path_aus = '/home/rnaura/mnt/md1/AUS/DAILY/PROCESSED/PARQUET_S2/'


# def reload_data_nldas(YYYY_str, MM_str):
#     cached_weather_files_ncep = {}
#     try:
#         start_time = time.time()
#         month_path = os.path.join(base_path_nldas, 's2_tokens_L5', 's2_tokens_L7', 's2_tokens_L9', f'Year={YYYY_str}', f'Month={MM_str}')
        
#         for day_dir in os.listdir(month_path):
#             day_path = os.path.join(month_path, day_dir)
#             if os.path.isdir(day_path):
#                 day_files = [os.path.join(day_path, f) for f in os.listdir(day_path) if f.endswith('.parquet')]
#                 if day_files:
#                     cached_weather_files_ncep[day_dir] = day_files

#         print(f"Directories reloaded in {time.time() - start_time:.4f} seconds")
#     except Exception as e:
#         print(f"Failed to reload directories: {e}")
#     return  cached_weather_files_ncep
# Configure Flask-Caching
cache = Cache(app, config={'CACHE_TYPE': 'simple'})

base_path_nldas = '/mnt/md1/NLDAS/PARQUETE_S2/s2_tokens_L5/s2_tokens_L7/s2_tokens_L9/'


# def load_all_data():
#     global cached_weather_df_nldas
#     all_tables = []

#     print("Starting to load all data into memory...")

#     try:
#         for year in os.listdir(base_path_nldas):
#             year_path = os.path.join(base_path_nldas, year)
#             if not year.startswith('Year='):
#                 print(f"Skipping {year_path}: Not a valid year directory.")
#                 continue

#             for month in os.listdir(year_path):
#                 month_path = os.path.join(year_path, month)
#                 if not month.startswith('Month='):
#                     print(f"Skipping {month_path}: Not a valid month directory.")
#                     continue

#                 for day_dir in os.listdir(month_path):
#                     day_path = os.path.join(month_path, day_dir)
#                     if not day_dir.startswith('Day='):
#                         print(f"Skipping {day_path}: Not a valid day directory.")
#                         continue

#                     day_files = [os.path.join(day_path, f) for f in os.listdir(day_path) if f.endswith('.parquet')]
#                     if day_files:
#                         print(f"Loading Parquet files from {day_path}.")
#                         with concurrent.futures.ThreadPoolExecutor() as executor:
#                             tables = list(executor.map(read_parquet_file, day_files))
#                         all_tables.extend(tables)

#         if all_tables:
#             combined_table = pa.concat_tables(all_tables)
#             cached_weather_df_nldas = combined_table.to_pandas()
#             cached_weather_df_nldas.dropna(inplace=True)
#             cached_weather_df_nldas.set_index(['s2_token_L9', 's2_token_L7', 's2_token_L5'], inplace=True)
#             print("All data successfully loaded into memory.")
#         else:
#             print("No Parquet tables were loaded.")
#         cache.set('cached_weather_df_nldas', cached_weather_df_nldas)

#     except Exception as e:
#         print(f"Failed to load data: {e}")

# def read_parquet_file(file_path):
#     try:
#         return pa.parquet.read_table(file_path)
#     except Exception as e:
#         print(f"Failed to read Parquet file {file_path}: {e}")
#         return None

# Load a specific file into the cache
# def load_file_into_cache(file_path):
#     try:
#         df = pd.read_parquet(file_path)
#         print(f"Loaded {file_path} with shape: {df.shape}")
#         cached_weather_df_nldas[file_path] = df
#         # Merge or append logic here
#     except Exception as e:
#         print(f"Failed to load {file_path}: {e}")


"""" Logic to monitor data """

WATCHED_DIR = "/mnt/md1/NLDAS/PARQUETE_S2/"

# Store the most recent file data from all subdirectories
recent_file_data = {}

# Custom event handler that reads files when updated
class DirectoryEventHandler(FileSystemEventHandler):
    def on_modified(self, event):
        if event.is_directory:
            return
        else:
            # Read the file that was modified
            file_path = event.src_path
            with open(file_path, 'r') as file:
                recent_file_data[file_path] = file.read()
            print(f"Updated {file_path}")


# Start monitoring the directory for changes
def start_monitoring():
    event_handler = DirectoryEventHandler()
    observer = Observer()
    observer.schedule(event_handler, WATCHED_DIR, recursive=True)
    observer.start()
    print(f"Monitoring started for directory: {WATCHED_DIR}")
    try:
        while True:
            time.sleep(1)
    except KeyboardInterrupt:
        observer.stop()
    observer.join()

cached_weather_df_nldas = {} 

def load_file_into_cache(file_path):
    try:
        df = pd.read_parquet(file_path)
        # print(f"File: {file_path} | Columns: {df.columns} | Shape: {df.shape}")
        cached_weather_df_nldas[file_path] = df  # Store the DataFrame in the cache
        print(f"Data cached for: {file_path}")
    except Exception as e:
        print(f"Failed to load {file_path}: {str(e)}")


# Load all Parquet files into the cache initially
def load_all_data():
    for root, dirs, files in os.walk(WATCHED_DIR):
        for file in files:
            if file.endswith(".parquet"):
                file_path = os.path.join(root, file)
                load_file_into_cache(file_path)
    print(f"All Parquet files loaded into cache")
    
def reload_data_ncep(YYYY_str, MM_str, Day_str):
    global cached_weather_files_ncep
    cached_weather_files_ncep = {}
    try:
        start_time = time.time()
        month_path = os.path.join(base_path_ncep, 's2_tokens_L5', 's2_tokens_L7', 's2_tokens_L9', f'Year={YYYY_str}', f'Month={MM_str}')
        df = pd.read_parquet('/mnt/md1/NLDAS/PARQUETE_S2/s2_tokens_L5/s2_tokens_L7/s2_tokens_L9/Year=2024/Month=06/Day=27/_20240627_016.parquet')
        df_new = pd.read_parquet(file_path) 
        print(f"Reloading data from: {month_path}")
        logging.info(f"Reloading data from: {month_path}")
        
        for day_dir in os.listdir(month_path):
            day_path = os.path.join(month_path, day_dir)
            if os.path.isdir(day_path):
                day_files = [os.path.join(day_path, f) for f in os.listdir(day_path) if f.endswith('.parquet')]
                if day_files:
                    cached_weather_files_ncep[day_dir] = day_files
        
        print(f"Directories reloaded in {time.time() - start_time:.4f} seconds")
        # logging.info(f"Directories reloaded in {time.time() - start_time:.4f} seconds")
        return cached_weather_files_ncep
    except Exception as e:
        print(f"Failed to reload directories: {e}")
        logging.error(f"Failed to reload directories: {e}")

cached_weather_df_ncep = pd.DataFrame() 


# def load_year_data(YYYY_str):
#     global cached_weather_df_ncep
#     cached_weather_df_ncep = pd.DataFrame()  # Initialize as an empty DataFrame
    
#     # Define the base path for the year's data
#     year_path = os.path.join(base_path_ncep, 's2_tokens_L5', 's2_tokens_L7', 's2_tokens_L9', f'Year={YYYY_str}')
#     print(f'Year path: {year_path}')
    
#     try:
#         # Loop through all subdirectories in the year directory that represent months
#         for month_dir in sorted(os.listdir(year_path)):
#             month_path = os.path.join(year_path, month_dir)
            
#             # Ensure the directory exists and is a month directory
#             if os.path.isdir(month_path) and month_dir.startswith("Month="):
#                 day_dfs = []
                
#                 # Loop through all day directories within the month
#                 for day_dir in sorted(os.listdir(month_path)):
#                     day_path = os.path.join(month_path, day_dir)
                    
#                     # Ensure the directory exists and is a day directory
#                     if os.path.isdir(day_path):
#                         day_files = [os.path.join(day_path, f) for f in os.listdir(day_path) if f.endswith('.parquet')]
                        
#                         # Load the data for each day and concatenate
#                         if day_files:
#                             day_dfs.extend([pd.read_parquet(f) for f in day_files])
                
#                 # Concatenate all day data for the month
#                 if day_dfs:
#                     month_df = pd.concat(day_dfs, ignore_index=True).dropna()
                    
#                     # Append the month's data to the year data
#                     cached_weather_df_ncep = pd.concat([cached_weather_df_ncep, month_df], ignore_index=True)
        
#         # Set the index with the specified columns if they exist in the DataFrame
#         if all(col in cached_weather_df_ncep.columns for col in ['s2_token_L9', 's2_token_L7', 's2_token_L5']):
#             cached_weather_df_ncep.set_index(['s2_token_L9', 's2_token_L7', 's2_token_L5'], inplace=True)
#             cached_weather_df_ncep.reset_index(inplace=True)
#         else:
#             print("Not all columns are available in the DataFrame for setting as index.")
            
#     except Exception as e:
#         print(f"Failed to load data for Year {YYYY_str}: {e}")
#         cached_weather_df_ncep = pd.DataFrame()

#     return cached_weather_df_ncep




# def load_day_data(day_files):
#     """Helper function to load and concatenate Parquet files for a single day."""
#     try:
#         tables = [read_parquet_file(f) for f in day_files]
#         dfs = [table.to_pandas() for table in tables]
#         return pd.concat(dfs).dropna()
#     except Exception as e:
#         print(f"Failed to load day data: {e}")
#         return pd.DataFrame()



# def load_all_months_data(YYYY_str, MM_str):
#     global cached_weather_df_ncep
#     cached_weather_df_ncep = pd.DataFrame()  # Initialize as an empty DataFrame
    
#     # Check if data for the given month is already cached
#     if cached_weather_df_ncep is not None and not cached_weather_df_ncep.empty:
#         return cached_weather_df_ncep

#     try:
#         # Define the base path for the month's data
#         month_path = os.path.join(base_path_ncep, 's2_tokens_L5', 's2_tokens_L7', 's2_tokens_L9', f'Year={YYYY_str}', f'Month={MM_str.zfill(2)}')
#         print(f'month-path---{month_path}')
        
#         # Ensure the month directory exists
#         if not os.path.isdir(month_path):
#             raise FileNotFoundError(f"Month directory {month_path} does not exist")

#         # Prepare list of tasks for multiprocessing
#         tasks = []
#         for day_dir in os.listdir(month_path):
#             day_path = os.path.join(month_path, day_dir)
#             if os.path.isdir(day_path):
#                 day_files = [os.path.join(day_path, f) for f in os.listdir(day_path) if f.endswith('.parquet')]
#                 if day_files:
#                     tasks.append(day_files)

#         nproc = max(1, multiprocessing.cpu_count() - 1)  # Use nproc-1 for safety

#         # Use multiprocessing to load data in parallel
#         with multiprocessing.Pool(processes=nproc) as pool:
#             all_dfs = pool.map(load_day_data, tasks)
        
#         # Concatenate all the DataFrames into a single DataFrame
#         if all_dfs:
#             cached_weather_df_ncep = pd.concat(all_dfs).dropna()
            
#             # Set the index with the specified columns
#             # if all(col in cached_weather_df_ncep.columns for col in ['s2_token_L9', 's2_token_L7', 's2_token_L5']):
#             #     cached_weather_df_ncep.set_index(['s2_token_L9', 's2_token_L7', 's2_token_L5'], inplace=True)
#         else:
#             cached_weather_df_ncep = pd.DataFrame()
#     except Exception as e:
#         print(f"Failed to load data for Year {YYYY_str}, Month {MM_str}: {e}")
#         # logging.error(f"Failed to load data for Year {YYYY_str}, Month {MM_str}: {e}")
#         cached_weather_df_ncep = pd.DataFrame()

#     return cached_weather_df_ncep


#NOAA
def reload_data_noaa(YYYY_str, MM_str, Day_str):
    global cached_weather_files_noaa
    cached_weather_files_noaa = {}
    try:
        start_time = time.time()
        month_path = os.path.join(base_path_noaa, 's2_tokens_L5', 's2_tokens_L7', 's2_tokens_L9', 's2_tokens_L20', f'Year={YYYY_str}', f'Month={MM_str}')
        
        print(f"Reloading data from: {month_path}")
        logging.info(f"Reloading data from: {month_path}")
        
        for day_dir in os.listdir(month_path):
            day_path = os.path.join(month_path, day_dir)
            if os.path.isdir(day_path):
                day_files = [os.path.join(day_path, f) for f in os.listdir(day_path) if f.endswith('.parquet')]
                if day_files:
                    cached_weather_files_noaa[day_dir] = day_files
        # print("=====cached_weather_files_noaa=====",cached_weather_files_noaa)
        print(f"Directories reloaded in {time.time() - start_time:.4f} seconds")
        # logging.info(f"Directories reloaded in {time.time() - start_time:.4f} seconds")
        return cached_weather_files_noaa
    except Exception as e:
        print(f"Failed to reload directories: {e}")
        logging.error(f"Failed to reload directories: {e}")

#NOAA FORECASTED
def reload_data_noaa_forecasted(YYYY_str, MM_str, Day_str):
    global cached_weather_files_noaa_forecasted
    cached_weather_files_noaa_forecasted = {}
    try:
        start_time = time.time()
        month_path = os.path.join(base_path_noaa_forecasted, 's2_tokens_L5', 's2_tokens_L7', 's2_tokens_L9', 's2_tokens_L20', f'Year={YYYY_str}', f'Month={MM_str}')
        
        print(f"Reloading data from: {month_path}")
        logging.info(f"Reloading data from: {month_path}")
        
        for day_dir in os.listdir(month_path):
            day_path = os.path.join(month_path, day_dir)
            if os.path.isdir(day_path):
                day_files = [os.path.join(day_path, f) for f in os.listdir(day_path) if f.endswith('.parquet')]
                if day_files:
                    cached_weather_files_noaa_forecasted[day_dir] = day_files
        # print("=====cached_weather_files_noaa_forecasted=====",cached_weather_files_noaa_forecasted)
        print(f"Directories reloaded in {time.time() - start_time:.4f} seconds")
        # logging.info(f"Directories reloaded in {time.time() - start_time:.4f} seconds")
        return cached_weather_files_noaa_forecasted
    except Exception as e:
        print(f"Failed to reload directories: {e}")
        logging.error(f"Failed to reload directories: {e}")

#GHCND
def reload_data_ghcnd(YYYY_str, MM_str, Day_str):
    global cached_weather_files_ghcnd
    cached_weather_files_ghcnd = {}
    try:
        start_time = time.time()
        month_path = os.path.join(base_path_ghcnd, 's2_tokens_L5', 's2_tokens_L7', 's2_tokens_L9', f'Year={YYYY_str}', f'Month={MM_str}')
        
        print(f"Reloading data from: {month_path}")
        logging.info(f"Reloading data from: {month_path}")
        
        for day_dir in os.listdir(month_path):
            day_path = os.path.join(month_path, day_dir)
            if os.path.isdir(day_path):
                day_files = [os.path.join(day_path, f) for f in os.listdir(day_path) if f.endswith('.parquet')]
                if day_files:
                    cached_weather_files_ghcnd[day_dir] = day_files
        # print("=====cached_weather_files_ghcnd=====",cached_weather_files_ghcnd)
        print(f"Directories reloaded in {time.time() - start_time:.4f} seconds")
        # logging.info(f"Directories reloaded in {time.time() - start_time:.4f} seconds")
        return cached_weather_files_ghcnd
    except Exception as e:
        print(f"Failed to reload directories: {e}")
        logging.error(f"Failed to reload directories: {e}")

#AUS FORECASTED
def reload_data_aus_forecasted(YYYY_str, MM_str, Day_str):
    global cached_weather_files_aus_forecasted
    cached_weather_files_aus_forecasted = {}
    try:
        start_time = time.time()
        month_path = os.path.join(base_path_aus_forecasted, 's2_tokens_L5', 's2_tokens_L7', 's2_tokens_L9', f'Year={YYYY_str}', f'Month={MM_str}')
        
        print(f"Reloading data from: {month_path}")
        logging.info(f"Reloading data from: {month_path}")
        
        for day_dir in os.listdir(month_path):
            day_path = os.path.join(month_path, day_dir)
            if os.path.isdir(day_path):
                day_files = [os.path.join(day_path, f) for f in os.listdir(day_path) if f.endswith('.parquet')]
                if day_files:
                    cached_weather_files_aus_forecasted[day_dir] = day_files
        print("=====cached_weather_files_aus_forecasted=====",cached_weather_files_aus_forecasted)
        print(f"Directories reloaded in {time.time() - start_time:.4f} seconds")
        # logging.info(f"Directories reloaded in {time.time() - start_time:.4f} seconds")
        return cached_weather_files_aus_forecasted
    except Exception as e:
        print(f"Failed to reload directories: {e}")
        logging.error(f"Failed to reload directories: {e}")

#AUS
def reload_data_aus(YYYY_str, MM_str, Day_str):
    global cached_weather_files_aus
    cached_weather_files_aus = {}
    try:
        start_time = time.time()
        month_path = os.path.join(base_path_aus, 's2_tokens_L5', 's2_tokens_L7', 's2_tokens_L9', f'Year={YYYY_str}', f'Month={MM_str}')
        
        print(f"Reloading data from: {month_path}")
        logging.info(f"Reloading data from: {month_path}")
        
        for day_dir in os.listdir(month_path):
            day_path = os.path.join(month_path, day_dir)
            if os.path.isdir(day_path):
                day_files = [os.path.join(day_path, f) for f in os.listdir(day_path) if f.endswith('.parquet')]
                if day_files:
                    cached_weather_files_aus[day_dir] = day_files
        print("=====cached_weather_files_aus=====",cached_weather_files_aus)
        print(f"Directories reloaded in {time.time() - start_time:.4f} seconds")
        # logging.info(f"Directories reloaded in {time.time() - start_time:.4f} seconds")
        return cached_weather_files_aus
    except Exception as e:
        print(f"Failed to reload directories: {e}")
        logging.error(f"Failed to reload directories: {e}")

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
        "time": "ISO 8601"
    }
    
    return metadata



def meta_data_with_description():
    metadata = {
        "CRAIN_surface": "Categorical Rain",
        "DLWRF_surface": "Downward Long-Wave Rad. Flux",
        "DPT_2maboveground": "Dew Point Temperature",
        "DSWRF_surface": "Downward Short-Wave Radiation Flux",
        "Day": "Day of the month ",
        "GUST_surface": "Wind Speed (Gust)",
        "Month": "Month of the year [integer",
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
        "Year": "Year of the record [integer]",
        "index": "Index of the record [integer]",
        "latitude": "Latitude [degrees]",
        "longitude": "Longitude [degrees]",
        "time": "Time of the record in ISO 8601 format [YYYY-MM-DDTHH:MM:SSZ]"
    }
    
    return metadata




def get_GHCND_metdata():
    metadata = {
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
        "TAVG": "Celsius",
        "TMAX": "Celsius",
        "TMIN": "Celsius",
        "WDF2": "degrees",
        "WDF5": "degrees",
        "WSF2": "km/h",
        "WSF5": "km/h",
        "Year": "integer",
        "aifstime_utc": "ISO 8601 format",
        "elev": "meters",
        "index": "integer",
        "lat": "degrees",
        "lon": "degrees"
    }
    return metadata


def get_noaa_metadata():
    metadata = {
        "Day": "integer",
        "ETo_FAO_inches": "inches",
        "ETo_HAR_inches": "inches",
        "ETo_average_inches": "inches",
        "Month": "integer",
        "P": "pascals",
        "RH_max": "percentage",
        "RH_mean": "percentage",
        "RH_min": "percentage",
        "R_n__MJpm2": "MJ / m²",
        "R_s__MJpm2": "MJ / m²",
        "TS": "ISO 8601 timestamp",
        "T_dew": "Celsius",
        "T_max": "Celsius",
        "T_mean": "Celsius",
        "T_min": "Celsius",
        "U_z": "m/s",
        "Year": "integer",
        "index": "integer",
        "lat": "degrees",
        "lon": "degrees",
        "z_msl": "meters"
    }

    return metadata



def get_noaa_description():
    metadata_noaa_description = {
        "Day": "Day of the month [integer]",
        "ETo_FAO_inches": "Evapotranspiration using FAO method [inches]",
        "ETo_HAR_inches": "Evapotranspiration using Hargreaves method [inches]",
        "ETo_average_inches": "Average evapotranspiration [inches]",
        "Month": "Month of the year [integer]",
        "P": "Pressure [pascals]",
        "RH_max": "Maximum relative humidity [percentage]",
        "RH_mean": "Average relative humidity [percentage]",
        "RH_min": "Minimum relative humidity [percentage]",
        "R_n__MJpm2": "Net radiation [MJ / m²]",
        "R_s__MJpm2": "Solar radiation [MJ / m²]",
        "TS": "Timestamp of the observation [ISO 8601 format]",
        "T_dew": "Dew point temperature [Celsius]",
        "T_max": "Maximum daily air temperature [Celsius]",
        "T_mean": "Mean daily air temperature [Celsius]",
        "T_min": "Minimum daily air temperature [Celsius]",
        "U_z": "Wind speed [m/s]",
        "Year": "Year of observation [integer]",
        "index": "Index or row number [integer]",
        "lat": "Latitude [degrees]",
        "lon": "Longitude [degrees]",
        "z_msl": "Elevation above mean sea level [meters]"
    }

    return metadata_noaa_description


def get_GHCND_description():
    metadata_global_description = {
        "AWND": "Average daily wind speed [meters per second]",
        "Day": "Day of the month [integer]",
        "ETo_AVG_IN": "Average evapotranspiration [inches]",
        "ETo_AVG_MM": "Average evapotranspiration [millimeters]",
        "ETo_FAO_IN": "Evapotranspiration using FAO method [inches]",
        "ETo_FAO_MM": "Evapotranspiration using FAO method [millimeters]",
        "ETo_HAR_IN": "Evapotranspiration using Hargreaves method [inches]",
        "ETo_HAR_MM": "Evapotranspiration using Hargreaves method [millimeters]",
        "Month": "Month of observation [integer]",
        "PRCP": "Precipitation [inches]",
        "TAVG": "Average temperature [Celsius]",
        "TMAX": "Maximum temperature [Celsius]",
        "TMIN": "Minimum temperature [Celsius]",
        "WDF2": "Direction of fastest 2-minute wind [degrees]",
        "WDF5": "Direction of fastest 5-minute wind [degrees]",
        "WSF2": "Fastest 2-minute wind speed [km/h]",
        "WSF5": "Fastest 5-minute wind speed [km/h]",
        "Year": "Year of observation [integer]",
        "aifstime_utc": "UTC timestamp [ISO 8601 format]",
        "elev": "Elevation [meters]",
        "index": "Index or row number [integer]",
        "lat": "Latitude [degrees]",
        "lon": "Longitude [degrees]"
    }

    return metadata_global_description

def get_aus_forcasted_metadata():
    metadata = {
        "Day": "integer",
        "ETo_FAO_inches": "inches",
        "ETo_HAR_inches": "inches",
        "ETo_average_inches": "inches",
        "Month": "integer",
        "R_n__MJpm2": "MJ/m^2",
        "R_s__MJpm2": "MJ/m^2",
        "STN": "string",
        "TS": "ISO 8601",
        "T_max": "Celsius",
        "T_mean": "Celsius",
        "T_min": "Celsius",
        "Year": "integer",
        "index": "integer",
        "lat": "degrees",
        "lon": "degrees",
        "sub_region": "string",
        "z_msl": "meters"
    }
    return metadata

def get_aus_forcasted_description():
    metadata_aus_forecasted_description = {
        "Day": "Day of the month [integer]",
        "ETo_FAO_inches": "Evapotranspiration using FAO method [inches]",
        "ETo_HAR_inches": "Evapotranspiration using Hargreaves method [inches]",
        "ETo_average_inches": "Average evapotranspiration [inches]",
        "Month": "Month of the year [integer]",
        "R_n__MJpm2": "Net radiation at crop surface",
        "R_s__MJpm2": "Solar radiation at crop surface",
        "STN": "Station identifier [string]",
        "TS": "Timestamp in ISO 8601 format",
        "T_max": "Maximum temperature of the day [Celsius]",
        "T_mean": "Mean temperature of the day [Celsius]",
        "T_min": "Minimum temperature of the day [Celsius]",
        "Year": "Year of the observation [integer]",
        "index": "Unique identifier for the record [integer]",
        "lat": "Latitude of the station [degrees]",
        "lon": "Longitude of the station [degrees]",
        "sub_region": "Geographic sub-region [string]",
        "z_msl": "Elevation above mean sea level [meters]"
    }
    return metadata_aus_forecasted_description

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
        "rain_trace": "mm",
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
        "Day": "Day of the month [integer]",
        "ETo_AVG_IN": "Average evapotranspiration [inches]",
        "ETo_FAO_IN": "Evapotranspiration using FAO method [inches]",
        "ETo_HAR_IN": "Evapotranspiration using Hargreaves method [inches]",
        "Month": "Month of the year [integer]",
        "RH_mean": "Mean relative humidity [%]",
        "T_dew": "Dew point temperature [Celsius]",
        "T_max": "Maximum temperature of the day [Celsius]",
        "T_mean": "Mean temperature of the day [Celsius]",
        "T_min": "Minimum temperature of the day [Celsius]",
        "U_z": "Wind speed at reference height [m/s]",
        "Year": "Year of the observation [integer]",
        "air_temp": "Air temperature [Celsius]",
        "apparent_t": "Apparent temperature [Celsius]",
        "cloud": "Cloud cover description [string]",
        "cloud_base_m": "Cloud base height [meters]",
        "cloud_oktas": "Cloud cover in oktas [oktas]",
        "cloud_type": "Type of clouds [string]",
        "cloud_type_id": "Cloud type identifier [integer]",
        "delta_t": "Temperature difference between dry and wet bulbs [Celsius]",
        "dewpt": "Dew point temperature [Celsius]",
        "gust_kmh": "Wind gust speed [km/h]",
        "gust_kt": "Wind gust speed [knots]",
        "history_product": "Historical weather product [string]",
        "index": "Unique identifier for the record [integer]",
        "lat": "Latitude of the location [degrees]",
        "local_date_time": "Local date and time [string]",
        "local_date_time_full": "Full local date and time in format YYYYMMDDHHMMSS",
        "lon": "Longitude of the location [degrees]",
        "name": "Station or location name [string]",
        "press": "Air pressure [hPa]",
        "press_msl": "Mean sea level air pressure [hPa]",
        "press_qnh": "Air pressure at the station [hPa]",
        "press_tend": "Pressure tendency [R for Rising, F for Falling]",
        "rain_trace": "Precipitation in the last 24 hours [mm]",
        "rel_hum": "Relative humidity [%]",
        "sea_state": "State of the sea [string]",
        "sort_order": "Sort order index [float]",
        "swell_dir_worded": "Swell direction [string]",
        "swell_height": "Swell height [float]",
        "swell_period": "Swell period [string]",
        "vis_km": "Visibility distance [km]",
        "weather": "Weather description [string]",
        "wind_dir": "Wind direction [string]",
        "wind_spd_kmh": "Wind speed [km/h]",
        "wind_spd_kt": "Wind speed [knots]",
        "wmo": "World Meteorological Organization station identifier [float]"
    }
    return metadata_aus_description

def clear_cache(cache):
    """Clear the specified cache."""
    cache.clear()


def memoize_ncep(func):
    cache = {}
    last_keys = {}

    @wraps(func)
    def wrapper(*args, **kwargs):
        clear_cache(cache)

        # Create a unique cache key based on all arguments
        key = ('NCEP', args, frozenset(kwargs.items()))
        if key not in cache:
            print(f"NCEP Cache miss for key: {key}")
            cache[key] = func(*args, **kwargs)
            last_keys[key] = time.time()
        else:
            print(f"NCEP Cache hit for key: {key}")
        return cache[key]

    return wrapper

def memoize_nldas(func):
    cache = {}
    last_keys = {}

    @wraps(func)
    def wrapper(*args, **kwargs):
        clear_cache(cache)

        # Create a unique cache key based on all arguments
        key = ('NLDAS', args, frozenset(kwargs.items()))
        if key not in cache:
            print(f"NLDAS Cache miss for key: {key}")
            cache[key] = func(*args, **kwargs)
            last_keys[key] = time.time()
        else:
            print(f"NLDAS Cache hit for key: {key}")
        return cache[key]

    return wrapper

#NOAA
def memoize_noaa(func):
    cache = {}
    last_keys = {}

    @wraps(func)
    def wrapper(*args, **kwargs):
        clear_cache(cache)

        # Create a unique cache key based on all arguments
        key = ('NOAA', args, frozenset(kwargs.items()))
        if key not in cache:
            print(f"NOAA Cache miss for key: {key}")
            cache[key] = func(*args, **kwargs)
            last_keys[key] = time.time()
        else:
            print(f"NOAA Cache hit for key: {key}")
        return cache[key]

    return wrapper

#NOAA FORECASTED
def memoize_noaa_forecasted(func):
    cache = {}
    last_keys = {}

    @wraps(func)
    def wrapper(*args, **kwargs):
        clear_cache(cache)

        # Create a unique cache key based on all arguments
        key = ('NOAA_FORECASTED', args, frozenset(kwargs.items()))
        if key not in cache:
            print(f"NOAA_FORECASTED Cache miss for key: {key}")
            cache[key] = func(*args, **kwargs)
            last_keys[key] = time.time()
        else:
            print(f"NOAA_FORECASTED Cache hit for key: {key}")
        return cache[key]

    return wrapper

#GHCND
def memoize_ghcnd(func):
    cache = {}
    last_keys = {}

    @wraps(func)
    def wrapper(*args, **kwargs):
        clear_cache(cache)

        # Create a unique cache key based on all arguments
        key = ('GHCND', args, frozenset(kwargs.items()))
        if key not in cache:
            print(f"GHCND Cache miss for key: {key}")
            cache[key] = func(*args, **kwargs)
            last_keys[key] = time.time()
        else:
            print(f"GHCND Cache hit for key: {key}")
        return cache[key]

    return wrapper

#AUS_FORECASTED
def memoize_aus_forecasted(func):
    cache = {}
    last_keys = {}

    @wraps(func)
    def wrapper(*args, **kwargs):
        clear_cache(cache)

        # Create a unique cache key based on all arguments
        key = ('AUS_FORECASTED', args, frozenset(kwargs.items()))
        if key not in cache:
            print(f"AUS FORECASTED Cache miss for key: {key}")
            cache[key] = func(*args, **kwargs)
            last_keys[key] = time.time()
        else:
            print(f"AUS FORECASTED Cache hit for key: {key}")
        return cache[key]

    return wrapper

#AUS
def memoize_aus(func):
    cache = {}
    last_keys = {}

    @wraps(func)
    def wrapper(*args, **kwargs):
        clear_cache(cache)

        # Create a unique cache key based on all arguments
        key = ('AUS', args, frozenset(kwargs.items()))
        if key not in cache:
            print(f"AUS Cache miss for key: {key}")
            cache[key] = func(*args, **kwargs)
            last_keys[key] = time.time()
        else:
            print(f"AUS Cache hit for key: {key}")
        return cache[key]

    return wrapper

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

def getWeatherFromNCEP(agstack_geoid, dtStr, end_date=None):
    filePath = '/mnt/md2/NCEP/PARQUET_S2/'
    tok = dtStr.split('-')
    YYYY_str = int(tok[0])
    MM_str = int(tok[1])
    DD_str = int(tok[2])

    # Create a datetime object for start date and localize it to UTC
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

    s1_time_start = time.time()
    w_ret = pd.DataFrame()
    try:
        # Fetch WKT polygon and extract latitude and longitude
        wkt_polygon = fetchWKT(agstack_geoid)
        lat, lon = extractLatLonFromWKT(wkt_polygon)
        start_time = time.time()
        
        # Get the list of S2 indices and CIDs for the data point
        s2_index__L5_list, L5_cids = get_s2_cellids_and_token_list(5, [lat], [lon])

    
        list_of_5_paths = [filePath + 's2_token_L5=' + x for x in s2_index__L5_list
                           if os.path.exists(filePath + 's2_token_L5=' + x)]
        if list_of_5_paths == []:
            return empty_response()

        weather_datasets = []
        for x in list_of_5_paths:
            weather_datasets.append(ds.dataset(x, format="parquet", partitioning="hive"))

        w_all = pd.DataFrame()
        for weatherDataset in weather_datasets:
            if end_date:
                # Filter between start and end date
                weather_df = weatherDataset.to_table(
                    filter=(
                        (ds.field('Year') >= YYYY_str) & (ds.field('Year') <= end_YYYY_str) &
                        (ds.field('Month') >= MM_str) & (ds.field('Month') <= end_MM_str) &
                        (ds.field('Day') >= DD_str) & (ds.field('Day') <= end_DD_str)
                    )
                ).to_pandas()
                w_all = pd.concat([w_all, weather_df], ignore_index=True)
            elif end_date is None:
                # Filter for start date only
                weather_df = weatherDataset.to_table(
                    filter=(
                        (ds.field('Year') == YYYY_str) &
                        (ds.field('Month') == MM_str) &
                        (ds.field('Day') == DD_str)
                    )
                ).to_pandas()

                w_all = pd.concat([w_all, weather_df], ignore_index=True)
            if not weather_df.empty:
                weather_df.dropna(inplace=True)

            # print(weather_df['time'].head())
            # print(weather_df['time'].tail())
            # print(weather_df.dtypes)
            # Concatenate to the main DataFrame

        if len(w_all) > 0:
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

    except Exception as e:
        print(e)
        w_ret = empty_response()

    end_time = time.time()
    time_elapsed = (end_time - start_time)

    return w_ret






def getWeatherFromNLDAS(agstack_geoid, dtStr):
    filePath = '/mnt/md2/NLDAS/PARQUET_S2/'
    tok = dtStr.split('-')
    YYYY_str = tok[0]
    MM_str = tok[1]
    DD_str = tok[2]
    


    # Create a datetime object and localize it to UTC
    local_dt = datetime(int(YYYY_str), int(MM_str), int(DD_str))
    utc_dt = pytz.utc.localize(local_dt)  # Convert to UTC

    # Generate the next 7 days
    date_range = [utc_dt + timedelta(days=i) for i in range(8)]
    
    s1_time_start = time.time()
    w_ret = pd.DataFrame()
    
        # Fetch WKT polygon and extract latitude and longitude
    try:

        
        wkt_polygon = fetchWKT(agstack_geoid)
        lat, lon = extractLatLonFromWKT(wkt_polygon)
        start_time = time.time()
        # Get the list of S2 indices and CIDs for the data point
        s2_index__L5_list, L5_cids = get_s2_cellids_and_token_list(5, [lat], [lon])
        list_of_5_paths = [filePath + 's2_token_L5=' + x for x in s2_index__L5_list
                        if os.path.exists(filePath + 's2_token_L5=' + x)]
    
        #     return empty_response()
        weather_datasets = []
        for x in list_of_5_paths:
            weather_datasets.append(ds.dataset(x, format="parquet", partitioning="hive"))

        w_all = pd.DataFrame()
        for weatherDataset in weather_datasets:
            
            for current_date in date_range:
                
                # Use UTC date components for filtering
                YYYY_list = [current_date.year]
                MM_list = [str(current_date.month).zfill(2)]
                DD_list = [str(current_date.day).zfill(2)]
            
            

                weather_df = weatherDataset.to_table(
                    filter=(
                        ds.field('Year').isin(YYYY_list) &
                        ds.field('Month').isin(MM_list) &
                        ds.field('Day').isin(DD_list)
                    )
                ).to_pandas()
                # print(weather_df)
                # if not weather_df.empty:
                #     weather_df.dropna(inplace=True)
                w_all = pd.concat([w_all, weather_df], ignore_index=True)
        
        if len(w_all) > 0:
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
        
    except Exception as e:
        print(f'Error -{e}')
        w_ret = empty_response()
    end_time = time.time()
    time_elapsed = (end_time - start_time)

    return w_ret
    

# @memoize_ncep
# def getWeatherFromNCEP(dtStr, agstack_geoid):
#     global cached_weather_df_ncep
#     start_time_total = time.time()

#     YYYY_str, MM_str, DD_str = dtStr.split('-')

#     # Load all available months' data for the specified year if not already loaded
#     if cached_weather_df_ncep.empty:
#         cached_weather_df_ncep = load_year_data(YYYY_str)

#     # Fetch WKT polygon and extract latitude and longitude
#     wkt_polygon = fetchWKT(agstack_geoid)
#     lat, lon = extractLatLonFromWKT(wkt_polygon)
#     s2_index__L9_list, _ = get_s2_cellids_and_token_list(9, [lat], [lon])
#     s2_index__L9_list = [str(token) for token in s2_index__L9_list]
#     weather_df = pd.DataFrame()

#     # Filtering by level 9, 7, and 5
#     try:
#         if not cached_weather_df_ncep.empty:
#             print('At level 9')
#             weather_df = cached_weather_df_ncep.loc[(s2_index__L9_list, slice(None), slice(None)), :]
#     except KeyError as e:
#         print(f'KeyError at level 9: {e}')
#         weather_df = pd.DataFrame()

#     # Filter data by level 7 tokens if level 9 data is empty
#     if weather_df.empty:
#         level_7_tokens = get_level_7_tokens(s2_index__L9_list)
#         try:
#             print('At level 7')
#             weather_df = cached_weather_df_ncep.loc[(slice(None), level_7_tokens, slice(None)), :]
#         except KeyError as e:
#             print(f'KeyError at level 7: {e}')
#             weather_df = pd.DataFrame()

#     # Filter data by level 5 tokens if level 7 data is empty
#     if weather_df.empty:
#         level_5_tokens = get_level_5_tokens(s2_index__L9_list)
#         try:
#             print('At level 5')
#             weather_df = cached_weather_df_ncep.loc[(slice(None), slice(None), level_5_tokens), :]
#         except KeyError as e:
#             print(f'KeyError at level 5: {e}')
#             weather_df = pd.DataFrame()

#     if not weather_df.empty and len(weather_df) > 1:
#         weather_df['time'] = pd.to_datetime(weather_df['time'], errors='coerce')
#         numeric_cols = weather_df.columns.drop('time')
#         weather_df[numeric_cols] = weather_df[numeric_cols].apply(pd.to_numeric, errors='coerce')
#         weather_df = weather_df.groupby('time').mean().reset_index()
#         weather_df['time'] = weather_df['time'].dt.strftime('%Y-%m-%dT%H:%M:%S')
#         weather_df = weather_df[
#             (weather_df['Year'].astype(int) == int(YYYY_str)) & 
#             (weather_df['Month'].astype(int) == int(MM_str)) & 
#             (weather_df['Day'].astype(int) == int(DD_str))
#         ]

#     if weather_df.empty:
#         return pd.DataFrame(empty_response())
    
#     end_time_total = time.time()
#     total_duration = end_time_total - start_time_total
#     print(f'Total time taken to complete the process: {total_duration}')
#     print("Reached end of function")

#     return weather_df


# @memoize_ncep
# def getWeatherFromNCEP(agstack_geoid,start_date,end_date=None):
    
#     # logging.info(f"Data Processing started for geoid: {agstack_geoid}, date: {dtStr}")

#     print(f'start_date---------{start_date}')

#     # Parse the date
#     # tok = start_date.split('-')
#     # YYYY_str = tok[0]
#     # MM_str = tok[1].zfill(2)
#     # DD_str = tok[2]
#     start_date = pd.to_datetime(start_date)
#     if end_date is not None:
#         end_date = pd.to_datetime(end_date)

#     YYYY_str = start_date.strftime('%Y')
#     MM_str = start_date.strftime('%m')
#     DD_str = start_date.strftime('%d')

#     global cached_weather_files_ncep, last_reload_time_ncep, cached_weather_df_ncep
#     start_time_total = time.time()

#     day_key = f'Day={DD_str}'

#     # Reload data if necessary
#     start_time_reload = time.time()
#     day_data = reload_data_ncep(YYYY_str, MM_str, DD_str)
#     reload_duration = time.time() - start_time_reload
#     # logging.info(f"Data loading into memory took {reload_duration:.2f} seconds")

    
#     # Read and combine Parquet files using threading
#     if end_date is not None:
#         date_range = pd.date_range(start=start_date, end=end_date, freq='D')
#         all_dfs = []
        
#         for current_date in date_range:
#             YYYY_str = current_date.strftime('%Y')
#             MM_str = current_date.strftime('%m')
#             DD_str = current_date.strftime('%d')
            
#             day_key = f'Day={DD_str}'

#             # Reload data for the current date
#             start_time_reload = time.time()
#             day_data = reload_data_ncep(YYYY_str, MM_str, DD_str)
#             reload_duration = time.time() - start_time_reload
            
#             try:
#                 with concurrent.futures.ThreadPoolExecutor() as executor:
#                     tables = list(executor.map(read_parquet_file, day_data[day_key]))
#                 dfs = [table.to_pandas() for table in tables]
#                 daily_df = pd.concat(dfs)
#                 daily_df.dropna(inplace=True)
#                 daily_df.set_index(['s2_token_L9', 's2_token_L7', 's2_token_L5'], inplace=True)
#                 all_dfs.append(daily_df)
#             except Exception as e:
#                 # logging.error(f"Error reading Parquet files for date {current_date}: {e}")
#                 continue  # Skip to the next day

#         # Concatenate all dataframes for all days in the range
#         if all_dfs:
#             cached_weather_df_ncep = pd.concat(all_dfs)
#         else:
#             cached_weather_df_ncep = pd.DataFrame()  # Empty dataframe if no data was fetched
#     else:
#         # Original logic for fetching only start_date data
#         day_key = f'Day={DD_str}'

#         start_time_reload = time.time()
#         day_data = reload_data_ncep(YYYY_str, MM_str, DD_str)
#         reload_duration = time.time() - start_time_reload

#         try:
#             with concurrent.futures.ThreadPoolExecutor() as executor:
#                 tables = list(executor.map(read_parquet_file, day_data[day_key]))
#             dfs = [table.to_pandas() for table in tables]
#             cached_weather_df_ncep = pd.concat(dfs)
#             cached_weather_df_ncep.dropna(inplace=True)
#             cached_weather_df_ncep.set_index(['s2_token_L9', 's2_token_L7', 's2_token_L5'], inplace=True)
#         except Exception as e:
#             # logging.error(f"Error reading Parquet files: {e}")
#             return empty_response()

#    # Fetch WKT polygon and extract lat/lon
#     # start_time_polygon = time.time()
#     wkt_polygon = fetchWKT(agstack_geoid)
#     lat, lon = extractLatLonFromWKT(wkt_polygon)
#     # polygon_duration = time.time() - start_time_polygon
#     # logging.info(f"Fetching and processing WKT polygon took {polygon_duration:.2f} seconds")

#     # Filter data by level 9 tokens
#     # start_time_filtering_level9 = time.time()
#     s2_index__L9_list, _ = get_s2_cellids_and_token_list(9, [lat], [lon])
#     s2_index__L9_list = [str(token) for token in s2_index__L9_list]
#     # print(s2_index__L9_list)
#     try:
#         weather_df = cached_weather_df_ncep.loc[(s2_index__L9_list, slice(None), slice(None))]
#         # filtering_level9_duration = time.time() - start_time_filtering_level9
#         # logging.info(f"Filtering data for level 9 tokens took {filtering_level9_duration:.2f} seconds")
#     except KeyError:
#         weather_df = pd.DataFrame()

#     # Filter data by level 7 tokens if level 9 data is empty
#     if weather_df.empty:
#         # start_time_filtering_level7 = time.time()
#         level_7_tokens = get_level_7_tokens(s2_index__L9_list)
#         try:
#             weather_df = cached_weather_df_ncep.loc[(slice(None), level_7_tokens, slice(None))]
#             # filtering_level7_duration = time.time() - start_time_filtering_level7
#             # logging.info(f"Filtering data for level 7 tokens took {filtering_level7_duration:.2f} seconds")
#         except KeyError as e:
#             print(f'{e} at level 7')
#             weather_df = pd.DataFrame()

    
#     print(weather_df.tail())
#     # Filter data by level 5 tokens if level 7 data is empty
#     if weather_df.empty:
#         # start_time_filtering_level5 = time.time()
#         level_5_tokens = get_level_5_tokens(s2_index__L9_list)
#         try:
#             weather_df = cached_weather_df_ncep.loc[(slice(None), slice(None), level_5_tokens)]
#             # filtering_level5_duration = time.time() - start_time_filtering_level5
#             # logging.info(f"Filtering data for level 5 tokens took {filtering_level5_duration:.2f} seconds")
#         except KeyError as e:
#             print(f'{e} at level 5')
#             weather_df = pd.DataFrame()

#     if not end_date == None:
#         print(f'weather_df["time"]--{weather_df["time"]}')
#         print(f'start_date--{start_date}')
#         print(f'end_date--{end_date}')
#             # Filter data by the date range
#         weather_df = weather_df[
#             (pd.to_datetime(weather_df['time']) >= start_date) &
#             (pd.to_datetime(weather_df['time']) <= end_date)
#         ]

        
#     elif end_date is None:

#         weather_df = weather_df[(weather_df['Year'].astype(int) == int(YYYY_str)) & 
#                         (weather_df['Month'].astype(int) == int(MM_str)) & 
#                         (weather_df['Day'].astype(int) == int(DD_str))]

#     # print(f'weather_df--{weather_df}')
#     if len(weather_df) > 1:
#         # Convert the time column to datetime first, if it isn't already
#         weather_df['time'] = pd.to_datetime(weather_df['time'], errors='coerce')

#         # Convert only the relevant columns to numeric, excluding 'time'
#         numeric_cols = weather_df.columns.drop('time')
#         weather_df[numeric_cols] = weather_df[numeric_cols].apply(pd.to_numeric, errors='coerce')

#         # Group by the 'time' column and take the mean of numeric columns
#         weather_df = weather_df.groupby('time').mean().reset_index()

#         # Ensure the 'time' column is a datetime object after grouping
#         weather_df['time'] = pd.to_datetime(weather_df['time'], errors='coerce')

#         # Convert the 'time' column to ISO format using strftime
#         weather_df['time'] = weather_df['time'].dt.strftime('%Y-%m-%dT%H:%M:%S')

#     # Return empty response if no data is found
#     if weather_df.empty:
#         weather_df = pd.DataFrame(empty_response())

#     # print(f"NCEP Data for geoid {agstack_geoid} on {dtStr}: {weather_df}")
#     end_time_total = time.time()
#     total_duration = end_time_total - start_time_total
#     # logging.info(f"Total function execution took {total_duration:.2f} seconds")
#     print(f'Total function execution took {total_duration:.2f} seconds')

#     return weather_df


@memoize_ncep
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

#NOAA
def getWeatherFromNOAA(agstack_geoid, start_date,end_date):
    # filePath = '/mnt/md1/NLDAS/PARQUETE_S2/'
    # filePath = '/home/rajat/Downloads/Rnaura_Work/mnt/md1/NOAA/DAILY/PROCESSED/PARQUET_S2/'
    filePath = '/mnt/md2/NOAA/DAILY/PROCESSED/PARQUET_S2/'
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
    try:


        w_ret = pd.DataFrame()
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

"""
@memoize_noaa
def getWeatherFromNOAA(dtStr, agstack_geoid):
    # Parse the date
    tok = dtStr.split('-')
    YYYY_str = tok[0]
    MM_str = tok[1].zfill(2)
    DD_str = tok[2]

    start_time_total = time.time()
    day_key = f'Day={DD_str}'
    cached_weather_df_noaa = pd.DataFrame()

    start_time_reload = time.time()
    day_data = reload_data_noaa(YYYY_str, MM_str, DD_str)
    reload_duration = time.time() - start_time_reload

    try:
        # Load and combine data from Parquet files
        with concurrent.futures.ThreadPoolExecutor() as executor:
            tables = list(executor.map(read_parquet_file, day_data[day_key]))

        combined_table = pa.concat_tables(tables)
        cached_weather_df_noaa = combined_table.to_pandas()

        # Remove NaN values and set the multi-level index
        cached_weather_df_noaa.dropna(inplace=True)
        cached_weather_df_noaa.set_index(['s2_token_L9', 's2_token_L7', 's2_token_L5'], inplace=True)

    except Exception as e:
        logging.error(f"Error reading Parquet files: {e}")
        return empty_response()

    # Extract lat/lon from the geoid
    wkt_polygon = fetchWKT(agstack_geoid)
    lat, lon = extractLatLonFromWKT(wkt_polygon)

    # Obtain level 9 tokens
    s2_index__L9_list, _ = get_s2_cellids_and_token_list(9, [lat], [lon])
    s2_index__L9_list = [str(token) for token in s2_index__L9_list]

    # Filter by level 9 tokens
    try:
        weather_df = cached_weather_df_noaa.loc[
            (slice(None), s2_index__L9_list, slice(None), slice(None))
        ]

    except KeyError:
        weather_df = pd.DataFrame()

    # If level 7 data is empty, filter by level 7 tokens
    if weather_df.empty:
        level_7_tokens = get_level_7_tokens(s2_index__L9_list)
        try:
            weather_df = cached_weather_df_noaa.loc[
                (slice(None), level_7_tokens, slice(None))
            ]
        except KeyError as e:
            print(f'{e} at level 7')
            weather_df = pd.DataFrame()
            

    # If level 5 data is empty, filter by level 5 tokens
    if weather_df.empty:
        level_5_tokens = get_level_5_tokens(s2_index__L9_list)
        try:
            weather_df = cached_weather_df_noaa.loc[
                (slice(None), slice(None), level_5_tokens)
            ]
        except KeyError as e:
            print(f'{e} at level 5')
            weather_df = pd.DataFrame()

    # Filter the DataFrame by date
    try:
        weather_df = weather_df[(weather_df['Year'].astype(int) == int(YYYY_str)) & 
                    (weather_df['Month'].astype(int) == int(MM_str)) & 
                    (weather_df['Day'].astype(int) == int(DD_str))]
    except Exception as e:
        print(e)

    # Handle time series data if more than one entry exists
    if len(weather_df) > 1:
        # Convert Timestamp column to datetime and calculate the mean of numeric columns
        weather_df['TS'] = pd.to_datetime(weather_df['TS'], errors='coerce')
        numeric_cols = weather_df.select_dtypes(include=['number']).columns
        non_numeric_cols = weather_df.select_dtypes(exclude=['number']).columns

        # Group by the TS column (now as part of the index), and aggregate as needed
        weather_df = weather_df.groupby('TS').agg({**{col: 'mean' for col in numeric_cols}, **{col: 'first' for col in non_numeric_cols if col != 'TS'}}).reset_index()

        # If the TS column is part of the index, it will automatically be added as a column after reset_index()
        weather_df['Year'] = weather_df['TS'].dt.year
        weather_df['Month'] = weather_df['TS'].dt.month
        weather_df['Day'] = weather_df['TS'].dt.day

        # Convert Year, Month, Day columns to integers
        weather_df['Year'] = weather_df['Year'].astype(int)
        weather_df['Month'] = weather_df['Month'].astype(int)
        weather_df['Day'] = weather_df['Day'].astype(int)

    elif weather_df.empty:
        weather_df = empty_response()
    try:

        weather_df = weather_df.drop('Timestamp',axis=1)
        weather_df = weather_df.drop('STN',axis=1)
    except:
        pass

    print(f"NCEP Data for geoid {agstack_geoid} on {dtStr}: {weather_df}")
    end_time_total = time.time()
    total_duration = end_time_total - start_time_total
    print(f'Total function execution took {total_duration:.2f} seconds')

    return weather_df
"""

# NOAA_FORECASTED
# @memoize_noaa_forecasted
# @memoize_noaa_forecasted
def getWeatherFromNOAAFORECASTED(dtStr, agstack_geoid):
    filePath = "/mnt/md2/NOAA/DAILY/FORECASTED/PARQUET_S2/"

    # Parse the date
    tok = dtStr.split('-')
    YYYY_str = tok[0]
    MM_str = tok[1].zfill(2)
    DD_str = tok[2]
    
    start_time = time.time()


    # Create a datetime object and localize it to UTC
    local_dt = datetime(int(YYYY_str), int(MM_str), int(DD_str))
    utc_dt = pytz.utc.localize(local_dt)  # Convert to UTC

    date_range = [utc_dt + timedelta(days=i) for i in range(8)]

    w_ret = pd.DataFrame()
    # Fetch WKT polygon and extract lat/lon
    wkt_polygon = fetchWKT(agstack_geoid)
    lat, lon = extractLatLonFromWKT(wkt_polygon)

    # Filter data by level 9 tokens
    s2_index__L5_list,cell_id = get_s2_cellids_and_token_list(5, [lat], [lon])
    
    list_of_5_paths = [filePath + 's2_token_L5=' + x for x in s2_index__L5_list
            if os.path.exists(filePath + 's2_token_L5=' + x)]

    if not list_of_5_paths:
        return empty_response()
    
    weather_datasets = []
    for x in list_of_5_paths:
        weather_datasets.append(ds.dataset(x, format="parquet", partitioning="hive"))
    w_all = pd.DataFrame()
    for weatherDataset in weather_datasets:
        for current_date in date_range:
            YYYY_list = [current_date.year]
            MM_list = [str(current_date.month).zfill(2)]
            DD_list = [str(current_date.day).zfill(2)]

            weather_df = weatherDataset.to_table(
                filter=(
                    ds.field('Year').isin(YYYY_list) &
                    ds.field('Month').isin(MM_list) &
                    ds.field('Day').isin(DD_list)
                )
            ).to_pandas()
            
            w_all = pd.concat([w_all, weather_df], ignore_index=True)
    
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
        
    end_time = time.time()
    time_elapsed = (end_time - start_time)

    print(f'total time elapsed --{time_elapsed}')
    
    return w_ret

# GHCND
def getWeatherFromGHCND(agstack_geoid, start_date,end_date):
    # filePath = '/mnt/md1/NLDAS/PARQUETE_S2/'
    # filePath = '/home/rajat/Downloads/Rnaura_Work/mnt/md1/GHCND/DAILY/PROCESSED/PARQUET_S2/'
    filePath = '/home/rnaura/mnt/md1/GHCND/DAILY/PROCESSED/PARQUET_S2/'
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
    #     w_ret.reset_index(inplace=True)
    # else:
    #     w_ret = pd.DataFrame()

    return w_ret
"""
@memoize_ghcnd
def getWeatherFromGHCND(dtStr, agstack_geoid):
    # logging.info(f"Data Processing started for geoid: {agstack_geoid}, date: {dtStr}")

    # Parse the date
    tok = dtStr.split('-')
    YYYY_str = tok[0]
    MM_str = tok[1].zfill(2)
    DD_str = tok[2]

    start_time_total = time.time()

    day_key = f'Day={DD_str}'

    # Reload data if necessary
    start_time_reload = time.time()
    day_data = reload_data_ghcnd(YYYY_str, MM_str, DD_str)
    reload_duration = time.time() - start_time_reload
    # logging.info(f"Data loading into memory took {reload_duration:.2f} seconds")

    # Read and combine Parquet files using threading
    try:
        with concurrent.futures.ThreadPoolExecutor() as executor:
            tables = list(executor.map(read_parquet_file, day_data[day_key]))

        combined_table = pa.concat_tables(tables)
        cached_weather_df_ghcnd = combined_table.to_pandas()
        # Set multi-level index for fast lookup
        cached_weather_df_ghcnd.set_index(['s2_token_L9', 's2_token_L7', 's2_token_L5'], inplace=True)
    except Exception as e:
        logging.error(f"Error reading Parquet files: {e}")
        return empty_response()
    
    try:

        # Fetch WKT polygon and extract lat/lon
        wkt_polygon = fetchWKT(agstack_geoid)
        lat, lon = extractLatLonFromWKT(wkt_polygon)
    except:
        pass

    try:

        # Filter data by level 9 tokens
        s2_index__L9_list, _ = get_s2_cellids_and_token_list(9, [lat], [lon])
        s2_index__L9_list = [str(token) for token in s2_index__L9_list]
    except:
        pass

    # Filter by level 9 tokens
    try:
        weather_df = cached_weather_df_ghcnd.loc[
            (s2_index__L9_list, slice(None), slice(None))
        ]
    except KeyError:
        weather_df = pd.DataFrame()

    # If level 9 data is empty, filter by level 7 tokens
    if weather_df.empty:
        level_7_tokens = get_level_7_tokens(s2_index__L9_list)
        try:
            weather_df = cached_weather_df_ghcnd.loc[
                (slice(None), level_7_tokens, slice(None))
            ]
        except KeyError as e:
            print(f'{e} at level 7')
            weather_df = pd.DataFrame()
            

    # If level 7 data is empty, filter by level 5 tokens
    if weather_df.empty:
        level_5_tokens = get_level_5_tokens(s2_index__L9_list)
        try:
            weather_df = cached_weather_df_ghcnd.loc[
                (slice(None), slice(None), level_5_tokens)
            ]
        except KeyError as e:
            print(f'{e} at level 5')
            weather_df = pd.DataFrame()

    # Filter the DataFrame by date
    try:
        weather_df = weather_df[(weather_df['Year'].astype(int) == int(YYYY_str)) & 
                    (weather_df['Month'].astype(int) == int(MM_str)) & 
                    (weather_df['Day'].astype(int) == int(DD_str))]
    except Exception as e:
        print(e)
    
    # print(weather_df['Date'])
    if len(weather_df) > 1:
   
        # Convert 'TS' to datetime
        weather_df['Date'] = pd.to_datetime(weather_df['Date'], format='%Y-%m-%dT%H:%M:%S', errors='coerce')
        # Identify numeric columns for aggregation
        numeric_cols = weather_df.select_dtypes(include=['number']).columns
        non_numeric_cols = weather_df.select_dtypes(exclude=['number']).columns

        # Aggregate numeric columns by 'TS'
        weather_df_numeric = weather_df[numeric_cols].groupby(weather_df['Date']).mean().reset_index()

        # Merge non-numeric columns
        weather_df_non_numeric = weather_df[non_numeric_cols].drop_duplicates(subset='Date')
        weather_df = pd.merge(weather_df_numeric, weather_df_non_numeric, on='Date', how='left')
        
        # Convert Year, Month, Day columns to integers
        weather_df['Year'] = weather_df['Year'].astype(int)
        weather_df['Month'] = weather_df['Month'].astype(int)
        weather_df['Day'] = weather_df['Day'].astype(int)

        # Format 'TS' column
        # weather_df['Date'] = weather_df['Date'].dt.strftime('%Y-%m-%dT%H:%M:%S')
        weather_df['Date'] = pd.to_datetime(weather_df['Date'], errors='coerce')
        weather_df['aifstime_utc'] = pd.to_datetime(weather_df['aifstime_utc'], errors='coerce')
        weather_df['local_date_time_full'] = pd.to_datetime(weather_df['local_date_time_full'], errors='coerce')
    # Return empty response if no data is found
    if weather_df.empty:
        weather_df = pd.DataFrame(empty_response())
    try:

        # Remove NaN values
        columns_to_drop = [
            'AWDR','SNOW', 'SNWD', 'TOBS', 'WESD', 'WSFI','StnID']
        weather_df = weather_df.drop(columns=columns_to_drop, errors='ignore')
        weather_df.dropna(inplace=True)
    except Exception as e:
        print(e)
    
    end_time_total = time.time()
    total_duration = end_time_total - start_time_total
    print(f'Total function execution took {total_duration:.2f} seconds')
    return weather_df
"""
from datetime import timedelta
# AUS_FORECASTED
def getWeatherFromAUSFORECASTED(agstack_geoid, dtStr):
    filePath = "/home/rnaura/mnt/md1/AUS/DAILY/FORECASTED/PARQUET_S2/"
    tok = dtStr.split('-')
    YYYY_str = tok[0]
    MM_str = tok[1]
    DD_str = tok[2]

    # Create a datetime object and localize it to UTC
    local_dt = datetime(int(YYYY_str), int(MM_str), int(DD_str))
    utc_dt = pytz.utc.localize(local_dt)  # Convert to UTC

    # Generate the next 7 days
    date_range = [utc_dt + timedelta(days=i) for i in range(8)]

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
        for current_date in date_range:
            YYYY_list = [current_date.year]
            MM_list = [str(current_date.month).zfill(2)]
            DD_list = [str(current_date.day).zfill(2)]

            weather_df = weatherDataset.to_table(
                filter=(
                    ds.field('Year').isin(YYYY_list) &
                    ds.field('Month').isin(MM_list) &
                    ds.field('Day').isin(DD_list)
                )
            ).to_pandas()

            # Ensure TS column is created
            if not weather_df.empty:
                weather_df['TS'] = pd.to_datetime(weather_df[['Year', 'Month', 'Day']].astype(str).agg('-'.join, axis=1)).dt.tz_localize('UTC')

                numeric_cols = weather_df.select_dtypes(include=[np.number]).columns
                w_df = pd.DataFrame(weather_df[numeric_cols].mean(axis=0)).T

                # Retain the first TS value or handle TS appropriately
                if 'TS' in weather_df.columns and not weather_df['TS'].empty:
                    w_df['TS'] = weather_df['TS'].iloc[0]
                else:
                    w_df['TS'] = pd.NaT  # or some other default value

                w_all = pd.concat([w_all, w_df], ignore_index=True)

    print("=======w_all======", w_all)
    if not w_all.empty:
        w_all = w_all.fillna(0)
        # Take the average of the columns
        w_ret = pd.DataFrame(w_all.groupby(['TS']).mean())
        w_ret['Year'] = w_ret['Year'].astype(int)
        w_ret['Month'] = w_ret['Month'].astype(int)
        w_ret['Day'] = w_ret['Day'].astype(int)
        w_ret.reset_index(inplace=True)
    else:
        w_ret = pd.DataFrame()

    end_time = time.time()
    time_elapsed = (end_time - start_time)
    print(f"Time elapsed: {time_elapsed} seconds")

    return w_ret
"""
@memoize_aus_forecasted
def getWeatherFromAUSFORECASTED(dtStr, agstack_geoid):
    # logging.info(f"Data Processing started for geoid: {agstack_geoid}, date: {dtStr}")

    # Parse the date
    tok = dtStr.split('-')
    YYYY_str = tok[0]
    MM_str = tok[1].zfill(2)
    DD_str = tok[2]

    global cached_weather_files_aus_forecasted, last_reload_time_aus_forecasted, cached_weather_df_aus_forecasted
    start_time_total = time.time()

    day_key = f'Day={DD_str}'

    # Reload data if necessary
    start_time_reload = time.time()
    day_data = reload_data_aus_forecasted(YYYY_str, MM_str, DD_str)
    reload_duration = time.time() - start_time_reload
    # logging.info(f"Data loading into memory took {reload_duration:.2f} seconds")

    # Read and combine Parquet files using threading
    try:
        with concurrent.futures.ThreadPoolExecutor() as executor:
            tables = list(executor.map(read_parquet_file, day_data[day_key]))

        combined_table = pa.concat_tables(tables)
        cached_weather_df_aus_forecasted = combined_table.to_pandas()

        # Remove NaN values
        cached_weather_df_aus_forecasted.dropna(inplace=True)

        # Set multi-level index for fast lookup
        cached_weather_df_aus_forecasted.set_index(['s2_token_L9', 's2_token_L7', 's2_token_L5'], inplace=True)
    except Exception as e:
        logging.error(f"Error reading Parquet files: {e}")
        return empty_response()
    
    # Fetch WKT polygon and extract lat/lon
    wkt_polygon = fetchWKT(agstack_geoid)
    lat, lon = extractLatLonFromWKT(wkt_polygon)

    # Filter data by level 9 tokens
    s2_index__L9_list, _ = get_s2_cellids_and_token_list(9, [lat], [lon])
    s2_index__L9_list = [str(token) for token in s2_index__L9_list]

    # Filter by level 9 tokens
    try:
        weather_df = cached_weather_df_aus_forecasted.loc[
            (s2_index__L9_list, slice(None), slice(None))
        ]
    except KeyError:
        weather_df = pd.DataFrame()

    # If level 9 data is empty, filter by level 7 tokens
    if weather_df.empty:
        level_7_tokens = get_level_7_tokens(s2_index__L9_list)
        try:
            weather_df = cached_weather_df_aus_forecasted.loc[
                (slice(None), level_7_tokens, slice(None))
            ]
        except KeyError as e:
            print(f'{e} at level 7')
            weather_df = pd.DataFrame()
            
    # If level 7 data is empty, filter by level 5 tokens
    if weather_df.empty:
        level_5_tokens = get_level_5_tokens(s2_index__L9_list)
        try:
            weather_df = cached_weather_df_aus_forecasted.loc[
                (slice(None), slice(None), level_5_tokens)
            ]
        except KeyError as e:
            print(f'{e} at level 5')
            weather_df = pd.DataFrame()

    # Filter the DataFrame by date
    try:
        weather_df = weather_df[(weather_df['Year'].astype(int) == int(YYYY_str)) & 
                    (weather_df['Month'].astype(int) == int(MM_str)) & 
                    (weather_df['Day'].astype(int) == int(DD_str))]
    except Exception as e:
        print(e)

    # Process time series data if multiple entries exist
    if len(weather_df) > 1:
        # Convert 'TS' to datetime
        weather_df['TS'] = pd.to_datetime(weather_df['TS'], errors='coerce')

        # Identify numeric columns for aggregation
        numeric_cols = weather_df.select_dtypes(include=['number']).columns
        non_numeric_cols = weather_df.select_dtypes(exclude=['number']).columns

        # Aggregate numeric columns by 'TS'
        weather_df_numeric = weather_df[numeric_cols].groupby(weather_df['TS']).mean().reset_index()

        # Merge non-numeric columns
        weather_df_non_numeric = weather_df[non_numeric_cols].drop_duplicates(subset='TS')
        weather_df = pd.merge(weather_df_numeric, weather_df_non_numeric, on='TS', how='left')

        # Ensure 'Year', 'Month', 'Day' columns are present after aggregation
        weather_df['Year'] = weather_df['TS'].dt.year
        weather_df['Month'] = weather_df['TS'].dt.month
        weather_df['Day'] = weather_df['TS'].dt.day

        # Convert Year, Month, Day columns to integers
        weather_df['Year'] = weather_df['Year'].astype(int)
        weather_df['Month'] = weather_df['Month'].astype(int)
        weather_df['Day'] = weather_df['Day'].astype(int)

        # Format 'TS' column
        # weather_df['TS'] = weather_df['TS'].dt.strftime('%Y-%m-%dT%H:%M:%S')

    # Return empty response if no data is found
    if weather_df.empty:
        weather_df = pd.DataFrame(empty_response())
    
    try:
        weather_df = weather_df.reset_index()
        weather_df = weather_df.drop(columns=['s2_token_L9', 's2_token_L7', 's2_token_L5'], errors='ignore')
        weather_df = weather_df.drop(columns=['geometry', 'region'], errors='ignore')
        print("Columns after dropping:", weather_df.columns)
    except Exception as e:
        print(f"Error dropping columns: {e}")

    logging.info(f"AUS FORECASTED Data for geoid {agstack_geoid} on {dtStr}: {weather_df}")
    end_time_total = time.time()
    total_duration = end_time_total - start_time_total
    logging.info(f'Total function execution took {total_duration:.2f} seconds')
    print(weather_df)
    return weather_df
"""
# AUS
# AUS
# AUS
def getWeatherFromAUS(agstack_geoid, dtStr,end_date=None):
    # filePath = '/mnt/md1/NLDAS/PARQUETE_S2/'
    # filePath = '/home/rajat/Downloads/Rnaura_Work/mnt/md1/AUS/DAILY/PROCESSED/PARQUET_S2_TEST/'
    filePath = "/home/rnaura/mnt/md1/AUS/DAILY/PROCESSED/PARQUET_S2/"

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
        w_all = pd.concat([w_all, weather_df], ignore_index=True)
    if len(w_all) > 0:
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

    print("=======w_all======", w_all)
    # if not w_all.empty:
    #     w_all = w_all.fillna(0)
    #     # Take the average of the columns
    #     w_ret = pd.DataFrame(w_all.groupby(['aifstime_utc']).mean())
    #     w_ret.reset_index(inplace=True)
    # else:
    #     w_ret = pd.DataFrame()

    end_time = time.time()
    time_elapsed = (end_time - start_time)

    return w_ret

"""
@memoize_aus
def getWeatherFromAUS(dtStr, agstack_geoid):
    # logging.info(f"Data Processing started for geoid: {agstack_geoid}, date: {dtStr}")

    # Parse the date
    tok = dtStr.split('-')
    YYYY_str = tok[0]
    MM_str = tok[1].zfill(2)
    DD_str = tok[2]

    global cached_weather_files_aus, last_reload_time_aus, cached_weather_df_aus
    start_time_total = time.time()

    day_key = f'Day={DD_str}'

    # Reload data if necessary
    start_time_reload = time.time()
    day_data = reload_data_aus(YYYY_str, MM_str, DD_str)
    reload_duration = time.time() - start_time_reload
    # logging.info(f"Data loading into memory took {reload_duration:.2f} seconds")

    # Read and combine Parquet files using threading
    try:
        with concurrent.futures.ThreadPoolExecutor() as executor:
            tables = list(executor.map(read_parquet_file, day_data[day_key]))

        combined_table = pa.concat_tables(tables)
        cached_weather_df_aus = combined_table.to_pandas()

        # Set multi-level index for fast lookup
        cached_weather_df_aus.set_index(['s2_token_L9', 's2_token_L7', 's2_token_L5'], inplace=True)
    except Exception as e:
        logging.error(f"Error reading Parquet files: {e}")
        return empty_response()
    
    # Fetch WKT polygon and extract lat/lon
    wkt_polygon = fetchWKT(agstack_geoid)
    lat, lon = extractLatLonFromWKT(wkt_polygon)

    # Filter data by level 9 tokens
    s2_index__L9_list, _ = get_s2_cellids_and_token_list(9, [lat], [lon])
    s2_index__L9_list = [str(token) for token in s2_index__L9_list]

    # Filter by level 9 tokens
    try:
        weather_df = cached_weather_df_aus.loc[
            (s2_index__L9_list, slice(None), slice(None))
        ]
    except KeyError:
        weather_df = pd.DataFrame()

    # If level 9 data is empty, filter by level 7 tokens
    if weather_df.empty:
        level_7_tokens = get_level_7_tokens(s2_index__L9_list)
        try:
            weather_df = cached_weather_df_aus.loc[
                (slice(None), level_7_tokens, slice(None))
            ]
        except KeyError as e:
            print(f'{e} at level 7')
            weather_df = pd.DataFrame()
            
    # If level 7 data is empty, filter by level 5 tokens
    if weather_df.empty:
        level_5_tokens = get_level_5_tokens(s2_index__L9_list)
        try:
            weather_df = cached_weather_df_aus.loc[
                (slice(None), slice(None), level_5_tokens)
            ]
        except KeyError as e:
            print(f'{e} at level 5')
            weather_df = pd.DataFrame()

    # Filter the DataFrame by date
    weather_df = weather_df[
        (weather_df['Year'].astype(int) == int(YYYY_str)) & 
        (weather_df['Month'].astype(int) == int(MM_str)) & 
        (weather_df['Day'].astype(int) == int(DD_str))
    ]

    # Check if there's data to process
    if not weather_df.empty:
        # Convert 'aifstime_utc' to datetime
        weather_df['aifstime_utc'] = pd.to_datetime(weather_df['aifstime_utc'], errors='coerce')

        # Extract date part (ignoring time)
        weather_df['date'] = weather_df['aifstime_utc'].dt.date

        # Aggregate numeric columns by date
        numeric_cols = weather_df.select_dtypes(include=['number']).columns
        weather_df_numeric = weather_df.groupby('date')[numeric_cols].mean().reset_index()

        # For non-numeric columns, take the first non-null value
        non_numeric_cols = weather_df.select_dtypes(exclude=['number']).columns.difference(['aifstime_utc', 'date'])
        weather_df_non_numeric = weather_df.groupby('date')[non_numeric_cols].first().reset_index()

        # Combine numeric and non-numeric columns
        weather_df_combined = pd.merge(weather_df_numeric, weather_df_non_numeric, on='date')

        # Filter out the result for the specific date
        weather_df_result = weather_df_combined[
            (weather_df_combined['date'] == pd.to_datetime(f'{YYYY_str}-{MM_str}-{DD_str}').date())
        ]

        # Drop the temporary 'date' column
        weather_df_result = weather_df_result.drop('date', axis=1)

        # Ensure the result is a single-row DataFrame
        if len(weather_df_result) == 0:
            weather_df_result = pd.DataFrame(empty_response())
    else:
        weather_df_result = pd.DataFrame(empty_response())

    try:
        weather_df_result = weather_df_result.drop(['geometry', 'region'], axis=1, errors='ignore')
    except:
        pass

    logging.info(f"AUS Data for geoid {agstack_geoid} on {dtStr}: {weather_df_result}")
    end_time_total = time.time()
    total_duration = end_time_total - start_time_total
    logging.info(f'Total function execution took {total_duration:.2f} seconds')

    print(weather_df_result)
    return weather_df_result
"""

def reload_gunicorn_workers():
    print("Reloading Gunicorn workers...")
    # Send the HUP signal to Gunicorn to gracefully reload the workers
    os.kill(os.getppid(), signal.SIGHUP)  # Get Gunicorn's parent process ID and send SIGHUP

def reload_data():
    print("Reloading data...")
    
    # Reload Gunicorn workers
    reload_gunicorn_workers()




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

# Function to start monitoring for the current day



cache_nldas = {}
@app.route('/')
#def hello_world():
 #   return 'Hello World!'
@app.route('/getNLDAS')
def getNLDASData():
    geoid = request.args['geoid']
    # s2_token = ["8094c20c","8094e9e4"]
    date = request.args.get('date', '')
    # Start timing for cache access
    start_time_cache_access = time.time()
    
    # Data not in cache, fetch from NCEP
    weather_df = getWeatherFromNLDAS(geoid,date)
    
    try:
        weather_df.reset_index(inplace=True)
    except:
        # logging.warning(f"Resetting index failed for dataframes with geoid: {geoid}")
        pass

    # Convert DataFrame to a dictionary
    json_data = weather_df.to_dict(orient='records')


    # Log the completion of the request
    # logging.info(f"Data for geoid: {geoid} fetched from NCEP and cached")
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

    # Start timing for cache access
    start_time_cache_access = time.time()
   

    if end_date is not None :

        weather_df = getWeatherFromNOAA(geoid, start_date, end_date)
    else:
        weather_df = getWeatherFromNOAA(geoid, start_date, end_date=None)
    
    try:
        weather_df.reset_index(inplace=True)
    except:
        # logging.warning(f"Resetting index failed for dataframe with geoid: {geoid}")
        pass

    # Convert DataFrame to a dictionary
    json_data = weather_df.to_dict(orient='records')


    response = {
        "data": json_data,
        "metadata": get_noaa_metadata(),
        "metadata-description":get_noaa_description()
    }

    # Log the completion of the request
    # logging.info(f"Data for geoid: {geoid} fetched from NOAA and cached")

    return jsonify(response)

#NOAA_FORECASTED
cache_noaa_forecasted = {}
@app.route('/getNOAAForecast')
def getNOAAFORECASTEDWeatherData():
    geoid = request.args['geoid']
    date = request.args.get('date', '')

    # Start timing for cache access
    start_time_cache_access = time.time()


    # Data not in cache, fetch from NOAA FORECASTED
    weather_df = getWeatherFromNOAAFORECASTED(date, geoid)
    
    try:
        weather_df.reset_index(inplace=True)
    except:
        logging.warning(f"Resetting index failed for dataframe with geoid: {geoid}")
        pass

    # Convert DataFrame to a dictionary
    json_data = weather_df.to_dict(orient='records')

    

    # Log the completion of the request
    # logging.info(f"Data for geoid: {geoid} fetched from NOAA FORECASTED and cached")
    response = {
        "data": json_data,
        "metadata": meta_data_response(),
        "metadata-description":meta_data_with_description()
    }
    return jsonify(response)

#GHCND
cache_ghcnd = {}
@app.route('/getGHCND')
def getGHCNDWeatherData():
    geoid = request.args['geoid']
    start_date = request.args.get('start_date', '')
    end_date = request.args.get('end_date', None)


    
    weather_df = getWeatherFromGHCND(geoid, start_date,end_date)
    
    try:
        weather_df.reset_index(inplace=True)
    except:
        # logging.warning(f"Resetting index failed for dataframe with geoid: {geoid}")
        pass

    # Convert DataFrame to a dictionary
    json_data = weather_df.to_dict(orient='records')

    # Cache the response
    # cache_ghcnd[cache_key_ghcnd] = json_data

    # Log the completion of the request
    # logging.info(f"Data for geoid: {geoid} fetched from GHCND and cached")

    response = {
        "data": json_data,
        "metadata": get_GHCND_metdata(),
        "metadata-description":get_GHCND_description()
    }

    return jsonify(response)

#AUS_FORECASTED
cache_aus_forecasted = {}
@app.route('/getAUSFORECASTED')
def getAUSFORECASTEDWeatherData():
    geoid = request.args['geoid']
    date = request.args.get('date', '')

    # Start timing for cache access
    start_time_cache_access = time.time()
    cache_key_aus_forecasted = ('AUS_FORECASTED',geoid, date)
    if cache_key_aus_forecasted in cache_aus_forecasted:
        # Log cache access time
        cache_access_duration = time.time() - start_time_cache_access
        logging.info(f"Returning cached data for geoid: {geoid} took {cache_access_duration:.2f} seconds")
        return jsonify(cache_aus_forecasted[cache_key_aus_forecasted])

    # Data not in cache, fetch from AUS FORECASTED
    weather_df = getWeatherFromAUSFORECASTED(geoid, date)
    
    try:
        weather_df.reset_index(inplace=True)
    except:
        logging.warning(f"Resetting index failed for dataframe with geoid: {geoid}")
        pass

    # Convert DataFrame to a dictionary
    json_data = weather_df.to_dict(orient='records')

    # Cache the response
    cache_aus_forecasted[cache_key_aus_forecasted] = json_data

    # Log the completion of the request
    logging.info(f"Data for geoid: {geoid} fetched from AUS FORECASTED and cached")

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

# @app.route('/getNDVIImg')
# def getNDVIImg():
#     dtStr =  request.args.get('dt', '')
#     agstack_geoid =  request.args.get('geoid', '')
#     final_res_gdf, bound_gdf, returnedDtStr = getNDVIgdfFromGeoid(agstack_geoid, dtStr)
#     imgDict = {}
#     imgDict['type']='Feature'
#     imgDict['image_geoDataFrameDict']=final_res_gdf.to_json()
#     imgDict['boundary_geoDataFrameDict']=bound_gdf.to_json()
#     imgDict['imageDt']=returnedDtStr

#     return imgDict

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
# base_path_noaa_reload = '/home/rnaura/Rajat/noaa/NOAA/DAILY/PROCESSED/S2_GEOJSON/PARQUET_S2/s2_tokens_L5/s2_tokens_L7/s2_tokens_L9/s2_tokens_L20/Year=2024/'
base_path_noaa_reload = '/home/rnaura/mnt/md1/NOAA/DAILY/PROCESSED/S2_GEOJSON/PARQUET_S2/s2_tokens_L5/s2_tokens_L7/s2_tokens_L9/s2_tokens_L20/Year=2024/'
# base_path_noaa_forecasted_reload = '/home/rajat/Downloads/Rnaura/mnt/md1/NOAA/DAILY/FORECASTED/PARQUET_S2/s2_tokens_L5/s2_tokens_L7/s2_tokens_L9/s2_tokens_L20/Year=2024/'
base_path_noaa_forecasted_reload = '/home/rnaura/mnt/md1/NOAA/DAILY/FORECASTED/PARQUET_S2/s2_tokens_L5/s2_tokens_L7/s2_tokens_L9/s2_tokens_L20/Year=2024/'
# base_path_ghcnd_reload = '/home/rajat/Downloads/Rnaura/mnt/md1/GHCND/PARQUET_S2/s2_tokens_L5/s2_tokens_L7/s2_tokens_L9/Year=2024/'
base_path_ghcnd_reload = '/home/rnaura/mnt/md1/GHCND/PARQUET_S2/s2_tokens_L5/s2_tokens_L7/s2_tokens_L9/Year=2024/'
# base_path_aus_forecasted_reload = '/home/rajat/Downloads/Rnaura_Work/mnt/md1/AUS/DAILY/FORECASTED/PARQUET_S2/s2_tokens_L5/s2_tokens_L7/s2_tokens_L9/Year=2024/'
base_path_aus_forecasted_reload = '/home/rnaura/mnt/md1/AUS/DAILY/FORECASTED/PARQUET_S2/s2_tokens_L5/s2_tokens_L7/s2_tokens_L9/Year=2024/'
# base_path_aus_reload = '/home/rajat/Downloads/Rnaura_Work/mnt/md1/AUS/DAILY/PROCESSED/PARQUET_S2/s2_tokens_L5/s2_tokens_L7/s2_tokens_L9/Year=2024/'
base_path_aus_reload = '/home/rnaura/mnt/md1/AUS/DAILY/PROCESSED/PARQUET_S2/s2_tokens_L5/s2_tokens_L7/s2_tokens_L9/Year=2024/'



extra_dirs = [base_path_nldas_reload, base_path_ncep_reload,base_path_noaa_reload,base_path_noaa_forecasted_reload,base_path_ghcnd_reload,base_path_aus_forecasted_reload,base_path_aus_reload]
extra_files = extra_dirs[:]

for extra_dir in extra_dirs:
    for dirname, dirs, files in walk(extra_dir):
        for filename in files:
            file_path = path.join(dirname, filename)
            if path.isfile(file_path):
                extra_files.append(file_path)


# @app.before_request
# def setup():
#     base_path_ncep = '/mnt/md1/NCEP/PARQUET_S2/'
#     YYYY_str = time.strftime('%Y')
#     MM_str = time.strftime('%m').zfill(2)
#     DD_str = time.strftime('%d').zfill(2)
#     month_path = os.path.join(base_path_ncep, 's2_tokens_L5', 's2_tokens_L7', 's2_tokens_L9', f'Year={YYYY_str}', f'Month={MM_str}')
    
#     # Start the monitoring thread
#     threading.Thread(target=start_monitoring, args=(month_path, reload_data)).start()


if __name__ == '__main__':
    # extra_files = [updated_data_available_file,]
    print("Loading data...")
    # from threading import Thread
    # monitoring_thread = Thread(target=start_monitoring)
    # monitoring_thread.daemon = True
    # monitoring_thread.start()

    # Load all Parquet files into cache initially
    # load_all_data()
    print("Starting Flask app...")
    # logging.basicConfig(level=logging.INFO)
    app.run(debug=True, extra_files=extra_files,port=5000)



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