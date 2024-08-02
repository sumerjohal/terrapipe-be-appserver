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

logging.basicConfig(
    filename='weather_data_processing.log',  # Log file name
    level=logging.INFO,  # Log level (INFO, DEBUG, ERROR, etc.)
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',  # Log message format
    datefmt='%Y-%m-%d %H:%M:%S'  # Date format
)# cache = Cache(app)




app = Flask(__name__)
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

last_update_time = None
# cached_weather_df = pd.DataFrame()  # Placeholder for cached data
last_reload_time = None

cached_weather_df_nldas = pd.DataFrame()
last_reload_time_nldas = None

cached_weather_df_ncep = pd.DataFrame()
cached_weather_dirs_ncep= {}
last_reload_time_ncep = None
cache_ncep = {}
region_weather_cache = {}
region_last_reload_time = {}
cached_day_files = {}
cached_weather_files_ncep = {}


base_path_nldas = '/mnt/md1/NLDAS/PARQUETE_S2/'
base_path_ncep = '/mnt/md1/NCEP/PARQUET_S2/'

class FileChangeHandler(FileSystemEventHandler):
    def __init__(self, reload_function):
        self.reload_function = reload_function

    def on_modified(self, event):
        if event.src_path.endswith('.parquet'):
            print(f"File modified: {event.src_path}")
            self.reload_function()

    def on_created(self, event):
        if event.src_path.endswith('.parquet'):
            print(f"File created: {event.src_path}")
            self.reload_function()

def start_file_monitoring(path, reload_function):
    event_handler = FileChangeHandler(reload_function)
    observer = Observer()
    observer.schedule(event_handler, path, recursive=True)
    observer_thread = threading.Thread(target=observer.start)
    observer_thread.daemon = True
    observer_thread.start()
    print(f"Started monitoring {path} for changes")

def reload_data_nldas(YYYY_str, MM_str):
    global cached_weather_df_nldas, last_reload_time_nldas
    try:
        start_time = time.time()
        month_path = os.path.join(base_path_nldas, 's2_tokens_L5', 's2_tokens_L7', 's2_tokens_L9', f'Year={YYYY_str}', f'Month={MM_str}')
        
        all_files = []
        for day_dir in os.listdir(month_path):
            day_path = os.path.join(month_path, day_dir)
            if os.path.isdir(day_path):
                all_files.extend([os.path.join(day_path, f) for f in os.listdir(day_path) if f.endswith('.parquet')])
        
        tables = [pq.read_table(file) for file in all_files]
        combined_table = pa.concat_tables(tables)
        cached_weather_df_nldas = combined_table.to_pandas()
        last_reload_time_nldas = time.time()
        
        print(f"Data reloaded in {time.time() - start_time:.4f} seconds")
    except Exception as e:
        print(f"Failed to reload data: {e}")

def reload_data_ncep(YYYY_str, MM_str, Day_str):
    cached_weather_files_ncep = {}
    try:
        start_time = time.time()
        month_path = os.path.join(base_path_ncep, 's2_tokens_L5', 's2_tokens_L7', 's2_tokens_L9', f'Year={YYYY_str}', f'Month={MM_str}')
        
        for day_dir in os.listdir(month_path):
            day_path = os.path.join(month_path, day_dir)
            if os.path.isdir(day_path):
                day_files = [os.path.join(day_path, f) for f in os.listdir(day_path) if f.endswith('.parquet')]
                if day_files:
                    cached_weather_files_ncep[day_dir] = day_files

        print(f"Directories reloaded in {time.time() - start_time:.4f} seconds")
    except Exception as e:
        print(f"Failed to reload directories: {e}")

    return cached_weather_files_ncep
world = gpd.read_file('/home/rnaura/terrapipe/ne_110m_admin_0_countries.shp')


last_region = None
def determine_region(lat, lon):
    point = Point(lon, lat)
    # Finding the row in the dataframe where the point lies within the country's boundary
    country = world[world.contains(point)]
    if not country.empty:
        # Extracting the first row as the result
        return country.iloc[0]['SOVEREIGNT']  # or another column name if different
    else:
        return 'Unknown'



def getWeatherFromNLDAS(dtStr, agstack_geoid):
    
    """    
    This API filters data from a PARQUET file based on S2 tokens and date. It retrieves polygons 
    from the asset registry and fetches the S2 tokens according to the specified level and polygons.
    params: date string , Geoid
    """

    configFile = '/home/rnaura/terrapipe/config/11SKA/11SKA.json'
    tok = dtStr.split('-')
    YYYY_str = tok[0]
    MM_str = tok[1].zfill(2)
    DD_str = tok[2]

    # Load and validate the JSON file
    with open(configFile, "r") as jsonfile:
        fieldJSON = json.load(jsonfile)
        field_geoid = fieldJSON['geoid']
        field_wkt = fieldJSON['wkt']

    # Load the WKT and calculate centroid
    fieldPoly = wkt.loads(field_wkt)
    c = fieldPoly.centroid
    lat = c.y
    lon = c.x
    global cached_weather_df_nldas, last_reload_time_nldas
    
    

    lats, lons = [lat], [lon]
    s2_index__L9_list, _ = get_s2_cellids_and_token_list(9, lats, lons) 

    start_time_total = time.time()

    if last_reload_time_nldas is None or time.time() - last_reload_time_nldas > 60 * 60* 5:
        reload_data_nldas(YYYY_str, MM_str)   
    
    weather_df_filtered = pd.DataFrame()
    if not cached_weather_df_ncep.empty:

        # Check if data needs reloading based on last_reload_time
        if 's2_token_L9' in cached_weather_df_nldas.columns and not cached_weather_df_nldas.empty:
            cached_weather_df_nldas = cached_weather_df_nldas[cached_weather_df_nldas['s2_token_L9'].isin(s2_index__L9_list)]
        
        if not cached_weather_df_nldas.empty:
            level_5_tokens = get_level_5_tokens(s2_index__L9_list)
            weather_df = cached_weather_df_nldas[cached_weather_df_nldas['s2_token_L5'].isin(level_5_tokens)]
            print(weather_df)
        if not cached_weather_df_nldas.empty:
            level_7_tokens = get_level_7_tokens(s2_index__L9_list)
            weather_df = cached_weather_df_nldas[cached_weather_df_nldas['s2_token_L7'].isin(level_7_tokens)]
            
        try:
            weather_df['Year'] = weather_df['Year'].astype(int)
            weather_df['Month'] = weather_df['Month'].astype(int)
            weather_df['Day'] = weather_df['Day'].astype(int)
            
            weather_df_filtered = weather_df[(weather_df['Year'] == int(YYYY_str)) & 
                                             (weather_df['Month'] == int(MM_str)) & 
                                             (weather_df['Day'] == int(DD_str))]
            
            print("After filtering:")
            print(weather_df_filtered.head())
            
        except KeyError as e:
            print(f"KeyError: {e}")
            # weather_df_filtered = pd.DataFrame(columns=weather_df.columns)  # Empty DataFrame with correct columns
        
        end_time_total = time.time()
        total_duration = end_time_total - start_time_total
        print(f"Total request handling time: {total_duration:.4f} seconds")

    else:
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


# logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

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

# # Custom cache key function
# def make_cache_key():
#     args = request.args
#     dtStr = args.get('date')
#     agstack_geoid = args.get('geoid')
#     key_str = f"{dtStr}_{agstack_geoid}"
#     print(f'Cache key generated: {hashlib.md5(key_str.encode()).hexdigest()}')
#     return hashlib.md5(key_str.encode()).hexdigest()

def make_key(*args, **kwargs):
   """A function which is called to derive the key for a computed value.
      The key in this case is the concat value of all the json request
      parameters. Other strategy could to use any hashing function.
   :returns: unique string for which the value should be cached.
   """
   user_data = request.args
   return ",".join([f"{key}={value}" for key, value in user_data.items()])


def memoize(func):
    cache = {}

    @wraps(func)
    def wrapper(*args,**kwargs):
        key = str(args) + str(kwargs)

        if key not in cache:
            cache[key] = func(*args,**kwargs)
        return cache[key]

    return wrapper


# @cache.cached(timeout=5000000, make_cache_key=make_key)
# @log_time_and_calls
@memoize
def getWeatherFromNCEP(dtStr, agstack_geoid):
    # Extract date and geoid from the request
    dtStr = request.args.get('date')
    agstack_geoid = request.args.get('geoid')

    logging.info(f"Data Processing started for geoid: {agstack_geoid}, date: {dtStr}")

    # Parse the date
    tok = dtStr.split('-')
    YYYY_str = tok[0]
    MM_str = tok[1].zfill(2)
    DD_str = tok[2]

    global cached_weather_files_ncep, last_reload_time_ncep, cached_weather_df_ncep
    start_time_total = time.time()

    day_key = f'Day={DD_str}'

    # Reload data if necessary
    start_time_reload = time.time()
    day_data = reload_data_ncep(YYYY_str, MM_str, DD_str)
    reload_duration = time.time() - start_time_reload
    logging.info(f"Data loading into memory took {reload_duration:.2f} seconds")

    try:
        # Read and combine Parquet files
        start_time_reading = time.time()
        tables = [pq.read_table(file) for file in day_data[day_key]]
        reading_duration = time.time() - start_time_reading
        logging.info(f"data table creation took {reading_duration:.2f} seconds")
        
        combined_table = pa.concat_tables(tables)
        count_table = len(combined_table)
        logging.info(f"parquet table count {count_table} ")
        start_time_combined_table= time.time()
        reading_duration = time.time() - start_time_combined_table
        logging.info(f"table concatenation took {reading_duration:.2f} seconds")

        start_time_table_to_dataframe= time.time()
        cached_weather_df_ncep = combined_table.to_pandas()

        reading_duration = time.time() - start_time_table_to_dataframe
        logging.info(f"tables dataframe conversion took {reading_duration:.2f} seconds")

        start_time_multi_indexing= time.time()
        cached_weather_df_ncep.set_index(['s2_token_L9', 's2_token_L7', 's2_token_L5'], inplace=True)
        reading_duration = time.time() - start_time_multi_indexing
        logging.info(f"Reading and combining Parquet files took {reading_duration:.2f} seconds")

    except Exception as e:
        logging.error(f"Error reading parquet files: {e}")
        return empty_response()

    # Fetch WKT polygon and extract lat/lon
    start_time_polygon = time.time()
    wkt_polygon = fetchWKT(agstack_geoid)
    lat, lon = extractLatLonFromWKT(wkt_polygon)
    polygon_duration = time.time() - start_time_polygon
    logging.info(f"Fetching and processing WKT polygon took {polygon_duration:.2f} seconds")

    # Filter data by level 9 tokens
    start_time_filtering_level9 = time.time()
    s2_index__L9_list, _ = get_s2_cellids_and_token_list(9, [lat], [lon])
    s2_index__L9_list = [str(token) for token in s2_index__L9_list]
    try:
        weather_df = cached_weather_df_ncep.loc[(s2_index__L9_list, slice(None), slice(None))]
        filtering_level9_duration = time.time() - start_time_filtering_level9
        logging.info(f"Filtering data for level 9 tokens took {filtering_level9_duration:.2f} seconds")
    except KeyError:
        weather_df = pd.DataFrame()

    # Filter data by level 7 tokens if level 9 data is empty
    if weather_df.empty:
        start_time_filtering_level7 = time.time()
        level_7_tokens = get_level_7_tokens(s2_index__L9_list)
        try:
            weather_df = cached_weather_df_ncep.loc[(slice(None), level_7_tokens, slice(None))]
            filtering_level7_duration = time.time() - start_time_filtering_level7
            logging.info(f"Filtering data for level 7 tokens took {filtering_level7_duration:.2f} seconds")
        except KeyError:
            weather_df = pd.DataFrame()

    # Filter data by level 5 tokens if level 7 data is empty
    if weather_df.empty:
        start_time_filtering_level5 = time.time()
        level_5_tokens = get_level_5_tokens(s2_index__L9_list)
        try:
            weather_df = cached_weather_df_ncep.loc[(slice(None), slice(None), level_5_tokens)]
            filtering_level5_duration = time.time() - start_time_filtering_level5
            logging.info(f"Filtering data for level 5 tokens took {filtering_level5_duration:.2f} seconds")
        except KeyError:
            weather_df = pd.DataFrame()

    # Return empty response if no data is found
    if weather_df.empty:
        weather_df = pd.DataFrame(empty_response())

    end_time_total = time.time()
    total_duration = end_time_total - start_time_total
    logging.info(f"Total function execution took {total_duration:.2f} seconds")

    return weather_df

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
@app.route('/')
#def hello_world():
 #   return 'Hello World!'
@app.route('/getDataFromNLDAS')
def getNLDASData():
    geoid = request.args['geoid']
    # s2_token = ["8094c20c","8094e9e4"]
    date = request.args.get('date', '')
    weather_df = getWeatherFromNLDAS(date,geoid)
    weather_df.reset_index(inplace=True)
    # Convert DataFrame to JSON
    json_data = weather_df.to_json(orient='records')
    return json_data

cache_ncep = {}

@app.route('/getDataFromNCEP')
def getNCEPWeatherData():
    geoid = request.args['geoid']
    date = request.args.get('date', '')

    # Start timing for cache access
    start_time_cache_access = time.time()

    if geoid in cache_ncep:
        # Log cache access time
        cache_access_duration = time.time() - start_time_cache_access
        # logging.info(f"Cache access for geoid: {geoid} took {cache_access_duration:.2f} seconds")

        logging.info(f"Returning cached data for geoid: {geoid} took {cache_access_duration:.2f} seconds")
        return jsonify(cache_ncep[geoid])

    # Data not in cache, fetch from NCEP
    weather_df = getWeatherFromNCEP(date, geoid)
    
    try:    
        weather_df.reset_index(inplace=True)
    except:
        logging.warning(f"Resetting index failed for dataframe with geoid: {geoid}")
        pass

    # Convert DataFrame to JSON
    json_data = weather_df.to_json(orient='records')

    # Cache the response
    cache_ncep[geoid] = json_data

    # Log the completion of the request
    logging.info(f"Data for geoid: {geoid} fetched from NCEP and cached")

    return json_data

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
extra_dirs = [base_path_nldas, base_path_ncep]
extra_files = extra_dirs[:]

for extra_dir in extra_dirs:
    for dirname, dirs, files in walk(extra_dir):
        for filename in files:
            file_path = path.join(dirname, filename)
            if path.isfile(file_path):
                extra_files.append(file_path)

if __name__ == '__main__':
    # extra_files = [updated_data_available_file,]
    logging.basicConfig(level=logging.INFO)
    app.run(debug=True, extra_files=extra_files)



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