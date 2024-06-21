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
from datetime import date
from datetime import date, timedelta, datetime
from dateutil.parser import parse
from shapely.geometry import mapping, Point, Polygon
import pyarrow.dataset as ds
import matplotlib.pyplot as plt
import matplotlib
warnings.filterwarnings('ignore')

################ GLOBALS
#global datasets


################  FUNCTIONS
################  FUNCTIONS
def getExtent(polygon):
    lats = list(polygon.exterior.xy[1])
    lons = list(polygon.exterior.xy[0])
    ext = [min(lons),max(lons),min(lats),max(lats)]
    
    return ext

def isValidPolygon(user_fieldWKT):
    try:
        poly = shapely.wkt.loads(user_fieldWKT)
        if poly.geom_type == 'Polygon':
            return True
        else:
            return False
    except Exception as e:
        print(e)
        return False

def getFieldOutline(user_fieldWKT, s2_max_resolution,maxResColName):
    poly = shapely.wkt.loads(user_fieldWKT)
    (ln, lt) = poly.boundary.xy
    lons=ln.tolist()
    lats=lt.tolist()
    
    min_level=s2_max_resolution
    max_level=s2_max_resolution
    r = s2.RegionCoverer()
    r.min_level = min_level
    r.max_level = max_level
    lb_lat = min(lats)
    ub_lat = max(lats)
    lb_lon = min(lons)
    ub_lon = max(lons)
    lb = s2.LatLng.from_degrees(lb_lat, lb_lon)
    ub = s2.LatLng.from_degrees(ub_lat, ub_lon)
    maxRes_cell_ids = r.get_covering(s2.LatLngRect.from_point_pair(lb,ub))


    s2_LmaxRes_token_list = []
    p_gdf = gpd.GeoDataFrame()
    idx=0
    for cellid in maxRes_cell_ids:
        c = s2.Cell(cellid)
        level = cellid.level()
        vertices = []
        for i in range(0,4):
            vertex = c.get_vertex(i)
            latlng = s2.LatLng.from_point(vertex)
            vertices.append((latlng.lng().degrees, latlng.lat().degrees))
        geo = shapely.geometry.Polygon(vertices)
        if poly.intersects(geo):
            s2_LmaxRes_token_list.append(cellid.to_token())
            #p_gdf.loc[idx,'level'] = level
            p_gdf.loc[idx,maxResColName] = cellid.to_token()
            p_gdf.loc[idx,'geometry']=geo
            idx+=1
        
    return p_gdf
    
def get_s2_cellids_and_token_list_for_wkt(resLevel, user_fieldWKT):
    poly = shapely.wkt.loads(user_fieldWKT)
    (ln, lt) = poly.boundary.xy
    lons=ln.tolist()
    lats=lt.tolist()

    min_level=resLevel
    max_level=resLevel
    r = s2.RegionCoverer()
    r.min_level = min_level
    r.max_level = max_level
    lb_lat = min(lats)
    ub_lat = max(lats)
    lb_lon = min(lons)
    ub_lon = max(lons)
    lb = s2.LatLng.from_degrees(lb_lat, lb_lon)
    ub = s2.LatLng.from_degrees(ub_lat, ub_lon)
    cids = r.get_covering(s2.LatLngRect.from_point_pair(lb,ub))
    s2_token_list = []
    for cellid in cids:
        s2_token_list.append(cellid.to_token())

    return s2_token_list, cids

def getBoundaryCoverageNew(cids, poly, maxResColName):
    s2_index__L19_list = []
    p_gdf = gpd.GeoDataFrame()
    idx=0
    for cellid in cids:
        c = s2.Cell(cellid)
        level = cellid.level()
        vertices = []
        for i in range(0,4):
            vertex = c.get_vertex(i)
            latlng = s2.LatLng.from_point(vertex)
            vertices.append((latlng.lng().degrees, latlng.lat().degrees))
        geo = shapely.geometry.Polygon(vertices)
        if poly.intersects(geo):
            s2_index__L19_list.append(cellid.to_token())
            p_gdf.loc[idx,maxResColName] = cellid.to_token()
            p_gdf.loc[idx,'geometry']=geo
        idx+=1
        
    return p_gdf

def getnewNDVIFromParquetForDataSets(datasets, NDays, s2_index__L10_list, s2_index__L15_list, s2_index__L19_list):
    start_time = time.time()
    subset_dirs = []
    endDt = datetime.now() #SET THE MARKET AS NOW

    startDt = endDt - timedelta(days=NDays) 
    dtStrs_list = [datetime.strftime(x, '%Y-%m-%d') for x in pd.date_range(startDt, periods=NDays).tolist()]
    
    p_df = pd.DataFrame()
    for dtset in datasets:

        p_subset_df = dtset.to_table(
            columns=['s2_index__L19', 'NDVI', 'UTC_DATETIME','YY_MM_DD'],
            filter=(
                ds.field('s2_index__L10').isin(s2_index__L10_list) & 
                ds.field('s2_index__L19').isin(s2_index__L19_list) &
                ds.field('YY_MM_DD').isin(dtStrs_list)
            )).to_pandas()

        p_df = p_df.append(p_subset_df, ignore_index=True)
    
    
    if len(p_df)==0:
        print('Empty Dataset')
    else:
        p_df.reset_index(drop=True, inplace=True)
        
    end_time = time.time()
    time_elapsed = (end_time - start_time)

    return p_df, time_elapsed


def getnewNDVIFromParquetForDataSets_maxRes(maxRes, datasets, NDays, s2_index__L10_list, s2_index__L15_list, s2_index__maxRes_list):
    start_time = time.time()
    subset_dirs = []
    endDt = datetime.now() #SET THE MARKET AS NOW
    
    maxResColName = 's2_index__L'+str(maxRes)

    startDt = endDt - timedelta(days=NDays) 
    dtStrs_list = [datetime.strftime(x, '%Y-%m-%d') for x in pd.date_range(startDt, periods=NDays).tolist()]
    
    p_df = pd.DataFrame()
    for dtset in datasets:

        p_subset_df = dtset.to_table(
            columns=[maxResColName, 'NDVI', 'UTC_DATETIME','YY_MM_DD'],
            filter=(
                ds.field('s2_index__L10').isin(s2_index__L10_list) & 
                ds.field(maxResColName).isin(s2_index__maxRes_list) &
                ds.field('YY_MM_DD').isin(dtStrs_list)
            )).to_pandas()

        p_df = p_df.append(p_subset_df, ignore_index=True)
    
    
    if len(p_df)==0:
        print('Empty Dataset')
    else:
        p_df.reset_index(drop=True, inplace=True)
        
    end_time = time.time()
    time_elapsed = (end_time - start_time)

    return p_df, time_elapsed


def getNDVIAndStatsForDtStr_maxRes(maxRes, p_df, p_gdf, dtStr):
    # For the first one
    maxResColName = 's2_index__L'+str(maxRes)
    keepCols = ['NDVI', maxResColName, 'YY_MM_DD']

    dtStr_df = p_df[p_df.YY_MM_DD==dtStr][keepCols]
    dtStr_df.reset_index(drop=True, inplace=True)

    #group by the index and take average
    n_df = dtStr_df.groupby(maxResColName).mean()
    n_df = n_df.dropna(subset=['NDVI'])
    n_df.reset_index(inplace=True)

    res_gdf = p_gdf.merge(n_df, on=maxResColName, how = 'inner')
    r_gdf = res_gdf.dropna(subset=['NDVI'])
    r_gdf.reset_index(drop=True, inplace=True)
    
    ndvi_mean = np.nanmean(r_gdf['NDVI'])
    ndvi_p10 = np.nanpercentile(r_gdf['NDVI'], 10)
    ndvi_p25 = np.nanpercentile(r_gdf['NDVI'], 25)
    ndvi_p50 = np.nanpercentile(r_gdf['NDVI'], 50)
    ndvi_p75 = np.nanpercentile(r_gdf['NDVI'], 75)
    ndvi_p90 = np.nanpercentile(r_gdf['NDVI'], 90)
    ndvi_min = np.min(r_gdf['NDVI'])
    ndvi_max = np.max(r_gdf['NDVI'])
    
    r_stats = [ndvi_mean,ndvi_min,ndvi_p10,ndvi_p25,ndvi_p50,ndvi_p75,ndvi_p90,ndvi_max]
    
    return r_gdf, r_stats
    
    
def getInterpolatedValueFromEdgeNeightbors(tkn):
    est_val = np.nan
   
    try:
        neighboring_cells = s2.CellId.from_token(tkn).get_edge_neighbors()
        arr_vals = [res_gdf[res_gdf.s2_index__L20==x.to_token()].NDVI.values.tolist() for x in neighboring_cells]
        est_val = np.nanmean([x[0] for x in arr_vals if len(x)>0])
    except:
        est_val = np.nan
       
    return est_val



app = Flask(__name__)

pa_filePath_s2_L19 = '/mnt/md1/SENTINEL/PARQUET_S2/'
pa_filePath_s2_L20 = '/mnt/md1/SENTINEL/PARQUET_S2_L20/'

resLevel = 19
pa_filePath = eval('pa_filePath_s2_L'+str(resLevel))
s2_max_resolution = resLevel
maxResColName = 's2_index__L'+str(resLevel)

################  ROUTES
@app.route('/')
#@app.route('/index')
#@app.route('/home')
#def index():
#    return render_template('index.html')
#
#@app.route('/')
def hello_world():
    return 'Hello World!'

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

