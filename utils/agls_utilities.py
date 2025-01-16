#!/usr/bin/env python
# coding: utf-8
from utils.imports import *
from utils.settings import *
import pytz 


def getArea(geom):
    geom_area = ops.transform(
    partial(
        pyproj.transform,
        pyproj.Proj(init='EPSG:4326'),
        pyproj.Proj(
            proj='aea',
            lat_1=geom.bounds[1],
            lat_2=geom.bounds[3])),
    geom)

    # Return the area in km^2
    return geom_area.area / 1000000


def fixFootprint(curFootprintStr):
    #if this is a list, then simply return a string of the list
    if type(curFootprintStr)==list:
        newFootprintStr = str(curFootprintStr)
    else:
        #First see if this has a "," in it, tokensize
        if not '[' in curFootprintStr:
            if ',' in curFootprintStr:
                tok = curFootprintStr.split(',')
                newFootprintStr = str(tok)
            else:
                newFootprintStr = '['+curFootprintStr+']'
        else:
            newFootprintStr=curFootprintStr        
            
    return newFootprintStr

def fixCustomerid(customerid):
    #First see if this has a "," in it, tokensize
    if type(customerid)==str:
        newCustomerid=int(customerid)
    else:
        newCustomerid=customerid
        
    return newCustomerid

def getSatTilesInfo(gdf):
    #find the Sentinel and LandSat tiles for each field
    sentinelOverlapTile_df = pd.DataFrame(columns=['TILE','geometry'])
    landsatOverlapTile_df = pd.DataFrame(columns=['TILE','geometry'])

    sTiles_df = pd.DataFrame(columns=['TILE','geometry'])
    lTiles_df = pd.DataFrame(columns=['TILE','geometry'])

    for i,field in gdf.iterrows():

        #print([field.customerid, field.fieldassetid])
        fieldSentinelTileList = []
        fieldLandsatTileList = []

        fieldPoly = field.geometry

        #Sentinel
        bool_list = gdf_S['geometry'].apply(lambda g: g.intersects(fieldPoly))
        gdf_S_subset = gdf_S.loc[bool_list]
        gdf_S_subset.reset_index(drop=True, inplace=True)
        tileList = pd.unique(gdf_S_subset.TILE).tolist()
        fieldSentinelTileList.append(tileList)
        sTiles_df = sTiles_df.append(gdf_S_subset, ignore_index=True)

        #Landsat
        bool_list = gdf_L['geometry'].apply(lambda g: g.intersects(fieldPoly))
        gdf_L_subset = gdf_L.loc[bool_list]
        gdf_L_subset.reset_index(drop=True, inplace=True)
        tileList = pd.unique(gdf_L_subset.TILE).tolist()
        fieldLandsatTileList.append(tileList)
        lTiles_df = lTiles_df.append(gdf_L_subset, ignore_index=True)


        #Transform the GeoDataFrame
        slist = fieldSentinelTileList[0]
        llist = fieldLandsatTileList[0]
        
        sfoot = fixFootprint(slist)
        #'['+ makeStrFromList(slist) + ']'
        lfoot = fixFootprint(llist)
        #'['+makeStrFromList(llist)+ ']'
        
        gdf.loc[i,'SentinelFootprint'] = sfoot
        gdf.loc[i,'LandsatFootprint'] = lfoot

        print([field.customerid, field.fieldassetid, sfoot, lfoot])
        
    gdf_updated = gdf.copy()
    
    return gdf_updated

def getFeatures(gdf):
    """Function to parse features from GeoDataFrame in such a manner that rasterio wants them"""
    return [json.loads(gdf.to_json())['features']]


# In[3]:


def getCCGeoDataFrame(gmlPath):

    gmlFileStr = 'MSK_CLOUDS_B00.gml'
    gmlFile = gmlPath + gmlFileStr
    
    f = open(gmlFile, 'r')
    lines = f.readlines()
    f.close()
    
    start=0
    name=0
    points = []
    data = []
    for l in lines:
        line = l.strip()
    
        if (line[0:13] =='<gml:Envelope'):
            l = line.strip().split(':')
            proj = l[-1][0:-2]
    
        if (line[0:12]=='<gml:Polygon'):
            l = line.strip()[13:-1]
            tok = l.split('=')
            name = tok[1].replace('\"','')

        if (line[0:12]=='<gml:posList'):
            start=1
        if (line=='</gml:LinearRing>'):
            start=0
        
        if start==1:
            #print(name)
            type = name.split('.')[0]
            s = line[30:-14]
            c = s.split(' ')
        
            r = int(len(c)/2)
            points=[]
            for i in range(r):
                points.append(Point([int(c[i*2]),int(c[i*2+1])]))
            #male polygon and Add to the geopandas    
            coords = sum(map(list, (p.coords for p in points)), [])
            poly = Polygon(coords)
            data.append(dict({'Cloud_Name': name, 'Mask_Type': type, 'geometry': poly}))
    
    df = pd.DataFrame(data)
    if (len(df)>0):
        e = 'epsg:'+proj
        crs = {'init': e}
        cloudsGDF = gpd.GeoDataFrame(df, crs=crs, geometry=df.geometry)
        return cloudsGDF
    else:
        nullDataFrame = gpd.GeoDataFrame()
        return nullDataFrame
    


# In[4]:


def getMostRecentFile(sentinel_tiff_files_array):
    dts = []
    #print(sentinel_tiff_files_array)
    for file in sentinel_tiff_files_array:
        d_str = file.split('_')[2][0:8]
        dt = parse(d_str)
        dts.append(dt)
    idx = [dts.index(x) for x in sorted(dts)]
    
    return sentinel_tiff_files_array[idx[-1]]


# In[5]:


def getAreaOfROI(roi_polygon):
    #Returns area in KM^2 and the fraction of the 100x100km Sentinel

    geom = shapely.geometry.Polygon(roi_polygon)
    #geom = Polygon([(0, 0), (0, 10), (10, 10), (10, 0), (0, 0)])
    geom_area = shapely.ops.transform(
        partial(
            pyproj.transform,
            pyproj.Proj(init='EPSG:4326'),
            pyproj.Proj(
                proj='aea',
                lat_1=geom.bounds[1],
                lat_2=geom.bounds[3])),
        geom)

    # return the area in Km^2
    a = (geom_area.area / 1000000)
    f = a / 12056.04 #109.8x 109.8km = 12,056.04 km ^2
    return [a,f]


# In[6]:


def getFootprint(XML_file):
    f = open(XML_file, "r")
    lines = f.readlines()
    s=0
    n=0
    crd = []

    for line in lines:
        n=n+1
        if (line.strip().startswith('<Global_Footprint>')):
            s=1
            n=0
        if (line.strip().startswith('</Global_Footprint>')):
            s=0
        if s==1 & n==1:
            l = line.strip()
            t = l.split('<EXT_POS_LIST>')
            c = t[1][:-15].strip()
            crds = c.split(' ')
        
            for i in range(int(len(crds)/2)):
                lat = float(crds[i*2])
                lon = float(crds[i*2+1])
                ll = [lon,lat]
                crd.append(ll)
    return crd


def getBounds(tileStr):
    tile_df = S_gdf[S_gdf.Name==tileStr].geometry.bounds
    tile_df.reset_index(drop=True,inplace=True)
    bounds = [
        tile_df.loc[0,'minx'],
        tile_df.loc[0,'miny'],
        tile_df.loc[0,'maxx'],
        tile_df.loc[0,'maxy']
    ]
    
    return bounds

def getWMSData(from_dt,to_dt,bounding_box):
    INSTANCE_ID = '5c53cc29-8f4a-43ed-b97b-cb27be196176'
    LAYER_NAME = '3_NDVI'
    bands_script = 'return [B01,B02,B04,B05,B08,B8A,B09,B10,B11,B12]'
    wms_bands_request = WmsRequest(layer=LAYER_NAME,
                               custom_url_params={
                                   CustomUrlParam.EVALSCRIPT: bands_script,
                                   CustomUrlParam.ATMFILTER: 'NONE'
                               },
                               bbox=bounding_box, 
                               time=(from_dt, to_dt), 
                               width=600, height=None,
                               image_format=MimeType.TIFF_d32f,
                               instance_id=INSTANCE_ID)
    wms_bands = wms_bands_request.get_data()
    return wms_bands

def plot_cloud_mask(mask, figsize=(15, 15), fig=None):
    """
    Utility function for plotting a binary cloud mask.
    """
    if fig == None:
        plt.figure(figsize=figsize)
    plt.imshow(mask, cmap=plt.cm.gray)
    
def getnir20(imagePath,imgFile):
    upscale_factor = 0.5
    with rasterio.open(imagePath+imgFile, driver='JP2OpenJPEG') as dataset:
        # resample data to target shape
        img = dataset.read(
            out_shape=(
                dataset.count,
                int(dataset.height * upscale_factor),
                int(dataset.width * upscale_factor)
            ),
            resampling=Resampling.bilinear
        )

        # scale image transform
        transform = dataset.transform * dataset.transform.scale(
            (dataset.width / img.shape[-1]),
            (dataset.height / img.shape[-2])
        )
        img = img.squeeze()

    img = img.astype('float')   
    nir_20m = img.copy()
    
    return nir_20m

def getCloudMask(imagePath,imgFile):
    upscale_factor = 2
    with rasterio.open(imagePath+imgFile, driver='JP2OpenJPEG') as dataset:
        # resample data to target shape
        img = dataset.read(
            out_shape=(
                dataset.count,
                int(dataset.height * upscale_factor),
                int(dataset.width * upscale_factor)
            ),
            resampling=Resampling.bilinear
        )

        # scale image transform
        transform = dataset.transform * dataset.transform.scale(
            (dataset.width / img.shape[-1]),
            (dataset.height / img.shape[-2])
        )
        img = img.squeeze()

    img = img.astype('float')   
    img[img>0]=np.nan
    img[img==0]=1.0    # Set all probabilitiy > 0 as 1, "MOST CONSERVATIVE"
    
    mask = img.copy()
    
    return mask

def getCloudMaskArchives(imagePath,imgFile):
    upscale_factor = 1
    with rasterio.open(imagePath+imgFile) as dataset:
        # resample data to target shape
        img = dataset.read(
            out_shape=(
                dataset.count,
                int(dataset.height * upscale_factor),
                int(dataset.width * upscale_factor)
            ),
            resampling=Resampling.bilinear
        )

        # scale image transform
        transform = dataset.transform * dataset.transform.scale(
            (dataset.width / img.shape[-1]),
            (dataset.height / img.shape[-2])
        )
        img = img.squeeze()

    img = img.astype('float')   
    img[img>0]=np.nan
    img[img==0]=1.0    # Set all probabilitiy > 0 as 1, "MOST CONSERVATIVE"
    
    mask = img.copy()
    
    return mask


def getCloudMask_20m(imagePath,imgFile):
    upscale_factor = 1 #no need to upsample since it is a 20m pixel
    with rasterio.open(imagePath+imgFile, driver='JP2OpenJPEG') as dataset:
        # resample data to target shape
        img = dataset.read(
            out_shape=(
                dataset.count,
                int(dataset.height * upscale_factor),
                int(dataset.width * upscale_factor)
            ),
            resampling=Resampling.bilinear
        )

        # scale image transform
        transform = dataset.transform * dataset.transform.scale(
            (dataset.width / img.shape[-1]),
            (dataset.height / img.shape[-2])
        )
        img = img.squeeze()

    img = img.astype('float')   
    img[img>0]=np.nan
    img[img==0]=1.0    # Set all probabilitiy > 0 as 1, "MOST CONSERVATIVE"
    
    mask = img.copy()
    
    return mask

def getTileWkt(tile):
    df = S_gdf[S_gdf.TILE==tile].geometry
    df.reset_index(drop=True,inplace=True)
    geom = df.geometry[0][0]
    wkt = geom.to_wkt()
    coords = list(geom.exterior.coords)
    
    return wkt


def getExtent(polygon):
    lats = list(polygon.exterior.xy[1])
    lons = list(polygon.exterior.xy[0])
    ext = [min(lons),max(lons),min(lats),max(lats)]
    
    return ext

def getImgDateStr(f):
    dtStr = os.path.split(f)[1].split('_')[2][0:8]
    imgDateStr = dtStr[0:4]+'-'+dtStr[4:6]+'-'+dtStr[6:8]
    
    return imgDateStr

def getImgDateStr2(f):
    #else if Landsat:
    #if f.split('/')[4]=='LANDSAT8':
    if f.split('/')[3]=='LANDSAT8':
        dtStr = os.path.split(f)[1].split('_')[3][0:8]
        imgDateStr = dtStr[0:4]+'-'+dtStr[4:6]+'-'+dtStr[6:8]
    else:
    #elif Sentinel:
        dtStr = os.path.split(f)[1].split('_')[2][0:8]
        imgDateStr = dtStr[0:4]+'-'+dtStr[4:6]+'-'+dtStr[6:8]
    
    try:
        datetime.strptime(imgDateStr, "%Y-%m-%d")
    except:
        print('Value Error: '+f+' '+dtStr+' '+imgDateStr)
        imgDateStr = 'NaT'
        
    return imgDateStr

#Generate the list of NDVI files - BOTH SENTINEL AND LANDSAT - for last Ndays
def getNDVIFilesWithDate(Ndays, S_tileList, L_tileList, dtStr):
    N = Ndays
    dt = datetime.strptime(dtStr,'%Y-%m-%d')
    endDay = dt
    startDay = dt - timedelta(days=N)
    
    #get the list of Tiles
    STiles = LTiles = []
    STiles7 = [os.path.abspath(os.path.join(sentinelProcFolder, p)) for p in os.listdir(sentinelProcFolder)]
    STiles = list(set(STiles7))
    
    LTiles = [x for x in os.listdir(landsatProcFolder)]
    
    SfileList = []
    for f in STiles:
        folder = f + '/'
        tiffFiles = [x for x in os.listdir(folder) if x.endswith('.tiff')]
        
        for fl in tiffFiles:
            fldtstr = fl.split('_')[2][0:8]
            fldt = datetime.strptime(fldtstr, '%Y%m%d')
            #print([dtStr, fldt, startDay, endDay])
            if ((fldt >= startDay) & (fldt <= endDay)):
                f = folder + fl
                SfileList.append(f)

    LfileList = []
    pattern = '*_2020*_T1'
    tifFolders = fnmatch.filter(os.listdir(landsatProcFolder), pattern)
    for fldr in tifFolders:
        fl = landsatProcFolder + fldr + '/' + fldr+'_NDVI_WGS84.tiff'
        fldtstr = fl.split('_')[3][0:8]
        fldt = datetime.strptime(fldtstr, '%Y%m%d')
        if ((fldt >= startDay) & (fldt <= endDay)):
            LfileList.append(fl)          
    

    #make the df
    ind=0
    ndvi_df = pd.DataFrame(columns = ['DATE','TILE','FILE'])
    for f in SfileList:
        s = f.split('/')

        if s[1]=='mnt':
            tile=s[5]
        else:
            tile=s[6]

        dtstr = s[-1].split('_')[2]
        ts = datetime.strptime(dtstr,'%Y%m%dT%H%M%S')
        #print([ts, tile, f])
        ndvi_df.loc[ind,'DATE']=ts
        ndvi_df.loc[ind,'TILE']=tile
        ndvi_df.loc[ind,'FILE']=f
        ind=ind+1

    #print(ndvi_df)


    for fl in LfileList:
        s = fl.split('/')
        tile = s[6].split('_')[2]
        dtstr = s[6].split('_')[3]
        ts = datetime.strptime(dtstr,'%Y%m%d')
        ndvi_df.loc[ind,'DATE']=ts
        ndvi_df.loc[ind,'TILE']=tile
        ndvi_df.loc[ind,'FILE']=fl
        ind=ind+1

  
    ndvi_df['DATE'] =pd.to_datetime(ndvi_df.DATE)
    

    #now filter this with the tileList provided for the field
    tileList = S_tileList + L_tileList
    
    ndf = ndvi_df[ndvi_df.TILE.isin(tileList)]
    df =  ndf.copy()
    #sort values
    df.sort_values(by=['DATE'],inplace=True, ascending=False)
    df.reset_index(drop = True)
    
    return df


#Generate the list of NII files - BOTH SENTINEL AND LANDSAT - for last 30 days
def getNIIFilesWithDate(Ndays, S_tileList, L_tileList, dtStr):
    N = Ndays
    dt = datetime.strptime(dtStr,'%Y-%m-%d')
    #print(dt)
    endDay = dt
    startDay = dt - timedelta(days=N)
    
    #get the list of Tiles
    STiles = LTiles = []
    STiles = [os.path.abspath(os.path.join(sentinel_EXTENSIONS_ProcFolder, p)) for p in os.listdir(sentinel_EXTENSIONS_ProcFolder)]
    STiles = list(set(STiles))

    #print(STiles)
    
    SfileList = []
    for f in STiles:
        folder = f + '/'
        tiffFiles = [x for x in os.listdir(folder) if x.endswith('.tiff')]
        #print(tiffFiles)

        for fl in tiffFiles:
            fldtstr = fl.split('_')[2][0:8]
            fldt = datetime.strptime(fldtstr, '%Y%m%d')
            if ((fldt >= startDay) & (fldt <= endDay)):
                f = folder + fl
                SfileList.append(f)   

    #print(SfileList)

    #make the df
    ind=0
    nii_df = pd.DataFrame(columns = ['DATE','TILE','FILE'])
    for f in SfileList:
        s = f.split('/')

        if s[1]=='mnt':
            tile=s[5]
        else:
            tile=s[6]

        dtstr = s[-1].split('_')[2]
        ts = datetime.strptime(dtstr,'%Y%m%dT%H%M%S')
        #print([ts, tile, f])
        nii_df.loc[ind,'DATE']=ts
        nii_df.loc[ind,'TILE']=tile
        nii_df.loc[ind,'FILE']=f
        ind=ind+1
  
    nii_df['DATE'] =pd.to_datetime(nii_df.DATE)

    #now filter this with the tileList provided for the field
    tileList = S_tileList
    
    ndf = nii_df[nii_df.TILE.isin(tileList)]
    df =  ndf.copy()
    #sort values
    df.sort_values(by=['DATE'],inplace=True, ascending=False)
    df.reset_index(drop = True)
    
    return df


#Generate the list of NDVI files - BOTH SENTINEL AND LANDSAT - for last 30 days
def getMultibandFilesWithDate(Ndays, S_tileList, fileType, dtStr):
    #fileType = '20m' or '10m'
    N = Ndays
    dt = datetime.strptime(dtStr,'%Y-%m-%d')
    endDay = dt
    startDay = dt - timedelta(days=N)
    
    #get the list of Tiles
    STiles = LTiles = []
    STiles7 = [os.path.abspath(os.path.join(sentinel_MULTIBAND_ProcFolder, p)) for p in os.listdir(sentinel_MULTIBAND_ProcFolder)]
    STiles = list(set(STiles7))
        
    SfileList = []
    for f in STiles:
        folder = f + '/'
        tiffFiles = [x for x in os.listdir(folder) if ((x.endswith('.tiff')) & (fileType in x))]
        
        for fl in tiffFiles:
            fldtstr = fl.split('_')[2][0:8]
            fldt = datetime.strptime(fldtstr, '%Y%m%d')
            #print([dtStr, fldt, startDay, endDay])
            if ((fldt >= startDay) & (fldt <= endDay)):
                f = folder + fl
                SfileList.append(f)

    #make the df
    ind=0
    mb_df = pd.DataFrame(columns = ['DATE','TILE','FILE'])
    for f in SfileList:
        s = f.split('/')
        tile=s[5]
        dtstr = s[-1].split('_')[2]
        ts = datetime.strptime(dtstr,'%Y%m%dT%H%M%S')
        #print([ts, tile, f])
        mb_df.loc[ind,'DATE']=ts
        mb_df.loc[ind,'TILE']=tile
        mb_df.loc[ind,'FILE']=f
        ind=ind+1

    mb_df['DATE'] =pd.to_datetime(mb_df.DATE)
    
    #now filter this with the tileList provided for the field
    ndf = mb_df[mb_df.TILE.isin(S_tileList)]
    df =  ndf.copy()
    #sort values
    df.sort_values(by=['DATE'],inplace=True, ascending=False)
    df.reset_index(drop = True)
    
    return df


def getNDVIValueWithStatsAndImg2(df,polygon):
    #Returns the image (2-d numpy array) and the extent
    
    for i,row in df.iterrows():
        
        tifFile = row.FILE
        f = tifFile
        ndvi_imgArr = []
        ndvi_extList = []
        imgDtStr = []
    
        mval_stats = pd.DataFrame(columns=['Mean','p90','p25','p75','Median','StdDev','Extent','ImgDateStr'])

        roi_polygon = list(polygon.exterior.coords)

        src = rasterio.open(f,'r', count=1) 
            
        roi_polygon_src_coords = warp.transform_geom({'init': 'epsg:4326'},
                                                     src.crs,
                                                     {'type': 'Polygon','coordinates': [roi_polygon]})
        imgDtStr = getImgDateStr2(f)
        ext = getExtent(polygon)
        extJSON = json.dumps(ext)
        
        #try:
        out_image,out_transform = mask.mask(src,[roi_polygon_src_coords],crop=True, all_touched=True)
        imm = out_image[0,::]

        if (np.nansum(imm) > 0):
            #Do the cutting and analysis
            ndvi = imm / 255
            ndvi = ndvi.astype('float')

            ndvi_imgArr = ndvi
            #replace illegitimate values

            ndvi[ndvi <= 0] = np.nan
            ndvi[ndvi > 1] = np.nan
            data = ndvi.flatten()
            fdata = data[~np.isnan(data)]
            #print('Field Img Carve Date: '+ imgDtStr)


            mval_stats.loc[0,'Mean']=np.mean(fdata)
            mval_stats.loc[0,'p90']=np.percentile(fdata,90)
            mval_stats.loc[0,'p25']=np.percentile(fdata,25)
            mval_stats.loc[0,'p75']=np.percentile(fdata,75)
            mval_stats.loc[0,'Median']=np.median(fdata)
            mval_stats.loc[0,'StdDev']=np.std(fdata)
            mval_stats.loc[0,'Extent']=extJSON
            mval_stats.loc[0,'ImgDateStr']=imgDtStr

            break
        else:
            mval_stats = []
            continue
        #except:
        #    #print('except...')
        #    continue
        
    return [mval_stats, ndvi_imgArr, extJSON , imgDtStr]


def getMBValueWithStatsAndImg(df,imgType,fileType,polygon):
    #fileType = '10m', '20m'
    #imgType = 'GNDVI', 'VARI', LAI, 'NDII', etc for either 10m or 20m

    #Returns the image (2-d numpy array) and the extent
    
    for i,row in df.iterrows():
        
        tifFile = row.FILE
        f = tifFile
        band_imgArr = []
        band_extList = []
        imgDtStr = []
    
        mval_stats = pd.DataFrame(columns=['Mean','p90','p25','p75','Median','StdDev','Extent','ImgDateStr'])

        roi_polygon = list(polygon.exterior.coords)

        try:
            if (fileType=='10m'):
                src = rasterio.open(f,'r', count=4)
            else:
                src = rasterio.open(f,'r', count=7)
        
                
            roi_polygon_src_coords = warp.transform_geom(
                {'init': 'epsg:4326'},
                src.crs,
                {'type': 'Polygon','coordinates': [roi_polygon]}
                )

            imgDtStr = getImgDateStr2(f)
            ext = getExtent(polygon)
            extJSON = json.dumps(ext)

            out_image,out_transform = mask.mask(src,[roi_polygon_src_coords],crop=True, all_touched=True)

            if (fileType=='10m'):
                #Define the 10m bands
                """
                gndvi_raw = out_image[0,::]
                vari_raw = out_image[1,::]
                savi_raw = out_image[2,::]
                evi_raw = out_image[3,::]
                """
                if (imgType=='GNDVI'):
                    imm = out_image[0,::]
                elif (imgType=='VARI'):
                    imm = out_image[1,::]
                elif (imgType=='SAVI'):
                    imm = out_image[2,::]
                elif (imgType=='EVI'):
                    imm = out_image[3,::]
                else: #(catch all)
                    imm = out_image[0,::]

            else:
                #define the 20m bands
                """
                ndre_raw = out_image[0,::]
                srre_raw = out_image[1,::]
                ndii_raw = out_image[2,::]
                lci_raw = out_image[3,::]
                lai_raw = out_image[4,::]
                cnu_raw = out_image[5,::]
                ch_raw = out_image[6,::]
                """
                if (imgType=='NDRE'):
                    imm = out_image[0,::]
                elif (imgType=='SRRE'):
                    imm = out_image[1,::]
                elif (imgType=='NDII'):
                    imm = out_image[2,::]
                elif (imgType=='LCI'):
                    imm = out_image[3,::]
                elif (imgType=='LAI'):
                    imm = out_image[4,::]
                elif (imgType=='CNU'):
                    imm = out_image[5,::]
                elif (imgType=='CH'):
                    imm = out_image[6,::]
                else: #(catch all)
                    imm = out_image[0,::]

            polyac = area(mapping(polygon)) * 0.000247
            maskac = np.count_nonzero(imm)*100*0.000247
            fac = maskac / polyac
            #print(fac)

            #if ((len(np.unique(imm[~np.isnan(imm)])) > 0) & (fac >= 0.75)):
            if len(np.unique(imm[~np.isnan(imm)])) > 0:
                #Do the cutting and analysis
                
                band = imm / 255
                band = band.astype('float')
                band_imgArr = band
                #replace illegitimate values
                #band[band <= 0] = np.nan
                #band[band > 1] = np.nan
                

                data = band.flatten()
                fdata = data[~np.isnan(data)]
                #print('Field Img Carve Date: '+ imgDtStr)


                mval_stats.loc[0,'Mean']=np.mean(fdata)
                mval_stats.loc[0,'p90']=np.percentile(fdata,90)
                mval_stats.loc[0,'p25']=np.percentile(fdata,25)
                mval_stats.loc[0,'p75']=np.percentile(fdata,75)
                mval_stats.loc[0,'Median']=np.median(fdata)
                mval_stats.loc[0,'StdDev']=np.std(fdata)
                mval_stats.loc[0,'Extent']=extJSON
                mval_stats.loc[0,'ImgDateStr']=imgDtStr

                break
            else:
                mval_stats = []
                continue
        except:
            continue
        
    return [mval_stats, band_imgArr, extJSON , imgDtStr]



def getNDVIValueWithStatsAndImg3(df,polygon):
    #Returns the image (2-d numpy array) and the extent
    mval_stats = pd.DataFrame()
    ndvi_imgArr = []
    ndvi_extList = []
    imgDtStr = ''
    extJSON = ''

    for i,row in df.iterrows():
        
        tifFile = row.FILE
        f = tifFile
        ndvi_imgArr = []
        ndvi_extList = []
        imgDtStr = ''
    
        mval_stats = pd.DataFrame(columns=['Mean','p90','p25','p75','Median','StdDev','Extent','ImgDateStr'])

        roi_polygon = list(polygon.exterior.coords)

        try:
            if (tifFile.split('/')[1]=='mnt'):
                src = rasterio.open(f,'r', count=1)
            else:
                src = rasterio.open(f,'r', count=7)
                
            roi_polygon_src_coords = warp.transform_geom(
                {'init': 'epsg:4326'},
                src.crs,
                {'type': 'Polygon','coordinates': [roi_polygon]}
                )

            imgDtStr = getImgDateStr2(f)
            ext = getExtent(polygon)
            extJSON = json.dumps(ext)

            out_image,out_transform = mask.mask(src,[roi_polygon_src_coords],crop=True, all_touched=True)

            if (tifFile.split('/')[1]=='mnt'):
                imm = out_image[0,::]
            else:
                imm = out_image[4,::]

            polyac = area(mapping(polygon)) * 0.000247
            maskac = np.count_nonzero(imm)*100*0.000247
            fac = maskac / polyac

            if ((np.nansum(imm) > 0) & (fac >= 0.75)):
                #Do the cutting and analysis
                ndvi = imm / 255
                ndvi = ndvi.astype('float')

                ndvi_imgArr = ndvi
                #replace illegitimate values

                ndvi[ndvi <= 0] = np.nan
                ndvi[ndvi > 1] = np.nan
                data = ndvi.flatten()
                fdata = data[~np.isnan(data)]
                #print('Field Img Carve Date: '+ imgDtStr)


                mval_stats.loc[0,'Mean']=np.mean(fdata)
                mval_stats.loc[0,'p90']=np.percentile(fdata,90)
                mval_stats.loc[0,'p25']=np.percentile(fdata,25)
                mval_stats.loc[0,'p75']=np.percentile(fdata,75)
                mval_stats.loc[0,'Median']=np.median(fdata)
                mval_stats.loc[0,'StdDev']=np.std(fdata)
                mval_stats.loc[0,'Extent']=extJSON
                mval_stats.loc[0,'ImgDateStr']=imgDtStr

                break
            else:
                mval_stats = []
                continue
        except:
            continue
        
    return [mval_stats, ndvi_imgArr, extJSON , imgDtStr]

def getNDVIValueWithStatsAndImg4(df,polygon):
    #Returns the image (2-d numpy array) and the extent
    #use epsg value of epsg:3857
    mval_stats = pd.DataFrame()
    ndvi_imgArr = []
    extJSON = ''
    imgDtStr = ''


    for i,row in df.iterrows():
        
        tifFile = row.FILE
        f = tifFile
        ndvi_imgArr = []
        ndvi_extList = []
        imgDtStr = []
    
        mval_stats = pd.DataFrame(columns=['Mean','p90','p25','p75','Median','StdDev','Extent','ImgDateStr'])

        roi_polygon = list(polygon.exterior.coords)

        #try:
        if (tifFile.split('/')[1]=='mnt'):
            src = rasterio.open(f,'r', count=1)
        else:
            src = rasterio.open(f,'r', count=7)
            

        #print(src.crs)

        roi_polygon_src_coords = warp.transform_geom(
            {'init': 'epsg:4326'},
            src.crs,
            {'type': 'Polygon','coordinates': [roi_polygon]}
            )


        imgDtStr = getImgDateStr2(f)
        ext = getExtent(polygon)
        extJSON = json.dumps(ext)

        out_image,out_transform = mask.mask(src,[roi_polygon_src_coords], crop=True, pad=True)

        if (tifFile.split('/')[1]=='mnt'):
            imm = out_image[0,::]
        else:
            imm = out_image[4,::]

        polyac = area(mapping(polygon)) * 0.000247
        maskac = np.count_nonzero(imm)*100*0.000247
        fac = maskac / polyac

        #print(fac)

        if ((np.nansum(imm) > 0) & (fac >= 0.75)):
            #Do the cutting and analysis
            ndvi = imm / 255
            ndvi = ndvi.astype('float')

            ndvi_imgArr = ndvi
            #replace illegitimate values

            ndvi[ndvi <= 0] = np.nan
            ndvi[ndvi > 1] = np.nan
            data = ndvi.flatten()
            fdata = data[~np.isnan(data)]
            #print('Field Img Carve Date: '+ imgDtStr)


            mval_stats.loc[0,'Mean']=np.mean(fdata)
            mval_stats.loc[0,'p90']=np.percentile(fdata,90)
            mval_stats.loc[0,'p25']=np.percentile(fdata,25)
            mval_stats.loc[0,'p75']=np.percentile(fdata,75)
            mval_stats.loc[0,'Median']=np.median(fdata)
            mval_stats.loc[0,'StdDev']=np.std(fdata)
            mval_stats.loc[0,'Extent']=extJSON
            mval_stats.loc[0,'ImgDateStr']=imgDtStr

            break
        else:
            mval_stats = []
            continue
        #except:
        #	continue
        
    return [mval_stats, ndvi_imgArr, extJSON , imgDtStr]

def getNDVIValueWithStatsAndImg5(df,polygon):
    #Returns the image (2-d numpy array) and the extent
    #use epsg value of epsg:3857
    mval_stats = pd.DataFrame()
    ndvi_imgArr = []
    extJSON = ''
    imgDtStr = ''


    for i,row in df.iterrows():
        
        tifFile = row.FILE
        f = tifFile
        ndvi_imgArr = []
        ndvi_extList = []
        imgDtStr = []
    
        mval_stats = pd.DataFrame(columns=['Mean','p90','p25','p75','Median','StdDev','Extent','ImgDateStr'])

        roi_polygon = list(polygon.exterior.coords)

        #try:
        if (tifFile.split('/')[1]=='mnt'):
            src = rasterio.open(f,'r', count=1)
        else:
            src = rasterio.open(f,'r', count=7)
            

        #print(src.crs)

        roi_polygon_src_coords = warp.transform_geom(
            {'init': 'epsg:4326'},
            src.crs,
            {'type': 'Polygon','coordinates': [roi_polygon]}
            )


        imgDtStr = getImgDateStr2(f)
        ext = getExtent(polygon)
        extJSON = json.dumps(ext)

        out_image,out_transform = mask.mask(src,[roi_polygon_src_coords], crop=True, pad=True)

        if (tifFile.split('/')[1]=='mnt'):
            imm = out_image[0,::]
        else:
            imm = out_image[4,::]

        polyac = area(mapping(polygon)) * 0.000247
        maskac = np.count_nonzero(imm)*100*0.000247
        fac = maskac / polyac

        #print(fac)

        if ((np.nansum(imm) > 0) & (fac >= 0.75)):
            #Do the cutting and analysis
            ndvi = imm / 255
            ndvi = ndvi.astype('float')

            ndvi_imgArr = ndvi
            #replace illegitimate values

            ndvi[ndvi <= 0] = np.nan
            ndvi[ndvi > 1] = np.nan
            data = ndvi.flatten()
            fdata = data[~np.isnan(data)]
            #print('Field Img Carve Date: '+ imgDtStr)


            mval_stats.loc[0,'Mean']=np.mean(fdata)
            mval_stats.loc[0,'mx']=np.max(fdata)
            mval_stats.loc[0,'p90']=np.percentile(fdata,90)
            mval_stats.loc[0,'p25']=np.percentile(fdata,25)
            mval_stats.loc[0,'p75']=np.percentile(fdata,75)
            mval_stats.loc[0,'p10']=np.percentile(fdata,10)
            mval_stats.loc[0,'mn']=np.min(fdata)
            mval_stats.loc[0,'Median']=np.median(fdata)
            mval_stats.loc[0,'StdDev']=np.std(fdata)
            mval_stats.loc[0,'Extent']=extJSON
            mval_stats.loc[0,'ImgDateStr']=imgDtStr

            break
        else:
            mval_stats = []
            continue
        #except:
        #	continue
        
    return [mval_stats, ndvi_imgArr, extJSON , imgDtStr]

def getNDVIValueWithStatsAndImg3_noscreen(df,polygon):
    #Returns the image (2-d numpy array) and the extent
    mval_stats = pd.DataFrame()
    ndvi_imgArr = []
    ndvi_extList = []
    imgDtStr = ''
    extJSON = ''

    for i,row in df.iterrows():
        
        tifFile = row.FILE
        f = tifFile
        ndvi_imgArr = []
        ndvi_extList = []
        imgDtStr = ''
    
        mval_stats = pd.DataFrame(columns=['Mean','p90','p25','p75','Median','StdDev','Extent','ImgDateStr'])

        roi_polygon = list(polygon.exterior.coords)

        try:
            if (tifFile.split('/')[1]=='mnt'):
                src = rasterio.open(f,'r', count=1)
            else:
                src = rasterio.open(f,'r', count=7)
                
            roi_polygon_src_coords = warp.transform_geom(
                {'init': 'epsg:4326'},
                src.crs,
                {'type': 'Polygon','coordinates': [roi_polygon]}
                )

            imgDtStr = getImgDateStr2(f)
            ext = getExtent(polygon)
            extJSON = json.dumps(ext)

            out_image,out_transform = mask.mask(src,[roi_polygon_src_coords],crop=True, all_touched=True)

            if (tifFile.split('/')[1]=='mnt'):
                imm = out_image[0,::]
            else:
                imm = out_image[4,::]

            polyac = area(mapping(polygon)) * 0.000247
            maskac = np.count_nonzero(imm)*100*0.000247
            fac = maskac / polyac

            if ((np.nansum(imm) > 0) & (fac >= 0.75)):
                #Do the cutting and analysis
                ndvi = imm / 255
                ndvi = ndvi.astype('float')

                ndvi_imgArr = ndvi
                #replace illegitimate values
                ######## <<<<<<<< ------ >>>>>> THIS IS THE CHANGE/SCREEN 
                #ndvi[ndvi <= 0] = np.nan
                #ndvi[ndvi > 1] = np.nan
                ########

                data = ndvi.flatten()
                fdata = data[~np.isnan(data)]
                #print('Field Img Carve Date: '+ imgDtStr)


                mval_stats.loc[0,'Mean']=np.mean(fdata)
                mval_stats.loc[0,'p90']=np.percentile(fdata,90)
                mval_stats.loc[0,'p25']=np.percentile(fdata,25)
                mval_stats.loc[0,'p75']=np.percentile(fdata,75)
                mval_stats.loc[0,'Median']=np.median(fdata)
                mval_stats.loc[0,'StdDev']=np.std(fdata)
                mval_stats.loc[0,'Extent']=extJSON
                mval_stats.loc[0,'ImgDateStr']=imgDtStr

                break
            else:
                mval_stats = []
                continue
        except:
            continue
        
    return [mval_stats, ndvi_imgArr, extJSON , imgDtStr]

def getNIIValueWithStatsAndImg3(df,polygon):
    
    mval_stats = pd.DataFrame(columns=['Mean','p90','p25','p75','Median','StdDev','Extent','ImgDateStr'])
    tifFile = ''
    f = ''
    nii_imgArr = []
    nii_extList = []
    imgDtStr = []
    extJSON=''

    for i,row in df.iterrows():
        mval_stats = pd.DataFrame(columns=['Mean','p90','p25','p75','Median','StdDev','Extent','ImgDateStr'])

        tifFile = row.FILE
        f = tifFile
        nii_imgArr = []
        nii_extList = []
        imgDtStr = []
        extJSON=''

        roi_polygon = list(polygon.exterior.coords)

        try:
            src = rasterio.open(f,'r', count=6)

                
            roi_polygon_src_coords = warp.transform_geom(
                {'init': 'epsg:4326'},
                src.crs,
                {'type': 'Polygon','coordinates': [roi_polygon]}
                )

            imgDtStr = getImgDateStr2(f)
            ext = getExtent(polygon)
            extJSON = json.dumps(ext)

            out_image,out_transform = mask.mask(src,[roi_polygon_src_coords],crop=True, all_touched=True)

            imm = out_image[5,::]


            polyac = area(mapping(polygon)) * 0.000247
            maskac = np.count_nonzero(imm)*100*0.000247
            fac = maskac / polyac

            if ((np.nansum(imm) > 0) & (fac >= 0.75)):
                #Do the cutting and analysis
                nii = imm / 255
                nii = nii.astype('float')

                nii_imgArr = nii

                #replace illegitimate values
                nii[nii <= 0] = np.nan
                #nii[nii > 1] = np.nan
                data = nii.flatten()
                fdata = data[~np.isnan(data)]

                mval_stats.loc[0,'Mean']=np.mean(fdata)
                mval_stats.loc[0,'p90']=np.percentile(fdata,90)
                mval_stats.loc[0,'p25']=np.percentile(fdata,25)
                mval_stats.loc[0,'p75']=np.percentile(fdata,75)
                mval_stats.loc[0,'Median']=np.median(fdata)
                mval_stats.loc[0,'StdDev']=np.std(fdata)
                mval_stats.loc[0,'Extent']=extJSON
                mval_stats.loc[0,'ImgDateStr']=imgDtStr

                break
            else:
                mval_stats = []
                continue
        except:
            continue
        
    return [mval_stats, nii_imgArr, extJSON , imgDtStr]


def getSinglePoly(r):
    if type(r)==shapely.geometry.polygon.Polygon:
        poly = r
    else:  
        n = len(r)
        poly = []
        if n>1:
            for i in range(n):
                p = r[i]
                if type(p)==shapely.geometry.polygon.Polygon:
                    poly=p
                else:
                    continue
        else:
            poly = []
            
    return poly

def removeIncompleteStns(weather_df,N):
    complete_weather_df = pd.DataFrame(columns = weather_df.columns)
    list_of_stns = pd.unique(weather_df.StationId)
    for stn in list_of_stns:
        stn_df = weather_df[weather_df.StationId==stn]
        if (len(stn_df)==N+1):
            complete_weather_df = complete_weather_df.append(stn_df)
        else:
            continue
    
    return complete_weather_df

def interpolateMissingVals(weather_df):
    #First drop the QA fields and index by Date and StantionId
    new_weather_df = pd.DataFrame()
    keepCols = ['StationId',
             'Date',
             'JulianDate',
             'Reference ETo(in)',
             'Precipitation(in)',
             'SolarRadiationAverage(Ly/day)',
             'AverageVaporPressure(mBars)',
             'Maximum Air Temperature(F)',
             'Minimum Air Temperature   (F)',
             'Average Air Temperature   (F)',
             'Maximum Relative Humidity   (%)',
             'Minimum Relative Humidity   (%)',
             'Average Relative Humidity   (%)',
             'Dew Point   (∞=F)',
             'Average Wind Speed   (mph)',
             'Wind Run   (miles)',
             'Average Soil Temperature   (F)'
            ]
    keepCols = [x.replace(' ','') for x in keepCols]
    keepCols = [x.replace('∞=','') for x in keepCols]
    df = weather_df[keepCols]
    list_of_stns = pd.unique(df.StationId)
    for stn in list_of_stns:
        stn_df = df[df.StationId==stn]
        stn_df = stn_df.replace('--',np.nan,method='bfill')
        stn_df = stn_df.astype({
            'StationId':'int',
            'Date':'datetime64[ns]',
            'JulianDate':'int64',
            'ReferenceETo(in)':'float64',
            'Precipitation(in)':'float64',
            'SolarRadiationAverage(Ly/day)':'float64',
            'AverageVaporPressure(mBars)':'float64',
            'MaximumAirTemperature(F)':'float64',
            'MinimumAirTemperature(F)':'float64',
            'AverageAirTemperature(F)':'float64',
            'MaximumRelativeHumidity(%)':'float64',
            'MinimumRelativeHumidity(%)':'float64',
            'AverageRelativeHumidity(%)':'float64',
            'DewPoint(F)':'float64',
            'AverageWindSpeed(mph)':'float64',
            'WindRun(miles)':'float64',
            'AverageSoilTemperature(F)':'float64'                      
                               })
        stn_df['Date'] = pd.to_datetime(stn_df['Date'])
        stn_df = stn_df.set_index('Date', drop=True)
        stn_df = stn_df.resample('D').mean()
        stn_df = stn_df.interpolate()
        stn_df.reset_index(inplace=True)
        
        new_weather_df = new_weather_df.append(stn_df)
    return new_weather_df

def getWeather2ForDates(p,list_of_dates,df):
    
    N=5 #default
    #gets the weather metrics for the METRIC calcs for a point for each of the dates based
    #of the closest N stns #given the geopandas geodataframe gpd
    
    #s = pt.split('POINT ')[1]
    #p = Point([float(s.split(' ')[0][1:-1]), float(s.split(' ')[1][0:-2])])
    
    #Screen out the null values

    df = df[df['Tavg']!='--']
    df = df[df['Rnet']!='--']
    df = df[df['Wavg']!='--']
    df = df[df['ReferenceETo(in)']!='--']
    df = df[df['Precipitation(in)']!='--'] 
    df = df[df['AverageSoilTemperature(F)']!='--']
    df = df[df['MaximumAirTemperature(F)']!='--']
    df = df[df['MinimumAirTemperature(F)']!='--']
    df = df[df['AverageAirTemperature(F)']!='--']
    df = df[df['MaximumRelativeHumidity(%)']!='--']
    df = df[df['MinimumRelativeHumidity(%)']!='--']
    df = df[df['AverageRelativeHumidity(%)']!='--']
    df = df[df['DewPoint(F)']!='--']
    df = df[df['WindRun(miles)']!='--']
    
    
    cols_list = ['Date',
                 'Tavg',
                 'Rnet',
                 'Wavg',
                 'ReferenceETo(in)',
                 'Precipitation(in)',
                 'AverageSoilTemperature(F)',
                 'MaximumAirTemperature(F)',
                 'MinimumAirTemperature(F)',
                 'AverageAirTemperature(F)',
                 'MaximumRelativeHumidity(%)',
                 'MinimumRelativeHumidity(%)',
                 'AverageRelativeHumidity(%)',
                 'DewPoint(F)',
                 'WindRun(miles)'
                ]
    df_weather_for_dates = pd.DataFrame(columns=cols_list)
    
    df.Date=pd.to_datetime(df.Date)
    
    
    
    idxx=0
    for dt in list_of_dates:
        df_dt = df[df.Date==dt]
        
        df_dt['StationId']=df_dt['StationId'].astype('uint16')
        stnPoints = list(pd.unique(df_dt.geometry))
        
        #change on 12-9-2019
        stns = pd.unique(df_dt.StationId)
        
        
        #print(stns)
        #first get the nearest stns and weights for the point provided
        [ids, pts, dis, wts] = getNearestNStnIDsInOrder2(p,stnPoints,5)
        #print([ids, pts, dis, wts])

        nearestStnIDs = stns[ids]
        
        weather=[]
        airTemp = []
        solarRad = []
        windSp = []
        eto=[]
        precip=[]
        avgSoilTemp=[]
        maxAirTemp=[]
        minAirTemp=[]
        avgAirTemp=[]
        maxRelHumPercent=[]
        minRelHumPercent=[]
        avgRelHumPercent=[]
        dewPoint=[]
        windRun=[]
                
        for stn in nearestStnIDs:
            airTemp.append(float(df_dt[df_dt.StationId==stn]['Tavg'].values.tolist()[0]))
            solarRad.append(float(df_dt[df_dt.StationId==stn]['Rnet'].values.tolist()[0]))
            windSp.append(float(df_dt[df_dt.StationId==stn]['Wavg'].values.tolist()[0]))
            eto.append(float(df_dt[df_dt.StationId==stn]['ReferenceETo(in)'].values.tolist()[0]))
            precip.append(float(df_dt[df_dt.StationId==stn]['Precipitation(in)'].values.tolist()[0]))
            avgSoilTemp.append(float(df_dt[df_dt.StationId==stn]['AverageSoilTemperature(F)'].values.tolist()[0]))
            maxAirTemp.append(float(df_dt[df_dt.StationId==stn]['MaximumAirTemperature(F)'].values.tolist()[0]))
            minAirTemp.append(float(df_dt[df_dt.StationId==stn]['MinimumAirTemperature(F)'].values.tolist()[0]))
            avgAirTemp.append(float(df_dt[df_dt.StationId==stn]['AverageAirTemperature(F)'].values.tolist()[0]))
            maxRelHumPercent.append(float(df_dt[df_dt.StationId==stn]['MaximumRelativeHumidity(%)'].values.tolist()[0]))
            minRelHumPercent.append(float(df_dt[df_dt.StationId==stn]['MinimumRelativeHumidity(%)'].values.tolist()[0]))
            avgRelHumPercent.append(float(df_dt[df_dt.StationId==stn]['AverageRelativeHumidity(%)'].values.tolist()[0]))
            dewPoint.append(float(df_dt[df_dt.StationId==stn]['DewPoint(F)'].values.tolist()[0]))
            windRun.append(float(df_dt[df_dt.StationId==stn]['WindRun(miles)'].values.tolist()[0]))
        
        weighted_airTemp = np.nansum(np.array(airTemp)*np.array(wts))
        weighted_solarRad = np.nansum(np.array(solarRad)*np.array(wts))
        weighted_windSp = np.nansum(np.array(windSp)*np.array(wts))
        weighted_eto = np.nansum(np.array(eto)*np.array(wts))
        weighted_precip = np.nansum(np.array(precip)*np.array(wts))
        weighted_avgSoilTemp = np.nansum(np.array(avgSoilTemp)*np.array(wts))
        weighted_maxAirTemp = np.nansum(np.array(maxAirTemp)*np.array(wts))
        weighted_minAirTemp = np.nansum(np.array(minAirTemp)*np.array(wts))
        weighted_avgAirTemp = np.nansum(np.array(avgAirTemp)*np.array(wts))
        weighted_maxRelHumPercent = np.nansum(np.array(maxRelHumPercent)*np.array(wts))
        weighted_minRelHumPercent = np.nansum(np.array(minRelHumPercent)*np.array(wts))
        weighted_avgRelHumPercent = np.nansum(np.array(avgRelHumPercent)*np.array(wts))
        weighted_dewPoint = np.nansum(np.array(dewPoint)*np.array(wts))
        weighted_windRun = np.nansum(np.array(windRun)*np.array(wts))
        
        df_weather_for_dates.loc[idxx,'Date']=dt
        df_weather_for_dates.loc[idxx,'Tavg']=weighted_airTemp
        df_weather_for_dates.loc[idxx,'Rnet']=weighted_solarRad
        df_weather_for_dates.loc[idxx,'Wavg']=weighted_windSp
        df_weather_for_dates.loc[idxx,'ReferenceETo(in)']=weighted_eto
        df_weather_for_dates.loc[idxx,'Precipitation(in)']=weighted_precip
        df_weather_for_dates.loc[idxx,'AverageSoilTemperature(F)']=weighted_avgSoilTemp
        df_weather_for_dates.loc[idxx,'MaximumAirTemperature(F)']=weighted_maxAirTemp
        df_weather_for_dates.loc[idxx,'MinimumAirTemperature(F)']=weighted_minAirTemp
        df_weather_for_dates.loc[idxx,'AverageAirTemperature(F)']=weighted_avgAirTemp
        df_weather_for_dates.loc[idxx,'MaximumRelativeHumidity(%)']=weighted_maxRelHumPercent
        df_weather_for_dates.loc[idxx,'MinimumRelativeHumidity(%)']=weighted_minRelHumPercent
        df_weather_for_dates.loc[idxx,'AverageRelativeHumidity(%)']=weighted_avgRelHumPercent
        df_weather_for_dates.loc[idxx,'DewPoint(F)']=weighted_dewPoint
        df_weather_for_dates.loc[idxx,'WindRun(miles)']=weighted_windRun
        
        idxx=idxx+1
          
    return df_weather_for_dates

def getEToForDates(pt,list_of_dates,df):
    N=5 #default
    #gets the ETo for a point for each of the dates based on the ETo of the closest N stns
    #given the geopandas geodataframe gpd
    
    p = pt
    df_eto_for_dates = pd.DataFrame(columns=['Date','ETo(in)'])
    
    df.Date=pd.to_datetime(df.Date)
    
    idxx=0
    for dt in list_of_dates:
        df_dt = df[df.Date==dt]
        df_dt['StationId']=df_dt['StationId'].astype('uint16')
        stnPoints = list(pd.unique(df_dt.geometry))
        stns = pd.unique(df_dt.StationId)
        
        #first get the nearest stns and weights for the point provided
        [ids, pts, dis, wts] = getNearestNStnIDsInOrder(p,stnPoints,5)
        #print([ids, pts, dis, wts])
        nearestStnIDs = stns[ids]
        
        etos=[]
        for stn in nearestStnIDs:
            
            eto = df_dt[df_dt.StationId==stn]['ReferenceETo(in)'].values.tolist()[0]
            #print(stn,eto)
            etos.append(float(eto))
            
        weighted_eto = sum(np.array(etos)*np.array(wts))
        
        df_eto_for_dates.loc[idxx,'Date']=dt
        df_eto_for_dates.loc[idxx,'ETo(in)']=weighted_eto
        idxx=idxx+1
    
    return df_eto_for_dates

def getFirstNonEmptyOrNonNaN_DF(ndvi_mList_df):
    df = ndvi_mList_df.copy()
    
    if len(df)>0:
        retRow = df.iloc[0]
    else:
        #just set everything to zero
        df = pd.DataFrame(columns=['Mean','p90','p25','p75','Median','StdDev'])
        for col in df.columns:
            df.loc[0,col] = np.nan
        retRow=df.iloc[0]
    
    return retRow

def getFirstNonEmptyOrNonNaN(mList):
    m = []
    
    #if the list is all NaNs, then return 0
    x = ~np.isnan(np.array(mList))*np.ones(len(mList))
    if (np.sum(x) == 0):
        return 0
    else:
        # remove the empty spaces
        mList = [x for x in mList if x]
    
        # remove NaNs
        idx_ok = ~np.isnan(mList)
        a = np.array(mList)
        b = a[idx_ok]
        mList = b.tolist()
     
        #now get the first on the top of the list
        m = next(s for s in mList if s)
        return m
    
def computeMetricET(Tavg,Rnet,Wavg,NDVI):
    #Tavg = Avg daily Temp in Deg-C
    #Rnet = Net Solar Radiation
    #Wavg = Avg Wind speed (m s-1)
    #NDVI is a single number which is the average of the field
    
    #Modifications as each is a vector of values (one for each point)
    
    #Constants:
    alpha=0.25
    lpd2wpm2= 0.484583
    
    #Tavg in C, already comes like that
    Rnet = Rnet*(lpd2wpm2/3)
    
    #Assumptions: Assuming Tsur = Tavg+1
    Tsur = Tavg+1
    
    #Filter out negative net radiation
    if Rnet < 0.0:
        Rnet=0.0
    
    #G is soil heat flux in W/m^2. alpha surface albedo. NDVI value per each pixel
    G=Rnet*Tsur*(0.0038+0.0074*alpha)*(1-0.98*NDVI**4)
    
    #Sensible heat flux:
    rhoa = 1.2 #air density (kg/m^3) http://www.cgd.ucar.edu/tss/aboutus/staff/bonan/ecoclim/1sted/Chapter07.pdf
    #Ts= air temperature in first reference height (surface of the crop)
    
    Cp = 1010 #Cp=heat capacity of air (Jkg-1K-1). http://www.cgd.ucar.edu/tss/aboutus/staff/bonan/ecoclim/1sted/Chapter07.pdf
    rah = 400/(1+0.536*Wavg) #rah
    
    H=rhoa*Cp*(Tsur-Tavg)/rah
    #Where, H=Sesible heat flux, rhoa= air density(W m-2)
    
    LE = Rnet-H-G
    A=LE/(Rnet-G)
    #LE= Latent heat flux (LE) [ LE=Rnet-H-G]
    
    #Actual ET
    gamma = 2270; #gamma= latent heat of vaporization (Jkg-1)
    Rho_w = 1000*(1-(Tavg+288.9414)/508929/(Tavg+68.12963)*(Tavg-3.9863)**2) #Rho_w= density of water (kg m-3)
    #Equation from http://calculator.tutorvista.com/water-density-calculator.html

    #As we have daily Rnet24 and G24 are same as Rnet and G respectively ????
    Rnet24=Rnet;
    G24=G;

    MetricET=86400*A*(Rnet24-G24)/(gamma*Rho_w)*0.03937 #in inches

    #where,
    #A=soil humidity constant
    #Rnet= daily net radiation (Wm-2)
    #G24 = daily soil heat flux(Wm-2)
    #gamma= latent heat of vaporization (Jkg-1)
    #Rho_w= density of water (kg m-3)

    return MetricET


def getNearestFile(sentinel_files_array,dtstr):
    
    dates_df = pd.DataFrame()
    res_df = pd.DataFrame()
    dts=[]
    fileNames =[]
    
    for file in sentinel_files_array:
        d_str = file.split('_')[2][0:8]
        dt = parse(d_str)
        dts.append(dt)
        fileNames.append(file)
    
    dates_df = pd.DataFrame(dts,columns=['Date']) 
    dates_df['Date']=pd.to_datetime(dates_df.Date)
    dates_df.sort_values(by='Date', inplace=True)
    old_idx = dates_df.index #get the old index
    new_df = dates_df.reset_index(drop=True).copy()
    
    res_df = new_df[new_df.Date <= pd.Timestamp(dtstr)]
    idx = res_df.iloc[-1].name
    filenm = fileNames[old_idx[idx]]
    
    return filenm

def getMostRecentFile(sentinel_files_array):
    dts = []
    for file in sentinel_files_array:
        d_str = file.split('_')[2][0:8]
        dt = parse(d_str)
        dts.append(dt)
    idx = [dts.index(x) for x in sorted(dts)]
    return sentinel_files_array[idx[-1]]

def getFootprint(XML_file):
    f = open(XML_file, "r")
    lines = f.readlines()
    s=0
    n=0
    crd = []

    for line in lines:
        n=n+1
        if (line.strip().startswith('<Global_Footprint>')):
            s=1
            n=0
        if (line.strip().startswith('</Global_Footprint>')):
            s=0
        if s==1 & n==1:
            l = line.strip()
            t = l.split('<EXT_POS_LIST>')
            c = t[1][:-15].strip()
            crds = c.split(' ')
        
            for i in range(int(len(crds)/2)):
                lat = float(crds[i*2])
                lon = float(crds[i*2+1])
                ll = [lon,lat]
                crd.append(ll)
    return crd

def getAreaOfROI(roi_polygon):
    #Returns area in KM^2 and the fraction of the 100x100km Sentinel
    import pyproj    
    import shapely
    import shapely.ops as ops
    from shapely.geometry.polygon import Polygon
    from functools import partial

    geom = shapely.geometry.Polygon(roi_polygon)
    #geom = Polygon([(0, 0), (0, 10), (10, 10), (10, 0), (0, 0)])
    geom_area = ops.transform(
        partial(
            pyproj.transform,
            pyproj.Proj(init='EPSG:4326'),
            pyproj.Proj(
                proj='aea',
                lat1=geom.bounds[1],
                lat2=geom.bounds[3])),
        geom)

    # return the area in Km^2
    a = (geom_area.area / 1000000)
    f = a / 12056.04 #109.8x109.8 km = 12,056.04 km ^2
    return [a,f]

def Convert(s): 
    li = list(eval(s))
    return li

def getNearestNStnIDsInOrder2(c,stnPoints,N):
    #returns the indexes, 
    #the sorted list of stns, 
    #the distances and 
    #the weights (w=1/dis)
    
    
    dist = []
    #stnPts = MultiPoint(stnPoints)
    #nearest_geoms = nearest_points(point, stnPoints) #ordered list of shapely points
    for p in stnPoints:
        d = float(c.distance(p))
        if d==0:
            d=0.000001
        else:
            dist.append(d)
    
    idx = np.array([i[0] for i in sorted(enumerate(dist), key=lambda x:x[1])])

    orderedStns = [stnPoints[i] for i in idx]
    
    orderedDist = [dist[i] for i in idx]
    
    M = min(N, len(stnPoints))

    orderedWghts=[]
    invDist = []
    for d in orderedDist[0:M]:
        invDist.append(1/d)
    den = sum(invDist)   
    for w in invDist:
        orderedWghts.append(w/den)
    
    return [idx[0:M],orderedStns[0:M],orderedDist[0:M],orderedWghts]

def getTileStr(tileListStr):
    tilestr = tileListStr.strip(']["')
    return tilestr

def getShapelyPolygon(PolyString):
    #cleanup any beginning or training parentheses
    newstr = PolyString.replace('(','').replace(')','')
    s = newstr.split('POLYGON Z ')[1]
    coords = s.split(',')
    lat = []
    lon = []
    ll = []
    for c in coords:
        n = c.strip().split(' ')
        lon = float(n[0])
        lat = float(n[1])
        ll.append([lon,lat])
        
    shapelyPolygon = shapely.geometry.Polygon(ll)
    
    return shapelyPolygon

def getAttributesJSON(pandasDataSeries):
    attrJSON = ''
    excluded_col_list = ['geometry', 'Sentinel2_TileList']
    ds = pandasDataSeries.drop(excluded_col_list)
    attrJSONstr = ds.to_json()
    
    return attrJSONstr

"""
def lookupAcctName(acct_gdf2,customerID):
    gdf = acct_gdf2.copy()
    catid = 'JAINLOGICCUSTOMERID_'+str(customerID)
    gdff = gdf[gdf.Category.str.upper()==catid]
    acctName = pd.unique(gdff.SubCategory)[0]
    #.replace(' ','').upper()
   
    return acctName
"""

def lookupAcctName(acct_gdf2,customerID):
    gdf = acct_gdf.copy()
    catid = str(customerID)
    gdff = gdf[gdf.customerid==int(catid)]
    acctName = pd.unique(gdff.SubCategory)[0].replace(' ','').upper()
   
    return acctName


def lookupRanchName(acct_gdf2,customerID):
    gdf = acct_gdf2.copy()
    #catid = str(customerID)
    gdff = gdf[gdf.customerid==customerID]
    if len(gdff)>0:
        ranchName = pd.unique(gdff.SubCategory)[0]
    else:
        ranchName = []
    
    return ranchName

def lookupFieldName(acct_gdf2,customerID,assetID):
    gdf = acct_gdf2.copy()
    gdff = gdf[(gdf.customerid==customerID)  & (gdf.fieldassetid==assetID)]
    
    if len(gdff)>0:
        fieldName = pd.unique(gdff.Name)[0]
    else:
        fieldName = []
        
    return fieldName

def lookupRanchID(acct_gdf2,ranchName):
    gdf = acct_gdf2.copy()
    gdff = gdf[gdf.SubCategory.str.replace(' ','').str.upper()==ranchName]
    ranchID = pd.unique(gdff.Category)[0].replace(' ','').upper()
    
    return ranchID

def lookupFieldID(acct_gdf2,ranchName,fieldName):
    gdf = acct_gdf2.copy()
    gdff = gdf[(gdf.SubCategory.str.replace(' ','').str.upper()==ranchName) & (gdf.Name.str.replace(' ','').str.upper()==fieldName)]
    fieldID = pd.unique(gdff.fieldassetid)[0]
    
    return fieldID


    
def getSmallerEvenNo(num):
    if (num % 2) == 0:
        retNum = num
    else:
        retNum = num - (num % 2)
    return retNum

def getTrimmedImg(imgArr):
    s = imgArr.shape
    x = tuple([getSmallerEvenNo(s[0]), getSmallerEvenNo(s[1])])
    imgArr_new = imgArr[0:x[0],0:x[1]]
    return imgArr_new

def getReTrimmedImg(img1,img2,img3):
    s1 = img1.shape
    s2 = img2.shape
    s3 = img3.shape
    
    row_min = min(s1[0],s2[0],s3[0])
    col_min = min(s1[1],s2[1],s3[1])
    
    img1_trimmed = img1[0:row_min,0:col_min]
    img2_trimmed = img2[0:row_min,0:col_min]
    img3_trimmed = img3[0:row_min,0:col_min]
    
    return [img1_trimmed,img2_trimmed,img3_trimmed]

# In[347]:
def getResizedImgs(img1, img2):
    s1 = img1.shape
    s2 = img2.shape
    
    row_1 = s1[0]
    row_2 = s2[0]
    col_1 = s1[1]
    col_2 = s2[1]
    
    delta_rows = delta_cols = 0
    
    if (row_1 > row_2):
        delta_rows = int((row_1 - row_2) / 2)
        #delta_cols = int((col_1 - col_2) / 2)
        #img2 = np.pad(img2,(delta_rows,delta_cols),'edge')
    if (row_2 > row_1):
        delta_rows = int((row_2 - row_1) / 2)
        #delta_cols = int((col_2 - col_1) / 2)
        #img1 = np.pad(img1,(delta_rows,delta_cols),'edge')
    if (col_1 > col_2):
        delta_cols = math.ceil((col_1 - col_2) / 2)
        #print(delta_cols)
    if (col_2 > col_1):
        delta_cols = math.ceil((col_2 - col_1) / 2)
        #print(delta_cols)

    img1 = np.pad(img1,(delta_rows,delta_cols),'edge')
    img2 = np.pad(img2,(delta_rows,delta_cols),'edge')
        
    return (img1, img2)

def dropCCRows(pkl_df):
    for i, row in pkl_df.iterrows():
        f = row.PKLFilePath
        ranch = row.Ranch
        field = row.Field
        loc = row.Location
    
        infile = open(f,'rb')
        img = pickle.load(infile)
        infile.close()
        #print([img.size, loc])
        pkl_df.loc[i,'Size']=img.size
    
    LL=pkl_df.groupby('Location')['Size'].unique()
    loc_list = LL.index.tolist()
    ind=0
    pklnew_df = pd.DataFrame()
    for i, row in pkl_df.iterrows():
        loc = row.Location
        val = max(LL[loc])
        pkl_df.loc[i,'MaxSize']=val
    pkl_df2 = pkl_df[pkl_df.Size == pkl_df.MaxSize]

    return pkl_df2

def color_above10_red(value):
    color='black'
    if value >= 10:
        color = 'red'
    return ('color: %s' % color)

def color_above10_green(value):
    color='black'
    if value >= 10:
        color = 'green'
    return ('color: %s' % color)


def getLatestDt(list_of_files):
    """
    returns the most recent file in this list of files
    """
    dts=[]
    for f in list_of_files:
        dtStr = f.split('.')[0][-10:]
        dt = datetime.strptime(dtStr,'%Y-%m-%d')
        dts.append(dt)
        
    dts.sort()
    latestDate = dts[-1]
    
    return latestDate


# In[114]:


def color_above10_red(value):
    color='black'
    if value >= 10:
        color = 'red'
    return ('color: %s' % color)


# In[116]:


def color_above10_green(value):
    color='black'
    if value >= 10:
        color = 'green'
    return ('color: %s' % color)


# In[117]:


def getLinkTo(pdfFile):
    
    file_path = Path(pdfFile)
    linkURL='https://www.agralogics.com/HyperView/'+file_path.name
    linkStr = '<a href=\"'+linkURL+'\">HyperView Report</a>'  
    
    return linkStr


def doHyperViewFTP(pdfFile):
    url = 'ftp.agralogics.com'
    username = 'agls@agralogics.com'
    passwd = 'AglsAdmin!2019'
    ftpFolder='HyperView/'

    file_path = Path(pdfFile)

    with FTP(url, username, passwd, timeout=100000) as ftp, open(file_path, 'rb') as file:
        ftp.storbinary(f'STOR {ftpFolder+file_path.name}', file)
    
    return


def doFTP(pdfFile):
    url = 'ftp.agralogics.com'
    username = 'agls@agralogics.com'
    passwd = 'AglsAdmin!2019'
    ftpFolder='HyperView/'

    file_path = Path(pdfFile)

    with FTP(url, username, passwd) as ftp, open(file_path, 'rb') as file:
        ftp.storbinary(f'STOR {ftpFolder+file_path.name}', file)
    
    return


# In[119]:


def getBase64EncodingForPDF(pdfFile):
    with open(pdfFile, "rb") as pdf_file:
        encoded_string = base64.b64encode(pdf_file.read())
    
    imgstr = '<img src=\"data:application/pdf;base64, '+encoded_string+'\">'
    
    retstr = '<a href=\"www.google.com\">'+imgstr+'</a>'
    return retstr


# In[120]:


def genPDFEmbededLink(pdfFilePath):
    #pre = '<div style=\"padding:20px\">'
    #embed = '<embed src=\"'+pdfFilePath+'\" type=\"application/pdf\" width=\"100%\", height=\"600px\"></embed>'
    #post = '</div>'
    #linkedStrAsCellValue = pre+embed+post
    
    pre = '<object data=\"'+pdfFilePath+'\" type=\"application/x-pdf\" title=\"HyperView\" width=\"500\" height=\"720\">'
    ref = '<a href=\"'+pdfFilePath+'\">HyperViewReport</a>'
    post='</object>'
    
    linkedStrAsCellValue = pre+ref+post
    
    return linkedStrAsCellValue


# In[121]:


def makeTrueRedFalseBlank(html_str):
    
    old_str = 'True'
    new_str = '<font color=\"red\">'+old_str+'</font>'
    
    old_str2 = 'False'
    new_str2 = ''
    
    ret_str = html_str.replace(old_str,new_str)
    ret_str2 = ret_str.replace(old_str2,new_str2)
    
    return ret_str2


# In[122]:


def parse_point(record):
    pieces = record.split() # splits each record into a list of 3
    x = float(pieces[2].lstrip('(')) # latitude
    y = float(pieces[3].rstrip('')) # longitude 
    z = float(pieces[4].rstrip(')')) # elevation 
    point = Point(x,y,z) # convert to Shapely Point
    return point

def unparse_point(record):
    pieces = record.split() # splits each record into a list of 3
    x = float(pieces[2].lstrip('(')) # latitude
    y = float(pieces[3].rstrip('')) # longitude 
    z = float(pieces[4].rstrip(')')) # elevation 
    point = [x,y,z] # convert to coordinates 
    return point

def add_basemap(ax, zoom, url):
    #url=ctx.sources.ST_TERRAIN_BACKGROUND):
    #url='http://tile.stamen.com/terrain/tileZ/tileX/tileY.png'):
    #url=ctx.tile_providers.OSM_A):
    xmin, xmax, ymin, ymax = ax.axis()
    basemap, extent = ctx.bounds2img(xmin, ymin, xmax, ymax, zoom=zoom, url=url)
    ax.imshow(basemap, extent=extent, interpolation='bilinear')
    ax.axis((xmin, xmax, ymin, ymax))


def getAllSentinelScenePolygons():
    sentinelFMLTileFile = data_dir + 'S2A_OPER_GIP_TILPAR_MPC__20151209T095117_V20150622T000000_21000101T000000_B00.kml'
    gdf = gpd.read_file(sentinelFMLTileFile, driver='KML')
    gdf.rename(columns={"Name": "TILE"}, inplace=True)
    sentinelTiles_gdf = gdf[['TILE','geometry']]

    
    return sentinelTiles_gdf

#Read ALL Landsat Tiles 
def getAllLandsatScenePolygons():
    landsat_shapefile  = data_dir + 'wrs2.shp'
    landsat_gdf = gpd.GeoDataFrame.from_file(landsat_shapefile)
    
    landsat_gdf['TILE'] = landsat_gdf['PATH'].astype('str').str.rjust(3,'0') + landsat_gdf['ROW'].astype('str').str.rjust(3,'0')
        
    l_gdf = landsat_gdf[['TILE','geometry']]
    l_gdf.reset_index(drop=True,inplace=True)
    
    return l_gdf


#Read US Sentinel Tiles
def getAllUSSentinelScenePolygons(usap):
    sentinelFMLTileFile = data_dir + 'S2A_OPER_GIP_TILPAR_MPC__20151209T095117_V20150622T000000_21000101T000000_B00.kml'
    sentinel_gdf = gpd.read_file(sentinelFMLTileFile, driver='KML')
    sentinel_gdf = sentinel_gdf.drop_duplicates(subset = 'Name')

    j=0
    SentinelUSTiles_df = pd.DataFrame(columns=['TILE','geometry']) 
    for i,tile in sentinel_gdf.iterrows():
        sentp = tile.geometry[0]
        tileName = tile.Name
        if (sentp.intersects(usap) | sentp.within(usap)):
            SentinelUSTiles_df.loc[j,'TILE']=tileName
            SentinelUSTiles_df.loc[j,'geometry']=sentp
            j=j+1 
                
    SentinelUSTiles_gdf = gpd.GeoDataFrame(SentinelUSTiles_df,geometry=SentinelUSTiles_df.geometry)
    s_gdf = SentinelUSTiles_gdf.copy()
    s_gdf.reset_index(drop=True,inplace=True)
    
    return s_gdf

#Read US Landsat Tiles 
def getAllUSLandsatScenePolygons(usap):
    landsat_shapefile  = data_dir + 'wrs2.shp'
    landsat_gdf = gpd.GeoDataFrame.from_file(landsat_shapefile)
    
    j=0
    LandsatUSTiles_df = pd.DataFrame(columns=['TILE','geometry']) 
    for i,tile in landsat_gdf.iterrows():
        landp = tile.geometry
        path = str(tile.PATH).rjust(3,'0')
        row = str(tile.ROW).rjust(3,'0')
        tileName = path+row
        if (landp.intersects(usap) | landp.within(usap)):
            LandsatUSTiles_df.loc[j,'TILE']=tileName
            LandsatUSTiles_df.loc[j,'geometry']=landp
            j=j+1
                
    LandsatUSTiles_gdf = gpd.GeoDataFrame(LandsatUSTiles_df,geometry=LandsatUSTiles_df.geometry)
    l_gdf = LandsatUSTiles_gdf.copy()
    l_gdf.reset_index(drop=True,inplace=True)
    
    return l_gdf

def makeJsonStr(key,val):
    jsonstr = []
    data = {}
    data[key] = val
    jsonstr = json.dumps(data)
    return jsonstr

def makeStrFromList(tileList):
    tileStr = ''
    i=0
    for tile in tileList:
        if i==0:
            tileStr = tile
        else:
            tileStr = tileStr + ',' + tile
        i=i+1
    return tileStr

def getAcres(WGS84_ShapelyPolygon):
    ac = []
    ac = area(mapping(g)) * 0.000247
    return ac


def makePoly(y):
    retList = []
    pairs = y.split(',').strip()
    for j in pairs:
        retList.append([j])
    return retList

def add_basemap(ax, zoom, url):
    xmin, xmax, ymin, ymax = ax.axis()
    basemap, extent = ctx.bounds2img(xmin, ymin, xmax, ymax, zoom=zoom, url=url)
    ax.imshow(basemap, extent=extent, interpolation='bilinear')
    ax.axis((xmin, xmax, ymin, ymax))
    return


def pullJainlogicData(NDays):
    
    #Get the Jainlogic token
    print('\nGetting Jainlogic Token')
    url = 'https://login.jainlogic.com/oauth/token'
    headers = {
        'Content-Type': 'application/json',
        'Cookie': 'did=s%3Av0%3A78554490-cc77-11ea-93fa-bb6dc2a87cd5.5aeNuB2QvAJUoB1aaB5RcQewIEcIgXZAAbDDP41%2BVDs; did_compat=s%3Av0%3A78554490-cc77-11ea-93fa-bb6dc2a87cd5.5aeNuB2QvAJUoB1aaB5RcQewIEcIgXZAAbDDP41%2BVDs'
    }
    data = {
        'grant_type': 'password',
        'username' : 'agralogicspipeline@jainsusa.com',
        'password' : 'P2UrEnwUq1UZzkPZzMWH',
        'audience': 'jain-logic-apis',
        'scope' : 'openid offline_access',
        'client_id': 'ebZs4I55JaJSyB3MD42RgqPxG1twdVOu',
        'client_secret': 'Yrbz2uEkCO_NO1CltbuhYISrxWdSPfUmhZ3TvOz5UmhwC4jDGdStgaLnxB8r-yre'
    }
    r = requests.post(url, data=json.dumps(data), headers=headers, verify=False)
    token = json.loads(r.text)
    print('Done!')

    #Get the data
    print('\nPulling Latest Jainlogic Accts and Assets')
    #url = 'https://apps.puresense.com/irrigationapi/assets/irrigationsets'
    url = 'https://api.jainlogic.com/irrigation/assets/irrigationsets'
    headers = {
        'Content-Type': 'application/json',
        'Authorization': token['token_type'] + ' ' + token['access_token']
    }

    r = requests.get(url, headers=headers)
    if (r.status_code == 200):
        print('\tPull Successful')
    else:
        print('\tPull Unsuccessful')
    print('Done!')

    #Read the field definitions
    print('\nReading Field Definitions...')
    fieldDefs_df = pd.DataFrame()
    for i in range(len(json.loads(r.text))):
        d = json.loads(r.text)[i]
        df = pd.DataFrame.from_dict(d, orient='index').T
        fieldDefs_df = fieldDefs_df.append(df, ignore_index=True)
    fieldDefs_df.geographylastupdated = pd.to_datetime(fieldDefs_df.geographylastupdated)
    for i,row in fieldDefs_df.iterrows():
        P = wkt.loads(row.wkt)
        fieldDefs_df.loc[i,'Acres']= area(mapping(P))*0.000247
    #Clean up the wkt    
    fieldDefs_df = fieldDefs_df[~fieldDefs_df['wkt'].isnull()] 
    fieldDefs_df['wkt'] = fieldDefs_df['wkt'].apply(wkt.loads)
    print('Done!')

    #Filter the ones updated in the last "NDays" days
    #NDays = 10
    """
    if (type(NDays)==int):
        updated_df = fieldDefs_df[
            fieldDefs_df.geographylastupdated>=datetime.today()-timedelta(days=NDays)
        ]
        updated_df.reset_index(drop=True, inplace=True)
    else: #get all
        updated_df = fieldDefs_df.copy()
    """
    updated_df = fieldDefs_df.copy()

    #Cleanup
    updated_df = updated_df.sort_values(["customerid", "fieldassetid"])
    updated_gdf = gpd.GeoDataFrame(updated_df, geometry=updated_df.wkt)
    updated_gdf.drop(columns=['wkt'],inplace=True)
    updated_gdf.sort_values(by=['customerid', 'fieldassetid'], inplace=True)
    updated_gdf.reset_index(drop=True,inplace=True)
    updated_gdf.crs = {'init' :'epsg:4326'} #initialize the CRS

    #Make copy for next operations
    gdf_acct = updated_gdf.copy()

    #print(pd.unique(gdf_acct.customerid).tolist())
    


    #Screen out anything **other than** Polygon or MultiPolygon
    gdf_acct = gdf_acct[(gdf_acct.geom_type=='Polygon') | (gdf_acct.geom_type=='MultiPolygon')]
    gdf_acct.reset_index(drop=True,inplace=True)

    #Get only the Polygons first
    gdf_onlyPolygons = gdf_acct[gdf_acct.geom_type == 'Polygon']
    gdf_onlyPolygons.reset_index(drop=True,inplace=True)

    #Explode out any MultiPolygons and select just their first Polygon
    gdf_multiPolygons = gdf_acct[gdf_acct.geom_type == 'MultiPolygon']
    for i,row in gdf_multiPolygons.iterrows():
        polygons = list(row.geometry)
        if len(polygons)>1:
            p = polygons[0]
            gdf_multiPolygons.loc[i,'geometry']=p
    gdf_multiPolygons.reset_index(drop=True,inplace=True)

    #Now concat both 
    gdf_UpdatedAccts = pd.concat([gdf_onlyPolygons, gdf_multiPolygons], axis=0)
    gdf_UpdatedAccts.sort_values(by=['customerid', 'fieldassetid'], inplace=True)
    gdf_UpdatedAccts.reset_index(drop=True,inplace=True)

    #update Acres
    for i,row in gdf_UpdatedAccts.iterrows():
        P = row.geometry
        gdf_UpdatedAccts.loc[i,'Acres']= area(mapping(P))*0.000247


    return gdf_UpdatedAccts

def setupAGLSDataFrame(gdf_UpdatedAccts):
    gdf_acct = gdf_UpdatedAccts.copy()

    print('\nConstructing AGLS data structures...')
    #Create a 'Name' (field name), SubCategory (Customer), and Location,column to identify
    gdf_acct['Name']=gdf_acct['fieldname']
    gdf_acct['SubCategory']=gdf_acct['customername']
    gdf_acct['Location']=gdf_acct['SubCategory'].map(str) + ': ' + gdf_acct['Name']
    #Set the TS as a string
    gdf_acct.geographylastupdated =  gdf_acct.geographylastupdated.astype('str')
    gdf_acct.reset_index(drop=True, inplace=True)
    print('Done!')

    #Screen for California polygons only
    print('\nCScreening for California ONLY ...')
    bool_list = gdf_acct.geometry.apply(lambda g: g.within(cap))
    gdf_acct = gdf_acct.loc[bool_list]
    gdf_acct.reset_index(drop=True, inplace=True)
    print('Done!')

    print('\nCScreening for MONITORED ONLY ...')
    #Screen for ONLY MONITORED, Accounts
    gdf_acct = gdf_acct[(gdf_acct.monitored==True)]
    gdf_acct.reset_index(drop=True, inplace=True)
    #& (gdf_acct.subscribedhypergrow==True) & (gdf_acct.subscribedhyperview==True)]
    print('Done!')

    return gdf_acct

def setupAGLSDataFrameForRegion(gdf_UpdatedAccts, region, ismonitored):
    #print(region)
    gdf_acct = gdf_UpdatedAccts.copy()

    print('\nConstructing AGLS data structures...')
    #Create a 'Name' (field name), SubCategory (Customer), and Location,column to identify
    gdf_acct['Name']=gdf_acct['fieldname']
    gdf_acct['SubCategory']=gdf_acct['customername']
    gdf_acct['Location']=gdf_acct['SubCategory'].map(str) + ': ' + gdf_acct['Name']
    #Set the TS as a string
    gdf_acct.geographylastupdated =  gdf_acct.geographylastupdated.astype('str')
    gdf_acct.reset_index(drop=True, inplace=True)
    print('Done!')

    list_of_valid_regions = ['CA','WA','FL','TX','AZ','AR', 'DE', 'MO', 'MT', 'OK', 'OR','AUS','US', 'IND', 'MEX','CAN']
    if not region in list_of_valid_regions:
        print('Please provide specific proper region...aborting')
        return []
    else:
        if region=='CA':
            #Screen for California polygons only
            print('\n\tScreening for CA ONLY ...')
            bool_series = gdf_acct.geometry.apply(lambda g: g.within(cap))
            gdf_acct = gdf_acct.loc[bool_series]
            gdf_acct.reset_index(drop=True, inplace=True)
            print('\tDone!')
        elif region=='AZ':
            #Screen for AZ polygons only
            print('\n\tScreening for AZ ONLY ...')
            bool_series = gdf_acct.geometry.apply(lambda g: g.within(azp))
            gdf_acct = gdf_acct.loc[bool_series]
            gdf_acct.reset_index(drop=True, inplace=True)
            print('\tDone!')
        elif region=='WA':
            #Screen for WA polygons only
            print('\n\tScreening for WA ONLY ...')
            bool_series = gdf_acct.geometry.apply(lambda g: g.within(wap))
            gdf_acct = gdf_acct.loc[bool_series]
            gdf_acct.reset_index(drop=True, inplace=True)
            print('\tDone!')
        elif region=='TX':
            #Screen for WA polygons only
            print('\n\tScreening for TX ONLY ...')
            bool_series = gdf_acct.geometry.apply(lambda g: g.within(txp))
            gdf_acct = gdf_acct.loc[bool_series]
            gdf_acct.reset_index(drop=True, inplace=True)
            print('\tDone!')
        elif region=='FL':
            #Screen for WA polygons only
            print('\n\tScreening for FL ONLY ...')
            bool_series = gdf_acct.geometry.apply(lambda g: g.within(flp))
            gdf_acct = gdf_acct.loc[bool_series]
            gdf_acct.reset_index(drop=True, inplace=True)
            print('\tDone!')
        elif region=='OR':
            print('\n\tScreening for OR ONLY ...')
            bool_series = gdf_acct.geometry.apply(lambda g: g.within(orp))
            gdf_acct = gdf_acct.loc[bool_series]
            gdf_acct.reset_index(drop=True, inplace=True)
            print('\tDone!')
        elif region=='OK':
            #Screen for WA polygons only
            print('\n\tScreening for OK ONLY ...')
            bool_series = gdf_acct.geometry.apply(lambda g: g.within(okp))
            gdf_acct = gdf_acct.loc[bool_series]
            gdf_acct.reset_index(drop=True, inplace=True)
            print('\tDone!')
        elif region=='MT':
            #Screen for WA polygons only
            print('\n\tScreening for MT ONLY ...')
            bool_series = gdf_acct.geometry.apply(lambda g: g.within(mtp))
            gdf_acct = gdf_acct.loc[bool_series]
            gdf_acct.reset_index(drop=True, inplace=True)
            print('\tDone!')
        elif region=='MO':
            #Screen for WA polygons only
            print('\n\tScreening for FL ONLY ...')
            bool_series = gdf_acct.geometry.apply(lambda g: g.within(mop))
            gdf_acct = gdf_acct.loc[bool_series]
            gdf_acct.reset_index(drop=True, inplace=True)
            print('\tDone!')
        elif region=='DE':
            #Screen for WA polygons only
            print('\n\tScreening for DE ONLY ...')
            bool_series = gdf_acct.geometry.apply(lambda g: g.within(dep))
            gdf_acct = gdf_acct.loc[bool_series]
            gdf_acct.reset_index(drop=True, inplace=True)
            print('\tDone!')
        elif region=='AR':
            #Screen for WA polygons only
            print('\n\tScreening for AR ONLY ...')
            bool_series = gdf_acct.geometry.apply(lambda g: g.within(arp))
            gdf_acct = gdf_acct.loc[bool_series]
            gdf_acct.reset_index(drop=True, inplace=True)
            print('\tDone!')
        elif region=='AUS':
            #Screen for AUS polygons only
            print('\n\tScreening for AUS ONLY ...')
            a_gdf = gpd.GeoDataFrame()
            acct_gdf = gpd.GeoDataFrame()
            for ausp in ausp_list:
                bool_series = gdf_acct.geometry.apply(lambda g: g.within(ausp))
                a_gdf = gdf_acct.loc[bool_series]
                acct_gdf = acct_gdf.append(a_gdf, ignore_index=True)
            gdf_acct = acct_gdf.copy()
            gdf_acct.reset_index(drop=True, inplace=True)
            print('\tDone!')
        elif region=='IND':
            #Screen for IND polygons only
            print('\n\tScreening for IND ONLY ...')
            a_gdf = gpd.GeoDataFrame()
            acct_gdf = gpd.GeoDataFrame()
            for indp in indp_list:
                bool_series = gdf_acct.geometry.apply(lambda g: g.within(indp))
                a_gdf = gdf_acct.loc[bool_series]
                acct_gdf = acct_gdf.append(a_gdf, ignore_index=True)
            gdf_acct = acct_gdf.copy()
            gdf_acct.reset_index(drop=True, inplace=True)
            print('\tDone!')

        elif region=='MEX':
            #Screen for IND polygons only
            print('\n\tScreening for MEX ONLY ...')
            a_gdf = gpd.GeoDataFrame()
            acct_gdf = gpd.GeoDataFrame()
            for mexp in mexp_list:
                bool_series = gdf_acct.geometry.apply(lambda g: g.within(mexp))
                a_gdf = gdf_acct.loc[bool_series]
                acct_gdf = acct_gdf.append(a_gdf, ignore_index=True)
            gdf_acct = acct_gdf.copy()
            gdf_acct.reset_index(drop=True, inplace=True)
            print('\tDone!')
        elif region=='CAN':
            #Screen for CAN polygons only
            print('\n\tScreening for CAN ONLY ...')
            a_gdf = gpd.GeoDataFrame()
            acct_gdf = gpd.GeoDataFrame()
            for canp in canp_list:
                bool_series = gdf_acct.geometry.apply(lambda g: g.within(canp))
                a_gdf = gdf_acct.loc[bool_series]
                acct_gdf = acct_gdf.append(a_gdf, ignore_index=True)
            gdf_acct = acct_gdf.copy()
            gdf_acct.reset_index(drop=True, inplace=True)
            print('\tDone!')

        if (ismonitored):
            print('\nScreening for MONITORED ONLY ...')
            #Screen for ONLY MONITORED, Accounts
            gdf_acct = gdf_acct[(gdf_acct.monitored==True)]

        gdf_acct.reset_index(drop=True, inplace=True)
        print('Done!')

        return gdf_acct

def getAcctNameForJLID(customerid):
    #load the acct name for an id 
    gdf_acct = readUpdatedAcctsConfig([customerid])
    acctNameForID = pd.unique(gdf_acct.customername).tolist()
    return acctNameForID

def getExistingPilotAcctListForRegion(region):

    list_of_existing_acctFolders = [f.path.split('/')[-1] for f in os.scandir(pilotAcctConfigPath) if f.is_dir()]
    #print(list_of_existing_acctFolders)

    list_of_existing_accts = []
    for f in list_of_existing_acctFolders:
        folder = pilotAcctConfigPath +f
        #print(folder)
        configFile = folder+ '/'+ f+ '.json'
        #print(configFile)
        if os.path.exists(configFile):
            list_of_existing_accts.append(f)

    list_of_existing_accts.sort()
    #print(list_of_existing_accts)

    #load the acct list into a dataframe
    gdf_acct = readPilotAcctsConfig(list_of_existing_accts)

    #print(gdf_acct)

    #Screen for the region
    list_of_valid_regions = ['CA','WA','FL','TX','AZ','AUS','US','IND', 'MEX','EGY', 'CAN']
    if not region in list_of_valid_regions:
        print('Please provide specific proper region...aborting')
        return []
    else:
        if region=='IND':
            #Screen for India polygons only
            print('\nScreening for IND ONLY ...')
            bool_series = gdf_acct.geometry.apply(lambda g: g.within(indp))
            gdf_acct = gdf_acct.loc[bool_series]
            gdf_acct.reset_index(drop=True, inplace=True)
            print('Done!')
        if region=='EGY':
            #Screen for India polygons only
            print('\nScreening for EGY ONLY ...')
            bool_series = gdf_acct.geometry.apply(lambda g: g.within(egyp))
            gdf_acct = gdf_acct.loc[bool_series]
            gdf_acct.reset_index(drop=True, inplace=True)
            print('Done!')
        if region=='MEX':
            #Screen for Mexico polygons only
            print('\nScreening for MEX ONLY ...')
            bool_series = gdf_acct.geometry.apply(lambda g: g.within(mexp))
            gdf_acct = gdf_acct.loc[bool_series]
            gdf_acct.reset_index(drop=True, inplace=True)
            print('Done!')
        if region=='CA':
            #Screen for California polygons only
            print('\nScreening for CA ONLY ...')
            bool_series = gdf_acct.geometry.apply(lambda g: g.within(cap))
            gdf_acct = gdf_acct.loc[bool_series]
            gdf_acct.reset_index(drop=True, inplace=True)
            print('Done!')
        elif region=='AZ':
            #Screen for AZ polygons only
            print('\nScreening for AZ ONLY ...')
            bool_series = gdf_acct.geometry.apply(lambda g: g.within(azp))
            gdf_acct = gdf_acct.loc[bool_series]
            gdf_acct.reset_index(drop=True, inplace=True)
            print('Done!')
        elif region=='WA':
            #Screen for WA polygons only
            print('\nCScreening for WA ONLY ...')
            bool_series = gdf_acct.geometry.apply(lambda g: g.within(wap))
            gdf_acct = gdf_acct.loc[bool_series]
            gdf_acct.reset_index(drop=True, inplace=True)
            print('Done!')
        elif region=='TX':
            #Screen for WA polygons only
            print('\nScreening for TX ONLY ...')
            bool_series = gdf_acct.geometry.apply(lambda g: g.within(txp))
            gdf_acct = gdf_acct.loc[bool_series]
            gdf_acct.reset_index(drop=True, inplace=True)
            print('Done!')
        elif region=='FL':
            #Screen for WA polygons only
            print('\nScreening for FL ONLY ...')
            bool_series = gdf_acct.geometry.apply(lambda g: g.within(flp))
            gdf_acct = gdf_acct.loc[bool_series]
            gdf_acct.reset_index(drop=True, inplace=True)
            print('Done!')
        elif region=='AUS':
            #Screen for WA polygons only
            print('\nScreening for AUS ONLY ...')
            a_gdf = gpd.GeoDataFrame()
            acct_gdf = gpd.GeoDataFrame()
            for ausp in ausp_list:
                bool_series = gdf_acct.geometry.apply(lambda g: g.within(ausp))
                a_gdf = gdf_acct.loc[bool_series]
                acct_gdf = acct_gdf.append(a_gdf, ignore_index=True)
            gdf_acct = acct_gdf.copy()
            gdf_acct.reset_index(drop=True, inplace=True)
        elif region=='CAN':
            #Screen for CAN polygons only
            print('\nScreening for CAN ONLY ...')
            a_gdf = gpd.GeoDataFrame()
            acct_gdf = gpd.GeoDataFrame()
            for canp in canp_list:
                bool_series = gdf_acct.geometry.apply(lambda g: g.within(canp))
                a_gdf = gdf_acct.loc[bool_series]
                acct_gdf = acct_gdf.append(a_gdf, ignore_index=True)
            gdf_acct = acct_gdf.copy()
            gdf_acct.reset_index(drop=True, inplace=True)
            print('Done!')

    list_of_valid_accts = pd.unique(gdf_acct.Category).tolist()
    list_of_valid_accts.sort()


    return list_of_valid_accts

def getExistingAcctListForRegion(region, acctName):
    if acctName=='ALL':
        list_of_existing_acctFolders = [x for x in os.listdir(agls_acctConfigPath)]
    else:
        list_of_existing_acctFolders = [x for x in os.listdir(agls_acctConfigPath) if x.startswith(acctName)]
    
    #print(list_of_existing_acctFolders)
    list_of_existing_accts = []
    for f in list_of_existing_acctFolders:
        folder = agls_acctConfigPath + f + '/'
        configFile = folder + f+'.json'
        if os.path.exists(configFile):
            list_of_existing_accts.append(f)


    list_of_existing_accts.sort()
    #load the acct list into a dataframe
    #print('List of accts: '+str(list_of_existing_accts))
    gdf_acct = readInternalAcctsConfig(list_of_existing_accts)
    #print(gdf_acct)
    #print(gdf_acct.columns)

    #Screen for the region
    list_of_valid_regions = ['CA','WA','FL','TX','AZ','AUS','US', 'IND', 'MEX','EGY', 'CAN']
    if not region in list_of_valid_regions:
        print('Please provide specific proper region...aborting')
        return []
    else:
        if region=='CA':
            #Screen for California polygons only
            print('\nScreening for CA ONLY ...')
            bool_series = gdf_acct.geometry.apply(lambda g: g.within(cap))
            gdf_acct = gdf_acct.loc[bool_series]
            gdf_acct.reset_index(drop=True, inplace=True)
            print('Done!')
        elif region=='IND':
            #Screen for IND polygons only
            print('\nScreening for IND ONLY ...')
            bool_series = gdf_acct.geometry.apply(lambda g: g.within(indp))
            gdf_acct = gdf_acct.loc[bool_series]
            gdf_acct.reset_index(drop=True, inplace=True)
            print('Done!')
        elif region=='MEX':
            #Screen for MEX polygons only
            print('\nScreening for MEX ONLY ...')
            bool_series = gdf_acct.geometry.apply(lambda g: g.within(mexp))
            gdf_acct = gdf_acct.loc[bool_series]
            gdf_acct.reset_index(drop=True, inplace=True)
            print('Done!')
        elif region=='EGY':
            #Screen for MEX polygons only
            print('\nScreening for EGY ONLY ...')
            bool_series = gdf_acct.geometry.apply(lambda g: g.within(egyp))
            gdf_acct = gdf_acct.loc[bool_series]
            gdf_acct.reset_index(drop=True, inplace=True)
            print('Done!')
        elif region=='AZ':
            #Screen for AZ polygons only
            print('\nScreening for AZ ONLY ...')
            bool_series = gdf_acct.geometry.apply(lambda g: g.within(azp))
            gdf_acct = gdf_acct.loc[bool_series]
            gdf_acct.reset_index(drop=True, inplace=True)
            print('Done!')
        elif region=='WA':
            #Screen for WA polygons only
            print('\nScreening for WA ONLY ...')
            bool_series = gdf_acct.geometry.apply(lambda g: g.within(wap))
            gdf_acct = gdf_acct.loc[bool_series]
            gdf_acct.reset_index(drop=True, inplace=True)
            print('Done!')
        elif region=='TX':
            #Screen for WA polygons only
            print('\nScreening for TX ONLY ...')
            bool_series = gdf_acct.geometry.apply(lambda g: g.within(txp))
            gdf_acct = gdf_acct.loc[bool_series]
            gdf_acct.reset_index(drop=True, inplace=True)
            print('Done!')
        elif region=='FL':
            #Screen for WA polygons only
            print('\nScreening for FL ONLY ...')
            bool_series = gdf_acct.geometry.apply(lambda g: g.within(flp))
            gdf_acct = gdf_acct.loc[bool_series]
            gdf_acct.reset_index(drop=True, inplace=True)
            print('Done!')
        elif region=='OR':
            #Screen for WA polygons only
            print('\nScreening for OR ONLY ...')
            bool_series = gdf_acct.geometry.apply(lambda g: g.within(orp))
            gdf_acct = gdf_acct.loc[bool_series]
            gdf_acct.reset_index(drop=True, inplace=True)
            print('Done!')
        elif region=='OK':
            #Screen for WA polygons only
            print('\nScreening for OK ONLY ...')
            bool_series = gdf_acct.geometry.apply(lambda g: g.within(okp))
            gdf_acct = gdf_acct.loc[bool_series]
            gdf_acct.reset_index(drop=True, inplace=True)
            print('Done!')
        elif region=='MT':
            #Screen for WA polygons only
            print('\nScreening for MT ONLY ...')
            bool_series = gdf_acct.geometry.apply(lambda g: g.within(mtp))
            gdf_acct = gdf_acct.loc[bool_series]
            gdf_acct.reset_index(drop=True, inplace=True)
            print('Done!')
        elif region=='MO':
            #Screen for WA polygons only
            print('\nScreening for MO ONLY ...')
            bool_series = gdf_acct.geometry.apply(lambda g: g.within(mop))
            gdf_acct = gdf_acct.loc[bool_series]
            gdf_acct.reset_index(drop=True, inplace=True)
            print('Done!')
        elif region=='DE':
            #Screen for WA polygons only
            print('\nScreening for DE ONLY ...')
            bool_series = gdf_acct.geometry.apply(lambda g: g.within(dep))
            gdf_acct = gdf_acct.loc[bool_series]
            gdf_acct.reset_index(drop=True, inplace=True)
            print('Done!')
        elif region=='AR':
            #Screen for WA polygons only
            print('\nScreening for AR ONLY ...')
            bool_series = gdf_acct.geometry.apply(lambda g: g.within(arp))
            gdf_acct = gdf_acct.loc[bool_series]
            gdf_acct.reset_index(drop=True, inplace=True)
            print('Done!')
        elif region=='AUS':
            #Screen for WA polygons only
            print('\nScreening for AUS ONLY ...')
            a_gdf = gpd.GeoDataFrame()
            acct_gdf = gpd.GeoDataFrame()
            for ausp in ausp_list:
                bool_series = gdf_acct.geometry.apply(lambda g: g.within(ausp))
                a_gdf = gdf_acct.loc[bool_series]
                acct_gdf = acct_gdf.append(a_gdf, ignore_index=True)
            gdf_acct = acct_gdf.copy()
            gdf_acct.reset_index(drop=True, inplace=True)
        elif region=='CAN':
            #Screen for CAN polygons only
            print('\nScreening for CAN ONLY ...')
            a_gdf = gpd.GeoDataFrame()
            acct_gdf = gpd.GeoDataFrame()
            for canp in canp_list:
                bool_series = gdf_acct.geometry.apply(lambda g: g.within(canp))
                a_gdf = gdf_acct.loc[bool_series]
                acct_gdf = acct_gdf.append(a_gdf, ignore_index=True)
            gdf_acct = acct_gdf.copy()
            gdf_acct.reset_index(drop=True, inplace=True)
            print('Done!')

    list_of_valid_accts = pd.unique(gdf_acct.Category).tolist()
    #print(list_of_valid_accts)

    return list_of_valid_accts

def getExistingAGLSAcctListForRegion(region):
    #list_of_existing_acctIDs = [x.split('_')[-1] for x in os.listdir(acctConfigPath) if (x.startswith('JainlogicCustomerid_') & ~(x.find('__')==-1))]
    list_of_existing_acctFolders = [x.split('_')[-1] for x in os.listdir(acctConfigPath) if x.startswith('JainlogicCustomerid_')]
    list_of_existing_acctIDs = []
    for f in list_of_existing_acctFolders:
        folder = acctConfigPath + 'JainlogicCustomerid_'+f
        configFile = folder+ '/JainlogicCustomerid_'+f+'.json'
        if os.path.exists(configFile):
            list_of_existing_acctIDs.append(f)


    list_of_existing_acctIDs.sort(key=int)
    list_of_existing_acctIDs_int = list(map(int, list_of_existing_acctIDs))
    list_of_existing_accts = ['JainlogicCustomerid_'+str(x) for x in list_of_existing_acctIDs_int]

    #load the acct list into a dataframe
    gdf_acct = readUpdatedAcctsConfig(list_of_existing_acctIDs_int)

    #Screen for the region
    list_of_valid_regions = ['CA','WA','FL','TX','AZ','OR','OK','MT','MO','DE','AR','AUS','US','IND','MEX','CAM']
    if not region in list_of_valid_regions:
        print('Please provide specific proper region...aborting')
        return []
    else:
        if region=='CA':
            #Screen for California polygons only
            print('\nScreening for CA ONLY ...')
            bool_series = gdf_acct.geometry.apply(lambda g: g.within(cap))
            gdf_acct = gdf_acct.loc[bool_series]
            gdf_acct.reset_index(drop=True, inplace=True)
            print('Done!')
        elif region=='AZ':
            #Screen for AZ polygons only
            print('\nScreening for AZ ONLY ...')
            bool_series = gdf_acct.geometry.apply(lambda g: g.within(azp))
            gdf_acct = gdf_acct.loc[bool_series]
            gdf_acct.reset_index(drop=True, inplace=True)
            print('Done!')
        elif region=='WA':
            #Screen for WA polygons only
            print('\nCScreening for WA ONLY ...')
            bool_series = gdf_acct.geometry.apply(lambda g: g.within(wap))
            gdf_acct = gdf_acct.loc[bool_series]
            gdf_acct.reset_index(drop=True, inplace=True)
            print('Done!')
        elif region=='TX':
            #Screen for WA polygons only
            print('\nScreening for TX ONLY ...')
            bool_series = gdf_acct.geometry.apply(lambda g: g.within(txp))
            gdf_acct = gdf_acct.loc[bool_series]
            gdf_acct.reset_index(drop=True, inplace=True)
            print('Done!')
        elif region=='FL':
            #Screen for WA polygons only
            print('\nScreening for FL ONLY ...')
            bool_series = gdf_acct.geometry.apply(lambda g: g.within(flp))
            gdf_acct = gdf_acct.loc[bool_series]
            gdf_acct.reset_index(drop=True, inplace=True)
            print('Done!')
        elif region=='OR':
            print('\nScreening for OR ONLY ...')
            bool_series = gdf_acct.geometry.apply(lambda g: g.within(orp))
            gdf_acct = gdf_acct.loc[bool_series]
            gdf_acct.reset_index(drop=True, inplace=True)
            print('Done!')
        elif region=='OK':
            print('\nScreening for OK ONLY ...')
            bool_series = gdf_acct.geometry.apply(lambda g: g.within(okp))
            gdf_acct = gdf_acct.loc[bool_series]
            gdf_acct.reset_index(drop=True, inplace=True)
            print('Done!')
        elif region=='MT':
            print('\nScreening for MT ONLY ...')
            bool_series = gdf_acct.geometry.apply(lambda g: g.within(mtp))
            gdf_acct = gdf_acct.loc[bool_series]
            gdf_acct.reset_index(drop=True, inplace=True)
            print('Done!')
        elif region=='MO':
            print('\nScreening for MO ONLY ...')
            bool_series = gdf_acct.geometry.apply(lambda g: g.within(mop))
            gdf_acct = gdf_acct.loc[bool_series]
            gdf_acct.reset_index(drop=True, inplace=True)
            print('Done!')
        elif region=='DE':
            print('\nScreening for DE ONLY ...')
            bool_series = gdf_acct.geometry.apply(lambda g: g.within(dep))
            gdf_acct = gdf_acct.loc[bool_series]
            gdf_acct.reset_index(drop=True, inplace=True)
            print('Done!')
        elif region=='AR':
            print('\nScreening for AR ONLY ...')
            bool_series = gdf_acct.geometry.apply(lambda g: g.within(arp))
            gdf_acct = gdf_acct.loc[bool_series]
            gdf_acct.reset_index(drop=True, inplace=True)
            print('Done!')
        elif region=='AUS':
            #Screen for WA polygons only
            print('\nScreening for AUS ONLY ...')
            a_gdf = gpd.GeoDataFrame()
            acct_gdf = gpd.GeoDataFrame()
            for ausp in ausp_list:
                bool_series = gdf_acct.geometry.apply(lambda g: g.within(ausp))
                a_gdf = gdf_acct.loc[bool_series]
                acct_gdf = acct_gdf.append(a_gdf, ignore_index=True)
            gdf_acct = acct_gdf.copy()
            gdf_acct.reset_index(drop=True, inplace=True)
        elif region=='IND':
            #Screen for IND polygons only
            print('\nScreening for IND ONLY ...')
            a_gdf = gpd.GeoDataFrame()
            acct_gdf = gpd.GeoDataFrame()
            for indp in indp_list:
                bool_series = gdf_acct.geometry.apply(lambda g: g.within(indp))
                a_gdf = gdf_acct.loc[bool_series]
                acct_gdf = acct_gdf.append(a_gdf, ignore_index=True)
            gdf_acct = acct_gdf.copy()
            gdf_acct.reset_index(drop=True, inplace=True)
        elif region=='MEX':
            #Screen for MEX polygons only
            print('\nScreening for MEX ONLY ...')
            a_gdf = gpd.GeoDataFrame()
            acct_gdf = gpd.GeoDataFrame()
            for mexp in mexp_list:
                bool_series = gdf_acct.geometry.apply(lambda g: g.within(mexp))
                a_gdf = gdf_acct.loc[bool_series]
                acct_gdf = acct_gdf.append(a_gdf, ignore_index=True)
            gdf_acct = acct_gdf.copy()
            gdf_acct.reset_index(drop=True, inplace=True)
        elif region=='CAN':
            #Screen for CAN polygons only
            print('\nScreening for CAN ONLY ...')
            a_gdf = gpd.GeoDataFrame()
            acct_gdf = gpd.GeoDataFrame()
            for canp in canp_list:
                bool_series = gdf_acct.geometry.apply(lambda g: g.within(canp))
                a_gdf = gdf_acct.loc[bool_series]
                acct_gdf = acct_gdf.append(a_gdf, ignore_index=True)
            gdf_acct = acct_gdf.copy()
            gdf_acct.reset_index(drop=True, inplace=True)
            print('Done!')

    list_of_valid_acctIDs = pd.unique(gdf_acct.customerid).tolist()
    #print(list_of_valid_acctIDs)

    list_of_valid_acctIDs_int = [int(x) for x in list_of_valid_acctIDs]
    list_of_valid_acctIDs_int.sort()

    #print(list_of_valid_acctIDs_int)

    return list_of_valid_acctIDs_int


def getExistingAGLSAcctsGDF():
    list_of_existing_acctIDs = [x.split('_')[-1] for x in os.listdir(acctConfigPath) if x.startswith('JainlogicCustomerid_')]
    list_of_existing_acctIDs.sort(key=int)
    list_of_existing_acctIDs_int = list(map(int, list_of_existing_acctIDs))
    list_of_existing_accts = ['JainlogicCustomerid_'+str(x) for x in list_of_existing_acctIDs_int]

    gdf_existingAccts = gpd.GeoDataFrame()
    for acct in list_of_existing_accts:
        configJsonFile = acctConfigPath + acct + '/'+acct +'.json'
        agdf = gpd.GeoDataFrame()

        #open and read this file if it exists
        if os.path.exists(configJsonFile):
            with open(configJsonFile, "r") as read_file:
                data = json.load(read_file)
            try:
                g = gpd.GeoDataFrame.from_features(data["features"])
                agdf = agdf.append(g,ignore_index=True)
            except:
                json_data = json.loads(data)
                g = gpd.GeoDataFrame.from_features(json_data["features"])
                agdf = agdf.append(g,ignore_index=True) 
            #Put the centroid of each field
            for i,row in agdf.iterrows():
                c = row.geometry.centroid
                agdf.loc[i,'CentroidPoint']=c.to_wkt()
            agdf.fieldassetid = agdf.fieldassetid.astype('int')
            gdf_existingAccts = gdf_existingAccts.append(agdf, ignore_index=True)

    return gdf_existingAccts



def getExistingAGLSAcctList():
    list_of_existing_acctIDs = [x.split('_')[-1] for x in os.listdir(acctConfigPath) if x.startswith('JainlogicCustomerid_')]
    list_of_existing_acctFolders = [x.split('_')[-1] for x in os.listdir(acctConfigPath) if x.startswith('JainlogicCustomerid_')]
    list_of_existing_acctIDs = []
    for f in list_of_existing_acctFolders:
        folder = acctConfigPath + 'JainlogicCustomerid_'+f
        configFile = folder+ '/JainlogicCustomerid_'+f+'.json'
        if os.path.exists(configFile):
            list_of_existing_acctIDs.append(f)

    list_of_existing_acctIDs.sort(key=int)
    list_of_existing_acctIDs_int = list(map(int, list_of_existing_acctIDs))
    list_of_existing_accts = ['JainlogicCustomerid_'+str(x) for x in list_of_existing_acctIDs_int]

    return list_of_existing_acctIDs_int

def createNewAGLSAccts(list_of_new_acctIDs_to_be_created, gdf_acct):
    print('\nCreating new accts for custmerids: '+str(list_of_new_acctIDs_to_be_created))
    #list_of_new_acctIDs_to_be_created
    gdf = gdf_acct[gdf_acct['customerid'].isin(list_of_new_acctIDs_to_be_created)]

    #get the satellite data
    gdf_newAccounts = getSatTilesInfo(gdf) #get the Satellite info populated

    #update the config files
    for acctID in list_of_new_acctIDs_to_be_created:
        acct = 'JainlogicCustomerid_'+ str(acctID)
        print('\tCreating acct config for acct: '+acct)
        html_ranch_map_dir = json_output_dir = outDir + acct + '/'

        #if the directory doesn't exist then make it
        if not os.path.exists(html_ranch_map_dir):
            os.makedirs(html_ranch_map_dir)

        #get the subset for the acct
        df = gdf_newAccounts[gdf_newAccounts.customerid == acctID ]
        df.reset_index(drop=True,inplace=True)

        #generate a ranch map and write the JSON
        first_polygon = df.loc[0,'geometry']
        c = [list(first_polygon.centroid.coords)[0][1],list(first_polygon.centroid.coords)[0][0]]
        m = folium.Map(location=c, zoom_start=10, tiles = 'stamentoner')

        for j in range(len(df)):
            gs = folium.GeoJson(df.iloc[j:j+1])
            p = df.loc[j,'geometry']
            ac = df.loc[j,'Acres']
            label = 'Name: '+df.loc[j,'Name']+', Size(Acres): '+str(round(ac,2))
            folium.Popup(label).add_to(gs)
            gs.add_to(m)
            uid = str(uuid.uuid4())
            df.loc[j,'Category']='Jainlogic Monitored'
            df.loc[j,'UUID']=uid
        html_file = html_ranch_map_dir + acct+'.html'
        #print(html_file)
        m.save(html_file)

        json_output_file = json_output_dir + acct + '.json'
        print(json_output_file)
        df.to_file(json_output_file, driver='GeoJSON')
        print('\tDone creating acct config!')

    return gdf_newAccounts

def updateExistingAGLSAcctsForNameChange(changedNames_gdf):
    gdf_updateAccts = changedNames_gdf.copy()
    list_of_acctIDs_to_be_updated = pd.unique(gdf_updateAccts.customerid).tolist()
    list_of_acctIDs_to_be_updated.sort(key=int)

    for a in list_of_acctIDs_to_be_updated:
        gdf_updateAccts.customerid=gdf_updateAccts.customerid.astype('str')
        acct_gdf = gdf_updateAccts[gdf_updateAccts.customerid==str(a)]
        acct_gdf.reset_index(drop=True, inplace=True)
        acct = 'JainlogicCustomerid_'+str(a)
        
        #generate a NEW config file
        json_output_dir = outDir + acct + '/'
        json_output_file = json_output_dir + acct + '.json'
        
        print('\tCreating New JSON: ')
        print(json_output_file)
        #Make sure the IDs are strings
        acct_gdf[['customerid','fieldassetid']] = acct_gdf[['customerid','fieldassetid']].astype('str')

        #Make sure that the footprints are strings
        acct_gdf[['SentinelFootprint','LandsatFootprint']] = acct_gdf[['SentinelFootprint','LandsatFootprint']].astype('str')

        acct_gdf.to_file(json_output_file, driver='GeoJSON')



    return

def updateExistingAGLSAccts(list_of_acctIDs_to_be_updated, gdf_acct):
    gdf_existingAccts = gpd.GeoDataFrame()
    gdf_updatedAccts = gpd.GeoDataFrame()
    gdf_updatedIncrAccts = gpd.GeoDataFrame()

    if len(list_of_acctIDs_to_be_updated)==0: #empty list
        gdf_existingAccts = gdf_acct.copy()
    else:
        gdf = gdf_acct[gdf_acct['customerid'].isin(list_of_acctIDs_to_be_updated)]
        gdf.reset_index(drop=True,inplace=True)
        gdf_newAssets = gdf.copy()

        for acctID in list_of_acctIDs_to_be_updated:

            acct = 'JainlogicCustomerid_'+str(acctID)
            print('\n\n*******************************************\nAcct: '+acct)

            #read the existing acct_gdf for this acctid
            configJsonFile = acctConfigPath + acct + '/'+acct +'.json'
            agdf = gpd.GeoDataFrame()

            if os.path.exists(configJsonFile):
                with open(configJsonFile, "r") as read_file:
                    data = json.load(read_file)
                try:
                    g = gpd.GeoDataFrame.from_features(data["features"])
                    agdf = agdf.append(g,ignore_index=True)
                except:
                    json_data = json.loads(data)
                    g = gpd.GeoDataFrame.from_features(json_data["features"])
                    agdf = agdf.append(g,ignore_index=True) 
                #Put the centroid of each field
                for i,row in agdf.iterrows():
                    c = row.geometry.centroid
                    agdf.loc[i,'CentroidPoint']=c.to_wkt()
                agdf.fieldassetid = agdf.fieldassetid.astype('int')
                gdf_existingAccts = gdf_existingAccts.append(agdf, ignore_index=True)

            print('\n\tExisting Acct Info: ')
            for j,row in agdf.iterrows():
                print(['CustomerID: '+str(row.customerid)+', AssetID: '+str(row.fieldassetid)])



            df = gdf_newAssets[gdf_newAssets.customerid == acctID ]
            df.reset_index(drop=True,inplace=True)
            dff = gpd.GeoDataFrame()
            k=0
            print('\n\tNew Info for Account to be Updated :')
            for j,row in df.iterrows():
                #PRINT HERE ONLY IF
                #a) IT IS A NEW ASSET **OR**
                #b) THE EXISTING POLYGON IS EXACTLY THE SAME (HAS BEEN UPDATED BEFORE)
                poly_new = row.geometry
                fid = row.fieldassetid
                poly_orig_df = agdf[agdf.fieldassetid == fid]
                poly_orig_df.reset_index(drop=True, inplace=True)

                new_geom=False
                if len(poly_orig_df) > 0:
                    new_poly=False
                    poly_orig = poly_orig_df.geometry[0]

                    if poly_new.equals(poly_orig):
                        new_geom = False
                    else:
                        new_geom = True
                else:
                    new_poly=True

                if ((new_poly) | (new_geom)):
                    print(['CustomerID: '+str(row.customerid)+', AssetID: '+str(row.fieldassetid)])
                    dff =  dff.append(row, ignore_index=True)

            concat_gdf = gpd.GeoDataFrame()
            df_new = gpd.GeoDataFrame()
            if len(dff)>0:
                for k,row in dff.iterrows():
                    uid = str(uuid.uuid4())
                    c = row.geometry.centroid
                    dff.loc[k,'Category']='Jainlogic Monitored'
                    dff.loc[k,'UUID']=uid
                    dff.loc[k,'AttributesJSON']=makeJsonStr('Acres',row.Acres)
                    dff.loc[k,'CentroidPoint']=c.to_wkt()
                #get the Sentinel and Landsat Footprint
                df_new = getSatTilesInfo(dff)

                #Drop the rows with this list of updated fieldassetids

                #concat the two datasets while dropping the redundanct rows from the original
                concat_gdf = pd.concat([agdf,df_new], axis=0, ignore_index=True, sort=False).drop_duplicates(['fieldassetid'], keep='last')
                concat_gdf.customerid = concat_gdf.customerid.astype('int')
                concat_gdf.fieldassetid = concat_gdf.fieldassetid.astype('int')
                concat_gdf = concat_gdf.sort_values(by=['fieldassetid'])
                concat_gdf.reset_index(drop=True, inplace=True)

                print('\n\tFinal Acct Config: ')
                for j,row in concat_gdf.iterrows():
                    print(['CustomerID: '+str(row.customerid)+', AssetID: '+str(row.fieldassetid)])


                #generate a NEW ranch map
                html_ranch_map_dir = json_output_dir = outDir + acct + '/'
                first_polygon = concat_gdf.loc[0,'geometry']
                c = [list(first_polygon.centroid.coords)[0][1],list(first_polygon.centroid.coords)[0][0]]
                m = folium.Map(location=c, zoom_start=10, tiles = 'stamentoner')
                html_file = html_ranch_map_dir + acct+'.html'
                #print(html_file)
                m.save(html_file)

                #generate a NEW config file
                json_output_file = json_output_dir + acct + '.json'
                print('\tCreating New JSON: ')
                print(json_output_file)
                #Make sure the IDs are strings
                concat_gdf[['customerid','fieldassetid']] = concat_gdf[['customerid','fieldassetid']].astype('str')

                #Make sure that the footprints are strings
                concat_gdf[['SentinelFootprint','LandsatFootprint']] = concat_gdf[['SentinelFootprint','LandsatFootprint']].astype('str')

                concat_gdf.to_file(json_output_file, driver='GeoJSON')
            else:
                print('NO CHANGES MADE')
            print('*******************************************\n')
            print('*******************************************\n')

            #Update the changes to both the updated accounts and added accounts
            if len(concat_gdf)>0:
                gdf_updatedAccts = gdf_updatedAccts.append(concat_gdf, ignore_index=True)
            if len(df_new)>0:
                df_new.customerid = df_new.customerid.astype('int')
                df_new.fieldassetid = df_new.fieldassetid.astype('int')
                df_new[['customerid','fieldassetid']] = df_new[['customerid','fieldassetid']].astype('str')
                gdf_updatedIncrAccts = gdf_updatedIncrAccts.append(df_new, ignore_index=True)

    return gdf_existingAccts, gdf_updatedAccts, gdf_updatedIncrAccts

def readAcctConfig(acctType, list_of_accts):
    acct_gdf = gpd.GeoDataFrame()
    #Now get the existing accounts for these accountIDs from their config files
    for acctStr in list_of_accts:
        print('\t\t\tAccount: '+str(acctStr))
        try:
            acct = str(acctStr)
            if acctType=='PILOT':
                configJsonFile =  pilotAcctConfigPath + acct + '/'+acct +'.json'
            elif acctType=='AGLS_ACCOUNTS':
                configJsonFile =  agls_acctConfigPath + acct + '/'+acct +'.json'
            else:
                print('Illegal Acct Type')
                exit()

            with open(configJsonFile, "r") as read_file:
                data = json.load(read_file)
            try:
                g = gpd.GeoDataFrame.from_features(data["features"])
                acct_gdf = acct_gdf.append(g,ignore_index=True)
            except:
                json_data = json.loads(data)
                g = gpd.GeoDataFrame.from_features(json_data["features"])
                acct_gdf = acct_gdf.append(g,ignore_index=True) 
        except:
            print('Exception in reading file for acct: '+str(acctID))
            continue	
    #Put the centroid of each field
    for i,row in acct_gdf.iterrows():
        c = row.geometry.centroid
        acct_gdf.loc[i,'CentroidPoint']=c.to_wkt()

    acct_gdf.reset_index(drop=True,inplace=True)

    return acct_gdf

def readInternalAcctsConfig(list_of_new_accts):
    acct_gdf = gpd.GeoDataFrame()
    #Now get the existing accounts for these accountIDs from their config files
    for a in list_of_new_accts:
        #print('\t\t\tAcct: '+str(acct))
        acct = a.upper()
        try:
            configJsonFile =  agls_acctConfigPath + acct+ '/' + acct +'.json'

            with open(configJsonFile, "r") as read_file:
                data = json.load(read_file)
            try:
                g = gpd.GeoDataFrame.from_features(data["features"])
                acct_gdf = acct_gdf.append(g,ignore_index=True)
            except:
                json_data = json.loads(data)
                g = gpd.GeoDataFrame.from_features(json_data["features"])
                acct_gdf = acct_gdf.append(g,ignore_index=True) 
        except:
            print('Exception in reading file for acct: '+str(acct))
            continue	
    #Put the centroid of each field
    for i,row in acct_gdf.iterrows():
        c = row.geometry.centroid
        acct_gdf.loc[i,'CentroidPoint']=c.to_wkt()

    acct_gdf.reset_index(drop=True,inplace=True)

    return acct_gdf

def readUpdatedAcctsConfig(list_of_new_acctIDs_to_be_updated):
    acct_gdf = gpd.GeoDataFrame()
    #Now get the existing accounts for these accountIDs from their config files
    for acctID in list_of_new_acctIDs_to_be_updated:
        print('\t\t\tcustomerid: '+str(acctID))
        try:
            acct = 'JainlogicCustomerid_'+str(acctID)
            configJsonFile =  acctConfigPath + acct + '/'+acct +'.json'

            with open(configJsonFile, "r") as read_file:
                data = json.load(read_file)
            try:
                g = gpd.GeoDataFrame.from_features(data["features"])
                acct_gdf = acct_gdf.append(g,ignore_index=True)
            except:
                json_data = json.loads(data)
                g = gpd.GeoDataFrame.from_features(json_data["features"])
                acct_gdf = acct_gdf.append(g,ignore_index=True) 
        except:
            print('Exception in reading file for acct: '+str(acctID))
            continue	
    #Put the centroid of each field
    for i,row in acct_gdf.iterrows():
        c = row.geometry.centroid
        acct_gdf.loc[i,'CentroidPoint']=c.to_wkt()

    acct_gdf.reset_index(drop=True,inplace=True)

    return acct_gdf

def readPilotAcctsConfig(list_of_new_accts):
    acct_gdf = gpd.GeoDataFrame()
    #Now get the existing accounts for these accountIDs from their config files
    for acctID in list_of_new_accts:
        print('\t\t\tPilot Acct: '+str(acctID))
        try:
            acct = str(acctID)
            configJsonFile =  pilotAcctConfigPath + acct + '/'+acct +'.json'

            with open(configJsonFile, "r") as read_file:
                data = json.load(read_file)
            try:
                g = gpd.GeoDataFrame.from_features(data["features"])
                acct_gdf = acct_gdf.append(g,ignore_index=True)
            except:
                json_data = json.loads(data)
                g = gpd.GeoDataFrame.from_features(json_data["features"])
                acct_gdf = acct_gdf.append(g,ignore_index=True) 
        except:
            print('Exception in reading file for acct: '+str(acctID))
            continue	
    #Put the centroid of each field
    for i,row in acct_gdf.iterrows():
        c = row.geometry.centroid
        acct_gdf.loc[i,'CentroidPoint']=c.to_wkt()

    acct_gdf.reset_index(drop=True,inplace=True)

    return acct_gdf


def writeCIMISDataToFile(configFile, startDtStr, endDtStr, stns_gdf):
    #Step 1: get the weather gdf ready
    download_folder = cimisRawPath
    dstFolder = download_folder + 'monthlyAllStns/'
    historic_cimis_loc = cimis_monthly_ftp_loc

    #get the daily weather 
    dailyCIMISdata = cimisProcPath + 'dailyWeather.xlsx'
    #read into a pandas dataframe
    dailyWeather_df = pd.read_excel(dailyCIMISdata)

    dailyCols = ['StationId',
                 'Date',
                 'JulianDate',
                 'Reference ETo(in)',
                 'QCforReference ETo',
                 'Precipitation(in)',
                 'QCforPrecipitation',
                 'SolarRadiationAverage(Ly/day)',
                 'QCforSolarRadiationAverage',
                 'AverageVaporPressure(mBars)',
                 'QCforAverageVaporPressure',
                 'Maximum Air Temperature(F)',
                 'QC for Maximum Air Temperature',
                 'Minimum Air Temperature   (F)',
                 'QC for Minimum Air Temperature',
                 'Average Air Temperature   (F)',
                 'QC for Average Air Temperature',
                 'Maximum Relative Humidity   (%)',
                 'QC for Maximum Relative Humidity',
                 'Minimum Relative Humidity   (%)',
                 'QC for Minimum Relative Humidity',
                 'Average Relative Humidity   (%)',
                 'QC for Average Relative Humidity',
                 'Dew Point   (∞=F)',
                 'QC for Dew Point',
                 'Average Wind Speed   (mph)',
                 'QC for Average Wind Speed',
                 'Wind Run   (miles)',
                 'QC for Wind Run',
                 'Average Soil Temperature   (F)',
                 'QC for Average Soil Temperature'
                ]
    dailyCols = [x.replace(' ','') for x in dailyCols]
    dailyCols = [x.replace('∞','') for x in dailyCols]
    dailyCols = [x.replace('=','') for x in dailyCols]


    dailyMonthlyWeatherPD = pd.DataFrame(columns=dailyCols)
    data = pd.DataFrame(columns = dailyCols)
    #read all the files into a big pandas dataframe
    allfiles = os.listdir(dstFolder)
    monthlyStnFiles = [fname for fname in allfiles if fname.endswith('.csv')]
    for fl in monthlyStnFiles:
        f = dstFolder + fl
        #print(f)
        mon = fl[0:3]
        stnNum = f.split('.')[0][-3:]
        if int(stnNum) != 259:
            if (int(stnNum) in stns_gdf.StationNumber.values): #if the stnnum is an active station
                #read each file
                data = pd.read_csv(f,header=None)
                data.columns = dailyCols
                #fix the Date to make it a Datetime obj
                data['Date'] =  pd.to_datetime(data['Date'],format='%m/%d/%Y')
                #dailyMonthlyWeatherPD = pd.concat([dailyMonthlyWeatherPD, data], ignore_index=True, sort=True)
                dailyMonthlyWeatherPD = dailyMonthlyWeatherPD.append(data, ignore_index=True)

    #Concatenate the two dataframes and filter on time
    dailyMonthlyWeatherPD = dailyMonthlyWeatherPD[dailyWeather_df.columns]
    dailyMonthlyWeatherPD_all = pd.concat([dailyMonthlyWeatherPD,dailyWeather_df], axis=0, ignore_index=True)
    dailyMonthlyWeatherPD_all.sort_values(by=['Date'], inplace=True)
    
    #startDtStr and endDtStr are defined earlier
    num_days_run = datetime.strptime(endDtStr,'%Y-%m-%d') - datetime.strptime(startDtStr,'%Y-%m-%d')

    dailyMonthlyWeatherPD = dailyMonthlyWeatherPD_all[(dailyMonthlyWeatherPD_all.Date >=startDtStr) & (dailyMonthlyWeatherPD_all.Date <= endDtStr)]
    dailyMonthlyWeatherPD.reset_index(drop=True,inplace=True) 

    """
    #startDtStr and endDtStr are defined earlier
    num_days_run = datetime.strptime(endDtStr,'%Y-%m-%d') - datetime.strptime(startDtStr,'%Y-%m-%d')

    dailyMonthlyWeatherPD = dailyWeather_df[(dailyWeather_df.Date >=startDtStr) & (dailyWeather_df.Date <= endDtStr)]
    dailyMonthlyWeatherPD.reset_index(drop=True,inplace=True)  
    """

    weather2020_df = dailyMonthlyWeatherPD.copy()
    new_weather_df = interpolateMissingVals(weather2020_df)
    nd = (num_days_run.days)
    complete_weather_df = new_weather_df.copy()
    #removeIncompleteStns(new_weather_df,nd)
    #Sort by Julian Date
    weather2020_df = complete_weather_df.iloc[complete_weather_df['JulianDate'].astype(int).argsort()]
    dmw_pd = weather2020_df.copy()

    #Screen out the null values
    dmw_pd['ReferenceETo(in)'] = dmw_pd['ReferenceETo(in)'].replace('--',np.nan, regex=True).astype('float')
    dmw_pd['AverageAirTemperature(F)'] = dmw_pd['AverageAirTemperature(F)'].replace('--',np.nan, regex=True).astype('float')
    dmw_pd['SolarRadiationAverage(Ly/day)'] = dmw_pd['SolarRadiationAverage(Ly/day)'].replace('--',np.nan, regex=True).astype('float')
    dmw_pd['AverageWindSpeed(mph)'] = dmw_pd['AverageWindSpeed(mph)'].replace('--',np.nan, regex=True).astype('float')

    #Put values for metric system calculations
    dmw_pd['Tavg']=(5/9)*(dmw_pd['AverageAirTemperature(F)'].values.astype('float')-32) #convert to DegC
    dmw_pd['Rnet']=dmw_pd['SolarRadiationAverage(Ly/day)'].values.astype('float') * 0.484583 #Convert to Wm-2
    dmw_pd['Wavg']=dmw_pd['AverageWindSpeed(mph)'].values.astype('float') * 0.44704 #Convert to ms-1


    keep_eto_cols_list = ['Date', 
                          'StationId', 
                          'ReferenceETo(in)',
                          'Precipitation(in)', 
                          'AverageSoilTemperature(F)',
                          'MaximumAirTemperature(F)',
                          'MinimumAirTemperature(F)',
                          'AverageAirTemperature(F)',
                          'MaximumRelativeHumidity(%)',
                          'MinimumRelativeHumidity(%)',
                          'AverageRelativeHumidity(%)',
                          'DewPoint(F)',
                          'WindRun(miles)'
                         ]
    eto = dmw_pd[keep_eto_cols_list]

    #print(eto)

    keep_metric_cols_list = ['Date','StationId','Tavg','Rnet','Wavg']
    metric = dmw_pd[keep_metric_cols_list] 

    #for each station, we have to interpolate the timeseries and fill in the values
    #first for ETo
    new_list_of_active_stns = []
    list_of_stns = pd.unique(eto.StationId)

    #print(list_of_stns)

    eto_interp = pd.DataFrame()
    for stn in list_of_stns:
        eto_stn = eto[eto.StationId==stn]
        eto_stn = eto_stn.set_index(['Date'])
        eto_stn = eto_stn.sort_index()
        n = len(eto_stn[eto_stn['ReferenceETo(in)'].isnull()==True])
        if n==0:
            eto_interp = eto_interp.append(eto_stn)
        elif (n<100):
            print('\tInterpolating for '+str(n)+' points for stn '+ str(stn))
            new_list_of_active_stns.append(stn)
            eto_stn_interpolate = eto_stn.interpolate(method='time')
            eto_interp = eto_interp.append(eto_stn_interpolate)
        else:
            print([n,str(stn)])

    stns_df = stns_gdf.copy()
    stns_df = stns_df.rename(columns={"StationNumber": "StationId"})
    eto_interp = eto_interp.reset_index()
    eto_interp.StationId = eto_interp.StationId.astype('int64')
    df = pd.merge(eto_interp,stns_df,on='StationId')
    #df_metric = pd.merge(metric,df,on='StationId')
    df_metric = pd.merge(df,metric,how='left', left_on=['StationId','Date'], right_on = ['StationId','Date'])
    #convert into a GeoDataFrame
    geom = [Point(xy) for xy in zip(df.Longitude, df.Latitude)]
    crs = {'init': 'epsg:4326'}
    df_metric.StationId = df_metric.StationId.astype('int64')
    gdf_metric = gpd.GeoDataFrame(df_metric, crs=crs, geometry=geom)

    print('\tWriting CIMIS data to json file: '+configFile)
    gdf_metric.to_file(configFile, driver="GeoJSON")
    print('\tDone!')
    
    return


def getCIMISStns():
    #Read the list of active stations
    activeStnsFile = cimisProcPath + 'ActiveStations.xlsx'
    stns_df = pd.read_excel(activeStnsFile)
    for i,row in stns_df.iterrows():
        stns_df.loc[i,'geometry']=Point([row.Longitude, row.Latitude, row.ELEV])
    stns_gdf =  gpd.GeoDataFrame(stns_df,geometry=stns_df.geometry)
    stnPoints = list(pd.unique(stns_gdf.geometry))
    stns = pd.unique(stns_gdf.StationNumber)

    return stns_gdf, stnPoints, stns

def generateAndMergeETo(acct_gdf,Weather2020_gdf,startDtStr, endDtStr):
    start = datetime.strptime(startDtStr,'%Y-%m-%d')
    end = datetime.strptime(endDtStr,'%Y-%m-%d')
    #-timedelta(days=1) #This date is NOT included

    list_of_dates = [start + timedelta(days=x+1) for x in range(0, (end-start).days)]

    startt = time.time()
    EToForAccounts_2020 = pd.DataFrame()
    for i,row in acct_gdf.iterrows():
        print('\tAcct: '+str(row.customerid)+', Asset: '+str(row.fieldassetid)+', Ranch:'+row.SubCategory+', Field:'+row.Name)
        start_time = time.time()
        c = row.CentroidPoint
        
        etv = getWeather2ForDates(c,list_of_dates,Weather2020_gdf)
        
        #now add the columns for the other fields to etv
        etv['UUID']=row.UUID
        etv['Category']=row.Category
        etv['SubCategory']=row.SubCategory
        etv['Name']=row.Name

        EToForAccounts_2020 = EToForAccounts_2020.append(etv, ignore_index = True)
    print("ETo Done --- %s seconds ---" % (time.time() - startt))
    EToForAccounts = EToForAccounts_2020.rename(columns={"ReferenceETo(in)": "ETo(in)"})
    merged_df = pd.merge(acct_gdf, EToForAccounts, on=['Category','SubCategory','Name', 'UUID'])
    merged_df.reset_index() #reset the index to avoid any export issues
    merged_df['Category'] = 'JainlogicCustomerid_' + merged_df['customerid'].astype(str)
    merged_df['Date'] = pd.to_datetime(merged_df['Date'])

    return merged_df


def generateAndMergeETo2(acct_gdf,Weather2020_gdf,startDtStr, endDtStr):
    start = datetime.strptime(startDtStr,'%Y-%m-%d')
    end = datetime.strptime(endDtStr,'%Y-%m-%d')
    #-timedelta(days=1) #This date is NOT included

    list_of_dates = [start + timedelta(days=x+1) for x in range(0, (end-start).days)]

    startt = time.time()
    EToForAccounts_2020 = pd.DataFrame()
    for i,row in acct_gdf.iterrows():
        print('\tAcct: '+str(row.customerid)+', Asset: '+str(row.fieldassetid)+', Ranch:'+row.SubCategory+', Field:'+row.Name)
        start_time = time.time()
        c = row.CentroidPoint
        
        etv = getWeather3ForDates(c,list_of_dates,Weather2020_gdf)
        
        #now add the columns for the other fields to etv
        etv['UUID']=row.UUID
        etv['Category']=row.Category
        etv['SubCategory']=row.SubCategory
        etv['Name']=row.Name

        EToForAccounts_2020 = EToForAccounts_2020.append(etv, ignore_index = True)
    print("ETo Done --- %s seconds ---" % (time.time() - startt))
    EToForAccounts = EToForAccounts_2020.rename(columns={"ReferenceETo(in)": "ETo(in)"})
    merged_df = pd.merge(acct_gdf, EToForAccounts, on=['Category','SubCategory','Name', 'UUID'])
    merged_df.reset_index() #reset the index to avoid any export issues
    merged_df['Category'] = 'JainlogicCustomerid_' + merged_df['customerid'].astype(str)
    merged_df['Date'] = pd.to_datetime(merged_df['Date'])

    return merged_df

def generateETc(merged_df):

    unique_accounts = pd.unique(merged_df.Category)
    startt = time.time()
    #Now to put this into the output folder for each account
    print("Start --- %s  ---" % datetime.now())

    for a in unique_accounts:
        try:
            starta = time.time()
            
            print('\n\tWriting...%s ' % a)
            print("\tStart --- %s  ---" % datetime.now())
            acctName = a.upper()
            acct_dir = acctOutputPath + a + '/data/'
            acct_dir_img = acctOutputPath + a + '/img/'
            acct_dir_pkl = acctOutputPath + a + '/pkl/'
            acct_dir_rep = acctOutputPath + a + '/reports/'
                
            #create the directory if one doesn't exist
            if not os.path.exists(acct_dir):
                os.makedirs(acct_dir)
            if not os.path.exists(acct_dir_img):
                os.makedirs(acct_dir_img)
            if not os.path.exists(acct_dir_pkl):
                os.makedirs(acct_dir_pkl)
            if not os.path.exists(acct_dir_rep):
                os.makedirs(acct_dir_rep)
                
            #For each account, take the subset of the acct_df 
            acctDF = merged_df[merged_df.Category == a]
            
            acctGDF = gpd.GeoDataFrame(acctDF,geometry=acctDF.geometry)
            
            #for every day, do the following analysis and printout
            acct_dates_list = acctGDF.Date.tolist()
            unique_dts = pd.unique(acct_dates_list) # so we don't multiple count for each field
            unique_dts.sort()
            
            for dt in unique_dts:
                
                endDay = dt        
                
                startd = time.time()
                #print("\n\t\tSTART --- %s  ---" % datetime.now())

                etoDateStr =  datetime.strftime(dt,'%Y-%m-%d')
                jsonFile = acct_dir + acctName+'__ET_'+etoDateStr+'.json'
                
                acctGDF_date = acctGDF[acctGDF.Date==dt]
                acctGDF_date.reset_index(drop=True, inplace=True)
                
                print('\t\t'+etoDateStr)
                for i,row in acctGDF_date.iterrows():
                    
                    acct = row.Category
                    ranch = row.SubCategory
                    field = row.Name
                    fieldassetid = 'ASSETID_'+str(row.fieldassetid)
                    
                    #ETo and weather values
                    eto = row['ETo(in)']
                    tavg = row['Tavg']
                    rnet = row['Rnet']
                    wavg = row['Wavg']
                    ############## add others Raw and Deried
                    precip = row['Precipitation(in)']
                    avgSoilTemp = row['AverageSoilTemperature(F)']
                    maxAirTemp = row['MaximumAirTemperature(F)']
                    minAirTemp = row['MinimumAirTemperature(F)']
                    avgAirTemp = row['AverageAirTemperature(F)']
                    maxRelHumPercent = row['MaximumRelativeHumidity(%)']
                    minRelHumPercent = row['MinimumRelativeHumidity(%)']
                    avgRelHumPercent = row['AverageRelativeHumidity(%)']
                    dewPoint = row['DewPoint(F)']
                    windRun = row['WindRun(miles)']
                    
                    #########################################
                
                    # NDVI value
                    P = getSinglePoly(row.geometry)

                    
                    ss = row.SentinelFootprint
                    ll = row.LandsatFootprint
                    
                    #do Sentinel
                    if not type(ss)==list:
                        if (ss.find(',')==-1): #There is only ONE tile
                            S = [ss.replace('[','').replace(']','').replace('\'','')]
                        else:
                            tok = ss.split(',')
                            l = len(tok)
                            S=[]
                            for ii in range(l):
                                tk = tok[ii].replace('[','').replace(']','').replace('\'','')
                                S.append(tk)
                    else:
                        S = ss
                    if (type(S)==str):
                        S = [S]
                
                    #do landsat
                    if not type(ll)==list:
                        if (ss.find(',')==-1): #There is only ONE tile
                            L = ll.replace('[','').replace(']','').replace('\'','')
                        else:
                            tok = ll.split(', ')
                            l = len(tok)
                            L=[]
                            for ii in range(l):
                                tk = tok[ii].replace('[','').replace(']','').replace('\'','')
                                L.append(tk.replace(' ',''))
                    else:
                        L = [str(x).zfill(6) for x in ll]

                    if (type(L)==str):
                        L = [L]
                    
                    ndvi_mList_df = pd.DataFrame()
                    
                    NDVIFiles_df = getNDVIFilesWithDate(15, S, L, etoDateStr)
                    if len(NDVIFiles_df)==0:
                        NDVIFiles_df = getNDVIFilesWithDate(32, S, L, etoDateStr)
                    if len(NDVIFiles_df)==0:
                        NDVIFiles_df = getNDVIFilesWithDate(62, S, L, etoDateStr)

                    NDVIFiles_df.reset_index(drop=True,inplace=True)
                    #print(NDVIFiles_df)
                    
                    [ndvi_df,ndvi_img,ndvi_ext,imgDtStr] = getNDVIValueWithStatsAndImg3(NDVIFiles_df,P)
                    
                    
                    if (len(ndvi_df)>0 & (len(ndvi_img)>0)): #The first time the df has valid values
                        ext = json.loads(ndvi_ext)
                        ndvi_mList_df = ndvi_df.copy()
                        
                        
                        acct_dir_img2 = acct_dir_img + fieldassetid + '/'
                        if not os.path.exists(acct_dir_img2):
                            os.makedirs(acct_dir_img2)
                        imgFile = acct_dir_img2 + acctName+'_'+fieldassetid+'__IMG_'+etoDateStr+'.png'
                        infoJsonFile = acct_dir_img2 + acctName+'_'+fieldassetid+'__INFO_'+etoDateStr+'.json'

                        acct_dir_pkl2 = acct_dir_pkl + fieldassetid + '/'
                        if not os.path.exists(acct_dir_pkl2):
                            os.makedirs(acct_dir_pkl2)
                        pklFile = acct_dir_pkl2 + acctName+'_'+fieldassetid+'__PKL_'+etoDateStr+'.pkl'
                        
                    
                        imgDict = {
                            'ImgFileName': os.path.split(imgFile)[1],
                            'ImgDate': imgDtStr,
                            'ImgType': 'ETc',
                            'ImgExtent': ndvi_ext, 
                            'ImgBoundary': P.to_wkt(),
                            'AssetID': row.fieldassetid,
                            'CustomerID':a.split('_')[1],
                            'Ranch': ranch,
                            'Field': field
                        }
                        with open(infoJsonFile, 'w', encoding='utf8') as outfile:
                            json.dump(imgDict, outfile, ensure_ascii=False)

                        plt.ioff();
                        with open(imgFile, 'wb') as outfile:
                            fig = plt.figure()
                            plt.imshow(ndvi_img, extent=ext, cmap=zoom1, vmin=0, vmax=1)
                            fig.gca().set_axis_off()
                            plt.savefig(outfile, dpi=600, bbox_inches='tight', pad_inches=0, transparent=True)
                            plt.close();
                            
                        with open(pklFile,'wb') as f:
                            pickle.dump(ndvi_img, f)
                    
                        retdf = ndvi_df.loc[0]
                        ndvi_mean = retdf.Mean
                        ndvi_p90 = retdf.p90
                        ndvi_p25 = retdf.p25
                        ndvi_p75 = retdf.p75
                        ndvi_median = retdf.Median
                        ndvi_std = retdf.StdDev
                        
                    else:
                        #set values if NaN
                        ndvi_mean = np.nan
                        ndvi_p90 = np.nan
                        ndvi_p25 = np.nan
                        ndvi_p75 = np.nan
                        ndvi_median = np.nan
                        ndvi_std = np.nan
                        
                    #derived
                    ndvi_unif = (ndvi_p25 / ndvi_median)*100
                
                    #compute the METRIC ETc
                    et_metric = computeMetricET(tavg,rnet,wavg,ndvi_mean)
                    
                    ##########if et_metric is NaN, make it zero
                    if (np.isnan(et_metric)):
                        et_metric=0.0
                
                    #the Kc
                    fac = 1.35 #somewhat arbitrary for now
                    k_c = 1.3558*ndvi_mean*fac + 0.05
                    k_c2 = 1.3558*ndvi_mean + 0.05
                
                    #compute the ETc
                    et_ndvi = eto*k_c
                    et_ndvi2 = eto*k_c2
                
                    #get the final ETc
                    # Include et_metric
                    k_a = 1 #We can reduce k_a to increase the effect of Metric (see above)
                    et_c =  (1.15 - k_a)*et_metric + ((k_a + 0.0)*et_ndvi)
                    et_c2 =  (1.15 - k_a)*et_metric + ((k_a + 0.0)*et_ndvi2)
                    
                    k_c_p90 = ((1.3558*ndvi_p90*fac + 0.05) + (1.3558*ndvi_p90 + 0.05))/2
                    et_p90 = eto*k_c_p90
                    
                    acctGDF_date.loc[i,'ETc']=et_c
                    acctGDF_date.loc[i,'ETc_low']=et_c2
                    acctGDF_date.loc[i,'Kc']=k_c
                    acctGDF_date.loc[i,'ETm']=et_metric
                    acctGDF_date.loc[i,'ET_p90']=et_p90
                    acctGDF_date.loc[i,'ET_unif']=ndvi_unif
                    acctGDF_date.loc[i,'Precipitation(in)']=precip
                    acctGDF_date.loc[i,'AverageSoilTemperature(F)']=avgSoilTemp
                    acctGDF_date.loc[i,'MaximumAirTemperature(F)']=maxAirTemp
                    acctGDF_date.loc[i,'MinimumAirTemperature(F)']=minAirTemp
                    acctGDF_date.loc[i,'AverageAirTemperature(F)']=avgAirTemp
                    acctGDF_date.loc[i,'MaximumRelativeHumidity(%)']=maxRelHumPercent
                    acctGDF_date.loc[i,'MinimumRelativeHumidity(%)']=minRelHumPercent
                    acctGDF_date.loc[i,'AverageRelativeHumidity(%)']=avgRelHumPercent
                    acctGDF_date.loc[i,'DewPoint(F)']=dewPoint
                    acctGDF_date.loc[i,'WindRun(miles)']=windRun
                    

                acctGDF_date['SentinelFootprint'] = acctGDF_date['SentinelFootprint'].astype('str')
                acctGDF_date['LandsatFootprint'] = acctGDF_date['LandsatFootprint'].astype('str')
                acctGDF_date['Date'] = acctGDF_date['Date'].astype('str')
                acctGDF_date.drop('CentroidPoint',inplace=True, axis=1)
                print(acctGDF_date.ETc.tolist())
                
                #save the values
                acctGDF_date.to_file(jsonFile, driver="GeoJSON")

            print("\tEND --- %s seconds ---" % (time.time() - starta))
            print('\tDone!')
        except:
            print('Exception in Acct: '+str(a))
            continue

    print("\nEND --- %s seconds ---" % (time.time() - startt))
    return


def multiprocessing_etc(merged_df, dt):
    time.sleep(2)
    dt_new = dt.astype('M8[D]').astype('O')
    generate_MP_ETc(merged_df,dt_new)
    return

def generate_MP_ETc(merged_df, dt):

    unique_accounts = pd.unique(merged_df.Category)
    startt = time.time()
    #Now to put this into the output folder for each account
    print("Start --- %s  ---" % datetime.now())

    for a in unique_accounts:
        #try:
        starta = time.time()
        
        print('\n\tWriting...%s ' % a)
        print("\tStart --- %s  ---" % datetime.now())
        acctName = a.upper()
        acct_dir = acctOutputPath + a + '/data/'
        acct_dir_img = acctOutputPath + a + '/img/'
        acct_dir_pkl = acctOutputPath + a + '/pkl/'
        acct_dir_rep = acctOutputPath + a + '/reports/'
            
        #create the directory if one doesn't exist
        if not os.path.exists(acct_dir):
            os.makedirs(acct_dir)
        if not os.path.exists(acct_dir_img):
            os.makedirs(acct_dir_img)
        if not os.path.exists(acct_dir_pkl):
            os.makedirs(acct_dir_pkl)
        if not os.path.exists(acct_dir_rep):
            os.makedirs(acct_dir_rep)
            
        #For each account, take the subset of the acct_df 
        acctDF = merged_df[merged_df.Category == a]
        
        acctGDF = gpd.GeoDataFrame(acctDF,geometry=acctDF.geometry)
        
        #for every day, do the following analysis and printout
        acct_dates_list = acctGDF.Date.tolist()
        unique_dts = pd.unique(acct_dates_list) # so we don't multiple count for each field
        unique_dts.sort()
        
        ####################### (COMMENT AND UN-INDENT)
        #for dt in unique_dts:
        ###############################################
            
        endDay = dt        
        
        startd = time.time()
        #print("\n\t\tSTART --- %s  ---" % datetime.now())

        etoDateStr =  datetime.strftime(dt,'%Y-%m-%d')
        jsonFile = acct_dir + acctName+'__ET_'+etoDateStr+'.json'
        
        acctGDF_date = acctGDF[acctGDF.Date==dt]
        acctGDF_date.reset_index(drop=True, inplace=True)
        
        print('\t\t'+etoDateStr)
        for i,row in acctGDF_date.iterrows():
            
            acct = row.Category
            ranch = row.SubCategory
            field = row.Name
            fieldassetid = 'ASSETID_'+str(row.fieldassetid)
            
            #ETo and weather values
            eto = row['ETo(in)']
            tavg = row['Tavg']
            rnet = row['Rnet']
            wavg = row['Wavg']
            ############## add others Raw and Deried
            precip = row['Precipitation(in)']
            avgSoilTemp = row['AverageSoilTemperature(F)']
            maxAirTemp = row['MaximumAirTemperature(F)']
            minAirTemp = row['MinimumAirTemperature(F)']
            avgAirTemp = row['AverageAirTemperature(F)']
            maxRelHumPercent = row['MaximumRelativeHumidity(%)']
            minRelHumPercent = row['MinimumRelativeHumidity(%)']
            avgRelHumPercent = row['AverageRelativeHumidity(%)']
            dewPoint = row['DewPoint(F)']
            windRun = row['WindRun(miles)']
            
            #########################################
        
            # NDVI value
            P = getSinglePoly(row.geometry)

            
            ss = row.SentinelFootprint
            ll = row.LandsatFootprint
            
            #do Sentinel
            if not type(ss)==list:
                if (ss.find(',')==-1): #There is only ONE tile
                    S = [ss.replace('[','').replace(']','').replace('\'','')]
                else:
                    tok = ss.split(',')
                    l = len(tok)
                    S=[]
                    for ii in range(l):
                        tk = tok[ii].replace('[','').replace(']','').replace('\'','')
                        S.append(tk)
            else:
                S = ss
            if (type(S)==str):
                S = [S]
        
            #do landsat
            if not type(ll)==list:
                if (ss.find(',')==-1): #There is only ONE tile
                    L = ll.replace('[','').replace(']','').replace('\'','')
                else:
                    tok = ll.split(', ')
                    l = len(tok)
                    L=[]
                    for ii in range(l):
                        tk = tok[ii].replace('[','').replace(']','').replace('\'','')
                        L.append(tk.replace(' ',''))
            else:
                L = [str(x).zfill(6) for x in ll]

            if (type(L)==str):
                L = [L]
            
            ndvi_mList_df = pd.DataFrame()
            
            NDVIFiles_df = getNDVIFilesWithDate(32, S, L, etoDateStr)
            NDVIFiles_df.reset_index(drop=True,inplace=True)
            #print(NDVIFiles_df)
            
            [ndvi_df,ndvi_img,ndvi_ext,imgDtStr] = getNDVIValueWithStatsAndImg3(NDVIFiles_df,P)
            
            
            if (len(ndvi_df)>0 & (len(ndvi_img)>0)): #The first time the df has valid values
                ext = json.loads(ndvi_ext)
                ndvi_mList_df = ndvi_df.copy()
                
                
                acct_dir_img2 = acct_dir_img + fieldassetid + '/'
                if not os.path.exists(acct_dir_img2):
                    os.makedirs(acct_dir_img2)
                imgFile = acct_dir_img2 + acctName+'_'+fieldassetid+'__IMG_'+etoDateStr+'.png'
                infoJsonFile = acct_dir_img2 + acctName+'_'+fieldassetid+'__INFO_'+etoDateStr+'.json'

                acct_dir_pkl2 = acct_dir_pkl + fieldassetid + '/'
                if not os.path.exists(acct_dir_pkl2):
                    os.makedirs(acct_dir_pkl2)
                pklFile = acct_dir_pkl2 + acctName+'_'+fieldassetid+'__PKL_'+etoDateStr+'.pkl'
                
            
                imgDict = {
                    'ImgFileName': os.path.split(imgFile)[1],
                    'ImgDate': imgDtStr,
                    'ImgType': 'ETc',
                    'ImgExtent': ndvi_ext, 
                    'ImgBoundary': P.to_wkt(),
                    'AssetID': row.fieldassetid,
                    'CustomerID':a.split('_')[1],
                    'Ranch': ranch,
                    'Field': field
                }
                with open(infoJsonFile, 'w', encoding='utf8') as outfile:
                    json.dump(imgDict, outfile, ensure_ascii=False)

                plt.ioff();
                with open(imgFile, 'wb') as outfile:
                    fig = plt.figure()
                    plt.imshow(ndvi_img, extent=ext, cmap=zoom1, vmin=0, vmax=1)
                    fig.gca().set_axis_off()
                    plt.savefig(outfile, dpi=600, bbox_inches='tight', pad_inches=0, transparent=True)
                    plt.close();
                    
                with open(pklFile,'wb') as f:
                    pickle.dump(ndvi_img, f)
            
                retdf = ndvi_df.loc[0]
                ndvi_mean = retdf.Mean
                ndvi_p90 = retdf.p90
                ndvi_p25 = retdf.p25
                ndvi_p75 = retdf.p75
                ndvi_median = retdf.Median
                ndvi_std = retdf.StdDev
                
            else:
                #set values if NaN
                ndvi_mean = np.nan
                ndvi_p90 = np.nan
                ndvi_p25 = np.nan
                ndvi_p75 = np.nan
                ndvi_median = np.nan
                ndvi_std = np.nan
                
                
            #derived
            ndvi_unif = (ndvi_p25 / ndvi_median)*100
        
            #compute the METRIC ETc
            et_metric = computeMetricET(tavg,rnet,wavg,ndvi_mean)
            
            ##########if et_metric is NaN, make it zero
            if (np.isnan(et_metric)):
                et_metric=0.0
        
            #the Kc
            fac = 1.35 #somewhat arbitrary for now
            k_c = 1.3558*ndvi_mean*fac + 0.05
            k_c2 = 1.3558*ndvi_mean + 0.05
        
            #compute the ETc
            et_ndvi = eto*k_c
            et_ndvi2 = eto*k_c2
        
            #get the final ETc
            # Include et_metric
            k_a = 1 #We can reduce k_a to increase the effect of Metric (see above)
            et_c =  (1.15 - k_a)*et_metric + ((k_a + 0.0)*et_ndvi)
            et_c2 =  (1.15 - k_a)*et_metric + ((k_a + 0.0)*et_ndvi2)
            
            k_c_p90 = ((1.3558*ndvi_p90*fac + 0.05) + (1.3558*ndvi_p90 + 0.05))/2
            et_p90 = eto*k_c_p90
            
            acctGDF_date.loc[i,'ETc']=et_c
            acctGDF_date.loc[i,'ETc_low']=et_c2
            acctGDF_date.loc[i,'Kc']=k_c
            acctGDF_date.loc[i,'ETm']=et_metric
            acctGDF_date.loc[i,'ET_p90']=et_p90
            acctGDF_date.loc[i,'ET_unif']=ndvi_unif
            acctGDF_date.loc[i,'Precipitation(in)']=precip
            acctGDF_date.loc[i,'AverageSoilTemperature(F)']=avgSoilTemp
            acctGDF_date.loc[i,'MaximumAirTemperature(F)']=maxAirTemp
            acctGDF_date.loc[i,'MinimumAirTemperature(F)']=minAirTemp
            acctGDF_date.loc[i,'AverageAirTemperature(F)']=avgAirTemp
            acctGDF_date.loc[i,'MaximumRelativeHumidity(%)']=maxRelHumPercent
            acctGDF_date.loc[i,'MinimumRelativeHumidity(%)']=minRelHumPercent
            acctGDF_date.loc[i,'AverageRelativeHumidity(%)']=avgRelHumPercent
            acctGDF_date.loc[i,'DewPoint(F)']=dewPoint
            acctGDF_date.loc[i,'WindRun(miles)']=windRun
                

            acctGDF_date['SentinelFootprint'] = acctGDF_date['SentinelFootprint'].astype('str')
            acctGDF_date['LandsatFootprint'] = acctGDF_date['LandsatFootprint'].astype('str')
            acctGDF_date['Date'] = acctGDF_date['Date'].astype('str')
            if 'CentroidPoint' in acctGDF_date.columns:
                acctGDF_date.drop('CentroidPoint',inplace=True, axis=1)
            print(acctGDF_date.ETc.tolist())
            
            #save the values
            acctGDF_date.to_file(jsonFile, driver="GeoJSON")

        #print("\tEND --- %s seconds ---" % (time.time() - starta))
        #print('\tDone!')
        #except:
        #	print('Exception in Acct: '+str(a))
        #	continue

    #print("\nEND --- %s seconds ---" % (time.time() - startt))
    return


def getFeatures(gdf):
    """Function to parse features from GeoDataFrame in such a manner that rasterio wants them"""
    return [json.loads(gdf.to_json())['features']]


# In[3]:


def getCCGeoDataFrame(gmlPath):

    gmlFileStr = 'MSK_CLOUDS_B00.gml'
    gmlFile = gmlPath + gmlFileStr
    
    f = open(gmlFile, 'r')
    lines = f.readlines()
    f.close()
    
    start=0
    name=0
    points = []
    data = []
    for l in lines:
        line = l.strip()
    
        if (line[0:13] =='<gml:Envelope'):
            l = line.strip().split(':')
            proj = l[-1][0:-2]
    
        if (line[0:12]=='<gml:Polygon'):
            l = line.strip()[13:-1]
            tok = l.split('=')
            name = tok[1].replace('\"','')

        if (line[0:12]=='<gml:posList'):
            start=1
        if (line=='</gml:LinearRing>'):
            start=0
        
        if start==1:
            #print(name)
            type = name.split('.')[0]
            s = line[30:-14]
            c = s.split(' ')
        
            r = int(len(c)/2)
            points=[]
            for i in range(r):
                points.append(Point([int(c[i*2]),int(c[i*2+1])]))
            #male polygon and Add to the geopandas    
            coords = sum(map(list, (p.coords for p in points)), [])
            poly = Polygon(coords)
            data.append(dict({'Cloud_Name': name, 'Mask_Type': type, 'geometry': poly}))
    
    df = pd.DataFrame(data)
    if (len(df)>0):
        e = 'epsg:'+proj
        crs = {'init': e}
        cloudsGDF = gpd.GeoDataFrame(df, crs=crs, geometry=df.geometry)
        return cloudsGDF
    else:
        nullDataFrame = gpd.GeoDataFrame()
        return nullDataFrame
    


# In[4]:


def getMostRecentFile(sentinel_tiff_files_array):
    dts = []
    #print(sentinel_tiff_files_array)
    for file in sentinel_tiff_files_array:
        d_str = file.split('_')[2][0:8]
        dt = parse(d_str)
        dts.append(dt)
    idx = [dts.index(x) for x in sorted(dts)]
    
    return sentinel_tiff_files_array[idx[-1]]


# In[5]:


def getAreaOfROI(roi_polygon):
    #Returns area in KM^2 and the fraction of the 100x100km Sentinel

    geom = shapely.geometry.Polygon(roi_polygon)
    #geom = Polygon([(0, 0), (0, 10), (10, 10), (10, 0), (0, 0)])
    geom_area = shapely.ops.transform(
        partial(
            pyproj.transform,
            pyproj.Proj(init='EPSG:4326'),
            pyproj.Proj(
                proj='aea',
                lat_1=geom.bounds[1],
                lat_2=geom.bounds[3])),
        geom)

    # return the area in Km^2
    a = (geom_area.area / 1000000)
    f = a / 12056.04 #109.8x 109.8km = 12,056.04 km ^2
    return [a,f]


# In[6]:


def getFootprint(XML_file):
    f = open(XML_file, "r")
    lines = f.readlines()
    s=0
    n=0
    crd = []

    for line in lines:
        n=n+1
        if (line.strip().startswith('<Global_Footprint>')):
            s=1
            n=0
        if (line.strip().startswith('</Global_Footprint>')):
            s=0
        if s==1 & n==1:
            l = line.strip()
            t = l.split('<EXT_POS_LIST>')
            c = t[1][:-15].strip()
            crds = c.split(' ')
        
            for i in range(int(len(crds)/2)):
                lat = float(crds[i*2])
                lon = float(crds[i*2+1])
                ll = [lon,lat]
                crd.append(ll)
    return crd

def getBounds(tileStr):
    tile_df = S_gdf[S_gdf.Name==tileStr].geometry.bounds
    tile_df.reset_index(drop=True,inplace=True)
    bounds = [
        tile_df.loc[0,'minx'],
        tile_df.loc[0,'miny'],
        tile_df.loc[0,'maxx'],
        tile_df.loc[0,'maxy']
    ]
    
    return bounds

def getWMSData(from_dt,to_dt,bounding_box):
    INSTANCE_ID = '5c53cc29-8f4a-43ed-b97b-cb27be196176'
    LAYER_NAME = '3_NDVI'
    bands_script = 'return [B01,B02,B04,B05,B08,B8A,B09,B10,B11,B12]'
    wms_bands_request = WmsRequest(layer=LAYER_NAME,
                               custom_url_params={
                                   CustomUrlParam.EVALSCRIPT: bands_script,
                                   CustomUrlParam.ATMFILTER: 'NONE'
                               },
                               bbox=bounding_box, 
                               time=(from_dt, to_dt), 
                               width=600, height=None,
                               image_format=MimeType.TIFF_d32f,
                               instance_id=INSTANCE_ID)
    wms_bands = wms_bands_request.get_data()
    return wms_bands

def plot_cloud_mask(mask, figsize=(15, 15), fig=None):
    """
    Utility function for plotting a binary cloud mask.
    """
    if fig == None:
        plt.figure(figsize=figsize)
    plt.imshow(mask, cmap=plt.cm.gray)
    

def getCloudMask(imagePath,imgFile):
    upscale_factor = 2
    with rasterio.open(imagePath+imgFile, driver='JP2OpenJPEG') as dataset:
        # resample data to target shape
        img = dataset.read(
            out_shape=(
                dataset.count,
                int(dataset.height * upscale_factor),
                int(dataset.width * upscale_factor)
            ),
            resampling=Resampling.bilinear
        )

        # scale image transform
        transform = dataset.transform * dataset.transform.scale(
            (dataset.width / img.shape[-1]),
            (dataset.height / img.shape[-2])
        )
        img = img.squeeze()

    img = img.astype('float')   
    img[img>0]=np.nan
    img[img==0]=1.0    # Set all probabilitiy > 0 as 1, "MOST CONSERVATIVE"
    
    mask = img.copy()
    
    return mask

def getTileWkt(tile):
    df = S_gdf[S_gdf.TILE==tile].geometry
    df.reset_index(drop=True,inplace=True)
    geom = df.geometry[0][0]
    wkt = geom.to_wkt()
    coords = list(geom.exterior.coords)
    
    return wkt


def interpolateNaNs(df, colStr):
    
    new_gdf = gpd.GeoDataFrame()
    df[colStr] = df[colStr].astype('float')
    
    list_of_unique_fieldassetids = pd.unique(df.fieldassetid).tolist()
    list_of_unique_fieldassetids = sorted(list_of_unique_fieldassetids, key = int)
    
    for assetid in list_of_unique_fieldassetids:
        dff = df[df.fieldassetid==assetid]
        dff.reset_index(drop=True, inplace=True)
        dff[colStr] = dff[colStr].interpolate(method='linear', limit_direction='both', axis=0)
        #dff[colStr] = dff[colStr].bfill() #back fill any trailing NaNs
        #dff[colStr] = dff[colStr].interpolate(method='pad') #forward fill any trailing NaNs
    

        if len(dff[dff[colStr].isnull()])==0:
            new_gdf = new_gdf.append(dff, ignore_index=True)
        else:
            continue #Ignore this assetID
        
    new_gdf.reset_index(drop=True, inplace=True)
    
    return new_gdf


def interpolateNaNs2(df, colStr):
    
    new_gdf = gpd.GeoDataFrame()
    df[colStr] = df[colStr].astype('float')
    
    list_of_unique_fieldassetids = pd.unique(df.Name).tolist()
    list_of_unique_fieldassetids.sort()
    
    for assetid in list_of_unique_fieldassetids:
        dff = df[df.Name==assetid]
        dff.reset_index(drop=True, inplace=True)
        dff[colStr] = dff[colStr].interpolate(method='linear', limit_direction='both', axis=0)

        if len(dff[dff[colStr].isnull()])==0:
            new_gdf = new_gdf.append(dff, ignore_index=True)
        else:
            continue 
        
    new_gdf.reset_index(drop=True, inplace=True)
    
    return new_gdf

def getTS_StartEndCols(acctData_gdf):
    #Pivot to get the daily ET for each assetID
    df = []
    df = pd.pivot_table(data=acctData_gdf,
                               values='ETc_average',
                               columns=['fieldassetid'],
                               index='Date', aggfunc=np.sum)
    df.reset_index(inplace=True)
    df['Date']=pd.to_datetime(df.Date)
    df = df.set_index('Date',drop=True)
    #interpolate any missing values
    df = df.interpolate(method='linear', limit_direction='both', axis=0)

    startDtStr = datetime.strftime(min(df.index),"%Y-%m-%dT%H:%M:%SZ")
    endDtStr = datetime.strftime(max(df.index),"%Y-%m-%dT%H:%M:%SZ")
    
    res_et = pd.DataFrame(columns = df.columns)
    #Write the sum to the result df
    res_et.loc[0,:]=df.sum(axis=0)
    et_dts = res_et.T
    et_dts.columns=['et']
    et_dts['start_ts'] = startDtStr
    et_dts['end_ts'] = endDtStr
    et_dts.drop('et',axis=1,inplace=True)
    
    et_dts.index = et_dts.index.map(int)
    
    return et_dts

def getTS_StartEndCols2(acctData_gdf):
    #Pivot to get the daily ET for each assetID
    df = []
    df = pd.pivot_table(data=acctData_gdf,
                               values='ETc_average',
                               columns=['Name'],
                               index='Date', aggfunc=np.sum)
    df.reset_index(inplace=True)
    df['Date']=pd.to_datetime(df.Date)
    df = df.set_index('Date',drop=True)
    #interpolate any missing values
    df = df.interpolate(method='linear', limit_direction='both', axis=0)

    startDtStr = datetime.strftime(min(df.index),"%Y-%m-%dT%H:%M:%SZ")
    endDtStr = datetime.strftime(max(df.index),"%Y-%m-%dT%H:%M:%SZ")
    
    res_et = pd.DataFrame(columns = df.columns)
    #Write the sum to the result df
    res_et.loc[0,:]=df.sum(axis=0)
    et_dts = res_et.T
    et_dts.columns=['et']
    et_dts['start_ts'] = startDtStr
    et_dts['end_ts'] = endDtStr
    et_dts.drop('et',axis=1,inplace=True)
    
    #et_dts.index = et_dts.index.map(int)
    
    return et_dts

def getWkt(acct_gdf):
    et_geoms = pd.DataFrame(columns = ['fieldassetid','assetWkt'])
    for i,row in acct_gdf.iterrows():
        et_geoms.loc[i,'fieldassetid'] = row.fieldassetid
        et_geoms.loc[i,'assetWkt'] = row.geometry.to_wkt()
    
    et_geoms.set_index('fieldassetid', inplace=True, drop=True)
    
    et_geoms.index = et_geoms.index.map(int)
    
    return et_geoms

def getWkt2(acct_gdf):
    et_geoms = pd.DataFrame(columns = ['Name','assetWkt'])
    for i,row in acct_gdf.iterrows():
        et_geoms.loc[i,'Name'] = row.Name
        et_geoms.loc[i,'assetWkt'] = row.geometry.to_wkt()
    
    et_geoms.set_index('Name', inplace=True, drop=True)
    #et_geoms.index = et_geoms.index.map(int)
    
    return et_geoms

def getETCol(acctData_gdf):
    #Pivot to get the daily ET for each assetID
    df = []

    acctData_gdf.ETc_average = acctData_gdf.ETc_average.astype('float')
    df = pd.pivot_table(data=acctData_gdf,
                               values='ETc_average',
                               columns=['fieldassetid'],
                               index='Date', aggfunc=np.sum)
    df.reset_index(inplace=True)
    df['Date']=pd.to_datetime(df.Date)
    df = df.set_index('Date',drop=True)
    #interpolate any missing values
    df = df.interpolate(method='linear', limit_direction='both', axis=0)

    #sum to get the sum for the entire week
    res_et = pd.DataFrame(columns = df.columns)
    #Write the sum to the result df
    res_et.loc[0,:]=df.sum(axis=0)
    et_avg = res_et.T
    et_avg.columns=['et']
    
    et_avg.index = et_avg.index.map(int)

    
    return et_avg

def getETCol2(acctData_gdf):
    #Pivot to get the daily ET for each assetID
    df = []

    acctData_gdf.ETc_average = acctData_gdf.ETc_average.astype('float')
    df = pd.pivot_table(data=acctData_gdf,
                               values='ETc_average',
                               columns=['Name'],
                               index='Date', aggfunc=np.sum)
    df.reset_index(inplace=True)
    df['Date']=pd.to_datetime(df.Date)
    df = df.set_index('Date',drop=True)
    #interpolate any missing values
    df = df.interpolate(method='linear', limit_direction='both', axis=0)

    #sum to get the sum for the entire week
    res_et = pd.DataFrame(columns = df.columns)
    #Write the sum to the result df
    res_et.loc[0,:]=df.sum(axis=0)
    et_avg = res_et.T
    et_avg.columns=['et']
    
    #et_avg.index = et_avg.index.map(int)

    
    return et_avg

def getET90Col(acctData_gdf):
    #Pivot to get the daily ET for each assetID
    df = []

    acctData_gdf.ET_p90 = acctData_gdf.ET_p90.astype('float')
    df = pd.pivot_table(data=acctData_gdf,
                               values='ET_p90',
                               columns=['fieldassetid'],
                               index='Date', aggfunc=np.sum)
    df.reset_index(inplace=True)
    df['Date']=pd.to_datetime(df.Date)
    df = df.set_index('Date',drop=True)
    #interpolate any missing values
    df = df.interpolate(method='linear', limit_direction='both', axis=0)

    #sum to get the sum for the entire week
    res_et = pd.DataFrame(columns = df.columns)
    #Write the sum to the result df
    res_et.loc[0,:]=df.sum(axis=0)
    et_90 = res_et.T
    et_90.columns=['et90']
    
    et_90.index = et_90.index.map(int)

    return et_90

def getET90Col2(acctData_gdf):
    #Pivot to get the daily ET for each assetID
    df = []

    acctData_gdf.ET_p90 = acctData_gdf.ET_p90.astype('float')
    df = pd.pivot_table(data=acctData_gdf,
                               values='ET_p90',
                               columns=['Name'],
                               index='Date', aggfunc=np.sum)
    df.reset_index(inplace=True)
    df['Date']=pd.to_datetime(df.Date)
    df = df.set_index('Date',drop=True)
    #interpolate any missing values
    df = df.interpolate(method='linear', limit_direction='both', axis=0)

    #sum to get the sum for the entire week
    res_et = pd.DataFrame(columns = df.columns)
    #Write the sum to the result df
    res_et.loc[0,:]=df.sum(axis=0)
    et_90 = res_et.T
    et_90.columns=['et90']
    
    #et_90.index = et_90.index.map(int)

    return et_90

def getETUnifCol(acctData_gdf):
    #Pivot to get the daily ET for each assetID
    df = []
    
    acctData_gdf.ET_unif = acctData_gdf.ET_unif.astype('float')
    df = pd.pivot_table(data=acctData_gdf,
                               values='ET_unif',
                               columns=['fieldassetid'],
                               index='Date', aggfunc=np.mean)
    df.reset_index(inplace=True)
    df['Date']=pd.to_datetime(df.Date)
    df = df.set_index('Date',drop=True)
    #interpolate any missing values
    df = df.interpolate(method='linear', limit_direction='both', axis=0)

    #sum to get the sum for the entire week
    res_et = pd.DataFrame(columns = df.columns)
    #Write the sum to the result df
    res_et.loc[0,:]=df.mean(axis=0)
    et_unif = res_et.T
    et_unif.columns=['etunif']
    
    et_unif.index = et_unif.index.map(int)

    return et_unif

def getETUnifCol2(acctData_gdf):
    #Pivot to get the daily ET for each assetID
    df = []
    
    acctData_gdf.ET_unif = acctData_gdf.ET_unif.astype('float')
    df = pd.pivot_table(data=acctData_gdf,
                               values='ET_unif',
                               columns=['Name'],
                               index='Date', aggfunc=np.mean)
    df.reset_index(inplace=True)
    df['Date']=pd.to_datetime(df.Date)
    df = df.set_index('Date',drop=True)
    #interpolate any missing values
    df = df.interpolate(method='linear', limit_direction='both', axis=0)

    #sum to get the sum for the entire week
    res_et = pd.DataFrame(columns = df.columns)
    #Write the sum to the result df
    res_et.loc[0,:]=df.mean(axis=0)
    et_unif = res_et.T
    et_unif.columns=['etunif']
    
    #et_unif.index = et_unif.index.map(int)

    return et_unif

def getPrecipCol(acctData_gdf):
    #Pivot to get the daily ET for each assetID
    df = []
    
    acctData_gdf.ET_unif = acctData_gdf.ET_unif.astype('float')
    df = pd.pivot_table(data=acctData_gdf,
                               values='Precipitation(in)',
                               columns=['fieldassetid'],
                               index='Date', aggfunc=np.mean)
    df.reset_index(inplace=True)
    df['Date']=pd.to_datetime(df.Date)
    df = df.set_index('Date',drop=True)
    #interpolate any missing values
    df = df.interpolate(method='linear', limit_direction='both', axis=0)

    #sum to get the sum for the entire week
    res_et = pd.DataFrame(columns = df.columns)
    #Write the sum to the result df
    res_et.loc[0,:]=df.mean(axis=0)
    precip = res_et.T
    precip.columns=['precip']
    
    precip.index = precip.index.map(int)

    return precip

def getPrecipCol2(acctData_gdf):
    #Pivot to get the daily ET for each assetID
    df = []
    
    acctData_gdf.ET_unif = acctData_gdf.ET_unif.astype('float')
    df = pd.pivot_table(data=acctData_gdf,
                               values='Precipitation(in)',
                               columns=['Name'],
                               index='Date', aggfunc=np.mean)
    df.reset_index(inplace=True)
    df['Date']=pd.to_datetime(df.Date)
    df = df.set_index('Date',drop=True)
    #interpolate any missing values
    df = df.interpolate(method='linear', limit_direction='both', axis=0)

    #sum to get the sum for the entire week
    res_et = pd.DataFrame(columns = df.columns)
    #Write the sum to the result df
    res_et.loc[0,:]=df.mean(axis=0)
    precip = res_et.T
    precip.columns=['precip']
    
    #precip.index = precip.index.map(int)

    return precip

def getETImgInfo(imgData_df):
    j=0
    etImgInfo = pd.DataFrame(columns = ['fieldassetid','ImageInfo'])
    
    list_of_unique_assetids  = pd.unique(imgData_df.fieldassetid)
    for asset in list_of_unique_assetids:
        df = imgData_df[imgData_df.fieldassetid==asset]
        df.reset_index(drop=True, inplace=True)
        jsonStr = df.to_json(orient='table')
        etImgInfo.loc[j,'fieldassetid'] = asset
        etImgInfo.loc[j,'ImageInfo'] = jsonStr
        j=j+1
    etImgInfo.set_index('fieldassetid',drop=True,inplace=True)
    
    return etImgInfo




def getETImgInfo2(imgData_df):
    j=0
    etImgInfo = pd.DataFrame(columns = ['fieldassetid','ImageInfo'])
    
    list_of_unique_assetids  = pd.unique(imgData_df.fieldassetid)
    for asset in list_of_unique_assetids:
        df = imgData_df[imgData_df.fieldassetid==asset]
        df.reset_index(drop=True, inplace=True)
        df = df.T
        jsonStr = df.to_json()    
        etImgInfo.loc[j,'fieldassetid'] = asset
        etImgInfo.loc[j,'ImageInfo'] = jsonStr
        j=j+1
    etImgInfo.fieldassetid = etImgInfo.fieldassetid.astype('int')
    etImgInfo.set_index('fieldassetid',drop=True,inplace=True)
    
    
    return etImgInfo


def getETImgInfo3(imgData_df):
    j=0
    etImgInfo = pd.DataFrame(columns = ['fieldassetid','ImageInfo'])
    

    list_of_unique_assetids  = pd.unique(imgData_df.fieldassetid).tolist()
    list_of_unique_assetids.sort()

    jsonDictList = []
    assetList = []
    for a in list_of_unique_assetids:
        df = imgData_df[imgData_df.fieldassetid==a]
        df.reset_index(drop=True, inplace=True)
        df = df.T
        jsonDict = df.to_dict()
        jsonDictList.append(jsonDict[0])
        assetList.append(int(a))
    etImgInfoSeries = pd.Series(jsonDictList, index=assetList)

    
    return etImgInfoSeries

def getETImgInfo4(imgData_df):

    j=0
    etImgInfo = pd.DataFrame(columns = ['Name','ImageInfo'])
    
    list_of_unique_assetids  = pd.unique(imgData_df.Name).tolist()
    list_of_unique_assetids.sort()

    jsonDictList = []
    assetList = []
    for a in list_of_unique_assetids:
        df = imgData_df[imgData_df.Name==a.upper()]
        df.reset_index(drop=True, inplace=True)
        df = df.T
        jsonDict = df.to_dict()
        jsonDictList.append(jsonDict[0])
        assetList.append(int(a))
    etImgInfoSeries = pd.Series(jsonDictList, index=assetList)

    
    return etImgInfoSeries

def getETImgInfo5(imgData_df):
    #print(imgData_df.columns)
    j=0
    etImgInfo = pd.DataFrame(columns = ['Field','ImageInfo'])
    
    list_of_unique_assetids  = pd.unique(imgData_df.Field).tolist()
    list_of_unique_assetids.sort()

    #print('getETImgInfo5: '+'\n')
    #print(list_of_unique_assetids)

    jsonDictList = []
    assetList = []
    for a in list_of_unique_assetids:
        df = imgData_df[imgData_df.Field.str.upper()==str(a).upper()]
        df.reset_index(drop=True, inplace=True)
        df = df.T
        jsonDict = df.to_dict()
        jsonDictList.append(jsonDict)
        assetList.append(a)
    etImgInfoSeries = pd.Series(jsonDictList, index=assetList)

    
    return etImgInfoSeries

def uploadHyperViewToAzure(SOURCE_FILE, acct, endDateStr):
    #base_url = 'https://api.jainlogic.com/irrigation/assets/irrigationsets/'
    #sas_url = base_url + 'agralogics-archives?sv=2019-02-02&si=agralogics-all&sr=c&sig=fqJ8dhIZ7wpn0xLL4tM5RzpumzJL5udWJY51iSjrNH8%3D'
    sas_url = 'https://jlsatelliteimagery.blob.core.windows.net/agralogics-archives?sv=2019-02-02&si=agralogics-all&sr=c&sig=fqJ8dhIZ7wpn0xLL4tM5RzpumzJL5udWJY51iSjrNH8%3D'
    container_client = ContainerClient.from_container_url(sas_url)

    uploadTimeStamp = datetime.strftime(datetime.today(),'%Y-%m-%dT%H:%M:%SZ')
    try:
        with open(SOURCE_FILE, "rb") as data:
            blob_client = container_client.upload_blob(name='Agralogics_HyperView__CustomerID_'+acct+'_onDate_'+endDateStr+ '_uploadTimePST_'+uploadTimeStamp, data=data)
        properties = blob_client.get_blob_properties()
        
    except Exception as ex:
        print('Exception in uploadToAzure:')
        print(ex)
        
    return properties


def uploadHyperGrowToAzure(SOURCE_FILE, acct, startDateStr, endDateStr):

    #base_url = 'https://api.jainlogic.com/irrigation/assets/irrigationsets/'
    #sas_url = base_url + 'agralogics-archives?sv=2019-02-02&si=agralogics-all&sr=c&sig=fqJ8dhIZ7wpn0xLL4tM5RzpumzJL5udWJY51iSjrNH8%3D'
    sas_url = 'https://jlsatelliteimagery.blob.core.windows.net/agralogics-archives?sv=2019-02-02&si=agralogics-all&sr=c&sig=fqJ8dhIZ7wpn0xLL4tM5RzpumzJL5udWJY51iSjrNH8%3D'
    container_client = ContainerClient.from_container_url(sas_url)

    uploadTimeStamp = datetime.strftime(datetime.today(),'%Y-%m-%dT%H:%M:%SZ')
    try:
        with open(SOURCE_FILE, "rb") as data:
            blob_client = container_client.upload_blob(name='Agralogics_HyperGrow__CustomerID_'+acct+'_startDate_'+startDateStr+'_endDate_'+endDateStr+ '_uploadTimePST_'+uploadTimeStamp, data=data)
        properties = blob_client.get_blob_properties()
        
    except Exception as ex:
        print('Exception in uploadToAzure:')
        print(ex)
        
    return properties

def uploadHyperViewToAzure_TEST(SOURCE_FILE, acct, endDateStr):

    #base_url = 'https://api.jainlogic.com/irrigation/assets/irrigationsets/'
    #sas_url = base_url + 'agralogics-archives?sv=2019-02-02&si=agralogics-all&sr=c&sig=4ph4aZuOAqWZGDuFIPU6tOKqFauAf9wiNf7vT9SmcCU%3D'
    sas_url = 'https://jlsatelliteimagerytest.blob.core.windows.net/agralogics-archives?sv=2019-02-02&si=agralogics-all&sr=c&sig=4ph4aZuOAqWZGDuFIPU6tOKqFauAf9wiNf7vT9SmcCU%3D'
    
    container_client = ContainerClient.from_container_url(sas_url)

    uploadTimeStamp = datetime.strftime(datetime.today(),'%Y-%m-%dT%H:%M:%SZ')
    try:
        with open(SOURCE_FILE, "rb") as data:
            blob_client = container_client.upload_blob(name='Agralogics_HyperView__CustomerID_'+acct+'_onDate_'+endDateStr+ '_uploadTimePST_'+uploadTimeStamp, data=data)
        properties = blob_client.get_blob_properties()
        
    except Exception as ex:
        print('Exception in uploadHyperViewToAzure (TEST):')
        print(ex)
        
    return properties

def uploadHyperGrowToAzure_TEST(SOURCE_FILE, acct, startDateStr, endDateStr):
    #print([SOURCE_FILE, acct, startDateStr, endDateStr])
    #base_url = 'https://api.jainlogic.com/irrigation/assets/irrigationsets/'
    #sas_url = base_url + 'agralogics-archives?sv=2019-02-02&si=agralogics-all&sr=c&sig=4ph4aZuOAqWZGDuFIPU6tOKqFauAf9wiNf7vT9SmcCU%3D'
    sas_url = 'https://jlsatelliteimagerytest.blob.core.windows.net/agralogics-archives?sv=2019-02-02&si=agralogics-all&sr=c&sig=4ph4aZuOAqWZGDuFIPU6tOKqFauAf9wiNf7vT9SmcCU%3D'
    
    properties = []
    
    container_client = ContainerClient.from_container_url(sas_url)

    uploadTimeStamp = datetime.strftime(datetime.today(),'%Y-%m-%dT%H:%M:%SZ')
    #try:
    with open(SOURCE_FILE, "rb") as data:
        blob_client = container_client.upload_blob(name='Agralogics_HyperGrow__CustomerID_'+acct+'_startDate_'+startDateStr+'_endDate_'+endDateStr+ '_uploadTimePST_'+uploadTimeStamp, data=data)
    properties = blob_client.get_blob_properties()
        
    #except Exception as ex:
    #    print('Exception in uploadToAzure (TEST):')
    #    print(ex)
        
    return properties


def dropNaNValueAssets(results_df):
    new_df = pd.DataFrame()
    df = results_df.copy()
    numeric_cols = [
    'Vigor','Uniformity',
    '1-Period-PDiff_PAcres',
    '1-Period-PDiff_NAcres',
    '1-Period-PDiff_PPercent',
    '1-Period-PDiff_NPercent',
    '4-Period-PDiff_PAcres',
    '4-Period-PDiff_NAcres',
    '4-Period-PDiff_PPercent',
    '4-Period-PDiff_NPercent',
    ]
    df[numeric_cols] = df[numeric_cols].astype('float')
    is_NaN = df[numeric_cols].isnull()
    row_has_NaN = is_NaN.any(axis=1)
    nan_df = df[row_has_NaN]
    
    
    new_df.reset_index(drop=True, inplace=True)
    
    if len(nan_df)==0:
        new_df = df.copy()
    else: #Drop these rows and return the rest 
        print('Dropped Rows: \n')
        print(nan_df)
        new_df = df[~row_has_NaN]
    
    new_df.reset_index(drop=True, inplace=True)
    
    return new_df


def saveHyperViewFig(img,cmp,vmin,vmax,imgOutFile):
    fig = plt.figure()
    plt.imshow(img,cmap=cmp, vmin=vmin,vmax=vmax)
    gca().set_xticklabels(['']*img.shape[0])
    gca().set_yticklabels(['']*img.shape[1])
    fig.gca().set_axis_off()
    fig.patch.set_facecolor('#ffffff')
    fig.patch.set_alpha(0)
    fig.savefig(imgOutFile, dpi=600, bbox_inches='tight', pad_inches=0, transparent=True)
    plt.close();
    return

def getImgDtStr(f):
    f_info = f.replace('PKL','INFO').replace('.pkl','.json')
    imgD = pklD.replace('/pkl/','/img/')
    infoFile = imgD + f_info
    with open(infoFile) as f:
        data = json.load(f)
        df = pd.DataFrame(data.items()).T
        df.columns = df.loc[0,:]
        df.drop(0,axis=0, inplace=True)
        df.reset_index(drop=True,inplace=True)
        
    imgDtStr = df.ImgDate[0]
    
    return imgDtStr

def addDateRateCols(results_df):
    dff = results_df.copy()
    
    list_of_dt_strs = results_df['1-Period-PDiff_DateRange'].tolist()
    i=0
    for dtRangeStr in list_of_dt_strs:
        startDtStr = dtRangeStr.split(' TO: ')[0]
        endDtStr = dtRangeStr.split(' TO: ')[1]
        startDtStrN = datetime.strftime(datetime.strptime(startDtStr,'%Y-%m-%d'),'%Y-%m-%dT%H:%M:%SZ')
        endDtStrN = datetime.strftime(datetime.strptime(endDtStr,'%Y-%m-%d'),'%Y-%m-%dT%H:%M:%SZ')
        dff.loc[i,'1-Period-PDiff_DateRange_Start'] = startDtStrN
        dff.loc[i,'1-Period-PDiff_DateRange_End'] = endDtStrN
        i=i+1
        
    list_of_dt_strs = results_df['4-Period-PDiff_DateRange'].tolist()
    i=0
    for dtRangeStr in list_of_dt_strs:
        startDtStr = dtRangeStr.split(' TO: ')[0]
        endDtStr = dtRangeStr.split(' TO: ')[1]
        startDtStrN = datetime.strftime(datetime.strptime(startDtStr,'%Y-%m-%d'),'%Y-%m-%dT%H:%M:%SZ')
        endDtStrN = datetime.strftime(datetime.strptime(endDtStr,'%Y-%m-%d'),'%Y-%m-%dT%H:%M:%SZ')
        dff.loc[i,'4-Period-PDiff_DateRange_Start'] = startDtStrN
        dff.loc[i,'4-Period-PDiff_DateRange_End'] = endDtStrN
        i=i+1
        
    return dff
############################ AUS
def getWeather4ForDates(p,list_of_dates,df):
    
    N=5 #default
    #gets the weather metrics for the METRIC calcs for a point for each of the dates based
    #of the closest N stns #given the geopandas geodataframe gpd
    
    #Screen out the null values
    df = df[df['ETo_average_inches']!='--']

    
    cols_list = ['Date',
                 'ETo_average_inches'
                ]
    df_weather_for_dates = pd.DataFrame(columns=cols_list)
    df.Date=pd.to_datetime(df.Date)
    

    idxx=0
    for dt in list_of_dates:
        df_dt = df[df.index==dt]
        
        df_dt['StationId']=df_dt['StationId'].astype('str')
        stnPoints = list(pd.unique(df_dt.geometry))
        stns = pd.unique(df_dt.StationId)

        [ids, pts, dis, wts] = getNearestNStnIDsInOrder2(p,stnPoints,5)

        ids_int = [int(x) for x in ids]
        nearestStnIDs = stns[ids_int]
        eto=[]
        for stn in nearestStnIDs:
            eto.append(float(df_dt[df_dt.StationId==stn]['ETo_average_inches'].values.tolist()[0]))
        
        weighted_eto = np.nansum(np.array(eto)*np.array(wts))
        
        df_weather_for_dates.loc[idxx,'Date']=dt
        df_weather_for_dates.loc[idxx,'ETo_average_inches']=weighted_eto
        idxx=idxx+1
          
    return df_weather_for_dates

def getWeather5ForDates(p,list_of_dates,df):
    
    N=5 #default
    #gets the weather metrics for the METRIC calcs for a point for each of the dates based
    #of the closest N stns #given the geopandas geodataframe gpd
    
    #Screen out the null values
    df = df[df['ETo_average_inches'] != '--']

    
    cols_list = ['Date',
                 'ETo_average_inches'
                ]
    df_weather_for_dates = pd.DataFrame(columns=cols_list)
    df.Date=pd.to_datetime(df.Date)
    
    idxx=0
    for dt in list_of_dates:
        df_dt = df[df.Date==dt]
        
        df_dt['StationId']=df_dt['StationId'].astype('str')
        stnPoints = list(pd.unique(df_dt.geometry))
        stns = pd.unique(df_dt.StationId)

        #print(df_dt)
        #print(stns)
        #print([p,stnPoints])

        [ids, pts, dis, wts] = getNearestNStnIDsInOrder2(p,stnPoints,5)

        #print([ids, pts, dis, wts])

        ids_int = [int(x) for x in ids]
        #print(ids_int)

        nearestStnIDs = stns[ids_int]

        eto=[]
        for stn in nearestStnIDs:
            eto.append(float(df_dt[df_dt.StationId==stn]['ETo_average_inches'].values.tolist()[0]))
        
        weighted_eto = np.nansum(np.array(eto)*np.array(wts))
        
        df_weather_for_dates.loc[idxx,'Date']=dt
        df_weather_for_dates.loc[idxx,'ETo_average_inches']=weighted_eto
        idxx=idxx+1
          
    return df_weather_for_dates


def getWeather6ForDates(p,list_of_dates,df):
    
    N=5 #default
    #gets the weather metrics for the METRIC calcs for a point for each of the dates based
    #of the closest N stns #given the geopandas geodataframe gpd
    
    #Screen out the null values
    df = df[df['ETo_average_inches'] != '--']

    
    cols_list = ['Date',
                 'ETo_average_inches'
                ]
    df_weather_for_dates = pd.DataFrame(columns=cols_list)
    df.Date=pd.to_datetime(df.Date)
    
    idxx=0
    for dt in list_of_dates:
        df_dt = df[df.Date==dt]
        
        df_dt['StationId']=df_dt['StationId'].astype('str')
        stnPoints = list(pd.unique(df_dt.geometry))
        stns = pd.unique(df_dt.StationId)

        #print(df_dt)
        #print(stns)
        #print([p,stnPoints])

        [ids, pts, dis, wts] = getNearestNStnIDsInOrder2(p,stnPoints,5)

        #print([ids, pts, dis, wts])

        ids_int = [int(x) for x in ids]
        #print(ids_int)

        nearestStnIDs = stns[ids_int]

        eto=[]
        for stn in nearestStnIDs:
            eto.append(float(df_dt[df_dt.StationId==stn]['ETo_average_inches'].values.tolist()[0]))
        
        weighted_eto = np.nansum(np.array(eto)*np.array(wts))
        
        df_weather_for_dates.loc[idxx,'Date']=dt
        df_weather_for_dates.loc[idxx,'ETo_average_inches']=weighted_eto
        idxx=idxx+1
          
    return df_weather_for_dates



def getWeather3ForDates(p,list_of_dates,df):
    
    N=5 #default
    #gets the weather metrics for the METRIC calcs for a point for each of the dates based
    #of the closest N stns #given the geopandas geodataframe gpd
    
    #Screen out the null values
    df = df[df['ETo_average_inches']!='--']

    
    cols_list = ['Date',
                 'ETo_average_inches'
                ]
    df_weather_for_dates = pd.DataFrame(columns=cols_list)
    df.Date=pd.to_datetime(df.Date)
    

    idxx=0
    for dt in list_of_dates:
        df_dt = df[df.index==dt]
        
        df_dt['StationId']=df_dt['StationId'].astype('int').astype('str')
        stnPoints = list(pd.unique(df_dt.geometry))
        stns = pd.unique(df_dt.StationId)

        [ids, pts, dis, wts] = getNearestNStnIDsInOrder2(p,stnPoints,5)

        ids_int = [int(x) for x in ids]
        nearestStnIDs = stns[ids_int]
        #eto=[]
        for stn in nearestStnIDs:
            eto.append(float(df_dt[df_dt.StationId==stn]['ETo_average_inches'].values.tolist()[0]))
        
        weighted_eto = np.nansum(np.array(eto)*np.array(wts))
        
        df_weather_for_dates.loc[idxx,'Date']=dt
        df_weather_for_dates.loc[idxx,'ETo_average_inches']=weighted_eto
        idxx=idxx+1
          
    return df_weather_for_dates

def generate_AUS_ETo(acct_gdf, stns_gdf, startDtStr, endDtStr):
    #print('hello')
    dailyCols = ['STN',
             'TS',
             'JulianDate',
             'ETo_average_inches',
             'Precip_in',
             'R_n__MJpm2',
             'P',
             'T_mean',
             'T_min',
             'T_max',
             'RH_mean',
             'RH_max',
             'RH_min',
             'T_dew',
             'Uz'
            ]
    dailyCols = [x.replace(' ','') for x in dailyCols]
    dailyCols = [x.replace('∞','') for x in dailyCols]
    dailyCols = [x.replace('=','') for x in dailyCols]

    dailyHistoricalWeatherPD = pd.DataFrame(columns=dailyCols)
    data = pd.DataFrame(columns = dailyCols)
    #read all the files into a big pandas dataframe
    dtFiles = []

    #years = ['2020']
    startYr = startDtStr.split('-')[0]
    endYr = endDtStr.split('-')[0]
    years = list(set([startYr,endYr]))
    years.sort()
    for yr in years:
        annualStnFiles = [f for f in os.listdir(aus_obsPath) if f.endswith('.geojson')]
        for f in annualStnFiles:
            try:
                y = f.split('__')[1].split('.')[0].startswith(yr)
                dtFiles.append(f)
            except:
                continue

    dtFiles.sort()
    #print(dtFiles)

    for fl in dtFiles:
        f = aus_obsPath + fl
        #print(f)
        data = gpd.read_file(f)
        #data.columns = dailyCols
        #print(str(len(data)))
        #fix the Date to make it a Datetime obj
        #################################################### 2021-01-12
        #data['Date'] =  pd.to_datetime(data['aifstime_utc'],format='%Y%m%d%H%M%S')
        data['aifstime_utc'] =  pd.to_datetime(data['aifstime_utc'],format='%Y%m%d%H%M%S')
        data['local_date_time_full'] =  pd.to_datetime(data['local_date_time_full'])
        #,format='%Y%m%d%H%M%S')
        data['Date'] = data.local_date_time_full

        #print(data.Date)
        #print(data)

        ####################################################
        dailyHistoricalWeatherPD = pd.concat([dailyHistoricalWeatherPD, data], ignore_index=True, sort=True)
        #print(str(len(dailyHistoricalWeatherPD)))

    dailyHistoricalWeatherPD['Date']=pd.to_datetime(dailyHistoricalWeatherPD['Date'])
    dailyHistoricalWeatherPD = dailyHistoricalWeatherPD.replace('-', np.nan)
    dailyHistoricalWeatherPD.reset_index(drop=True, inplace=True)
    float_cols = ['air_temp', 'dewpt', 'rain_trace', 'press','rel_hum','wind_spd_kmh', 'ETo_FAO_IN', 'ETo_HAR_IN', 'ETo_AVG_IN']
    dailyHistoricalWeatherPD[float_cols] = dailyHistoricalWeatherPD[float_cols].astype('float')

    #print(dailyHistoricalWeatherPD)

    #resample to get Daily data for each station
    weather_df = pd.DataFrame()
    weather_all_df = pd.DataFrame()
    list_of_unique_stns = list(set(dailyHistoricalWeatherPD['wmo']))
    i=0
    for stn in list_of_unique_stns:
        df = dailyHistoricalWeatherPD[dailyHistoricalWeatherPD.wmo==stn]
        df.set_index('Date', drop=True, inplace=True)
        if not np.isnan(int(stn)):
            weather_df['WMO']=int(stn)
            weather_df['T_mean__degC'] = df['air_temp'].resample('D').mean().squeeze()
            weather_df['T_max__degC'] = df['air_temp'].resample('D').max().squeeze()
            weather_df['T_min__degC'] = df['air_temp'].resample('D').min().squeeze()
            weather_df['RH_mean__pct'] = df['rel_hum'].resample('D').mean().squeeze()
            weather_df['RH_max__pct'] = df['rel_hum'].resample('D').max().squeeze()
            weather_df['RH_min__pct'] = df['rel_hum'].resample('D').min().squeeze()
            weather_df['T_dew__degC'] = df['dewpt'].resample('D').mean().squeeze()
            weather_df['Precip_in'] = df['rain_trace'].resample('D').sum().squeeze()*0.0393701 #convert to inches
            weather_df['P_mean__mbar'] = df['press'].resample('D').mean().squeeze()
            weather_df['WindSpeed_max__mps'] = df['wind_spd_kmh'].resample('D').max().squeeze()*0.277778
            weather_df['WindSpeed_min__mps'] = df['wind_spd_kmh'].resample('D').min().squeeze()*0.277778
            weather_df['WindSpeed_mean__mps'] = df['wind_spd_kmh'].resample('D').mean().squeeze()*0.277778
            weather_df['ETo_FAO_avg__in'] = df['ETo_FAO_IN'].resample('D').mean().squeeze()
            weather_df['ETo_HAR_avg__in'] = df['ETo_HAR_IN'].resample('D').mean().squeeze()
            weather_df['ETo_average_inches'] = df['ETo_AVG_IN'].resample('D').mean().squeeze()
            weather_df['Date']=datetime.strftime(df.index[0],'%Y-%m-%d')
            weather_df['Lat']=df['lat'].resample('D').mean().squeeze()
            weather_df['Lon']=df['lon'].resample('D').mean().squeeze()

            weather_all_df =  weather_all_df.append(weather_df)
    geom = [Point(xy) for xy in zip(weather_all_df['Lon'], weather_all_df['Lat'])]
    # Coordinate reference system : WGS84
    crs = {'init': 'epsg:4326'}
    # Creating a Geographic data frame 
    weather_gdf = gpd.GeoDataFrame(weather_all_df, crs=crs, geometry=geom)
    weather_gdf = weather_gdf.dropna()

    #print(weather_gdf)

    stns = pd.unique(weather_gdf.WMO).tolist()
    stns.sort()
    stnPoints=[]
    for stn in stns:
        stnStr = str(int(stn))
        p = stns_gdf[stns_gdf.WMO==stnStr].geometry.tolist()[0]
        stnPoints.append(p)

    dmw_gpd = weather_gdf.copy()

    #Screen out the null values
    dmw_gpd['ETo_average_inches'] = dmw_gpd['ETo_average_inches'].replace('--',np.nan, regex=True).astype('float')
    dmw_gpd['T_mean__degC'] = dmw_gpd['T_mean__degC'].replace('--',np.nan, regex=True).astype('float')
    #dmw_gpd['R_n__MJpm2'] = dmw_gpd['R_n__MJpm2'].replace('--',np.nan, regex=True).astype('float')
    dmw_gpd['WindSpeed_mean__mps'] = dmw_gpd['WindSpeed_mean__mps'].replace('--',np.nan, regex=True).astype('float')

    #Put values for metric system calculations
    dmw_gpd['Tavg']=dmw_gpd['T_mean__degC'].astype('float') 
    #dmw_gpd['Rnet']=dmw_gpd['R_n__MJpm2'].values.astype('float') * 11.574 #Convert to Wm-2
    dmw_gpd['Wavg']=dmw_gpd['WindSpeed_mean__mps'].astype('float')
    #dmw_gpd['Date']=dmw_gpd['TS']
    dmw_gpd['StationId']=dmw_gpd['WMO']
    dmw_gpd['T_mean']=dmw_gpd['T_mean__degC'].astype('float')
    dmw_gpd['T_min']=dmw_gpd['T_min__degC'].astype('float')
    dmw_gpd['T_max']=dmw_gpd['T_max__degC'].astype('float')

    dmw_gpd['RH_mean']=dmw_gpd['RH_mean__pct'].astype('float')
    dmw_gpd['RH_min']=dmw_gpd['RH_min__pct'].astype('float')
    dmw_gpd['RH_max']=dmw_gpd['RH_max__pct'].astype('float')
    dmw_gpd['T_dew']=dmw_gpd['T_dew__degC'].astype('float')
    #
    weather_gdf = dmw_gpd.copy()
    keep_metric_cols_list = ['Date','StationId','Tavg','Wavg']
    #keep_metric_cols_list = ['Date','StationId','Tavg','Rnet','Wavg']
    metric = dmw_gpd[keep_metric_cols_list] 


    #Set the horizon from the inputs
    startDt = datetime.strptime(startDtStr, '%Y-%m-%d')
    endDt = datetime.strptime(endDtStr,'%Y-%m-%d')

    #print([startDt, endDt])

    list_of_dates = [datetime.strptime(d.strftime('%Y-%m-%d'),'%Y-%m-%d') for d in pd.date_range(startDt, endDt)]

    startt = time.time()
    EToForAccounts = pd.DataFrame()
    for i,row in acct_gdf.iterrows():
        try:
            print('\t\tCalculating ETo for'+
                ' Customerid: '+str(row.customerid)+
                ', Assetid: '+str(row.fieldassetid)+
                ', Acct: ' + row.Category+
                ', Ranch: '+row.SubCategory+
                ', Field: '+row.Name)
            start_time = time.time()
            p = loads(row.CentroidPoint) #shapely Point obj
            
            etv = getWeather3ForDates(p,list_of_dates,weather_gdf)
            
            #now add the columns for the other fields to etv
            etv['UUID']=row.UUID
            etv['Category']=row.Category
            etv['SubCategory']=row.SubCategory
            etv['Name']=row.Name

            EToForAccounts = EToForAccounts.append(etv, ignore_index = True)
            print("\t\tEND ETo--- %s seconds ---" % (time.time() - start_time))
        except:
            continue
    print("\tETo Done --- %s seconds ---" % (time.time() - startt))
    print('\tDone!')

    EToForAccounts = EToForAccounts.rename(columns={"ETo_average_inches": "ETo(in)"})


    #Merge the acct
    merged_df = pd.DataFrame()
    if len(EToForAccounts)>0:
        merged_df = pd.merge(acct_gdf, EToForAccounts, how='inner', on=['UUID'])
        #on=['Category','SubCategory','Name', 'UUID'])
        merged_df.reset_index() #reset the index to avoid any export issues
    
    return merged_df


def generate_ISR_ETo2(acct_gdf, stns_gdf, startDtStr, endDtStr):
    #print('hello')
    dailyCols = ['STN',
             'TS',
             'JulianDate',
             'ETo_average_inches',
            ]
    dailyCols = [x.replace(' ','') for x in dailyCols]
    dailyCols = [x.replace('∞','') for x in dailyCols]
    dailyCols = [x.replace('=','') for x in dailyCols]

    dailyHistoricalWeatherPD = pd.DataFrame(columns=dailyCols)
    data = pd.DataFrame(columns = dailyCols)
    dtFiles = []

    startYr = startDtStr.split('-')[0]
    endYr = endDtStr.split('-')[0]
    years = list(set([startYr,endYr]))
    years.sort()


    for yr in years:
        annualStnFiles = [f for f in os.listdir(isr_obsPath) if (f.startswith('ISR') & f.endswith('.geojson'))]
        for f in annualStnFiles:
            try:
                y = f.split('__')[1].split('.')[0].startswith(yr)
                dtFiles.append(f)
            except:
                continue
    dtFiles.sort()

    for fl in dtFiles:
        f = isr_obsPath_new + fl
        data = gpd.read_file(f)
        #################################################### 2021-01-12
        data['aifstime_utc'] =  pd.to_datetime(data['aifstime_utc'],format='%Y%m%d%H%M%S')
        data['local_date_time_full'] =  pd.to_datetime(data['local_date_time_full'])
        data['Date'] = data.local_date_time_full

        ####################################################
        dailyHistoricalWeatherPD = pd.concat([dailyHistoricalWeatherPD, data], ignore_index=True, sort=True)

    dailyHistoricalWeatherPD['Date']=pd.to_datetime(dailyHistoricalWeatherPD['Date'], format='%Y%m%d%H%M%S', errors='ignore')

    dailyHistoricalWeatherPD = dailyHistoricalWeatherPD.replace('-', np.nan)
    dailyHistoricalWeatherPD.reset_index(drop=True, inplace=True)
    float_cols = ['PRCP','TAVG','TMAX','TMIN','ETo_FAO_IN', 'ETo_HAR_IN', 'ETo_AVG_IN']
    dailyHistoricalWeatherPD[float_cols] = dailyHistoricalWeatherPD[float_cols].astype('float')

    #print(dailyHistoricalWeatherPD)

    #resample to get Daily data for each station
    weather_df = pd.DataFrame()
    weather_all_df = pd.DataFrame()
    list_of_unique_stns = list(set(dailyHistoricalWeatherPD['StnID']))

    #print(list_of_unique_stns)

    i=0
    for stn in list_of_unique_stns:
        df = dailyHistoricalWeatherPD[dailyHistoricalWeatherPD.StnID==stn]
        df.set_index('Date', drop=True, inplace=True)
        
        weather_df['StnID']=stn
        weather_df['T_mean__degC'] = df['TAVG'].resample('D').mean().squeeze()
        weather_df['T_max__degC'] = df['TMAX'].resample('D').max().squeeze()
        weather_df['T_min__degC'] = df['TMIN'].resample('D').min().squeeze()
        weather_df['Precip_in'] = df['PRCP'].resample('D').sum().squeeze()*0.0393701 #convert to inches
        weather_df['ETo_FAO_avg__in'] = df['ETo_FAO_IN'].resample('D').mean().squeeze()
        weather_df['ETo_HAR_avg__in'] = df['ETo_HAR_IN'].resample('D').mean().squeeze()
        weather_df['ETo_average_inches'] = df['ETo_AVG_IN'].resample('D').mean().squeeze()
        weather_df['Date']=datetime.strftime(df.index[0],'%Y-%m-%d')
        weather_df['Lat']=df['lat'].resample('D').mean().squeeze()
        weather_df['Lon']=df['lon'].resample('D').mean().squeeze()

        weather_all_df =  weather_all_df.append(weather_df)

    #print(weather_all_df.columns)
    geom = [Point(xy) for xy in zip(weather_all_df['Lon'], weather_all_df['Lat'])]
    # Coordinate reference system : WGS84
    crs = {'init': 'epsg:4326'}
    # Creating a Geographic data frame 
    weather_gdf = gpd.GeoDataFrame(weather_all_df, crs=crs, geometry=geom)
    weather_gdf = weather_gdf.dropna()

    #print(weather_gdf)

    stns = pd.unique(weather_gdf.StnID).tolist()
    stns.sort()
    #print(stns)
    stnPoints=[]

    #print(stns_gdf.columns)
    #print(len(stns_gdf))
    for stn in stns:
        stnStr = str(stn)
        p = stns_gdf[stns_gdf.stnid==stnStr].geometry
        stnPoints.append(p)

    dmw_gpd = weather_gdf.copy()

    #Screen out the null values
    dmw_gpd['ETo_average_inches'] = dmw_gpd['ETo_average_inches'].replace('--',np.nan, regex=True).astype('float')
    dmw_gpd['StationId']=dmw_gpd['StnID']
    dmw_gpd['T_mean__degC'] = dmw_gpd['T_mean__degC'].replace('--',np.nan, regex=True).astype('float')	
    
    #Put values for metric system calculations
    dmw_gpd['Tavg']=dmw_gpd['T_mean__degC'].astype('float') 
    dmw_gpd['Wavg']=0.0

    dmw_gpd['T_mean']=dmw_gpd['T_mean__degC'].astype('float')
    dmw_gpd['T_min']=dmw_gpd['T_min__degC'].astype('float')
    dmw_gpd['T_max']=dmw_gpd['T_max__degC'].astype('float')

    #
    weather_gdf = dmw_gpd.copy()
    keep_metric_cols_list = ['Date','StationId','Tavg','Wavg']
    metric = dmw_gpd[keep_metric_cols_list] 

    #Set the horizon from the inputs
    startDt = datetime.strptime(startDtStr, '%Y-%m-%d')
    endDt = datetime.strptime(endDtStr,'%Y-%m-%d')

    #print([startDt, endDt])

    list_of_dates = [datetime.strptime(d.strftime('%Y-%m-%d'),'%Y-%m-%d') for d in pd.date_range(startDt, endDt)]

    startt = time.time()
    EToForAccounts = pd.DataFrame()
    for i,row in acct_gdf.iterrows():
        print('\t\tCalculating ETo for'+
            ', Acct: ' + row.Category+
            ', Ranch: '+row.SubCategory+
            ', Field: '+row.Name)
        start_time = time.time()
        p = loads(row.CentroidPoint) #shapely Point obj
        
        etv = getWeather4ForDates(p,list_of_dates,weather_gdf)

        #now add the columns for the other fields to etv
        etv['UUID']=row.UUID
        etv['Category']=row.Category
        etv['SubCategory']=row.SubCategory
        etv['Name']=row.Name

        EToForAccounts = EToForAccounts.append(etv, ignore_index = True)
        print("\t\tEND ETo--- %s seconds ---" % (time.time() - start_time))
    print("\tETo Done --- %s seconds ---" % (time.time() - startt))
    print('\tDone!')

    EToForAccounts = EToForAccounts.rename(columns={"ETo_average_inches": "ETo(in)"})


    #Merge the acct
    merged_df = pd.merge(acct_gdf, EToForAccounts, on=['Category','SubCategory','Name', 'UUID'])
    merged_df.reset_index() #reset the index to avoid any export issues
    
    return merged_df

def generate_EGY_ETo2(acct_gdf, stns_gdf, startDtStr, endDtStr):
    #print('hello')
    dailyCols = ['STN',
             'TS',
             'JulianDate',
             'ETo_average_inches',
            ]
    dailyCols = [x.replace(' ','') for x in dailyCols]
    dailyCols = [x.replace('∞','') for x in dailyCols]
    dailyCols = [x.replace('=','') for x in dailyCols]

    dailyHistoricalWeatherPD = pd.DataFrame(columns=dailyCols)
    data = pd.DataFrame(columns = dailyCols)
    dtFiles = []

    startYr = startDtStr.split('-')[0]
    endYr = endDtStr.split('-')[0]
    years = list(set([startYr,endYr]))
    years.sort()


    for yr in years:
        annualStnFiles = [f for f in os.listdir(isr_obsPath) if (f.startswith('EGY') & f.endswith('.geojson'))]
        for f in annualStnFiles:
            try:
                y = f.split('__')[1].split('.')[0].startswith(yr)
                dtFiles.append(f)
            except:
                continue
    dtFiles.sort()
    #print(annualStnFiles)
    annualStnFiles.sort()

    for fl in dtFiles:
        f = egy_obsPath + fl
        data = gpd.read_file(f)
        #################################################### 2021-01-12
        data['aifstime_utc'] =  pd.to_datetime(data['aifstime_utc'],format='%Y%m%d%H%M%S')
        data['local_date_time_full'] =  pd.to_datetime(data['local_date_time_full'])
        data['Date'] = data.local_date_time_full

        ####################################################
        dailyHistoricalWeatherPD = pd.concat([dailyHistoricalWeatherPD, data], ignore_index=True, sort=True)

    dailyHistoricalWeatherPD['Date']=pd.to_datetime(dailyHistoricalWeatherPD['Date'], format='%Y%m%d%H%M%S', errors='ignore')

    dailyHistoricalWeatherPD = dailyHistoricalWeatherPD.replace('-', np.nan)
    dailyHistoricalWeatherPD.reset_index(drop=True, inplace=True)
    float_cols = ['PRCP','TAVG','TMAX','TMIN','ETo_FAO_IN', 'ETo_HAR_IN', 'ETo_AVG_IN']
    dailyHistoricalWeatherPD[float_cols] = dailyHistoricalWeatherPD[float_cols].astype('float')

    #print(dailyHistoricalWeatherPD)

    #resample to get Daily data for each station
    weather_df = pd.DataFrame()
    weather_all_df = pd.DataFrame()
    list_of_unique_stns = list(set(dailyHistoricalWeatherPD['StnID']))

    #print(list_of_unique_stns)

    i=0
    for stn in list_of_unique_stns:
        df = dailyHistoricalWeatherPD[dailyHistoricalWeatherPD.StnID==stn]
        df.set_index('Date', drop=True, inplace=True)
        
        weather_df['StnID']=stn
        weather_df['T_mean__degC'] = df['TAVG'].resample('D').mean().squeeze()
        weather_df['T_max__degC'] = df['TMAX'].resample('D').max().squeeze()
        weather_df['T_min__degC'] = df['TMIN'].resample('D').min().squeeze()
        weather_df['Precip_in'] = df['PRCP'].resample('D').sum().squeeze()*0.0393701 #convert to inches
        weather_df['ETo_FAO_avg__in'] = df['ETo_FAO_IN'].resample('D').mean().squeeze()
        weather_df['ETo_HAR_avg__in'] = df['ETo_HAR_IN'].resample('D').mean().squeeze()
        weather_df['ETo_average_inches'] = df['ETo_AVG_IN'].resample('D').mean().squeeze()
        weather_df['Date']=datetime.strftime(df.index[0],'%Y-%m-%d')
        weather_df['Lat']=df['lat'].resample('D').mean().squeeze()
        weather_df['Lon']=df['lon'].resample('D').mean().squeeze()

        weather_all_df =  weather_all_df.append(weather_df)

    #print(weather_all_df.columns)
    geom = [Point(xy) for xy in zip(weather_all_df['Lon'], weather_all_df['Lat'])]
    # Coordinate reference system : WGS84
    crs = {'init': 'epsg:4326'}
    # Creating a Geographic data frame 
    weather_gdf = gpd.GeoDataFrame(weather_all_df, crs=crs, geometry=geom)
    weather_gdf = weather_gdf.dropna()

    #print(weather_gdf)

    stns = pd.unique(weather_gdf.StnID).tolist()
    stns.sort()
    #print(stns)
    stnPoints=[]

    #print(stns_gdf.columns)
    #print(len(stns_gdf))
    for stn in stns:
        stnStr = str(stn)
        p = stns_gdf[stns_gdf.stnid==stnStr].geometry
        stnPoints.append(p)

    dmw_gpd = weather_gdf.copy()

    #Screen out the null values
    dmw_gpd['ETo_average_inches'] = dmw_gpd['ETo_average_inches'].replace('--',np.nan, regex=True).astype('float')
    dmw_gpd['StationId']=dmw_gpd['StnID']
    dmw_gpd['T_mean__degC'] = dmw_gpd['T_mean__degC'].replace('--',np.nan, regex=True).astype('float')	
    
    #Put values for metric system calculations
    dmw_gpd['Tavg']=dmw_gpd['T_mean__degC'].astype('float') 
    dmw_gpd['Wavg']=0.0

    dmw_gpd['T_mean']=dmw_gpd['T_mean__degC'].astype('float')
    dmw_gpd['T_min']=dmw_gpd['T_min__degC'].astype('float')
    dmw_gpd['T_max']=dmw_gpd['T_max__degC'].astype('float')

    #
    weather_gdf = dmw_gpd.copy()
    keep_metric_cols_list = ['Date','StationId','Tavg','Wavg']
    metric = dmw_gpd[keep_metric_cols_list] 

    #Set the horizon from the inputs
    startDt = datetime.strptime(startDtStr, '%Y-%m-%d')
    endDt = datetime.strptime(endDtStr,'%Y-%m-%d')

    #print([startDt, endDt])

    list_of_dates = [datetime.strptime(d.strftime('%Y-%m-%d'),'%Y-%m-%d') for d in pd.date_range(startDt, endDt)]

    startt = time.time()
    EToForAccounts = pd.DataFrame()
    for i,row in acct_gdf.iterrows():
        print('\t\tCalculating ETo for'+
            ', Acct: ' + row.Category+
            ', Ranch: '+row.SubCategory+
            ', Field: '+row.Name)
        start_time = time.time()
        p = loads(row.CentroidPoint) #shapely Point obj
        
        etv = getWeather4ForDates(p,list_of_dates,weather_gdf)

        #now add the columns for the other fields to etv
        etv['UUID']=row.UUID
        etv['Category']=row.Category
        etv['SubCategory']=row.SubCategory
        etv['Name']=row.Name

        EToForAccounts = EToForAccounts.append(etv, ignore_index = True)
        print("\t\tEND ETo--- %s seconds ---" % (time.time() - start_time))
    print("\tETo Done --- %s seconds ---" % (time.time() - startt))
    print('\tDone!')

    EToForAccounts = EToForAccounts.rename(columns={"ETo_average_inches": "ETo(in)"})


    #Merge the acct
    merged_df = pd.merge(acct_gdf, EToForAccounts, on=['Category','SubCategory','Name', 'UUID'])
    merged_df.reset_index() #reset the index to avoid any export issues
    
    return merged_df

def generate_GHCND_ETo(acct_gdf, stns_gdf, startDtStr, endDtStr, region):
    #print('hello')
    dailyCols = ['STN',
             'TS',
             'JulianDate',
             'ETo_average_inches',
            ]
    dailyCols = [x.replace(' ','') for x in dailyCols]
    dailyCols = [x.replace('∞','') for x in dailyCols]
    dailyCols = [x.replace('=','') for x in dailyCols]

    dailyHistoricalWeatherPD = pd.DataFrame(columns=dailyCols)
    data = pd.DataFrame(columns = dailyCols)
    dtFiles = []

    startYr = startDtStr.split('-')[0]
    endYr = endDtStr.split('-')[0]
    years = list(set([startYr,endYr]))
    years.sort()


    for yr in years:
        annualStnFiles = [f for f in os.listdir(isr_obsPath) if (f.startswith(region) & f.endswith('.geojson'))]
        for f in annualStnFiles:
            try:
                y = f.split('__')[1].split('.')[0].startswith(yr)
                dtFiles.append(f)
            except:
                continue
    dtFiles.sort()
    #print(annualStnFiles)
    annualStnFiles.sort()

    for fl in dtFiles:
        f = egy_obsPath + fl
        data = gpd.read_file(f)
        #################################################### 2021-01-12
        data['aifstime_utc'] =  pd.to_datetime(data['aifstime_utc'],format='%Y%m%d%H%M%S')
        data['local_date_time_full'] =  pd.to_datetime(data['local_date_time_full'])
        data['Date'] = data.local_date_time_full

        ####################################################
        dailyHistoricalWeatherPD = pd.concat([dailyHistoricalWeatherPD, data], ignore_index=True, sort=True)

    dailyHistoricalWeatherPD['Date']=pd.to_datetime(dailyHistoricalWeatherPD['Date'], format='%Y%m%d%H%M%S', errors='ignore')

    dailyHistoricalWeatherPD = dailyHistoricalWeatherPD.replace('-', np.nan)
    dailyHistoricalWeatherPD.reset_index(drop=True, inplace=True)
    float_cols = ['PRCP','TAVG','TMAX','TMIN','ETo_FAO_IN', 'ETo_HAR_IN', 'ETo_AVG_IN']
    dailyHistoricalWeatherPD[float_cols] = dailyHistoricalWeatherPD[float_cols].astype('float')

    #print(dailyHistoricalWeatherPD)

    #resample to get Daily data for each station
    weather_df = pd.DataFrame()
    weather_all_df = pd.DataFrame()
    list_of_unique_stns = list(set(dailyHistoricalWeatherPD['StnID']))

    #print(list_of_unique_stns)

    i=0
    for stn in list_of_unique_stns:
        df = dailyHistoricalWeatherPD[dailyHistoricalWeatherPD.StnID==stn]
        df.set_index('Date', drop=True, inplace=True)
        
        weather_df['StnID']=stn
        weather_df['T_mean__degC'] = df['TAVG'].resample('D').mean().squeeze()
        weather_df['T_max__degC'] = df['TMAX'].resample('D').max().squeeze()
        weather_df['T_min__degC'] = df['TMIN'].resample('D').min().squeeze()
        weather_df['Precip_in'] = df['PRCP'].resample('D').sum().squeeze()*0.0393701 #convert to inches
        weather_df['ETo_FAO_avg__in'] = df['ETo_FAO_IN'].resample('D').mean().squeeze()
        weather_df['ETo_HAR_avg__in'] = df['ETo_HAR_IN'].resample('D').mean().squeeze()
        weather_df['ETo_average_inches'] = df['ETo_AVG_IN'].resample('D').mean().squeeze()
        weather_df['Date']=datetime.strftime(df.index[0],'%Y-%m-%d')
        weather_df['Lat']=df['lat'].resample('D').mean().squeeze()
        weather_df['Lon']=df['lon'].resample('D').mean().squeeze()

        weather_all_df =  weather_all_df.append(weather_df)

    #print(weather_all_df.columns)
    geom = [Point(xy) for xy in zip(weather_all_df['Lon'], weather_all_df['Lat'])]
    # Coordinate reference system : WGS84
    crs = {'init': 'epsg:4326'}
    # Creating a Geographic data frame 
    weather_gdf = gpd.GeoDataFrame(weather_all_df, crs=crs, geometry=geom)
    weather_gdf = weather_gdf.dropna()

    #print(weather_gdf)

    stns = pd.unique(weather_gdf.StnID).tolist()
    stns.sort()
    #print(stns)
    stnPoints=[]

    #print(stns_gdf.columns)
    #print(len(stns_gdf))
    for stn in stns:
        stnStr = str(stn)
        p = stns_gdf[stns_gdf.stnid==stnStr].geometry
        stnPoints.append(p)

    dmw_gpd = weather_gdf.copy()

    #Screen out the null values
    dmw_gpd['ETo_average_inches'] = dmw_gpd['ETo_average_inches'].replace('--',np.nan, regex=True).astype('float')
    dmw_gpd['StationId']=dmw_gpd['StnID']
    dmw_gpd['T_mean__degC'] = dmw_gpd['T_mean__degC'].replace('--',np.nan, regex=True).astype('float')	
    
    #Put values for metric system calculations
    dmw_gpd['Tavg']=dmw_gpd['T_mean__degC'].astype('float') 
    dmw_gpd['Wavg']=0.0

    dmw_gpd['T_mean']=dmw_gpd['T_mean__degC'].astype('float')
    dmw_gpd['T_min']=dmw_gpd['T_min__degC'].astype('float')
    dmw_gpd['T_max']=dmw_gpd['T_max__degC'].astype('float')

    #
    weather_gdf = dmw_gpd.copy()
    keep_metric_cols_list = ['Date','StationId','Tavg','Wavg']
    metric = dmw_gpd[keep_metric_cols_list] 

    #Set the horizon from the inputs
    startDt = datetime.strptime(startDtStr, '%Y-%m-%d')
    endDt = datetime.strptime(endDtStr,'%Y-%m-%d')

    #print([startDt, endDt])

    list_of_dates = [datetime.strptime(d.strftime('%Y-%m-%d'),'%Y-%m-%d') for d in pd.date_range(startDt, endDt)]

    startt = time.time()
    EToForAccounts = pd.DataFrame()
    for i,row in acct_gdf.iterrows():
        print('\t\tCalculating ETo for'+
            ', Acct: ' + row.Category+
            ', Ranch: '+row.SubCategory+
            ', Field: '+row.Name)
        start_time = time.time()
        p = loads(row.CentroidPoint) #shapely Point obj
        
        etv = getWeather4ForDates(p,list_of_dates,weather_gdf)

        #now add the columns for the other fields to etv
        etv['UUID']=row.UUID
        etv['Category']=row.Category
        etv['SubCategory']=row.SubCategory
        etv['Name']=row.Name

        EToForAccounts = EToForAccounts.append(etv, ignore_index = True)
        print("\t\tEND ETo--- %s seconds ---" % (time.time() - start_time))
    print("\tETo Done --- %s seconds ---" % (time.time() - startt))
    print('\tDone!')

    EToForAccounts = EToForAccounts.rename(columns={"ETo_average_inches": "ETo(in)"})


    #Merge the acct
    merged_df = pd.merge(acct_gdf, EToForAccounts, on=['Category','SubCategory','Name', 'UUID'])
    merged_df.reset_index() #reset the index to avoid any export issues
    
    return merged_df

def generate_AUS_ETo2(acct_gdf, stns_gdf, startDtStr, endDtStr):
    #print('hello')
    dailyCols = ['STN',
             'TS',
             'JulianDate',
             'ETo_average_inches',
            ]
    dailyCols = [x.replace(' ','') for x in dailyCols]
    dailyCols = [x.replace('∞','') for x in dailyCols]
    dailyCols = [x.replace('=','') for x in dailyCols]

    dailyHistoricalWeatherPD = pd.DataFrame(columns=dailyCols)
    data = pd.DataFrame(columns = dailyCols)
    dtFiles = []

    startYr = startDtStr.split('-')[0]
    endYr = endDtStr.split('-')[0]
    years = list(set([startYr,endYr]))
    years.sort()


    for yr in years:
        annualStnFiles = [f for f in os.listdir(aus_obsPath_new) if (f.startswith('AUS') & f.endswith('.geojson'))]
        for f in annualStnFiles:
            try:
                y = f.split('__')[1].split('.')[0].startswith(yr)
                dtFiles.append(f)
            except:
                continue
    dtFiles.sort()

    for fl in dtFiles:
        f = aus_obsPath_new + fl
        data = gpd.read_file(f)
        #################################################### 2021-01-12
        data['aifstime_utc'] =  pd.to_datetime(data['aifstime_utc'],format='%Y%m%d%H%M%S')
        data['local_date_time_full'] =  pd.to_datetime(data['local_date_time_full'])
        data['Date'] = data.local_date_time_full

        ####################################################
        dailyHistoricalWeatherPD = pd.concat([dailyHistoricalWeatherPD, data], ignore_index=True, sort=True)

    #print(dailyHistoricalWeatherPD)
    #print(dailyHistoricalWeatherPD['Date'])
    #print(type(dailyHistoricalWeatherPD['Date']))

    dailyHistoricalWeatherPD['Date']=pd.to_datetime(dailyHistoricalWeatherPD['Date'], format='%Y%m%d%H%M%S', errors='ignore')
    #infer_datetime_format=True)
    #

    dailyHistoricalWeatherPD = dailyHistoricalWeatherPD.replace('-', np.nan)
    dailyHistoricalWeatherPD.reset_index(drop=True, inplace=True)
    float_cols = ['PRCP','TAVG','TMAX','TMIN','ETo_FAO_IN', 'ETo_HAR_IN', 'ETo_AVG_IN']
    dailyHistoricalWeatherPD[float_cols] = dailyHistoricalWeatherPD[float_cols].astype('float')

    #print(dailyHistoricalWeatherPD)

    #resample to get Daily data for each station
    weather_df = pd.DataFrame()
    weather_all_df = pd.DataFrame()
    list_of_unique_stns = list(set(dailyHistoricalWeatherPD['StnID']))

    #print(list_of_unique_stns)

    i=0
    for stn in list_of_unique_stns:
        df = dailyHistoricalWeatherPD[dailyHistoricalWeatherPD.StnID==stn]
        df.set_index('Date', drop=True, inplace=True)
        #
        #if (stn):
        #print(df)
        #print(df.columns)
        #print(len(weather_df))
        #(str(stn).replace(' ','').strip()):
        weather_df['StnID']=stn
        weather_df['T_mean__degC'] = df['TAVG'].resample('D').mean().squeeze()
        weather_df['T_max__degC'] = df['TMAX'].resample('D').max().squeeze()
        weather_df['T_min__degC'] = df['TMIN'].resample('D').min().squeeze()
        weather_df['Precip_in'] = df['PRCP'].resample('D').sum().squeeze()*0.0393701 #convert to inches

        weather_df['ETo_FAO_avg__in'] = df['ETo_FAO_IN'].resample('D').mean().squeeze()
        weather_df['ETo_HAR_avg__in'] = df['ETo_HAR_IN'].resample('D').mean().squeeze()
        weather_df['ETo_average_inches'] = df['ETo_AVG_IN'].resample('D').mean().squeeze()
        weather_df['Date']=datetime.strftime(df.index[0],'%Y-%m-%d')
        weather_df['Lat']=df['lat'].resample('D').mean().squeeze()
        weather_df['Lon']=df['lon'].resample('D').mean().squeeze()

        weather_all_df =  weather_all_df.append(weather_df)

    #print(weather_all_df.columns)
    geom = [Point(xy) for xy in zip(weather_all_df['Lon'], weather_all_df['Lat'])]
    # Coordinate reference system : WGS84
    crs = {'init': 'epsg:4326'}
    # Creating a Geographic data frame 
    weather_gdf = gpd.GeoDataFrame(weather_all_df, crs=crs, geometry=geom)
    weather_gdf = weather_gdf.dropna()

    #print(weather_gdf)

    stns = pd.unique(weather_gdf.StnID).tolist()
    stns.sort()
    stnPoints=[]

    #print(stns_gdf.columns)
    #print(len(stns_gdf))
    for stn in stns:
        stnStr = str(stn)
        p = stns_gdf[stns_gdf.stnid==stn].geometry.tolist()[0]
        stnPoints.append(p)

    dmw_gpd = weather_gdf.copy()

    #Screen out the null values
    dmw_gpd['ETo_average_inches'] = dmw_gpd['ETo_average_inches'].replace('--',np.nan, regex=True).astype('float')
    dmw_gpd['StationId']=dmw_gpd['StnID']
    dmw_gpd['T_mean__degC'] = dmw_gpd['T_mean__degC'].replace('--',np.nan, regex=True).astype('float')	
    
    #Put values for metric system calculations
    dmw_gpd['Tavg']=dmw_gpd['T_mean__degC'].astype('float') 
    dmw_gpd['Wavg']=0.0

    dmw_gpd['T_mean']=dmw_gpd['T_mean__degC'].astype('float')
    dmw_gpd['T_min']=dmw_gpd['T_min__degC'].astype('float')
    dmw_gpd['T_max']=dmw_gpd['T_max__degC'].astype('float')

    #
    weather_gdf = dmw_gpd.copy()
    keep_metric_cols_list = ['Date','StationId','Tavg','Wavg']
    metric = dmw_gpd[keep_metric_cols_list] 

    #Set the horizon from the inputs
    startDt = datetime.strptime(startDtStr, '%Y-%m-%d')
    endDt = datetime.strptime(endDtStr,'%Y-%m-%d')

    #print([startDt, endDt])

    list_of_dates = [datetime.strptime(d.strftime('%Y-%m-%d'),'%Y-%m-%d') for d in pd.date_range(startDt, endDt)]

    startt = time.time()
    EToForAccounts = pd.DataFrame()
    for i,row in acct_gdf.iterrows():
        print('\t\tCalculating ETo for'+
            ' Customerid: '+str(row.customerid)+
            ', Assetid: '+str(row.fieldassetid)+
            ', Acct: ' + row.Category+
            ', Ranch: '+row.SubCategory+
            ', Field: '+row.Name)
        start_time = time.time()
        p = loads(row.CentroidPoint) #shapely Point obj
        
        etv = getWeather4ForDates(p,list_of_dates,weather_gdf)
        #print('etv in generate_AUS_ETo2 = '+str(etv))

        #now add the columns for the other fields to etv
        etv['UUID']=row.UUID
        etv['Category']=row.Category
        etv['SubCategory']=row.SubCategory
        etv['Name']=row.Name

        EToForAccounts = EToForAccounts.append(etv, ignore_index = True)
        print("\t\tEND ETo--- %s seconds ---" % (time.time() - start_time))
    print("\tETo Done --- %s seconds ---" % (time.time() - startt))
    print('\tDone!')

    EToForAccounts = EToForAccounts.rename(columns={"ETo_average_inches": "ETo(in)"})


    #Merge the acct
    merged_df = pd.merge(acct_gdf, EToForAccounts, on=['Category','SubCategory','Name', 'UUID'])
    merged_df.reset_index() #reset the index to avoid any export issues
    
    return merged_df

def getForecastedAUSEToForPtAndDtStr(dtStr,p):
    try:
        forecastDir = '/mnt/md0/AUS/FORECAST/'
        jsonFile = 'AUS_ETO__'+dtStr+'.geojson'
        forecast_gdf = gpd.read_file(forecastDir+jsonFile)
        
        for i,row in forecast_gdf.iterrows():
            forecast_gdf.loc[i,'DIST']=p.distance(row.geometry)
        forecast_gdf = forecast_gdf.sort_values(by=['DIST'])
        forecast_gdf.reset_index(drop=True, inplace=True)
        eto_forecast = forecast_gdf.ETo_average_inches[0]
        #print('\nGetting ETo '+dtStr+' from Forecast')
    except:
        eto_forecast = np.nan

    return eto_forecast

def getDailySummaryAUSEToForPtAndDtStr(dtStr,p):
    try:
        dailySummaryDir = '/mnt/md0/GHCND/DAILY/DAILYSUMMARY_PROCESSED/'
        jsonFile = 'AUS_ETO__'+dtStr+'.geojson'
        dailysum_gdf = gpd.read_file(dailySummaryDir+jsonFile)

        for i,row in dailysum_gdf.iterrows():
            dailysum_gdf.loc[i,'DIST']=p.distance(row.geometry)
        dailysum_gdf = dailysum_gdf.sort_values(by=['DIST'])
        dailysum_gdf.reset_index(drop=True, inplace=True)
        eto_dailysummary = dailysum_gdf.ETo_average_inches[0]
        #print('\nGetting ETo for '+dtStr+' from Daily Summary')
    except:
        eto_dailysummary = np.nan
    return eto_dailysummary


def generate_AUS_ETo3(acct_gdf, stns_gdf, startDtStr, endDtStr):
    #Insert logic for computing the ETO from forecasted data in case the obs data is missing

    dailyCols = ['STN',
             'TS',
             'JulianDate',
             'ETo_average_inches',
            ]
    dailyCols = [x.replace(' ','') for x in dailyCols]
    dailyCols = [x.replace('∞','') for x in dailyCols]
    dailyCols = [x.replace('=','') for x in dailyCols]

    dailyHistoricalWeatherPD = pd.DataFrame(columns=dailyCols)
    data = pd.DataFrame(columns = dailyCols)
    dtFiles = []

    startYr = startDtStr.split('-')[0]
    endYr = endDtStr.split('-')[0]
    years = list(set([startYr,endYr]))
    years.sort()


    for yr in years:
        annualStnFiles = [f for f in os.listdir(aus_obsPath_new) if (f.startswith('AUS') & f.endswith('.geojson'))]
        for f in annualStnFiles:
            try:
                y = f.split('__')[1].split('.')[0].startswith(yr)
                dtFiles.append(f)
            except:
                continue
    dtFiles.sort()

    for fl in dtFiles:
        f = aus_obsPath_new + fl
        data = gpd.read_file(f)
        #################################################### 2021-01-12
        data['aifstime_utc'] =  pd.to_datetime(data['aifstime_utc'],format='%Y%m%d%H%M%S')
        data['local_date_time_full'] =  pd.to_datetime(data['local_date_time_full'])
        data['Date'] = data.local_date_time_full

        ####################################################
        dailyHistoricalWeatherPD = pd.concat([dailyHistoricalWeatherPD, data], ignore_index=True, sort=True)


    dailyHistoricalWeatherPD['Date']=pd.to_datetime(dailyHistoricalWeatherPD['Date'], format='%Y%m%d%H%M%S', errors='ignore')


    dailyHistoricalWeatherPD = dailyHistoricalWeatherPD.replace('-', np.nan)
    dailyHistoricalWeatherPD.reset_index(drop=True, inplace=True)
    float_cols = ['PRCP','TAVG','TMAX','TMIN','ETo_FAO_IN', 'ETo_HAR_IN', 'ETo_AVG_IN']
    dailyHistoricalWeatherPD[float_cols] = dailyHistoricalWeatherPD[float_cols].astype('float')


    #resample to get Daily data for each station
    weather_df = pd.DataFrame()
    weather_all_df = pd.DataFrame()
    list_of_unique_stns = list(set(dailyHistoricalWeatherPD['StnID']))


    i=0
    for stn in list_of_unique_stns:
        df = dailyHistoricalWeatherPD[dailyHistoricalWeatherPD.StnID==stn]
        df.set_index('Date', drop=True, inplace=True)

        weather_df['StnID']=stn
        weather_df['T_mean__degC'] = df['TAVG'].resample('D').mean().squeeze()
        weather_df['T_max__degC'] = df['TMAX'].resample('D').max().squeeze()
        weather_df['T_min__degC'] = df['TMIN'].resample('D').min().squeeze()
        weather_df['Precip_in'] = df['PRCP'].resample('D').sum().squeeze()*0.0393701 #convert to inches

        weather_df['ETo_FAO_avg__in'] = df['ETo_FAO_IN'].resample('D').mean().squeeze()
        weather_df['ETo_HAR_avg__in'] = df['ETo_HAR_IN'].resample('D').mean().squeeze()
        weather_df['ETo_average_inches'] = df['ETo_AVG_IN'].resample('D').mean().squeeze()
        #weather_df['Date']=datetime.strftime(df.index[0],'%Y-%m-%d')
        weather_df['Lat']=df['lat'].resample('D').mean().squeeze()
        weather_df['Lon']=df['lon'].resample('D').mean().squeeze()


        weather_all_df =  weather_all_df.append(weather_df)


    #fix any missing dates
    ######################
    #Set the horizon from the inputs
    startDt = datetime.strptime(startDtStr, '%Y-%m-%d')
    endDt = datetime.strptime(endDtStr,'%Y-%m-%d')

    weather_all_df = getEToForMissingDates2(weather_all_df, startDt, endDt)
    #print(weather_all_df.columns)
    #print(weather_all_df.head())

    geom = [Point(xy) for xy in zip(weather_all_df['Lon'], weather_all_df['Lat'])]
    # Coordinate reference system : WGS84
    crs = {'init': 'epsg:4326'}
    # Creating a Geographic data frame 
    weather_gdf = gpd.GeoDataFrame(weather_all_df, crs=crs, geometry=geom)
    weather_gdf = weather_gdf.dropna(subset=['StnID'])

    #rename the stns
    weather_gdf['StationId']=weather_gdf['StnID']
    

    list_of_dates = [datetime.strptime(d.strftime('%Y-%m-%d'),'%Y-%m-%d') for d in pd.date_range(startDt, endDt)]

    startt = time.time()
    EToForAccounts = pd.DataFrame()
    for i,row in acct_gdf.iterrows():
        print('\t\tCalculating ETo for'+
            ' Customerid: '+str(row.customerid)+
            ', Assetid: '+str(row.fieldassetid)+
            ', Acct: ' + row.Category+
            ', Ranch: '+row.SubCategory+
            ', Field: '+row.Name)
        start_time = time.time()
        p = loads(row.CentroidPoint) #shapely Point obj

        #print(weather_gdf)
        etv = getWeather6ForDates(p,list_of_dates,weather_gdf)
        #print(etv)

        #now add the columns for the other fields to etv
        etv['UUID']=row.UUID
        etv['Category']=row.Category
        etv['SubCategory']=row.SubCategory
        etv['Name']=row.Name

        EToForAccounts = EToForAccounts.append(etv, ignore_index = True)
        print("\t\tEND ETo--- %s seconds ---" % (time.time() - start_time))
    print("\tETo Done --- %s seconds ---" % (time.time() - startt))
    print('\tDone!')

    EToForAccounts = EToForAccounts.rename(columns={"ETo_average_inches": "ETo(in)"})


    #Merge the acct
    merged_df = pd.merge(acct_gdf, EToForAccounts, on=['Category','SubCategory','Name', 'UUID'])
    merged_df.reset_index() #reset the index to avoid any export issues

    return merged_df

def getEToForMissingDates(weather_all_df, startDt, endDt):
    #only look at the dates that are in the range
    weather_all_df = weather_all_df[((weather_all_df.index>=startDt) & (weather_all_df.index<=endDt))]
    #drop the rows that have Lat or Lon as NaN
    weather_all_df = weather_all_df.dropna(subset = ['Lon', 'Lat'])

    list_of_stn_LonLat_points = [(row.Lon, row.Lat) for i,row in weather_all_df.iterrows()]
    list_of_stn_LonLat_points = list(set(list_of_stn_LonLat_points))
    weather_all_df_fixed = pd.DataFrame()
    for pt in list_of_stn_LonLat_points:
        #print(pt)
        lon=pt[0]
        df = weather_all_df[((weather_all_df.Lat==pt[1]) & (weather_all_df.Lon==pt[0]))]
        #find the missing dates
        currDtList = [datetime.strftime(x,'%Y-%m-%d') for x in df.index.tolist()]
        reqdDtList = [datetime.strftime(x,'%Y-%m-%d') for x in pd.date_range(startDt+timedelta(days=1),endDt,freq='d')]
        newDtList = [x for x in reqdDtList if x not in currDtList]
        #print(newDtList)

        if len(newDtList)>0:
            #add these datapoints to the df
            new_df=df.copy()
            new_df.reset_index(inplace=True)
            stnp = Point(pt)
            idx=len(new_df)
            stnid=new_df.StnID[0]
            for dtStr in newDtList:
                eto = getForecastedAUSEToForPtAndDtStr(dtStr,stnp)
                #print(eto)        
                new_df.loc[idx,'Date']=datetime.strptime(dtStr,'%Y-%m-%d')
                new_df.loc[idx,'StnID']=stnid
                new_df.loc[idx,'ETo_average_inches']=eto
                new_df.loc[idx,'Lon']=pt[0]
                new_df.loc[idx,'Lat']=pt[1]
                idx+=1
            new_df.Date = pd.to_datetime(new_df.Date)
            new_df['aifstime_utc']=new_df.Date
            new_df['local_date_time_full']=getLocalTime(new_df.Date.tolist(), lon)
            new_df= new_df.sort_values(by=['Date'])
            #print(new_df)

            weather_all_df_fixed = weather_all_df_fixed.append(new_df, ignore_index=True)
        else:
            weather_all_df_fixed = weather_all_df_fixed.append(df, ignore_index=True)
        
        #print(weather_all_df_fixed)	

    weather_all_df_fixed.reset_index(drop=True, inplace=True)
    
    return weather_all_df_fixed


def getEToForMissingDates2(weather_all_df, startDt, endDt):
    #drop the rows that have Lat or Lon as NaN
    weather_all_df = weather_all_df.dropna(subset = ['Lon', 'Lat'])

    list_of_stn_LonLat_points = [(row.Lon, row.Lat) for i,row in weather_all_df.iterrows()]
    list_of_stn_LonLat_points = list(set(list_of_stn_LonLat_points))

    #only look at the dates that are in the range
    weather_all_df = weather_all_df[((weather_all_df.index>=startDt) & (weather_all_df.index<=endDt))]

    cols = weather_all_df.columns
    
    weather_all_df_fixed = pd.DataFrame()
    for pt in list_of_stn_LonLat_points:
        #print(pt)
        lon=pt[0]
        lat=pt[1]

        df = weather_all_df[((weather_all_df.Lat==pt[1]) & (weather_all_df.Lon==pt[0]))]
        if len(df)==0: #all dates are missing
            df = pd.DataFrame(columns=cols)
            newDtList = [datetime.strftime(x,'%Y-%m-%d') for x in pd.date_range(startDt+timedelta(days=1),endDt,freq='d')]
            stnid=str(s2.CellId.from_lat_lng(s2.LatLng.from_degrees(float(lat), float(lon))).to_token())
        else:
            #find the missing dates
            currDtList = [datetime.strftime(x,'%Y-%m-%d') for x in df.index.tolist()]
            reqdDtList = [datetime.strftime(x,'%Y-%m-%d') for x in pd.date_range(startDt+timedelta(days=1),endDt,freq='d')]
            newDtList = [x for x in reqdDtList if x not in currDtList]
            stnid=df.StnID[0]
        

        if len(newDtList)>0:
            #add these datapoints to the df
            new_df=df.copy()
            new_df.reset_index(inplace=True)
            stnp = Point(pt)
            idx=len(new_df)
            
            for dtStr in newDtList:
                #First get this from the Daily Summary
                eto = getDailySummaryAUSEToForPtAndDtStr(dtStr,stnp)
                if np.isnan(eto):
                    #If eto is 0 or nan for that dtStr, then get it from forecast
                    eto = getForecastedAUSEToForPtAndDtStr(dtStr,stnp)
                    
                    #print(eto)        
                new_df.loc[idx,'Date']=datetime.strptime(dtStr,'%Y-%m-%d')
                new_df.loc[idx,'StnID']=stnid
                new_df.loc[idx,'ETo_average_inches']=eto
                new_df.loc[idx,'Lon']=pt[0]
                new_df.loc[idx,'Lat']=pt[1]
                idx+=1
            new_df.Date = pd.to_datetime(new_df.Date)
            new_df['aifstime_utc']=new_df.Date
            #print(lon)
            new_df['local_date_time_full']=getLocalTime(new_df.Date.tolist(), lon)
            new_df= new_df.sort_values(by=['Date'])
            #print(new_df)

            weather_all_df_fixed = weather_all_df_fixed.append(new_df, ignore_index=True)
        else:
            weather_all_df_fixed = weather_all_df_fixed.append(df, ignore_index=True)
        
        #print(weather_all_df_fixed)	

    weather_all_df_fixed.reset_index(drop=True, inplace=True)
    
    return weather_all_df_fixed


def multiprocessing_etc_aus(merged_df,dt):
    time.sleep(2)
    dt_new = dt
    #dt.astype('M8[D]').astype('O')
    generate_MP_ETc_AUS(merged_df,dt_new)
    return

def generate_MP_ETc_AUS(merged_df,dt):
    unique_accounts = pd.unique(merged_df.customerid)
    startt = time.time()
    #Now to put this into the output folder for each account
    print("Start --- %s  ---" % datetime.now())

    for cid in unique_accounts:
        a = 'JainlogicCustomerid_'+str(cid)
        #try:
        starta = time.time()
        
        print('\n\tWriting...%s ' % a)
        print("\tStart --- %s  ---" % datetime.now())

        acctName = a.replace(' ','').upper()

        acct_dir = acctOutputPath + a + '/data/'
        acct_dir_img = acctOutputPath + a + '/img/'
        acct_dir_pkl = acctOutputPath + a + '/pkl/'
        acct_dir_rep = acctOutputPath + a + '/reports/'
            
        #create the directory if one doesn't exist
        if not os.path.exists(acct_dir):
            os.makedirs(acct_dir)
        if not os.path.exists(acct_dir_img):
            os.makedirs(acct_dir_img)
        if not os.path.exists(acct_dir_pkl):
            os.makedirs(acct_dir_pkl)
        if not os.path.exists(acct_dir_rep):
            os.makedirs(acct_dir_rep)
            
        #For each account, take the subset of the acct_df 
        #acctDF = merged_df[merged_df.Category == a]
        acctDF = merged_df[merged_df.customerid == cid]
        
        acctGDF = gpd.GeoDataFrame(acctDF,geometry=acctDF.geometry)
        
        #for every day, do the following analysis and printout
        acct_dates_list = acctGDF.Date.tolist()
        unique_dts = pd.unique(acct_dates_list) # so we don't multiple count for each field
        unique_dts.sort()
        
        #for dt in unique_dts:
            
        endDay = dt        
        
        startd = time.time()
        #print("\n\t\tSTART --- %s  ---" % datetime.now())

        etoDateStr =  datetime.strftime(dt,'%Y-%m-%d')
        jsonFile = acct_dir + acctName+'__ET_'+etoDateStr+'.json'
        
        acctGDF_date = acctGDF[acctGDF.Date==dt]
        acctGDF_date.reset_index(drop=True, inplace=True)
        
        print('\t\t'+etoDateStr)
        for i,row in acctGDF_date.iterrows():
            
            customerid = str(row.customerid)
            acct = row.Category
            ranch = row.SubCategory
            field = row.Name
            fieldassetid = 'ASSETID_'+str(row.fieldassetid)
            
            #ETo and weather values
            eto = row['ETo(in)']
            tavg = row['Tavg']
            #rnet = row['Rnet']
            wavg = row['Wavg']
            ############## add others Raw and Deried
            precip = row['Precip_in']
            maxAirTemp = row['T_max']
            minAirTemp = row['T_min']
            avgAirTemp = row['T_mean']
            maxRelHumPercent = row['RH_max']
            minRelHumPercent = row['RH_min']
            avgRelHumPercent = row['RH_mean']
            dewPoint = row['T_dew']
            
            #########################################
        
            # NDVI value
            P = getSinglePoly(row.geometry)

            
            ss = row.SentinelFootprint
            ll = row.LandsatFootprint
            
            #do Sentinel
            if not type(ss)==list:
                if (ss.find(',')==-1): #There is only ONE tile
                    S = [ss.replace('[','').replace(']','').replace('\'','')]
                else:
                    tok = ss.split(',')
                    l = len(tok)
                    S=[]
                    for ii in range(l):
                        tk = tok[ii].replace('[','').replace(']','').replace('\'','')
                        S.append(tk)
            else:
                S = ss
            if (type(S)==str):
                S = [S]
        
            #do landsat
            if not type(ll)==list:
                if (ss.find(',')==-1): #There is only ONE tile
                    L = ll.replace('[','').replace(']','').replace('\'','')
                else:
                    tok = ll.split(', ')
                    l = len(tok)
                    L=[]
                    for ii in range(l):
                        tk = tok[ii].replace('[','').replace(']','').replace('\'','')
                        L.append(tk.replace(' ',''))
            else:
                L = [str(x).zfill(6) for x in ll]

            if (type(L)==str):
                L = [L]
            
            ndvi_mList_df = pd.DataFrame()
            
            NDVIFiles_df = getNDVIFilesWithDate(32, S, L, etoDateStr)
            NDVIFiles_df.reset_index(drop=True,inplace=True)
            #print(NDVIFiles_df)
            
            [ndvi_df,ndvi_img,ndvi_ext,imgDtStr] = getNDVIValueWithStatsAndImg3(NDVIFiles_df,P)
            
            
            if (len(ndvi_df)>0 & (len(ndvi_img)>0)): #The first time the df has valid values
                ext = json.loads(ndvi_ext)
                ndvi_mList_df = ndvi_df.copy()
                
                
                acct_dir_img2 = acct_dir_img + fieldassetid + '/'
                if not os.path.exists(acct_dir_img2):
                    os.makedirs(acct_dir_img2)
                imgFile = acct_dir_img2 + acctName+'_'+fieldassetid+'__IMG_'+etoDateStr+'.png'
                infoJsonFile = acct_dir_img2 + acctName+'_'+fieldassetid+'__INFO_'+etoDateStr+'.json'

                acct_dir_pkl2 = acct_dir_pkl + fieldassetid + '/'
                if not os.path.exists(acct_dir_pkl2):
                    os.makedirs(acct_dir_pkl2)
                pklFile = acct_dir_pkl2 + acctName+'_'+fieldassetid+'__PKL_'+etoDateStr+'.pkl'
                
            
                imgDict = {
                    'ImgFileName': os.path.split(imgFile)[1],
                    'ImgDate': imgDtStr,
                    'ImgType': 'ETc',
                    'ImgExtent': ndvi_ext, 
                    'ImgBoundary': P.to_wkt(),
                    'AssetID': row.fieldassetid,
                    'CustomerID': row.customerid,
                    'Ranch': ranch,
                    'Field': field
                }
                #'ImgType': 'NDVI',
                with open(infoJsonFile, 'w', encoding='utf8') as outfile:
                    json.dump(imgDict, outfile, ensure_ascii=False)

                plt.ioff();
                with open(imgFile, 'wb') as outfile:
                    fig = plt.figure()
                    plt.imshow(ndvi_img, extent=ext, cmap=zoom1, vmin=0, vmax=1)
                    fig.gca().set_axis_off()
                    plt.savefig(outfile, dpi=600, bbox_inches='tight', pad_inches=0, transparent=True)
                    plt.close();
                    
                with open(pklFile,'wb') as f:
                    pickle.dump(ndvi_img, f)
            
                retdf = ndvi_df.loc[0]

                ndvi_mean = retdf.Mean
                ndvi_p90 = retdf.p90
                ndvi_p25 = retdf.p25
                ndvi_p75 = retdf.p75
                ndvi_median = retdf.Median
                ndvi_std = retdf.StdDev
                
            else:
                #set values if NaN
                ndvi_mean = np.nan
                ndvi_p90 = np.nan
                ndvi_p25 = np.nan
                ndvi_p75 = np.nan
                ndvi_median = np.nan
                ndvi_std = np.nan
                
                
            #derived
            ndvi_unif = (ndvi_p25 / ndvi_median)*100
        
            #compute the METRIC ETc
            et_metric = 0.0 #Don't compute METRIC ET due to no rnet
            #computeMetricET(tavg,rnet,wavg,ndvi_mean)
            
            ##########if et_metric is NaN, make it zero
            if (np.isnan(et_metric)):
                et_metric=0.0
        
            #the Kc
            fac = 1.35 #somewhat arbitrary for now
            k_c = 1.3558*ndvi_mean*fac + 0.05
            k_c2 = 1.3558*ndvi_mean + 0.05
        
            #compute the ETc
            et_ndvi = eto*k_c
            et_ndvi2 = eto*k_c2
        
            #get the final ETc
            # Include et_metric
            k_a = 1 #We can reduce k_a to increase the effect of Metric (see above)
            et_c =  (1.15 - k_a)*et_metric + ((k_a + 0.0)*et_ndvi)
            et_c2 =  (1.15 - k_a)*et_metric + ((k_a + 0.0)*et_ndvi2)
            
            k_c_p90 = ((1.3558*ndvi_p90*fac + 0.05) + (1.3558*ndvi_p90 + 0.05))/2
            et_p90 = eto*k_c_p90
            
            acctGDF_date.loc[i,'ETc']=et_c
            acctGDF_date.loc[i,'ETc_low']=et_c2
            acctGDF_date.loc[i,'Kc']=k_c
            acctGDF_date.loc[i,'ETm']=et_metric
            acctGDF_date.loc[i,'ET_p90']=et_p90
            acctGDF_date.loc[i,'ET_unif']=ndvi_unif
            acctGDF_date.loc[i,'Precipitation(in)']=precip
            acctGDF_date.loc[i,'MaximumAirTemperature(F)']=maxAirTemp
            acctGDF_date.loc[i,'MinimumAirTemperature(F)']=minAirTemp
            acctGDF_date.loc[i,'AverageAirTemperature(F)']=avgAirTemp
            acctGDF_date.loc[i,'MaximumRelativeHumidity(%)']=maxRelHumPercent
            acctGDF_date.loc[i,'MinimumRelativeHumidity(%)']=minRelHumPercent
            acctGDF_date.loc[i,'AverageRelativeHumidity(%)']=avgRelHumPercent
            acctGDF_date.loc[i,'DewPoint(F)']=dewPoint
            

        acctGDF_date['SentinelFootprint'] = acctGDF_date['SentinelFootprint'].astype('str')
        acctGDF_date['LandsatFootprint'] = acctGDF_date['LandsatFootprint'].astype('str')
        acctGDF_date['Date'] = acctGDF_date['Date'].astype('str')
        if 'CentroidPoint' in acctGDF_date.columns:
            acctGDF_date.drop('CentroidPoint',inplace=True, axis=1)
        
        
        #save the values
        acctGDF_date.to_file(jsonFile, driver="GeoJSON")

        #Print to screen
        print(acctGDF_date.ETc.tolist())
        print("\t\tElapsed --- %s seconds ---" % (time.time() - startd))

        print("\tEND --- %s seconds ---" % (time.time() - starta))
        print('\tDone!')
        #except:
            #print('Exception in Acct: '+str(Exception))
            #continue

    print("\nEND --- %s seconds ---" % (time.time() - startt))
    return


def generateETc_NOAA3(merged_df):

    unique_accounts = pd.unique(merged_df.customerid)
    startt = time.time()
    #Now to put this into the output folder for each account
    print("Start --- %s  ---" % datetime.now())

    for cid in unique_accounts:
        a = 'JainlogicCustomerid_'+str(cid)
        #try:
        starta = time.time()
        
        print('\n\tWriting...%s ' % a)
        print("\tStart --- %s  ---" % datetime.now())

        acctName = a.replace(' ','').upper()

        acct_dir = acctOutputPath + a + '/data/'
        acct_dir_img = acctOutputPath + a + '/img/'
        acct_dir_pkl = acctOutputPath + a + '/pkl/'
        acct_dir_rep = acctOutputPath + a + '/reports/'
            
        #create the directory if one doesn't exist
        if not os.path.exists(acct_dir):
            os.makedirs(acct_dir)
        if not os.path.exists(acct_dir_img):
            os.makedirs(acct_dir_img)
        if not os.path.exists(acct_dir_pkl):
            os.makedirs(acct_dir_pkl)
        if not os.path.exists(acct_dir_rep):
            os.makedirs(acct_dir_rep)
            
        #For each account, take the subset of the acct_df 
        #acctDF = merged_df[merged_df.Category == a]
        acctDF = merged_df[merged_df.customerid == cid]
        
        acctGDF = gpd.GeoDataFrame(acctDF,geometry=acctDF.geometry)
        
        #for every day, do the following analysis and printout
        acct_dates_list = acctGDF.Date.tolist()
        unique_dts = pd.unique(acct_dates_list) # so we don't multiple count for each field
        unique_dts.sort()
        
        for dt in unique_dts:
            
            endDay = dt        
            
            startd = time.time()
            #print("\n\t\tSTART --- %s  ---" % datetime.now())

            etoDateStr =  datetime.strftime(dt,'%Y-%m-%d')
            jsonFile = acct_dir + acctName+'__ET_'+etoDateStr+'.json'
            
            acctGDF_date = acctGDF[acctGDF.Date==dt]
            acctGDF_date.reset_index(drop=True, inplace=True)
            
            print('\t\t'+etoDateStr)
            for i,row in acctGDF_date.iterrows():
                
                customerid = str(row.customerid)
                acct = row.Category
                ranch = row.SubCategory
                field = row.Name
                fieldassetid = 'ASSETID_'+str(row.fieldassetid)
                
                #ETo and weather values
                eto = row['ETo(in)']
                """
                tavg = row['Tavg']
                #rnet = row['Rnet']
                wavg = row['Wavg']
                ############## add others Raw and Deried
                precip = row['Precip_in']
                maxAirTemp = row['T_max']
                minAirTemp = row['T_min']
                avgAirTemp = row['T_mean']
                maxRelHumPercent = row['RH_max']
                minRelHumPercent = row['RH_min']
                avgRelHumPercent = row['RH_mean']
                dewPoint = row['T_dew']
                """
                
                #########################################
            
                # NDVI value
                P = getSinglePoly(row.geometry)

                
                ss = row.SentinelFootprint
                ll = row.LandsatFootprint
                
                #do Sentinel
                if not type(ss)==list:
                    if (ss.find(',')==-1): #There is only ONE tile
                        S = [ss.replace('[','').replace(']','').replace('\'','')]
                    else:
                        tok = ss.split(',')
                        l = len(tok)
                        S=[]
                        for ii in range(l):
                            tk = tok[ii].replace('[','').replace(']','').replace('\'','')
                            S.append(tk)
                else:
                    S = ss
                if (type(S)==str):
                    S = [S]
            
                #do landsat
                if not type(ll)==list:
                    if (ss.find(',')==-1): #There is only ONE tile
                        L = ll.replace('[','').replace(']','').replace('\'','')
                    else:
                        tok = ll.split(', ')
                        l = len(tok)
                        L=[]
                        for ii in range(l):
                            tk = tok[ii].replace('[','').replace(']','').replace('\'','')
                            L.append(tk.replace(' ',''))
                else:
                    L = [str(x).zfill(6) for x in ll]

                if (type(L)==str):
                    L = [L]
                
                ndvi_mList_df = pd.DataFrame()
                
                NDVIFiles_df = getNDVIFilesWithDate(32, S, L, etoDateStr)
                NDVIFiles_df.reset_index(drop=True,inplace=True)
                #print(S, L, NDVIFiles_df)
                
                [ndvi_df,ndvi_img,ndvi_ext,imgDtStr] = getNDVIValueWithStatsAndImg3(NDVIFiles_df,P)
                
                
                if (len(ndvi_df)>0 & (len(ndvi_img)>0)): #The first time the df has valid values
                    ext = json.loads(ndvi_ext)
                    ndvi_mList_df = ndvi_df.copy()
                    
                    
                    acct_dir_img2 = acct_dir_img + fieldassetid + '/'
                    if not os.path.exists(acct_dir_img2):
                        os.makedirs(acct_dir_img2)
                    imgFile = acct_dir_img2 + acctName+'_'+fieldassetid+'__IMG_'+etoDateStr+'.png'
                    infoJsonFile = acct_dir_img2 + acctName+'_'+fieldassetid+'__INFO_'+etoDateStr+'.json'

                    acct_dir_pkl2 = acct_dir_pkl + fieldassetid + '/'
                    if not os.path.exists(acct_dir_pkl2):
                        os.makedirs(acct_dir_pkl2)
                    pklFile = acct_dir_pkl2 + acctName+'_'+fieldassetid+'__PKL_'+etoDateStr+'.pkl'
                    
                
                    imgDict = {
                        'ImgFileName': os.path.split(imgFile)[1],
                        'ImgDate': imgDtStr,
                        'ImgType': 'ETc',
                        'ImgExtent': ndvi_ext, 
                        'ImgBoundary': P.to_wkt(),
                        'AssetID': row.fieldassetid,
                        'CustomerID': row.customerid,
                        'Ranch': ranch,
                        'Field': field
                    }
                    with open(infoJsonFile, 'w', encoding='utf8') as outfile:
                        json.dump(imgDict, outfile, ensure_ascii=False)

                    plt.ioff();
                    with open(imgFile, 'wb') as outfile:
                        fig = plt.figure()
                        plt.imshow(ndvi_img, extent=ext, cmap=zoom1, vmin=0, vmax=1)
                        fig.gca().set_axis_off()
                        plt.savefig(outfile, dpi=600, bbox_inches='tight', pad_inches=0, transparent=True)
                        plt.close();
                        
                    with open(pklFile,'wb') as f:
                        pickle.dump(ndvi_img, f)
                
                    retdf = ndvi_df.loc[0]

                    ndvi_mean = retdf.Mean
                    ndvi_p90 = retdf.p90
                    ndvi_p25 = retdf.p25
                    ndvi_p75 = retdf.p75
                    ndvi_median = retdf.Median
                    ndvi_std = retdf.StdDev
                    
                else:
                    #set values if NaN
                    ndvi_mean = np.nan
                    ndvi_p90 = np.nan
                    ndvi_p25 = np.nan
                    ndvi_p75 = np.nan
                    ndvi_median = np.nan
                    ndvi_std = np.nan
                    
                    
                #derived
                ndvi_unif = (ndvi_p25 / ndvi_median)*100
            
                #compute the METRIC ETc
                et_metric = 0.0 #Don't compute METRIC ET due to no rnet
                #computeMetricET(tavg,rnet,wavg,ndvi_mean)
                
                ##########if et_metric is NaN, make it zero
                if (np.isnan(et_metric)):
                    et_metric=0.0
            
                #the Kc
                fac = 1.35 #somewhat arbitrary for now
                k_c = 1.3558*ndvi_mean*fac + 0.05
                k_c2 = 1.3558*ndvi_mean + 0.05
            
                #compute the ETc
                et_ndvi = eto*k_c
                et_ndvi2 = eto*k_c2
            
                #get the final ETc
                # Include et_metric
                k_a = 1 #We can reduce k_a to increase the effect of Metric (see above)
                et_c =  (1.15 - k_a)*et_metric + ((k_a + 0.0)*et_ndvi)
                et_c2 =  (1.15 - k_a)*et_metric + ((k_a + 0.0)*et_ndvi2)
                
                k_c_p90 = ((1.3558*ndvi_p90*fac + 0.05) + (1.3558*ndvi_p90 + 0.05))/2
                et_p90 = eto*k_c_p90
                
                acctGDF_date.loc[i,'ETc']=et_c
                acctGDF_date.loc[i,'ETc_low']=et_c2
                acctGDF_date.loc[i,'Kc']=k_c
                acctGDF_date.loc[i,'ETm']=et_metric
                acctGDF_date.loc[i,'ET_p90']=et_p90
                acctGDF_date.loc[i,'ET_unif']=ndvi_unif
                """
                acctGDF_date.loc[i,'Precipitation(in)']=precip
                acctGDF_date.loc[i,'MaximumAirTemperature(F)']=maxAirTemp
                acctGDF_date.loc[i,'MinimumAirTemperature(F)']=minAirTemp
                acctGDF_date.loc[i,'AverageAirTemperature(F)']=avgAirTemp
                acctGDF_date.loc[i,'MaximumRelativeHumidity(%)']=maxRelHumPercent
                acctGDF_date.loc[i,'MinimumRelativeHumidity(%)']=minRelHumPercent
                acctGDF_date.loc[i,'AverageRelativeHumidity(%)']=avgRelHumPercent
                acctGDF_date.loc[i,'DewPoint(F)']=dewPoint
                """
                

            acctGDF_date['SentinelFootprint'] = acctGDF_date['SentinelFootprint'].astype('str')
            acctGDF_date['LandsatFootprint'] = acctGDF_date['LandsatFootprint'].astype('str')
            acctGDF_date['Date'] = acctGDF_date['Date'].astype('str')
            if 'CentroidPoint' in acctGDF_date.columns:
                acctGDF_date.drop('CentroidPoint',inplace=True, axis=1)
            
            
            #save the values
            acctGDF_date.to_file(jsonFile, driver="GeoJSON")

            #Print to screen
            print(acctGDF_date.ETc.tolist())
            print("\tElapsed --- %s seconds ---" % (time.time() - startd))

        print("\tEND --- %s seconds ---" % (time.time() - starta))
        print('\tDone!')
        #except:
        #	print('Exception in Acct: '+str(Exception))
        #	continue

    print("\nEND --- %s seconds ---" % (time.time() - startt))
    return

def generate_ISR_ETc(merged_df):

    unique_accounts = pd.unique(merged_df.Category)
    startt = time.time()
    #Now to put this into the output folder for each account
    print("Start --- %s  ---" % datetime.now())

    for cid in unique_accounts:
        a = str(cid)
        #try:
        starta = time.time()
        
        print('\n\tWriting...%s ' % a)
        print("\tStart --- %s  ---" % datetime.now())

        acctName = a.replace(' ','').upper()

        acct_dir = pilotAcctOutputPath + a + '/data/'
        acct_dir_img = pilotAcctOutputPath + a + '/img/'
        acct_dir_pkl = pilotAcctOutputPath + a + '/pkl/'
        acct_dir_rep = pilotAcctOutputPath + a + '/reports/'
            
        #create the directory if one doesn't exist
        if not os.path.exists(acct_dir):
            os.makedirs(acct_dir)
        if not os.path.exists(acct_dir_img):
            os.makedirs(acct_dir_img)
        if not os.path.exists(acct_dir_pkl):
            os.makedirs(acct_dir_pkl)
        if not os.path.exists(acct_dir_rep):
            os.makedirs(acct_dir_rep)
            
        #For each account, take the subset of the acct_df 
        #acctDF = merged_df[merged_df.Category == a]
        acctDF = merged_df[merged_df.Category == cid]
        
        acctGDF = gpd.GeoDataFrame(acctDF,geometry=acctDF.geometry)
        
        #for every day, do the following analysis and printout
        acct_dates_list = acctGDF.Date.tolist()
        unique_dts = pd.unique(acct_dates_list) # so we don't multiple count for each field
        unique_dts.sort()
        
        for dt in unique_dts:
            
            endDay = dt        
            
            startd = time.time()
            #print("\n\t\tSTART --- %s  ---" % datetime.now())

            etoDateStr =  datetime.strftime(dt,'%Y-%m-%d')
            jsonFile = acct_dir + acctName+'__ET_'+etoDateStr+'.json'
            
            acctGDF_date = acctGDF[acctGDF.Date==dt]
            acctGDF_date.reset_index(drop=True, inplace=True)
            
            print('\t\t'+etoDateStr)
            for i,row in acctGDF_date.iterrows():
                
                acct = row.Category
                ranch = row.SubCategory
                field = row.Name.replace(' ','')
                
                #ETo and weather values
                eto = row['ETo(in)']
                
                            
                # NDVI value
                P = getSinglePoly(row.geometry)

                ss = row.SentinelFootprint
                ll = row.LandsatFootprint
                
                #do Sentinel
                if not type(ss)==list:
                    if (ss.find(',')==-1): #There is only ONE tile
                        S = [ss.replace('[','').replace(']','').replace('\'','')]
                    else:
                        tok = ss.split(',')
                        l = len(tok)
                        S=[]
                        for ii in range(l):
                            tk = tok[ii].replace('[','').replace(']','').replace('\'','')
                            S.append(tk)
                else:
                    S = ss
                if (type(S)==str):
                    S = [S]
            
                #do landsat
                if not type(ll)==list:
                    if (ss.find(',')==-1): #There is only ONE tile
                        L = ll.replace('[','').replace(']','').replace('\'','')
                    else:
                        tok = ll.split(', ')
                        l = len(tok)
                        L=[]
                        for ii in range(l):
                            tk = tok[ii].replace('[','').replace(']','').replace('\'','')
                            L.append(tk.replace(' ',''))
                else:
                    L = [str(x).zfill(6) for x in ll]

                if (type(L)==str):
                    L = [L]
                
                ndvi_mList_df = pd.DataFrame()
                
                NDVIFiles_df = getNDVIFilesWithDate(32, S, L, etoDateStr)
                NDVIFiles_df.reset_index(drop=True,inplace=True)
                #print(S, L, NDVIFiles_df)
                
                [ndvi_df,ndvi_img,ndvi_ext,imgDtStr] = getNDVIValueWithStatsAndImg4(NDVIFiles_df,P)
                
                
                if (len(ndvi_df)>0 & (len(ndvi_img)>0)): #The first time the df has valid values
                    ext = json.loads(ndvi_ext)
                    ndvi_mList_df = ndvi_df.copy()
                    
                    
                    acct_dir_img2 = acct_dir_img 
                    if not os.path.exists(acct_dir_img2):
                        os.makedirs(acct_dir_img2)
                    imgFile = acct_dir_img2 + acctName+'_'+field+'__IMG_'+etoDateStr+'.png'
                    infoJsonFile = acct_dir_img2 + acctName+'_'+field+'__INFO_'+etoDateStr+'.json'

                    acct_dir_pkl2 = acct_dir_pkl 
                    if not os.path.exists(acct_dir_pkl2):
                        os.makedirs(acct_dir_pkl2)
                    pklFile = acct_dir_pkl2 + acctName+'_'+field+'__PKL_'+etoDateStr+'.pkl'
                    
                
                    imgDict = {
                        'ImgFileName': os.path.split(imgFile)[1],
                        'ImgDate': imgDtStr,
                        'ImgType': 'ETc',
                        'ImgExtent': ndvi_ext, 
                        'ImgBoundary': P.to_wkt(),
                        'SubCategory': row.SubCategory,
                        'Category': row.Category,
                        'Ranch': ranch,
                        'Field': field
                    }
                    #print(infoJsonFile)
                    with open(infoJsonFile, 'w', encoding='utf8') as outfile:
                        json.dump(imgDict, outfile, ensure_ascii=False)

                    plt.ioff();
                    #print(imgFile)
                    with open(imgFile, 'wb') as outfile:
                        fig = plt.figure()
                        plt.imshow(ndvi_img, extent=ext, cmap=zoom1, vmin=0, vmax=1)
                        fig.gca().set_axis_off()
                        plt.savefig(outfile, dpi=600, bbox_inches='tight', pad_inches=0, transparent=True)
                        plt.close();
                        
                    #print(pklFile)
                    with open(pklFile,'wb') as f:
                        pickle.dump(ndvi_img, f)
                
                    retdf = ndvi_df.loc[0]

                    ndvi_mean = retdf.Mean
                    ndvi_p90 = retdf.p90
                    ndvi_p25 = retdf.p25
                    ndvi_p75 = retdf.p75
                    ndvi_median = retdf.Median
                    ndvi_std = retdf.StdDev
                    
                else:
                    #set values if NaN
                    ndvi_mean = np.nan
                    ndvi_p90 = np.nan
                    ndvi_p25 = np.nan
                    ndvi_p75 = np.nan
                    ndvi_median = np.nan
                    ndvi_std = np.nan
                    
                    
                #derived
                ndvi_unif = (ndvi_p25 / ndvi_median)*100
            
                #compute the METRIC ETc
                et_metric = 0.0 #Don't compute METRIC ET due to no rnet
                #computeMetricET(tavg,rnet,wavg,ndvi_mean)
                
                ##########if et_metric is NaN, make it zero
                if (np.isnan(et_metric)):
                    et_metric=0.0
            
                #the Kc
                fac = 1.35 #somewhat arbitrary for now
                k_c = 1.3558*ndvi_mean*fac + 0.05
                k_c2 = 1.3558*ndvi_mean + 0.05
            
                #compute the ETc
                et_ndvi = eto*k_c
                et_ndvi2 = eto*k_c2
            
                #get the final ETc
                # Include et_metric
                k_a = 1 #We can reduce k_a to increase the effect of Metric (see above)
                et_c =  (1.15 - k_a)*et_metric + ((k_a + 0.0)*et_ndvi)
                et_c2 =  (1.15 - k_a)*et_metric + ((k_a + 0.0)*et_ndvi2)
                
                k_c_p90 = ((1.3558*ndvi_p90*fac + 0.05) + (1.3558*ndvi_p90 + 0.05))/2
                et_p90 = eto*k_c_p90
                
                acctGDF_date.loc[i,'ETc']=et_c
                acctGDF_date.loc[i,'ETc_low']=et_c2
                acctGDF_date.loc[i,'Kc']=k_c
                acctGDF_date.loc[i,'ETm']=et_metric
                acctGDF_date.loc[i,'ET_p90']=et_p90
                acctGDF_date.loc[i,'ET_unif']=ndvi_unif			

            acctGDF_date['SentinelFootprint'] = acctGDF_date['SentinelFootprint'].astype('str')
            acctGDF_date['LandsatFootprint'] = acctGDF_date['LandsatFootprint'].astype('str')
            acctGDF_date['Date'] = acctGDF_date['Date'].astype('str')
            if 'CentroidPoint' in acctGDF_date.columns:
                acctGDF_date.drop('CentroidPoint',inplace=True, axis=1)
            
            
            #save the values
            #print(jsonFile)
            acctGDF_date.to_file(jsonFile, driver="GeoJSON")

            #Print to screen
            print(acctGDF_date.ETc.tolist())
            print("\tElapsed --- %s seconds ---" % (time.time() - startd))

        print("\tEND --- %s seconds ---" % (time.time() - starta))
        print('\tDone!')
        #except:
        #	print('Exception in Acct: '+str(Exception))
        #	continue

    print("\nEND --- %s seconds ---" % (time.time() - startt))
    return


def generate_ISR_ETc2(merged_df, measureType):

    unique_accounts = pd.unique(merged_df.Category)
    startt = time.time()
    #Now to put this into the output folder for each account
    print("Start --- %s  ---" % datetime.now())

    for cid in unique_accounts:
        a = str(cid)
        try:
            starta = time.time()
            
            print('\n\tWriting...%s ' % a)
            print("\tStart --- %s  ---" % datetime.now())

            acctName = a.replace(' ','').upper()

            acct_dir = pilotAcctOutputPath + a + '/data/'
            acct_dir_img = pilotAcctOutputPath + a + '/img/'
            acct_dir_pkl = pilotAcctOutputPath + a + '/pkl/'
            acct_dir_rep = pilotAcctOutputPath + a + '/reports/'
                
            #create the directory if one doesn't exist
            if not os.path.exists(acct_dir):
                os.makedirs(acct_dir)
            if not os.path.exists(acct_dir_img):
                os.makedirs(acct_dir_img)
            if not os.path.exists(acct_dir_pkl):
                os.makedirs(acct_dir_pkl)
            if not os.path.exists(acct_dir_rep):
                os.makedirs(acct_dir_rep)
                
            #For each account, take the subset of the acct_df 
            #acctDF = merged_df[merged_df.Category == a]
            acctDF = merged_df[merged_df.Category == cid]
            
            acctGDF = gpd.GeoDataFrame(acctDF,geometry=acctDF.geometry)
            
            #for every day, do the following analysis and printout
            acct_dates_list = acctGDF.Date.tolist()
            unique_dts = pd.unique(acct_dates_list) # so we don't multiple count for each field
            unique_dts.sort()
            
            for dt in unique_dts:
                try:
                    endDay = dt        
                    
                    startd = time.time()
                    #print("\n\t\tSTART --- %s  ---" % datetime.now())

                    etoDateStr =  datetime.strftime(dt,'%Y-%m-%d')
                    jsonFile = acct_dir + acctName+'__ET_'+etoDateStr+'.json'
                    
                    acctGDF_date = acctGDF[acctGDF.Date==dt]
                    acctGDF_date.reset_index(drop=True, inplace=True)
                    
                    print('\t\t'+etoDateStr)
                    for i,row in acctGDF_date.iterrows():
                        
                        acct = row.Category
                        ranch = row.SubCategory
                        field = row.Name.replace(' ','')
                        
                        #ETo and weather values
                        eto = row['ETo(in)']
                        
                                    
                        # NDVI value
                        P = getSinglePoly(row.geometry)

                        ss = row.SentinelFootprint
                        ll = row.LandsatFootprint
                        
                        #do Sentinel
                        if not type(ss)==list:
                            if (ss.find(',')==-1): #There is only ONE tile
                                S = [ss.replace('[','').replace(']','').replace('\'','')]
                            else:
                                tok = ss.split(',')
                                l = len(tok)
                                S=[]
                                for ii in range(l):
                                    tk = tok[ii].replace('[','').replace(']','').replace('\'','')
                                    S.append(tk)
                        else:
                            S = ss
                        if (type(S)==str):
                            S = [S]
                    
                        #do landsat
                        if not type(ll)==list:
                            if (ss.find(',')==-1): #There is only ONE tile
                                L = ll.replace('[','').replace(']','').replace('\'','')
                            else:
                                tok = ll.split(', ')
                                l = len(tok)
                                L=[]
                                for ii in range(l):
                                    tk = tok[ii].replace('[','').replace(']','').replace('\'','')
                                    L.append(tk.replace(' ',''))
                        else:
                            L = [str(x).zfill(6) for x in ll]

                        if (type(L)==str):
                            L = [L]
                        
                        ndvi_mList_df = pd.DataFrame()
                        
                        NDVIFiles_df = getNDVIFilesWithDate(32, S, L, etoDateStr)
                        NDVIFiles_df.reset_index(drop=True,inplace=True)
                        #print(S, L, NDVIFiles_df)
                        
                        [ndvi_df,ndvi_img,ndvi_ext,imgDtStr] = getNDVIValueWithStatsAndImg5(NDVIFiles_df,P)
                        
                        
                        if (len(ndvi_df)>0 & (len(ndvi_img)>0)): #The first time the df has valid values
                            ext = json.loads(ndvi_ext)
                            ndvi_mList_df = ndvi_df.copy()
                            
                            
                            acct_dir_img2 = acct_dir_img 
                            if not os.path.exists(acct_dir_img2):
                                os.makedirs(acct_dir_img2)
                            imgFile = acct_dir_img2 + acctName+'_'+field+'__IMG_'+etoDateStr+'.png'
                            infoJsonFile = acct_dir_img2 + acctName+'_'+field+'__INFO_'+etoDateStr+'.json'

                            acct_dir_pkl2 = acct_dir_pkl 
                            if not os.path.exists(acct_dir_pkl2):
                                os.makedirs(acct_dir_pkl2)
                            pklFile = acct_dir_pkl2 + acctName+'_'+field+'__PKL_'+etoDateStr+'.pkl'
                            
                        
                            imgDict = {
                                'ImgFileName': os.path.split(imgFile)[1],
                                'ImgDate': imgDtStr,
                                'ImgType': 'ETc',
                                'ImgExtent': ndvi_ext, 
                                'ImgBoundary': P.to_wkt(),
                                'SubCategory': row.SubCategory,
                                'Category': row.Category,
                                'Ranch': ranch,
                                'Field': field
                            }
                            #print(infoJsonFile)
                            with open(infoJsonFile, 'w', encoding='utf8') as outfile:
                                json.dump(imgDict, outfile, ensure_ascii=False)

                            plt.ioff();
                            #print(imgFile)
                            with open(imgFile, 'wb') as outfile:
                                fig = plt.figure()
                                plt.imshow(ndvi_img, extent=ext, cmap=zoom1, vmin=0, vmax=1)
                                fig.gca().set_axis_off()
                                plt.savefig(outfile, dpi=600, bbox_inches='tight', pad_inches=0, transparent=True)
                                plt.close();
                                
                            #print(pklFile)
                            with open(pklFile,'wb') as f:
                                pickle.dump(ndvi_img, f)
                        
                            retdf = ndvi_df.loc[0]

                            ndvi_mean = retdf.Mean
                            ndvi_p90 = retdf.p90
                            ndvi_p25 = retdf.p25
                            ndvi_p75 = retdf.p75
                            ndvi_p10 = retdf.p10
                            ndvi_mx = retdf.mx
                            ndvi_mn = retdf.mn
                            ndvi_median = retdf.Median
                            ndvi_std = retdf.StdDev
                            
                        else:
                            #set values if NaN
                            ndvi_mean = np.nan
                            ndvi_mx = np.nan
                            ndvi_mn = np.nan
                            ndvi_p90 = np.nan
                            ndvi_p25 = np.nan
                            ndvi_p75 = np.nan
                            ndvi_p10 = np.nan
                            ndvi_median = np.nan
                            ndvi_std = np.nan
                            
                            
                        #derived
                        ndvi_unif = (ndvi_p25 / ndvi_median)*100
                    
                        #compute the METRIC ETc
                        et_metric = 0.0 #Don't compute METRIC ET due to no rnet
                        #computeMetricET(tavg,rnet,wavg,ndvi_mean)
                        
                        ##########if et_metric is NaN, make it zero
                        if (np.isnan(et_metric)):
                            et_metric=0.0
                    
                        #the Kc
                        fac = 1.35 #somewhat arbitrary for now
                        k_c = 1.3558*ndvi_mean*fac + 0.05
                        k_c2 = 1.3558*ndvi_mean + 0.05
                    
                        #compute the ETc
                        et_ndvi = eto*k_c
                        et_ndvi2 = eto*k_c2
                    
                        #get the final ETc
                        # Include et_metric
                        k_a = 1 #We can reduce k_a to increase the effect of Metric (see above)
                        et_c =  (1.15 - k_a)*et_metric + ((k_a + 0.0)*et_ndvi)
                        et_c2 =  (1.15 - k_a)*et_metric + ((k_a + 0.0)*et_ndvi2)
                        
                        k_c_p90 = ((1.3558*ndvi_p90*fac + 0.05) + (1.3558*ndvi_p90 + 0.05))/2
                        k_c_p75 = ((1.3558*ndvi_p75*fac + 0.05) + (1.3558*ndvi_p75 + 0.05))/2
                        k_c_p25 = ((1.3558*ndvi_p25*fac + 0.05) + (1.3558*ndvi_p25 + 0.05))/2
                        k_c_p10 = ((1.3558*ndvi_p10*fac + 0.05) + (1.3558*ndvi_p10 + 0.05))/2
                        k_c_max = ((1.3558*ndvi_mx*fac + 0.05) + (1.3558*ndvi_mx + 0.05))/2
                        k_c_min = ((1.3558*ndvi_mn*fac + 0.05) + (1.3558*ndvi_mn + 0.05))/2


                        et_p90 = eto*k_c_p90
                        et_p75 = eto*k_c_p75
                        et_p25 = eto*k_c_p25
                        et_p10 = eto*k_c_p10
                        et_max = eto*k_c_max
                        et_min = eto*k_c_min
                        
                        acctGDF_date.loc[i,'ETc']=et_c
                        acctGDF_date.loc[i,'ETc_low']=et_c2
                        acctGDF_date.loc[i,'Kc']=k_c
                        acctGDF_date.loc[i,'ETm']=et_metric
                        acctGDF_date.loc[i,'ET_p90']=et_p90
                        acctGDF_date.loc[i,'ET_unif']=ndvi_unif
                        acctGDF_date.loc[i,'ET_p75']=et_p75
                        acctGDF_date.loc[i,'ET_p25']=et_p25
                        acctGDF_date.loc[i,'ET_p10']=et_p10
                        acctGDF_date.loc[i,'ET_max']=et_max
                        acctGDF_date.loc[i,'ET_min']=et_min

                    acctGDF_date['SentinelFootprint'] = acctGDF_date['SentinelFootprint'].astype('str')
                    acctGDF_date['LandsatFootprint'] = acctGDF_date['LandsatFootprint'].astype('str')
                    acctGDF_date['Date'] = acctGDF_date['Date'].astype('str')
                    if 'CentroidPoint' in acctGDF_date.columns:
                        acctGDF_date.drop('CentroidPoint',inplace=True, axis=1)
                    
                    
                    #save the values
                    #print(jsonFile)
                    acctGDF_date.to_file(jsonFile, driver="GeoJSON")

                    #Print to screen
                    print(acctGDF_date.ETc.tolist())
                    print("\tElapsed --- %s seconds ---" % (time.time() - startd))

                except:
                    print('Exception in Acct: '+str(Exception))
                    continue 

            print("\tEND --- %s seconds ---" % (time.time() - starta))
            print('\tFinally Done!')

        except:
            print('Exception in Acct: '+str(Exception))
            continue

    print("\nEND --- %s seconds ---" % (time.time() - startt))
    return

def generate_EGY_ETc2(merged_df, measureType):

    unique_accounts = pd.unique(merged_df.Category)
    startt = time.time()
    #Now to put this into the output folder for each account
    print("Start --- %s  ---" % datetime.now())

    for cid in unique_accounts:
        a = str(cid)
        try:
            starta = time.time()
            
            print('\n\tWriting...%s ' % a)
            print("\tStart --- %s  ---" % datetime.now())

            acctName = a.replace(' ','').upper()

            acct_dir = pilotAcctOutputPath + a + '/data/'
            acct_dir_img = pilotAcctOutputPath + a + '/img/'
            acct_dir_pkl = pilotAcctOutputPath + a + '/pkl/'
            acct_dir_rep = pilotAcctOutputPath + a + '/reports/'
                
            #create the directory if one doesn't exist
            if not os.path.exists(acct_dir):
                os.makedirs(acct_dir)
            if not os.path.exists(acct_dir_img):
                os.makedirs(acct_dir_img)
            if not os.path.exists(acct_dir_pkl):
                os.makedirs(acct_dir_pkl)
            if not os.path.exists(acct_dir_rep):
                os.makedirs(acct_dir_rep)
                
            #For each account, take the subset of the acct_df 
            #acctDF = merged_df[merged_df.Category == a]
            acctDF = merged_df[merged_df.Category == cid]
            
            acctGDF = gpd.GeoDataFrame(acctDF,geometry=acctDF.geometry)
            
            #for every day, do the following analysis and printout
            acct_dates_list = acctGDF.Date.tolist()
            unique_dts = pd.unique(acct_dates_list) # so we don't multiple count for each field
            unique_dts.sort()
            
            for dt in unique_dts:
                try:
                    endDay = dt        
                    
                    startd = time.time()
                    #print("\n\t\tSTART --- %s  ---" % datetime.now())

                    etoDateStr =  datetime.strftime(dt,'%Y-%m-%d')
                    jsonFile = acct_dir + acctName+'__ET_'+etoDateStr+'.json'
                    
                    acctGDF_date = acctGDF[acctGDF.Date==dt]
                    acctGDF_date.reset_index(drop=True, inplace=True)
                    
                    print('\t\t'+etoDateStr)
                    for i,row in acctGDF_date.iterrows():
                        
                        acct = row.Category
                        ranch = row.SubCategory
                        field = row.Name.replace(' ','')
                        
                        #ETo and weather values
                        eto = row['ETo(in)']
                        
                                    
                        # NDVI value
                        P = getSinglePoly(row.geometry)

                        ss = row.SentinelFootprint
                        ll = row.LandsatFootprint
                        
                        #do Sentinel
                        if not type(ss)==list:
                            if (ss.find(',')==-1): #There is only ONE tile
                                S = [ss.replace('[','').replace(']','').replace('\'','')]
                            else:
                                tok = ss.split(',')
                                l = len(tok)
                                S=[]
                                for ii in range(l):
                                    tk = tok[ii].replace('[','').replace(']','').replace('\'','')
                                    S.append(tk)
                        else:
                            S = ss
                        if (type(S)==str):
                            S = [S]
                    
                        #do landsat
                        if not type(ll)==list:
                            if (ss.find(',')==-1): #There is only ONE tile
                                L = ll.replace('[','').replace(']','').replace('\'','')
                            else:
                                tok = ll.split(', ')
                                l = len(tok)
                                L=[]
                                for ii in range(l):
                                    tk = tok[ii].replace('[','').replace(']','').replace('\'','')
                                    L.append(tk.replace(' ',''))
                        else:
                            L = [str(x).zfill(6) for x in ll]

                        if (type(L)==str):
                            L = [L]
                        
                        ndvi_mList_df = pd.DataFrame()
                        
                        NDVIFiles_df = getNDVIFilesWithDate(32, S, L, etoDateStr)
                        NDVIFiles_df.reset_index(drop=True,inplace=True)
                        #print(S, L, NDVIFiles_df)
                        
                        [ndvi_df,ndvi_img,ndvi_ext,imgDtStr] = getNDVIValueWithStatsAndImg5(NDVIFiles_df,P)
                        
                        
                        if (len(ndvi_df)>0 & (len(ndvi_img)>0)): #The first time the df has valid values
                            ext = json.loads(ndvi_ext)
                            ndvi_mList_df = ndvi_df.copy()
                            
                            
                            acct_dir_img2 = acct_dir_img 
                            if not os.path.exists(acct_dir_img2):
                                os.makedirs(acct_dir_img2)
                            imgFile = acct_dir_img2 + acctName+'_'+field+'__IMG_'+etoDateStr+'.png'
                            infoJsonFile = acct_dir_img2 + acctName+'_'+field+'__INFO_'+etoDateStr+'.json'

                            acct_dir_pkl2 = acct_dir_pkl 
                            if not os.path.exists(acct_dir_pkl2):
                                os.makedirs(acct_dir_pkl2)
                            pklFile = acct_dir_pkl2 + acctName+'_'+field+'__PKL_'+etoDateStr+'.pkl'
                            
                        
                            imgDict = {
                                'ImgFileName': os.path.split(imgFile)[1],
                                'ImgDate': imgDtStr,
                                'ImgType': 'ETc',
                                'ImgExtent': ndvi_ext, 
                                'ImgBoundary': P.to_wkt(),
                                'SubCategory': row.SubCategory,
                                'Category': row.Category,
                                'Ranch': ranch,
                                'Field': field
                            }
                            #print(infoJsonFile)
                            with open(infoJsonFile, 'w', encoding='utf8') as outfile:
                                json.dump(imgDict, outfile, ensure_ascii=False)

                            plt.ioff();
                            #print(imgFile)
                            with open(imgFile, 'wb') as outfile:
                                fig = plt.figure()
                                plt.imshow(ndvi_img, extent=ext, cmap=zoom1, vmin=0, vmax=1)
                                fig.gca().set_axis_off()
                                plt.savefig(outfile, dpi=600, bbox_inches='tight', pad_inches=0, transparent=True)
                                plt.close();
                                
                            #print(pklFile)
                            with open(pklFile,'wb') as f:
                                pickle.dump(ndvi_img, f)
                        
                            retdf = ndvi_df.loc[0]

                            ndvi_mean = retdf.Mean
                            ndvi_p90 = retdf.p90
                            ndvi_p25 = retdf.p25
                            ndvi_p75 = retdf.p75
                            ndvi_p10 = retdf.p10
                            ndvi_mx = retdf.mx
                            ndvi_mn = retdf.mn
                            ndvi_median = retdf.Median
                            ndvi_std = retdf.StdDev
                            
                        else:
                            #set values if NaN
                            ndvi_mean = np.nan
                            ndvi_mx = np.nan
                            ndvi_mn = np.nan
                            ndvi_p90 = np.nan
                            ndvi_p25 = np.nan
                            ndvi_p75 = np.nan
                            ndvi_p10 = np.nan
                            ndvi_median = np.nan
                            ndvi_std = np.nan
                            
                        #print(ndvi_mean)

                        #derived
                        ndvi_unif = (ndvi_p25 / ndvi_median)*100
                    
                        #compute the METRIC ETc
                        et_metric = 0.0 #Don't compute METRIC ET due to no rnet
                        #computeMetricET(tavg,rnet,wavg,ndvi_mean)
                        
                        ##########if et_metric is NaN, make it zero
                        if (np.isnan(et_metric)):
                            et_metric=0.0
                    
                        #the Kc
                        fac = 1.35 #somewhat arbitrary for now
                        k_c = 1.3558*ndvi_mean*fac + 0.05
                        k_c2 = 1.3558*ndvi_mean + 0.05
                    
                        #compute the ETc
                        et_ndvi = eto*k_c
                        et_ndvi2 = eto*k_c2
                    
                        #get the final ETc
                        # Include et_metric
                        k_a = 1 #We can reduce k_a to increase the effect of Metric (see above)
                        et_c =  (1.15 - k_a)*et_metric + ((k_a + 0.0)*et_ndvi)
                        et_c2 =  (1.15 - k_a)*et_metric + ((k_a + 0.0)*et_ndvi2)
                        
                        k_c_p90 = ((1.3558*ndvi_p90*fac + 0.05) + (1.3558*ndvi_p90 + 0.05))/2
                        k_c_p75 = ((1.3558*ndvi_p75*fac + 0.05) + (1.3558*ndvi_p75 + 0.05))/2
                        k_c_p25 = ((1.3558*ndvi_p25*fac + 0.05) + (1.3558*ndvi_p25 + 0.05))/2
                        k_c_p10 = ((1.3558*ndvi_p10*fac + 0.05) + (1.3558*ndvi_p10 + 0.05))/2
                        k_c_max = ((1.3558*ndvi_mx*fac + 0.05) + (1.3558*ndvi_mx + 0.05))/2
                        k_c_min = ((1.3558*ndvi_mn*fac + 0.05) + (1.3558*ndvi_mn + 0.05))/2


                        et_p90 = eto*k_c_p90
                        et_p75 = eto*k_c_p75
                        et_p25 = eto*k_c_p25
                        et_p10 = eto*k_c_p10
                        et_max = eto*k_c_max
                        et_min = eto*k_c_min
                        
                        acctGDF_date.loc[i,'ETc']=et_c
                        acctGDF_date.loc[i,'ETc_low']=et_c2
                        acctGDF_date.loc[i,'Kc']=k_c
                        acctGDF_date.loc[i,'ETm']=et_metric
                        acctGDF_date.loc[i,'ET_p90']=et_p90
                        acctGDF_date.loc[i,'ET_unif']=ndvi_unif
                        acctGDF_date.loc[i,'ET_p75']=et_p75
                        acctGDF_date.loc[i,'ET_p25']=et_p25
                        acctGDF_date.loc[i,'ET_p10']=et_p10
                        acctGDF_date.loc[i,'ET_max']=et_max
                        acctGDF_date.loc[i,'ET_min']=et_min

                    acctGDF_date['SentinelFootprint'] = acctGDF_date['SentinelFootprint'].astype('str')
                    acctGDF_date['LandsatFootprint'] = acctGDF_date['LandsatFootprint'].astype('str')
                    acctGDF_date['Date'] = acctGDF_date['Date'].astype('str')
                    if 'CentroidPoint' in acctGDF_date.columns:
                        acctGDF_date.drop('CentroidPoint',inplace=True, axis=1)
                    
                    
                    #save the values
                    #print(jsonFile)
                    acctGDF_date.to_file(jsonFile, driver="GeoJSON")

                    #Print to screen
                    print(acctGDF_date.ETc.tolist())
                    print("\tElapsed --- %s seconds ---" % (time.time() - startd))

                except:
                    print('Exception in Acct: '+str(Exception))
                    continue 

            print("\tEND --- %s seconds ---" % (time.time() - starta))
            print('\tFinally Done!')

        except:
            print('Exception in Acct: '+str(Exception))
            continue

    print("\nEND --- %s seconds ---" % (time.time() - startt))
    return

def generate_IND_KPIs(acct_gdf, kpi, unique_dts):

    unique_accounts = pd.unique(acct_gdf.Category)
    startt = time.time()
    #Now to put this into the output folder for each account
    print("Start --- %s  ---" % datetime.now())

    for acctName in unique_accounts:
        
        starta = time.time()

        acct = acctName.upper().replace(' ','')
        
        print('\n\tWriting...%s ' % acct)
        print("\tStart --- %s  ---" % datetime.now())

        acct_dir = pilotAcctOutputPath + acct + '/data/'
        acct_dir_img = pilotAcctOutputPath + acct + '/img/'
        acct_dir_pkl = pilotAcctOutputPath + acct + '/pkl/'
        acct_dir_rep = pilotAcctOutputPath + acct + '/reports/'
            
        #create the directory if one doesn't exist
        if not os.path.exists(acct_dir):
            os.makedirs(acct_dir)
        if not os.path.exists(acct_dir_img):
            os.makedirs(acct_dir_img)
        if not os.path.exists(acct_dir_pkl):
            os.makedirs(acct_dir_pkl)
        if not os.path.exists(acct_dir_rep):
            os.makedirs(acct_dir_rep)
            
        acctDF = acct_gdf[acct_gdf.Category == acctName]
        acctDF.reset_index(drop=True, inplace=True)
        acctGDF = gpd.GeoDataFrame(acctDF,geometry=acctDF.geometry)
        
        
        for dt in unique_dts:
            
            endDay = dt        
            startd = time.time()
            dtStr =  datetime.strftime(dt,'%Y-%m-%d')
            jsonFile = acct_dir + acct+'__NII_'+dtStr+'.json'
            
            #acctGDF_date = acctGDF[acctGDF.Date==dt]
            #acctGDF_date.reset_index(drop=True, inplace=True)
            
            print('\t\t'+dtStr)
            for i,row in acctGDF.iterrows():
                ranch = row.SubCategory
                field = row.Name				
            
                # NII value
                P = getSinglePoly(row.geometry)
                
                ss = row.SentinelFootprint
                ll = row.LandsatFootprint
                
                #do Sentinel
                if not type(ss)==list:
                    if (ss.find(',')==-1): #There is only ONE tile
                        S = [ss.replace('[','').replace(']','').replace('\'','')]
                    else:
                        tok = ss.split(',')
                        l = len(tok)
                        S=[]
                        for ii in range(l):
                            tk = tok[ii].replace('[','').replace(']','').replace('\'','')
                            S.append(tk)
                else:
                    S = ss
                if (type(S)==str):
                    S = [S]
            
                #do landsat
                if not type(ll)==list:
                    if (ss.find(',')==-1): #There is only ONE tile
                        L = ll.replace('[','').replace(']','').replace('\'','')
                    else:
                        tok = ll.split(', ')
                        l = len(tok)
                        L=[]
                        for ii in range(l):
                            tk = tok[ii].replace('[','').replace(']','').replace('\'','')
                            L.append(tk.replace(' ',''))
                else:
                    L = [str(x).zfill(6) for x in ll]

                if (type(L)==str):
                    L = [L]
                
                nii_mList_df = pd.DataFrame()
                
                NIIFiles_df = getNIIFilesWithDate(32, S, L, dtStr)
                NIIFiles_df.reset_index(drop=True,inplace=True)
                #print(NIIFiles_df)
                
                [nii_df,nii_img,nii_ext,imgDtStr] = getNIIValueWithStatsAndImg3(NIIFiles_df,P)
                
                
                if (len(nii_df)>0 & (len(nii_img)>0)): #The first time the df has valid values
                    ext = json.loads(nii_ext)
                    nii_mList_df = nii_df.copy()
                    
                    
                    acct_dir_img2 = acct_dir_img
                    if not os.path.exists(acct_dir_img2):
                        os.makedirs(acct_dir_img2)

                    imgFile = acct_dir_img2 + acctName+'_'+field+'__IMG_'+dtStr+'.png'
                    infoJsonFile = acct_dir_img2 + acctName+'_'+field+'__INFO_'+dtStr+'.json'

                    acct_dir_pkl2 = acct_dir_pkl
                    if not os.path.exists(acct_dir_pkl2):
                        os.makedirs(acct_dir_pkl2)
                    pklFile = acct_dir_pkl2 + acctName+'_'+field+'__PKL_'+dtStr+'.pkl'
                    
                
                    imgDict = {
                        'ImgFileName': os.path.split(imgFile)[1],
                        'ImgDate': imgDtStr,
                        'ImgType': 'NII',
                        'ImgExtent': nii_ext, 
                        'ImgBoundary': P.to_wkt(),
                        'Acct': acct,
                        'Ranch': ranch,
                        'Field': field
                    }
                    with open(infoJsonFile, 'w', encoding='utf8') as outfile:
                        json.dump(imgDict, outfile, ensure_ascii=False)

                    plt.ioff();
                    with open(imgFile, 'wb') as outfile:
                        fig = plt.figure()
                        plt.imshow(nii_img, extent=ext, cmap=zoom1)
                        fig.gca().set_axis_off()
                        plt.savefig(outfile, dpi=600, bbox_inches='tight', pad_inches=0, transparent=True)
                        plt.close();
                        
                    with open(pklFile,'wb') as f:
                        pickle.dump(nii_img, f)
                
                    retdf = nii_df.loc[0]

                    nii_mean = retdf.Mean
                    nii_p90 = retdf.p90
                    nii_p25 = retdf.p25
                    nii_p75 = retdf.p75
                    nii_median = retdf.Median
                    nii_std = retdf.StdDev
                    
                else:
                    #set values if NaN
                    nii_mean = np.nan
                    nii_p90 = np.nan
                    nii_p25 = np.nan
                    nii_p75 = np.nan
                    nii_median = np.nan
                    nii_std = np.nan
                    
                    
                #derived
                nii_unif = (nii_p25 / nii_median)*100
                acctGDF.loc[i,'Date']=dtStr
                acctGDF.loc[i,'NII_mean']=nii_mean
                acctGDF.loc[i,'NII_p25']=nii_p25
                acctGDF.loc[i,'NII_p75']=nii_p75
                acctGDF.loc[i,'NII_p90']=nii_p90
                acctGDF.loc[i,'NII_median']=nii_median
                acctGDF.loc[i,'NII_std']=nii_std
                acctGDF.loc[i,'NII_unif']=nii_unif
            
    
            acctGDF['SentinelFootprint'] = acctGDF['SentinelFootprint'].astype('str')
            acctGDF['Date'] = acctGDF['Date'].astype('str')
            if 'CentroidPoint' in acctGDF.columns:
                acctGDF.drop('CentroidPoint',inplace=True, axis=1)
            
            
            #save the values
            acctGDF.to_file(jsonFile, driver="GeoJSON")

            #Print to screen
            print(acctGDF.NII_mean.tolist())
            print("\tElapsed --- %s seconds ---" % (time.time() - startd))

        print("\tEND --- %s seconds ---" % (time.time() - starta))
        print('\tDone!')
        #except:
        #	print('Exception in Acct: '+str(Exception))
        #	continue

    print("\nEND --- %s seconds ---" % (time.time() - startt))
    return


########################  AZ Functions
def getWeatherForDates_NOAA(p,list_of_dates,df):
    
    N=5 #default
    #gets the weather metrics for the METRIC calcs for a point for each of the dates based
    #of the closest N stns #given the geopandas geodataframe gpd
    
    #Screen out the null values

    df = df[df['Tavg']!='--']
    df = df[df['Rnet']!='--']
    df = df[df['Wavg']!='--']
    df = df[df['ETo_average_inches']!='--']
    df = df[df['Precip_in']!='--'] 
    df = df[df['T_max']!='--']
    df = df[df['T_min']!='--']
    df = df[df['T_mean']!='--']
    df = df[df['RH_max']!='--']
    df = df[df['RH_min']!='--']
    df = df[df['RH_mean']!='--']
    df = df[df['T_dew']!='--']
    
    
    cols_list = ['Date',
                 'Tavg',
                 'Wavg',
                 'ETo_average_inches',
                 'Precip_in',
                 'T_max',
                 'T_min',
                 'T_mean',
                 'RH_max',
                 'RH_min',
                 'RH_mean',
                 'T_dew'
                ]
    df_weather_for_dates = pd.DataFrame(columns=cols_list)
    df.Date=pd.to_datetime(df.Date)
    

    idxx=0
    for dt in list_of_dates:
        #df_dt = df[df.index==dt]
        df_dt = df[df.Date==dt]
        
        #df_dt['StationId']=df_dt['StationId'].astype('str')
        stnPoints = list(pd.unique(df_dt.geometry))

        stns = pd.unique(df_dt.StationId).tolist()
        stns.sort()
        
        #first get the nearest stns and weights for the point provided
        [ids, pts, dis, wts] = getNearestNStnIDsInOrder2(p,stnPoints,5)
        
        #print('getWeatherForDates_NOAA ids:')
        #print(ids)
        nearestStnIDs = [stns[i] for i in ids]
        #print(nearestStnIDs)

        airTemp = []
        solarRad = []
        windSp = []
        eto=[]
        precip=[]
        maxAirTemp=[]
        minAirTemp=[]
        avgAirTemp=[]
        maxRelHumPercent=[]
        minRelHumPercent=[]
        avgRelHumPercent=[]
        dewPoint=[]
                
        for stn in nearestStnIDs:
            airTemp.append(float(df_dt[df_dt.StationId==stn]['Tavg'].values.tolist()[0]))
            solarRad.append(float(df_dt[df_dt.StationId==stn]['Rnet'].values.tolist()[0]))
            windSp.append(float(df_dt[df_dt.StationId==stn]['Wavg'].values.tolist()[0]))
            eto.append(float(df_dt[df_dt.StationId==stn]['ETo_average_inches'].values.tolist()[0]))
            precip.append(float(df_dt[df_dt.StationId==stn]['Precip_in'].values.tolist()[0]))
            maxAirTemp.append(float(df_dt[df_dt.StationId==stn]['T_max'].values.tolist()[0]))
            minAirTemp.append(float(df_dt[df_dt.StationId==stn]['T_min'].values.tolist()[0]))
            avgAirTemp.append(float(df_dt[df_dt.StationId==stn]['T_mean'].values.tolist()[0]))
            maxRelHumPercent.append(float(df_dt[df_dt.StationId==stn]['RH_max'].values.tolist()[0]))
            minRelHumPercent.append(float(df_dt[df_dt.StationId==stn]['RH_min'].values.tolist()[0]))
            avgRelHumPercent.append(float(df_dt[df_dt.StationId==stn]['RH_mean'].values.tolist()[0]))
            dewPoint.append(float(df_dt[df_dt.StationId==stn]['T_dew'].values.tolist()[0]))
        
        weighted_airTemp = np.nansum(np.array(airTemp)*np.array(wts))
        weighted_solarRad = np.nansum(np.array(solarRad)*np.array(wts))
        weighted_windSp = np.nansum(np.array(windSp)*np.array(wts))
        weighted_eto = np.nansum(np.array(eto)*np.array(wts))
        weighted_precip = np.nansum(np.array(precip)*np.array(wts))
        weighted_maxAirTemp = np.nansum(np.array(maxAirTemp)*np.array(wts))
        weighted_minAirTemp = np.nansum(np.array(minAirTemp)*np.array(wts))
        weighted_avgAirTemp = np.nansum(np.array(avgAirTemp)*np.array(wts))
        weighted_maxRelHumPercent = np.nansum(np.array(maxRelHumPercent)*np.array(wts))
        weighted_minRelHumPercent = np.nansum(np.array(minRelHumPercent)*np.array(wts))
        weighted_avgRelHumPercent = np.nansum(np.array(avgRelHumPercent)*np.array(wts))
        weighted_dewPoint = np.nansum(np.array(dewPoint)*np.array(wts))
        
        df_weather_for_dates.loc[idxx,'Date']=dt
        df_weather_for_dates.loc[idxx,'Tavg']=weighted_airTemp
        df_weather_for_dates.loc[idxx,'Rnet']=weighted_solarRad
        df_weather_for_dates.loc[idxx,'Wavg']=weighted_windSp
        df_weather_for_dates.loc[idxx,'ETo_average_inches']=weighted_eto
        df_weather_for_dates.loc[idxx,'Precip_in']=weighted_precip
        df_weather_for_dates.loc[idxx,'T_max']=weighted_maxAirTemp
        df_weather_for_dates.loc[idxx,'T_min']=weighted_minAirTemp
        df_weather_for_dates.loc[idxx,'T_mean']=weighted_avgAirTemp
        df_weather_for_dates.loc[idxx,'RH_max']=weighted_maxRelHumPercent
        df_weather_for_dates.loc[idxx,'RH_min']=weighted_minRelHumPercent
        df_weather_for_dates.loc[idxx,'RH_mean']=weighted_avgRelHumPercent
        df_weather_for_dates.loc[idxx,'T_dew']=weighted_dewPoint
        
        idxx=idxx+1
          
    return df_weather_for_dates

def getWeatherForDates_NOAA2(p,list_of_dates,df):
    N=5 #default
    #of the closest N stns #given the geopandas geodataframe gpd

    df = df[df['ReferenceETo(in)']!='--']

    cols_list = ['Date','ReferenceETo(in)']
    df_weather_for_dates = pd.DataFrame(columns=cols_list)
    df.Date=pd.to_datetime(df.Date)

    idxx=0
    for dt in list_of_dates:
        df_dt = df[df.Date==dt]

        stnPoints = list(pd.unique(df_dt.geometry))

        stns = pd.unique(df_dt.StationId).tolist()
        stns.sort()

        #first get the nearest stns and weights for the point provided
        [ids, pts, dis, wts] = getNearestNStnIDsInOrder2(p,stnPoints,5)
        
        nearestStnIDs = [stns[i] for i in ids]

        eto=[]

        for stn in nearestStnIDs:
            eto.append(float(df_dt[df_dt.StationId==stn]['ReferenceETo(in)'].values.tolist()[0]))

       
        weighted_eto = np.nansum(np.array(eto)*np.array(wts))

        df_weather_for_dates.loc[idxx,'Date']=dt
        df_weather_for_dates.loc[idxx,'ReferenceETo(in)']=weighted_eto
        
        idxx=idxx+1
        
    
    df_weather_for_dates = df_weather_for_dates.rename(columns={'ReferenceETo(in)':'ETo_average_inches'})

    return df_weather_for_dates
    

def generate_NOAA_ETo3(acct_gdf, stns_gdf, startDtStr, endDtStr, obsPathStr, region):
    
    dailyCols = ['STN',
             'TS',
             'JulianDate',
             'ETo_average_inches',
             'Precip_in',
             'R_n__MJpm2',
             'P',
             'T_mean',
             'T_min',
             'T_max',
             'RH_mean',
             'RH_max',
             'RH_min',
             'T_dew',
             'Uz'
            ]
    dailyCols = [x.replace(' ','') for x in dailyCols]
    dailyCols = [x.replace('∞','') for x in dailyCols]
    dailyCols = [x.replace('=','') for x in dailyCols]

    dailyHistoricalWeatherPD = pd.DataFrame(columns=dailyCols)
    data = pd.DataFrame(columns = dailyCols)
    #read all the files into a big pandas dataframe
    dtFiles = []


    startYr = startDtStr.split('-')[0]
    endYr = endDtStr.split('-')[0]
    years = list(set([startYr,endYr]))
    years.sort()

    #print(years)

    for yr in years:
        annualStnFiles = [fname for fname in os.listdir(obsPathStr) if (fname.endswith('.geojson') & (fname.split('__')[1].split('.')[0].startswith(yr)) )]
        annualStnFiles.sort()

        #print(annualStnFiles)

        for fl in annualStnFiles:
            f = obsPathStr + fl
            #print(f)
            data = gpd.read_file(f)
            flDtStr = fl.split('.')[0].split('_')[-1]
            flDt = datetime.strptime(flDtStr,'%Y-%m-%d')
            #data.columns = dailyCols
            #fix the Date to make it a Datetime obj
            #print(data.columns)

            if 'TS' in data.columns:
                data['Date'] =  pd.to_datetime(data['TS'],format='%Y-%m-%d')
            else:
                data['Date'] =  pd.to_datetime(data['Date'],format='%Y-%m-%d')

            if not data['Date'].dt==flDt:
                data['Date']=flDt

            if not 'ETo_average_inches' in data.columns:
                data['ETo_average_inches'] =  data['ETo_AVG_IN']
            dailyHistoricalWeatherPD = dailyHistoricalWeatherPD.append(data, ignore_index=True)
    
    #print(dailyHistoricalWeatherPD.columns)
    #print(dailyHistoricalWeatherPD[['Date', 'ETo_average_inches']])

    #Screen these dates
    dailyHistoricalWeatherPD.Date = pd.to_datetime(dailyHistoricalWeatherPD.Date)
    #print([startDtStr, endDtStr])
    dailyHistoricalWeatherPD = dailyHistoricalWeatherPD[(dailyHistoricalWeatherPD.Date >= startDtStr) & (dailyHistoricalWeatherPD.Date <= endDtStr)]
    dailyHistoricalWeatherPD.reset_index(drop=True, inplace=True)

    #print(dailyHistoricalWeatherPD)

    #Sort by Julian Date
    dailyHistoricalWeatherPD['JulianDate'] = dailyHistoricalWeatherPD.Date.dt.dayofyear
    weather_df =  dailyHistoricalWeatherPD.iloc[dailyHistoricalWeatherPD['JulianDate'].astype(int).argsort()]
    weather_df.reset_index(inplace=True,drop=True)
    dmw_pd = weather_df.copy()

    
    #Put values for metric system calculations
    dmw_pd['Tavg']=dmw_pd['T_mean'].values.astype('float') #convert to DegC
    dmw_pd['Rnet']=dmw_pd['R_n__MJpm2'].values.astype('float')
    dmw_pd['Wavg']=dmw_pd['Uz'].values.astype('float')

    #Screen out the null values
    dmw_pd['ReferenceETo(in)'] = dmw_pd['ETo_average_inches'].replace('--',np.nan, regex=True).astype('float')
    dmw_pd['AverageAirTemperature(F)'] = dmw_pd['T_mean'].replace('--',np.nan, regex=True).astype('float')
    dmw_pd['SolarRadiationAverage(Ly/day)'] = dmw_pd['R_n__MJpm2'].replace('--',np.nan, regex=True).astype('float') * 23.900574 #convert to Mangley /day
    dmw_pd['AverageWindSpeed(mph)'] = dmw_pd['Uz'].replace('--',np.nan, regex=True).astype('float') * 2.24 #convert to Mph

    dmw_pd['StationId'] = dmw_pd['STN']
    dmw_pd['Precipitation(in)'] = dmw_pd['Precip_in'].replace('--',np.nan, regex=True).astype('float')
    dmw_pd['AverageAirTemperature(F)'] = dmw_pd['T_mean'].replace('--',np.nan, regex=True).astype('float') * (9/5) + 32
    dmw_pd['MinimumAirTemperature(F)'] = dmw_pd['T_min'].replace('--',np.nan, regex=True).astype('float')* (9/5) + 32
    dmw_pd['MaximumAirTemperature(F)'] = dmw_pd['T_max'].replace('--',np.nan, regex=True).astype('float')* (9/5) + 32
    dmw_pd['AverageRelativeHumidity(%)'] = dmw_pd['RH_mean'].replace('--',np.nan, regex=True).astype('float')
    dmw_pd['DewPoint(F)'] = dmw_pd['T_dew'].replace('--',np.nan, regex=True).astype('float')* (9/5) + 32
    #print(dmw_pd)
    #print(dmw_pd.columns)

    keep_eto_cols_list = ['Date', 
                          'StationId', 
                          'ReferenceETo(in)',
                          'Precipitation(in)', 
                          'MaximumAirTemperature(F)',
                          'MinimumAirTemperature(F)',
                          'AverageAirTemperature(F)',
                          'AverageRelativeHumidity(%)',
                          'DewPoint(F)',
                          'ETo_average_inches',
                          'Precip_in','R_n__MJpm2','P','T_mean','T_min','T_max','RH_mean','RH_max','RH_min','T_dew','Uz',
                         ]
    eto = dmw_pd[keep_eto_cols_list]

    

    keep_metric_cols_list = ['Date','StationId','Tavg','Rnet','Wavg']
    metric = dmw_pd[keep_metric_cols_list] 

    #for each station, we have to interpolate the timeseries and fill in the values
    #first for ETo
    new_list_of_active_stns = []
    list_of_stns = pd.unique(eto.StationId)
    eto_interp = pd.DataFrame()
    #print(list_of_stns)
    for stn in list_of_stns:
        eto_stn = eto[eto.StationId==stn]
        eto_stn = eto_stn.set_index(['Date'])
        eto_stn = eto_stn.sort_index()
        n = len(eto_stn[eto_stn['ReferenceETo(in)'].isnull()==True])
        if n==0:
            eto_interp = eto_interp.append(eto_stn)
        else:
            if (n>=100):
                print([n,stn])
                #dropping these
            else:
                print('\tInterpolating for '+str(n)+' points for stn '+ str(stn))
                new_list_of_active_stns.append(stn)
                eto_stn_interpolate = eto_stn.interpolate(method='time')
                eto_interp = eto_interp.append(eto_stn_interpolate)

    stns_df = stns_gdf.copy()
    stns_df = stns_df.rename(columns={"stationIdentifier": "StationId"})
    eto_interp = eto_interp.reset_index()

    #print(eto_interp)
    #print(stns_df)

    if 'stnid' in stns_df.columns:
        stns_df = stns_df.rename(columns={'stnid':'StationId'})

    df = pd.merge(eto_interp,stns_df,on='StationId')
    df_metric = pd.merge(df,metric,how='left', left_on=['StationId','Date'], right_on = ['StationId','Date'])
    #convert into a GeoDataFrame
    crs = {'init': 'epsg:4326'}
    gdf_metric = gpd.GeoDataFrame(df_metric, crs=crs, geometry=df_metric.geometry)

    #print(gdf_metric)

    historicWeather_gdf =  gdf_metric.copy()
    historicWeather_gdf.reset_index(drop=True, inplace=True)
    historicWeather_gdf['Wavg']=historicWeather_gdf['Wavg'].values.astype('float')
    historicWeather_gdf['Rnet']=historicWeather_gdf['Rnet'].values.astype('float')
    historicWeather_gdf['Uz']=historicWeather_gdf['Uz'].values.astype('float')
    historicWeather_gdf['Precip_in']=historicWeather_gdf['Precip_in'].values.astype('float')

    #print(historicWeather_gdf)
    #print(historicWeather_gdf.columns)

    #Set the horizon from the inputs
    deltaT = datetime.strptime(startDtStr,'%Y-%m-%d') - datetime.strptime(endDtStr,'%Y-%m-%d')
    list_of_dates = pd.date_range(pd.datetime.strptime(startDtStr,'%Y-%m-%d'), periods=-deltaT.days+1).to_pydatetime().tolist()

    startt = time.time()
    EToForAccounts = pd.DataFrame()
    for i,row in acct_gdf.iterrows():
        print('\t\tCalculating ETo for'+
            ', Acct: ' + row.Category+
            ', Ranch: '+row.SubCategory+
            ', Field: '+row.Name)
        start_time = time.time()
        p = loads(row.CentroidPoint) #shapely Point obj

        #print(pd.unique(historicWeather_gdf.Date))
        
        etv = getWeatherForDates_NOAA(p,list_of_dates,historicWeather_gdf)

        #print(etv)

        #now add the columns for the other fields to etv
        etv['UUID']=row.UUID
        etv['Category']=row.Category
        etv['SubCategory']=row.SubCategory
        etv['Name']=row.Name

        #print(etv)

        EToForAccounts = EToForAccounts.append(etv, ignore_index = True)
        #print(EToForAccounts)
        print("\t\tEND ETo--- %s seconds ---" % (time.time() - start_time))
    print("\tETo Done --- %s seconds ---" % (time.time() - startt))
    print('\tDone!')
    

    EToForAccounts = EToForAccounts.rename(columns={'ETo_average_inches': 'ETo(in)', 
                                                        'Precip_in': 'Precipitation(in)',
                                                        'T_max': 'MaximumAirTemperature(F)',
                                                        'T_min': 'MinimumAirTemperature(F)',
                                                        'T_mean': 'AverageAirTemperature(F)',
                                                        'RH_max': 'MaximumRelativeHumidity(%)',
                                                        'RH_min': 'MinimumRelativeHumidity(%)',
                                                        'RH_mean': 'AverageRelativeHumidity(%)',
                                                        'T_dew': 'DewPoint(F)'
                                                       })


    #Merge the acct
    #print('EToForAccounts.columns '+EToForAccounts.columns)
    merged_df = pd.merge(acct_gdf, EToForAccounts, on=['Category','SubCategory','Name','UUID'])
    merged_df.reset_index() #reset the index to avoid any export issues
    merged_df['Date'] = pd.to_datetime(merged_df['Date'])
    
    #print(merged_df)

    return merged_df


def generate_NOAA_ETo2(acct_gdf, stns_gdf, startDtStr, endDtStr, obsPathStr, region):
    dailyCols = ['Date', 'TAVG', 'TMAX', 'TMIN', 'StnID', 'lat', 'lon', 'elev', 'desc',
           'ETo_FAO_MM', 'ETo_HAR_MM', 'ETo_AVG_MM', 'ETo_FAO_IN', 'ETo_HAR_IN',
           'ETo_AVG_IN', 'aifstime_utc', 'local_date_time_full', 'PRCP',
           'geometry', 'ETo_average_inches']
    dailyCols = [x.replace(' ','') for x in dailyCols]
    dailyCols = [x.replace('∞','') for x in dailyCols]
    dailyCols = [x.replace('=','') for x in dailyCols]

    dailyHistoricalWeatherPD = pd.DataFrame(columns=dailyCols)
    data = pd.DataFrame(columns = dailyCols)
    #read all the files into a big pandas dataframe
    dtFiles = []


    startYr = startDtStr.split('-')[0]
    endYr = endDtStr.split('-')[0]
    years = list(set([startYr,endYr]))
    years.sort()

    #print(years)

    for yr in years:
        annualStnFiles = [fname for fname in os.listdir(obsPathStr) if (fname.endswith('.geojson') & (fname.split('__')[1].split('.')[0].startswith(yr)) & (fname.split('_')[0]==region.upper()) )]
        annualStnFiles.sort()

        #print(annualStnFiles)

        for fl in annualStnFiles:
            f = obsPathStr + fl
            #print(f)
            data = gpd.read_file(f)
            flDtStr = fl.split('.')[0].split('_')[-1]
            flDt = datetime.strptime(flDtStr,'%Y-%m-%d')
            #data.columns = dailyCols
            #fix the Date to make it a Datetime obj
            #print(data.columns)

            if 'TS' in data.columns:
                data['Date'] =  pd.to_datetime(data['TS'],format='%Y-%m-%d')
            else:
                data['Date'] =  pd.to_datetime(data['Date'],format='%Y-%m-%d')

            if not data['Date'].dt==flDt:
                data['Date']=flDt

            if not 'ETo_average_inches' in data.columns:
                data['ETo_average_inches'] =  data['ETo_AVG_IN']
            dailyHistoricalWeatherPD = dailyHistoricalWeatherPD.append(data, ignore_index=True)


    #Screen these dates
    dailyHistoricalWeatherPD.Date = pd.to_datetime(dailyHistoricalWeatherPD.Date)
    startDt = datetime.strptime(startDtStr, '%Y-%m-%d')
    endDt = datetime.strptime(endDtStr, '%Y-%m-%d')
    dailyHistoricalWeatherPD = dailyHistoricalWeatherPD[(dailyHistoricalWeatherPD.Date >= startDt) & (dailyHistoricalWeatherPD.Date <= endDt)]
    dailyHistoricalWeatherPD.reset_index(drop=True, inplace=True)

    #Sort by Julian Date
    dailyHistoricalWeatherPD['JulianDate'] = dailyHistoricalWeatherPD.Date.dt.dayofyear
    weather_df =  dailyHistoricalWeatherPD.iloc[dailyHistoricalWeatherPD['JulianDate'].astype(int).argsort()]
    weather_df.reset_index(inplace=True,drop=True)
    dmw_pd = weather_df.copy()

    #Screen out the null values
    dmw_pd['ReferenceETo(in)'] = dmw_pd['ETo_average_inches'].replace('--',np.nan, regex=True).astype('float')
    dmw_pd['StationId'] = dmw_pd['StnID']

    keep_eto_cols_list = ['Date', 
                          'StationId', 
                          'ReferenceETo(in)',
                         ]
    eto = dmw_pd[keep_eto_cols_list]

    #for each station, we have to interpolate the timeseries and fill in the values
    #first for ETo
    new_list_of_active_stns = []
    list_of_stns = pd.unique(eto.StationId)
    eto_interp = pd.DataFrame()
    #print(list_of_stns)
    for stn in list_of_stns:
        eto_stn = eto[eto.StationId==stn]
        eto_stn = eto_stn.set_index(['Date'])
        eto_stn = eto_stn.sort_index()
        n = len(eto_stn[eto_stn['ReferenceETo(in)'].isnull()==True])
        if n==0:
            eto_interp = eto_interp.append(eto_stn)
        else:
            if (n>=100):
                print([n,stn])
                #dropping these
            else:
                print('\tInterpolating for '+str(n)+' points for stn '+ str(stn))
                new_list_of_active_stns.append(stn)
                eto_stn_interpolate = eto_stn.interpolate(method='time')
                eto_interp = eto_interp.append(eto_stn_interpolate)

    stns_df = stns_gdf.copy()
    eto_interp = eto_interp.reset_index()


    if 'stnid' in stns_df.columns:
        stns_df = stns_df.rename(columns={'stnid':'StationId'})

    df = pd.merge(eto_interp,stns_df,on='StationId')
    #convert into a GeoDataFrame
    crs = {'init': 'epsg:4326'}
    gdf_metric = gpd.GeoDataFrame(df, crs=crs, geometry=df.geometry)

    #print(gdf_metric)

    historicWeather_gdf =  gdf_metric.copy()

    #print(historicWeather_gdf)
    #print(historicWeather_gdf.columns)

    #Set the horizon from the inputs
    deltaT = datetime.strptime(startDtStr,'%Y-%m-%d') - datetime.strptime(endDtStr,'%Y-%m-%d')
    list_of_dates = pd.date_range(pd.datetime.strptime(startDtStr,'%Y-%m-%d'), periods=-deltaT.days+1).to_pydatetime().tolist()

    startt = time.time()
    EToForAccounts = pd.DataFrame()
    for i,row in acct_gdf.iterrows():
        print('\t\tCalculating ETo for'+
            ', Acct: ' + row.Category+
            ', Ranch: '+row.SubCategory+
            ', Field: '+row.Name)
        start_time = time.time()
        p = loads(row.CentroidPoint) #shapely Point obj

        #print(pd.unique(historicWeather_gdf.Date))

        etv = getWeatherForDates_NOAA2(p,list_of_dates,historicWeather_gdf)

        #print(etv)

        #now add the columns for the other fields to etv
        etv['UUID']=row.UUID
        etv['Category']=row.Category
        etv['SubCategory']=row.SubCategory
        etv['Name']=row.Name

        #print(etv)

        EToForAccounts = EToForAccounts.append(etv, ignore_index = True)
        #print(EToForAccounts)
        print("\t\tEND ETo--- %s seconds ---" % (time.time() - start_time))
    print("\tETo Done --- %s seconds ---" % (time.time() - startt))
    print('\tDone!')


    EToForAccounts = EToForAccounts.rename(columns={'ETo_average_inches': 'ETo(in)'})


    #Merge the acct
    #print('EToForAccounts.columns '+EToForAccounts.columns)
    merged_df = pd.merge(acct_gdf, EToForAccounts, on=['Category','SubCategory','Name','UUID'])
    merged_df.reset_index() #reset the index to avoid any export issues
    merged_df['Date'] = pd.to_datetime(merged_df['Date'])

    return merged_df

def generate_NOAA_ETo(acct_gdf, stns_gdf, startDtStr, endDtStr, obsPathStr):
    
    dailyCols = ['STN',
             'TS',
             'JulianDate',
             'ETo_average_inches',
             'Precip_in',
             'R_n__MJpm2',
             'P',
             'T_mean',
             'T_min',
             'T_max',
             'RH_mean',
             'RH_max',
             'RH_min',
             'T_dew',
             'Uz'
            ]
    dailyCols = [x.replace(' ','') for x in dailyCols]
    dailyCols = [x.replace('∞','') for x in dailyCols]
    dailyCols = [x.replace('=','') for x in dailyCols]

    dailyHistoricalWeatherPD = pd.DataFrame(columns=dailyCols)
    data = pd.DataFrame(columns = dailyCols)
    #read all the files into a big pandas dataframe
    dtFiles = []


    startYr = startDtStr.split('-')[0]
    endYr = endDtStr.split('-')[0]
    years = list(set([startYr,endYr]))
    years.sort()
    for yr in years:
        annualStnFiles = [fname for fname in os.listdir(obsPathStr) if (fname.endswith('.geojson') & fname.split('__')[1].split('.')[0].startswith(yr))]
        for fl in annualStnFiles:
            f = obsPathStr + fl
            #print(f)
            data = gpd.read_file(f)
            #data.columns = dailyCols
            #fix the Date to make it a Datetime obj
            data['Date'] =  pd.to_datetime(data['TS'],format='%Y-%m-%d')
            dailyHistoricalWeatherPD = pd.concat([dailyHistoricalWeatherPD, data], ignore_index=True, sort=True)
    #Screen these dates
    dailyHistoricalWeatherPD = dailyHistoricalWeatherPD[(dailyHistoricalWeatherPD.Date >=startDtStr) & (dailyHistoricalWeatherPD.Date <= endDtStr)]
    dailyHistoricalWeatherPD.reset_index(drop=True, inplace=True)


    #Sort by Julian Date
    dailyHistoricalWeatherPD['JulianDate'] = dailyHistoricalWeatherPD.Date.dt.dayofyear
    weather_df =  dailyHistoricalWeatherPD.iloc[dailyHistoricalWeatherPD['JulianDate'].astype(int).argsort()]
    weather_df.reset_index(inplace=True,drop=True)
    dmw_pd = weather_df.copy()

    
    #Put values for metric system calculations
    dmw_pd['Tavg']=dmw_pd['T_mean'].values.astype('float') #convert to DegC
    dmw_pd['Rnet']=dmw_pd['R_n__MJpm2'].values.astype('float')
    dmw_pd['Wavg']=dmw_pd['Uz'].values.astype('float')

    #Screen out the null values
    dmw_pd['ReferenceETo(in)'] = dmw_pd['ETo_average_inches'].replace('--',np.nan, regex=True).astype('float')
    dmw_pd['AverageAirTemperature(F)'] = dmw_pd['T_mean'].replace('--',np.nan, regex=True).astype('float')
    dmw_pd['SolarRadiationAverage(Ly/day)'] = dmw_pd['R_n__MJpm2'].replace('--',np.nan, regex=True).astype('float') * 23.900574 #convert to Mangley /day
    dmw_pd['AverageWindSpeed(mph)'] = dmw_pd['Uz'].replace('--',np.nan, regex=True).astype('float') * 2.24 #convert to Mph

    dmw_pd['StationId'] = dmw_pd['STN']
    dmw_pd['Precipitation(in)'] = dmw_pd['Precip_in'].replace('--',np.nan, regex=True).astype('float')
    dmw_pd['AverageAirTemperature(F)'] = dmw_pd['T_mean'].replace('--',np.nan, regex=True).astype('float') * (9/5) + 32
    dmw_pd['MinimumAirTemperature(F)'] = dmw_pd['T_min'].replace('--',np.nan, regex=True).astype('float')* (9/5) + 32
    dmw_pd['MaximumAirTemperature(F)'] = dmw_pd['T_max'].replace('--',np.nan, regex=True).astype('float')* (9/5) + 32
    dmw_pd['AverageRelativeHumidity(%)'] = dmw_pd['RH_mean'].replace('--',np.nan, regex=True).astype('float')
    dmw_pd['DewPoint(F)'] = dmw_pd['T_dew'].replace('--',np.nan, regex=True).astype('float')* (9/5) + 32


    ######## If one of the Etc's is null, then recompute the average for that
    dmw_pd['ETo_average_inches'] = dmw_pd[['ETo_HAR_inches','ETo_FAO_inches']].mean(axis=1)
    ########
    keep_eto_cols_list = ['Date', 
                          'StationId', 
                          'ReferenceETo(in)',
                          'Precipitation(in)', 
                          'MaximumAirTemperature(F)',
                          'MinimumAirTemperature(F)',
                          'AverageAirTemperature(F)',
                          'AverageRelativeHumidity(%)',
                          'DewPoint(F)',
                          'ETo_average_inches',
                          'Precip_in','R_n__MJpm2','P','T_mean','T_min','T_max','RH_mean','RH_max','RH_min','T_dew','Uz'
                         ]
    eto = dmw_pd[keep_eto_cols_list]

    #print(eto)

    keep_metric_cols_list = ['Date','StationId','Tavg','Rnet','Wavg']
    metric = dmw_pd[keep_metric_cols_list] 

    #for each station, we have to interpolate the timeseries and fill in the values
    #first for ETo
    new_list_of_active_stns = []
    list_of_stns = pd.unique(eto.StationId)
    eto_interp = pd.DataFrame()
    #print(list_of_stns)
    for stn in list_of_stns:
        eto_stn = eto[eto.StationId==stn]
        eto_stn = eto_stn.set_index(['Date'])
        eto_stn = eto_stn.sort_index()
        n = len(eto_stn[eto_stn['ReferenceETo(in)'].isnull()==True])
        if n==0:
            eto_interp = eto_interp.append(eto_stn)
        else:
            if (n>=100):
                print([n,stn])
                #dropping these
            else:
                print('\tInterpolating for '+str(n)+' points for stn '+ str(stn))
                new_list_of_active_stns.append(stn)
                eto_stn_interpolate = eto_stn.interpolate(method='time')
                eto_interp = eto_interp.append(eto_stn_interpolate)

    stns_df = stns_gdf.copy()
    stns_df = stns_df.rename(columns={"stationIdentifier": "StationId"})
    eto_interp = eto_interp.reset_index()

    #print(eto_interp.columns)
    #print(stns_df.columns)

    df = pd.merge(eto_interp,stns_df,on='StationId')
    df_metric = pd.merge(df,metric,how='left', left_on=['StationId','Date'], right_on = ['StationId','Date'])
    #convert into a GeoDataFrame
    crs = {'init': 'epsg:4326'}
    gdf_metric = gpd.GeoDataFrame(df_metric, crs=crs, geometry=df_metric.geometry)

    #print(gdf_metric)

    historicWeather_gdf =  gdf_metric.copy()
    historicWeather_gdf.reset_index(drop=True, inplace=True)
    historicWeather_gdf['Wavg']=historicWeather_gdf['Wavg'].values.astype('float')
    historicWeather_gdf['Rnet']=historicWeather_gdf['Rnet'].values.astype('float')
    historicWeather_gdf['Uz']=historicWeather_gdf['Uz'].values.astype('float')
    historicWeather_gdf['Precip_in']=historicWeather_gdf['Precip_in'].values.astype('float')

    #print(historicWeather_gdf.index)

    #Set the horizon from the inputs
    deltaT = datetime.strptime(startDtStr,'%Y-%m-%d') - datetime.strptime(endDtStr,'%Y-%m-%d')
    list_of_dates = pd.date_range(pd.datetime.strptime(startDtStr,'%Y-%m-%d'), periods=-deltaT.days+1).to_pydatetime().tolist()

    #print(acct_gdf.columns)

    startt = time.time()
    EToForAccounts = pd.DataFrame()
    for i,row in acct_gdf.iterrows():
        print('\t\tCalculating ETo for'+
            #' Customerid: '+str(row.customerid)+
            #', Assetid: '+str(row.fieldassetid)+
            ', Acct: ' + row.Category+
            ', Ranch: '+row.SubCategory+
            ', Field: '+row.Name)
        start_time = time.time()
        p = loads(row.CentroidPoint) #shapely Point obj
        
        etv = getWeatherForDates_NOAA(p,list_of_dates,historicWeather_gdf)

        #print(etv)

        #now add the columns for the other fields to etv
        etv['UUID']=row.UUID
        #etv['customerid']=row.customerid
        #etv['fieldassetid']=row.fieldassetid
        etv['Category']=row.Category
        etv['SubCategory']=row.SubCategory
        etv['Name']=row.Name

        #print(etv)

        EToForAccounts = EToForAccounts.append(etv, ignore_index = True)
        #print(EToForAccounts)
        print("\t\tEND ETo--- %s seconds ---" % (time.time() - start_time))
    print("\tETo Done --- %s seconds ---" % (time.time() - startt))
    print('\tDone!')
    

    EToForAccounts = EToForAccounts.rename(columns={'ETo_average_inches': 'ETo(in)', 
                                                        'Precip_in': 'Precipitation(in)',
                                                        'T_max': 'MaximumAirTemperature(F)',
                                                        'T_min': 'MinimumAirTemperature(F)',
                                                        'T_mean': 'AverageAirTemperature(F)',
                                                        'RH_max': 'MaximumRelativeHumidity(%)',
                                                        'RH_min': 'MinimumRelativeHumidity(%)',
                                                        'RH_mean': 'AverageRelativeHumidity(%)',
                                                        'T_dew': 'DewPoint(F)'
                                                       })


    #Merge the acct
    #print('EToForAccounts.columns '+EToForAccounts.columns)
    #merged_df = pd.merge(acct_gdf, EToForAccounts, on=['Category','SubCategory','Name','customerid','fieldassetid','UUID'])
    merged_df = pd.merge(acct_gdf, EToForAccounts, on=['Category','SubCategory','Name','UUID'])
    #merged_df = pd.merge(acct_gdf, EToForAccounts, on=['customerid','fieldassetid', 'UUID'])
    merged_df.reset_index() #reset the index to avoid any export issues
    merged_df['Date'] = pd.to_datetime(merged_df['Date'])
    
    #print(merged_df)

    return merged_df

def generateETc_NOAA(merged_df):

    unique_accounts = pd.unique(merged_df.customerid)
    startt = time.time()
    #Now to put this into the output folder for each account
    print("Start --- %s  ---" % datetime.now())

    for acct in unique_accounts:
        a = 'JainlogicCustomerid_'+str(acct)
        #try:
        starta = time.time()
        
        print('\n\tWriting...%s ' % a)
        print("\tStart --- %s  ---" % datetime.now())
        acctName = a.upper()
        acct_dir = acctOutputPath + a + '/data/'
        acct_dir_img = acctOutputPath + a + '/img/'
        acct_dir_pkl = acctOutputPath + a + '/pkl/'
        acct_dir_rep = acctOutputPath + a + '/reports/'
            
        #create the directory if one doesn't exist
        if not os.path.exists(acct_dir):
            os.makedirs(acct_dir)
        if not os.path.exists(acct_dir_img):
            os.makedirs(acct_dir_img)
        if not os.path.exists(acct_dir_pkl):
            os.makedirs(acct_dir_pkl)
        if not os.path.exists(acct_dir_rep):
            os.makedirs(acct_dir_rep)
            
        #For each account, take the subset of the acct_df 
        acctDF = merged_df[merged_df.customerid == acct]
        
        acctGDF = gpd.GeoDataFrame(acctDF,geometry=acctDF.geometry)
        
        #for every day, do the following analysis and printout
        acct_dates_list = acctGDF.Date.tolist()
        unique_dts = pd.unique(acct_dates_list) # so we don't multiple count for each field
        unique_dts.sort()
        
        for dt in unique_dts:
            
            endDay = dt        
            
            startd = time.time()
            #print("\n\t\tSTART --- %s  ---" % datetime.now())

            etoDateStr =  datetime.strftime(dt,'%Y-%m-%d')
            jsonFile = acct_dir + acctName+'__ET_'+etoDateStr+'.json'
            
            acctGDF_date = acctGDF[acctGDF.Date==dt]
            acctGDF_date.reset_index(drop=True, inplace=True)
            
            print('\t\t'+etoDateStr)
            #print(acctGDF_date.columns)
            for i,row in acctGDF_date.iterrows():

                customerid = str(row.customerid)
                acct = row.Category
                ranch = row.SubCategory
                field = row.Name
                fieldassetid = 'ASSETID_'+str(row.fieldassetid)
                
                #ETo and weather values
                #print('row: '+str(row))

                eto = row['ETo(in)']
                tavg = row['Tavg']
                rnet = row['Rnet']
                wavg = row['Wavg']
                ############## add others Raw and Deried
                precip = row['Precipitation(in)']
                maxAirTemp = row['MaximumAirTemperature(F)']
                minAirTemp = row['MinimumAirTemperature(F)']
                avgAirTemp = row['AverageAirTemperature(F)']
                maxRelHumPercent = row['MaximumRelativeHumidity(%)']
                minRelHumPercent = row['MinimumRelativeHumidity(%)']
                avgRelHumPercent = row['AverageRelativeHumidity(%)']
                dewPoint = row['DewPoint(F)']
                
                #########################################
            
                # NDVI value
                P = getSinglePoly(row.geometry)

                
                ss = row.SentinelFootprint
                ll = row.LandsatFootprint
                
                #do Sentinel
                if not type(ss)==list:
                    if (ss.find(',')==-1): #There is only ONE tile
                        S = [ss.replace('[','').replace(']','').replace('\'','')]
                    else:
                        tok = ss.split(',')
                        l = len(tok)
                        S=[]
                        for ii in range(l):
                            tk = tok[ii].replace('[','').replace(']','').replace('\'','')
                            S.append(tk)
                else:
                    S = ss
                if (type(S)==str):
                    S = [S]
            
                #do landsat
                if not type(ll)==list:
                    if (ss.find(',')==-1): #There is only ONE tile
                        L = ll.replace('[','').replace(']','').replace('\'','')
                    else:
                        tok = ll.split(', ')
                        l = len(tok)
                        L=[]
                        for ii in range(l):
                            tk = tok[ii].replace('[','').replace(']','').replace('\'','')
                            L.append(tk.replace(' ',''))
                else:
                    L = [str(x).zfill(6) for x in ll]

                if (type(L)==str):
                    L = [L]
                

                #print([S,L,etoDateStr])

                ndvi_mList_df = pd.DataFrame()
                
                NDVIFiles_df = getNDVIFilesWithDate(32, S, L, etoDateStr)
                NDVIFiles_df.reset_index(drop=True,inplace=True)

                #print(NDVIFiles_df)
                
                [ndvi_df,ndvi_img,ndvi_ext,imgDtStr] = getNDVIValueWithStatsAndImg3(NDVIFiles_df,P)
                
                #print([ndvi_df,ndvi_ext,imgDtStr])
                
                if (len(ndvi_df)>0 & (len(ndvi_img)>0)): #The first time the df has valid values
                    ext = json.loads(ndvi_ext)
                    ndvi_mList_df = ndvi_df.copy()
                    
                    
                    acct_dir_img2 = acct_dir_img + fieldassetid + '/'
                    if not os.path.exists(acct_dir_img2):
                        os.makedirs(acct_dir_img2)
                    imgFile = acct_dir_img2 + acctName+'_'+fieldassetid+'__IMG_'+etoDateStr+'.png'
                    infoJsonFile = acct_dir_img2 + acctName+'_'+fieldassetid+'__INFO_'+etoDateStr+'.json'

                    acct_dir_pkl2 = acct_dir_pkl + fieldassetid + '/'
                    if not os.path.exists(acct_dir_pkl2):
                        os.makedirs(acct_dir_pkl2)
                    pklFile = acct_dir_pkl2 + acctName+'_'+fieldassetid+'__PKL_'+etoDateStr+'.pkl'
                    
                
                    imgDict = {
                        'ImgFileName': os.path.split(imgFile)[1],
                        'ImgDate': imgDtStr,
                        'ImgType': 'ETc',
                        'ImgExtent': ndvi_ext, 
                        'ImgBoundary': P.to_wkt(),
                        'AssetID': row.fieldassetid,
                        'CustomerID':a.split('_')[1],
                        'Ranch': ranch,
                        'Field': field
                    }
                    with open(infoJsonFile, 'w', encoding='utf8') as outfile:
                        json.dump(imgDict, outfile, ensure_ascii=False)

                    plt.ioff();
                    with open(imgFile, 'wb') as outfile:
                        fig = plt.figure()
                        plt.imshow(ndvi_img, extent=ext, cmap=zoom1, vmin=0, vmax=1)
                        fig.gca().set_axis_off()
                        plt.savefig(outfile, dpi=600, bbox_inches='tight', pad_inches=0, transparent=True)
                        plt.close();
                        
                    with open(pklFile,'wb') as f:
                        pickle.dump(ndvi_img, f)
                
                    retdf = ndvi_df.loc[0]
                    ndvi_mean = retdf.Mean
                    ndvi_p90 = retdf.p90
                    ndvi_p25 = retdf.p25
                    ndvi_p75 = retdf.p75
                    ndvi_median = retdf.Median
                    ndvi_std = retdf.StdDev
                    
                else:
                    #set values if NaN
                    ndvi_mean = np.nan
                    ndvi_p90 = np.nan
                    ndvi_p25 = np.nan
                    ndvi_p75 = np.nan
                    ndvi_median = np.nan
                    ndvi_std = np.nan
                    
                #print(['NDVI_mean = '+str(ndvi_mean)])
                #derived
                ndvi_unif = (ndvi_p25 / ndvi_median)*100
            
                #compute the METRIC ETc
                #print('Metric inputs:')
                #print([tavg,rnet,wavg,ndvi_mean])
                et_metric = computeMetricET(tavg,rnet,wavg,ndvi_mean)
                
                ##########if et_metric is NaN, make it zero
                if (np.isnan(et_metric)):
                    et_metric=0.0

                #print('ETc metric='+str(et_metric))
                #the Kc
                fac = 1.35 #somewhat arbitrary for now
                k_c = 1.3558*ndvi_mean*fac + 0.05
                k_c2 = 1.3558*ndvi_mean + 0.05
            
                #compute the ETc
                et_ndvi = eto*k_c
                et_ndvi2 = eto*k_c2
            
                #get the final ETc
                # Include et_metric
                k_a = 1 #We can reduce k_a to increase the effect of Metric (see above)
                et_c =  (1.15 - k_a)*et_metric + ((k_a + 0.0)*et_ndvi)
                et_c2 =  (1.15 - k_a)*et_metric + ((k_a + 0.0)*et_ndvi2)
                
                k_c_p90 = ((1.3558*ndvi_p90*fac + 0.05) + (1.3558*ndvi_p90 + 0.05))/2
                et_p90 = eto*k_c_p90
                
                #print([eto,k_c,k_c2,et_metric,et_c, et_c2, et_p90,ndvi_unif,precip])

                acctGDF_date.loc[i,'ETc']=et_c
                acctGDF_date.loc[i,'ETc_low']=et_c2
                acctGDF_date.loc[i,'Kc']=k_c
                acctGDF_date.loc[i,'ETm']=et_metric
                acctGDF_date.loc[i,'ET_p90']=et_p90
                acctGDF_date.loc[i,'ET_unif']=ndvi_unif
                acctGDF_date.loc[i,'Precipitation(in)']=precip
                acctGDF_date.loc[i,'MaximumAirTemperature(F)']=maxAirTemp
                acctGDF_date.loc[i,'MinimumAirTemperature(F)']=minAirTemp
                acctGDF_date.loc[i,'AverageAirTemperature(F)']=avgAirTemp
                acctGDF_date.loc[i,'MaximumRelativeHumidity(%)']=maxRelHumPercent
                acctGDF_date.loc[i,'MinimumRelativeHumidity(%)']=minRelHumPercent
                acctGDF_date.loc[i,'AverageRelativeHumidity(%)']=avgRelHumPercent
                acctGDF_date.loc[i,'DewPoint(F)']=dewPoint
                

            acctGDF_date['SentinelFootprint'] = acctGDF_date['SentinelFootprint'].astype('str')
            acctGDF_date['LandsatFootprint'] = acctGDF_date['LandsatFootprint'].astype('str')
            acctGDF_date['Date'] = acctGDF_date['Date'].astype('str')
            acctGDF_date.drop('CentroidPoint',inplace=True, axis=1)
            print(acctGDF_date.ETc.tolist())
            
            #save the values
            acctGDF_date.to_file(jsonFile, driver="GeoJSON")

        print("\tEND --- %s seconds ---" % (time.time() - starta))
        print('\tDone!')
        #except:
        #	print('Exception in Acct: '+str(a))
        #	continue

    print("\nEND --- %s seconds ---" % (time.time() - startt))
    return

def generateETc_NOAA2(merged_df):

    unique_accounts = pd.unique(merged_df.Category)
    startt = time.time()
    #Now to put this into the output folder for each account
    print("Start --- %s  ---" % datetime.now())

    for acct in unique_accounts:
        a = str(acct)
        try:
            starta = time.time()
            
            print('\n\tWriting...%s ' % a)
            print("\tStart --- %s  ---" % datetime.now())
            acctName = a.upper()
            acct_dir = agls_acctOutputPath + a + '/data/'
            acct_dir_img = agls_acctOutputPath + a + '/images/'
            acct_dir_pkl = agls_acctOutputPath + a + '/pkl/'
            acct_dir_rep = agls_acctOutputPath + a + '/reports/'
                
            #create the directory if one doesn't exist
            if not os.path.exists(acct_dir):
                os.makedirs(acct_dir)
            if not os.path.exists(acct_dir_img):
                os.makedirs(acct_dir_img)
            if not os.path.exists(acct_dir_pkl):
                os.makedirs(acct_dir_pkl)
            if not os.path.exists(acct_dir_rep):
                os.makedirs(acct_dir_rep)
                
            #For each account, take the subset of the acct_df 
            acctDF = merged_df[merged_df.Category == acct]
            
            acctGDF = gpd.GeoDataFrame(acctDF,geometry=acctDF.geometry)

            #print(type(acctGDF))
            
            #for every day, do the following analysis and printout
            acct_dates_list = acctGDF.Date.tolist()
            unique_dts = pd.unique(acct_dates_list) # so we don't multiple count for each field
            unique_dts.sort()
            
            for dt in unique_dts:
                
                endDay = dt        
                
                startd = time.time()
                #print("\n\t\tSTART --- %s  ---" % datetime.now())

                etoDateStr =  datetime.strftime(dt,'%Y-%m-%d')
                jsonFile = acct_dir + acctName+'_ET_'+etoDateStr+'.json'
                
                acctGDF_date = acctGDF[acctGDF.Date==dt]
                acctGDF_date.reset_index(drop=True, inplace=True)

                #Convert back to GPD
                #acctGDF_date = gpd.GeoDataFrame(acctGDF_date, geometry=acctDF.geometry)
                
                print('\t\t'+etoDateStr)
                #print(acctGDF_date.columns)
                for i,row in acctGDF_date.iterrows():

                    acct = row.Category
                    ranch = row.SubCategory
                    field = row.Name
                    fieldName = field.upper()
                    
                    #ETo and weather values
                    #print('row: '+str(row))

                    eto = row['ETo(in)']
                    tavg = row['Tavg']
                    rnet = row['Rnet']
                    wavg = row['Wavg']
                    ############## add others Raw and Deried
                    precip = row['Precipitation(in)']
                    maxAirTemp = row['MaximumAirTemperature(F)']
                    minAirTemp = row['MinimumAirTemperature(F)']
                    avgAirTemp = row['AverageAirTemperature(F)']
                    maxRelHumPercent = row['MaximumRelativeHumidity(%)']
                    minRelHumPercent = row['MinimumRelativeHumidity(%)']
                    avgRelHumPercent = row['AverageRelativeHumidity(%)']
                    dewPoint = row['DewPoint(F)']
                    
                    #########################################
                
                    # NDVI value
                    #print(row.geometry)
                    #print(type(row.geometry))
                    P = getSinglePoly(row.geometry)

                    
                    ss = row.SentinelFootprint
                    ll = row.LandsatFootprint
                    
                    #do Sentinel
                    if not type(ss)==list:
                        if (ss.find(',')==-1): #There is only ONE tile
                            S = [ss.replace('[','').replace(']','').replace('\'','')]
                        else:
                            tok = ss.split(',')
                            l = len(tok)
                            S=[]
                            for ii in range(l):
                                tk = tok[ii].replace('[','').replace(']','').replace('\'','')
                                S.append(tk)
                    else:
                        S = ss
                    if (type(S)==str):
                        S = [S]
                
                    #do landsat
                    if not type(ll)==list:
                        if (ss.find(',')==-1): #There is only ONE tile
                            L = ll.replace('[','').replace(']','').replace('\'','')
                        else:
                            tok = ll.split(', ')
                            l = len(tok)
                            L=[]
                            for ii in range(l):
                                tk = tok[ii].replace('[','').replace(']','').replace('\'','')
                                L.append(tk.replace(' ',''))
                    else:
                        L = [str(x).zfill(6) for x in ll]

                    if (type(L)==str):
                        L = [L]
                    

                    #print([S,L,etoDateStr])

                    ndvi_mList_df = pd.DataFrame()
                    
                    NDVIFiles_df = getNDVIFilesWithDate(32, S, L, etoDateStr)
                    NDVIFiles_df.reset_index(drop=True,inplace=True)

                    #print(NDVIFiles_df)
                    
                    [ndvi_df,ndvi_img,ndvi_ext,imgDtStr] = getNDVIValueWithStatsAndImg3(NDVIFiles_df,P)
                    
                    #print([ndvi_df,ndvi_ext,imgDtStr])
                    
                    if (len(ndvi_df)>0 & (len(ndvi_img)>0)): #The first time the df has valid values
                        ext = json.loads(ndvi_ext)
                        ndvi_mList_df = ndvi_df.copy()
                        
                        
                        acct_dir_img2 = acct_dir_img 
                        if not os.path.exists(acct_dir_img2):
                            os.makedirs(acct_dir_img2)
                        imgFile = acct_dir_img2 + acctName+'_'+field+'__IMG_'+etoDateStr+'.png'
                        infoJsonFile = acct_dir_img2 + acctName+'_'+field+'__INFO_'+etoDateStr+'.json'

                        acct_dir_pkl2 = acct_dir_pkl 
                        if not os.path.exists(acct_dir_pkl2):
                            os.makedirs(acct_dir_pkl2)
                        pklFile = acct_dir_pkl2 + acctName+'_'+field+'__PKL_'+etoDateStr+'.pkl'
                        
                    
                        imgDict = {
                            'ImgFileName': os.path.split(imgFile)[1],
                            'ImgDate': imgDtStr,
                            'ImgType': 'ETc',
                            'ImgExtent': ndvi_ext, 
                            'ImgBoundary': P.to_wkt(),
                            'Ranch': ranch,
                            'Field': field
                        }
                        with open(infoJsonFile, 'w', encoding='utf8') as outfile:
                            json.dump(imgDict, outfile, ensure_ascii=False)

                        plt.ioff();
                        with open(imgFile, 'wb') as outfile:
                            fig = plt.figure()
                            plt.imshow(ndvi_img, extent=ext, cmap=zoom1, vmin=0, vmax=1)
                            fig.gca().set_axis_off()
                            plt.savefig(outfile, dpi=600, bbox_inches='tight', pad_inches=0, transparent=True)
                            plt.close();
                            
                        with open(pklFile,'wb') as f:
                            pickle.dump(ndvi_img, f)
                    
                        retdf = ndvi_df.loc[0]
                        ndvi_mean = retdf.Mean
                        ndvi_p90 = retdf.p90
                        ndvi_p25 = retdf.p25
                        ndvi_p75 = retdf.p75
                        ndvi_median = retdf.Median
                        ndvi_std = retdf.StdDev
                        
                    else:
                        #set values if NaN
                        ndvi_mean = np.nan
                        ndvi_p90 = np.nan
                        ndvi_p25 = np.nan
                        ndvi_p75 = np.nan
                        ndvi_median = np.nan
                        ndvi_std = np.nan
                        
                    #print(['NDVI_mean = '+str(ndvi_mean)])
                    #derived
                    ndvi_unif = (ndvi_p25 / ndvi_median)*100
                
                    #compute the METRIC ETc
                    #print('Metric inputs:')
                    #print([tavg,rnet,wavg,ndvi_mean])
                    et_metric = computeMetricET(tavg,rnet,wavg,ndvi_mean)
                    
                    ##########if et_metric is NaN, make it zero
                    if (np.isnan(et_metric)):
                        et_metric=0.0

                    #print('ETc metric='+str(et_metric))
                    #the Kc
                    fac = 1.35 #somewhat arbitrary for now
                    k_c = 1.3558*ndvi_mean*fac + 0.05
                    k_c2 = 1.3558*ndvi_mean + 0.05
                
                    #compute the ETc
                    et_ndvi = eto*k_c
                    et_ndvi2 = eto*k_c2
                
                    #get the final ETc
                    # Include et_metric
                    k_a = 1 #We can reduce k_a to increase the effect of Metric (see above)
                    et_c =  (1.15 - k_a)*et_metric + ((k_a + 0.0)*et_ndvi)
                    et_c2 =  (1.15 - k_a)*et_metric + ((k_a + 0.0)*et_ndvi2)
                    
                    k_c_p90 = ((1.3558*ndvi_p90*fac + 0.05) + (1.3558*ndvi_p90 + 0.05))/2
                    et_p90 = eto*k_c_p90
                    
                    #print([eto,k_c,k_c2,et_metric,et_c, et_c2, et_p90,ndvi_unif,precip])

                    acctGDF_date.loc[i,'ETc']=et_c
                    acctGDF_date.loc[i,'ETc_low']=et_c2
                    acctGDF_date.loc[i,'Kc']=k_c
                    acctGDF_date.loc[i,'ETm']=et_metric
                    acctGDF_date.loc[i,'ET_p90']=et_p90
                    acctGDF_date.loc[i,'ET_unif']=ndvi_unif
                    acctGDF_date.loc[i,'Precipitation(in)']=precip
                    acctGDF_date.loc[i,'MaximumAirTemperature(F)']=maxAirTemp
                    acctGDF_date.loc[i,'MinimumAirTemperature(F)']=minAirTemp
                    acctGDF_date.loc[i,'AverageAirTemperature(F)']=avgAirTemp
                    acctGDF_date.loc[i,'MaximumRelativeHumidity(%)']=maxRelHumPercent
                    acctGDF_date.loc[i,'MinimumRelativeHumidity(%)']=minRelHumPercent
                    acctGDF_date.loc[i,'AverageRelativeHumidity(%)']=avgRelHumPercent
                    acctGDF_date.loc[i,'DewPoint(F)']=dewPoint
                    

                acctGDF_date['SentinelFootprint'] = acctGDF_date['SentinelFootprint'].astype('str')
                acctGDF_date['LandsatFootprint'] = acctGDF_date['LandsatFootprint'].astype('str')
                acctGDF_date['to_email_list'] = acctGDF_date['to_email_list'].astype('str')
                acctGDF_date['cc_email_list'] = acctGDF_date['cc_email_list'].astype('str')
                acctGDF_date['bcc_email_list'] = acctGDF_date['bcc_email_list'].astype('str')
                acctGDF_date['Date'] = acctGDF_date['Date'].astype('str')
                acctGDF_date.drop('CentroidPoint',inplace=True, axis=1)

                #print(type(acctGDF_date))
                print(acctGDF_date.ETc.tolist())
                
                #save the values
                #acctGDF_date = acctGDF_date.astype('str')
                acctGDF_date.to_file(jsonFile, driver="GeoJSON")

            print("\tEND --- %s seconds ---" % (time.time() - starta))
            print('\tDone!')
        except:
            print('Exception in Acct: '+str(a))
            continue

    print("\nEND --- %s seconds ---" % (time.time() - startt))
    return

def generateETc_NOAA3(merged_df):

    unique_accounts = pd.unique(merged_df.customerid)
    startt = time.time()
    #Now to put this into the output folder for each account
    print("Start --- %s  ---" % datetime.now())

    for acct in unique_accounts:
        a = 'JainlogicCustomerid_'+str(acct)
        #try:
        starta = time.time()
        
        print('\n\tWriting...%s ' % a)
        print("\tStart --- %s  ---" % datetime.now())
        acctName = a.upper()
        acct_dir = acctOutputPath + a + '/data/'
        acct_dir_img = acctOutputPath + a + '/img/'
        acct_dir_pkl = acctOutputPath + a + '/pkl/'
        acct_dir_rep = acctOutputPath + a + '/reports/'
            
        #create the directory if one doesn't exist
        if not os.path.exists(acct_dir):
            os.makedirs(acct_dir)
        if not os.path.exists(acct_dir_img):
            os.makedirs(acct_dir_img)
        if not os.path.exists(acct_dir_pkl):
            os.makedirs(acct_dir_pkl)
        if not os.path.exists(acct_dir_rep):
            os.makedirs(acct_dir_rep)
            
        #For each account, take the subset of the acct_df 
        acctDF = merged_df[merged_df.customerid == acct]
        
        acctGDF = gpd.GeoDataFrame(acctDF,geometry=acctDF.geometry)
        
        #for every day, do the following analysis and printout
        acct_dates_list = acctGDF.Date.tolist()
        unique_dts = pd.unique(acct_dates_list) # so we don't multiple count for each field
        unique_dts.sort()

        for dt in unique_dts:
            
            endDay = dt        
            
            startd = time.time()
            #print("\n\t\tSTART --- %s  ---" % datetime.now())

            etoDateStr =  datetime.strftime(dt,'%Y-%m-%d')
            jsonFile = acct_dir + acctName+'__ET_'+etoDateStr+'.json'
            
            acctGDF_date = acctGDF[acctGDF.Date==dt]
            acctGDF_date.reset_index(drop=True, inplace=True)
            
            print('\t\t'+etoDateStr)
            #print(acctGDF_date.columns)
            for i,row in acctGDF_date.iterrows():

                customerid = str(row.customerid)
                acct = row.Category
                ranch = row.SubCategory
                field = row.Name
                fieldassetid = 'ASSETID_'+str(row.fieldassetid)

                #ETo and weather values
                #print('row: '+str(row))
                if 'Tavg' in row.index.tolist():
                    eto = row['ETo(in)']
                    tavg = row['Tavg']
                    rnet = row['Rnet']
                    wavg = row['Wavg']
                    ############## add others Raw and Deried
                    precip = row['Precipitation(in)']
                    maxAirTemp = row['MaximumAirTemperature(F)']
                    minAirTemp = row['MinimumAirTemperature(F)']
                    avgAirTemp = row['AverageAirTemperature(F)']
                    maxRelHumPercent = row['MaximumRelativeHumidity(%)']
                    minRelHumPercent = row['MinimumRelativeHumidity(%)']
                    avgRelHumPercent = row['AverageRelativeHumidity(%)']
                    dewPoint = row['DewPoint(F)']
                else:
                    eto = row['ETo(in)']

                #########################################

                P = getSinglePoly(row.geometry)


                ss = row.SentinelFootprint
                ll = row.LandsatFootprint

                #do Sentinel
                if not type(ss)==list:
                    if (ss.find(',')==-1): #There is only ONE tile
                        S = [ss.replace('[','').replace(']','').replace('\'','')]
                    else:
                        tok = ss.split(',')
                        l = len(tok)
                        S=[]
                        for ii in range(l):
                            tk = tok[ii].replace('[','').replace(']','').replace('\'','')
                            S.append(tk)
                else:
                    S = ss
                if (type(S)==str):
                    S = [S]

                #do landsat
                if not type(ll)==list:
                    if (ss.find(',')==-1): #There is only ONE tile
                        L = ll.replace('[','').replace(']','').replace('\'','')
                    else:
                        tok = ll.split(', ')
                        l = len(tok)
                        L=[]
                        for ii in range(l):
                            tk = tok[ii].replace('[','').replace(']','').replace('\'','')
                            L.append(tk.replace(' ',''))
                else:
                    L = [str(x).zfill(6) for x in ll]

                if (type(L)==str):
                    L = [L]


                #print([S,L,etoDateStr])

                ndvi_mList_df = pd.DataFrame()

                NDVIFiles_df = getNDVIFilesWithDate(62, S, L, etoDateStr)
                NDVIFiles_df.reset_index(drop=True,inplace=True)

                #print(NDVIFiles_df)

                [ndvi_df,ndvi_img,ndvi_ext,imgDtStr] = getNDVIValueWithStatsAndImg3(NDVIFiles_df,P)

                #print([ndvi_df,ndvi_ext,imgDtStr])

                if (len(ndvi_df)>0 & (len(ndvi_img)>0)): #The first time the df has valid values
                    ext = json.loads(ndvi_ext)
                    ndvi_mList_df = ndvi_df.copy()
                    
                    
                    acct_dir_img2 = acct_dir_img + fieldassetid + '/'
                    if not os.path.exists(acct_dir_img2):
                        os.makedirs(acct_dir_img2)
                    imgFile = acct_dir_img2 + acctName+'_'+fieldassetid+'__IMG_'+etoDateStr+'.png'
                    infoJsonFile = acct_dir_img2 + acctName+'_'+fieldassetid+'__INFO_'+etoDateStr+'.json'

                    acct_dir_pkl2 = acct_dir_pkl + fieldassetid + '/'
                    if not os.path.exists(acct_dir_pkl2):
                        os.makedirs(acct_dir_pkl2)
                    pklFile = acct_dir_pkl2 + acctName+'_'+fieldassetid+'__PKL_'+etoDateStr+'.pkl'
                    
                
                    imgDict = {
                        'ImgFileName': os.path.split(imgFile)[1],
                        'ImgDate': imgDtStr,
                        'ImgType': 'ETc',
                        'ImgExtent': ndvi_ext, 
                        'ImgBoundary': P.to_wkt(),
                        'AssetID': row.fieldassetid,
                        'CustomerID':a.split('_')[1],
                        'Ranch': ranch,
                        'Field': field
                    }
                    with open(infoJsonFile, 'w', encoding='utf8') as outfile:
                        json.dump(imgDict, outfile, ensure_ascii=False)

                    plt.ioff();
                    with open(imgFile, 'wb') as outfile:
                        fig = plt.figure()
                        plt.imshow(ndvi_img, extent=ext, cmap=zoom1, vmin=0, vmax=1)
                        fig.gca().set_axis_off()
                        plt.savefig(outfile, dpi=600, bbox_inches='tight', pad_inches=0, transparent=True)
                        plt.close();

                    with open(pklFile,'wb') as f:
                        pickle.dump(ndvi_img, f)

                    retdf = ndvi_df.loc[0]
                    ndvi_mean = retdf.Mean
                    ndvi_p90 = retdf.p90
                    ndvi_p25 = retdf.p25
                    ndvi_p75 = retdf.p75
                    ndvi_median = retdf.Median
                    ndvi_std = retdf.StdDev

                else:
                    #set values if NaN
                    ndvi_mean = np.nan
                    ndvi_p90 = np.nan
                    ndvi_p25 = np.nan
                    ndvi_p75 = np.nan
                    ndvi_median = np.nan
                    ndvi_std = np.nan

                #print(['NDVI_mean = '+str(ndvi_mean)])
                #derived
                ndvi_unif = (ndvi_p25 / ndvi_median)*100

                #compute the METRIC ETc
                if 'Tavg' in row.index.tolist():
                    et_metric = computeMetricET(tavg,rnet,wavg,ndvi_mean)
                else:
                    et_metric=0.0

                ##########if et_metric is NaN, make it zero
                if (np.isnan(et_metric)):
                    et_metric=0.0

                #print('ETc metric='+str(et_metric))
                #the Kc
                fac = 1.35 #somewhat arbitrary for now
                k_c = 1.3558*ndvi_mean*fac + 0.05
                k_c2 = 1.3558*ndvi_mean + 0.05

                #compute the ETc
                et_ndvi = eto*k_c
                et_ndvi2 = eto*k_c2

                #get the final ETc
                # Include et_metric
                k_a = 1 #We can reduce k_a to increase the effect of Metric (see above)
                et_c =  (1.15 - k_a)*et_metric + ((k_a + 0.0)*et_ndvi)
                et_c2 =  (1.15 - k_a)*et_metric + ((k_a + 0.0)*et_ndvi2)

                k_c_p90 = ((1.3558*ndvi_p90*fac + 0.05) + (1.3558*ndvi_p90 + 0.05))/2
                et_p90 = eto*k_c_p90

                #print([eto,k_c,k_c2,et_metric,et_c, et_c2, et_p90,ndvi_unif,precip])

                acctGDF_date.loc[i,'ETc']=et_c
                acctGDF_date.loc[i,'ETc_low']=et_c2
                acctGDF_date.loc[i,'Kc']=k_c
                acctGDF_date.loc[i,'ETm']=et_metric
                acctGDF_date.loc[i,'ET_p90']=et_p90
                acctGDF_date.loc[i,'ET_p25']=ndvi_p25
                acctGDF_date.loc[i,'ET_p75']=ndvi_p75
                acctGDF_date.loc[i,'ET_median']=ndvi_median
                acctGDF_date.loc[i,'ET_StdDev']=ndvi_std
                acctGDF_date.loc[i,'ET_unif']=ndvi_unif
                if 'Tavg' in row.index.tolist():
                    acctGDF_date.loc[i,'Precipitation(in)']=precip
                    acctGDF_date.loc[i,'MaximumAirTemperature(F)']=maxAirTemp
                    acctGDF_date.loc[i,'MinimumAirTemperature(F)']=minAirTemp
                    acctGDF_date.loc[i,'AverageAirTemperature(F)']=avgAirTemp
                    acctGDF_date.loc[i,'MaximumRelativeHumidity(%)']=maxRelHumPercent
                    acctGDF_date.loc[i,'MinimumRelativeHumidity(%)']=minRelHumPercent
                    acctGDF_date.loc[i,'AverageRelativeHumidity(%)']=avgRelHumPercent
                    acctGDF_date.loc[i,'DewPoint(F)']=dewPoint


            acctGDF_date['SentinelFootprint'] = acctGDF_date['SentinelFootprint'].astype('str')
            acctGDF_date['LandsatFootprint'] = acctGDF_date['LandsatFootprint'].astype('str')
            acctGDF_date['Date'] = acctGDF_date['Date'].astype('str')
            acctGDF_date.drop('CentroidPoint',inplace=True, axis=1)

            #print(type(acctGDF_date))
            print(acctGDF_date.ETc.tolist())

            #save the values
            #acctGDF_date = acctGDF_date.astype('str')
            #print(jsonFile)
            acctGDF_date.to_file(jsonFile, driver="GeoJSON")

        print("\tEND --- %s seconds ---" % (time.time() - starta))
        print('\tDone!')
        #except:
        #	print('Exception in Acct: '+str(a))
        #	continue

    print("\nEND --- %s seconds ---" % (time.time() - startt))
    return


def generate_AUS_ETc(merged_df):

    unique_accounts = pd.unique(merged_df.customerid)
    startt = time.time()
    #Now to put this into the output folder for each account
    print("Start --- %s  ---" % datetime.now())

    for acct in unique_accounts:
        a = 'JainlogicCustomerid_'+str(acct)
        #try:
        starta = time.time()

        print('\n\tWriting...%s ' % a)
        print("\tStart --- %s  ---" % datetime.now())
        acctName = a.upper()
        acct_dir = acctOutputPath + a + '/data/'
        acct_dir_img = acctOutputPath + a + '/img/'
        acct_dir_pkl = acctOutputPath + a + '/pkl/'
        acct_dir_rep = acctOutputPath + a + '/reports/'

        #create the directory if one doesn't exist
        if not os.path.exists(acct_dir):
            os.makedirs(acct_dir)
        if not os.path.exists(acct_dir_img):
            os.makedirs(acct_dir_img)
        if not os.path.exists(acct_dir_pkl):
            os.makedirs(acct_dir_pkl)
        if not os.path.exists(acct_dir_rep):
            os.makedirs(acct_dir_rep)

        #For each account, take the subset of the acct_df 
        acctDF = merged_df[merged_df.customerid == acct]

        acctGDF = gpd.GeoDataFrame(acctDF,geometry=acctDF.geometry)

        #for every day, do the following analysis and printout
        acct_dates_list = acctGDF.Date.tolist()
        unique_dts = pd.unique(acct_dates_list) # so we don't multiple count for each field
        unique_dts.sort()

        #print(unique_dts)
        #print(acctGDF.columns)
        #print(acctGDF.head())

        for dt in unique_dts:

            endDay = dt        

            startd = time.time()
            #print("\n\t\tSTART --- %s  ---" % datetime.now())

            etoDateStr =  datetime.strftime(dt,'%Y-%m-%d')
            jsonFile = acct_dir + acctName+'__ET_'+etoDateStr+'.json'

            acctGDF_date = acctGDF[acctGDF.Date==dt]
            acctGDF_date.reset_index(drop=True, inplace=True)

            print('\t\t'+etoDateStr)
            #print(acctGDF_date.columns)
            for i,row in acctGDF_date.iterrows():

                customerid = str(row.customerid)
                acct = row.Category
                ranch = row.SubCategory
                field = row.Name
                fieldassetid = 'ASSETID_'+str(row.fieldassetid)

                #ETo and weather values
                #print('row: '+str(row))
                if 'Tavg' in row.index.tolist():
                    eto = row['ETo(in)']
                    tavg = row['Tavg']
                    rnet = row['Rnet']
                    wavg = row['Wavg']
                    ############## add others Raw and Deried
                    precip = row['Precipitation(in)']
                    maxAirTemp = row['MaximumAirTemperature(F)']
                    minAirTemp = row['MinimumAirTemperature(F)']
                    avgAirTemp = row['AverageAirTemperature(F)']
                    maxRelHumPercent = row['MaximumRelativeHumidity(%)']
                    minRelHumPercent = row['MinimumRelativeHumidity(%)']
                    avgRelHumPercent = row['AverageRelativeHumidity(%)']
                    dewPoint = row['DewPoint(F)']
                else:
                    eto = row['ETo(in)']

                #########################################

                P = getSinglePoly(row.geometry)


                ss = row.SentinelFootprint
                ll = row.LandsatFootprint

                #do Sentinel
                if not type(ss)==list:
                    if (ss.find(',')==-1): #There is only ONE tile
                        S = [ss.replace('[','').replace(']','').replace('\'','')]
                    else:
                        tok = ss.split(',')
                        l = len(tok)
                        S=[]
                        for ii in range(l):
                            tk = tok[ii].replace('[','').replace(']','').replace('\'','')
                            S.append(tk)
                else:
                    S = ss
                if (type(S)==str):
                    S = [S]

                #do landsat
                if not type(ll)==list:
                    if (ss.find(',')==-1): #There is only ONE tile
                        L = ll.replace('[','').replace(']','').replace('\'','')
                    else:
                        tok = ll.split(', ')
                        l = len(tok)
                        L=[]
                        for ii in range(l):
                            tk = tok[ii].replace('[','').replace(']','').replace('\'','')
                            L.append(tk.replace(' ',''))
                else:
                    L = [str(x).zfill(6) for x in ll]

                if (type(L)==str):
                    L = [L]


                #print([S,L,etoDateStr])

                ndvi_mList_df = pd.DataFrame()

                NDVIFiles_df = getNDVIFilesWithDate(62, S, L, etoDateStr)
                NDVIFiles_df.reset_index(drop=True,inplace=True)

                #print(NDVIFiles_df)

                [ndvi_df,ndvi_img,ndvi_ext,imgDtStr] = getNDVIValueWithStatsAndImg3(NDVIFiles_df,P)
                
                #print(ndvi_df)

                #print([ndvi_df,ndvi_ext,imgDtStr])

                if (len(ndvi_df)>0 & (len(ndvi_img)>0)): #The first time the df has valid values
                    ext = json.loads(ndvi_ext)
                    ndvi_mList_df = ndvi_df.copy()


                    acct_dir_img2 = acct_dir_img + fieldassetid + '/'
                    if not os.path.exists(acct_dir_img2):
                        os.makedirs(acct_dir_img2)
                    imgFile = acct_dir_img2 + acctName+'_'+fieldassetid+'__IMG_'+etoDateStr+'.png'
                    infoJsonFile = acct_dir_img2 + acctName+'_'+fieldassetid+'__INFO_'+etoDateStr+'.json'

                    acct_dir_pkl2 = acct_dir_pkl + fieldassetid + '/'
                    if not os.path.exists(acct_dir_pkl2):
                        os.makedirs(acct_dir_pkl2)
                    pklFile = acct_dir_pkl2 + acctName+'_'+fieldassetid+'__PKL_'+etoDateStr+'.pkl'


                    imgDict = {
                        'ImgFileName': os.path.split(imgFile)[1],
                        'ImgDate': imgDtStr,
                        'ImgType': 'ETc',
                        'ImgExtent': ndvi_ext, 
                        'ImgBoundary': P.to_wkt(),
                        'AssetID': row.fieldassetid,
                        'CustomerID':a.split('_')[1],
                        'Ranch': ranch,
                        'Field': field
                    }
                    with open(infoJsonFile, 'w', encoding='utf8') as outfile:
                        json.dump(imgDict, outfile, ensure_ascii=False)

                    plt.ioff();
                    with open(imgFile, 'wb') as outfile:
                        fig = plt.figure()
                        plt.imshow(ndvi_img, extent=ext, cmap=zoom1, vmin=0, vmax=1)
                        fig.gca().set_axis_off()
                        plt.savefig(outfile, dpi=600, bbox_inches='tight', pad_inches=0, transparent=True)
                        plt.close();

                    with open(pklFile,'wb') as f:
                        pickle.dump(ndvi_img, f)

                    retdf = ndvi_df.loc[0]
                    ndvi_mean = retdf.Mean
                    ndvi_p90 = retdf.p90
                    ndvi_p25 = retdf.p25
                    ndvi_p75 = retdf.p75
                    ndvi_median = retdf.Median
                    ndvi_std = retdf.StdDev

                else:
                    #set values if NaN
                    ndvi_mean = np.nan
                    ndvi_p90 = np.nan
                    ndvi_p25 = np.nan
                    ndvi_p75 = np.nan
                    ndvi_median = np.nan
                    ndvi_std = np.nan

                #print(['NDVI_mean = '+str(ndvi_mean)])
                #derived
                ndvi_unif = (ndvi_p25 / ndvi_median)*100

                #compute the METRIC ETc
                if 'Tavg' in row.index.tolist():
                    et_metric = computeMetricET(tavg,rnet,wavg,ndvi_mean)
                else:
                    et_metric=0.0

                ##########if et_metric is NaN, make it zero
                if (np.isnan(et_metric)):
                    et_metric=0.0

                #print('ETc metric='+str(et_metric))
                #the Kc
                fac = 1.35 #somewhat arbitrary for now
                k_c = 1.3558*ndvi_mean*fac + 0.05
                k_c2 = 1.3558*ndvi_mean + 0.05

                #compute the ETc
                et_ndvi = eto*k_c
                et_ndvi2 = eto*k_c2

                #get the final ETc
                # Include et_metric
                k_a = 1 #We can reduce k_a to increase the effect of Metric (see above)
                et_c =  (1.15 - k_a)*et_metric + ((k_a + 0.0)*et_ndvi)
                et_c2 =  (1.15 - k_a)*et_metric + ((k_a + 0.0)*et_ndvi2)

                k_c_p90 = ((1.3558*ndvi_p90*fac + 0.05) + (1.3558*ndvi_p90 + 0.05))/2
                et_p90 = eto*k_c_p90

                #print([eto,k_c,k_c2,et_metric,et_c, et_c2, et_p90,ndvi_unif,precip])

                acctGDF_date.loc[i,'ETc']=et_c
                acctGDF_date.loc[i,'ETc_low']=et_c2
                acctGDF_date.loc[i,'Kc']=k_c
                acctGDF_date.loc[i,'ETm']=et_metric
                acctGDF_date.loc[i,'ET_p90']=et_p90
                acctGDF_date.loc[i,'ET_p25']=ndvi_p25
                acctGDF_date.loc[i,'ET_p75']=ndvi_p75
                acctGDF_date.loc[i,'ET_median']=ndvi_median
                acctGDF_date.loc[i,'ET_StdDev']=ndvi_std
                acctGDF_date.loc[i,'ET_unif']=ndvi_unif
                if 'Tavg' in row.index.tolist():
                    acctGDF_date.loc[i,'Precipitation(in)']=precip
                    acctGDF_date.loc[i,'MaximumAirTemperature(F)']=maxAirTemp
                    acctGDF_date.loc[i,'MinimumAirTemperature(F)']=minAirTemp
                    acctGDF_date.loc[i,'AverageAirTemperature(F)']=avgAirTemp
                    acctGDF_date.loc[i,'MaximumRelativeHumidity(%)']=maxRelHumPercent
                    acctGDF_date.loc[i,'MinimumRelativeHumidity(%)']=minRelHumPercent
                    acctGDF_date.loc[i,'AverageRelativeHumidity(%)']=avgRelHumPercent
                    acctGDF_date.loc[i,'DewPoint(F)']=dewPoint


            acctGDF_date['SentinelFootprint'] = acctGDF_date['SentinelFootprint'].astype('str')
            acctGDF_date['LandsatFootprint'] = acctGDF_date['LandsatFootprint'].astype('str')
            acctGDF_date['Date'] = acctGDF_date['Date'].astype('str')
            acctGDF_date.drop('CentroidPoint',inplace=True, axis=1)

            #print(type(acctGDF_date))
            print(acctGDF_date.ETc.tolist())

            #save the values
            #acctGDF_date = acctGDF_date.astype('str')
            #print(jsonFile)
            acctGDF_date.to_file(jsonFile, driver="GeoJSON")

        print("\tEND --- %s seconds ---" % (time.time() - starta))
        print('\tDone!')
        #except:
        #	print('Exception in Acct: '+str(a))
        #	continue

    print("\nEND --- %s seconds ---" % (time.time() - startt))
    return
"""
def generate_AUS_ETc(merged_df):

    unique_accounts = pd.unique(merged_df.customerid)
    startt = time.time()
    #Now to put this into the output folder for each account
    print("Start --- %s  ---" % datetime.now())

    for acct in unique_accounts:
        a = 'JainlogicCustomerid_'+str(acct)
        #try:
        starta = time.time()

        print('\n\tWriting...%s ' % a)
        print("\tStart --- %s  ---" % datetime.now())
        acctName = a.replace(' ','').upper()
        acct_dir = acctOutputPath + a + '/data/'
        acct_dir_img = acctOutputPath + a + '/img/'
        acct_dir_pkl = acctOutputPath + a + '/pkl/'
        acct_dir_rep = acctOutputPath + a + '/reports/'

        #create the directory if one doesn't exist
        if not os.path.exists(acct_dir):
            os.makedirs(acct_dir)
        if not os.path.exists(acct_dir_img):
            os.makedirs(acct_dir_img)
        if not os.path.exists(acct_dir_pkl):
            os.makedirs(acct_dir_pkl)
        if not os.path.exists(acct_dir_rep):
            os.makedirs(acct_dir_rep)

        #For each account, take the subset of the acct_df 
        acctDF = merged_df[merged_df.Category == a]

        acctGDF = gpd.GeoDataFrame(acctDF,geometry=acctDF.geometry)

        #for every day, do the following analysis and printout
        acct_dates_list = acctGDF.Date.tolist()
        unique_dts = pd.unique(acct_dates_list) # so we don't multiple count for each field
        unique_dts.sort()


        #print(unique_dts)

        for dt in unique_dts:

            endDay = dt        

            startd = time.time()
            #print("\n\t\tSTART --- %s  ---" % datetime.now())

            etoDateStr =  datetime.strftime(dt,'%Y-%m-%d')
            jsonFile = acct_dir + acctName+'__ET_'+etoDateStr+'.json'

            acctGDF_date = acctGDF[acctGDF.Date==dt]
            acctGDF_date.reset_index(drop=True, inplace=True)

            print('\t\t'+etoDateStr)
            for i,row in acctGDF_date.iterrows():

                customerid = str(row.customerid)
                acct = row.Category
                ranch = row.SubCategory
                field = row.Name
                fieldassetid = 'ASSETID_'+str(row.fieldassetid)

                #ETo and weather values
                eto = row['ETo(in)']
                tavg = row['Tavg']
                #rnet = row['Rnet']
                wavg = row['Wavg']
                ############## add others Raw and Deried
                precip = row['Precip_in']
                maxAirTemp = row['T_max']
                minAirTemp = row['T_min']
                avgAirTemp = row['T_mean']
                maxRelHumPercent = row['RH_max']
                minRelHumPercent = row['RH_min']
                avgRelHumPercent = row['RH_mean']
                dewPoint = row['T_dew']

                #########################################

                # NDVI value
                P = getSinglePoly(row.geometry)


                ss = row.SentinelFootprint
                ll = row.LandsatFootprint

                #do Sentinel
                if not type(ss)==list:
                    if (ss.find(',')==-1): #There is only ONE tile
                        S = [ss.replace('[','').replace(']','').replace('\'','')]
                    else:
                        tok = ss.split(',')
                        l = len(tok)
                        S=[]
                        for ii in range(l):
                            tk = tok[ii].replace('[','').replace(']','').replace('\'','')
                            S.append(tk)
                else:
                    S = ss
                if (type(S)==str):
                    S = [S]

                #do landsat
                if not type(ll)==list:
                    if (ss.find(',')==-1): #There is only ONE tile
                        L = ll.replace('[','').replace(']','').replace('\'','')
                    else:
                        tok = ll.split(', ')
                        l = len(tok)
                        L=[]
                        for ii in range(l):
                            tk = tok[ii].replace('[','').replace(']','').replace('\'','')
                            L.append(tk.replace(' ',''))
                else:
                    L = [str(x).zfill(6) for x in ll]

                if (type(L)==str):
                    L = [L]

                ndvi_mList_df = pd.DataFrame()

                NDVIFiles_df = getNDVIFilesWithDate(32, S, L, etoDateStr)
                NDVIFiles_df.reset_index(drop=True,inplace=True)
                #print(NDVIFiles_df)

                [ndvi_df,ndvi_img,ndvi_ext,imgDtStr] = getNDVIValueWithStatsAndImg3(NDVIFiles_df,P)


                if (len(ndvi_df)>0 & (len(ndvi_img)>0)): #The first time the df has valid values
                    ext = json.loads(ndvi_ext)
                    ndvi_mList_df = ndvi_df.copy()


                    acct_dir_img2 = acct_dir_img + fieldassetid + '/'
                    if not os.path.exists(acct_dir_img2):
                        os.makedirs(acct_dir_img2)
                    imgFile = acct_dir_img2 + acctName+'_'+fieldassetid+'__IMG_'+etoDateStr+'.png'
                    infoJsonFile = acct_dir_img2 + acctName+'_'+fieldassetid+'__INFO_'+etoDateStr+'.json'

                    acct_dir_pkl2 = acct_dir_pkl + fieldassetid + '/'
                    if not os.path.exists(acct_dir_pkl2):
                        os.makedirs(acct_dir_pkl2)
                    pklFile = acct_dir_pkl2 + acctName+'_'+fieldassetid+'__PKL_'+etoDateStr+'.pkl'


                    imgDict = {
                        'ImgFileName': os.path.split(imgFile)[1],
                        'ImgDate': imgDtStr,
                        'ImgType': 'NDVI',
                        'ImgExtent': ndvi_ext, 
                        'ImgBoundary': P.to_wkt(),
                        'AssetID': row.fieldassetid,
                        'CustomerID': row.customerid,
                        'Ranch': ranch,
                        'Field': field
                    }
                    with open(infoJsonFile, 'w', encoding='utf8') as outfile:
                        json.dump(imgDict, outfile, ensure_ascii=False)

                    plt.ioff();
                    with open(imgFile, 'wb') as outfile:
                        fig = plt.figure()
                        plt.imshow(ndvi_img, extent=ext, cmap=zoom1, vmin=0, vmax=1)
                        fig.gca().set_axis_off()
                        plt.savefig(outfile, dpi=600, bbox_inches='tight', pad_inches=0, transparent=True)
                        plt.close();

                    with open(pklFile,'wb') as f:
                        pickle.dump(ndvi_img, f)

                    retdf = ndvi_df.loc[0]

                    ndvi_mean = retdf.Mean
                    ndvi_p90 = retdf.p90
                    ndvi_p25 = retdf.p25
                    ndvi_p75 = retdf.p75
                    ndvi_median = retdf.Median
                    ndvi_std = retdf.StdDev

                else:
                    #set values if NaN
                    ndvi_mean = np.nan
                    ndvi_p90 = np.nan
                    ndvi_p25 = np.nan
                    ndvi_p75 = np.nan
                    ndvi_median = np.nan
                    ndvi_std = np.nan


                #derived
                ndvi_unif = (ndvi_p25 / ndvi_median)*100

                #compute the METRIC ETc
                et_metric = 0.0 #Don't compute METRIC ET due to no rnet
                #computeMetricET(tavg,rnet,wavg,ndvi_mean)

                ##########if et_metric is NaN, make it zero
                if (np.isnan(et_metric)):
                    et_metric=0.0

                #the Kc
                fac = 1.35 #somewhat arbitrary for now
                k_c = 1.3558*ndvi_mean*fac + 0.05
                k_c2 = 1.3558*ndvi_mean + 0.05

                #compute the ETc
                et_ndvi = eto*k_c
                et_ndvi2 = eto*k_c2

                #get the final ETc
                # Include et_metric
                k_a = 1 #We can reduce k_a to increase the effect of Metric (see above)
                et_c =  (1.15 - k_a)*et_metric + ((k_a + 0.0)*et_ndvi)
                et_c2 =  (1.15 - k_a)*et_metric + ((k_a + 0.0)*et_ndvi2)

                k_c_p90 = ((1.3558*ndvi_p90*fac + 0.05) + (1.3558*ndvi_p90 + 0.05))/2
                et_p90 = eto*k_c_p90

                acctGDF_date.loc[i,'ETc']=et_c
                acctGDF_date.loc[i,'ETc_low']=et_c2
                acctGDF_date.loc[i,'Kc']=k_c
                acctGDF_date.loc[i,'ETm']=et_metric
                acctGDF_date.loc[i,'ET_p90']=et_p90
                acctGDF_date.loc[i,'ET_unif']=ndvi_unif
                acctGDF_date.loc[i,'Precipitation(in)']=precip
                acctGDF_date.loc[i,'MaximumAirTemperature(F)']=maxAirTemp
                acctGDF_date.loc[i,'MinimumAirTemperature(F)']=minAirTemp
                acctGDF_date.loc[i,'AverageAirTemperature(F)']=avgAirTemp
                acctGDF_date.loc[i,'MaximumRelativeHumidity(%)']=maxRelHumPercent
                acctGDF_date.loc[i,'MinimumRelativeHumidity(%)']=minRelHumPercent
                acctGDF_date.loc[i,'AverageRelativeHumidity(%)']=avgRelHumPercent
                acctGDF_date.loc[i,'DewPoint(F)']=dewPoint


            acctGDF_date['SentinelFootprint'] = acctGDF_date['SentinelFootprint'].astype('str')
            #acctGDF_date['LandsatFootprint'] = acctGDF_date['LandsatFootprint'].astype('str')
            acctGDF_date['Date'] = acctGDF_date['Date'].astype('str')
            if 'CentroidPoint' in acctGDF_date.columns:
                acctGDF_date.drop('CentroidPoint',inplace=True, axis=1)
            print(acctGDF_date.ETc.tolist())

            #save the values
            acctGDF_date.to_file(jsonFile, driver="GeoJSON")

        print("\tEND --- %s seconds ---" % (time.time() - starta))
        print('\tDone!')
        #except:
            #print('Exception in Acct: '+str(Exception))
            #continue

    print("\nEND --- %s seconds ---" % (time.time() - startt))
    return
"""

def generateETc_NOAA3_noscreen(merged_df):

    unique_accounts = pd.unique(merged_df.customerid)
    startt = time.time()
    #Now to put this into the output folder for each account
    print("Start --- %s  ---" % datetime.now())

    for acct in unique_accounts:
        a = 'JainlogicCustomerid_'+str(acct)
        #try:
        starta = time.time()
        
        print('\n\tWriting...%s ' % a)
        print("\tStart --- %s  ---" % datetime.now())
        acctName = a.upper()
        acct_dir = acctOutputPath + a + '/data/'
        acct_dir_img = acctOutputPath + a + '/img/'
        acct_dir_pkl = acctOutputPath + a + '/pkl/'
        acct_dir_rep = acctOutputPath + a + '/reports/'
            
        #create the directory if one doesn't exist
        if not os.path.exists(acct_dir):
            os.makedirs(acct_dir)
        if not os.path.exists(acct_dir_img):
            os.makedirs(acct_dir_img)
        if not os.path.exists(acct_dir_pkl):
            os.makedirs(acct_dir_pkl)
        if not os.path.exists(acct_dir_rep):
            os.makedirs(acct_dir_rep)
            
        #For each account, take the subset of the acct_df 
        acctDF = merged_df[merged_df.customerid == acct]
        
        acctGDF = gpd.GeoDataFrame(acctDF,geometry=acctDF.geometry)
        
        #for every day, do the following analysis and printout
        acct_dates_list = acctGDF.Date.tolist()
        unique_dts = pd.unique(acct_dates_list) # so we don't multiple count for each field
        unique_dts.sort()

        for dt in unique_dts:
            
            endDay = dt        
            
            startd = time.time()
            #print("\n\t\tSTART --- %s  ---" % datetime.now())

            etoDateStr =  datetime.strftime(dt,'%Y-%m-%d')
            jsonFile = acct_dir + acctName+'__ET_'+etoDateStr+'.json'
            
            acctGDF_date = acctGDF[acctGDF.Date==dt]
            acctGDF_date.reset_index(drop=True, inplace=True)
            
            print('\t\t'+etoDateStr)
            #print(acctGDF_date.columns)
            for i,row in acctGDF_date.iterrows():

                customerid = str(row.customerid)
                acct = row.Category
                ranch = row.SubCategory
                field = row.Name
                fieldassetid = 'ASSETID_'+str(row.fieldassetid)

                #ETo and weather values
                #print('row: '+str(row))
                if 'Tavg' in row.index.tolist():
                    eto = row['ETo(in)']
                    tavg = row['Tavg']
                    rnet = row['Rnet']
                    wavg = row['Wavg']
                    ############## add others Raw and Deried
                    precip = row['Precipitation(in)']
                    maxAirTemp = row['MaximumAirTemperature(F)']
                    minAirTemp = row['MinimumAirTemperature(F)']
                    avgAirTemp = row['AverageAirTemperature(F)']
                    maxRelHumPercent = row['MaximumRelativeHumidity(%)']
                    minRelHumPercent = row['MinimumRelativeHumidity(%)']
                    avgRelHumPercent = row['AverageRelativeHumidity(%)']
                    dewPoint = row['DewPoint(F)']
                else:
                    eto = row['ETo(in)']

                #########################################

                P = getSinglePoly(row.geometry)


                ss = row.SentinelFootprint
                ll = row.LandsatFootprint

                #do Sentinel
                if not type(ss)==list:
                    if (ss.find(',')==-1): #There is only ONE tile
                        S = [ss.replace('[','').replace(']','').replace('\'','')]
                    else:
                        tok = ss.split(',')
                        l = len(tok)
                        S=[]
                        for ii in range(l):
                            tk = tok[ii].replace('[','').replace(']','').replace('\'','')
                            S.append(tk)
                else:
                    S = ss
                if (type(S)==str):
                    S = [S]

                #do landsat
                if not type(ll)==list:
                    if (ss.find(',')==-1): #There is only ONE tile
                        L = ll.replace('[','').replace(']','').replace('\'','')
                    else:
                        tok = ll.split(', ')
                        l = len(tok)
                        L=[]
                        for ii in range(l):
                            tk = tok[ii].replace('[','').replace(']','').replace('\'','')
                            L.append(tk.replace(' ',''))
                else:
                    L = [str(x).zfill(6) for x in ll]

                if (type(L)==str):
                    L = [L]


                #print([S,L,etoDateStr])

                ndvi_mList_df = pd.DataFrame()

                NDVIFiles_df = getNDVIFilesWithDate(62, S, L, etoDateStr)
                NDVIFiles_df.reset_index(drop=True,inplace=True)

                #print(NDVIFiles_df)

                [ndvi_df,ndvi_img,ndvi_ext,imgDtStr] = getNDVIValueWithStatsAndImg3_noscreen(NDVIFiles_df,P)

                #print([ndvi_df,ndvi_ext,imgDtStr])

                if (len(ndvi_df)>0 & (len(ndvi_img)>0)): #The first time the df has valid values
                    ext = json.loads(ndvi_ext)
                    ndvi_mList_df = ndvi_df.copy()
                    
                    
                    acct_dir_img2 = acct_dir_img + fieldassetid + '/'
                    if not os.path.exists(acct_dir_img2):
                        os.makedirs(acct_dir_img2)
                    imgFile = acct_dir_img2 + acctName+'_'+fieldassetid+'__IMG_'+etoDateStr+'.png'
                    infoJsonFile = acct_dir_img2 + acctName+'_'+fieldassetid+'__INFO_'+etoDateStr+'.json'

                    acct_dir_pkl2 = acct_dir_pkl + fieldassetid + '/'
                    if not os.path.exists(acct_dir_pkl2):
                        os.makedirs(acct_dir_pkl2)
                    pklFile = acct_dir_pkl2 + acctName+'_'+fieldassetid+'__PKL_'+etoDateStr+'.pkl'
                    
                
                    imgDict = {
                        'ImgFileName': os.path.split(imgFile)[1],
                        'ImgDate': imgDtStr,
                        'ImgType': 'ETc',
                        'ImgExtent': ndvi_ext, 
                        'ImgBoundary': P.to_wkt(),
                        'AssetID': row.fieldassetid,
                        'CustomerID':a.split('_')[1],
                        'Ranch': ranch,
                        'Field': field
                    }
                    with open(infoJsonFile, 'w', encoding='utf8') as outfile:
                        json.dump(imgDict, outfile, ensure_ascii=False)

                    plt.ioff();
                    with open(imgFile, 'wb') as outfile:
                        fig = plt.figure()
                        plt.imshow(ndvi_img, extent=ext, cmap=zoom1, vmin=0, vmax=1)
                        fig.gca().set_axis_off()
                        plt.savefig(outfile, dpi=600, bbox_inches='tight', pad_inches=0, transparent=True)
                        plt.close();

                    with open(pklFile,'wb') as f:
                        pickle.dump(ndvi_img, f)

                    retdf = ndvi_df.loc[0]
                    ndvi_mean = retdf.Mean
                    ndvi_p90 = retdf.p90
                    ndvi_p25 = retdf.p25
                    ndvi_p75 = retdf.p75
                    ndvi_median = retdf.Median
                    ndvi_std = retdf.StdDev

                else:
                    #set values if NaN
                    ndvi_mean = np.nan
                    ndvi_p90 = np.nan
                    ndvi_p25 = np.nan
                    ndvi_p75 = np.nan
                    ndvi_median = np.nan
                    ndvi_std = np.nan

                #print(['NDVI_mean = '+str(ndvi_mean)])
                #derived
                ndvi_unif = (ndvi_p25 / ndvi_median)*100

                #compute the METRIC ETc
                if 'Tavg' in row.index.tolist():
                    et_metric = computeMetricET(tavg,rnet,wavg,ndvi_mean)
                else:
                    et_metric=0.0

                ##########if et_metric is NaN, make it zero
                if (np.isnan(et_metric)):
                    et_metric=0.0

                #print('ETc metric='+str(et_metric))
                #the Kc
                fac = 1.35 #somewhat arbitrary for now
                k_c = 1.3558*ndvi_mean*fac + 0.05
                k_c2 = 1.3558*ndvi_mean + 0.05

                #compute the ETc
                et_ndvi = eto*k_c
                et_ndvi2 = eto*k_c2

                #get the final ETc
                # Include et_metric
                k_a = 1 #We can reduce k_a to increase the effect of Metric (see above)
                et_c =  (1.15 - k_a)*et_metric + ((k_a + 0.0)*et_ndvi)
                et_c2 =  (1.15 - k_a)*et_metric + ((k_a + 0.0)*et_ndvi2)

                k_c_p90 = ((1.3558*ndvi_p90*fac + 0.05) + (1.3558*ndvi_p90 + 0.05))/2
                et_p90 = eto*k_c_p90

                #print([eto,k_c,k_c2,et_metric,et_c, et_c2, et_p90,ndvi_unif,precip])

                acctGDF_date.loc[i,'ETc']=et_c
                acctGDF_date.loc[i,'ETc_low']=et_c2
                acctGDF_date.loc[i,'Kc']=k_c
                acctGDF_date.loc[i,'ETm']=et_metric
                acctGDF_date.loc[i,'ET_p90']=et_p90
                acctGDF_date.loc[i,'ET_p25']=ndvi_p25
                acctGDF_date.loc[i,'ET_p75']=ndvi_p75
                acctGDF_date.loc[i,'ET_median']=ndvi_median
                acctGDF_date.loc[i,'ET_StdDev']=ndvi_std
                acctGDF_date.loc[i,'ET_unif']=ndvi_unif
                if 'Tavg' in row.index.tolist():
                    acctGDF_date.loc[i,'Precipitation(in)']=precip
                    acctGDF_date.loc[i,'MaximumAirTemperature(F)']=maxAirTemp
                    acctGDF_date.loc[i,'MinimumAirTemperature(F)']=minAirTemp
                    acctGDF_date.loc[i,'AverageAirTemperature(F)']=avgAirTemp
                    acctGDF_date.loc[i,'MaximumRelativeHumidity(%)']=maxRelHumPercent
                    acctGDF_date.loc[i,'MinimumRelativeHumidity(%)']=minRelHumPercent
                    acctGDF_date.loc[i,'AverageRelativeHumidity(%)']=avgRelHumPercent
                    acctGDF_date.loc[i,'DewPoint(F)']=dewPoint


            acctGDF_date['SentinelFootprint'] = acctGDF_date['SentinelFootprint'].astype('str')
            acctGDF_date['LandsatFootprint'] = acctGDF_date['LandsatFootprint'].astype('str')
            acctGDF_date['Date'] = acctGDF_date['Date'].astype('str')
            acctGDF_date.drop('CentroidPoint',inplace=True, axis=1)

            #print(type(acctGDF_date))
            print(acctGDF_date.ETc.tolist())

            #save the values
            #acctGDF_date = acctGDF_date.astype('str')
            #print(jsonFile)
            acctGDF_date.to_file(jsonFile, driver="GeoJSON")

        print("\tEND --- %s seconds ---" % (time.time() - starta))
        print('\tDone!')
        #except:
        #	print('Exception in Acct: '+str(a))
        #	continue

    print("\nEND --- %s seconds ---" % (time.time() - startt))
    return

def multiprocessing_etc_noaa(merged_df, dt):
    time.sleep(2)
    dt_new = dt.astype('M8[D]').astype('O')
    generate_MP_ETc_NOAA(merged_df,dt_new)
    return

def generate_MP_ETc_NOAA(merged_df,dt):

    unique_accounts = pd.unique(merged_df.customerid)
    startt = time.time()
    #Now to put this into the output folder for each account
    print("Start MP ETc--- %s  ---" % datetime.now())

    for acct in unique_accounts:
        a = 'JainlogicCustomerid_'+str(acct)
        #try:
        starta = time.time()
        
        print('\n\tWriting MP ETc...%s ' % a)
        print("\tStart --- %s  ---" % datetime.now())
        acctName = a.upper()
        acct_dir = acctOutputPath + a + '/data/'
        acct_dir_img = acctOutputPath + a + '/img/'
        acct_dir_pkl = acctOutputPath + a + '/pkl/'
        acct_dir_rep = acctOutputPath + a + '/reports/'
            
        #create the directory if one doesn't exist
        if not os.path.exists(acct_dir):
            os.makedirs(acct_dir)
        if not os.path.exists(acct_dir_img):
            os.makedirs(acct_dir_img)
        if not os.path.exists(acct_dir_pkl):
            os.makedirs(acct_dir_pkl)
        if not os.path.exists(acct_dir_rep):
            os.makedirs(acct_dir_rep)
            
        #For each account, take the subset of the acct_df 
        acctDF = merged_df[merged_df.customerid == acct]
        
        acctGDF = gpd.GeoDataFrame(acctDF,geometry=acctDF.geometry)
        
        #for every day, do the following analysis and printout
        acct_dates_list = acctGDF.Date.tolist()
        unique_dts = pd.unique(acct_dates_list) # so we don't multiple count for each field
        unique_dts.sort()
        
        #for dt in unique_dts:
            
        endDay = dt        
        
        startd = time.time()
        #print("\n\t\tSTART --- %s  ---" % datetime.now())

        etoDateStr =  datetime.strftime(dt,'%Y-%m-%d')
        jsonFile = acct_dir + acctName+'__ET_'+etoDateStr+'.json'
        
        acctGDF_date = acctGDF[acctGDF.Date==dt]
        acctGDF_date.reset_index(drop=True, inplace=True)
        
        print('\t\t'+etoDateStr)
        #print(acctGDF_date.columns)
        for i,row in acctGDF_date.iterrows():

            customerid = str(row.customerid)
            acct = row.Category
            ranch = row.SubCategory
            field = row.Name
            fieldassetid = 'ASSETID_'+str(row.fieldassetid)
            
            #ETo and weather values
            #print('row: '+str(row))

            eto = row['ETo(in)']
            tavg = row['Tavg']
            rnet = row['Rnet']
            wavg = row['Wavg']
            ############## add others Raw and Deried
            precip = row['Precipitation(in)']
            maxAirTemp = row['MaximumAirTemperature(F)']
            minAirTemp = row['MinimumAirTemperature(F)']
            avgAirTemp = row['AverageAirTemperature(F)']
            maxRelHumPercent = row['MaximumRelativeHumidity(%)']
            minRelHumPercent = row['MinimumRelativeHumidity(%)']
            avgRelHumPercent = row['AverageRelativeHumidity(%)']
            dewPoint = row['DewPoint(F)']
            
            #########################################
        
            # NDVI value
            P = getSinglePoly(row.geometry)

            
            ss = row.SentinelFootprint
            ll = row.LandsatFootprint
            
            #do Sentinel
            if not type(ss)==list:
                if (ss.find(',')==-1): #There is only ONE tile
                    S = [ss.replace('[','').replace(']','').replace('\'','')]
                else:
                    tok = ss.split(',')
                    l = len(tok)
                    S=[]
                    for ii in range(l):
                        tk = tok[ii].replace('[','').replace(']','').replace('\'','')
                        S.append(tk)
            else:
                S = ss
            if (type(S)==str):
                S = [S]
        
            #do landsat
            if not type(ll)==list:
                if (ss.find(',')==-1): #There is only ONE tile
                    L = ll.replace('[','').replace(']','').replace('\'','')
                else:
                    tok = ll.split(', ')
                    l = len(tok)
                    L=[]
                    for ii in range(l):
                        tk = tok[ii].replace('[','').replace(']','').replace('\'','')
                        L.append(tk.replace(' ',''))
            else:
                L = [str(x).zfill(6) for x in ll]

            if (type(L)==str):
                L = [L]
            

            #print([S,L,etoDateStr])

            ndvi_mList_df = pd.DataFrame()
            
            NDVIFiles_df = getNDVIFilesWithDate(32, S, L, etoDateStr)
            NDVIFiles_df.reset_index(drop=True,inplace=True)

            #print(NDVIFiles_df)
            
            [ndvi_df,ndvi_img,ndvi_ext,imgDtStr] = getNDVIValueWithStatsAndImg3(NDVIFiles_df,P)
            
            #print([ndvi_df,ndvi_ext,imgDtStr])
            
            if (len(ndvi_df)>0 & (len(ndvi_img)>0)): #The first time the df has valid values
                ext = json.loads(ndvi_ext)
                ndvi_mList_df = ndvi_df.copy()
                
                
                acct_dir_img2 = acct_dir_img + fieldassetid + '/'
                if not os.path.exists(acct_dir_img2):
                    os.makedirs(acct_dir_img2)
                imgFile = acct_dir_img2 + acctName+'_'+fieldassetid+'__IMG_'+etoDateStr+'.png'
                infoJsonFile = acct_dir_img2 + acctName+'_'+fieldassetid+'__INFO_'+etoDateStr+'.json'

                acct_dir_pkl2 = acct_dir_pkl + fieldassetid + '/'
                if not os.path.exists(acct_dir_pkl2):
                    os.makedirs(acct_dir_pkl2)
                pklFile = acct_dir_pkl2 + acctName+'_'+fieldassetid+'__PKL_'+etoDateStr+'.pkl'
                
            
                imgDict = {
                    'ImgFileName': os.path.split(imgFile)[1],
                    'ImgDate': imgDtStr,
                    'ImgType': 'ETc',
                    'ImgExtent': ndvi_ext, 
                    'ImgBoundary': P.to_wkt(),
                    'AssetID': row.fieldassetid,
                    'CustomerID':a.split('_')[1],
                    'Ranch': ranch,
                    'Field': field
                }
                with open(infoJsonFile, 'w', encoding='utf8') as outfile:
                    json.dump(imgDict, outfile, ensure_ascii=False)

                plt.ioff();
                with open(imgFile, 'wb') as outfile:
                    fig = plt.figure()
                    plt.imshow(ndvi_img, extent=ext, cmap=zoom1, vmin=0, vmax=1)
                    fig.gca().set_axis_off()
                    plt.savefig(outfile, dpi=600, bbox_inches='tight', pad_inches=0, transparent=True)
                    plt.close();
                    
                with open(pklFile,'wb') as f:
                    pickle.dump(ndvi_img, f)
            
                retdf = ndvi_df.loc[0]
                ndvi_mean = retdf.Mean
                ndvi_p90 = retdf.p90
                ndvi_p25 = retdf.p25
                ndvi_p75 = retdf.p75
                ndvi_median = retdf.Median
                ndvi_std = retdf.StdDev
                
            else:
                #set values if NaN
                ndvi_mean = np.nan
                ndvi_p90 = np.nan
                ndvi_p25 = np.nan
                ndvi_p75 = np.nan
                ndvi_median = np.nan
                ndvi_std = np.nan
                
            #print(['NDVI_mean = '+str(ndvi_mean)])
            #derived
            ndvi_unif = (ndvi_p25 / ndvi_median)*100
        
            #compute the METRIC ETc
            #print('Metric inputs:')
            #print([tavg,rnet,wavg,ndvi_mean])
            et_metric = computeMetricET(tavg,rnet,wavg,ndvi_mean)
            
            ##########if et_metric is NaN, make it zero
            if (np.isnan(et_metric)):
                et_metric=0.0

            #print('ETc metric='+str(et_metric))
            #the Kc
            fac = 1.35 #somewhat arbitrary for now
            k_c = 1.3558*ndvi_mean*fac + 0.05
            k_c2 = 1.3558*ndvi_mean + 0.05
        
            #compute the ETc
            et_ndvi = eto*k_c
            et_ndvi2 = eto*k_c2
        
            #get the final ETc
            # Include et_metric
            k_a = 1 #We can reduce k_a to increase the effect of Metric (see above)
            et_c =  (1.15 - k_a)*et_metric + ((k_a + 0.0)*et_ndvi)
            et_c2 =  (1.15 - k_a)*et_metric + ((k_a + 0.0)*et_ndvi2)
            
            k_c_p90 = ((1.3558*ndvi_p90*fac + 0.05) + (1.3558*ndvi_p90 + 0.05))/2
            et_p90 = eto*k_c_p90
            
            #print([eto,k_c,k_c2,et_metric,et_c, et_c2, et_p90,ndvi_unif,precip])

            acctGDF_date.loc[i,'ETc']=et_c
            acctGDF_date.loc[i,'ETc_low']=et_c2
            acctGDF_date.loc[i,'Kc']=k_c
            acctGDF_date.loc[i,'ETm']=et_metric
            acctGDF_date.loc[i,'ET_p90']=et_p90
            acctGDF_date.loc[i,'ET_unif']=ndvi_unif
            acctGDF_date.loc[i,'Precipitation(in)']=precip
            acctGDF_date.loc[i,'MaximumAirTemperature(F)']=maxAirTemp
            acctGDF_date.loc[i,'MinimumAirTemperature(F)']=minAirTemp
            acctGDF_date.loc[i,'AverageAirTemperature(F)']=avgAirTemp
            acctGDF_date.loc[i,'MaximumRelativeHumidity(%)']=maxRelHumPercent
            acctGDF_date.loc[i,'MinimumRelativeHumidity(%)']=minRelHumPercent
            acctGDF_date.loc[i,'AverageRelativeHumidity(%)']=avgRelHumPercent
            acctGDF_date.loc[i,'DewPoint(F)']=dewPoint
            

        acctGDF_date['SentinelFootprint'] = acctGDF_date['SentinelFootprint'].astype('str')
        acctGDF_date['LandsatFootprint'] = acctGDF_date['LandsatFootprint'].astype('str')
        acctGDF_date['Date'] = acctGDF_date['Date'].astype('str')
        acctGDF_date.drop('CentroidPoint',inplace=True, axis=1)
        print(acctGDF_date.ETc.tolist())
        
        #save the values
        acctGDF_date.to_file(jsonFile, driver="GeoJSON")

        #print("\tEND --- %s seconds ---" % (time.time() - starta))
        #print('\tDone!')
        #except:
        #	print('Exception in Acct: '+str(a))
        #	continue

    print("\nEND --- %s seconds ---" % (time.time() - startt))
    return

def downloadNewSentinelTiles(newSentinelTileList, N):
    #make directories if there aren't any

    init_list_tiles = newSentinelTileList
    #Create directories if not there
    for tile in init_list_tiles:
        rawTilePath = sentinelRawFolder + tile
        procTilePath = sentinelProcFolder + tile
        if not os.path.exists(rawTilePath):
            os.mkdir(rawTilePath)
        if not os.path.exists(procTilePath):
            os.mkdir(procTilePath)

    #Set the dates (LATER) from parameters
    yrStr = datetime.strftime(datetime.now(),'%Y')

    endDt = datetime.now() - timedelta(days=2)
    startDt = endDt - timedelta(days=int(N))

    list_of_ts = pd.date_range(startDt, endDt, freq='m')
    list_of_dates = [datetime.strptime(datetime.strftime(x,'%Y-%m-%d'), '%Y-%m-%d') for x in list_of_ts]
    for dt in list_of_dates:
        startDt = dt
        endDt = startDt + timedelta(days=30)

        ################ Run the numbers
        startDtStr = datetime.strftime(startDt, '%Y-%m-%d')
        endDtStr = datetime.strftime(endDt, '%Y-%m-%d')
        print('\nDownload Horizon from '+startDtStr+' to '+endDtStr)

        #Setup the query
        api = SentinelAPI('sumer_johal2019', 'AgraLogics!2019', 'https://scihub.copernicus.eu/dhus')
        #DISCARD tiles with more than 20% Cloud Cover
        query_kwargs = {
            'platformname': 'Sentinel-2',
            'producttype': 'S2MSI2A',
            'cloudcoverpercentage': (0, 75) #ignore if Cloud Cover > 75%
        }

        for tile in init_list_tiles: 
            try:
                print('\n\nProcessing Tile: '+tile)
                products = OrderedDict()
                kw = query_kwargs.copy()
                kw['date'] = (startDt, endDt)
                pp = api.query(getTileWkt(tile),**kw)
                products.update(pp)
                download_path = sentinelRawFolder + tile + '/' 
                pdf = pd.DataFrame(products).T
                pdf.reset_index(inplace=True)
                pdf['TILE'] = pdf.title.apply(lambda x: x.split('_')[5][-5:])    
                pdf = pdf[pdf.TILE==tile]
                pdf.set_index('index',drop=True, inplace=True)
                pdf.drop('TILE',inplace=True, axis=1)
                del pdf.index.name
                products_filtered = pdf.to_dict(into=OrderedDict, orient='index')

                #1. Download the Raw Data for the tile
                api.download_all(products_filtered, directory_path=download_path)

                #2. Save info to a pandas dataframe
                pdf = api.to_dataframe(products_filtered)

                #Now unzip
                for i,j in pdf.iterrows(): 
                    productID = pdf.identifier[i]
                    f = pdf.filename[i]

                    #3. Unzip the file
                    zipFile = download_path + pdf.title[i]
                    zf = zipFile + '.zip'
                    print("\n*******************\nExtracting file ...",zf)
                    with ZipFile(zf, 'r') as zipObj:
                        #Extract all the contents of zip file in current directory
                        zipObj.extractall(path = download_path)
                    print('...' + 'Done.')  

                    #4. Now delete the zip file
                    print('Deleting Zip file ...' + zf)
                    os.remove(zf)
                    print('...' + 'Done!')


                    ############ CHECK TO SEE IF THE FOOTPRINT IS VALID/BIG ENOUGH
                    footprint_path = download_path + f + '/'       
                    footprint_file = footprint_path + 'MTD_MSIL2A.xml'

                    fp = getFootprint(footprint_file)
                    area_footprint, coverageFraction = getAreaOfROI(fp)

                    if coverageFraction >= 0.25: #at least 25% the raster is there
                        print('\nCoverage is >= 25%: Coverage%= ' + str(coverageFraction*100))
                        #5. Convert to GEOTIFF & process for CC
                        tiffFile_path = download_path.replace("/RAW", "/PROCESSED")
                        #construct the entire path
                        dir_1 = download_path + f + '/GRANULE/'
                        d1 = os.listdir(dir_1)

                        gmlPath = dir_1 + d1[0] + '/QI_DATA/'

                        dir_2 = dir_1 + d1[0] + '/IMG_DATA/R10m/'
                        bands = os.listdir(dir_2)

                        #sort the list
                        bands.sort()
                        imagePath = dir_2

                        gmlPath = dir_1 + d1[0] + '/QI_DATA/'
                        gmlFile = 'MSK_CLDPRB_20m.jp2'
                        cloudMask = getCloudMask(gmlPath,gmlFile)


                        band2_src = rasterio.open(imagePath+bands[1], driver='JP2OpenJPEG')
                        blue = band2_src.read(1) * cloudMask #blue

                        band3_src = rasterio.open(imagePath+bands[2], driver='JP2OpenJPEG')
                        green = band3_src.read(1) * cloudMask #green

                        band4_src = rasterio.open(imagePath+bands[3], driver='JP2OpenJPEG')
                        red = band4_src.read(1) * cloudMask #red

                        band8_src = rasterio.open(imagePath+bands[4], driver='JP2OpenJPEG')
                        nir = band8_src.read(1) * cloudMask #nir
                        band8_src.close()

                        out_meta = band2_src.meta.copy()
                        epsg_code = int(band2_src.crs.data['init'][5:])

                        #Derived bands
                        ndvi = ((nir - red) / (nir + red))
                        ndvi_scaled = ndvi*255
                        ndvi_uint16 = ndvi_scaled.astype(np.uint16)

                        gndvi = ((nir - green) / (nir + green))
                        gndvi_scaled = gndvi*255
                        gndvi_uint16 = gndvi_scaled.astype(np.uint16)

                        vari = ((green - red) / (green + red - blue))
                        vari_scaled = vari*255
                        vari_uint16 = vari_scaled.astype(np.uint16)

                        savi = ((nir - red) / (nir + red + 0.5)) * 1.5
                        savi_scaled = savi*255
                        savi_uint16 = savi_scaled.astype(np.uint16)

                        evi = 2.5 * ((nir - red) / (nir + 6*red - 7.5*blue + 1))
                        evi_scaled = evi*255
                        evi_uint16 = evi_scaled.astype(np.uint16)





                        pcoverPercent = (np.count_nonzero(~np.isnan(ndvi)) / np.count_nonzero(~np.isnan(band2_src.read(1))))*100
                        pcPercent = str('{0:.2f}%'.format(pcoverPercent))

                        band2_src.close()
                        band3_src.close()
                        band4_src.close()
                        band8_src.close()

                        print('**********DATE: '+str(date.today())+', Percent CloudThrough:'+pcPercent)

                        #write the NDVI to a Geotiff for that scene
                        tifname = tiffFile_path + f.replace(".SAFE","") + '_NDVI_WGS84.tiff'    
                        print('Writing TIFF file: ' + tifname)
                        ndviTIFF = rasterio.open(tifname,'w',driver='Gtiff',
                                         width=band4_src.width, height=band4_src.height,
                                         count=1,
                                         crs=band4_src.crs,
                                         transform=band4_src.transform,
                                         dtype=band4_src.dtypes[0]
                                         )       
                        ndviTIFF.write(ndvi_uint16,1) #NDVI
                        ndviTIFF.close()

                        print('...Done!\n*******************\n\n')

                    else:
                        print('Coverage is < 25%: Coverage%= ' + str(coverageFraction*100))
                        #Delete the raw files
                        rmFolderName = download_path + f
                        print('Deleting .SAFE file folder ...' + f)
                        try:
                            shutil.rmtree(rmFolderName)
                        except OSError as e:
                            print ("Error: %s - %s." % (e.filename, e.strerror))
                        print('...' + 'Done!')
            except:
                print('Exception for tile: '+tile)
                continue

    return

def getUniqueSentinelTiles2(gdf_acct):
    list_of_sentinel_tiles = []
    sentinel_tiles = []
    for i,row in gdf_acct.iterrows():
        tileList = []
        fp = row.SentinelFootprint

        if type(fp)==list:
            sentinel_tiles =  sentinel_tiles + fp
        elif type(fp)==str:
            toks = fp.split(',')
            for i in range(len(toks)):
                tileList.append(str(toks[i]))
            sentinel_tiles =  sentinel_tiles + tileList

    list_of_sentinel_tiles =  list(set(sentinel_tiles))
    list_of_sentinel_tiles.sort()

    return list_of_sentinel_tiles


def getUniqueSentinelTiles(gdf_acct):
    list_of_sentinel_tiles = []
    sentinel_tiles = []
    for i,row in gdf_acct.iterrows():
        tileList = []
        fp = row.SentinelFootprint

        if type(fp)==list:
            tileList + fp
        elif type(fp)==str:
            toks = fp.split(',')
            for i in range(len(toks)):
                tileList.append(str(toks[i]))
            
        sentinel_tiles =  sentinel_tiles + tileList

    list_of_sentinel_tiles =  list(set(sentinel_tiles))
    list_of_sentinel_tiles.sort()

    """
    sentinel_tiles = []
    for i,row in gdf_acct.iterrows():
        tilelist = ast.literal_eval(str(row.SentinelFootprint))
        sentinel_tiles =  sentinel_tiles + tilelist

    list_of_sentinel_tiles =  list(set(sentinel_tiles))
    list_of_sentinel_tiles.sort()
    """
    return list_of_sentinel_tiles


#Make a list of the sentinel directories that already exist
existing_sentinel_list = [x for x in os.listdir(sentinelProcFolder)]
existing_landsat_list = [x for x in os.listdir(landsatProcFolder)]




def validate(email): 
    match=re.search(r"(^[a-zA-Z0-9_.+-]+@[a-zA-Z0-9]+\.[a-zA-Z0-9.]*\.*[com|org|edu]{3}$)",email)
    if match:
        return True
    else:
        return False


def fitBlankETValsInDateRange(DF, TypeStr, AcctName, RanchName,FieldName, fromDtstr, toDtStr):
    """
    fits a linear fit to de-NaN an ETc array of values
    """
    gdf = DF[(DF.SubCategory==RanchName) &
            (DF.Name==FieldName) &
            (DF.Category==AcctName)
           ]
    df = pd.DataFrame(gdf)
    #first make the Datetime searachable
    df.Date = pd.to_datetime(df.Date)
    
    #Try and fit the values with the week in question
    df = df[(df.Date >= fromDtStr) & (df.Date <= toDtStr)]
    df['ETc_average'] = df[['ETc', 'ETc_low']].mean(axis=1)
    ret_df = df.loc[:,['Date','ETc_average','ETc','ETc_low']]
    ret_df.sort_values(by='Date', inplace=True)
    ret_df.set_index('Date',inplace=True)
    
    if TypeStr=='ETc_low':
        mm='ETc_low'
    elif TypeStr=='ETc_high':
        mm = 'ETc'
    elif TypeStr=='ETc_mid':
        mm='ETc_average'
    elif TypeStr=='':
        mm='ETc_average'
    else:
        mm='ETc_average'
        
    
    ret_df['ETc_interp'] = ret_df[mm].interpolate(method='linear', limit_direction = 'both')
    ret_df.ETc_interp = ret_df.ETc_interp.bfill() #back fill any trailing NaNs
    ret_df.ETc_interp = ret_df.ETc_interp.interpolate(method='pad') #forward fill any trailing NaNs
    
    return ret_df

def annualETBudgetForAcct(acct_gdf):
    """
    makes a pivot table of the monthly budget and returns the table
    """
    #First get last year's dates
    curYrStr = datetime.strftime(datetime.today(),'%Y')
    prevYrStr = str(int(curYrStr)-1)
    
    prevYr_StartDateStr = prevYrStr+'-01-01'
    prevYr_StartDate = datetime.strptime(prevYr_StartDateStr,'%Y-%m-%d')
    
    prevYr_EndDateStr = prevYrStr+'-12-31'
    prevYr_EndDate = datetime.strptime(prevYr_EndDateStr,'%Y-%m-%d')
    
    
    tmp_df = acct_gdf.copy()
    tmp_df.Date = pd.to_datetime(tmp_df.Date, errors='coerse')
    tmp_df = tmp_df[(tmp_df.Date >=prevYr_StartDate) & (tmp_df.Date<=prevYr_EndDate)]
    tmp_df['Location'] = tmp_df['SubCategory'].map(str) + ': ' + tmp_df['Name'].map(str)
    tmp_df['ETc_average'] = tmp_df[['ETc', 'ETc_low']].mean(axis=1)
    tmp_df.set_index([tmp_df.index, tmp_df.Date], inplace=True)
    
    tmpResult = pd.pivot_table(data=tmp_df, 
                           columns=pd.Grouper(freq='M', level='Date', closed='right'), 
                           values=["ETc_average"],
                           index='Location', 
                           aggfunc=np.sum,
                          ).rename_axis(None).reset_index()
    tmpResult.columns = tmpResult.columns.droplevel(0)
    tmpResult.reset_index(drop=True, inplace=True)
    colHeadingsDates=tmpResult.columns[1:].tolist()
    colHeadingsMonth =[]
    
    for col in colHeadingsDates:
        colHeadingsMonth.append(col.strftime("%b"))
    colHeadingsPad = ['Location'] +colHeadingsMonth
    tmpResult.columns = colHeadingsPad
    tmpResult['Annual Total'] = tmpResult[colHeadingsMonth].sum(axis=1)
    
    return tmpResult


def uniformityTrend8Week(acct_gdf):
    #Get the uniformity trend for last 6 weeks
    #First get last year's dates
    curYrStr = datetime.strftime(datetime.today(),'%Y')
    
    currYr_StartDateStr = curYrStr+'-01-01'
    currYr_StartDate = datetime.strptime(currYr_StartDateStr,'%Y-%m-%d')
    
    currYr_EndDate = datetime.strftime(datetime.today(),'%Y-%m-%d')
    
    
    tmp_df = acct_gdf.copy()
    tmp_df.Date = pd.to_datetime(tmp_df.Date, errors='coerse')
    tmp_df = tmp_df[(tmp_df.Date >=currYr_StartDate) & (tmp_df.Date<=currYr_EndDate)]
    tmp_df['Location'] = tmp_df['SubCategory'].map(str) + ': ' + tmp_df['Name'].map(str)
    tmp_df.set_index([tmp_df.index, tmp_df.Date], inplace=True)
    
    tmpResult = pd.pivot_table(data=tmp_df, 
                           columns=pd.Grouper(freq='W', level='Date', closed='right'), 
                           values=["ET_unif"],
                           index='Location', 
                           aggfunc=np.mean,
                          ).rename_axis(None).reset_index()
    tmpResult.columns = tmpResult.columns.droplevel(0)
    tmpResult.columns = ['Location'] + tmpResult.columns[1:].tolist()
    tmpResult.reset_index(drop=True, inplace=True)
    tmpResult = tmpResult.set_index('Location')

    tr = tmpResult.iloc[:,-9:-1]
    
    slope = []
    for i, row in tr.iterrows():
        x = list(range(len(row)))
        y = row.tolist()
        slope.append(scipy.stats.linregress(x,y)[0])
    
    tr['8-Wk Average']=tr.mean(axis=1)
    tr['4-Wk Slope']=slope
    
    return tr


def waterAccounting(acct_gdf):
    #Get the uniformity trend for last 6 weeks
    #First get last year's dates
    curYrStr = datetime.strftime(datetime.today(),'%Y')
    
    currYr_StartDateStr = curYrStr+'-01-01'
    currYr_StartDate = datetime.strptime(currYr_StartDateStr,'%Y-%m-%d')
    
    currYr_EndDate = datetime.strftime(datetime.today(),'%Y-%m-%d')
    
    
    tmp_df = acct_gdf.copy()
    tmp_df.Date = pd.to_datetime(tmp_df.Date, errors='coerse')
    tmp_df = tmp_df[(tmp_df.Date >=currYr_StartDate) & (tmp_df.Date<=currYr_EndDate)]
    tmp_df['Location'] = tmp_df['SubCategory'].map(str) + ': ' + tmp_df['Name'].map(str)
    tmp_df.set_index([tmp_df.index, tmp_df.Date], inplace=True)
    
    tmpResult = pd.pivot_table(data=tmp_df, 
                           columns=pd.Grouper(freq='W', level='Date', closed='right'), 
                           values=["ETc_low"],
                           index='Location', 
                           aggfunc=np.sum,
                          ).rename_axis(None).reset_index()
    tmpResult.columns = tmpResult.columns.droplevel(0)
    tmpResult.columns = ['Location'] + tmpResult.columns[1:].tolist()
    tmpResult.reset_index(drop=True, inplace=True)
    tmpResult = tmpResult.set_index('Location')

    tr = tmpResult.iloc[:,-9:-1]
    #Convert Columns to Text Strings
    cols= list(tr.columns.values)
    for col in cols:
        d = datetime.strptime(str(col),'%Y-%m-%d %H:%M:%S')
        wkno = d.isocalendar()[1]
        colStr = 'Wk'+str(wkno)+' '+datetime.strftime(d,'%Y-%m-%d')
        tr.rename(columns={col:colStr}, inplace=True)
    
    tr_new = tr.copy().round(2)
    tr_new['YTDminus8W']= tmpResult.iloc[:, 0:-10].sum(axis=1)
    #reorder the columns so the last column is the first
    cols= list(tr_new.columns.values)
    nCols = [cols[-1],cols[0],cols[1],cols[2],cols[3],cols[4],cols[5],cols[6],cols[7]]
    tr_new = tr_new[nCols].round(2)
    
    
    
    tr_cumsum = tr_new.cumsum(axis=1).round(2)
    tr_cumsum = tr_cumsum.rename(columns={"YTDminus8W": "Wkminus9"})
    cols= list(tr_cumsum.columns.values)
    for col in cols:
        nCol = str(col)+'_YTD'
        tr_cumsum.rename(columns={col:nCol}, inplace=True)
        
    tr_res = tr_new.join(tr_cumsum, how='outer')
    
    
    #tr_res.assign(C=" ") #Empty Column Space
    #tr_res.rename(columns = {'C':' '}, inplace = True) #So this shows up as empty
    
    #NOW GET PRECIP
    precip = pd.pivot_table(data=tmp_df, 
                           columns=pd.Grouper(freq='A', level='Date', closed='right'), 
                           values=['Precipitation(in)'],
                           index='Location', 
                           aggfunc=np.sum,
                          ).rename_axis(None).reset_index()
    precip.columns = precip.columns.droplevel(0)
    precip.columns = ['Location'] + ['YTD Precipitation(in)']
    #precip.columns[1:].tolist()
    precip.reset_index(drop=True, inplace=True)
    precip.set_index('Location', inplace=True)
    precip=precip.round(2)
    
    tr_res = tr_res.join(precip, how='outer')
    
    #GET YTD ETc
    et = pd.pivot_table(data=tmp_df, 
                           columns=pd.Grouper(freq='A', level='Date', closed='right'), 
                           values=["ETc_low"],
                           index='Location', 
                           aggfunc=np.sum,
                          ).rename_axis(None).reset_index()
    et.columns = et.columns.droplevel(0)
    et.columns = ['Location'] + ['YTD ETc(in)']
    et.reset_index(drop=True, inplace=True)
    et.set_index('Location', inplace=True)
    #et.assign(C=" ") #Empty Column Space
    #et.rename(columns = {'C':' '}, inplace = True) #So this shows up as empty
    
    tr_res = tr_res.join(et, how='outer')
    
    
    return tr_res


def makeTrueRedFalseBlank(html_str):
    """
    Colors elements in a dateframe
    green if positive and red if
    negative. Does not color NaN
    values.
    """
    
    old_str = '-'
    new_str = '<font color=\"red\">'+old_str+'</font>'
    
    old_str2 = 'False'
    new_str2 = ''
    
    ret_str = html_str.replace(old_str,new_str)
    ret_str2 = ret_str.replace(old_str2,new_str2)
    
    return ret_str2


def color_negative_red(value):
    """
    Colors elements in a dateframe
    green if positive and red if
    negative. Does not color NaN
    values.
    """
    if value < 0:
        color = 'red'
    elif value > 0:
        color = 'green'
    else:
        color = 'black'

    return 'color: %s' % color


def getAcresForField(Ranch, Field, acctProfile_df):
    acres = acctProfile_df[(acctProfile_df.SubCategory==Ranch) & (acctProfile_df.Name==Field)].Acres.tolist()[0]
    return acres

def getAcresForField2(fieldassetid, acct_gdf):
    #returns area in Acres
    acct_gdf.crs = {'init' :'epsg:4326'}
    f = acct_gdf[acct_gdf.fieldassetid==str(fieldassetid)].reset_index(drop=True).geometry
    area_ac = area(mapping(f[0]))* 0.000247105
    
    return area_ac

def getAcresForField3(fieldassetid, acct_gdf):
    #returns area in Acres
    acct_gdf.Name = acct_gdf.Name.str.upper()
    fieldassetid = fieldassetid.upper()
    acct_gdf.crs = {'init' :'epsg:4326'}
    f = acct_gdf[acct_gdf.Name==str(fieldassetid)].reset_index(drop=True).geometry
    area_ac = area(mapping(f[0]))* 0.000247105
    
    return area_ac

def doHyperGrowFTP(File):
    url = 'ftp.agralogics.com'
    username = 'agls@agralogics.com'
    passwd = 'AglsAdmin!2019'
    ftpFolder='HyperGrow/'
    file_path = Path(File)

    with FTP(url, username, passwd, timeout=100000) as ftp, open(file_path, 'rb') as file:
        ftp.storbinary(f'STOR {ftpFolder+file_path.name}', file)
    return

def doHyperGrowFTPToFolder(folder, File):
    url = 'ftp.agralogics.com'
    username = 'agls@agralogics.com'
    passwd = 'AglsAdmin!2019'
    ftpFolder=folder+'/'

    file_path = Path(File)

    with FTP(url, username, passwd) as ftp, open(file_path, 'rb') as file:
        ftp.storbinary(f'STOR {ftpFolder+file_path.name}', file)
    
    return

# In[830]:


def getLinkTo(imgFile):
    
    file_path = Path(imgFile)
    linkURL='https://www.agralogics.com/HyperGrow/'+file_path.name
    linkStr = '<a href=\"'+linkURL+'\">HyperGrow</a>'  
    
    return linkStr

def getHTMLLinkTo(htmlFile):
    
    file_path = Path(htmlFile)
    linkURL='https://www.agralogics.com/HyperGrow/'+file_path.name
    linkStr = '<a href=\"'+linkURL+'\">FieldMap</a>'  
    
    return linkStr

# In[831]:


def genEmbededLink(filePath):

    pre = '<object data=\"'+filePath+'\" type=\"image/jpeg\" title=\"HyperGrow\" width=\"500\" height=\"720\">'
    ref = '<a href=\"'+filePath+'\">HyperGrow</a>'
    post='</object>'
    
    linkedStrAsCellValue = pre+ref+post
    
    return linkedStrAsCellValue


# In[832]:


#find the most recent pickle files for each field
def getNearestImageFile(dt,ranch,field):
    N=30
    to_dt = dt
    from_dt = to_dt-timedelta(days=N)

    subset_df = imgfiles_df[(imgfiles_df.Date>=from_dt) &
                      (imgfiles_df.Date<=to_dt) &
                      (imgfiles_df.SubCategory==ranch) &
                      (imgfiles_df.Name==field)]
    subset_df.sort_values(by=['Date'], inplace=True)
    subset_df = subset_df.reset_index(drop=True)
    fp = subset_df.loc[(len(subset_df)-1),'filePath'] 
    
    return fp


# In[833]:


#find the most recent pickle files for each field
def getNearestPklFile(dt,ranch,field):
    N=30
    to_dt = dt
    from_dt = to_dt-timedelta(days=N)

    subset_df = pklfiles_df[(pklfiles_df.Date>=from_dt) &
                      (pklfiles_df.Date<=to_dt) &
                      (pklfiles_df.SubCategory==ranch) &
                      (pklfiles_df.Name==field)]
    subset_df.sort_values(by=['Date'], inplace=True)
    subset_df = subset_df.reset_index(drop=True)
    #print(pklfiles_df)
    fp = subset_df.loc[(len(subset_df)-1),'filePath'] 
    
    return fp


# In[56]:


def getPrevYrStartEndDates(curYrStartDtStr,curYrEndDtStr):
    startDtStr_prev = endDtStr_prev = ''
    curStartDt = datetime.strptime(curYrStartDtStr,'%Y-%m-%d')
    curEndDt = datetime.strptime(curYrEndDtStr,'%Y-%m-%d')

    #startDtStr_prev
    if ((curStartDt.month==2) & (curStartDt.day==29)):
        day='28'
        startDtStr_prev = str(curStartDt.year-1) + '-' + str(curStartDt.month).rjust(2,'0') + '-' + day
    else:
        startDtStr_prev = str(curStartDt.year-1) + '-' + str(curStartDt.month).rjust(2,'0') + '-' + str(curStartDt.day).rjust(2,'0')
    
    #endDtStr_prev
    if ((curEndDt.month==2) & (curEndDt.day==29)):
        day='28'
        endDtStr_prev = str(curEndDt.year-1) + '-' + str(curEndDt.month).rjust(2,'0') + '-' + day
    else:
        endDtStr_prev = str(curEndDt.year-1) + '-' + str(curEndDt.month).rjust(2,'0') + '-' + str(curEndDt.day).rjust(2,'0')

    return [startDtStr_prev,endDtStr_prev]

def getCumStartEndDates(currDtStr,reportEndDate):
    startDtStr = endDtStr = ''
    currDt = datetime.strptime(currDtStr,'%Y-%m-%d')

    startDtStr = str(currDt.year) + '-01-01'
    endDtStr = datetime.strftime(reportEndDate,'%Y-%m-%d')

    return [startDtStr,endDtStr]


def getPrevYrCumStartEndDates(currDtStr,reportEndDate):
    startDtStr = endDtStr = ''
    currDt = datetime.strptime(currDtStr,'%Y-%m-%d')

    if ((currDt.month==2) & (currDt.day==29)):
        day='28'
    elif ((currDt.month==2) & (currDt.day==28) & (((currDt.year-1)%4)==0)):
        day='29'
    else:
        day = str(reportEndDate.day)

    startDtStr = str(currDt.year-1) + '-01-01'
    endDtStr = str(reportEndDate.year-1) + '-' + str(reportEndDate.month).rjust(2,'0') + '-' + day.rjust(2,'0')

    return [startDtStr,endDtStr]

def getPrev7DayPeriod(currDtStr):
    startDtStr = endDtStr = ''
    currDt = datetime.strptime(currDtStr,'%Y-%m-%d')
    startDt = currDt - timedelta(days=7)
    endDt = currDt - timedelta(days=1)

    startDtStr = datetime.strftime(startDt,'%Y-%m-%d')
    endDtStr = datetime.strftime(endDt,'%Y-%m-%d')
    return [startDtStr, endDtStr]

def getPrev15DayPeriod(currDtStr):
    startDtStr = endDtStr = ''
    currDt = datetime.strptime(currDtStr,'%Y-%m-%d')
    prev = currDt.replace(day=1) - timedelta(days=1)
    if (currDt.day > 15):
        first_day_of_current_month = currDt.replace(day=1)
        startDtStr = datetime.strftime(first_day_of_current_month,'%Y-%m-%d')
        endDtStr = datetime.strftime(first_day_of_current_month + timedelta(days=14),'%Y-%m-%d')
    else:
        startDtStr = datetime.strftime(prev.replace(day=1)+ timedelta(days=15),'%Y-%m-%d')
        endDtStr = datetime.strftime(prev,'%Y-%m-%d')

    return [startDtStr, endDtStr]

def getPrevYr15DayPeriod(currDtStr):
    startDtStrPrevYr = endDtStrPrevYr = ''
    currDt = datetime.strptime(currDtStr,'%Y-%m-%d')

    [startDtStr_currYr, endDtStr_currYr] = getPrev15DayPeriod(currDtStr)
    startDt_currYr = datetime.strptime(startDtStr_currYr,'%Y-%m-%d')
    endDt_currYr = datetime.strptime(endDtStr_currYr,'%Y-%m-%d')

    #print(currDt)

    if (currDt.day > 15):
        startDtStr = str(startDt_currYr.year-1) + '-' + str(startDt_currYr.month).rjust(2,'0') + '-' + '01'
        endDtStr = str(endDt_currYr.year-1) + '-' + str(endDt_currYr.month).rjust(2,'0') + '-' + '15'
        startDtStrPrevYr = startDtStr
        endDtStrPrevYr = endDtStr
    else:
        yr = startDt_currYr.year-1
        mo = startDt_currYr.month
        startDtStr = (str(yr) + '-' + str(mo).rjust(2,'0') + '-' + '16')
        startDtStrPrevYr = startDtStr

        firstDayOfMonStr = str(currDt.year-1) + '-' + str(currDt.month).rjust(2,'0') + '-' + '01'
        first_day_of_month = datetime.strptime(firstDayOfMonStr,'%Y-%m-%d')
        last_day_of_month = first_day_of_month - timedelta(days=1)
        endDtStrPrevYr = datetime.strftime(last_day_of_month,'%Y-%m-%d')

    return [startDtStrPrevYr, endDtStrPrevYr]


# In[118]:
def averageETos(df):
    eto_df = df[['ETo(in)','ReferenceETo(in)']].astype('float')
    eto_df['ReferenceETo(in)'] = eto_df.mean(skipna = True, axis=1)
    
    new_df = df.copy()
    new_df['ReferenceETo(in)'] = eto_df['ReferenceETo(in)']
    
    return new_df

def averageETs(currYr_df,prevYr_df):
    et_currYr_df = currYr_df[['ET_p90','ET_unif','ETc','ETc_low','ETm']].astype('float')
    et_prevYr_df = prevYr_df[['ET_p90','ET_unif','ETc','ETc_low','ETm']].astype('float')
    
    et_currYr_df[et_currYr_df.isnull()] = et_prevYr_df #update the NaN values with prevYr's values
    
    et_currYr_df = et_currYr_df.interpolate() #Interpolate any NaNs still there.
    
    new_df = currYr_df.copy()
    new_df['ET_p90'] = et_currYr_df['ET_p90']
    new_df['ET_unif'] = et_currYr_df['ET_unif']
    new_df['ETc'] = et_currYr_df['ETc']
    new_df['ETc_low'] = et_currYr_df['ETc_low']
    new_df['ETm'] = et_currYr_df['ETm']
    
    return new_df

def interpolateAndFill(currYr_df,prevYr_df,reportStartDateStr,reportEndDateStr):
    reportStartDt = datetime.strptime(reportStartDateStr,'%Y-%m-%d')
    reportEndDt = datetime.strptime(reportEndDateStr,'%Y-%m-%d')
    
    [prevYrStartDateStr,preYrEndDateStr] = getPrevYrStartEndDates(reportStartDateStr,reportEndDateStr)
    prevYrReportStartDt = datetime.strptime(prevYrStartDateStr,'%Y-%m-%d')
    prevYrReportEndDt = datetime.strptime(preYrEndDateStr,'%Y-%m-%d')
    
    field_df = averageETos(currYr_df) #fix the ETo
    field_df = averageETs(field_df,prevYr_df) #fix the ETc numbers with prev Yr if NaN
    
    return field_df


def getStartEndDates(adf):
    startDtStr = datetime.strftime(min(adf.Date),'%Y-%m-%d')
    endDtStr = datetime.strftime(max(adf.Date),'%Y-%m-%d')
    
    return [startDtStr, endDtStr]


#IF THE ET IS NULL, THEN USE (ET FROM LAST YEAR * RATIO OF ETO) AS THESE ARE PERMANENT ORCHARDS
def fixDF(adf):
    [reportStartDateStr,reportEndDateStr]=getStartEndDates(adf)
    
    [startDtStr_prev,endDtStr_prev] = getPrevYrStartEndDates(reportStartDateStr,reportEndDateStr)
    prevYr_adf = acct_gdf[(acct_gdf.Date >= startDtStr_prev) & (acct_gdf.Date <= endDtStr_prev)]
    prevYr_adf.reset_index(drop=True,inplace=True)
    
    list_of_unique_fields = pd.unique(adf.Name)
    a_new_gdf = pd.DataFrame() 
    for field in list_of_unique_fields:
        currYr_df = adf[adf.Name==field]
        currYr_df.reset_index(drop=True,inplace=True)
    
        prevYr_df = prevYr_adf[prevYr_adf.Name==field]
        prevYr_df.reset_index(drop=True,inplace=True)
    
        df = interpolateAndFill(currYr_df,prevYr_df,reportStartDateStr,reportEndDateStr)

        a_new_gdf = a_new_gdf.append(df, ignore_index=True)
    
    return a_new_gdf

def getAGLS_emailText():
    
    it = 'Dear Agralogics Customer, <br><br>Please see your HyperGrow* report below.'
    et = 'Email us with any questions or concerns you might have support@jagralogics.com .' + '<br><br>' + 'Sincerely,' + '<br><br>' + 'The Agralogics Team'
    
    return [it,et]
def getJAIN_emailText():
    
    it = 'Dear Jain Irrigation Customer, <br><br>Please see your HyperGrow* report below powered by Agralogics,Inc.'
    et = 'Email us with any questions or concerns you might have jiisupport@jainsusa.com .' + '<br><br>' + 'Sincerely,' + '<br><br>' + 'Jain Agricultural Services'
    
    return [it,et]

def getPAG_emailText():
    it = 'Dear Principle Ag Customer, <br><br>Please see your HyperGrow* report below, powered by Agralogics,Inc.'
    et = 'Email us with any questions or concerns you might have.'+ '<br><br>' + 'Jens 509-760-2572 jens@principleag.com'+ '<br>' +'Mike 509-989-5955 mike@principleag.com'+ '<br><br>' + 'Sincerely,' + '<br><br>' + 'The Principle Ag Team'
    #et = 'Email us with any questions or concerns you might have support@jagralogics.com .' + '<br><br>' + 'Sincerely,' + '<br><br>' + 'The Agralogics Team'
    return [it,et]

def getLinkToLogo(pngFile):
    file_path = Path(pngFile)
    linkURL='https://www.agralogics.com/HyperGrow/'+file_path.name
    linkStr = '<img src=\"'+linkURL+'\" alt="Image" height="100" width="300">'
    return linkStr

def getLinkToAGLSLogo(pngFile):
    file_path = Path(pngFile)
    linkURL='https://www.agralogics.com/HyperGrow/'+file_path.name
    linkStr = '<img src=\"'+linkURL+'\" alt="Image">'
    return linkStr

def getLinkToPAGLogo(pngFile):
    file_path = Path(pngFile)
    linkURL='https://www.agralogics.com/HyperGrow/'+file_path.name
    linkStr = '<img src=\"'+linkURL+'\" alt="Image">'
    return linkStr

def saveHyperGrowFig(img,cmp,vmin,vmax,imgOutFile):
    fig = plt.figure()
    plt.imshow(img,cmap=cmp, vmin=vmin,vmax=vmax)
    gca().set_xticklabels(['']*img.shape[0])
    gca().set_yticklabels(['']*img.shape[1])
    fig.gca().set_axis_off()
    fig.patch.set_facecolor('#ffffff')
    fig.patch.set_alpha(0)
    fig.savefig(imgOutFile, dpi=600, bbox_inches='tight', pad_inches=0, transparent=True)
    return

# In[13]:
def makeFieldMap(imgFile2,acctProfile_df):
    #try:
    ranchNameUpper = imgFile2.split('/')[8].split('_')[0]
    fieldNameUpper = imgFile2.split('/')[8].split('_')[1]
    imgName = imgFile2.split('/')[8].split('__')[0]
    dt = imgFile2.split('.')[0].split('_')[-1]
    #print(ranchNameUpper, fieldNameUpper, dt)

    #get the bounding box of the img
    gdf = acctProfile_df[(acctProfile_df.Name.str.replace(' ','').str.upper()==fieldNameUpper) & 
                (acctProfile_df.SubCategory.str.replace(' ','').str.upper()==ranchNameUpper)]
    gdf.reset_index(inplace=True,drop=True)
    gdf.crs = {'init': 'epsg:4326', 'no_defs': True}

    #gdf['geometry'] = gdf['geometry'].to_crs(epsg=3857)
    #gdf = gdf.to_crs({'init': 'epsg:3857'})

    fieldPoly = gdf.geometry[0]
    #set the extent
    min_lon = fieldPoly.bounds[0]
    max_lon = fieldPoly.bounds[2]
    min_lat = fieldPoly.bounds[1]
    max_lat = fieldPoly.bounds[3]

    center_lat = fieldPoly.centroid.coords.xy[1][0]
    center_lon = fieldPoly.centroid.coords.xy[0][0]
    # create the map
    m = folium.Map(location=[center_lat, center_lon],
           tiles='https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}', 
           projection=4326,
           attr='Rendered by AGRALOGICS using ESRI & GeoEye Basemaps',
           zoom_start = 15)

    folium.raster_layers.ImageOverlay(
    image=imgFile2,
    name=imgName,
    bounds=[[min_lat, min_lon], [max_lat, max_lon]],
    opacity=1,
    interactive=False,
    cross_origin=False,
    zindex=1,
    alt=imgName
    ).add_to(m)

    #folium.LayerControl().add_to(m)
    outFile = '/'.join(imgFile2.replace('/images/','/reports/').split('/')[0:-1])+'/'+imgName+'_HyperGrow_'+dt+'.html'
    m.save(outFile)


    of = outFile
    
    #except:
    #    print('Exception in makeFieldMap()')
    #    of = []
   
    return of


def makeFieldMapJL(imgFile2,acctProfile_df):
    #try:
    customerid = imgFile2.split('/')[-1].split('_')[1]
    fieldassetid = imgFile2.split('/')[-1].split('_')[3]
    imgName = imgFile2.split('/')[-1].split('__')[0]
    dt = imgFile2.split('.')[0].split('_')[-1]

    #get the bounding box of the img
    gdf = acctProfile_df[(acctProfile_df.customerid.str.replace(' ','')==customerid) & 
                (acctProfile_df.fieldassetid.str.replace(' ','')==fieldassetid)]
    gdf.reset_index(inplace=True,drop=True)
    gdf.crs = {'init': 'epsg:4326', 'no_defs': True}

    fieldPoly = gdf.geometry[0]
    #set the extent
    min_lon = fieldPoly.bounds[0]
    max_lon = fieldPoly.bounds[2]
    min_lat = fieldPoly.bounds[1]
    max_lat = fieldPoly.bounds[3]

    center_lat = fieldPoly.centroid.coords.xy[1][0]
    center_lon = fieldPoly.centroid.coords.xy[0][0]
    # create the map
    m = folium.Map(location=[center_lat, center_lon],
           tiles='https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}', 
           projection=4326,
           attr='Rendered by AGRALOGICS, Inc. using ESRI & GeoEye Basemaps',
           zoom_start = 15)

    folium.raster_layers.ImageOverlay(
    image=imgFile2,
    name=imgName,
    bounds=[[min_lat, min_lon], [max_lat, max_lon]],
    opacity=1,
    interactive=False,
    cross_origin=False,
    zindex=1,
    alt=imgName
    ).add_to(m)

    #folium.LayerControl().add_to(m)
    outFile = tmpDir+imgName+'_HyperGrow_'+dt+'.html'
    m.save(outFile)


    of = outFile
    
    #except:
    #    print('Exception in makeFieldMap()')
    #    of = []
   
    return of

def makeFieldMapPAG(imgFile2,acctProfile_df):
    #try:
    customerid = imgFile2.split('/')[-1].split('_')[0]
    fieldassetid = '_'.join(imgFile2.split('/')[-1].split('__')[0].split('_')[1:])
    imgName = imgFile2.split('/')[-1].split('__')[0].upper()
    dt = imgFile2.split('.')[0].split('_')[-1]

    #print([customerid, fieldassetid, imgName, dt])
    #print(acctProfile_df.iloc[0])
    #print(acctProfile_df.columns)

    #get the bounding box of the img

    gdf = acctProfile_df[(acctProfile_df.Category.str.replace(' ','')==customerid.replace(' ','')) & 
                (acctProfile_df.Name.str.replace(' ','')==fieldassetid.replace(' ',''))]
    gdf.reset_index(inplace=True,drop=True)
    gdf.crs = {'init': 'epsg:4326', 'no_defs': True}

    #print(gdf)

    fieldPoly = gdf.geometry[0]
    #set the extent
    min_lon = fieldPoly.bounds[0]
    max_lon = fieldPoly.bounds[2]
    min_lat = fieldPoly.bounds[1]
    max_lat = fieldPoly.bounds[3]

    center_lat = fieldPoly.centroid.coords.xy[1][0]
    center_lon = fieldPoly.centroid.coords.xy[0][0]
    # create the map
    m = folium.Map(location=[center_lat, center_lon],
           tiles='https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}', 
           projection=4326,
           attr='Rendered by AGRALOGICS, Inc. using ESRI & GeoEye Basemaps',
           zoom_start = 15)

    folium.raster_layers.ImageOverlay(
    image=imgFile2,
    name=imgName,
    bounds=[[min_lat, min_lon], [max_lat, max_lon]],
    opacity=1,
    interactive=False,
    cross_origin=False,
    zindex=1,
    alt=imgName
    ).add_to(m)

    #folium.LayerControl().add_to(m)
    outFile = tmpDir+imgName+'_HyperGrow_'+dt+'.html'
    m.save(outFile)


    of = outFile
    
    #except:
    #    print('Exception in makeFieldMap()')
    #    of = []
   
    return of

def getHGReportDates(acct):
    rpt_endStr = '_HyperGrow.json'
    rptPath = acctReportRoot + 'JainlogicCustomerid_'+str(acct)+'/reports/'
    list_of_HG_endDtStrs = [datetime.strptime(x.split('_')[3], '%Y-%m-%d')  for x in os.listdir(rptPath) if x.endswith(rpt_endStr)]
    latestEndDt = max(list_of_HG_endDtStrs)
    latestStartDt = latestEndDt - timedelta(days=6)
    endDtStr = datetime.strftime(latestEndDt, '%Y-%m-%d')
    startDtStr = datetime.strftime(latestStartDt, '%Y-%m-%d')
    
    return startDtStr, endDtStr

def getHVReportDate(acct):
    rpt_endStr = '_HyperView.json'
    rptPath = acctReportRoot + 'JainlogicCustomerid_'+str(acct)+'/reports/'
    list_of_HV_endDtStrs = [datetime.strptime(x.split('_')[2], '%Y-%m-%d')  for x in os.listdir(rptPath) if x.endswith(rpt_endStr)]
    latestEndDt = max(list_of_HV_endDtStrs)
    latestStartDt = latestEndDt - timedelta(days=6)
    endDtStr = datetime.strftime(latestEndDt, '%Y-%m-%d')
    
    return endDtStr


def getAGLSHGReportDates(acct):
    rpt_endStr = '_HyperGrow.json'
    rptPath = agls_acctOutputPath + str(acct)+'/reports/'
    list_of_HG_endDtStrs = [datetime.strptime(x.split('_')[3], '%Y-%m-%d')  for x in os.listdir(rptPath) if x.endswith(rpt_endStr)]
    latestEndDt = max(list_of_HG_endDtStrs)
    latestStartDt = latestEndDt - timedelta(days=6)
    endDtStr = datetime.strftime(latestEndDt, '%Y-%m-%d')
    startDtStr = datetime.strftime(latestStartDt, '%Y-%m-%d')
    
    return startDtStr, endDtStr

def getAGLSHVReportDate(acct):
    rpt_endStr = '_HyperView.json'
    rptPath = agls_acctOutputPath +str(acct)+'/reports/'
    list_of_HV_endDtStrs = [datetime.strptime(x.split('_')[2], '%Y-%m-%d')  for x in os.listdir(rptPath) if x.endswith(rpt_endStr)]
    latestEndDt = max(list_of_HV_endDtStrs)
    latestStartDt = latestEndDt - timedelta(days=6)
    endDtStr = datetime.strftime(latestEndDt, '%Y-%m-%d')
    
    return endDtStr	

def getRegionForJLAcctid(p):
    region=''
    
    if p.within(cap):
        region='CA'
    elif p.within(azp):
        region='AZ'
    elif p.within(wap):
        region='WA'
    elif p.within(txp):
        region='TX'
    elif p.within(flp):
        region='FL'
    else:
        for ausp in ausp_list:
            if p.within(ausp):
                region='AUS'
                break
            else:
                continue
    if not len(region)>0:
        region='NA'

    return region

def generateAndMergeForecastETo2(acct_gdf,Weather2020_gdf,startDtStr, endDtStr):
    start = datetime.strptime(startDtStr,'%Y-%m-%d')
    end = datetime.strptime(endDtStr,'%Y-%m-%d')
    #-timedelta(days=1) #This date is NOT included

    list_of_dates = [start + timedelta(days=x) for x in range(0, (end-start).days+1)]

    startt = time.time()
    EToForAccounts_2020 = pd.DataFrame()
    for i,row in acct_gdf.iterrows():
        print('\tAcct: '+str(row.Category)+', Ranch:'+row.SubCategory+', Field:'+row.Name)
        start_time = time.time()
        c = row.CentroidPoint
        if (type(c)==str):
            c = wkt.loads(row.CentroidPoint)

        etv = getForecastEToForDates(c,list_of_dates,Weather2020_gdf)

        #now add the columns for the other fields to etv
        etv['UUID']=row.UUID
        etv['Category']=row.Category
        etv['SubCategory']=row.SubCategory
        etv['Name']=row.Name

        EToForAccounts_2020 = EToForAccounts_2020.append(etv, ignore_index = True)
    print("ETo Done --- %s seconds ---" % (time.time() - startt))
    EToForAccounts = EToForAccounts_2020.rename(columns={"ReferenceETo(in)": "ETo(in)"})
    merged_df = pd.merge(acct_gdf, EToForAccounts, on=['Category','SubCategory','Name', 'UUID'])
    merged_df.reset_index() #reset the index to avoid any export issues
    #merged_df['Category'] = 'JainlogicCustomerid_' + merged_df['customerid'].astype(str)
    merged_df['Date'] = pd.to_datetime(merged_df['Date'])

    return merged_df


def generateAndMergeForecastETo(acct_gdf,Weather2020_gdf,startDtStr, endDtStr):
    start = datetime.strptime(startDtStr,'%Y-%m-%d')
    end = datetime.strptime(endDtStr,'%Y-%m-%d')
    #-timedelta(days=1) #This date is NOT included

    list_of_dates = [start + timedelta(days=x) for x in range(0, (end-start).days+1)]

    startt = time.time()
    EToForAccounts_2020 = pd.DataFrame()
    for i,row in acct_gdf.iterrows():
        print('\tAcct: '+str(row.customerid)+', Asset: '+str(row.fieldassetid)+', Ranch:'+row.SubCategory+', Field:'+row.Name)
        start_time = time.time()
        c = row.CentroidPoint
        if (type(c)==str):
            c = wkt.loads(row.CentroidPoint)

        etv = getForecastEToForDates(c,list_of_dates,Weather2020_gdf)

        #now add the columns for the other fields to etv
        etv['UUID']=row.UUID
        etv['Category']=row.Category
        etv['SubCategory']=row.SubCategory
        etv['Name']=row.Name

        EToForAccounts_2020 = EToForAccounts_2020.append(etv, ignore_index = True)
    print("ETo Done --- %s seconds ---" % (time.time() - startt))
    EToForAccounts = EToForAccounts_2020.rename(columns={"ReferenceETo(in)": "ETo(in)"})
    merged_df = pd.merge(acct_gdf, EToForAccounts, on=['Category','SubCategory','Name', 'UUID'])
    merged_df.reset_index() #reset the index to avoid any export issues
    merged_df['Category'] = 'JainlogicCustomerid_' + merged_df['customerid'].astype(str)
    merged_df['Date'] = pd.to_datetime(merged_df['Date'])

    return merged_df

def getForecastEToForDates(p,list_of_dates,df):
    N=5 #default
    #gets the weather metrics for the METRIC calcs for a point for each of the dates based
    #of the closest N stns #given the geopandas geodataframe gpd

    #Screen out the null values
    df = df[df['ReferenceETo(in)']!='--']

    cols_list = ['Date',
                 'ReferenceETo(in)'
                ]
    
    df_eto_for_dates = pd.DataFrame(columns=cols_list)
    try:
        df.Date=pd.to_datetime(df.Date)
    except:
        df.Date=pd.to_datetime(df.Date, utc=True)
        
    idxx=0
    for dt in list_of_dates:
        df_dt = df[df.Date==dt]
        stnPoints = list(pd.unique(df_dt.geometry))
        stns = pd.unique(df_dt.StationId)

        #first get the nearest stns and weights for the point provided
        [ids, pts, dis, wts] = getNearestNStnIDsInOrder_Forecast(p,stnPoints,5)

        nearestStnIDs = stns[ids]
        eto=[]
        for stn in nearestStnIDs:
            eto.append(float(df_dt[df_dt.StationId==stn]['ReferenceETo(in)'].values.tolist()[0]))
    
        weighted_eto = np.nansum(np.array(eto)*np.array(wts))

        df_eto_for_dates.loc[idxx,'Date']=dt
        df_eto_for_dates.loc[idxx,'ReferenceETo(in)']=weighted_eto
        idxx=idxx+1

    return df_eto_for_dates

def getNearestNStnIDsInOrder_Forecast(c,stnPoints,N):
    #returns the indexes, 
    #the sorted list of stns, 
    #the distances and 
    #the weights (w=1/dis)

    dist = []
    for p in stnPoints:
        d = float(c.distance(p))
        if d==0:
            d=0.000001
        else:
            dist.append(d)

    idx = np.array([i[0] for i in sorted(enumerate(dist), key=lambda x:x[1])])
    orderedStns = [stnPoints[i] for i in idx]
    orderedDist = [dist[i] for i in idx]

    M = min(N, len(stnPoints))

    orderedWghts=[]
    invDist = []
    for d in orderedDist[0:M]:
        invDist.append(1/d)
    den = sum(invDist)   
    for w in invDist:
        orderedWghts.append(w/den)

    return [idx[0:M],orderedStns[0:M],orderedDist[0:M],orderedWghts]

def generateForecastedETc2(merged_df):
    unique_accounts = pd.unique(merged_df.Category)
    print('\nRunning ETc for accts'+str(unique_accounts))
    startt = time.time()
    #Now to put this into the output folder for each account
    print("Start --- %s  ---" % datetime.now())

    for a in unique_accounts:
        #try:
        starta = time.time()

        print('\n\tWriting...%s ' % a)
        print("\tStart --- %s  ---" % datetime.now())
        acctName = a.upper()
        acct_dir = agls_acctForecastOutputPath + a + '/data/'
        #create the directory if one doesn't exist
        if not os.path.exists(acct_dir):
            os.makedirs(acct_dir)

        #For each account, take the subset of the acct_df 
        acctDF = merged_df[merged_df.Category == a]

        acctGDF = gpd.GeoDataFrame(acctDF,geometry=acctDF.geometry)

        #for every day, do the following analysis and print
        acct_dates_list = acctGDF.Date.tolist()
        unique_dts = pd.unique(acct_dates_list) # so we don't multiple count for each field
        unique_dts.sort()


        #read the ETc file for the last 7 days and file the latest date
        startDt = unique_dts[0] #get the first date as the startDt
        etc_acct_dir = agls_acctOutputPath + a + '/data/'
        #etc_acct_dir = acctOutputPath + 'JainlogicCustomerid_'+ a + '/data/'

        jsonFileList = [x for x in os.listdir(etc_acct_dir) if x.endswith('.json')]
        jsonFileList.sort()

        #print(jsonFileList)

        list_of_recent_etcfiles = [x for x in jsonFileList if datetime.strptime(x.split('.')[0].split('_')[-1], '%Y-%m-%d') >= (startDt-timedelta(days=60)) ]
        list_of_recent_etcfiles.sort()
        print('\t\nTotal files in last 60 days files = '+str(len(list_of_recent_etcfiles)))

        
        if len(list_of_recent_etcfiles)>0:
            etc_gdf = gpd.GeoDataFrame()
            for jsonFile in list_of_recent_etcfiles:
                etFile = etc_acct_dir+jsonFile
                with open(etFile, "r") as read_file:
                    data = json.load(read_file)
                    try:
                        g = gpd.GeoDataFrame.from_features(data["features"])
                        etc_gdf = etc_gdf.append(g,ignore_index=True)
                    except:
                        json_data = json.loads(data)
                        g = gpd.GeoDataFrame.from_features(json_data["features"])
                        etc_gdf = etc_gdf.append(g,ignore_index=True)
            currentETraw_df = etc_gdf[['Date','Category','SubCategory','Name','ETo(in)','ETc','ETc_low']]
            currentETraw_df.Date = pd.to_datetime(currentETraw_df.Date)
            floatCols = ['ETo(in)','ETc','ETc_low']
            currentETraw_df[floatCols] = currentETraw_df[floatCols].astype('float')
            currentETraw_df[floatCols] = currentETraw_df[floatCols].replace(0,np.nan)

            #Now interpolate for each assetid
            list_of_et_fields = pd.unique(currentETraw_df.Name).tolist()
            list_of_et_fields.sort()

            currentET_df = gpd.GeoDataFrame()
            for field in list_of_et_fields:
                gdf = currentETraw_df[currentETraw_df.Name==field]
                gdf.reset_index(drop=True, inplace=True)
                gdf.set_index('Date', inplace=True)
                #gdf = gdf.resample('D').mean()
                #print('\nBefore interpolation')
                #print(gdf)
                gdf[floatCols] = gdf[floatCols].interpolate(method='linear', limit_direction='both', axis=0)
                #print('\nAfter interpolation')
                #print(gdf)
                gdf.reset_index(inplace=True)
                gdf = gdf.sort_values(by=['Date'])
                gdf = gdf.iloc[-1:]
                #tail(1)
                gdf.reset_index(drop=True, inplace=True)
                currentET_df = currentET_df.append(gdf, ignore_index=True)

            currentET_df.dropna(subset=floatCols)
            #currentET_df = currentET_df.sort_values(by=['Date'])
            #get the last row (latest value)

            """
            latest_etc_file = list_of_recent_etcfiles[-1]
            latest_etc_dt = datetime.strptime(latest_etc_file.split('__')[1].split('_')[1].split('.')[0], '%Y-%m-%d')
            latest_etc_dtStr = datetime.strftime(latest_etc_dt, '%Y-%m-%d')

            #Read the ETc File
            etc_gdf = gpd.GeoDataFrame()
            with open(etc_acct_dir+latest_etc_file, "r") as read_file:
                data = json.load(read_file)
                try:
                    g = gpd.GeoDataFrame.from_features(data["features"])
                    etc_gdf = etc_gdf.append(g,ignore_index=True)
                except:
                    json_data = json.loads(data)
                    g = gpd.GeoDataFrame.from_features(json_data["features"])
                    etc_gdf = etc_gdf.append(g,ignore_index=True)
            currentET_df = etc_gdf[['Date','customerid', 'fieldassetid','ETo(in)','ETc','ETc_low']]
            currentET_df.Date = pd.to_datetime(currentET_df.Date)
            floatCols = ['ETo(in)','ETc','ETc_low']
            currentET_df[floatCols] = currentET_df[floatCols].astype('float')
            """
            ### DEBUG
            #print('\n\nFinal currentET_df:')
            #print(len(currentET_df))
            #print(currentET_df)
            ##############


            for dt in unique_dts:
                endDay = dt        
                startd = time.time()

                etoDateStr =  datetime.strftime(dt,'%Y-%m-%d')
                jsonFile = acct_dir + acctName+'__ForecastedET_'+etoDateStr+'.json'

                acctGDF_date = acctGDF[acctGDF.Date==dt]
                acctGDF_date.reset_index(drop=True, inplace=True)

                print('\t\t'+etoDateStr)
                for i,row in acctGDF_date.iterrows():

                    acct = row.Category
                    ranch = row.SubCategory
                    field = row.Name

                    #ETo 
                    eto = row['ETo(in)']

                    #get the last ETo, ETc and ETc_low for this assetid
                    etdf = currentET_df[(currentET_df.SubCategory==ranch) &
                                        (currentET_df.Name==field)]
                    etdf.reset_index(drop=True, inplace=True)

                    #print(etdf)


                    lastETo = etdf['ETo(in)'][0]
                    lastETc = etdf.ETc[0]
                    lastETc_low = etdf.ETc_low[0]

                    eto_ratio = eto / lastETo
                    ETc_prorated = eto_ratio * lastETc
                    ETc_low_prorated = eto_ratio * lastETc_low

                    #print([customerid, fieldassetid, eto, lastETo, ETc_prorated, ETc_low_prorated])
                    acctGDF_date.loc[i,'EToRatio']=eto_ratio
                    acctGDF_date.loc[i,'ETc']=ETc_prorated
                    acctGDF_date.loc[i,'ETc_low']=ETc_low_prorated

                """
                acctGDF_date['Date'] = acctGDF_date['Date'].astype('str')
                acctGDF_date['Acres'] = acctGDF_date['Acres'].astype('str')
                acctGDF_date['EToRatio'] = acctGDF_date['EToRatio'].astype('str')
                acctGDF_date['ETc'] = acctGDF_date['ETc'].astype('str')
                acctGDF_date['ETc_low'] = acctGDF_date['ETc_low'].astype('str')
                print(acctGDF_date.ETc.tolist())

                #save the values
                acctGDF_date['SentinelFootprint'] = acctGDF_date['SentinelFootprint'].astype('str')
                acctGDF_date['LandsatFootprint'] = acctGDF_date['LandsatFootprint'].astype('str')
                if 'CentroidPoint' in acctGDF_date.columns:
                    acctGDF_date.drop('CentroidPoint',inplace=True, axis=1)
                #print(type(acctGDF_date))
                """
                acctGDF_date['SentinelFootprint'] = acctGDF_date['SentinelFootprint'].astype('str')
                acctGDF_date['LandsatFootprint'] = acctGDF_date['LandsatFootprint'].astype('str')
                acctGDF_date['to_email_list'] = acctGDF_date['to_email_list'].astype('str')
                acctGDF_date['cc_email_list'] = acctGDF_date['cc_email_list'].astype('str')
                acctGDF_date['bcc_email_list'] = acctGDF_date['bcc_email_list'].astype('str')
                acctGDF_date['Date'] = acctGDF_date['Date'].astype('str')
                acctGDF_date.drop('CentroidPoint',inplace=True, axis=1)

                #print(type(acctGDF_date))
                print(acctGDF_date.ETc.tolist())
                
                #save the values

                #for idx, roww in acctGDF_date.iterrows():
                #	print(roww)


                acctGDF_date.to_file(jsonFile, driver = 'GeoJSON', encoding='utf-8')

            print("\tEND --- %s seconds ---" % (time.time() - starta))
            print('\tDone!')
        else:
            print('No historical data for Acct: '+str(a))
            continue
        #except:
        #	print('Exception in Acct: '+str(a))
        #	continue

    print("\nEND --- %s seconds ---" % (time.time() - startt))
    return


def generateForecastedETc(merged_df):
    unique_accounts = pd.unique(merged_df.Category)
    print('\nRunning ETc for accts'+str(unique_accounts))
    startt = time.time()
    #Now to put this into the output folder for each account
    print("Start --- %s  ---" % datetime.now())

    for a in unique_accounts:
        #try:
        starta = time.time()

        print('\n\tWriting...%s ' % a)
        print("\tStart --- %s  ---" % datetime.now())
        acctName = a.upper()
        acct_dir = acctForecastOutputPath + a + '/data/'
        #create the directory if one doesn't exist
        if not os.path.exists(acct_dir):
            os.makedirs(acct_dir)

        #For each account, take the subset of the acct_df 
        acctDF = merged_df[merged_df.Category == a]

        acctGDF = gpd.GeoDataFrame(acctDF,geometry=acctDF.geometry)

        #for every day, do the following analysis and print
        acct_dates_list = acctGDF.Date.tolist()
        unique_dts = pd.unique(acct_dates_list) # so we don't multiple count for each field
        unique_dts.sort()


        #read the ETc file for the last 7 days and file the latest date
        startDt = unique_dts[0] #get the first date as the startDt
        etc_acct_dir = acctOutputPath + a + '/data/'
        #etc_acct_dir = acctOutputPath + 'JainlogicCustomerid_'+ a + '/data/'

        jsonFileList = [x for x in os.listdir(etc_acct_dir) if x.endswith('.json')]
        jsonFileList.sort()

        #print(jsonFileList)

        list_of_recent_etcfiles = [x for x in jsonFileList if datetime.strptime(x.split('__')[1].split('_')[1].split('.')[0], '%Y-%m-%d') >= (startDt-timedelta(days=60)) ]
        list_of_recent_etcfiles.sort()
        print('\t\nTotal files in last 30 days files = '+str(len(list_of_recent_etcfiles)))

        
        if len(list_of_recent_etcfiles)>0:
            etc_gdf = gpd.GeoDataFrame()
            for jsonFile in list_of_recent_etcfiles:
                etFile = etc_acct_dir+jsonFile
                with open(etFile, "r") as read_file:
                    data = json.load(read_file)
                    try:
                        g = gpd.GeoDataFrame.from_features(data["features"])
                        etc_gdf = etc_gdf.append(g,ignore_index=True)
                    except:
                        json_data = json.loads(data)
                        g = gpd.GeoDataFrame.from_features(json_data["features"])
                        etc_gdf = etc_gdf.append(g,ignore_index=True)
            currentETraw_df = etc_gdf[['Date','customerid','fieldassetid','ETo(in)','ETc','ETc_low']]
            currentETraw_df.Date = pd.to_datetime(currentETraw_df.Date)
            floatCols = ['ETo(in)','ETc','ETc_low']
            currentETraw_df[floatCols] = currentETraw_df[floatCols].astype('float')
            currentETraw_df[floatCols] = currentETraw_df[floatCols].replace(0,np.nan)

            #Now interpolate for each assetid
            list_of_et_assetsids = pd.unique(currentETraw_df.fieldassetid).tolist()
            list_of_et_assetsids.sort()

            currentET_df = gpd.GeoDataFrame()
            for assetid in list_of_et_assetsids:
                gdf = currentETraw_df[currentETraw_df.fieldassetid==assetid]
                gdf.reset_index(drop=True, inplace=True)
                gdf.set_index('Date', inplace=True)
                #gdf = gdf.resample('D').mean()
                #print('\nBefore interpolation')
                #print(gdf)
                gdf[floatCols] = gdf[floatCols].interpolate(method='linear', limit_direction='both', axis=0)
                #print('\nAfter interpolation')
                #print(gdf)
                gdf.reset_index(inplace=True)
                gdf = gdf.sort_values(by=['Date'])
                gdf = gdf.iloc[-1:]
                #tail(1)
                gdf.reset_index(drop=True, inplace=True)
                currentET_df = currentET_df.append(gdf, ignore_index=True)

            currentET_df.dropna(subset=floatCols)
            #currentET_df = currentET_df.sort_values(by=['Date'])
            #get the last row (latest value)

            """
            latest_etc_file = list_of_recent_etcfiles[-1]
            latest_etc_dt = datetime.strptime(latest_etc_file.split('__')[1].split('_')[1].split('.')[0], '%Y-%m-%d')
            latest_etc_dtStr = datetime.strftime(latest_etc_dt, '%Y-%m-%d')

            #Read the ETc File
            etc_gdf = gpd.GeoDataFrame()
            with open(etc_acct_dir+latest_etc_file, "r") as read_file:
                data = json.load(read_file)
                try:
                    g = gpd.GeoDataFrame.from_features(data["features"])
                    etc_gdf = etc_gdf.append(g,ignore_index=True)
                except:
                    json_data = json.loads(data)
                    g = gpd.GeoDataFrame.from_features(json_data["features"])
                    etc_gdf = etc_gdf.append(g,ignore_index=True)
            currentET_df = etc_gdf[['Date','customerid', 'fieldassetid','ETo(in)','ETc','ETc_low']]
            currentET_df.Date = pd.to_datetime(currentET_df.Date)
            floatCols = ['ETo(in)','ETc','ETc_low']
            currentET_df[floatCols] = currentET_df[floatCols].astype('float')
            """
            ### DEBUG
            #print('\n\nFinal currentET_df:')
            #print(len(currentET_df))
            #print(currentET_df)
            ##############


            for dt in unique_dts:
                endDay = dt        
                startd = time.time()

                etoDateStr =  datetime.strftime(dt,'%Y-%m-%d')
                jsonFile = acct_dir + acctName+'__ForecastedET_'+etoDateStr+'.json'

                acctGDF_date = acctGDF[acctGDF.Date==dt]
                acctGDF_date.reset_index(drop=True, inplace=True)

                print('\t\t'+etoDateStr)
                for i,row in acctGDF_date.iterrows():

                    acct = row.Category
                    ranch = row.SubCategory
                    field = row.Name
                    customerid = str(row.customerid)
                    fieldassetid = str(row.fieldassetid)

                    #ETo 
                    eto = row['ETo(in)']

                    #get the last ETo, ETc and ETc_low for this assetid
                    etdf = currentET_df[(currentET_df.customerid==customerid) &
                                        (currentET_df.fieldassetid==fieldassetid)]
                    etdf.reset_index(drop=True, inplace=True)
                    lastETo = etdf['ETo(in)'][0]
                    lastETc = etdf.ETc[0]
                    lastETc_low = etdf.ETc_low[0]

                    eto_ratio = eto / lastETo
                    ETc_prorated = eto_ratio * lastETc
                    ETc_low_prorated = eto_ratio * lastETc_low

                    #print([customerid, fieldassetid, eto, lastETo, ETc_prorated, ETc_low_prorated])
                    acctGDF_date.loc[i,'EToRatio']=eto_ratio
                    acctGDF_date.loc[i,'ETc']=ETc_prorated
                    acctGDF_date.loc[i,'ETc_low']=ETc_low_prorated

                acctGDF_date['Date'] = acctGDF_date['Date'].astype('str')
                print(acctGDF_date.ETc.tolist())

                #save the values
                acctGDF_date['SentinelFootprint'] = acctGDF_date['SentinelFootprint'].astype('str')
                acctGDF_date['LandsatFootprint'] = acctGDF_date['LandsatFootprint'].astype('str')
                acctGDF_date['Date'] = acctGDF_date['Date'].astype('str')
                acctGDF_date.drop('CentroidPoint',inplace=True, axis=1)
                acctGDF_date.to_file(jsonFile, driver = 'GeoJSON')

            print("\tEND --- %s seconds ---" % (time.time() - starta))
            print('\tDone!')
        else:
            print('No historical data for Acct: '+str(a))
            continue
        #except:
        #	print('Exception in Acct: '+str(a))
        #	continue

    print("\nEND --- %s seconds ---" % (time.time() - startt))
    return

def makeTrendChart_Unif(a_gdf):
    #a_gdf is the valid dataframe for a SINGLE assetid
    #try:

    customerid = pd.unique(a_gdf.customerid).tolist()[0]
    fieldassetid = pd.unique(a_gdf.fieldassetid).tolist()[0]
    repDtStr = datetime.strftime(datetime.strptime(max(a_gdf.Date), '%Y-%m-%d %H:%M:%S'), '%Y-%m-%d')

    unif_df = a_gdf[['Date','ET_unif']]
    unif_df.Date = pd.to_datetime(unif_df.Date)
    unif_df.set_index('Date', drop=True, inplace=True)

    unif_weekly_df = unif_df.resample('W').mean()
    unif_weekly_df.reset_index(inplace=True)
    slope, intercept, r_value, p_value, std_err = stats.linregress(unif_weekly_df.index,unif_weekly_df.ET_unif)
    chgPercent = ((((((slope)*(len(unif_weekly_df)-1) + intercept) - intercept)/(intercept))*100)/(len(unif_weekly_df)-1)).round(2)
    #slopePercent = (slope*100).round(2)
    slopePercent = round(slope*100,2)
    r_value = round(r_value,2)
    std_err = round(std_err,2)
    chgPercent = round(chgPercent,2)
    titleStr = 'Percent Slope:'+str(slopePercent)+'%'+', R'+'2'.translate(trans)+' Value:'+str(r_value)+', Std. Error:'+str(std_err)+', Avg. Weekly %Change:'+str(chgPercent)+'%'
    fig = px.scatter(unif_weekly_df, x="Date", y="ET_unif", trendline="ols")
    fig.update_yaxes(title_text='<b>ETc Uniformity%</b>')
    fig.update_layout(title=titleStr)

    outFile = tmpDir+str(customerid)+'_'+str(fieldassetid)+'__UniformityTrend_'+repDtStr+'.html'
    plotly.offline.plot(fig, filename=outFile, auto_open=False)

    of = outFile
    kpi = chgPercent

    #except:
    #    print('Exception in makeTrendChart_Unif()')
    #    of = []

    return of, kpi

def makeTrendChart_Unif_PAG(a_gdf):
    #a_gdf is the valid dataframe for a SINGLE assetid
    #try:

    customerid = pd.unique(a_gdf.Category).tolist()[0]
    fieldassetid = pd.unique(a_gdf.Name).tolist()[0]

    #print(max(a_gdf.Date))
    repDtStr = datetime.strftime(datetime.strptime(max(a_gdf.Date), '%Y-%m-%d'), '%Y-%m-%d')

    unif_df = a_gdf[['Date','ET_unif']]
    unif_df.Date = pd.to_datetime(unif_df.Date)
    unif_df.set_index('Date', drop=True, inplace=True)

    unif_weekly_df = unif_df.resample('W').mean()
    unif_weekly_df.reset_index(inplace=True)
    slope, intercept, r_value, p_value, std_err = stats.linregress(unif_weekly_df.index,unif_weekly_df.ET_unif)
    chgPercent = ((((((slope)*(len(unif_weekly_df)-1) + intercept) - intercept)/(intercept))*100)/(len(unif_weekly_df)-1)).round(2)
    #slopePercent = (slope*100).round(2)
    slopePercent = round(slope*100,2)
    r_value = round(r_value,2)
    std_err = round(std_err,2)
    chgPercent = round(chgPercent,2)
    titleStr = 'Percent Slope:'+str(slopePercent)+'%'+', R'+'2'.translate(trans)+' Value:'+str(r_value)+', Std. Error:'+str(std_err)+', Avg. Weekly %Change:'+str(chgPercent)+'%'
    fig = px.scatter(unif_weekly_df, x="Date", y="ET_unif", trendline="ols")
    fig.update_yaxes(title_text='<b>ETc Uniformity%</b>')
    fig.update_layout(title=titleStr)

    outFile = tmpDir+str(customerid)+'_'+str(fieldassetid)+'__UniformityTrend_'+repDtStr+'.html'
    plotly.offline.plot(fig, filename=outFile, auto_open=False)

    of = outFile
    kpi = chgPercent

    #except:
    #    print('Exception in makeTrendChart_Unif()')
    #    of = []

    return of, kpi

"""
def makeTrendChart_Unif(a_gdf):
    #a_gdf is the valid dataframe for a SINGLE assetid
    #try:

    print(a_gdf.end_ts)

    customerid = pd.unique(a_gdf.customerid).tolist()[0]
    fieldassetid = pd.unique(a_gdf.fieldassetid).tolist()[0]
    repDtStr = datetime.strftime(datetime.strptime(max(a_gdf.end_ts), '%Y-%m-%dT%H:%M:%SZ'), '%Y-%m-%d')

    unif_df = a_gdf[['end_ts','etunif']]
    unif_df.rename(columns={'end_ts':'Date','etunif':'ET_unif'}, inplace=True)
    unif_df.Date = pd.to_datetime(unif_df.Date)
    unif_df.ET_unif = unif_df.ET_unif.astype('float')
    unif_df.set_index('Date', drop=True, inplace=True)
    unif_weekly_df = unif_df.resample('W').mean()
    unif_weekly_df.reset_index(inplace=True)

    #unif_df = a_gdf[['Date','ET_unif']]
    #unif_df.Date = pd.to_datetime(unif_df.Date)
    #unif_df.set_index('Date', drop=True, inplace=True)
    #unif_weekly_df = unif_df.resample('W').mean()

    
    slope, intercept, r_value, p_value, std_err = stats.linregress(unif_weekly_df.index,unif_weekly_df.ET_unif)
    chgPercent = ((((((slope)*(len(unif_weekly_df)-1) + intercept) - intercept)/(intercept))*100)/(len(unif_weekly_df)-1)).round(2)

    titleStr = 'Percent Slope:'+str(slopePercent)+'%'+', R'+'2'.translate(trans)+' Value:'+str((r_value).round(2))+', Std. Error:'+str(std_err.round(2))+', Avg. Weekly %Change:'+str(chgPercent)+'%'
    fig = px.scatter(unif_weekly_df, x="Date", y="ET_unif", trendline="ols")
    fig.update_layout(title=titleStr)
    
    outFile = tmpDir+str(customerid)+'_'+str(fieldassetid)+'__UniformityTrend_'+repDtStr+'.html'
    plotly.offline.plot(fig, filename=outFile)

    of = outFile
    kpi = chgPercent
    
    #except:
    #    print('Exception in makeFieldMap()')
    #    of = []
   
    return of, kpi
"""
def makeTrendChart_et(at_gdf, region):
    #a_gdf is the valid dataframe for a SINGLE assetid
    #try:
    
    if region=='AUS':
        units = ' (mm) '
        at_gdf.ETc_average = at_gdf.ETc_average * 25.4
    else:
        units = ' (inches) '

    customerid = pd.unique(at_gdf.customerid).tolist()[0]
    fieldassetid = pd.unique(at_gdf.fieldassetid).tolist()[0]
    
    customername = pd.unique(at_gdf.customername).tolist()[0]
    fieldname = pd.unique(at_gdf.fieldname).tolist()[0]
    

    et_df = at_gdf[['Date','ETc_average']]
    et_df.Date = pd.to_datetime(et_df.Date)
    et_df.ETc_average = et_df.ETc_average.astype('float')
    repDtStr = datetime.strftime(max(et_df.Date),'%Y-%m-%d')
    et_df.set_index('Date', drop=True, inplace=True)
    
    ########### Interpolate Missing Values
    et_df = et_df.resample('D').interpolate()
    ######################################

    et_weekly_df = et_df.resample('W').sum()
    et_weekly_df.reset_index(inplace=True)
    et_weekly_df['ETc_cum'] = et_weekly_df['ETc_average'].cumsum()
    
    
    fig = make_subplots(specs=[[{"secondary_y": True}]])
    
    y1_titleStr = 'ETc (average)'+units
    y2_titleStr = 'ETc (cumulative)'+units
    
    meanETc = go.Scatter(
        name=y1_titleStr,
        x=et_weekly_df.Date,
        y=et_weekly_df.ETc_average,
        mode='lines+markers',
        line=dict(color='green'),
    )
    cumETc = go.Bar(
        name=y2_titleStr,
        x=et_weekly_df.Date,
        y=et_weekly_df.ETc_cum,
        opacity=0.3,
        marker_color='blue'
    )
    
    fig.add_trace(meanETc,secondary_y=False)
    fig.add_trace(cumETc,secondary_y=True)

    fig.update_yaxes(title_text='<b>'+y1_titleStr+'</b>', secondary_y=False)
    fig.update_yaxes(title_text='<b>'+y2_titleStr+'</b>', secondary_y=True)
    fig.update_layout(title='Weekly Average and Cumulative ETc Trend For Field: '+str(fieldname))

    fig.update_layout(
    yaxis=dict(
        titlefont=dict(color="green"),
        tick0=0
    ),
    yaxis2=dict(
        titlefont=dict(color="blue"),
        tick0=0
    ))
    """
    cumETc = go.Scatter(
        name=y2_titleStr,
        x=et_weekly_df.Date,
        y=et_weekly_df.ETc_cum,
        mode='lines+markers',
        line=dict(color='black',dash='dashdot'),
    )
    
    fig.add_trace(meanETc,secondary_y=False)
    fig.add_trace(cumETc,secondary_y=True)

    fig.update_yaxes(title_text='<b>'+y1_titleStr+'</b>', secondary_y=False)
    fig.update_yaxes(title_text='<b>'+y2_titleStr+'</b>', secondary_y=True)
    fig.update_layout(title='Weekly Average and Cumulative ETc Trend For Field: '+str(fieldname))

    fig.update_layout(
    yaxis=dict(
        titlefont=dict(color="green"),
        tick0=0
    ),
    yaxis2=dict(
        titlefont=dict(color="black"),
        tick0=0
    ))
    """
        
    outFile = tmpDir+str(customerid)+'_'+str(fieldassetid)+'__ETcTrend_'+repDtStr+'.html'
    plotly.offline.plot(fig, filename=outFile, auto_open=False)

    of = outFile
                  
    #except:
    #    print('Exception in makeTrendChart_et()')
    #    of = []

    return of


def makeTrendChart_et_PAG(at_gdf, region):
    #a_gdf is the valid dataframe for a SINGLE assetid
    #try:
    
    if region=='AUS':
        units = ' (mm) '
        at_gdf.ETc_average = at_gdf.ETc_average * 25.4
    else:
        units = ' (inches) '

    #print(at_gdf)

    customerid = pd.unique(at_gdf.Category).tolist()[0]
    fieldassetid = pd.unique(at_gdf.Name).tolist()[0]
    
    customername = pd.unique(at_gdf.Category).tolist()[0]
    fieldname = pd.unique(at_gdf.Name).tolist()[0]
    

    et_df = at_gdf[['Date','ETc_average']]
    et_df.Date = pd.to_datetime(et_df.Date)
    et_df.ETc_average = et_df.ETc_average.astype('float')
    repDtStr = datetime.strftime(max(et_df.Date),'%Y-%m-%d')
    et_df.set_index('Date', drop=True, inplace=True)
    
    ########### Interpolate Missing Values
    et_df = et_df.resample('D').interpolate()
    ######################################

    et_weekly_df = et_df.resample('W').sum()
    et_weekly_df.reset_index(inplace=True)
    et_weekly_df['ETc_cum'] = et_weekly_df['ETc_average'].cumsum()
    
    
    fig = make_subplots(specs=[[{"secondary_y": True}]])
    
    y1_titleStr = 'ETc (average)'+units
    y2_titleStr = 'ETc (cumulative)'+units
    
    meanETc = go.Scatter(
        name=y1_titleStr,
        x=et_weekly_df.Date,
        y=et_weekly_df.ETc_average,
        mode='lines+markers',
        line=dict(color='green'),
    )
    cumETc = go.Bar(
        name=y2_titleStr,
        x=et_weekly_df.Date,
        y=et_weekly_df.ETc_cum,
        opacity=0.3,
        marker_color='blue'
    )
    
    fig.add_trace(meanETc,secondary_y=False)
    fig.add_trace(cumETc,secondary_y=True)

    fig.update_yaxes(title_text='<b>'+y1_titleStr+'</b>', secondary_y=False)
    fig.update_yaxes(title_text='<b>'+y2_titleStr+'</b>', secondary_y=True)
    fig.update_layout(title='Weekly Average and Cumulative ETc Trend For Field: '+str(fieldname))

    fig.update_layout(
    yaxis=dict(
        titlefont=dict(color="green"),
        tick0=0
    ),
    yaxis2=dict(
        titlefont=dict(color="blue"),
        tick0=0
    ))
    """
    cumETc = go.Scatter(
        name=y2_titleStr,
        x=et_weekly_df.Date,
        y=et_weekly_df.ETc_cum,
        mode='lines+markers',
        line=dict(color='black',dash='dashdot'),
    )
    
    fig.add_trace(meanETc,secondary_y=False)
    fig.add_trace(cumETc,secondary_y=True)

    fig.update_yaxes(title_text='<b>'+y1_titleStr+'</b>', secondary_y=False)
    fig.update_yaxes(title_text='<b>'+y2_titleStr+'</b>', secondary_y=True)
    fig.update_layout(title='Weekly Average and Cumulative ETc Trend For Field: '+str(fieldname))

    fig.update_layout(
    yaxis=dict(
        titlefont=dict(color="green"),
        tick0=0
    ),
    yaxis2=dict(
        titlefont=dict(color="black"),
        tick0=0
    ))
    """
        
    outFile = tmpDir+str(customerid)+'_'+str(fieldassetid)+'__ETcTrend_'+repDtStr+'.html'
    plotly.offline.plot(fig, filename=outFile, auto_open=False)

    of = outFile
                  
    #except:
    #    print('Exception in makeTrendChart_et()')
    #    of = []

    return of

def getHTMLLinkToChart(htmlFile, kpi):
    
    file_path = Path(htmlFile)

    if float(kpi) < 0:
        stl = ' style=\"color:#FF0000;\" '
        sn = ''
    else:
        stl = ' style=\"color:#000000;\" '
        sn='+'

    linkURL='https://www.agralogics.com/HyperGrow/'+file_path.name
    linkStr = '<a href=\"'+linkURL+'\"'+stl+'>'+str(sn)+str(kpi)+'%</a>'  
    
    return linkStr


def getHTMLLinkToETcChart(htmlFile, kpi):
    
    file_path = Path(htmlFile)

    if float(kpi) < 0:
        stl = ' style=\"color:#FF0000;\" '
    else:
        stl = ' style=\"color:#000000;\" '

    linkURL='https://www.agralogics.com/HyperGrow/'+file_path.name
    linkStr = '<a href=\"'+linkURL+'\"'+stl+'>'+str(kpi)+'</a>'  
    
    return linkStr

def writeActiveFilesToJSONConfigFiles(gdf_AllJLAccts):

    activeFields_df = gdf_AllJLAccts[['customerid','customername','fieldassetid',
           'fieldname', 'monitored', 'subscribedhypergrow', 'subscribedhyperview']]
    list_of_HG_activeCustomerids = pd.unique(activeFields_df[activeFields_df.subscribedhypergrow==True].customerid).tolist()
    list_of_HG_activeCustomerids.sort()

    list_of_HV_activeCustomerids = pd.unique(activeFields_df[activeFields_df.subscribedhyperview==True].customerid).tolist()
    list_of_HV_activeCustomerids.sort()

    for customerid in list_of_HG_activeCustomerids:
        activeJsonFile = acctConfigPath + 'JainlogicCustomerid_'+str(customerid)+'/JainlogicCustomerid_'+str(customerid)+'__HG_activeFieldsConfig.json'

        df = activeFields_df[activeFields_df.customerid==customerid]
        dff = df[df.subscribedhypergrow==True]
        dff.reset_index(drop=True, inplace=True)
        #Write to the json file if the dff is > 0
        if len(dff)>0:
            #Make the directory of it doesn't exists
            if not os.path.exists(acctConfigPath + 'JainlogicCustomerid_'+str(customerid)):
                os.makedirs(acctConfigPath + 'JainlogicCustomerid_'+str(customerid));
            dff.to_json(activeJsonFile)
        else:
            print('\nNo active HG fields in customerid = '+str(customerid))

        
    for customerid in list_of_HV_activeCustomerids:
        activeJsonFile = acctConfigPath + 'JainlogicCustomerid_'+str(customerid)+'/JainlogicCustomerid_'+str(customerid)+'__HV_activeFieldsConfig.json'

        df = activeFields_df[activeFields_df.customerid==customerid]
        dff = df[df.subscribedhyperview==True]
        dff.reset_index(drop=True, inplace=True)

        #Write to the json file if the dff is > 0
        if len(dff)>0:
            #Make the directory of it doesn't exists
            if not os.path.exists(acctConfigPath + 'JainlogicCustomerid_'+str(customerid)):
                os.makedirs(acctConfigPath + 'JainlogicCustomerid_'+str(customerid));
            dff.to_json(activeJsonFile)
        else:
            print('\nNo active HV fields in customerid = '+str(customerid))

    return

def normalizeCustomerName(acct_gdf):

    normalized_acct_gdf = acct_gdf.copy()

    if "customername" in normalized_acct_gdf.columns:
        list_of_customernames = pd.unique(normalized_acct_gdf.customername)
        num_of_customernames = len(list_of_customernames)
        max_len=0;
        cname=''
        if num_of_customernames >= 2:
            for c in list_of_customernames:
                a_gdf = normalized_acct_gdf[normalized_acct_gdf.customername==c]
                if len(a_gdf)>=max_len:
                    max_len=len(a_gdf)
                    cname = c
                else:
                    continue
            normalized_acct_gdf.customername=cname
    elif "Ranch" in normalized_acct_gdf.columns:
        list_of_customernames = pd.unique(normalized_acct_gdf.Ranch)
        num_of_customernames = len(list_of_customernames)
        max_len=0;
        cname=''
        if num_of_customernames >= 2:
            for c in list_of_customernames:
                a_gdf = normalized_acct_gdf[normalized_acct_gdf.Ranch==c]
                if len(a_gdf)>=max_len:
                    max_len=len(a_gdf)
                    cname = c
                else:
                    continue
            normalized_acct_gdf.Ranch=cname

    return normalized_acct_gdf


def interpolateAndOverWriteNaNs(a_gdf, startDtStr, endDtStr):
    #returns a corrected DataFrame if NaNs are detected in any of the fields data
    #returns empty dataframe if no NaNs are detected

    list_of_accts = pd.unique(a_gdf.customerid).tolist()
    list_of_accts.sort()
    acct = list_of_accts[0] # Just do the first acct in this dataframe

    acct = str(acct)
    flag=0

    acct_gdf = a_gdf.copy()

    #Now read the ET Files
    startDt = datetime.strptime(startDtStr, '%Y-%m-%d')
    endDt = datetime.strptime(endDtStr, '%Y-%m-%d')
    etDataDir = acctOutputPath + 'JainlogicCustomerid_'+acct+'/data/'

    acctData_gdf = gpd.GeoDataFrame()
    list_of_etFiles=[]
    list_of_etFilesWithoutRoot = []
    list_of_dtStrs = [datetime.strftime(x, '%Y-%m-%d') for x in pd.date_range(startDt,endDt-timedelta(days=1),freq='d')]
    for dtStr in list_of_dtStrs:
        etFile = etDataDir + 'JAINLOGICCUSTOMERID_'+acct+'__ET_'+dtStr+'.json'
        list_of_etFiles.append(etFile)
        gdf = gpd.read_file(etFile)
        acctData_gdf = acctData_gdf.append(gdf,ignore_index=True)


    list_of_fids = pd.unique(acctData_gdf.fieldassetid)
    list_of_fids.sort()

    new_acctData_gdf = gpd.GeoDataFrame()

    for fid in list_of_fids:
        a_gdf = acctData_gdf[acctData_gdf.fieldassetid==fid]
        a_gdf.reset_index(drop=True, inplace=True)
        keepCols=['Date','ETc', 'ETc_low', 'Kc', 'ETm', 'ET_p90', 'ET_unif']
        asub_gdf = a_gdf[keepCols]

        # Print the number of NaN entries
        nan_count = asub_gdf['ETc'].isna().sum()
        if nan_count>0:
            print('\n\tFound %s NaNs in Acct %s, field: %s' % (str(nan_count), str(acct), str(fid)))
            asub_gdf.Date = pd.to_datetime(asub_gdf.Date)
            asub_gdf = asub_gdf.set_index('Date')
            aint_gdf = asub_gdf.interpolate(method='linear', limit_direction='both', axis=0)
            #resample('D').interpolate(method='linear', limit_direction='both', axis=0)
            aint_gdf.reset_index(inplace=True)
            a_gdf[keepCols]=aint_gdf.copy()
            print('\tFixed. New DF has %s NaNs' % (a_gdf['ETc'].isna().sum()))
            flag=1
        new_acctData_gdf = new_acctData_gdf.append(a_gdf, ignore_index=True)

    if flag==1:
        return new_acctData_gdf, list_of_etFiles
    else:
        return pd.DataFrame(), []
def getEToForStnList(obs_latest_gdf, list_of_stns):

    allETo_gdf = gpd.GeoDataFrame()
    minElements = ['TT','TN1','TX1']
    for stn in list_of_stns:
        obs_gdf = obs_latest_gdf[obs_latest_gdf.StationNumber==stn]
        obs_gdf.reset_index(drop=True, inplace=True)
        obs_gdf.ParameterValue=obs_gdf.ParameterValue.astype('float')
        list_of_elements = pd.unique(obs_gdf.ParameterShortName).tolist()
    
        if set(minElements).issubset(list_of_elements): #make sure they are in there
            ######### Process the ETo
            obs_gdf = obs_gdf.dropna()
            obs_gdf.reset_index(drop=True, inplace=True)
            et1 = ETo()
            freq = 'D'
            lat=float(obs_gdf.StationLatitude[0])
            lon=float(obs_gdf.StationLongitude[0])
            z_msl=float(float(obs_gdf.StationHeight[0]))

            #Pivot for a flat Datastructure
            allCols = ['Date', 'StationName', 'StationNumber','ParameterShortName','ParameterFullName', 'ParameterValue', 'ParameterUnits']
            obs_df = obs_gdf[allCols]
            obs_df['Date'] = pd.to_datetime(obs_df['Date']).apply(lambda x: x.strftime('%Y-%m-%dT%H%M%S'))
            print('\tComputing ETo for stn '+str(stn))
            tt = pd.pivot_table(obs_df, values='ParameterValue', index=['Date'], columns=['ParameterShortName'], aggfunc=np.mean)
            tt.index = pd.to_datetime(tt.index)
            o_df = tt.resample('H').interpolate(method='linear', limit_direction='both')

            colsToDivideBy10 = tt.columns.tolist()
            colsToDivideBy10.remove('RH')
            o_df[colsToDivideBy10] = o_df[colsToDivideBy10]/10
            
            
            o_df['StnID']=obs_gdf['StationNumber'][0]
            o_df['lat']=lat
            o_df['lon']=lon
            o_df['elev']=z_msl
            o_df['desc']=obs_gdf['StationName'][0]
            o_df['geometry']=obs_gdf.geometry[0]
            
            set0 = ['TT','TN1','TX1','RH','FF','TD', 'P0']
            set1 = ['TT','TN1','TX1','RH','FF','TD']
            set2 = ['TT','TN1','TX1','RH','TD']
            set3 = ['TT','TN1','TX1']

            if all(elem in o_df.columns.tolist()  for elem in set0):
                p_df = o_df[set0]
                p_df.reset_index(drop=False, inplace=True)
                p_df['P0']=p_df['P0']/10
                
                eto_cols = ['TS','T_mean','T_min','T_max', 'RH_mean','U_z','T_dew', 'P']
                p_df.columns=eto_cols
                p_df.set_index('TS',drop=True, inplace=True)
                p_df.index.name='date'
                p_df = p_df.astype('float')

            if all(elem in o_df.columns.tolist()  for elem in set1):
                p_df = o_df[set1]
                p_df.reset_index(drop=False, inplace=True)

                eto_cols = ['TS','T_mean','T_min','T_max', 'RH_mean','U_z','T_dew']
                p_df.columns=eto_cols
                p_df.set_index('TS',drop=True, inplace=True)
                p_df.index.name='date'
                p_df = p_df.astype('float')

            elif all(elem in o_df.columns.tolist()  for elem in set2):
                p_df = o_df[set2]
                p_df.reset_index(drop=False, inplace=True)

                eto_cols = ['TS','T_mean','T_min','T_max', 'RH_mean','T_dew']
                p_df.columns=eto_cols
                p_df.set_index('TS',drop=True, inplace=True)
                p_df.index.name='date'
                p_df = p_df.astype('float')
            elif all(elem in o_df.columns.tolist()  for elem in set3):
                p_df = o_df[set3]
                p_df.reset_index(drop=False, inplace=True)

                eto_cols = ['TS','T_mean','T_min','T_max']
                p_df.columns=eto_cols
                p_df.set_index('TS',drop=True, inplace=True)
                p_df.index.name='date'
                p_df = p_df.astype('float')

            else:
                break
            
            et1.param_est(p_df,freq,z_msl,lat,lon)
            res_df = et1.ts_param
            Rn=float(res_df.R_n[0]) #Net radiation (MJ/m2)
            Rs=float(res_df.R_s[0]) #Incoming shortwave radiation (MJ/m2)

            eto1 = et1.eto_fao()
            eto2 = et1.eto_hargreaves()

            o_df['ETo_FAO_MM'] = eto1
            o_df['ETo_HAR_MM'] = eto2
            o_df['ETo_AVG_MM'] = (eto1 + eto2)/2
            
            #print([eto1, eto2])
            
            # set other parameters
            
            utc_time_list = [x for x in eto1.index.tolist()]
            local_time_list = getLocalTime(utc_time_list, lon)
            aifstime_utc = [datetime.strftime(x, '%Y%m%d%H%M%S') for x in utc_time_list]
            local_date_time_full = [datetime.strftime(x, '%Y%m%d%H%M%S') for x in local_time_list]
            o_df["aifstime_utc"] = aifstime_utc
            o_df['local_date_time_full'] = local_date_time_full
            
        
            o_gdf = gpd.GeoDataFrame(o_df, geometry=o_df.geometry)

            o_gdf.reset_index(inplace=True, drop=True)
            
            allETo_gdf = allETo_gdf.append(o_gdf, ignore_index=True)
            
    return allETo_gdf

def getLocalTime(utc, longitude):
    """
    :param utc list: string Ex. '2008-12-2'
    :param longitude: longitude
    :return: Local Mean Time Timestamp
    """
    lmt = [x + timedelta(seconds=round(4*60*longitude)) for x in utc]
    lmt = [x.replace(tzinfo=None) for x in lmt]
    return lmt   




def type_py2sql(pytype):
    '''Return the closest sql type for a given python type'''
    if pytype in _type_py2sql_dict:
        return _type_py2sql_dict[pytype]
    else:
        raise NotImplementedError(
            "You may add custom `sqltype` to `"+str(pytype)+"` assignment in `_type_py2sql_dict`.")


def insertDF(gdf, dbname):
    #connect to the PostgreSQL DB
    hostname = 'localhost'
    username = 'agls_admin'
    passwd = 'Pass1234'
    database = dbname
    #'agls'
    port = '5432'
    host='127.0.0.1'

    #2. Create a table and insert into that table
    #TO PUT DATA
    cmdStr = 'postgresql+psycopg2://'+username+':'+passwd+'@'+hostname+':'+port+'/'+database
    engine = create_engine(cmdStr)
    conn = engine.raw_connection()
    cur = conn.cursor()

    writeType = 'replace'
    indexType=False
    geomType='POLYGON'
    geomCol = 'geom'

    # make the insert
    gdf.to_sql(name=dbname, 
                    con=engine, if_exists=writeType, index=indexType,
                    dtype={geomCol: Geometry(geomType, srid= 4326)})
    conn.close()

    return
def parse_point(record):
    pieces = record.split() # splits each record into a list of 3
    x = float(pieces[2].lstrip('(')) # latitude
    y = float(pieces[3].rstrip('')) # longitude 
    z = float(pieces[4].rstrip(')')) # elevation 
    point = Point(x,y,z) # convert to Shapely Point
    return point

def unparse_point(record):
    pieces = record.split() # splits each record into a list of 3
    x = float(pieces[2].lstrip('(')) # latitude
    y = float(pieces[3].rstrip('')) # longitude 
    z = float(pieces[4].rstrip(')')) # elevation 
    point = [x,y,z] # convert to coordinates 
    return point

def getNDates(N):
    startDtStr = datetime.strftime(datetime.now()-timedelta(days=N),'%Y-%m-%dT00:00:00+00:00')
    endDtStr = datetime.strftime(datetime.now()-timedelta(days=0), '%Y-%m-%dT00:00:00+00:00')

    #Put the URL replacements
    startDtStrN = startDtStr.replace(':','%3A').replace('+','%2B')
    endDtStrN = endDtStr.replace(':','%3A').replace('+','%2B')
    
    return [startDtStrN,endDtStrN]

def getDates():
    startDtStr = datetime.strftime(datetime.now()-timedelta(days=1),'%Y-%m-%dT00:00:00+00:00')
    endDtStr = datetime.strftime(datetime.now()-timedelta(days=0), '%Y-%m-%dT00:00:00+00:00')

    #Put the URL replacements
    startDtStrN = startDtStr.replace(':','%3A').replace('+','%2B')
    endDtStrN = endDtStr.replace(':','%3A').replace('+','%2B')
    
    return [startDtStrN,endDtStrN]

def getWeatherCols(obsData_gdf, weather_cols):
    for i,row in obsData_gdf.iterrows():
        for metric in weather_cols:
            colDict = row[metric]
            if type(colDict)==float:
                v=np.nan
                continue
            else:
                val = colDict['value']
                if val:
                    v = val
                else:
                    v = np.nan          

                colHeader = metric+'__'+colDict['unitCode'].split(':')[-1]
                obsData_gdf.loc[i,colHeader] = v
                
    return obsData_gdf

def getIFRAMELinkTo(htmlFile,metricStr):
    file_path = Path(htmlFile)
    linkURL='https://www.agralogics.com/HyperGrow/'+file_path.name
    return linkURL

def getIFRAMELinkToWithFolder(folder, htmlFile,metricStr):
    file_path = Path(htmlFile)
    linkURL='https://www.agralogics.com/'+folder+'/'+file_path.name
    return linkURL


def fixZinPolygon(poly):
    coordsList = [[poly.exterior.xy[0][x],poly.exterior.xy[1][x]] for x in range(0,len(poly.exterior.xy[0]))]
    new_poly = shapely.geometry.Polygon(coordsList)
    return new_poly


def logDtMISMatchErrors(acct, startDtStr, endDtStr, listOfJsonFiles, listOfJsonDts, errLogFile):
    
    startDt = datetime.strptime(startDtStr,'%Y-%m-%d')
    endDt = datetime.strptime(endDtStr,'%Y-%m-%d')
    
    errStr = '\n\n'+str(datetime.now())
    if ((min(listOfJsonDts) != startDt) | (max(listOfJsonDts) !=endDt) | (len(listOfJsonFiles) != 7)):
        errStr += '\nET Error in Acct: '+acct
        if (min(listOfJsonDts) != startDt):
            errStr += '\nReport Start Date and ET JSON Start Date are not identifcal !!'
            errStr += '\nReport Start Date: '+startDtStr+', ET JSON start date: '+datetime.strftime(min(listOfJsonDts), '%Y-%m-%d')
        elif (max(listOfJsonDts) != endDt):
            errStr += '\nReport End Date and ET JSON End Date are not identifcal !!'
            errStr += '\nReport End Date: '+endDtStr+', ET JSON End Date: '+datetime.strftime(max(listOfJsonDts), '%Y-%m-%d')
        elif len(listOfJsonFiles) != 7:
            errStr += '\nnot enough Jsons for a whole week'
            errStr += '\nNo of Json Files: '+str(len(listOfJsonFiles))
    else:
        errStr += '\nET Data Dt Check for Acct '+str(a)+'for Period=['+startDtStr+' - '+endDtStr+'] ...OK'
        
        
    with open(errLogFile, 'a') as errLog:
            errLog.write(errStr)
            
    return

def get_s2_cellids_and_token_list(resLevel, lats, lons):
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

def getEToForArchiveObs(region_gdf):
    
    #get a list of all stns
    list_of_stns = pd.unique(region_gdf.STATION).tolist()
    list_of_stns.sort()
    
    allETo_gdf = gpd.GeoDataFrame()
    
    for StnID in list_of_stns:
        print(StnID)
        obs_gdf = region_gdf[region_gdf.STATION==StnID]
        obs_gdf.reset_index(drop=True, inplace=True)
        
        obs_df = pd.DataFrame(obs_gdf)
        
        #Select the columns
        keepCols = ['DATE', 'TEMP', 'MIN', 'MAX', 'DEWP', 'STP', 'WDSP']
        o_df = obs_df[keepCols]
        
        #Rename the columns
        o_df = o_df.rename(columns={
            'DATE': 'date',
            'TEMP': 'T_mean',
            'MIN': 'T_min',
            'MAX': 'T_max',
            'DEWP': 'T_dew',
            'STP': 'P',
            'WDSP': 'U_z',
        })
        #Correct the units
        o_df['T_mean'] = (5/9) * (o_df['T_mean'] - 32) 
        o_df['T_min'] = (5/9) * (o_df['T_min'] - 32)
        o_df['T_max'] = (5/9) * (o_df['T_max'] - 32)
        o_df['T_dew'] = (5/9) * (o_df['T_dew'] - 32)
        o_df['P'] = 0.1 * o_df['P'] #convert mbar to kpa
        o_df['U_z'] = 0.514444 * o_df['U_z'] #convert to mps
        
        #set the date
        o_df.set_index('date',drop=True, inplace=True)
        o_df = o_df.astype('float')
            
        #setup the array
        et1 = ETo()
        freq = 'D'
        lat = obs_df.LATITUDE[0]
        lon = obs_df.LONGITUDE[0]
        z_msl = obs_df.ELEVATION[0]
        geom = obs_df.geometry[0]
        
        
        et1.param_est(o_df,freq,z_msl,lat,lon)
        res_df = et1.ts_param
        Rn=float(res_df.R_n[0]) #Net radiation (MJ/m2)
        Rs=float(res_df.R_s[0]) #Incoming shortwave radiation (MJ/m2)

        eto1 = et1.eto_fao()
        eto2 = et1.eto_hargreaves()

        o_df['ETo_FAO_MM'] = eto1
        o_df['ETo_HAR_MM'] = eto2
        o_df['ETo_AVG_MM'] = (eto1 + eto2)/2


        o_df['ETo_FAO_IN'] = eto1 / 25.4
        o_df['ETo_HAR_IN'] = eto2 / 25.4
        o_df['ETo_AVG_IN'] = ((eto1 + eto2)/2) / 25.4
        
        o_df['geometry']=geom
        o_gdf = gpd.GeoDataFrame(o_df, geometry=o_df.geometry)

        o_gdf.reset_index(inplace=True)
        allETo_gdf = allETo_gdf.append(o_gdf, ignore_index=True)
        
        
    return allETo_gdf

import time

def getEtoFromNCEP(lat, lon, filePath, start_date, end_date):
    lats=[lat]
    lons=[lon]
    start_time = time.time()
    tok = start_date.split('-')
    YYYY_str = tok[0]
    MM_str = tok[1]
    DD_str = tok[2]
    
    if end_date:
        # Parse end date and localize it to UTC if provided
        end_tok = end_date.split('-')
        end_YYYY_str = int(end_tok[0])
        end_MM_str = int(end_tok[1])
        end_DD_str = int(end_tok[2])
        end_local_dt = datetime(end_YYYY_str, end_MM_str, end_DD_str)
        utc_end_dt = pytz.utc.localize(end_local_dt)
    
    
    #get the list of S2 indeces and CIDs for the data point
    s2_index__L5_list, L5_cids = get_s2_cellids_and_token_list(5, lats, lons)
    # s2_index__L7_list, L7_cids = get_s2_cellids_and_token_list(7, lats, lons)
    # s2_index__L9_list, L9_cids = get_s2_cellids_and_token_list(9, lats, lons)
    
    # print(s2_index__L5_list)
    # print(s2_index__L7_list)
    # print(s2_index__L9_list)

    list_of_L5_paths = []
    list_of_L7_paths = []
    list_of_L9_paths = []
    
    list_of_L5_paths = [filePath+'s2_token_L5='+x for x in s2_index__L5_list 
                        if os.path.exists(filePath+'s2_token_L5='+x)]
    
    weather_datasets = []
        
    for x in list_of_L5_paths:
        weather_datasets.append(ds.dataset(x, format="parquet", partitioning="hive"))
    
    #get the Data
    w_all = pd.DataFrame()
    for weatherDataset in weather_datasets:
        """
        dataPath = filePath+'s2_index__L3='+s2_index__L3_list[0]+'/s2_index__L5='+s2_index__L5_list[0]+'/s2_index__L8='+s2_index__L8_list[0] 
        if os.path.exists(dataPath):
            dataDir = dataPath
        else:
            dataDir = []
            return []
        
        #get the parquet dataset from the directory
        weatherDataset = ds.dataset(dataDir,format="parquet", partitioning="hive")
        """
        
        dt = datetime.strptime(YYYY_str+'-'+MM_str+'-'+DD_str, '%Y-%m-%d')
        try:
            
            weather_df = weatherDataset.to_table(
                columns=['Year','Month','Day','ETo_AVG_IN'],
                filter=(
                    (ds.field('Year') >= int(YYYY_str)) & (ds.field('Year') <= end_YYYY_str) &
                    (ds.field('Month') >= int(MM_str)) & (ds.field('Month') <= end_MM_str) &
                    (ds.field('Day') >= int(DD_str)) & (ds.field('Day') <= end_DD_str)
                )
            ).to_pandas()

            #get the DataFrame for the average of that day
            weather_df['YYYY']=weather_df['Year'].astype(str)
            weather_df['MM']=weather_df['Month'].astype(str).str.zfill(2)
            weather_df['DD']=weather_df['Day'].astype(str).str.zfill(2)

            weather_df['Date']=weather_df.YYYY+'-'+weather_df.MM+'-'+weather_df.DD
            weather_df.Date = pd.to_datetime(weather_df.Date)

            weather_df = weather_df.drop(columns=['YYYY','MM','DD'], axis=1)
            w_df = pd.DataFrame(weather_df.mean(axis=0)).T
            w_df['Date']=dt

            cols = ['Date','ETo_AVG_IN']
            w_df = w_df[cols]
            w_all = pd.concat([w_all, w_df], ignore_index=True)
            
        except Exception as e:
            print(e)
    
    if len(w_all) > 0:
        w_all['ETo_AVG_IN'] = pd.to_numeric(w_all['ETo_AVG_IN'], errors='coerce')
        w_all = w_all.dropna(subset=['ETo_AVG_IN'])  # Drop rows with NaN in ETo_AVG_IN

        if not w_all.empty:
            # Take the average of the columns
            w_ret = pd.DataFrame(w_all.groupby(["Date"])["ETo_AVG_IN"].mean())
            w_ret.reset_index(inplace=True)
        else:
            w_ret = pd.DataFrame()
    else:
        w_ret = pd.DataFrame()

    
    end_time = time.time()
    time_elapsed = (end_time - start_time)
    
    return  w_ret

def getEtoFromNLDAS(lat, lon, filePath, start_date, end_date=None):
    lats=[lat]
    lons=[lon]
    start_time = time.time()
    tok = start_date.split('-')
    YYYY_str = tok[0]
    MM_str = tok[1]
    DD_str = tok[2]
    
    if end_date:
        # Parse end date and localize it to UTC if provided
        end_tok = end_date.split('-')
        end_YYYY_str = int(end_tok[0])
        end_MM_str = int(end_tok[1])
        end_DD_str = int(end_tok[2])
        end_local_dt = datetime(end_YYYY_str, end_MM_str, end_DD_str)
        utc_end_dt = pytz.utc.localize(end_local_dt)
    
    
    #get the list of S2 indeces and CIDs for the data point
    s2_index__L5_list, L5_cids = get_s2_cellids_and_token_list(5, lats, lons)
    list_of_L5_paths = []

    list_of_L5_paths = [filePath+'s2_token_L5='+x for x in s2_index__L5_list 
                        if os.path.exists(filePath+'s2_token_L5='+x)]
    
    weather_datasets = []
    for x in list_of_L5_paths:
        weather_datasets.append(ds.dataset(x, format="parquet", partitioning="hive"))
    
    #get the Data
    w_all = pd.DataFrame()
    try:
        
        for weatherDataset in weather_datasets:
            """
            dataPath = filePath+'s2_index__L3='+s2_index__L3_list[0]+'/s2_index__L5='+s2_index__L5_list[0]+'/s2_index__L8='+s2_index__L8_list[0] 
            if os.path.exists(dataPath):
                dataDir = dataPath
            else:
                dataDir = []
                return []
            
            #get the parquet dataset from the directory
            weatherDataset = ds.dataset(dataDir,format="parquet", partitioning="hive")
            """
            
            dt = datetime.strptime(YYYY_str+'-'+MM_str+'-'+DD_str, '%Y-%m-%d')

            weather_df = weatherDataset.to_table(
                columns=['Year','Month','Day','ETo_AVG_IN'],
                filter=(
                    (ds.field('Year') >= int(YYYY_str)) & (ds.field('Year') <= end_YYYY_str) &
                    (ds.field('Month') >= int(MM_str)) & (ds.field('Month') <= end_MM_str) &
                    (ds.field('Day') >= int(DD_str)) & (ds.field('Day') <= end_DD_str)
                )
            ).to_pandas()      
            
            #get the DataFrame for the average of that day
            weather_df['YYYY']=weather_df['Year'].astype(str)
            weather_df['MM']=weather_df['Month'].astype(str).str.zfill(2)
            weather_df['DD']=weather_df['Day'].astype(str).str.zfill(2)

            weather_df['Date']=weather_df.YYYY+'-'+weather_df.MM+'-'+weather_df.DD
            weather_df.Date = pd.to_datetime(weather_df.Date)

            weather_df = weather_df.drop(columns=['YYYY','MM','DD'], axis=1)
            w_df = pd.DataFrame(weather_df.mean(axis=0)).T
            w_df['Date']=dt

            cols = ['Date','ETo_AVG_IN']
            w_df = w_df[cols]
            w_all = pd.concat([w_all, w_df], ignore_index=True)
    except Exception as e:
        print(e)
            
        
    if len(w_all) > 0:
        w_all['ETo_AVG_IN'] = pd.to_numeric(w_all['ETo_AVG_IN'], errors='coerce')
        w_all = w_all.dropna(subset=['ETo_AVG_IN'])  # Drop rows with NaN in ETo_AVG_IN

        if not w_all.empty:
            # Take the average of the columns
            w_ret = pd.DataFrame(w_all.groupby(["Date"])["ETo_AVG_IN"].mean())
            w_ret.reset_index(inplace=True)
        else:
            w_ret = pd.DataFrame()
    else:
        w_ret = pd.DataFrame()

    
    end_time = time.time()
    time_elapsed = (end_time - start_time)
    
    return  w_ret

def getEtoFromAUS(lat, lon, filePath, start_date, end_date):
    lats = [lat]
    lons = [lon]
    start_time = time.time()

    # Get the list of S2 indices and CIDs for the data point
    s2_index__L5_list, L5_cids = get_s2_cellids_and_token_list(5, lats, lons)

    tok = start_date.split('-')
    YYYY_str = tok[0]
    MM_str = tok[1]
    DD_str = tok[2]

    if end_date:
        # Parse end date and localize it to UTC if provided
        end_tok = end_date.split('-')
        end_YYYY_str = int(end_tok[0])
        end_MM_str = int(end_tok[1])
        end_DD_str = int(end_tok[2])
        end_local_dt = datetime(end_YYYY_str, end_MM_str, end_DD_str)
        utc_end_dt = pytz.utc.localize(end_local_dt)
        
    list_of_L5_paths = [filePath + 's2_token_L5=' + x for x in s2_index__L5_list if os.path.exists(filePath + 's2_token_L5=' + x)]

    weather_datasets = []
    for x in list_of_L5_paths:
        weather_datasets.append(ds.dataset(x, format="parquet", partitioning="hive"))

    # Get the data
    w_all = pd.DataFrame()
    try:
        for weatherDataset in weather_datasets:
            # List the columns available in the dataset to identify the correct field name
            # print(f"Columns in dataset: {weatherDataset.schema}")
            
            yrStr = YYYY_str
            moStr = MM_str
            dtStr = DD_str

            dt = datetime.strptime(yrStr + '-' + moStr + '-' + dtStr, '%Y-%m-%d')

            # Adjust the field name to match the actual column names in the dataset
            weather_df = weatherDataset.to_table(
                columns=['Year', 'Month', 'Day', 'ETo_AVG_IN'],  # Correct column name here
                filter=(
                    (ds.field('Year') >= int(YYYY_str)) & (ds.field('Year') <= end_YYYY_str) &
                    (ds.field('Month') >= int(MM_str)) & (ds.field('Month') <= end_MM_str) &
                    (ds.field('Day') >= int(DD_str)) & (ds.field('Day') <= end_DD_str)
                )
            ).to_pandas()

            #get the DataFrame for the average of that day
            weather_df['YYYY']=weather_df['Year'].astype(str)
            weather_df['MM']=weather_df['Month'].astype(str).str.zfill(2)
            weather_df['DD']=weather_df['Day'].astype(str).str.zfill(2)

            weather_df['Date']=weather_df.YYYY+'-'+weather_df.MM+'-'+weather_df.DD
            weather_df.Date = pd.to_datetime(weather_df.Date)

            weather_df = weather_df.drop(columns=['YYYY','MM','DD'], axis=1)
                    # weather_df = weather_df.drop(columns=['YYYY','MM','DD'], axis=1)
            numeric_cols = weather_df.select_dtypes(include=['number']).columns
            
            if not numeric_cols.empty:
                w_df = pd.DataFrame(weather_df[numeric_cols].mean(axis=0)).T
            else:
                print("No numeric columns found for mean calculation.")
                continue
            
            w_df['Date']=dt

            cols = ['Date','ETo_AVG_IN']
            w_df = w_df[cols]
            w_all = pd.concat([w_all, w_df], ignore_index=True)
    except Exception as e:
        w_ret=pd.DataFrame()    
    
    if len(w_all) > 0:
        w_all['ETo_AVG_IN'] = pd.to_numeric(w_all['ETo_AVG_IN'], errors='coerce')
        w_all = w_all.dropna(subset=['ETo_AVG_IN'])  # Drop rows with NaN in ETo_AVG_IN

        if not w_all.empty:
            # Take the average of the columns
            w_ret = pd.DataFrame(w_all.groupby(["Date"])["ETo_AVG_IN"].mean())
            w_ret.reset_index(inplace=True)
        else:
            w_ret = pd.DataFrame()
    else:
        w_ret = pd.DataFrame()

    
    end_time = time.time()
    time_elapsed = (end_time - start_time)
    
    return w_ret

def getEtoFromAUSForecast(lat, lon, filePath,start_date,end_date):
    lats=[lat]
    lons=[lon]
    start_time = time.time()
    #get the list of S2 indeces and CIDs for the data point
    s2_index__L5_list, L5_cids = get_s2_cellids_and_token_list(5, lats, lons)
    
    #print(s2_index__L3_list)
    # print(f'AUS Forecast--{s2_index__L5_list}')
    #print(s2_index__L8_list)
    
    tok = start_date.split('-')
    YYYY_str = tok[0]
    MM_str = tok[1]
    DD_str = tok[2]
    
    if end_date:
        # Parse end date and localize it to UTC if provided
        end_tok = end_date.split('-')
        end_YYYY_str = int(end_tok[0])
        end_MM_str = int(end_tok[1])
        end_DD_str = int(end_tok[2])
        end_local_dt = datetime(end_YYYY_str, end_MM_str, end_DD_str)
        utc_end_dt = pytz.utc.localize(end_local_dt)

    list_of_L5_paths = []
    
    # s2_index__L7_list = ['8099c']
    # s2_index__L9_list = ['8099bc']

    list_of_L5_paths = [filePath+'s2_token_L5='+x for x in s2_index__L5_list 
                        if os.path.exists(filePath+'s2_token_L5='+x)]
    
    
    weather_datasets = []
        
    for x in list_of_L5_paths:
        weather_datasets.append(ds.dataset(x, format="parquet", partitioning="hive"))
    #get the Data
    try:
        w_all = pd.DataFrame()
        for weatherDataset in weather_datasets:
            """
            dataPath = filePath+'s2_index__L3='+s2_index__L3_list[0]+'/s2_index__L5='+s2_index__L5_list[0]+'/s2_index__L8='+s2_index__L8_list[0] 
            if os.path.exists(dataPath):
                dataDir = dataPath
            else:
                dataDir = []
                return []
            
            #get the parquet dataset from the directory
            weatherDataset = ds.dataset(dataDir,format="parquet", partitioning="hive")
            """
            yrStr=  YYYY_str
            moStr = MM_str
            dtStr = DD_str
            
            dt = datetime.strptime(yrStr+'-'+moStr+'-'+dtStr, '%Y-%m-%d')
                
            weather_df = weatherDataset.to_table(
                columns=['Year','Month','Day','ETo_average_inches'],
                filter=(
                    (ds.field('Year') >= int(YYYY_str)) & (ds.field('Year') <= end_YYYY_str) &
                    (ds.field('Month') >= int(MM_str)) & (ds.field('Month') <= end_MM_str) &
                    (ds.field('Day') >= int(DD_str)) & (ds.field('Day') <= end_DD_str)
                )
            ).to_pandas()

            #get the DataFrame for the average of that day
            weather_df['YYYY']=weather_df['Year'].astype(str)
            weather_df['MM']=weather_df['Month'].astype(str).str.zfill(2)
            weather_df['DD']=weather_df['Day'].astype(str).str.zfill(2)

            weather_df['Date']=weather_df.YYYY+'-'+weather_df.MM+'-'+weather_df.DD
            weather_df.Date = pd.to_datetime(weather_df.Date)

            weather_df = weather_df.drop(columns=['YYYY','MM','DD'], axis=1)
            w_df = pd.DataFrame(weather_df.mean(axis=0)).T
            w_df['Date']=dt

            cols = ['Date','ETo_average_inches']
            w_df = w_df[cols]
            w_all = pd.concat([w_all, w_df], ignore_index=True)
    except Exception as e:
        w_ret=pd.DataFrame()

        
    if len(w_all) > 0:
        w_all['ETo_AVG_IN'] = pd.to_numeric(w_all['ETo_average_inches'], errors='coerce')
        w_all = w_all.dropna(subset=['ETo_AVG_IN'])  # Drop rows with NaN in ETo_AVG_IN

        if not w_all.empty:
            # Take the average of the columns
            w_ret = pd.DataFrame(w_all.groupby(["Date"])["ETo_AVG_IN"].mean())
            w_ret.reset_index(inplace=True)
        else:
            w_ret = pd.DataFrame()
    else:
        w_ret = pd.DataFrame()

    
    end_time = time.time()
    time_elapsed = (end_time - start_time)
    
    return w_ret

def getEtoFromNOAA(lat, lon, filePath,start_date,end_date):
    lats=[lat]
    lons=[lon]
    start_time = time.time()
    #get the list of S2 indeces and CIDs for the data point
    s2_index__L3_list, L3_cids = get_s2_cellids_and_token_list(3, lats, lons)
    
    #print(s2_index__L3_list)
    # print(f'NOAA DAily-{s2_index__L3_list}')
    #print(s2_index__L8_list)
    
    tok = start_date.split('-')
    YYYY_str = tok[0]
    MM_str = tok[1]
    DD_str = tok[2]
    
    if end_date:
        # Parse end date and localize it to UTC if provided
        end_tok = end_date.split('-')
        end_YYYY_str = int(end_tok[0])
        end_MM_str = int(end_tok[1])
        end_DD_str = int(end_tok[2])
        end_local_dt = datetime(end_YYYY_str, end_MM_str, end_DD_str)
        utc_end_dt = pytz.utc.localize(end_local_dt)

    list_of_L3_paths = []
    
    # s2_index__L7_list = ['8099c']
    # s2_index__L9_list = ['8099bc']

    list_of_L3_paths = [filePath+'s2_index__L3='+x for x in s2_index__L3_list 
                        if os.path.exists(filePath+'s2_index__L3='+x)]
    
    
    weather_datasets = []
        
    for x in list_of_L3_paths:
        weather_datasets.append(ds.dataset(x, format="parquet", partitioning="hive"))
    #get the Data
    w_all = pd.DataFrame()
    try:
        
        for weatherDataset in weather_datasets:
            """
            dataPath = filePath+'s2_index__L3='+s2_index__L3_list[0]+'/s2_index__L5='+s2_index__L5_list[0]+'/s2_index__L8='+s2_index__L8_list[0] 
            if os.path.exists(dataPath):
                dataDir = dataPath
            else:
                dataDir = []
                return []
            
            #get the parquet dataset from the directory
            weatherDataset = ds.dataset(dataDir,format="parquet", partitioning="hive")
            """
            yrStr=  YYYY_str
            moStr = MM_str
            dtStr = DD_str
            
            dt = datetime.strptime(yrStr+'-'+moStr+'-'+dtStr, '%Y-%m-%d')

                
            weather_df = weatherDataset.to_table(
                columns=['YYYY', 'MM', 'DD', 'ETo_AVG_IN'],  # Use the correct column name
                filter=(
                    (ds.field('YYYY') >= int(YYYY_str)) & (ds.field('YYYY') <= int(end_YYYY_str)) &
                    (ds.field('MM') >= int(MM_str)) & (ds.field('MM') <= int(end_MM_str)) &
                    (ds.field('DD') >= int(DD_str)) & (ds.field('DD') <= int(end_DD_str))
                )
            ).to_pandas()


            # #get the DataFrame for the average of that day
            weather_df['YYYY']=weather_df['YYYY'].astype(str)
            weather_df['MM']=weather_df['MM'].astype(str).str.zfill(2)
            weather_df['DD']=weather_df['DD'].astype(str).str.zfill(2)

            weather_df['Date']=weather_df.YYYY+'-'+weather_df.MM+'-'+weather_df.DD
            weather_df.Date = pd.to_datetime(weather_df.Date)

            # weather_df = weather_df.drop(columns=['YYYY','MM','DD'], axis=1)
            numeric_cols = weather_df.select_dtypes(include=['number']).columns
            
            if not numeric_cols.empty:
                w_df = pd.DataFrame(weather_df[numeric_cols].mean(axis=0)).T
            else:
                print("No numeric columns found for mean calculation.")
                continue
            w_df['Date']=dt

            cols = ['Date','ETo_AVG_IN']
            w_df = w_df[cols]
            w_all = pd.concat([w_all, w_df], ignore_index=True)
    except Exception as e:
        w_ret = pd.DataFrame()
    
    if len(w_all) > 0:
        w_all['ETo_AVG_IN'] = pd.to_numeric(w_all['ETo_AVG_IN'], errors='coerce')
        w_all = w_all.dropna(subset=['ETo_AVG_IN'])  # Drop rows with NaN in ETo_AVG_IN

        if not w_all.empty:
            # Take the average of the columns
            w_ret = pd.DataFrame(w_all.groupby(["Date"])["ETo_AVG_IN"].mean())
            w_ret.reset_index(inplace=True)
        else:
            w_ret = pd.DataFrame()
    else:
        w_ret = pd.DataFrame()

    
    end_time = time.time()
    time_elapsed = (end_time - start_time)
    
    return w_ret
    

def getEtoFromNOAAForecast(lat, lon, filePath,start_date,end_date):
    lats=[lat]
    lons=[lon]
    start_time = time.time()
    #get the list of S2 indeces and CIDs for the data point
    s2_index__L5_list, L5_cids = get_s2_cellids_and_token_list(5, lats, lons)

    #print(s2_index__L3_list)
    # print(f'NOAA Forecast-{s2_index__L5_list}')
    #print(s2_index__L8_list)
    
    tok = start_date.split('-')
    YYYY_str = tok[0]
    MM_str = tok[1]
    DD_str = tok[2]
    
    if end_date:
        # Parse end date and localize it to UTC if provided
        end_tok = end_date.split('-')
        end_YYYY_str = int(end_tok[0])
        end_MM_str = int(end_tok[1])
        end_DD_str = int(end_tok[2])
        end_local_dt = datetime(end_YYYY_str, end_MM_str, end_DD_str)
        utc_end_dt = pytz.utc.localize(end_local_dt)


    list_of_L5_paths = []    
    # s2_index__L7_list = ['8099c']
    # s2_index__L9_list = ['8099bc']

    list_of_L5_paths = [filePath+'s2_tokens_l5='+x for x in s2_index__L5_list 
                        if os.path.exists(filePath+'s2_tokens_l5='+x)]
    
    weather_datasets = []
    for x in list_of_L5_paths:
        weather_datasets.append(ds.dataset(x, format="parquet", partitioning="hive"))
    #get the Data
    w_all = pd.DataFrame()
    try:
        for weatherDataset in weather_datasets:
            """
            dataPath = filePath+'s2_index__L3='+s2_index__L3_list[0]+'/s2_index__L5='+s2_index__L5_list[0]+'/s2_index__L8='+s2_index__L8_list[0] 
            if os.path.exists(dataPath):
                dataDir = dataPath
            else:
                dataDir = []
                return []
            
            #get the parquet dataset from the directory
            weatherDataset = ds.dataset(dataDir,format="parquet", partitioning="hive")
            """
            yrStr=  YYYY_str
            moStr = MM_str
            dtStr = DD_str
            
            dt = datetime.strptime(yrStr+'-'+moStr+'-'+dtStr, '%Y-%m-%d')

                
            weather_df = weatherDataset.to_table(
                columns=['Year', 'Month', 'Day', 'ETo_average_inches'],  # Use the correct column name
                filter=(
                    (ds.field('Year') >= int(YYYY_str)) & (ds.field('Year') <= int(end_YYYY_str)) &
                    (ds.field('Month') >= int(MM_str)) & (ds.field('Month') <= int(end_MM_str)) &
                    (ds.field('Day') >= int(DD_str)) & (ds.field('Day') <= int(end_DD_str))
                )
            ).to_pandas()



            # #get the DataFrame for the average of that day
            weather_df['YYYY']=weather_df['Year'].astype(str)
            weather_df['MM']=weather_df['Month'].astype(str).str.zfill(2)
            weather_df['DD']=weather_df['Day'].astype(str).str.zfill(2)

            weather_df['Date']=weather_df.YYYY+'-'+weather_df.MM+'-'+weather_df.DD
            weather_df.Date = pd.to_datetime(weather_df.Date)

            weather_df = weather_df.drop(columns=['YYYY','MM','DD'], axis=1)
            w_df = pd.DataFrame(weather_df.mean(axis=0)).T
            w_df['Date']=dt
            # print(w_df.columns)
            cols = ['Date','ETo_average_inches']
            w_df = w_df[cols]
            w_all = pd.concat([w_all, w_df], ignore_index=True)

    except Exception as e:
        w_ret = pd.DataFrame()     
        
        
    if len(w_all) > 0:
        w_all['ETo_AVG_IN'] = pd.to_numeric(w_all['ETo_average_inches'], errors='coerce')
        w_all = w_all.dropna(subset=['ETo_AVG_IN'])  # Drop rows with NaN in ETo_AVG_IN

        if not w_all.empty:
            # Take the average of the columns
            w_ret = pd.DataFrame(w_all.groupby(["Date"])["ETo_AVG_IN"].mean())
            w_ret.reset_index(inplace=True)
        else:
            w_ret = pd.DataFrame()
    else:
        w_ret = pd.DataFrame()

    
    end_time = time.time()
    time_elapsed = (end_time - start_time)

    return w_ret


def getEtoFromGHCND(lat, lon, filePath, start_date,end_date):
    lats=[lat]
    lons=[lon]
    
    #Parse the date
    
    tok = start_date.split('-')
    YYYY_str = tok[0]
    MM_str = tok[1]
    DD_str = tok[2]
    
    if end_date:
        # Parse end date and localize it to UTC if provided
        end_tok = end_date.split('-')
        end_YYYY_str = int(end_tok[0])
        end_MM_str = int(end_tok[1])
        end_DD_str = int(end_tok[2])
        end_local_dt = datetime(end_YYYY_str, end_MM_str, end_DD_str)
        utc_end_dt = pytz.utc.localize(end_local_dt)
    start_time = time.time()
    #get the list of S2 indeces and CIDs for the data point
    s2_index__L5_list, L5_cids = get_s2_cellids_and_token_list(5, lats, lons)
    

    # print(s2_index__L5_list)
    # print(s2_index__L7_list)
    # print(s2_index__L9_list)
    # s2_index__L7_list = ['8099c']
    # s2_index__L9_list = ['8099bc']


    list_of_L5_paths = []
    list_of_L5_paths = [filePath+'s2_token_L5='+x for x in s2_index__L5_list 
                        if os.path.exists(filePath+'s2_token_L5='+x)]

    
    weather_datasets = []
    for x in list_of_L5_paths:
        weather_datasets.append(ds.dataset(x, format="parquet", partitioning="hive"))
    #get the Data
    w_all = pd.DataFrame()
    for weatherDataset in weather_datasets:
        """
        dataPath = filePath+'s2_index__L3='+s2_index__L3_list[0]+'/s2_index__L5='+s2_index__L5_list[0]+'/s2_index__L8='+s2_index__L8_list[0] 
        if os.path.exists(dataPath):
            dataDir = dataPath
        else:
            dataDir = []
            return []
        
        #get the parquet dataset from the directory
        weatherDataset = ds.dataset(dataDir,format="parquet", partitioning="hive")
        """
        yrStr=  YYYY_str
        moStr = MM_str
        dtStr = DD_str
        
        dt = datetime.strptime(yrStr+'-'+moStr+'-'+dtStr, '%Y-%m-%d')
        try:
            
            weather_df = weatherDataset.to_table(
                columns=['Year','Month','Day','ETo_AVG_IN'],
                filter=(
                    (ds.field('Year') >= int(YYYY_str)) & (ds.field('Year') <= end_YYYY_str) &
                    (ds.field('Month') >= int(MM_str)) & (ds.field('Month') <= end_MM_str) &
                    (ds.field('Day') >= int(DD_str)) & (ds.field('Day') <= end_DD_str)
                )
            ).to_pandas()

            #get the DataFrame for the average of that day
            weather_df['YYYY']=weather_df['Year'].astype(str)
            weather_df['MM']=weather_df['Month'].astype(str).str.zfill(2)
            weather_df['DD']=weather_df['Day'].astype(str).str.zfill(2)

            weather_df['Date']=weather_df.YYYY+'-'+weather_df.MM+'-'+weather_df.DD
            weather_df.Date = pd.to_datetime(weather_df.Date)

            weather_df = weather_df.drop(columns=['YYYY','MM','DD'], axis=1)
            w_df = pd.DataFrame(weather_df.mean(axis=0)).T
            w_df['Date']=dt

            cols = ['Date','ETo_AVG_IN']
            w_df = w_df[cols]
            w_all = pd.concat([w_all, w_df], ignore_index=True)
        except Exception as e:
            print(e)
    
    if len(w_all) > 0:
        w_all['ETo_AVG_IN'] = pd.to_numeric(w_all['ETo_AVG_IN'], errors='coerce')
        w_all = w_all.dropna(subset=['ETo_AVG_IN'])  # Drop rows with NaN in ETo_AVG_IN

        if not w_all.empty:
            # Take the average of the columns
            w_ret = pd.DataFrame(w_all.groupby(["Date"])["ETo_AVG_IN"].mean())
            w_ret.reset_index(inplace=True)
        else:
            w_ret = pd.DataFrame()
    else:
        w_ret = pd.DataFrame()

    
    end_time = time.time()
    time_elapsed = (end_time - start_time)
    
    return w_ret

def getEToFromCIMIS(lat, lon, filePath, start_date,end_date):
          
    #Parse the date
    
    tok = start_date.split('-')
    YYYY_str = tok[0]
    MM_str = tok[1]
    DD_str = tok[2]
    
    if end_date:
        # Parse end date and localize it to UTC if provided
        end_tok = end_date.split('-')
        end_YYYY_str = int(end_tok[0])
        end_MM_str = int(end_tok[1])
        end_DD_str = int(end_tok[2])
        end_local_dt = datetime(end_YYYY_str, end_MM_str, end_DD_str)
        utc_end_dt = pytz.utc.localize(end_local_dt)
    
    
    ## Get S2 indices and CIDs
    s2_index__L5_list, L5_cids = get_s2_cellids_and_token_list(5, [lat], [lon])
    
    # Generate paths for available datasets
    list_of_L5_paths = [
        os.path.join(filePath, f"s2_token_L5={x}")
        for x in s2_index__L5_list
        if os.path.exists(os.path.join(filePath, f"s2_token_L5={x}"))
    ]
    
    # Load weather datasets
    weather_datasets = [
        ds.dataset(path, format="parquet", partitioning="hive")
        for path in list_of_L5_paths
    ]
    
    
    w_all = pd.DataFrame()
    w_ret=pd.DataFrame()
    for weatherDataset in weather_datasets:
        
        """
            dataPath = filePath+'s2_index__L3='+s2_index__L3_list[0]+'/s2_index__L5='+s2_index__L5_list[0]+'/s2_index__L8='+s2_index__L8_list[0] 
            if os.path.exists(dataPath):
                dataDir = dataPath
            else:
                dataDir = []
                return []
            
            #get the parquet dataset from the directory
            weatherDataset = ds.dataset(dataDir,format="parquet", partitioning="hive")
        """
        yrStr=  YYYY_str
        moStr = MM_str
        dtStr = DD_str
        
        dt = datetime.strptime(yrStr+'-'+moStr+'-'+dtStr, '%Y-%m-%d')
        
        weather_df = weatherDataset.to_table(
            columns=["Year", "Month", "Day", "HlyEto"],
            filter=(
                (ds.field('Year') >= int(YYYY_str)) & (ds.field('Year') <= end_YYYY_str) &
                (ds.field('Month') >= int(MM_str)) & (ds.field('Month') <= end_MM_str) &
                (ds.field('Day') >= int(DD_str)) & (ds.field('Day') <= end_DD_str)
            )
        ).to_pandas()
        
                #get the DataFrame for the average of that day
        weather_df['YYYY']=weather_df['Year'].astype(str)
        weather_df['MM']=weather_df['Month'].astype(str).str.zfill(2)
        weather_df['DD']=weather_df['Day'].astype(str).str.zfill(2)

        weather_df['Date']=weather_df.YYYY+'-'+weather_df.MM+'-'+weather_df.DD
        weather_df.Date = pd.to_datetime(weather_df.Date)

        weather_df = weather_df.drop(columns=['YYYY','MM','DD'], axis=1)
        # Convert HlyEto to numeric (handle coercion)
        weather_df['HlyEto'] = pd.to_numeric(weather_df['HlyEto'], errors='coerce')

        # Calculate mean for numeric columns
        numeric_cols = weather_df.select_dtypes(include=['number']).columns
        if not numeric_cols.empty:
            w_df = pd.DataFrame(weather_df[numeric_cols].mean(axis=0)).T
        else:
            print("No numeric columns found for mean calculation.")
            continue

        w_df['Date']=dt
        
        cols = ['Date','HlyEto']
        w_df = w_df[cols]
        w_all = pd.concat([w_all, w_df], ignore_index=True)         
            
    if len(w_all) > 0:
        w_all['ETo_AVG_IN'] = pd.to_numeric(w_all['HlyEto'], errors='coerce')
        w_all = w_all.dropna(subset=['ETo_AVG_IN'])  # Drop rows with NaN in ETo_AVG_IN

        if not w_all.empty:
            # Take the average of the columns
            w_ret = pd.DataFrame(w_all.groupby(["Date"])["ETo_AVG_IN"].mean())
            w_ret.reset_index(inplace=True)
        else:
            w_ret = pd.DataFrame()
    else:
        w_ret = pd.DataFrame()


    return w_ret


# from bisect import bisect_left , bisect_right

# def filter_by_date_range(df,start_date,end_date):
    
#     # ensure date column is sorted
    
#     df = df.sort_values(by='Date')
#     date_list = df['Date'].tolist()
    
#     # Use binary search for fast filtering
#     start_idx = bisect_left(date_list , start_date)
#     end_idx = bisect_right(date_list , end_date)
    
#     return df.iloc[start_idx:end_idx]



# def getEToFromCIMIS(lat, lon, filePath, start_date,end_date):
    
#     # parse the start date and end date
#     start_dt = datetime.strptime(start_date,'%Y-%m-%d')
#     utc_end_dt = None
#     if end_date:
#         end_dt = datetime.strptime(end_date,'%Y-%m-%d')
#         utc_end_dt = pytz.utc.localize(end_dt)
    
#     # Get S2 indices and CIDs
#     s2_index__L5_list, L5_cids = get_s2_cellids_and_token_list(5, [lat], [lon])
    
#     # Generate paths for available datasets
#     list_of_L5_paths = [
#         os.path.join(filePath, f"s2_token_L5={x}")
#         for x in s2_index__L5_list
#         if os.path.exists(os.path.join(filePath, f"s2_token_L5={x}"))
#     ]
    
#     # Load weather datasets
#     weather_datasets = [
#         ds.dataset(path, format="parquet", partitioning="hive")
#         for path in list_of_L5_paths
#     ]
    
#     # Prepare filters for dataset
#     filters = (
#         (ds.field("Year") >= start_dt.year) & (ds.field("Year") <= end_dt.year) &
#         (ds.field("Month") >= start_dt.month) & (ds.field("Month") <= end_dt.month) &
#         (ds.field("Day") >= start_dt.day) & (ds.field("Day") <= end_dt.day)
#     )
    
#     weather_frames = []
#     for weather_dataset in weather_datasets:
#         weather_df = weather_dataset.to_table(
#             columns=["Year", "Month", "Day", "HlyEto"], filter=filters
#         ).to_pandas()

#         if not weather_df.empty:
#             weather_df["Date"] = pd.to_datetime(
#                 weather_df[["Year", "Month", "Day"]]
#                     .astype(str)
#                     .agg("-".join, axis=1)
#             )
#             weather_df = weather_df.drop(columns=["Year", "Month", "Day"])
#             weather_df = filter_by_date_range(weather_df, start_dt, end_dt)  # Binary search filtering
#             print(weather_df.columns)
#             weather_frames.append(weather_df)

#      # Concatenate and aggregate data
#     if weather_frames:
#         w_all = pd.concat(weather_frames, ignore_index=True)
#         print(w_all["HlyEto"])
#         w_all["ETo_AVG_IN"] = w_all["HlyEto"]
#         w_ret = w_all.groupby("Date", as_index=False)["ETo_AVG_IN"].sum()
#     else:
#         w_ret = pd.DataFrame(columns=["Date", "ETo_AVG_IN"])
    
#     return w_ret


################################################################
################################################################	
#********** DATA
#get the USA Mainland Polygon
# print('\nBuilding Boundary Polygons...')
# print('\n\tBuilding Contiguous USA Mainland Boundary Polygon...')
# uspolys = pd.DataFrame(columns=['Area','geometry'])
# us_boundary = data_dir + 'cb_2018_us_nation_5m.shp'
# us_gdf = gpd.read_file(us_boundary, driver='SHAPEFILE')
# us_polygon_list = list(us_gdf.geometry[0])
# i=1
# for geom in us_polygon_list:
#     a = getArea(geom)
#     uspolys.loc[i,'Area']=a
#     uspolys.loc[i,'geometry']=geom
#     i=i+1
# uspolys.sort_values(by='Area',inplace=True)
# uspolys.reset_index(drop=True,inplace=True)
# usap = uspolys.loc[len(uspolys)-1,'geometry']    

# print('\n\tBuilding USA States Boundary Polygons...')
# usstates_boundary = data_dir + 'tl_2019_us_state.shp'
# usstates_gdf = gpd.read_file(usstates_boundary, driver='SHAPEFILE')
# usstates_gdf = usstates_gdf.to_crs(4326)

# print('\n\tBuilding CA Boundary Polygon...')
# ca_gdf = usstates_gdf[usstates_gdf.NAME=='California']
# ca_gdf.reset_index(drop=True, inplace=True)
# cap = ca_gdf.geometry[0]

# print('\n\tBuilding AZ Boundary Polygon...')
# az_gdf = usstates_gdf[usstates_gdf.NAME=='Arizona']
# az_gdf.reset_index(drop=True, inplace=True)
# azp = az_gdf.geometry[0]

# print('\n\tBuilding ID Boundary Polygon...')
# id_gdf = usstates_gdf[usstates_gdf.NAME=='Idaho']
# id_gdf.reset_index(drop=True, inplace=True)
# idp = id_gdf.geometry[0]

# print('\n\tBuilding IA Boundary Polygon...')
# ia_gdf = usstates_gdf[usstates_gdf.NAME=='Iowa']
# ia_gdf.reset_index(drop=True, inplace=True)
# iap = ia_gdf.geometry[0]

# print('\n\tBuilding WA Boundary Polygon...')
# wa_gdf = usstates_gdf[usstates_gdf.NAME=='Washington']
# wa_gdf.reset_index(drop=True, inplace=True)
# wap = wa_gdf.geometry[0]

# print('\n\tBuilding TX Boundary Polygon...')
# tx_gdf = usstates_gdf[usstates_gdf.NAME=='Texas']
# tx_gdf.reset_index(drop=True, inplace=True)
# txp = tx_gdf.geometry[0]

# print('\n\tBuilding FL Boundary Polygon...')
# fl_gdf = usstates_gdf[usstates_gdf.NAME=='Florida']
# fl_gdf.reset_index(drop=True, inplace=True)
# flp = fl_gdf.geometry[0]

# print('\n\tBuilding OR Boundary Polygon...')
# or_gdf = usstates_gdf[usstates_gdf.NAME=='Oregon']
# or_gdf.reset_index(drop=True, inplace=True)
# orp = or_gdf.geometry[0]


# print('\n\tBuilding OK Boundary Polygon...')
# ok_gdf = usstates_gdf[usstates_gdf.NAME=='Oklahoma']
# ok_gdf.reset_index(drop=True, inplace=True)
# okp = ok_gdf.geometry[0]


# print('\n\tBuilding MT Boundary Polygon...')
# mt_gdf = usstates_gdf[usstates_gdf.NAME=='Montana']
# mt_gdf.reset_index(drop=True, inplace=True)
# mtp = mt_gdf.geometry[0]


# print('\n\tBuilding MO Boundary Polygon...')
# mo_gdf = usstates_gdf[usstates_gdf.NAME=='Missouri']
# mo_gdf.reset_index(drop=True, inplace=True)
# mop = mo_gdf.geometry[0]


# print('\n\tBuilding DE Boundary Polygon...')
# de_gdf = usstates_gdf[usstates_gdf.NAME=='Delaware']
# de_gdf.reset_index(drop=True, inplace=True)
# dep = de_gdf.geometry[0]


# print('\n\tBuilding AR Boundary Polygon...')
# ar_gdf = usstates_gdf[usstates_gdf.NAME=='Arkansas']
# ar_gdf.reset_index(drop=True, inplace=True)
# arp = ar_gdf.geometry[0]


# print('\nBuilding Global Boundary Polygons...')
# wrs_gdf = gpd.read_file(worldShpFile)
# wrs_gdf = wrs_gdf.to_crs(4326)


# print('\n\tBuilding ISR Boundary Polygon...')
# isr_gdf = wrs_gdf[wrs_gdf.CNTRY_NAME=='Israel']
# isr_gdf.reset_index(drop=True, inplace=True)
# isrp = isr_gdf.geometry[0]

# print('\n\tBuilding IND Boundary Polygon...')
# ind_gdf = wrs_gdf[wrs_gdf.CNTRY_NAME=='India']
# ind_gdf.reset_index(drop=True, inplace=True)
# ind_gdf = ind_gdf.to_crs(4326)
# ind_poly_list = [row.geometry for i,row in ind_gdf.iterrows()]
# indp_list = []
# for p in ind_poly_list:
#     if type(p)==shapely.geometry.multipolygon.MultiPolygon:
#         indp_list=indp_list+list(p) #Break up the polygons itn
#     else:
#         indp_list.append(p)
# print('Done!')
# #indp = ind_gdf.geometry[0]

# print('\n\tBuilding MEX Boundary Polygon...')
# mex_gdf = wrs_gdf[wrs_gdf.CNTRY_NAME=='Mexico']
# mex_gdf.reset_index(drop=True, inplace=True)
# mex_gdf = mex_gdf.to_crs(4326)
# mex_poly_list = [row.geometry for i,row in mex_gdf.iterrows()]
# mexp_list = []
# for p in mex_poly_list:
#     if type(p)==shapely.geometry.multipolygon.MultiPolygon:
#         mexp_list=mexp_list+list(p) #Break up the polygons itn
#     else:
#         mexp_list.append(p)
# print('Done!')

# print('\n\tBuilding EGY Boundary Polygon...')
# egy_gdf = wrs_gdf[wrs_gdf.CNTRY_NAME=='Egypt']
# egy_gdf.reset_index(drop=True, inplace=True)
# egy_gdf = egy_gdf.to_crs(4326)
# egy_poly_list = [row.geometry for i,row in egy_gdf.iterrows()]
# egyp_list = []
# for p in egy_poly_list:
#     if type(p)==shapely.geometry.multipolygon.MultiPolygon:
#         egyp_list=egyp_list+list(p) #Break up the polygons itn
#     else:
#         egyp_list.append(p)
# print('Done!')

# print('\n\tBuilding CAN Boundary Polygon...')
# can_gdf = wrs_gdf[wrs_gdf.CNTRY_NAME=='Canada']
# can_gdf.reset_index(drop=True, inplace=True)
# can_gdf = can_gdf.to_crs(4326)
# can_poly_list = [row.geometry for i,row in can_gdf.iterrows()]
# canp_list = []
# for p in can_poly_list:
#     if type(p)==shapely.geometry.multipolygon.MultiPolygon:
#         canp_list=canp_list+list(p) #Break up the polygons itn
#     else:
#         canp_list.append(p)
# print('Done!')

# print('\n\tBuilding PAK Boundary Polygon...')
# pak_gdf = wrs_gdf[wrs_gdf.CNTRY_NAME=='Pakistan']
# pak_gdf.reset_index(drop=True, inplace=True)
# pak_gdf = pak_gdf.to_crs(4326)
# pak_poly_list = [row.geometry for i,row in pak_gdf.iterrows()]
# pakp_list = []
# for p in pak_poly_list:
#     if type(p)==shapely.geometry.multipolygon.MultiPolygon:
#         pakp_list=pakp_list+list(p) #Break up the polygons itn
#     else:
#         pakp_list.append(p)
# print('Done!')


# print('\n\tBuilding UK Boundary Polygon...')
# uk_gdf = wrs_gdf[wrs_gdf.CNTRY_NAME=='United Kingdom']
# uk_gdf.reset_index(drop=True, inplace=True)
# ukp = uk_gdf.geometry[0]


# print('\n\tBuilding AUS Boundary Polygons...')
# """
# list_of_ausstate_shapefiles = [aus_data_dir+x for x in os.listdir(aus_data_dir) if x.endswith('.shp')]
# aus_gdf = gpd.GeoDataFrame()
# for stateShpFile in list_of_ausstate_shapefiles:
#     ausstates_gdf = gpd.read_file(stateShpFile, driver='SHAPEFILE')
#     aus_gdf = aus_gdf.append(ausstates_gdf, ignore_index=True)
# aus_gdf = aus_gdf.to_crs(4326)
# ausregions_gdf = aus_gdf.dissolve(by='STATE_PID')
# ausregions_gdf.reset_index(drop=True, inplace=True)
# ausregions_gdf = ausregions_gdf.to_crs(4326)
# aus_poly_list = [row.geometry for i,row in ausregions_gdf.iterrows()]


# aus_shapefile = aus_data_dir + 'IDM00013.shp'
# aus_gdf = gpd.read_file(aus_shapefile, driver='SHAPEFILE')
# aus_gdf = aus_gdf.to_crs(4326)
# """
# aus_gdf = wrs_gdf[wrs_gdf.CNTRY_NAME=='Australia']
# aus_gdf.reset_index(drop=True, inplace=True)
# aus_gdf = aus_gdf.to_crs(4326)
# aus_poly_list = [row.geometry for i,row in aus_gdf.iterrows()]


# ausp_list = []
# for p in aus_poly_list:
#     if type(p)==shapely.geometry.multipolygon.MultiPolygon:
#         ausp_list=ausp_list+list(p) #Break up the polygons itn
#     else:
#         ausp_list.append(p)
# print('Done!')


# print('\nGetting All Satellite Scenes...')
# gdf_S = getAllSentinelScenePolygons()
# gdf_L = getAllLandsatScenePolygons()
# S_gdf = gdf_S.copy() #for same call, different name
# print('Done!')


os.environ['ASSET_REGISTRY_BASE_URL'] = 'https://api-ar.agstack.org/'
from shapely import wkt

def fetchWKT(geoid):
    """
    Fetch the Well-Known Text (WKT) representation of a polygon for a given geoid.
    
    Args:
    geoid (str): The geoid for which to fetch the WKT polygon.
    
    Returns:
    str: The WKT polygon as a string.
    """
    # Construct the full URL for the API request
    url = os.getenv('ASSET_REGISTRY_BASE_URL') + "/fetch-field-wkt/" + str(geoid)
    # Make the GET request to the API
    res = requests.get(url)
    # Extract the WKT from the JSON response
    field_data = res.json()['WKT']
    # print(field_data)
    # Access the WKT string
    fieldWKT = field_data
    return fieldWKT

def extractLatLonFromWKT(wkt_polygon):
    """
    Extract the centroid latitude and longitude from the WKT polygon string using Shapely.
    
    Args:
    wkt_polygon (str): The WKT polygon as a string.
    
    Returns:
    tuple: A tuple containing the centroid's latitude and longitude.
    """
    try:
        # Parse the WKT polygon string into a Shapely geometry object
        fieldPoly = wkt.loads(wkt_polygon)
        # print(fieldPoly)
        # Calculate the centroid of the polygon
        c = fieldPoly.centroid
        # Extract the latitude and longitude from the centroid
        lat = c.y
        lon = c.x
        print(lat, lon)
        return lat, lon
    except Exception as e:
        print(f"Error processing WKT with Shapely: {e}")
        return None, None
    
    
session = requests.session()
session.headers = headers = {
    'Accept': 'application/json',
    'Content-Type': 'application/json'

}

asset_registry_base = "https://api-ar.agstack.org"

def polygon_to_s2_tokens(geoid, level):
    import s2sphere
    """
    Convert a geoid into a Polygon WKT string to a GeoDataFrame containing all S2 tokens 
    that cover the polygon, along with the geometry of each token as a Polygon.

    Parameters:
        polygon_wkt (str): The WKT string representing the polygon.
        level (int): The S2 cell level. Default is 10.

    Returns:
        geopandas.GeoDataFrame: A GeoDataFrame with two columns: 's2_token' and 'geometry'.
    """

    if not level:
        level=10
    #print('level='+str(level))

        # Define the request body (adjust based on API requirements)
    req_body = {
        "level": level  # Example key-value pair; modify as per the API documentation
    }

    res = session.get(asset_registry_base +"/fetch-field-wkt/"+geoid, json=req_body)
    polygon_wkt = res.json().get('WKT')
       
    # Convert WKT to Shapely Polygon
    polygon = wkt.loads(polygon_wkt)

    
    # Extract polygon vertices as LatLng points
    vertices = [
        s2sphere.LatLng.from_degrees(lat, lon)
        for lon, lat in polygon.exterior.coords
    ]
    #print('vertices='+str(vertices))
    
    # Use S2RegionCoverer to find covering cells
    region_coverer = s2sphere.RegionCoverer()
    region_coverer.min_level = level
    region_coverer.max_level = level

    #print(vertices[0], vertices[2])
    # Cover the bounding box of the polygon
    rect = s2sphere.LatLngRect.from_point_pair(vertices[0], vertices[2])
    cell_ids = region_coverer.get_covering(rect)
    #print(cell_ids)
    
    # Collect S2 tokens and geometries
    rows = []
    for cell_id in cell_ids:
        #print('level='+str(cell_id.level()))
        
        cell = s2sphere.Cell(cell_id)
        cell_vertices = [
            s2sphere.LatLng.from_point(cell.get_vertex(i))
            for i in range(4)
        ]
        coords = [
            (v.lng().degrees, v.lat().degrees)
            for v in cell_vertices
        ]
        coords.append(coords[0])  # Close the polygon
        
        # Create a Shapely Polygon
        cell_polygon = wkt.loads(
            f"POLYGON(({','.join([f'{lon} {lat}' for lon, lat in coords])}))"
        )
        
        # Append the S2 token and geometry
        try:
            rows.append({
                's2_index__L8': cell_id.parent(8).to_token(),
                's2_index__L10': cell_id.parent(10).to_token(),  # This can fail if level is too low
                's2_index__L'+str(level): cell_id.to_token(),
                'geometry': cell_polygon
            })
        except AssertionError:
            # Handle the error gracefully, maybe skip this row or log a warning
            pass


    # Create a GeoDataFrame for the Polygon
    field_gdf = gpd.GeoDataFrame({'geometry': [polygon]})
    
    # Add additional columns
    field_gdf['geoid'] = [geoid]
    
    # Set a CRS (Coordinate Reference System) if needed
    field_gdf.set_crs(epsg=4326, inplace=True)  # WGS84
    
    # Create and return both the GeoDataFrames
    gdf = gpd.GeoDataFrame(rows, columns=['s2_index__L8', 's2_index__L10', 's2_index__L'+str(level), 'geometry'])

    #convert to the same CRS 4326
    field_gdf = field_gdf.set_crs(epsg=4326)
    gdf = gdf.set_crs(epsg=4326)
    
    return gdf, field_gdf


def polygon_to_s2_tokens_tmf(geoid, level):
    import s2sphere

    if not level:
        level=8
    #print('level='+str(level))

        # Define the request body (adjust based on API requirements)
    req_body = {
        "level": level  # Example key-value pair; modify as per the API documentation
    }

    res = session.get(asset_registry_base +"/fetch-field-wkt/"+geoid, json=req_body)
    polygon_wkt = res.json().get('WKT')
       
    # Convert WKT to Shapely Polygon
    polygon = wkt.loads(polygon_wkt)

    
    # Extract polygon vertices as LatLng points
    vertices = [
        s2sphere.LatLng.from_degrees(lat, lon)
        for lon, lat in polygon.exterior.coords
    ]
    #print('vertices='+str(vertices))
    
    # Use S2RegionCoverer to find covering cells
    region_coverer = s2sphere.RegionCoverer()
    region_coverer.min_level = level
    region_coverer.max_level = level

    #print(vertices[0], vertices[2])
    # Cover the bounding box of the polygon
    rect = s2sphere.LatLngRect.from_point_pair(vertices[0], vertices[2])
    cell_ids = region_coverer.get_covering(rect)
    #print(cell_ids)
    
    # Collect S2 tokens and geometries
    rows = []
    for cell_id in cell_ids:
        #print('level='+str(cell_id.level()))
        
        cell = s2sphere.Cell(cell_id)
        cell_vertices = [
            s2sphere.LatLng.from_point(cell.get_vertex(i))
            for i in range(4)
        ]
        coords = [
            (v.lng().degrees, v.lat().degrees)
            for v in cell_vertices
        ]
        coords.append(coords[0])  # Close the polygon
        
        # Create a Shapely Polygon
        cell_polygon = wkt.loads(
            f"POLYGON(({','.join([f'{lon} {lat}' for lon, lat in coords])}))"
        )
        
        # Append the S2 token and geometry
        try:
            rows.append({
                's2_index__L3': cell_id.parent(3).to_token(),
                's2_index__L8': cell_id.parent(8).to_token(),
                's2_index__L13': cell_id.parent(13).to_token(),
                # 's2_index__L10': cell_id.parent(10).to_token(),  # This can fail if level is too low
                's2_index__L'+str(level): cell_id.to_token(),
                'geometry': cell_polygon
            })
        except AssertionError:
            # Handle the error gracefully, maybe skip this row or log a warning
            pass


    # Create a GeoDataFrame for the Polygon
    field_gdf = gpd.GeoDataFrame({'geometry': [polygon]})
    
    # Add additional columns
    field_gdf['geoid'] = [geoid]
    
    # Set a CRS (Coordinate Reference System) if needed
    field_gdf.set_crs(epsg=4326, inplace=True)  # WGS84
    
    # Create and return both the GeoDataFrames
    gdf = gpd.GeoDataFrame(rows, columns=['s2_index__L3','s2_index__L8','s2_index__L13','s2_index__L'+str(level), 'geometry'])

    #convert to the same CRS 4326
    field_gdf = field_gdf.set_crs(epsg=4326)
    gdf = gdf.set_crs(epsg=4326)
    
    return gdf, field_gdf
