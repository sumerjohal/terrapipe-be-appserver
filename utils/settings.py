from utils.imports import *

ST_TONER = 'http://tile.stamen.com/toner/{z}/{x}/{y}.png'
ST_TONER_HYBRID = 'http://tile.stamen.com/toner-hybrid/{z}/{x}/{y}.png'
ST_TONER_LABELS = 'http://tile.stamen.com/toner-labels/{z}/{x}/{y}.png'
ST_ST_TONER_LINES = 'http://tile.stamen.com/toner-lines/{z}/{x}/{y}.png'
ST_TONER_BACKGROUND = 'http://tile.stamen.com/toner-background/{z}/{x}/{y}.png'
ST_TONER_LITE = 'http://tile.stamen.com/toner-lite/{z}/{x}/{y}.png'

ST_TERRAIN = 'http://tile.stamen.com/terrain/{z}/{x}/{y}.png'
ST_TERRAIN_LABELS = 'http://tile.stamen.com/terrain-labels/{z}/{x}/{y}.png'
ST_TERRAIN_LINES = 'http://tile.stamen.com/terrain-lines/{z}/{x}/{y}.png'
ST_TERRAIN_BACKGROUND = 'http://tile.stamen.com/terrain-background/{z}/{x}/{y}.png'

URL_BACKGROUND = ctx.sources.ST_TERRAIN
URL_BACKGROUND2 = ctx.sources.ST_TERRAIN_BACKGROUND
URL_LINES = ctx.sources.ST_TONER_LINES
URL_LABELS = ctx.sources.ST_TERRAIN_LABELS

#colormap customization - 1
top = []
bottom=[]
top = cm.get_cmap('Oranges_r', 64)
bottom = cm.get_cmap('Blues', 256)
newcolors = np.vstack((top(np.linspace(0, 1, 64)),bottom(np.linspace(0, 1, 256))))
zoom1 = ListedColormap(newcolors, name='Zoom1')

#colormap customization - 2
t2 = []
b2 =[]
t2 = cm.get_cmap('Reds_r', 128)
b2 = cm.get_cmap('Blues', 128)
newcolors2 = np.vstack((t2(np.linspace(0, 1, 128)),b2(np.linspace(0, 1, 128))))
zoom2 = ListedColormap(newcolors2, name='Zoom2')

#colormap customization - 3
t3 = []
b3=[]
t3 = cm.get_cmap('Reds_r', 128)
b3 = cm.get_cmap('Greens', 128)
newcolors3 = np.vstack((t3(np.linspace(0, 1, 128)),b3(np.linspace(0, 1, 128))))
zoom3 = ListedColormap(newcolors3, name='Zoom3')

#colormap customization - 5
top = []
bottom=[]
top = cm.get_cmap('Oranges_r', 64)
bottom = cm.get_cmap('Greens', 256)
newcolors = np.vstack((top(np.linspace(0, 1, 64)),bottom(np.linspace(0, 1, 256))))
zoom5 = ListedColormap(newcolors, name='Zoom5')

#********* CHANGE FOR VISHNU
# data_dir = '/home/rnaura/terrapipe/shapefiles/'
data_dir = '/home/rajat/Downloads/Rnaura_Work/terrapipe-dataserver/shapefiles/'
isr_shapefile = data_dir + 'ISR_adm0.shp'

aus_data_dir = '/home/rnaura/terrapipe/shapefiles/AUS_ADMIN_BOUNDARIES/State Boundaries FEBRUARY 2020/Standard/'
worldShpFile = '/home/rnaura/terrapipe/data/99bfd9e7-bb42-4728-87b5-07f8c8ac631c2020328-1-1vef4ev.lu5nk.shp'


outDir = '/mnt/md0/CONFIG/JAINLOGIC_MONITORED/'
tmpDir = '/home/sumer/my_project_dir/tmp/'

sentinelRawFolder = '/mnt/md1/SENTINEL/RAW/'
# sentinelProcFolder = '/mnt/md0/SENTINEL/PROCESSED/'
sentinelProcFolder = '/home/rajat/Downloads/Rnaura_Work/mnt/md0/SENTINEL/PROCESSED/'

sentinelProcFolder_NIROG = '/mnt/md1/SENTINEL/NIROG/'
sentinel_MULTIBAND_ProcFolder = '/mnt/md1/SENTINEL/MULTIBAND/'
sentinel_EXTENSIONS_ProcFolder = '/mnt/md0/SENTINEL_EXTENSIONS/PROCESSED'


landsatRawFolder = '/mnt/md1/LANDSAT8/RAW/'
# landsatProcFolder = '/mnt/md1/LANDSAT8/PROCESSED/'
landsatProcFolder = '/home/rajat/Downloads/Rnaura_Work/mnt/md1/LANDSAT8/PROCESSED/'

landsatProcFolder_NIROG = '/mnt/md1/LANDSAT8/NIROG/'

configPath = '/mnt/md0/CONFIG/'

cimisRawPath = '/mnt/md0/CIMIS/RAW/'
cimisProcPath = '/mnt/md0/CIMIS/PROCESSED/'
cimisConfigPath = configPath + 'CIMIS/'
cimis_monthly_ftp_loc = 'ftp://ftpcimis.water.ca.gov/pub2/monthly/'


acctConfigPath = configPath + 'JAINLOGIC_MONITORED/'
jain_acctConfigPath = configPath + 'JAINLOGIC_MONITORED/'
agls_acctConfigPath = configPath + 'AGLS_ACCOUNTS/'
pilotAcctConfigPath = configPath + 'PILOT/'

acctOutputPath = '/mnt/md0/OUTPUT/JAINLOGIC_MONITORED/'
jain_acctOutputPath = '/mnt/md0/OUTPUT/JLV2/JAINLOGIC_MONITORED/'
agls_acctOutputPath = '/mnt/md0/OUTPUT/AGLS_ACCOUNTS/'
pilotAcctOutputPath = '/mnt/md0/OUTPUT/PILOT/'


NOAA_StnRoot_tmp = '/mnt/md0/NOAA/TMP/'

NOAA_ProcessedRoot = '/mnt/md0/NOAA/DAILY/PROCESSED/'

NOAA_ForecastRoot = '/mnt/md0/NOAA/DAILY/FORECASTED/'
NOAA_ForecastConfigPath = '/mnt/md0/CONFIG/NOAA/FORECAST/'

acctForecastOutputPath = '/mnt/md0/OUTPUT/FORECASTED_JAINLOGIC_MONITORED/'
agls_acctForecastOutputPath = '/mnt/md0/OUTPUT/FORECASTED_AGLS_ACCOUNTS/'

AUS_ForecastPath = '/mnt/md0/AUS/FORECAST/'
IND_ForecastPath = '/mnt/md0/IND/FORECAST/'
MEX_ForecastPath = '/mnt/md0/MEX/FORECAST/'

WA_ForecastPath = NOAA_ForecastRoot + 'WA/'
CA_ForecastPath = NOAA_ForecastRoot + 'CA/'
IA_ForecastPath = NOAA_ForecastRoot + 'IA/'

IMS_OutputPath_hrly = '/mnt/md0/IMS/DAILY/RAW/'
IMS_OutputPath_dly = '/mnt/md0/IMS/DAILY/PROCESSED/'

noaa_parquet_dir = '/mnt/md1/NOAA/PARQUET/'

#S2 Paths for weather
nldas_filePath = '/mnt/md1/NLDAS/PARQUET_S2/'
noaa_filePath = '/mnt/md1/NOAA/PARQUET_S2/'
ncep_filePath = '/mnt/md1/NCEP/PARQUET_S2/'


configRoot = acctConfigPath
acctOutputRoot = acctOutputPath
acctReportRoot = acctOutputPath
jain_acctOutputRoot = jain_acctOutputPath
agls_acctReportRoot = agls_acctOutputPath
jain_acctReportRoot = jain_acctOutputPath

stnProcPath = '/mnt/md0/'
aus_obsPath = stnProcPath + 'AUS/'
aus_stnFile = stnProcPath +'AUS/AUS_STNS.geojson'

aus_obsPath_new = stnProcPath + 'GHCND/DAILY/PROCESSED/'
aus_stnFile_new = stnProcPath +'GHCND/GHCND_STNS__AUS.geojson'

ind_obsPath = stnProcPath + 'GHCND/DAILY/PROCESSED/'
ind_stnFile = stnProcPath +'GHCND/GHCND_STNS__IND.geojson'

isr_obsPath = stnProcPath + 'GHCND/DAILY/PROCESSED/'
isr_stnFile = stnProcPath +'GHCND/GHCND_STNS__ISR.geojson'

mex_obsPath = stnProcPath + 'GHCND/DAILY/PROCESSED/'
mex_stnFile = stnProcPath +'GHCND/GHCND_STNS__MEX.geojson'

egy_obsPath = stnProcPath + 'GHCND/DAILY/PROCESSED/'
egy_stnFile = stnProcPath +'GHCND/GHCND_STNS__EGY.geojson'

can_obsPath = stnProcPath + 'GHCND/DAILY/PROCESSED/'
can_stnFile = stnProcPath +'GHCND/GHCND_STNS__CAN.geojson'

pak_obsPath = stnProcPath + 'GHCND/DAILY/PROCESSED/'
pak_stnFile = stnProcPath +'GHCND/GHCND_STNS__PAK.geojson'

noaa_stnProcPath = stnProcPath + 'NOAA/DAILY/PROCESSED/'

az_obsPath = noaa_stnProcPath + 'AZ/'
az_stnFile = stnProcPath +'NOAA/NOAA_STNS__AZ.geojson'

id_obsPath = noaa_stnProcPath + 'ID/'
id_stnFile = stnProcPath +'NOAA/NOAA_STNS__ID.geojson'

ia_obsPath = noaa_stnProcPath + 'IA/'
ia_stnFile = stnProcPath +'NOAA/NOAA_STNS__IA.geojson'

ca_obsPath = noaa_stnProcPath + 'CA/'
ca_stnFile = stnProcPath +'NOAA/NOAA_STNS__CA.geojson'

wa_obsPath = noaa_stnProcPath + 'WA/'
wa_stnFile = stnProcPath +'NOAA/NOAA_STNS__WA.geojson'

fl_obsPath = noaa_stnProcPath + 'FL/'
fl_stnFile = stnProcPath +'NOAA/NOAA_STNS__FL.geojson'

tx_obsPath = noaa_stnProcPath + 'TX/'
tx_stnFile = stnProcPath +'NOAA/NOAA_STNS__TX.geojson'

or_obsPath = noaa_stnProcPath + 'OR/'
or_stnFile = stnProcPath +'NOAA/NOAA_STNS__OR.geojson'

ok_obsPath = noaa_stnProcPath + 'OK/'
ok_stnFile = stnProcPath +'NOAA/NOAA_STNS__OK.geojson'

mt_obsPath = noaa_stnProcPath + 'MT/'
mt_stnFile = stnProcPath +'NOAA/NOAA_STNS__MT.geojson'

mo_obsPath = noaa_stnProcPath + 'MO/'
mo_stnFile = stnProcPath +'NOAA/NOAA_STNS__MO.geojson'

de_obsPath = noaa_stnProcPath + 'DE/'
de_stnFile = stnProcPath +'NOAA/NOAA_STNS__DE.geojson'

ar_obsPath = noaa_stnProcPath + 'AR/'
ar_stnFile = stnProcPath +'NOAA/NOAA_STNS__AR.geojson'


logo_root = '/home/sumer/my_project_dir/logos/'
AGLS_logoFile = logo_root + 'AgralogicsLogo.png'
JAIN_logoFile = logo_root + 'jainlogo.png'
PAG_logoFile = logo_root + 'PrincipleAg_Logo_FullColor.png'

shiva_acctConfigPath = '/media/sumer/E6/CONFIG/'
acctMapGetID = dict({'ARBUSTO': 1513, 'THURLA':1556, 'JALISCO2':1207, 'CAVAOUR': 1624})
acctMapGetName = dict({1513:'ARBUSTO', 1556:'THURLA', 1207:'JALISCO2', 'CAVAOUR': 1624})


errLogFileRoot = '/home/sumer/my_project_dir/err/'


bounds = np.array([-50,-40, -30, -20, -10, 0, 10, 20, 30, 40, 50])
norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)

superscript_map = {
    "0": "⁰", "1": "¹", "2": "²", "3": "³", "4": "⁴", "5": "⁵", "6": "⁶",
    "7": "⁷", "8": "⁸", "9": "⁹", "a": "ᵃ", "b": "ᵇ", "c": "ᶜ", "d": "ᵈ",
    "e": "ᵉ", "f": "ᶠ", "g": "ᵍ", "h": "ʰ", "i": "ᶦ", "j": "ʲ", "k": "ᵏ",
    "l": "ˡ", "m": "ᵐ", "n": "ⁿ", "o": "ᵒ", "p": "ᵖ", "q": "۹", "r": "ʳ",
    "s": "ˢ", "t": "ᵗ", "u": "ᵘ", "v": "ᵛ", "w": "ʷ", "x": "ˣ", "y": "ʸ",
    "z": "ᶻ", "A": "ᴬ", "B": "ᴮ", "C": "ᶜ", "D": "ᴰ", "E": "ᴱ", "F": "ᶠ",
    "G": "ᴳ", "H": "ᴴ", "I": "ᴵ", "J": "ᴶ", "K": "ᴷ", "L": "ᴸ", "M": "ᴹ",
    "N": "ᴺ", "O": "ᴼ", "P": "ᴾ", "Q": "Q", "R": "ᴿ", "S": "ˢ", "T": "ᵀ",
    "U": "ᵁ", "V": "ⱽ", "W": "ᵂ", "X": "ˣ", "Y": "ʸ", "Z": "ᶻ", "+": "⁺",
    "-": "⁻", "=": "⁼", "(": "⁽", ")": "⁾"}

subscript_map = {
    "0": "₀", "1": "₁", "2": "₂", "3": "₃", "4": "₄", "5": "₅", "6": "₆",
    "7": "₇", "8": "₈", "9": "₉", "a": "ₐ", "b": "♭", "c": "꜀", "d": "ᑯ",
    "e": "ₑ", "f": "բ", "g": "₉", "h": "ₕ", "i": "ᵢ", "j": "ⱼ", "k": "ₖ",
    "l": "ₗ", "m": "ₘ", "n": "ₙ", "o": "ₒ", "p": "ₚ", "q": "૧", "r": "ᵣ",
    "s": "ₛ", "t": "ₜ", "u": "ᵤ", "v": "ᵥ", "w": "w", "x": "ₓ", "y": "ᵧ",
    "z": "₂", "A": "ₐ", "B": "₈", "C": "C", "D": "D", "E": "ₑ", "F": "բ",
    "G": "G", "H": "ₕ", "I": "ᵢ", "J": "ⱼ", "K": "ₖ", "L": "ₗ", "M": "ₘ",
    "N": "ₙ", "O": "ₒ", "P": "ₚ", "Q": "Q", "R": "ᵣ", "S": "ₛ", "T": "ₜ",
    "U": "ᵤ", "V": "ᵥ", "W": "w", "X": "ₓ", "Y": "ᵧ", "Z": "Z", "+": "₊",
    "-": "₋", "=": "₌", "(": "₍", ")": "₎"}
trans = str.maketrans(
    ''.join(superscript_map.keys()),
    ''.join(superscript_map.values()))
sub_trans = str.maketrans(
    ''.join(subscript_map.keys()),
    ''.join(subscript_map.values()))