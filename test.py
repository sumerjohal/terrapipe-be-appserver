from flask import Flask, request, jsonify
import os
import json
import pandas as pd
import pyarrow.parquet as pq
from shapely import wkt
import time
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler
import threading
import s2sphere as s2
import pyarrow as pa


app = Flask(__name__)

base_path = '/mnt/md1/NCEP/PARQUET_S2/'
configFile = '/home/rnaura/terrapipe/config/11SKA/11SKA.json'
last_update_time = None
cached_weather_df = pd.DataFrame()  # Placeholder for cached data
last_reload_time = None

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

def configure_reload_on_write(reload_function):
    event_handler = FileChangeHandler(reload_function)
    observer = Observer()
    observer.schedule(event_handler, base_path, recursive=True)
    observer.start()
    print(f"Started monitoring {base_path} for changes")

    # Ensure observer runs indefinitely in the background
    try:
        while True:
            observer.join(1)
    except KeyboardInterrupt:
        observer.stop()
    observer.join()

    
def reload_data(YYYY_str, MM_str, DD_str):
    global cached_weather_df, last_reload_time
    
    try:
        start_time = time.time()
        paths = "/mnt/md1/NCEP/PARQUET_S2/s2_tokens_L5/s2_tokens_L7/s2_tokens_L9/Year=2024/Month=06/Day=25/"
        all_files = [os.path.join(paths, f) for f in os.listdir(paths) if f.endswith('.parquet')]
        
        tables = [pq.read_table(file) for file in all_files]
        combined_table = pa.concat_tables(tables)
        cached_weather_df = combined_table.to_pandas()
        last_reload_time = time.time()
        
        print(f"Data reloaded in {time.time() - start_time:.4f} seconds")
        
    except Exception as e:
        print(f"Error reloading data: {e}")

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

@app.route('/getDataFromNLDAS', methods=['GET'])
def get_weather_from_nldas():
    global cached_weather_df, last_reload_time
    
    start_time_total = time.time()

    agstack_geoid = request.args.get('geoid')
    dtStr = request.args.get('date')
    
    try:
        tok = dtStr.split('-')
        YYYY_str = tok[0]
        MM_str = tok[1].zfill(2)
        DD_str = tok[2]
        
        # Check if data needs reloading based on last_reload_time
        if last_reload_time is None or time.time() - last_reload_time > 60 * 5:  # Reload every 5 minutes
            reload_data(YYYY_str, MM_str, DD_str)
        
        with open(configFile, "r") as jsonfile:
            fieldJSON = json.load(jsonfile)
            field_geoid = fieldJSON['geoid']
            field_wkt = fieldJSON['wkt']

        if agstack_geoid != field_geoid:
            raise ValueError("Geoid does not match the geoid in the JSON file.")
        
        fieldPoly = wkt.loads(field_wkt)
        c = fieldPoly.centroid
        lat, lon = c.y, c.x
        
        lats, lons = [lat], [lon]
        s2_index__L9_list, _ = get_s2_cellids_and_token_list(9, lats, lons)
        
        # print(cached_weather_df.head())  # Debugging: Print the head of the cached DataFrame
        
        # Filter DataFrame by s2_token_L9 values
        if 's2_token_L9' in cached_weather_df.columns:
            filtered_df = cached_weather_df[cached_weather_df['s2_token_L9'].isin(s2_index__L9_list)]
        else:
            filtered_df = pd.DataFrame()

        # Calculate total request handling time
        end_time_total = time.time()
        total_duration = end_time_total - start_time_total
        
        # Print durations to screen
        print(f"Total request handling time: {total_duration:.4f} seconds")
        
        return jsonify(filtered_df.to_dict(orient='records'))
    
    except Exception as e:
        print(f"Error handling request: {e}")
        return jsonify({"error": str(e)}), 500
    

if __name__ == '__main__':
    start_file_monitoring(base_path, reload_data)
    app.run(debug=True)
