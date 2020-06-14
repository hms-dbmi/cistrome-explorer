import requests
import yaml
from pprint import pprint

CISTROME_DB_INDEX_URL = "http://dc2.cistrome.org/api/main_filter_ng?allqc=false&cellinfos=all&completed=false&curated=false&factors=all&keyword=&page=1&peakqc=false&run=false&species=all"

def generate_config(species=None, biological_sources=None, factors=None):



if __name__ == "__main__":

    config = generate_config()
    with open('config.yml', 'w') as f:
        yaml.dump(config, f, default_flow_style=False)