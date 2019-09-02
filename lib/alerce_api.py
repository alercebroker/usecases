import requests
import pandas as pd
from pandas.io.json import json_normalize
from IPython.display import HTML
#from xml.etree import ElementTree
#import numpy as np
#from io import BytesIO
#from PIL import Image
#import base64


class alerce_api(object):
    'ALeRCE API python wrapper class'
    
    # init
    def __init__(self, **kwargs):
                
        self.ztf_url = "http://ztf.alerce.online"
        if "ztf_url" in kwargs.keys():
            self.ztf_url = kwargs["ztf_url"]
            
        self.catsHTM_url = "http://catshtm.alerce.online:5000"
        if "catsHTM_url" in kwargs.keys():
            self.catsHTM_url = kwargs["catsHTM_url"]
            
        self.oid = ""


    def query(self, params):
        'do general query given json query parameters'
        
        # show api results
        r = requests.post(url = "%s/query" % self.ztf_url, json = params) 
        df = pd.DataFrame(r.json())
        query_results = json_normalize(df.result)
        query_results.set_index('oid', inplace=True)
        return query_results

    def get_sql(self, params):
        'get sql from given json query parameters'
        
        r = requests.post(url = "%s/get_sql" % self.ztf_url, json = params)
        return r.content

    def get_detections(self, oid):
        'get detections given oid as pandas dataframe'
        
        #oid
        params = {
            "oid": oid
        }
        
        # show api results
        r = requests.post(url = "%s/get_detections" % self.ztf_url, json = params) 
        df = pd.DataFrame(r.json())
        detections = json_normalize(df.result.detections)
        detections.set_index('candid', inplace=True)
        return detections

    
    def get_non_detections(self, oid):
        '# get non detections given oid as pandas  dataframe'
        
        #oid
        params = {
            "oid": oid
        }
        
        # show api results
        r = requests.post(url = "%s/get_non_detections" % self.ztf_url, json = params) 
        df = pd.DataFrame(r.json())
        non_detections = json_normalize(df.result.non_detections)
        if non_detections.shape[0] > 0:
            non_detections.set_index('mjd', inplace=True)
        return non_detections

    def get_stats(self, oid, verbose=False):
        'get stats given oid as pandas dataframe'
        
        #oid
        params = {
            "oid": oid
        }
        
        # show api results
        r = requests.post(url = "%s/get_stats" % self.ztf_url, json = params) 
        df = pd.DataFrame(r.json())
        stats = json_normalize(df.result.stats)
        stats.set_index('oid', inplace=True)
        
        # save some features
        if verbose:
            print("Setting ra, dec, firstMJD")
        self.oid = oid
        self.meanra = stats.meanra
        self.meandec = stats.meandec
        self.firstMJD = stats.firstmjd
        
        return stats

    def get_probabilities(self, oid, dolate=True):
        'get probabilities given oid as pandas dataframe (late or early)'
        
        #oid
        params = {
            "oid": oid
        }
        
        # show api results
        r = requests.post(url = "%s/get_probabilities" % self.ztf_url, json = params) 
        early = json_normalize(r.json()["result"]["probabilities"]["early_classifier"])
        early.set_index("oid", inplace=True)
        if dolate:
            late = json_normalize(r.json()["result"]["probabilities"]["random_forest"])
            late.set_index("oid", inplace=True)
            return early, late
        else:
            return early

    def get_features(self, oid):
        'get features given oid as pandas dataframe'
    
        #oid
        params = {
            "oid": oid
        }
        
        # show api results
        r = requests.post(url = "%s/get_features" % self.ztf_url, json = params) 
        features = json_normalize(r.json())
        features.set_index('oid', inplace=True)
        return features
    
    def catsHTM(self, oid, catalog_name, radius):
        'get catsHTM crossmatch given oid, catalog_name (all is also allowed) and search radius in arcsec'
        
        # get ra, dec
        if oid != self.oid:
            self.get_stats(oid)
        
        params = {
            "catalog": "%s" % catalog_name,
            "ra": "%f" % self.meanra,
            "dec": "%f" % self.meandec,
            "radius": "%f" % radius
        }
        if catalog_name != "all":
            r = requests.get(url = "%s/crossmatch" % self.catsHTM_url, params = params)
        else:
            r = requests.get(url = "%s/crossmatch_all" % self.catsHTM_url, params = params)
        return r.json()

    def catsHTM_redshift(self, oid, radius):
        'get redshift given oid using catsHTM crossmatch'

        # get ra, dec
        if oid != self.oid:
            self.get_stats(oid)
                
        # search data in catsHTM
        xmatches = self.catsHTM(oid, "all", radius)
        
        # check whether redshift exists
        for xmatch in xmatches:
            for catname in xmatch.keys():
                for col in xmatch[catname].keys():
                    if col == 'z':
                        print("Found redshift in catalog %s: %f" % (catname, xmatch[catname][col]["value"]))
                        return float(xmatch[catname][col]["value"])

        # check whether photometric redshift exists
        for xmatch in xmatches:
            for catname in xmatch.keys():
                for col in xmatch[catname].keys():
                    if col in ['zphot',  'z_phot']:
                        print("Found photometric redshift in catalog %s: %f" % (catname, xmatch[catname][col]["value"]))
                        return float(xmatch[catname][col]["value"])
                    
        return

    def plot_stamp(self, oid, candid=None):
        'plot stamp in a notebook given oid. It uses IPython HTML.'
        
        # if candid is None, get minimum candid
        if candid is None:
            candid = min(self.get_detections(oid).index)
            
        print(oid, candid)
        science = "http://avro.alerce.online/get_stamp?oid=%s&candid=%s&type=science&format=png" % (oid, candid)
        images="""
        &emsp;&emsp;&emsp;&emsp;&emsp;
        Science
        &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp; 
        Template
        &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp; 
        Difference
        <div class="container">
        <div style="float:left;width:20%%"><img src="%s"></div>
        <div style="float:left;width:20%%"><img src="%s"></div>
        <div style="float:left;width:20%%"><img src="%s"></div>
        </div>
        """ % (science, science.replace("science", "template"), science.replace("science", "difference"))
        display(HTML(images))
