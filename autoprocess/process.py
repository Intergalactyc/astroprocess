import sys
import os
import glob
import re
import warnings
import gc

from pathlib import PurePath
from tqdm import tqdm

from astropy import units as u
from astropy.nddata import CCDData
from astropy.stats import mad_std
from astropy import log
from astropy.utils.exceptions import AstropyWarning

import ccdproc as ccdp
import numpy as np

warnings.simplefilter('ignore', category=AstropyWarning)
log.setLevel('WARNING')

ramlimit = 70e8
sclipthresh = 4

def msigcombine(frames): return ccdp.combine(frames,method="average",dtype="float32",sigma_clip=True,sigma_clip_low_thresh=sclipthresh,sigma_clip_high_thresh=sclipthresh,sigma_clip_func=np.ma.median,sigma_clip_dev_func=mad_std,mem_limit=ramlimit,unit="adu")
def medcombine(frames): return ccdp.combine(frames,method="median",dtype="float32",mem_limit=ramlimit,unit="adu")

def clearfolder(folder):
    for filename in os.listdir(folder):
        file_path = os.path.join(folder, filename)
        try:
            os.unlink(file_path)
        except Exception as e:
            print("Failed to delete %s. Reason: %s" % (file_path,e))

args = sys.argv
directory = ""
simplecombine = False
toplevel = False
useexisting = False

if len(args) == 2:
    directory = args[1].rstrip("/")
elif len(args) == 3:
    directory = args[1].rstrip("/")
    if "simple" in args[2] or "median" in args[2]: simplecombine = True
    if "top" in args[2]: toplevel = True
    if "use" in args[2] or "exist" in args[2]: useexisting = True
else:
    print("Invalid number of arguments. Pass the directory to process, as well as optionally the value of the boolean simplecombine.")
    quit()

subfolders = []
subpipes = ["masters","intermediate","calibrated_science","calibrated_standards","aligned_science","aligned_standards"]

if toplevel:
    subfolders.append("!TOPLEVEL!")
else:
    scan = os.scandir(directory)
    for item in scan: 
        if item.is_dir(): subfolders.append(item.name)

for folder in subfolders:
    if folder == "!TOPLEVEL!":
        base_path = directory
    else:
        base_path = os.path.join(directory,folder)
    
    out_path = os.path.join(base_path,"pipelineout")

    for sp in subpipes:
        path = os.path.join(out_path,sp)
        os.makedirs(path, exist_ok = True)

    files = glob.glob(base_path+"/**/*.fit", recursive=True) # FITS files as .fit
    files.extend(glob.glob(base_path+"/**/*.fits",recursive=True)) # FITS files as .fits

    biases = [] # list
    darks = [] # list
    flats = {} # dictionary of lists
    science = {} # dictionary of dictionaries of lists
    standards = {} # dictionary of dictionaries of lists
    masters = {"master_bias": None, "master_dark": None}

    print("Sorting data from " + base_path)

    for file in tqdm(files, desc="Sorting data from " + base_path):
        filepath = PurePath(file)
        classified = False
        if useexisting:
            if "master_bias" in filepath.name:
                masters["master_bias"] = file
                classified = True
            elif "master_dark" in filepath.name:
                masters["master_dark"] = file
                classified = True
            elif "master_flat" in filepath.name:
                try:
                    masters[filepath.name.split(".")[0]] = file
                except:
                    print("Master flat found with name " + filepath.name +" but could not be classified.")
                classified = True
        if not (classified or "pipelineout" in file):
            if "bias" in filepath.parent.name.lower():
                biases.append(file)
            elif "dark" in filepath.parent.name.lower():
                darks.append(file)
            elif "flat" in filepath.parent.parent.name.lower():
                filt = filepath.parent.name
                if filt not in flats.keys(): flats[filt]=[]
                flats[filt].append(file)
            elif "science" in filepath.parent.parent.name.lower() or "light" in filepath.parent.parent.name.lower(): # If arranged as science/filter/*
                filt = filepath.parent.name
                if "obj" not in science.keys(): science["obj"] = {}
                if filt not in science["obj"].keys(): science["obj"][filt]=[]
                science["obj"][filt].append(file)
            elif "science" in filepath.parent.parent.parent.name.lower() or "light" in filepath.parent.parent.parent.name.lower(): # If arranged as science/obj/filter/*
                obj = filepath.parent.parent.name
                filt = filepath.parent.name
                if obj not in science.keys(): science[obj] = {}
                if filt not in science[obj].keys(): science[obj][filt]=[]
                science[obj][filt].append(file)
            elif "standard" in filepath.parent.parent.parent.name.lower(): # must be arranged as standards/sobj/filter/*
                obj = filepath.parent.parent.name
                filt = filepath.parent.name
                if obj not in standards.keys(): standards[obj] = {}
                if filt not in standards[obj].keys(): standards[obj][filt]=[]
                science[obj][filt].append(file)
            else:
                print(file+" could not be classified.")

    print("Data sorted.")

    if masters.get("master_bias",None) is not None:
        master_bias = CCDData.read(masters.get("master_bias"),unit="adu")
        print("Master bias already exists for " + base_path)
    else:
        print("Forming master bias for " + base_path)
        if simplecombine:
            master_bias = medcombine(biases)
        else:
            master_bias = msigcombine(biases)
        master_bias.meta["combined"] = True
        print(str(len(biases)) + " biases combined to form " + out_path + "\\masters\\master_bias.fit")
    master_bias.data = master_bias.data.astype("float32")
    master_bias.write(out_path + "\\masters\\master_bias.fit",overwrite=True)
    master_bias = CCDData.read(out_path + "\\masters\\master_bias.fit",unit="adu")

    if masters.get("master_dark",None) is not None:
        master_dark = CCDData.read(masters.get("master_dark"),unit="adu")
        print("Master dark already exists for " + base_path)
    else:
        for dark in tqdm(darks, desc = "Debiasing darks for " + base_path): # Need to try considering different exposure times and making master darks for each exposure time if there are multiple
            name = PurePath(dark).name
            ccd = CCDData.read(dark,unit="adu")
            ccd = ccdp.subtract_bias(ccd,master_bias)
            ccd.data = ccd.data.astype("float32")
            ccd.write(out_path + "\\intermediate\\debiasdark_"+name+".fit",overwrite=True)
        print(str(len(darks)) + " darks debiased.")
        print("Forming master dark for " + base_path)
        ddarks = glob.glob(out_path + "\\intermediate\\debiasdark_*.fit")
        if simplecombine:
            master_dark = medcombine(ddarks)
        else:
            master_dark = msigcombine(ddarks)
        master_dark.meta["combined"] = True
        print(str(len(ddarks)) + " darks combined to form " + out_path + "\\masters\\master_dark.fit")
        clearfolder(out_path + "\\intermediate")
    master_dark.data = master_dark.data.astype("float32")
    master_dark.write(out_path + "\\masters\\master_dark.fit",overwrite=True)
    master_dark = CCDData.read(out_path + "\\masters\\master_dark.fit",unit="adu")

    for filt in flats.keys():
        if masters.get("master_flat_"+filt,None) is not None:
            master_flat = CCDData.read(masters.get("master_flat_"+filt),unit="adu")
            print("Master flat already exists for filter " + filt + " in " + base_path)
        else:
            for flat in tqdm(flats[filt], desc = "Calibrating flats for filter " + filt + " in " + base_path):
                name = PurePath(flat).name
                ccd = CCDData.read(flat,unit="adu")
                ccd.data = ccd.data.astype("float32")
                ccd = ccdp.subtract_bias(ccd,master_bias)
                ccd = ccdp.subtract_dark(ccd,master_dark,exposure_time="EXPTIME",exposure_unit=u.second,scale=True)
                ccd.write(out_path + "\\intermediate\\calibratedflat-"+filt+"_"+name+".fit",overwrite=True)
            cflats = glob.glob(out_path + "\\intermediate\\calibratedflat-"+filt+"_*.fit")
            print(str(len(cflats)) + " flats calibrated.")
            print("Forming master flat for filter " + filt + " in " + base_path)
            if simplecombine:
                master_flat = medcombine(cflats)
            else:
                master_flat = msigcombine(cflats)
            master_flat.meta["combined"] = True
            master_flat.write(out_path + "\\masters\\master_flat_"+filt+".fit",overwrite=True)
            print(str(len(cflats)) + " flats combined to form " + out_path + "\\masters\\master_flat_"+filt+".fit")
            clearfolder(out_path + "\\intermediate")

    for obj in science.keys():
        for filt in science[obj].keys():
            path = os.path.join(out_path,"calibrated_science",obj,filt)
            os.makedirs(path, exist_ok = True)
            master_flat = CCDData.read(out_path + "\\masters\\master_flat_"+filt+".fit", unit="adu")
            master_flat.data = master_flat.data.astype("float32")
            for frame in tqdm(science[obj][filt], desc="Calibrating frames of " + obj + " for filter " + filt + " in " + base_path):
                name = PurePath(frame).name
                ccd = CCDData.read(frame,unit="adu")
                ccd.data = ccd.data.astype("float32")
                ccd = ccdp.subtract_bias(ccd,master_bias)
                ccd = ccdp.subtract_dark(ccd,master_dark,exposure_time="EXPTIME",exposure_unit=u.second,scale=True)
                ccd = ccdp.flat_correct(ccd,master_flat)
                ccd.uncertainty = None
                ccd.mask = None
                ccd.write(path + "\\calibrated_"+name+".fit",overwrite=True)
            print(str(len(science[obj][filt])) + " frames calibrated.")

    for obj in standards.keys():
        for filt in standards[obj].keys():
            print("Calibrating frames of " + obj + " for filter " + filt + " in " + base_path)
            path = os.path.join(out_path,"calibrated_standards",obj,filt)
            os.makedirs(path, exist_ok = True)
            master_flat = CCDData.read(out_path + "\\masters\\master_flat_"+filt+".fit", unit="adu")
            for frame in standards[obj][filt]:
                name = PurePath(frame).name
                ccd = CCDData.read(frame,unit="adu")
                ccd.data = ccd.data.astype("float32")
                ccd = ccdp.subtract_bias(ccd,master_bias)
                ccd = ccdp.subtract_dark(ccd,master_dark,exposure_time="EXPTIME",exposure_unit=u.second,scale=True)
                ccd = ccdp.flat_correct(ccd,master_flat)
                ccd.write(path + "\\calibrated_"+name+".fit",overwrite=True)
            print(str(len(standards[obj][filt])) + " frames calibrated.")
    
    # delete "intermediate" subdir
    # are all of the castings good? is there any way to tell it to use this dtype during intermediate operations?
        # also can I make it so that the uncertainty and mask frames aren't even formed in the first place, rather than having to delete them at the end?
    # create a loading bar
        # how to measure for the combine operations...?
    # time estimate? (could go along with loading bar)
        # rough note: a run over testdata (~5 GB) took about 13 minutes, which is an average processing rate of ~6-7 MBps
            # how does AIJ do it so fing fast?
    # control RAM dedication and stuff?
        # is there a way to do this over the internet from a dedicated machine somewhere?
    # generally look at how I can reduce the processing time / optimize this
        # what if I just convert everything into numpy arrays and have numpy zoom through things?
            # look into fitsio library
    # comment your damn code, me
    # add a useexisting option: if selected, then if a master frame already exists, don't create a new one and just use what's already there

    # future:
        # artifact detection
        # automated alignment
            # RA/Dec assignment