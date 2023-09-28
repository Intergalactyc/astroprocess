import sys
import os

args = sys.argv
directory = ""
name = ""
objs = []
filts = []
standards = []
if len(args) == 3:
    name = args[1]
    directory = args[2]
    objinput = input("Input names of objects, separated by commas: ")
    objsraw = objinput.split(",")
    for obj in objsraw:
        objs.append(obj.lstrip().rstrip())
    filtinput = input("Input filters to use, separated by commas: ")
    filtsraw = filtinput.split(",")
    for filt in filtsraw:
        filts.append(filt.lstrip().rstrip())
    standardinput = input("If using standard stars, input standard stars to observe, separated by commas (RETURN if not): ")
    standardsraw = standardinput.split(",")
    for standard in standardsraw:
        standards.append(standard.lstrip().rstrip())
    if standards == [""]: standards = []
elif len(args) == 6:
    directory = args[1]
    name = args[2]
    objs = args[3]
    filts = args[4]
    standards = args[5]
else:
    print("Invalid number of arguments. Either input folder name and parent directory, or folder name, parent directory, object list, filter list, and standard star list.")
    quit()

try:
    os.makedirs(os.path.join(directory,name,"Bias"))
    os.makedirs(os.path.join(directory,name,"Dark"),exist_ok=True)
    for filt in filts:
        os.makedirs(os.path.join(directory,name,"Flat",filt),exist_ok=True)
    for obj in objs:
        for filt in filts:
            os.makedirs(os.path.join(directory,name,"Science",obj,filt),exist_ok=True)
    for standard in standards:
        for filt in filts:
            os.makedirs(os.path.join(directory,name,"StandardStars",standard,filt),exist_ok=True)
except:
    print("Failed.")
else:
    print("Completed.")