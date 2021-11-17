import arcpy
import os
import sys
# Dongdong Kong, CUG-atmos, 2021-03-22

# from arcpy.sa import *
arcpy.CheckOutExtension("spatial") # necessary!
arcpy.env.workspace = "."          # necessary!
arcpy.env.overwriteOutput = True
arcpy.env.parallelProcessingFactor = "40%"


def filldem(dem, outfile = "demfill.tif"):
    print("fill dem ...")
    outFill = arcpy.sa.Fill(dem)
    outFill.save(outfile)

def cal_flowdir(dem = "demfill.tif", outfile = "flowdir.tif"):
    print("running flowdir ...")
    print(dem)
    flowdir = arcpy.sa.FlowDirection(dem, "FORCE")
    flowdir.save(outfile)

def cal_flowaccu(flowdir = "flowdir.tif", outfile = "flowaccu.tif"):
    print("running flowaccu ...")
    outFlowAccumulation = arcpy.sa.FlowAccumulation(flowdir)
    outFlowAccumulation.save(outfile)

def cal_watershed(flowdir, pour, outfile = "watershed"):
    outWatershed = arcpy.sa.Watershed(flowdir, pour)
    outWatershed.save(outfile)

def stream_order(accu):
    ""

if __name__ == '__main__':
    # dem = "N:/ChinaWater/project/dem_guanshan.tif"
    nargs = len(sys.argv) 
    if nargs >= 2:
        # print 'Number of arguments:', len(sys.argv), 'arguments.'
        # print 'Argument List:', str(sys.argv)
        # example: 
        #	fastbasin dem.tif
        # 	fastbasin demfill.tif flowdir
        # 	fastbasin flowdir.tif flowaccu
        infile = sys.argv[1]
        arcpy.env.workspace = os.path.dirname(os.path.abspath(infile))
                
        if (nargs > 2):
            begin = sys.argv[2] # fill, flowdir, flowaccu

            options = {"fill": 1, "flowdir": 2, "flowaccu": 3}
            Ibegin = options[begin]

            if Ibegin <= 1: 
                filldem(infile)
                cal_flowdir()
                cal_flowaccu()	
            if Ibegin <= 2: 
                cal_flowdir(infile)
                cal_flowaccu()
            if Ibegin <= 3: 
                cal_flowaccu(infile)
        else:
            filldem(infile)
            cal_flowdir()
            cal_flowaccu()
