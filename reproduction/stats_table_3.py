#!/usr/bin/python

import sys
import numpy as np
import re
import OpenEXR
import Imath

def loadEXRImage(inputFilePath):
    img = OpenEXR.InputFile(inputFilePath)
    dw = img.header()['dataWindow']

    size = (dw.max.x - dw.min.x + 1, dw.max.y - dw.min.y + 1)
    (rc, gc, bc) = img.channels("RGB", Imath.PixelType(Imath.PixelType.FLOAT))

    r = np.frombuffer(rc, dtype=np.float32)
    g = np.frombuffer(gc, dtype=np.float32)
    b = np.frombuffer(bc, dtype=np.float32)

    size = (dw.max.x - dw.min.x + 1, dw.max.y - dw.min.y + 1)

    return [r, g, b, size]

def computeRelMSE(resultImage, referenceImage):
    try:
        [resR, resG, resB, _] = loadEXRImage(resultImage);
        [refR, refG, refB, _] = loadEXRImage(referenceImage);
    except:
        print("failed to open the EXR images");
        exit(1);

    # "relMSE" (as used by our guiding papers, normalization applied to the absolute difference, error added rather than averaged over color channels)
    # as an intuition: normalize the absolute error (with epsilon) and square the result
    epsilon = 0.001

    r = (refR - resR) * (refR - resR) / (refR * (refR + (2 * epsilon)) + (epsilon * epsilon))
    g = (refG - resG) * (refG - resG) / (refG * (refG + (2 * epsilon)) + (epsilon * epsilon))
    b = (refB - resB) * (refB - resB) / (refB * (refB + (2 * epsilon)) + (epsilon * epsilon))

    perPixelError = (r + g + b)

    # discard the 0.1% highest values
    n = int(0.999*perPixelError.shape[0]);
    partitioned = np.partition(perPixelError, n-1);
    mean = np.mean(partitioned[0:n])

    return mean;

def parseMitsubaTime(timeString):
    if timeString.endswith('s'):
        return float(timeString.rstrip('s'))
    if timeString.endswith('m'):
        return float(timeString.rstrip('m'))*60.0
    if timeString.endswith('h'):
        return float(timeString.rstrip('h'))*60.0*60.0
    if timeString.endswith('d'):
        #the Mitsuba day has only 12 hours :) (I don't know why...)
        return float(timeString.rstrip('d'))*60.0*60.0*12.0
    else:
        return float(timeString);

def parseMitsubaNumber(numberString):
    if numberString.endswith('K'):
        return float(numberString.rstrip('K'))*1000.0
    if numberString.endswith('M'):
        return float(numberString.rstrip('M'))*1000000.0
    if numberString.endswith('G'):
        return float(numberString.rstrip('G'))*1000000000.0
    if numberString.endswith('T'):
        return float(numberString.rstrip('T'))*1000000000000.0
    else:
        return float(numberString);

def parseMitsubaMemory(memoryString):
    if memoryString.endswith(' B'):
        return float(memoryString.rstrip(' B'))
    if memoryString.endswith(' KiB'):
        return float(memoryString.rstrip(' KiB'))*1024.0
    if memoryString.endswith(' MiB'):
        return float(memoryString.rstrip(' MiB'))*1024.0*1024.0
    if memoryString.endswith(' GiB'):
        return float(memoryString.rstrip(' GiB'))*1024.0*1024.0*1024.0
    if memoryString.endswith(' TiB'):
        return float(memoryString.rstrip(' TiB'))*1024.0*1024.0*1024.0*1024.0
    if memoryString.endswith(' PiB'):
        return float(memoryString.rstrip(' PiB'))*1024.0*1024.0*1024.0*1024.0*1024.0
    else:
        return float(memoryString);

def computeStats(resultImage, resultLog, referenceImage):
    trainingTime = None;
    renderTime = None;
    trainingSamples = None;
    renderSamples = None;
    avgPathLength = None;
    avgNumComponents = None;
    numRegions = None;
    memoryGuiding = None;
    memorySamples = None;
    relMSE = computeRelMSE(resultImage, referenceImage);

    memoryGuidingTree = None

    logFile = open(resultLog, 'r')
    for line in logFile:
        trainingSamplesAndTimeRegEx = re.search(r"training phase finished after ([0-9]+) samples per pixel in ([0-9]+\.?[0-9]*[smhd]).", line);
        if trainingSamplesAndTimeRegEx:
            trainingSamples = int(trainingSamplesAndTimeRegEx.group(1));
            trainingTime = parseMitsubaTime(trainingSamplesAndTimeRegEx.group(2));
        renderSamplesAndTimeRegEx = re.search(r"rendered ([0-9]+) samples per pixel in ([0-9]+\.?[0-9]*[smhd]).", line);
        if renderSamplesAndTimeRegEx:
            renderSamples = int(renderSamplesAndTimeRegEx.group(1));
            renderTime = parseMitsubaTime(renderSamplesAndTimeRegEx.group(2));
        avgPathLengthRegEx = re.search(r"Average path length\s*:\s*([0-9]+\.?[0-9]*[Ee]?-?[0-9]*)", line);
        if avgPathLengthRegEx:
            avgPathLength = float(avgPathLengthRegEx.group(1));
        avgNumComponentsRegEx = re.search(r"avgNumComponents\s*=\s*([0-9]+\.?[0-9]*[Ee]?-?[0-9]*)", line);
        if avgNumComponentsRegEx:
            avgNumComponents = float(avgNumComponentsRegEx.group(1));
        numRegionsRegEx = re.search(r"numRegions\s*=\s*([0-9]+)", line);
        if numRegionsRegEx:
            numRegions = int(numRegionsRegEx.group(1))
        memoryGuidingTreeRegEx = re.search(r"usedNodeStorage\s*=\s*([0-9]+\.?[0-9]*[Ee]?-?[0-9]* [KMGTPiB]+)", line);
        if memoryGuidingTreeRegEx:
            memoryGuidingTree = parseMitsubaMemory(memoryGuidingTreeRegEx.group(1))/(1024.0*1024.0);
        memoryGuidingRegEx = re.search(r"usedRegionStorage\s*=\s*([0-9]+\.?[0-9]*[Ee]?-?[0-9]* [KMGTPiB]+)", line);
        if memoryGuidingRegEx:
            memoryGuiding = parseMitsubaMemory(memoryGuidingRegEx.group(1))/(1024.0*1024.0)+memoryGuidingTree;
        memorySamplesRegEx = re.search(r"usedTrainingSampleStorage\s*=\s*([0-9]+\.?[0-9]*[Ee]?-?[0-9]* [KMGTPiB]+)", line);
        if memorySamplesRegEx:
            memorySamples = parseMitsubaMemory(memorySamplesRegEx.group(1))/(1024.0*1024.0);

    return str(trainingTime)+';'+str(renderTime)+';'+str(trainingSamples)+';'+str(renderSamples)+';'+str(avgPathLength)+';'+str(avgNumComponents)+';'+str(numRegions)+';'+str(memoryGuiding)+';'+str(memorySamples)+';'+str(relMSE);

if __name__ == '__main__':
    if (len(sys.argv) >= 2 and sys.argv[1] == '--header'):
        print("trainingTime(s);renderTime(s);trainingSamples;renderSamples;avgPathLength;avgNumComponents;numRegions;memoryGuiding(MB);memorySamples(MB);relMSE");
        exit(0);

    if len(sys.argv) < 3:
        print("usage: "+sys.argv[0]+" <result.exr> <result.log> <reference.exr>");
        exit(1);

    print(computeStats(sys.argv[1], sys.argv[2], sys.argv[3]));
