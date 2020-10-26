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

def computeStats(resultImage, resultLog, referenceImage):
    renderSamples = None;
    relMSE = computeRelMSE(resultImage, referenceImage);

    memoryGuidingTree = None

    logFile = open(resultLog, 'r')
    for line in logFile:
        renderSamplesAndTimeRegEx = re.search(r"rendered ([0-9]+) samples per pixel in ([0-9]+\.?[0-9]*[smhd]).", line);
        if renderSamplesAndTimeRegEx:
            renderSamples = int(renderSamplesAndTimeRegEx.group(1));

    return str(renderSamples)+';'+str(relMSE);

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("usage: "+sys.argv[0]+" <result.exr> <result.log> [<result2.exr> <result2.log> ...] <reference.exr>");
        exit(1);

    result = '';

    i = 1;
    while (i+2 < len(sys.argv)):
        result += computeStats(sys.argv[i], sys.argv[i+1], sys.argv[len(sys.argv)-1])+';'
        i += 2;

    print(result.rstrip(';'));
