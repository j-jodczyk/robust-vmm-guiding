import sys, os

mtsdir = '/your/path/to/mitsuba/dist' # TODO: put your Mitsuba path here
sys.path.append(mtsdir + '/python/'+str(sys.version_info.major)+'.'+str(sys.version_info.minor))
os.environ['PATH'] = mtsdir + os.pathsep + os.environ['PATH'];

import mitsuba
import mitsuba.core as mtsCore
import mitsuba.render as mtsRender
# Pathguiding Module
import mitsuba.guiding as mtsGuiding

import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

import multiprocessing

scheduler = mtsCore.Scheduler.getInstance()

# Start up the scheduling system with one worker per local core
for i in range(0, multiprocessing.cpu_count()):
    scheduler.registerWorker(mtsCore.LocalWorker(i, 'wrk%i' % i))
scheduler.start()

fileResolver = mtsCore.Thread.getThread().getFileResolver();
fileResolver.appendPath("/your/path/to/scenes/torus/"); # TODO: put your scene path here

paramMap = mtsCore.StringMap();
scene = mtsRender.SceneHandler.loadScene(fileResolver.resolve("torus_auto.xml"), paramMap);

scene.configure();
scene.initialize();

# TODO: choose whatever pixel you like
pixel = mtsCore.Point2(640,435);
its = mtsGuiding.PathGuidingUtilities.getFirstSmoothSurfaceInteraction(scene, pixel);

print(its.p)
print(its.shFrame.n)

# if you like, you can render the scene using a spherical camera with a different integrator here
#sphericalView = mtsGuiding.PathGuidingUtilities.renderSphericalView(scene, its.p, "path", 32);
#sphericalView.plot();
#sphericalView.write('torus_640_435.exr');


# TODO: choose what file you want to look at
guidingFieldStream = mtsCore.FileStream('training/32spp_total_guiding_field.serialized', mtsCore.FileStream.EReadOnly);
guidingField = mtsGuiding.GuidingField(guidingFieldStream);

# choose the size of the rendered image of the distribution here
mtsGuiding.PathGuidingUtilities.renderSize = mtsCore.Vector2i(512, 256);

region = guidingField.guidingTree.getRegion(its.p)
# optionally load the matching samples (make sure to export them beforehand - its not enabled by default)
#samples = mtsGuiding.PathGuidingUtilities.loadSampleRange('training/32spp_total_samples.serialized', region.dStart, region.nSamples);

vmm = region.distribution;
vmmPDF = mtsGuiding.PathGuidingUtilities.renderVMMPDF(vmm);

def getGridAnd2DArrayFromBitmap(bitmap):
    gridRes = bitmap.getSize();
    gridX = np.linspace(math.pi, -math.pi, gridRes.x);
    gridY = np.linspace(0, math.pi, gridRes.y);

    greyscaleBitmap = bitmap.convert(mtsCore.Bitmap.ELuminance, mtsCore.Bitmap.EFloat);
    data = np.array(greyscaleBitmap.buffer());

    return [gridX, gridY, data];

class FalsecolorNormalize(colors.Normalize):
    def __init__(self, vmin=None, vmax=None, clip=False):
        colors.Normalize.__init__(self, vmin, vmax, clip);

    def __call__(self, value, clip=None):
        with np.errstate(divide='ignore'):
            return np.where(value > 0.0, np.log2(value), -float("inf"))/10.0+0.5;

def plotBitmapInFalsecolor(bitmap, alpha=1.0):
    if len(plt.get_fignums()) == 0:
        plt.figure(dpi=100, figsize=[20.0, 6.0]);
    plt.xlabel(r'$\varphi$');
    plt.ylabel(r'$\theta$');
    #this is the way Mitsuba expects envmaps to be laid out
    plt.xlim(math.pi, -math.pi);
    plt.ylim(math.pi, 0.0);

    [gridX, gridY, data] = getGridAnd2DArrayFromBitmap(bitmap);

    plt.imshow(data, extent=[math.pi, -math.pi, math.pi, 0], alpha=alpha, cmap="viridis", norm=FalsecolorNormalize(), resample=False, interpolation='none');

plotBitmapInFalsecolor(vmmPDF);
plt.show();
# alternatively you can also write this bitmap to an exr file
vmmPDF.write('vmm_pdf.exr');
