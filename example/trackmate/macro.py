# @String input_dir
import os
import glob

from ij import IJ
from fiji.plugin.trackmate import Logger, Model, SelectionModel, Settings, TrackMate
from fiji.plugin.trackmate.detection import LabelImageDetectorFactory
from fiji.plugin.trackmate.features import FeatureFilter
from fiji.plugin.trackmate.gui.displaysettings import DisplaySettingsIO
from fiji.plugin.trackmate.io import CSVExporter, TmXmlWriter
from fiji.plugin.trackmate.tracking.jaqaman import SparseLAPTrackerFactory
from java.io import File


def run_trackmate(image_path):
    imp = IJ.openImage(image_path)
    if imp is None: return
    logger = Logger.IJ_LOGGER
    model = Model(); model.setLogger(Logger.IJ_LOGGER)
    settings = Settings(imp)
    settings.detectorFactory = LabelImageDetectorFactory()
    settings.detectorSettings["TARGET_CHANNEL"] = 1
    settings.detectorSettings["SIMPLIFY_CONTOURS"] = True
    settings.addSpotFilter(FeatureFilter('QUALITY', 0.0, True))
    settings.trackerFactory = SparseLAPTrackerFactory()
    settings.trackerSettings = settings.trackerFactory.getDefaultSettings()
    settings.trackerSettings.update({
        "LINKING_MAX_DISTANCE": 200.0,
        "GAP_CLOSING_MAX_DISTANCE": 200.0,
        "MAX_FRAME_GAP": 2,
        "ALLOW_TRACK_SPLITTING": True,
        "SPLITTING_MAX_DISTANCE": 15.0,
        "ALLOW_TRACK_MERGING": False,
        "MERGING_MAX_DISTANCE": 0.0
    })
    settings.addAllAnalyzers()
    trackmate = TrackMate(model, settings)
    if not trackmate.checkInput() or not trackmate.process(): imp.close(); return
    base_name = os.path.splitext(os.path.basename(image_path))[0]
    folder = os.path.dirname(image_path)
    ds = DisplaySettingsIO.readUserDefault()
    
    out_xml = File(os.path.join(folder, base_name+".xml"))
    writer = TmXmlWriter(out_xml)
    writer.appendSettings(settings)
    writer.appendModel(model)		
    writer.appendDisplaySettings(ds)	
    writer.writeToFile()
    CSVExporter.exportSpots(os.path.join(folder, base_name+"_spots.csv"), model, True)
    imp.close()


input_folder = "example/trackmate"
if os.path.exists(input_folder):
    for tiff_path in sorted(glob.glob(os.path.join(input_folder, "*merge.tif"))):
        name = os.path.splitext(os.path.basename(tiff_path))[0]
        print("NAME " + name)
        run_trackmate(tiff_path)
        print("SUCCESS")
else:
    print("Error: folder not found")
