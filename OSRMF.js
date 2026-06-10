/**
 * @title Optical and SAR Data Fusion for Paddy Rice Mapping
 * @description Extracts paddy rice extent using Sentinel-1 (VH) and Sentinel-2 (NDVI, MNDWI) 
 * based on phenological windows (transplanting and peak growth phases).
 */

// ==============================================================================
// 1. CONFIGURATION & PARAMETERS
// ==============================================================================
var CONFIG = {
  region: {
    country: 'China',
    province: 'Jilin Sheng' // Note: FAO/GAUL uses Wade-Giles spelling for Jilin
  },
  dates: {
    start: '2019-01-01',
    end: '2019-12-31'
  },
  // Day of Year (DOY) Windows for Rice Phenology
  doy: {
    ndviTransplanting: [110, 180],
    ndviPeakGrowth: [200, 280],
    mndwiFlooding: [145, 180],
    sarMinEarly: [120, 190],
    sarMaxLate: [200, 270]
  },
  // SAR Backscatter Empirical Thresholds
  sarThresholds: {
    v: ee.Image.constant(-15.56163934),
    w: ee.Image.constant(-22.81511617)
  },
  vis: {
    min: 0, max: 1, scale: 10, palette: ["black", "green"]
  }
};

// ==============================================================================
// 2. REGION OF INTEREST (ROI)
// ==============================================================================
var roi = ee.FeatureCollection("FAO/GAUL/2015/level1")
  .filter(ee.Filter.eq('ADM0_NAME', CONFIG.region.country))
  .filter(ee.Filter.eq('ADM1_NAME', CONFIG.region.province));

Map.centerObject(roi, 7);
Map.addLayer(roi.style({color: 'red', fillColor: '00000000', width: 2}), {}, 'Study Area Boundary');

// ==============================================================================
// 3. FUNCTIONS
// ==============================================================================

/**
 * Applies a custom Lee-Sigma filter to smooth SAR speckle noise
 */
function applyLeeSigmaFilter(image) {
  var size = 7;
  var kernel = ee.Kernel.square(size);
  
  var count = image.reduceNeighborhood({reducer: ee.Reducer.count(), kernel: kernel});
  var mean = image.reduceNeighborhood({reducer: ee.Reducer.mean(), kernel: kernel});
  var variance = image.reduceNeighborhood({reducer: ee.Reducer.variance(), kernel: kernel});
  
  var snr = variance.divide(mean.multiply(mean));
  var filtered = mean.updateMask(snr.lt(1.0)).convolve(ee.Kernel.gaussian(2.0));
  
  return filtered.copyProperties(image, ['system:time_start', 'system:time_end']);
}

/**
 * Calculates and adds NDVI and MNDWI bands
 */
function addSpectralIndices(image) {
  var ndvi = image.normalizedDifference(['B8', 'B4']).rename('NDVI');
  var mndwi = image.normalizedDifference(['B3', 'B11']).rename('MNDWI');
  return image.addBands([ndvi, mndwi]).float();
}

// ==============================================================================
// 4. DATA PREPARATION (Sentinel-1 & Sentinel-2)
// ==============================================================================

// --- Sentinel-1 SAR Processing ---
var s1Filtered = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterBounds(roi)
  .filterDate(CONFIG.dates.start, CONFIG.dates.end)
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
  .filter(ee.Filter.eq('instrumentMode', 'IW'))
  .map(function(image) {
    var validMask = image.mask().and(image.lt(-30.0).not());
    return image.updateMask(validMask);
  })
  .select('VH')
  .map(applyLeeSigmaFilter); 

// --- Sentinel-2 Optical Processing ---
var csPlus = ee.ImageCollection('GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED');
var s2Clouds = ee.ImageCollection('COPERNICUS/S2_CLOUD_PROBABILITY');

var s2Processed = ee.ImageCollection("COPERNICUS/S2_HARMONIZED")
  .filterBounds(roi)
  .filterDate(CONFIG.dates.start, CONFIG.dates.end)
  .linkCollection(csPlus, ['cs_cdf'])
  .linkCollection(s2Clouds, ['probability'])
  .map(addSpectralIndices)
  .map(function(img) {
    // Cloud Score+ Masking
    return img.updateMask(img.select('cs_cdf').gte(0.5));
  });

// ==============================================================================
// 5. PHENOLOGICAL METRICS EXTRACTION
// ==============================================================================

var ndviMinTrans = s2Processed.filter(ee.Filter.calendarRange(CONFIG.doy.ndviTransplanting[0], CONFIG.doy.ndviTransplanting[1], 'day_of_year')).select('NDVI').min();
var ndviMaxPeak  = s2Processed.filter(ee.Filter.calendarRange(CONFIG.doy.ndviPeakGrowth[0], CONFIG.doy.ndviPeakGrowth[1], 'day_of_year')).select('NDVI').max();
var mndwiMeanFlood = s2Processed.filter(ee.Filter.calendarRange(CONFIG.doy.mndwiFlooding[0], CONFIG.doy.mndwiFlooding[1], 'day_of_year')).select('MNDWI').mean();

// Note: applyLeeSigmaFilter renames the output band to 'VH_mean'
var vhMinEarly = s1Filtered.filter(ee.Filter.calendarRange(CONFIG.doy.sarMinEarly[0], CONFIG.doy.sarMinEarly[1], 'day_of_year')).select('VH_mean').min();
var vhMaxLate  = s1Filtered.filter(ee.Filter.calendarRange(CONFIG.doy.sarMaxLate[0], CONFIG.doy.sarMaxLate[1], 'day_of_year')).select('VH_mean').max();

// ==============================================================================
// 6. RICE MAPPING LOGIC & DECISION RULES
// ==============================================================================

var ONE = ee.Image.constant(1);
var V = CONFIG.sarThresholds.v;
var W = CONFIG.sarThresholds.w;

// Optical Rules
var f1 = ONE.subtract(ndviMinTrans.multiply(ndviMinTrans));
var f2 = ONE.subtract(ONE.subtract(ndviMaxPeak).pow(2));
    f2 = f2.updateMask(f2.gt(0)).unmask(0);
var f3 = mndwiMeanFlood.gte(0);

// SAR Phenomenological Rules
var exponent = V.subtract(W).divide(2).subtract(vhMaxLate.subtract(vhMinEarly));
var fH = ONE.divide(ONE.add(exponent.exp()));

var W_img = ee.Image(0)
  .where(vhMinEarly.gte(V), 1)
  .where(vhMinEarly.gte(W).and(vhMinEarly.lt(V)), vhMinEarly.subtract(W).divide(V.subtract(W)));
var fW = ONE.subtract(W_img.multiply(W_img));

var V_img = ee.Image(0)
  .where(vhMaxLate.lte(W), 1)
  .where(vhMaxLate.gt(W).and(vhMaxLate.lte(V)), V.subtract(vhMaxLate).divide(V.subtract(W)));
var fV = ONE.subtract(V_img.multiply(V_img));

// Synthesis & Thresholding
var fw = fW.gt(0.6).unmask(0);
var fh = fH.gt(0.6).unmask(0);
var fv = fV.gt(0.6).unmask(0);

// Calculate Final Index
var OSDT = f1.multiply(f2).multiply(f3).multiply(fw).multiply(fh).multiply(fv);
var riceMask = OSDT.gte(0.7);

// Apply final clipping only to the output result
var riceFinal = riceMask.updateMask(riceMask.neq(0)).clip(roi);

Map.addLayer(riceFinal, CONFIG.vis, "Paddy Rice Map");
