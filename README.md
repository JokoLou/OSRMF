# Optical–SAR Rice Mapping Framework (OSRMF)

OSRMF is a phenology-based paddy rice mapping framework that integrates Sentinel-1 SAR and Sentinel-2 optical observations within Google Earth Engine (GEE). The framework combines vegetation, water, and radar backscatter dynamics during key rice growth stages to identify paddy rice fields under varying environmental conditions.

## Features

* Integration of Sentinel-1 SAR and Sentinel-2 optical imagery
* Cloud-robust mapping through adaptive optical–SAR fusion
* Phenology-based rice identification
* Automatic extraction of flooding/transplanting and peak-growth signals
* Implementation in Google Earth Engine (JavaScript)

## Methodology

The framework utilizes three categories of indicators:

### Optical Indicators

* NDVI minimum during the transplanting stage
* NDVI maximum during the peak-growth stage
* Mean MNDWI during the flooding stage

### SAR Indicators

* Minimum VH backscatter during early growth
* Maximum VH backscatter during late growth
* Logistic transformation of temporal VH dynamics

### Fusion Strategy

The final Optical–SAR Decision Index (OSDT) is calculated by combining optical and SAR indicators:

```text
OSDT = f1 × f2 × f3 × fW × fH × fV
```

Pixels with OSDT ≥ 0.8 are classified as paddy rice.

## Data Sources

### Sentinel-1 SAR

* Dataset: COPERNICUS/S1_GRD
* Polarization: VH
* Mode: IW

### Sentinel-2 Optical

* Dataset: COPERNICUS/S2_HARMONIZED
* Indices:

  * NDVI
  * MNDWI

### Cloud Mask

* GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED
* COPERNICUS/S2_CLOUD_PROBABILITY

### Administrative Boundaries

* FAO/GAUL/2015/level1

## Workflow

1. Define study area.
2. Acquire Sentinel-1 and Sentinel-2 imagery.
3. Apply SAR speckle filtering.
4. Generate NDVI and MNDWI time series.
5. Extract phenological metrics.
6. Calculate optical indicators.
7. Calculate SAR indicators.
8. Generate OSDT index.
9. Extract paddy rice map.

## Example

The example script maps paddy rice in Jilin Province, China, for the 2019 growing season.

```javascript
region: China - Jilin Sheng
year: 2019
spatial resolution: 10 m
```

## Requirements

* Google Earth Engine account
* GEE JavaScript Code Editor

## Citation

If you use this code in scientific research, please cite:

Lou, J., et al. Optical–SAR Rice Mapping Framework (OSRMF): A phenology-based approach for large-scale paddy rice mapping using Sentinel-1 and Sentinel-2 observations.

## Author

Joko Lou

Ningbo University / Zhejiang University

Research Area: Agricultural Remote Sensing and Paddy Rice Mapping
