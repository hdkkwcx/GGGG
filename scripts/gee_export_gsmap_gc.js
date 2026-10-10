// Run in the GEE Code Editor, once per YEAR (2022, 2023, 2024).
// Exports hourly gauge-calibrated GSMaP (hourlyPrecipRateGC) at the stations, one CSV per month.
//
// Files keep the names and value column of the uncalibrated export
// (gsmap_hubei_hourly_long_YYYYMM.csv, gsmap_mm_h) so prepare_satellite_wide reads them unchanged;
// only the Drive folder differs (GEE_GSMAP_GC_HUBEI). Download it to data/GEE_GSMAP_GC_HUBEI.
//
// For 2022-2024 every image should be 'permanent': the provisional (NRT) gauge adjustment is replaced
// by the one fitted to that day's NOAA/CPC daily gauge analysis. The print below checks this, and
// gauge_quality (daily mean number of gauges used in the pixel) is exported alongside.

var YEAR = 2022;
var MONTHS = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12];
var START = ee.Date.fromYMD(YEAR, 1, 1);
var END = START.advance(1, 'year'); // half-open

var AOI = ee.Geometry.Rectangle(
  [109.4, 31.2, 111.6, 33.4], null, false);

var stationsRaw = ee.FeatureCollection(
  'projects/booming-list-489917-p4/assets/Hubei-STATIONS');

// Tolerate station_id/site/id naming differences in the station asset.
function pickKeyContains(feature, keyword) {
  var keys = ee.List(feature.propertyNames());
  var matched = keys.map(function(key) {
    key = ee.String(key);
    return ee.Algorithms.If(
      key.replace('_', '', 'g').toLowerCase().index(keyword).gte(0), key, null);
  }).removeAll([null]);
  return ee.Algorithms.If(
    matched.size().gt(0), feature.get(ee.String(matched.get(0))), null);
}

function normalizeStation(feature) {
  var hiddenSite = pickKeyContains(feature, 'site');
  var stationId = ee.Algorithms.If(
    feature.get('station_id'), feature.get('station_id'),
    ee.Algorithms.If(
      feature.get('site'), feature.get('site'),
      ee.Algorithms.If(
        hiddenSite, hiddenSite,
        ee.Algorithms.If(
          feature.get('id'), feature.get('id'), feature.get('system:index')))));
  var lon = ee.Algorithms.If(
    feature.get('lon'), feature.get('lon'), pickKeyContains(feature, 'lon'));
  var lat = ee.Algorithms.If(
    feature.get('lat'), feature.get('lat'), pickKeyContains(feature, 'lat'));

  return ee.Feature(
    ee.Geometry.Point([ee.Number(lon), ee.Number(lat)]), {
      station_id: ee.String(stationId),
      lon: ee.Number(lon),
      lat: ee.Number(lat)
    });
}

// Numeric closed-interval filter, avoiding geometry-boundary float error.
var stations = stationsRaw
  .map(normalizeStation)
  .filter(ee.Filter.notNull(['station_id', 'lon', 'lat']))
  .filter(ee.Filter.gte('lon', 109.4))
  .filter(ee.Filter.lte('lon', 111.6))
  .filter(ee.Filter.gte('lat', 31.2))
  .filter(ee.Filter.lte('lat', 33.4));

// GSMaP v8, gauge-calibrated hourly rate in mm/hr.
var gsmapAll = ee.ImageCollection('JAXA/GPM_L3/GSMaP/v8/operational')
  .filterDate(START, END);
var gsmap = gsmapAll.select(['hourlyPrecipRateGC', 'gaugeQualityInfo']);

print('Images in ' + YEAR + ':', gsmapAll.size());
print('Provisional images (should be 0):',
  gsmapAll.filter(ee.Filter.neq('status', 'permanent')).size());

function hourlyLongTable(monthStart, monthEnd) {
  var hourCount = monthEnd.difference(monthStart, 'hour');
  var hours = ee.List.sequence(0, hourCount.subtract(1));

  return ee.FeatureCollection(hours.map(function(offset) {
    var hourStart = monthStart.advance(ee.Number(offset), 'hour');
    var hourEnd = hourStart.advance(1, 'hour');
    // One GSMaP image per hour; mean keeps the mm/hr value.
    var image = gsmap.filterDate(hourStart, hourEnd).mean();
    var sampled = image.reduceRegions({
      collection: stations,
      reducer: ee.Reducer.mean(),
      scale: 11132,
      tileScale: 4
    });

    return sampled.map(function(feature) {
      return ee.Feature(feature.geometry(), {
        gsmap_mm_h: feature.get('hourlyPrecipRateGC'),
        gauge_quality: feature.get('gaugeQualityInfo'),
        station_id: feature.get('station_id'),
        time: hourStart.format("yyyy-MM-dd'T'HH:mm:ss'Z'")
      });
    }).filter(ee.Filter.notNull(['station_id', 'gsmap_mm_h']));
  })).flatten();
}

function exportMonth(year, month) {
  var monthStart = ee.Date.fromYMD(year, month, 1);
  var monthEnd = monthStart.advance(1, 'month');
  var table = hourlyLongTable(monthStart, monthEnd);
  var suffix = String(year) + (month < 10 ? '0' + month : String(month));

  print('Created export task: gsmap_gc_hubei_hourly_long_' + suffix);
  Export.table.toDrive({
    collection: table,
    description: 'gsmap_gc_hubei_hourly_long_' + suffix,
    folder: 'GEE_GSMAP_GC_HUBEI',
    fileNamePrefix: 'gsmap_hubei_hourly_long_' + suffix,
    fileFormat: 'CSV',
    selectors: ['system:index', 'gsmap_mm_h', 'gauge_quality', 'station_id', 'time', '.geo']
  });
}

MONTHS.forEach(function(month) {
  exportMonth(YEAR, month);
});

Map.centerObject(AOI, 7);
Map.addLayer(stations, {color: 'yellow'}, 'Hubei stations');
