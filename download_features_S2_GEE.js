var table = Table users/dgaso/Field_001

var aoi = table001.geometry();
var Field = 'Field001'
var proj = 'EPSG:32721' 

var start_date = '2019-11-15';
var end_date = '2020-04-30';

// Name of the export folder [string]
var exp_fol = 'S2_Features';

// Spatial resolution of the NDMI stack [numeric]
var sp_res = 20;
var tiles = ['15TTF'];

// sequence starts beggining of season in DOY and end -15 of end of season:

// Setting for US:
//growing season in US: Apr30 to Oct15
//var biweeks = ee.List.sequence(120,270,15); 
//var start_date = '2019-04-30';
//var end_date = '2019-10-15';

/////  Setting for UY:
//growing season in UY single crop: Nov15 to Apr30:
var biweeks = ee.List.sequence(320,350,15); 
var biweeks2 = ee.List.sequence(0,105,15); 

//growing season in UY doble crop: Dec15 to May30:
//var biweeks = ee.List.sequence(350,350,15); 
//var biweeks2 = ee.List.sequence(0,105,15); //135
//var start_date = '2018-12-15';
//var end_date = '2019-04-30';

print(biweeks,biweeks2)


// Vegetation Index list to download: 'CIG_mean', 'WDRVI_mean', 'TCARI_mean','MCARI_mean', 'NDVI_mean','GNDVI_mean','NDMI_mean'
// Spectral Bands list to download: 'B3_mean', 'B4_mean','B5_mean','B6_mean','B7_mean','B8A_mean', 'B11_mean', 'B12_mean'
// Select 1 index and 1 band
var band = 'B2_mean'
var index = 'CIRE_mean'

// Params could mask
var CLOUD_FILTER = 100
var CLD_PRB_THRESH = 10
var NIR_DRK_THRESH = 0.15
var CLD_PRJ_DIST = 1
var BUFFER = 50
 
// SPECTRAL MODULE: connect to awesone spectral
var spectral = require("users/dmlmont/spectral:spectral");

// ---------------- DO NOT WRITE AFTER THIS LINE ---------------------------
// function to merge bands in multi-image
var mergeBands = function mergeBands(image, previous) {
  return ee.Image(previous).addBands(image);
};

/////////////////////  Compute vegetation indexes by using awesome spectral ////////////////

// DATASET TO USE: SENTINEL-2 SR
var dataset = 'COPERNICUS/S2_SR';

// FUNCTION TO MAP OVER AN IMAGE COLLECTION
function addIndices(img) {
  
// REQUIRED PARAMETERS ACCORDING TO THE REQUIRED BANDS:

  var parameters = {
    "A": img.select("B1"),
    "B": img.select("B2"),
    "G": img.select("B3"),
    "R": img.select("B4"),
    "RE1": img.select("B5"),
    "RE2": img.select("B6"),
    "RE3": img.select("B7"),
    "RE4": img.select("B8A"),
    "N": img.select("B8"),
    "S1": img.select("B11"),
    "S2": img.select("B12"),
    "alpha": 0.5,
  };
  
  // SCALE THE IMAGE
  img = spectral.scale(img,dataset);
  
// COMPUTE INDEXES: list:https://awesome-ee-spectral-indices.readthedocs.io/en/latest/list.html#kernel
  return spectral.computeIndex(img,["CIRE","CIG","WDRVI","TCARI","MCARI","NDVI","GNDVI","NDMI"],parameters);
}

print("MCARI",spectral.indices.MCARI.formula);
print("TCARI",spectral.indices.TCARI.formula);
print("CIG",spectral.indices.CIG.formula);
print("WDRVI",spectral.indices.WDRVI.formula);

//////////////////////////////  Cloud mask  //////////////////////////////////////////////////////////////////

function get_s2_sr_cld_col(aoi, start_date, end_date) {
    // # Import and filter S2 SR.
    var s2_sr_col = (ee.ImageCollection('COPERNICUS/S2_SR')
        //.filter(ee.Filter.inList('MGRS_TILE', tiles))
        .filterBounds(aoi)
        .filterDate(start_date, end_date)
        .filter(ee.Filter.lte('CLOUDY_PIXEL_PERCENTAGE', CLOUD_FILTER)))

    // # Import and filter s2cloudless.
    var s2_cloudless_col = (ee.ImageCollection('COPERNICUS/S2_CLOUD_PROBABILITY')
        //.filter(ee.Filter.inList('MGRS_TILE', tiles))
        .filterBounds(aoi)
        .filterDate(start_date, end_date))

    // # Join the filtered s2cloudless collection to the SR collection by the 'system:index' property.
    return ee.ImageCollection(ee.Join.saveFirst('s2cloudless').apply({
        'primary': s2_sr_col,
        'secondary': s2_cloudless_col,
        'condition': ee.Filter.equals({
            'leftField': 'system:index',
            'rightField': 'system:index'
        })
    }))
}

function add_cloud_bands(img) {
    // # Get s2cloudless image, subset the probability band.
    var cld_prb = ee.Image(img.get('s2cloudless')).select('probability')

    // # Condition s2cloudless by the probability threshold value.
    var is_cloud = cld_prb.gt(CLD_PRB_THRESH).rename('clouds')

    // # Add the cloud probability layer and cloud mask as image bands.
    return img.addBands(ee.Image([cld_prb, is_cloud]))
    
}


function add_shadow_bands(img) {
    // # Identify water pixels from the SCL band.
    var not_water = img.select('SCL').neq(6)

    // # Identify dark NIR pixels that are not water (potential cloud shadow pixels).
    var SR_BAND_SCALE = 1e4
    var dark_pixels = img.select('B8').lt(NIR_DRK_THRESH*SR_BAND_SCALE).multiply(not_water).rename('dark_pixels')

    // # Determine the direction to project cloud shadow from clouds (assumes UTM projection).
    var shadow_azimuth = ee.Number(90).subtract(ee.Number(img.get('MEAN_SOLAR_AZIMUTH_ANGLE')));

    // # Project shadows from clouds for the distance specified by the CLD_PRJ_DIST input.
    var cld_proj = (img.select('clouds').directionalDistanceTransform(shadow_azimuth, CLD_PRJ_DIST*10)
        .reproject({'crs': img.select(0).projection(), 'scale': 100})
        .select('distance')
        .mask()
        .rename('cloud_transform'))

    // # Identify the intersection of dark pixels with cloud shadow projection.
    var shadows = cld_proj.multiply(dark_pixels).rename('shadows')

    // # Add dark pixels, cloud projection, and identified shadows as image bands.
    return img.addBands(ee.Image([dark_pixels, cld_proj, shadows]))
}

function add_cld_shdw_mask(img) {
    // # Add cloud component bands.
    var img_cloud = add_cloud_bands(img)

    // # Add cloud shadow component bands.
    var img_cloud_shadow = add_shadow_bands(img_cloud)

    // # Combine cloud and shadow mask, set cloud and shadow as value 1, else 0.
    var is_cld_shdw = img_cloud_shadow.select('clouds').add(img_cloud_shadow.select('shadows')).gt(0)

    // # Remove small cloud-shadow patches and dilate remaining pixels by BUFFER input.
    // # 20 m scale is for speed, and assumes clouds don't require 10 m precision.
    is_cld_shdw = (is_cld_shdw.focal_min(2).focal_max(BUFFER*2/20)
        .reproject({'crs': img.select([0]).projection(), 'scale': 20})
        .rename('cloudmask'))
    // # Add the final cloud-shadow mask to the image.
    return img_cloud_shadow.addBands(is_cld_shdw)
}


function apply_cld_shdw_mask(img) {
    // # Subset the cloudmask band and invert it so clouds/shadow are 0, else 1.
    var not_cld_shdw = img.select('cloudmask').not()

    // # Subset reflectance bands and update their masks, return the result.
    return img.select('B.*').updateMask(not_cld_shdw)
}

var s2_sr_cld_col = get_s2_sr_cld_col(aoi, start_date, end_date)

print('s2_sr_cld_col',s2_sr_cld_col)

var s2_sr_median = (s2_sr_cld_col.map(add_cld_shdw_mask)
                             .map(apply_cld_shdw_mask)
                             )
print('s2_sr_median',s2_sr_median)                             

// convert to float
var S2 = s2_sr_median.map(function(image) {return image.toFloat()});
var withVI = S2.map(addIndices);

////////////////////////     Biweekly means    //////////////////////////////////////////////

var Biweekly_mean = ee.ImageCollection.fromImages(
biweeks.map(function(m) {
m = ee.Number(m)
return withVI.filter(ee.Filter.calendarRange(m, m.add(15), 'day_of_year'))
.filterBounds(aoi)
.reduce(ee.Reducer.mean())
.set('day_of_year', m);
}));
//Map.addLayer(Monthly_mean.first().select('CIG_mean').clip(roi), {min: 0, max: 3}, 'S2 monthly average')
print('Biweekly_mean',Biweekly_mean);


var Biweekly_mean2 = ee.ImageCollection.fromImages(
biweeks2.map(function(m) {
m = ee.Number(m)
return withVI.filter(ee.Filter.calendarRange(m, m.add(15), 'day_of_year'))
.filterBounds(aoi)
.reduce(ee.Reducer.mean())
.set('day_of_year', m);
}));
print('Biweekly_mean2',Biweekly_mean2);

var Biweekly_mean = Biweekly_mean.merge(Biweekly_mean2).sort('system:index')
print('Biweekly_mean FIANL',Biweekly_mean)

var n = Biweekly_mean.size();
var colList = Biweekly_mean.toList(n);
var img = ee.Image(colList.get(5));
Map.addLayer(img.clip(aoi), {min: 0, max: 0.5}, 'sub image')

//////////////////////////  clip collection and download    /////////////////////////// 

//get crs and transform to download
var n = withVI.size();
var withVIList = withVI.toList(n);
var pr = ee.Image(withVIList.get(1)).select('B5').projection().crs();
var transform = ee.Image(withVIList.get(1)).select('B5').projection().transform().getInfo();
print('transf',transform)

//clip collection by AOI
var clip_aoi = Biweekly_mean
    .map(function(image){return image.clip(aoi)
    });
print('clip_aoi',clip_aoi);

// select index and stack images to download (per index or spectral band)
var download_index = clip_aoi
  .map(function(image) {return image.select(index)
  })
.iterate(mergeBands, ee.Image([]));
print('download_aoi',download_aoi);

Export.image.toDrive({
              image: download_index,
              description: index + '_' + Field,
              scale: sp_res,
              crs: proj,
              //crsTransform:[20,0,399960,0,-20,6300040], 
              fileFormat: 'GeoTIFF',
              folder: exp_fol,
              region: aoi,
              maxPixels: 10e10});

//////////////////   download bands from S2  ///////////////////////
var download_band = clip_aoi
  .map(function(image) {return image.select(band)
  })
  .iterate(mergeBands, ee.Image([]));
 
// export to drive
Export.image.toDrive({
              image: download_band,
              description: B2 + '_' + Field,
              scale: sp_res,
              crs: proj,
              //crsTransform:[20,0,399960,0,-20,6300040],
              fileFormat: 'GeoTIFF',
              folder: exp_fol,
              region: aoi,
              maxPixels: 10e10});

  .iterate(mergeBands, ee.Image([]));
  

Map.addLayer(aoi,null,'Region of interest');
Map.centerObject(aoi);
print('withVI',withVI)

//Map.centerObject(testPoint, 10)
var chart = ui.Chart.image.series({
    imageCollection: withVI.select('NDVI'),
    reducer:ee.Reducer.mean(),
    scale:20,
    region: aoi
    }).setOptions({
      interpolateNulls: true,
      lineWidth: 1,
      pointSize: 3,
      title: ' Seasonality',
      vAxis: {title: 'NDVI'},
      hAxis: {title: 'Date', format: 'YYYY-MMM', gridlines: {count: 12}}
    })
print(chart)
