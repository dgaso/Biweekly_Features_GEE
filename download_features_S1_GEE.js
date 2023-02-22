// This Script generates and downloads biweekly averages of vegetation indexes 
// from Sentinel-1, based on  https://doi.org/10.1016/j.isprsjprs.2021.05.013

var table = Table users/dgaso/Field_001
var geometry1 = table.geometry();
var Field = 'Field001'
var proj = 'EPSG:32721' 

var Start_date = '2019-11-15';
var Finish_date = '2020-04-30';

// Name of the export folder [string]
var exp_fold = 'S1_Features';
var tiles = ['15TTF']; //tile of S2 to match pixels with S1

// sequence starts beggining of season in DOY and end -15 of end of season:

/////  Setting for US:
//growing season in US: Apr30 to Oct15
//var biweeks = ee.List.sequence(120,270,15); 
//var start_date = '2019-04-30';
//var end_date = '2019-10-15';

/////  Setting for UY:
//it will merge two sequences (biweeks + biweeks2) when working with the summer season of the South Hemisphere
//growing season in UY single crop: Nov15 to Apr30:
var biweeks = ee.List.sequence(320,350,15); 
var biweeks2 = ee.List.sequence(0,105,15); 

//growing season in UY doble crop: Dec15 to May30:
//var biweeks = ee.List.sequence(350,350,15); 
//var biweeks2 = ee.List.sequence(0,105,15); //135
//var start_date = '2018-12-15';
//var end_date = '2019-04-30';

print(biweeks,biweeks2)

// Vegetation Index list to download:'q_mean','theta_mean','H_mean'
var index = 'm_mean';

//var cloud_threshold = 10;

// Select S1 GRD ImageCollection
var s1 = ee.ImageCollection('COPERNICUS/S1_GRD_FLOAT')
.filter(ee.Filter.eq('instrumentMode', 'IW'))
.filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
.filter(ee.Filter.eq('resolution_meters', 10))
.filterDate(Start_date, Finish_date)
//.filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'))
//.filter(ee.Filter.eq('platform_number', 'B'))
//.filter(ee.Filter.eq('relativeOrbitNumber_start', 141))
.filterBounds(geometry1);
print(s1)
var bandNames = s1.first().bandNames().remove('angle');

///////////////// Compute Indices from S1 ///////////////////////////
// indexes from : https://doi.org/10.1016/j.isprsjprs.2021.05.013
var ratio_param = function(image) {
  var ratio = image.expression(
    'VH/VV',
    {'VH': image.select('VH'),
    'VV' : image.select('VV')}).rename('q')
    return image.addBands(ratio)
}

var purity_param = function(image) {
  var purity = image.expression(
    '(1-q)/(1+q)',
    {'q': image.select('q')}).rename('m')
    return image.addBands(purity)
}

var theta_param = function(image) {
  var theta = image.expression(
    'atan(pow((1-q),2)/(1+pow(q,2)-q))',
    {'q': image.select('q')}).rename('theta')
    return image.addBands(theta)
}


var entropy_param = function(image) {
  var p1 = image.expression(
    '1/(1+q)',
    {'q': image.select('q')})
    
  var p2 = image.expression(
    'q/(1+q)',
    {'q': image.select('q')})
    
  var h1 = p1.expression(
    'p1*(log(p1)/log(2))',
    {'p1':p1})
    
  var h2 = p2.expression(
    'p2*(log(p2)/log(2))',
    {'p2':p2})
    
  var H = image.expression(
    '-1*(h1+h2)',
    {'h1':h1, 'h2':h2}).rename('H')
    return image.addBands(H)
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
// map indices
var s1_params = s1.select(bandNames).map(ratio_param).map(purity_param).map(theta_param).map(entropy_param);
print('s1_params',s1_params);
var img = s1_params.first();

var n = s1_params.size()
var colList = s1_params.toList(n);
print('colList',colList.get(8))

/*
//with replaced: to be used in case of missing S1 images
var currentDate = ee.Date('2020-08-20');
var replaced_image = s1_params.filterDate(currentDate.advance(-12, 'day'), currentDate.advance(12, 'day')).mean();
var replaced_image = replaced_image.set({'system:time_start': currentDate.millis()})
//var new_coll = ee.ImageCollection([s1_params,replaced_image])
print('replaced_image',replaced_image)
var col2 = ee.ImageCollection([replaced_image])
print('col2',col2)
var s1_params = s1_params.merge(col2).sort('system:index')
print('coll',s1_params)

var replaced = s1_params.map(function(image){
//var currentDate = ee.Date(image.get('system:time_start'));
var currentDate = ee.Date('2020-08-20');
var meanVal = s1_params.filterDate(
currentDate.advance(-12, 'day'), currentDate.advance(12, 'day')).mean();
return meanVal.where(image, image)
});
print('replaced',replaced)
var col =replaced.select('q').toBands();
print('col',col)
//.set({'system:time_start': image.get('segmentStartTime')})
*/

////////  Monthly means  ////////

var Biweekly_mean = ee.ImageCollection.fromImages(
biweeks.map(function(m) {
m = ee.Number(m)
return s1_params.filter(ee.Filter.calendarRange(m, m.add(15), 'day_of_year'))
.filterBounds(geometry1)
.reduce(ee.Reducer.mean())
.set('day_of_year', m);
}));
print('Monthly_mean',Biweekly_mean);

var Biweekly_mean2 = ee.ImageCollection.fromImages(
biweeks2.map(function(m) {
m = ee.Number(m)
return s1_params.filter(ee.Filter.calendarRange(m, m.add(15), 'day_of_year'))
.filterBounds(geometry1)
.reduce(ee.Reducer.mean())
.set('day_of_year', m);
}));
print('Biweekly_mean2',Biweekly_mean2);

var Biweekly_mean = Biweekly_mean.merge(Biweekly_mean2).sort('system:index');
print('Biweekly_mean FIANL',Biweekly_mean);


var n = Biweekly_mean.size();
print(n)
var colList = Biweekly_mean.toList(n);
print('collist',colList)

//choose 1 img from colList with get
var img = ee.Image(colList.get(2));
//print(img)
Map.addLayer(img.clip(geometry1), {min: 0, max: 1}, 'sub image')

//clip collection by roi
var clip_roi = Biweekly_mean
.map(function(image){return image.clip(geometry1)

});
print('clip_roi',clip_roi);

//////////////////// resampling S1 to S2 ///////////////////////////////////
var s2=ee.ImageCollection('COPERNICUS/S2_SR')
//.filter(ee.Filter.inList('MGRS_TILE',tiles))
.filterDate(Start_date, Finish_date)
.filterBounds(geometry1);
print('collection S2', s2);

//Loop to choose one images from s2 and get Projection from this band
var n = s2.size().getInfo();
var colList = s2.toList(n);
for (var i = 0; i < n; i++) {
var img = ee.Image(colList.get(1));
var img = img.select('SCL')
}
var Projection = img.projection()
print('proj_S2',Projection)

var img_repo = function(image){
var s1_res = image
// Request the data at the scale and projection of the s2 image.
.reproject({
crs: Projection
})
return s1_res
}

var reproj_s1 = clip_roi.map(img_repo)
print('test',reproj_s1)

////////////////////////////////////////////////////////////////////////////////////
// function to merge bands in multi-image
var mergeBands = function mergeBands(image, previous) {
return ee.Image(previous).addBands(image);
};

// select index and stack images to download
var download_roi = reproj_s1
.map(function(image) {return image.select(index)
})
.iterate(mergeBands, ee.Image([]));
print('download_roi',download_roi);

//////////////////////////download //////////////////////////////////////////////////////

Export.image.toDrive({
image: download_roi,
description: index + '_' + Field,
crs: proj,
crsTransform:[20,0,199980,0,-20,4600020], //copy from projection of S2
scale: 20,
fileFormat: 'GeoTIFF',
folder: exp_fold,
region: geometry1,
maxPixels: 10e10});

//// Plot S1 

//clip collection by roi
var clip_geometry1 = s1_params
.map(function(image){return image.clip(geometry1) });

// select bands with VI
var vi_s1 = clip_geometry1
.map(function(image) {return image.select('m','q','theta')});
print('Index_S1',vi_s1)

//plots
var chart_s1 =
ui.Chart.image.series({imageCollection: vi_s1, region: geometry1, reducer: ee.Reducer.mean(), scale: 10})
//.setSeriesNames(['VH'])
.setOptions({
title: 'S1-GRD VI mean',
hAxis: {
title: 'Time',
titleTextStyle: {italic: false, bold: true},
},
vAxis:
{title: 'VI', viewWindow: {min: 0, max: 1},titleTextStyle: {italic: false, bold: true}},
colors: ['cf513e', '1d6b99', 'f0af07']
});
print(chart_s1);

var chart_roi =
ui.Chart.image.series({imageCollection: vi_s1, region: geometry1,reducer: ee.Reducer.stdDev(), scale: 10})
//.setSeriesNames(['VH'])
.setOptions({
title: 'S1-GRD VI std',
hAxis: {
title: 'Time',
titleTextStyle: {italic: false, bold: true},
},
vAxis:
{title: 'VI', viewWindow: {min: 0, max: 0.4},titleTextStyle: {italic: false, bold: true}},
colors: ['cf513e', '1d6b99', 'f0af07']
});
print(chart_roi);


Map.centerObject(geometry1, 12);
Map.addLayer(vi_s1,{bands:'q', min:0,max:1}, 'q');
Map.addLayer(vi_s1,{bands:'m', min:0,max:1}, 'm');
Map.addLayer(vi_s1,{bands:'theta', min:0,max:1}, 'theta');
